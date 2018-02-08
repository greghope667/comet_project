from astropy.io import fits
from astropy.table import Table
from scipy.optimize import curve_fit
from gatspy.periodic import LombScargleFast
import numpy as np
cimport numpy as np
import math
import sys,os

model = LombScargleFast(fit_period=True, silence_warnings=True)

def import_lightcurve(file_path):
    """Returns (N by 2) table, columns are (time, flux)."""

    try:
        hdulist = fits.open(file_path)
    except FileNotFoundError:
        print("Import failed: file not found")
        return

    scidata = hdulist[1].data
    table = Table(scidata)['TIME','PDCSAP_FLUX']

    # Delete rows containing NaN values.
    nan_rows = [ i for i in range(len(table)) if
            math.isnan(table[i][1]) or math.isnan(table[i][0]) ]

    table.remove_rows(nan_rows)

    # Smooth data by deleting overly 'spikey' points.
    spikes = [ i for i in range(1,len(table)-1) if \
            abs(table[i][1] - 0.5*(table[i-1][1]+table[i+1][1])) \
            > 3*abs(table[i+1][1] - table[i-1][1])]

    for i in spikes:
        table[i][1] = 0.5*(table[i-1][1] + table[i+1][1])

    return table


def calculate_timestep(table):
    """Returns median value of time differences between data points,
    estimate of time delta data points."""

    dt = [ table[i+1][0] - table[i][0] for i in range(len(table)-1) ]
    dt.sort()
    return dt[int(len(dt)/2)]


def clean_data(table):
    """Interpolates missing data points, so we have equal time gaps
    between points. Returns three numpy arrays, time, flux, real.
    real is 0 if data point interpolated, 1 otherwise."""

    t = []
    x = []
    r = []
    timestep = calculate_timestep(table)

    for row in table:
        ti, xi = row

        if len(t) > 0:
            steps = int(round( (ti - t[-1])/timestep ))

            if steps > 1:
                fluxstep = (xi - x[-1])/steps

                for _ in range(steps-1):
                    t.append(timestep + t[-1])
                    x.append(fluxstep + x[-1])
                    r.append(0)

        t.append(ti)
        x.append(xi)
        r.append(1)

    return np.array(t),np.array(x),np.array(r)


def normalise_flux(flux):
    """Requires flux to be a numpy array.
    Normalisation is x --> (x/mean(x)) - 1"""

    return flux/flux.mean() - np.ones(len(flux))


def fourier_filter(flux,freq_count):
    """Attempt to remove periodic noise by finding and subtracting
    freq_count number of peaks in (discrete) fourier transform."""

    A = np.fft.rfft(flux)
    A_mag = np.abs(A)

    # Find frequencies with largest amplitudes.
    freq_index = np.argsort(-A_mag)[0:freq_count]

    # Mult by 1j so numpy knows we are using complex numbers
    B = np.zeros(len(A)) * 1j
    for i in freq_index:
        B[i] = A[i]

    # Fitted flux is our periodic approximation to the flux
    fitted_flux = np.fft.irfft(B,len(flux))

    return flux - fitted_flux


def lombscargle_filter(time,flux,real,min_score):
    """Also removes periodic noise, using lomb scargle methods."""
    time_real = time[real == 1]

    period = time[-1]-time[0]
    N = len(time)
    nyquist_period = (2*period)/N
    model.optimizer.period_range = (nyquist_period,period)
    model.optimizer.quiet = True

    try:
        for _ in range(30):
            flux_real = flux[real == 1]
            model.fit(time_real,flux_real)

            if model.score(model.best_period) < min_score:
                break

            flux -= model.predict(time)
    except:
        pass


def test_statistic_array(np.ndarray[np.float64_t,ndim=1] flux, int max_half_width):
    cdef int N = flux.shape[0]
    cdef int n = max_half_width

    cdef int i, m, j
    cdef float mu,sigma,norm_factor
    sigma = flux.std()

    cdef np.ndarray[dtype=np.float64_t,ndim=2] t_test = np.zeros([2*n,N])
#    cdef np.ndarray[dtype=np.float64_t,ndim=1] flux_points = np.zeros(2*n)
    for m in range(1,2*n):

        m1 = math.floor((m-1)/2)
        m2 = (m-1) - m1

        norm_factor = 1 / (m**0.5 * sigma)

        mu = flux[0:m].sum()
        t_test[m][m1] = mu * norm_factor

        for i in range(m1+1,N-m2-1):

            ##t_test[m][i] = flux[(i-m1):(i+m2+1)].sum() * norm_factor
            mu += (flux[i+m2] - flux[i-m1-1])
            t_test[m][i] = mu * norm_factor

    return t_test


def gauss(x,A,mu,sigma):
    return abs(A)*np.exp( -(x - mu)**2 / (2 * sigma**2) )

def bimodal(x,A1,mu1,sigma1,A2,mu2,sigma2):
    return gauss(x,A1,mu1,sigma1)+gauss(x,A2,mu2,sigma2)

def skewed_gauss(x,A,mu,sigma1,sigma2):
    y = np.zeros(len(x))
    for i in range(len(x)):
        if x[i] < mu:
            y[i] = gauss(x[i],A,mu,sigma1)
        else:
            y[i] = gauss(x[i],A,mu,sigma2)
    return y


def comet_curve(x,A,mu,sigma,tail):
    y = np.zeros(len(x))
    for i in range(len(x)):
        if x[i] < mu:
            y[i] = gauss(x[i],A,mu,sigma)
        else:
            y[i] = A*math.exp(-abs(x[i]-mu)/tail)
    return y


def single_gaussian_curve_fit(x,y):
    # Initial parameters guess
    i = np.argmax(y)
    A0 = y[i]
    mu0 = x[i]
    sigma0 = 1

    params_bounds = [[0,x[0],0], [np.inf,x[-1],np.inf]]
    params,cov = curve_fit(gauss,x,y,[A0,mu0,sigma0],bounds=params_bounds)
    return params


def nonzero(T):
    """Returns a 1d array of the nonzero elements of the array T"""
    return np.array([i for i in T.flat if i != 0])


def skewed_gaussian_curve_fit(x,y):
    # Initial parameters guess
    i = np.argmax(y)
    A0 = y[i]
    mu0 = x[i]
    sigma0 = 1

    params_bounds = [[0,x[0],0,0], [np.inf,x[-1],np.inf,np.inf]]
    params,cov = curve_fit(skewed_gauss,x,y,[A0,mu0,sigma0,sigma0],
                           bounds=params_bounds)
    return params


def comet_curve_fit(x,y):
    # Initial parameters guess
    i = np.argmax(y)
    A0 = y[i]
    mu0 = x[i]
    sigma0 = 1
    tail0 = 1

    width = (x[-1]-x[0])

    params_bounds = [[0,x[0],0,0], [np.inf,x[-1],width,width]]
    params,cov = curve_fit(comet_curve,x,y,[A0,mu0,sigma0,tail0],
                           bounds=params_bounds)
    return params


def double_gaussian_curve_fit(T):
    """Fit two normal distributions to a test statistic vector T.
    Returns (A1,mu1,sigma1,A2,mu2,sigma2)"""

    data = nonzero(T)
    N = len(data)

    T_min = data.min()
    T_max = data.max()

    # Split data into 100 bins, so we can approximate pdf.
    bins = np.linspace(T_min,T_max,101)
    y,bins = np.histogram(data,bins)
    x = (bins[1:] + bins[:-1])/2


    # We fit the two gaussians one by one, as this is more
    #  sensitive to small outlying bumps.
    params1 = single_gaussian_curve_fit(x,y)
    y1_fit = np.maximum(gauss(x,*params1),1)

    y2 = y/y1_fit
    params2 = single_gaussian_curve_fit(x,y2)

    params = [*params1,*params2]

    return params


def score_fit(y,fit):
    return sum(((y[i]-fit[i])**2 for i in range(len(y))))


def interpret(params):
    # Choose A1,mu1,sigma1 to be stats for larger peak
    if params[0]>params[3]:
        A1,mu1,sigma1,A2,mu2,sigma2 = params
    else:
        A2,mu2,sigma2,A1,mu1,sigma1 = params

    height_ratio = A2/A1
    separation = (mu2 - mu1)/sigma1

    return height_ratio,separation


def classify(m,n,real):
    N = len(real)
    if n-2*m < 0 or n+2*m > N:
        return "end"
    elif m < 3:
        return "point"
    elif real[(n-2*m):(n-m)].sum() < 0.5*m:
        return "artefact"
    else:
        return "maybeTransit"


def calc_asymmetry(m,n,time,flux):
    """Fit both symmetric and comet-like transit profiles and compare fit.

    Returns ratio of (errors squared)
    Possible errors and return values:
    -1 : Divide by zero as comet profile is exact fit
    -2 : Too close to end of light curve to fit profile
    -3 : Unable to fit model (e.g. timeout)
    """
    if n-3*m >= 0 and n+3*m < len(time):
        t = time[n-3*m:n+3*m]
        x = flux[n-3*m:n+3*m]
        background_level = (sum(x[:m]) + sum(x[5*m:]))/(2*m)
        x -= background_level

        try:
            params1 = single_gaussian_curve_fit(t,-x)
            params2 = comet_curve_fit(t,-x)
        except:
            return -3

        fit1 = -gauss(t,*params1)
        fit2 = -comet_curve(t,*params2)

        scores = [score_fit(x,fit) for fit in [fit1,fit2]]
        if scores[1] > 0:
            return scores[0]/scores[1]
        else:
            return -1
    else:
        return -2

