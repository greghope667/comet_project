from astropy.io import fits
from astropy.table import Table
from scipy.optimize import curve_fit
import numpy as np
cimport numpy as np
import math
##import matplotlib.pyplot as plt



def import_lightcurve(file_path):
    """Returns (N by 2) table, columns are (time, flux)."""

    try:
        hdulist = fits.open(file_path)
    except FileNotFoundError:
        print("Import failed: file not found")

        return (0, 0)

    scidata = hdulist[1].data
    table = Table(scidata)['TIME','PDCSAP_FLUX']

    # Delete rows containing NaN values.
    nan_rows = [ i for i in range(len(table)) if math.isnan(table[i][1]) or math.isnan(table[i][0]) ]

    table.remove_rows(nan_rows)
    
    return table


def calculate_timestep(table):
    """Returns median value of time differences between data points,
    estimate of time delta data points."""

    dt = [ table[i+1][0] - table[i][0] for i in range(len(table)-1) ]
    dt.sort()
    return dt[int(len(dt)/2)]


def clean_data(table):
    """Interpolates missing data points, so we have equal time gaps 
    between points. Returns two numpy arrays, time, flux"""

    t = []
    x = []
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

        t.append(ti)
        x.append(xi)

    return np.array(t),np.array(x)


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

    B = np.zeros(len(A))
    for i in freq_index:
        B[i] = A[i]
    
    # Fitted flux is our periodic approximation to the flux
    fitted_flux = np.fft.irfft(B,len(flux))

    return flux - fitted_flux


def test_statistic_array(np.ndarray[np.float64_t,ndim=1]  flux, int max_half_width):
    cdef int N = flux.shape[0]
    cdef int n = max_half_width

    cdef int i, m, j
    cdef float mu,sigma,norm_factor
    sigma = flux.std()

    cdef np.ndarray[dtype=np.float64_t,ndim=2] t_test = np.zeros([n,N])
#    cdef np.ndarray[dtype=np.float64_t,ndim=1] flux_points = np.zeros(2*n)
    for m in range(1,n): 

        norm_factor = 1 / ((2*m)**0.5 * sigma)

        mu = flux[0:(2*m)].sum()
        t_test[m][m] = mu * norm_factor

        for i in range(m+1,N-m):

            ##mu = flux[(i-m):(i+m)].sum()
            mu += (flux[i+m-1] - flux[i-m-1])
            t_test[m][i] = mu * norm_factor

    return t_test


def gauss(x,A,mu,sigma):
    return abs(A)*np.exp( -(x - mu)**2 / (2 * sigma**2) )

def bimodal(x,A1,mu1,sigma1,A2,mu2,sigma2):
    return gauss(x,A1,mu1,sigma1)+gauss(x,A2,mu2,sigma2)

def double_gaussian_curve_fit(T):
    """Fit two normal distributions to a test statistic vector T.
    Returns (A1,mu1,sigma1,A2,mu2,sigma2)"""

    data = np.array([i for i in T.flat if i != 0.0])
    N = len(data)

    T_min = data.min()
    T_max = data.max()

    # Split data into 100 bins, so we can approximate pdf.
    bins = np.linspace(T_min,T_max,101)
    y,bins = np.histogram(data,bins)
    x = (bins[1:] + bins[:-1])/2

    # Initial guess at parameters
    A1 = N/2
    A2 = N/2

    mu1 = (2*T_min + T_max)/3
    mu2 = (T_min + 2*T_max)/3

    sigma1 = 1
    sigma2 = 1

    start_params = (A1,mu1,sigma1,A2,mu2,sigma2)
    params,cov=curve_fit(bimodal,x,y,start_params)
    
    for i in [0,2,3,5]:
        params[i] = abs(params[i])

    return params


def interpret(params):
    # Choose A1,mu1,sigma1 to be stats for larger peak
    if params[0]>params[3]:
        A1,mu1,sigma1,A2,mu2,sigma2 = params
    else:
        A2,mu2,sigma2,A1,mu1,sigma1 = params

    height_ratio = A2/A1
    separation = (mu1 - mu2)/sigma1

    return height_ratio,separation

    



