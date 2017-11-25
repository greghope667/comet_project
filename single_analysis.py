#!/usr/bin/env python3
from analysis_tools_cython import *
from scipy.stats import skew
import sys
import os
import matplotlib.pyplot as plt

if len(sys.argv) > 1:
    fits_file = sys.argv[1]
    table=import_lightcurve(fits_file)
else:
    print("Missing argument.")
    print("Usage:",sys.argv[0],"[FILENAME]")
    sys.exit()


timestep = calculate_timestep(table)
t,flux = clean_data(table)
N = len(t)
ones = np.ones(N)

flux = normalise_flux(flux)
filteredflux = fourier_filter(flux,8)
A_mag = np.abs(np.fft.rfft(flux))
periodicnoise = flux-filteredflux
sigma = flux.std()

T = test_statistic_array(filteredflux,60)
data = nonzero(T)

# Find minimum test statistic value, and its location.
minT_pos = np.unravel_index(T.argmin(),T.shape)
minT = T[minT_pos]
minT_time = t[minT_pos[1]]
minT_duration = 2*minT_pos[0]*timestep
print("Maximum transit chance:")
print("   Time =",round(minT_time,2),"days.")
print("   Duration =",round(minT_duration,2),"days.")
print("   T =",round(minT,1))
print("Skew =",round(skew(data),2))
print("Transit depth =",round(1+flux[minT_pos[1]],6))


fig1,axarr = plt.subplots(4)
axarr[0].plot(A_mag)
axarr[1].plot(t,flux+ones,t,periodicnoise+ones)
axarr[2].plot(t,filteredflux+ones)
cax = axarr[3].imshow(T)
axarr[3].set_aspect('auto')
fig1.colorbar(cax)

params = double_gaussian_curve_fit(T)
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
T_test_nonzero = np.array(data)
_,bins,_ = ax2.hist(T_test_nonzero,bins=100,log=True)
y = np.maximum(bimodal(bins,*params),10)
ax2.plot(bins,y)

plt.show()

