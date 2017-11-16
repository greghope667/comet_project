from analysis_tools_cython import *
import os,sys

if len(sys.argv) > 0:
    directory = sys.argv[1]
else:
    print("Usage:",sys.argv[0],"[directory]")

if not directory.endswith("/"):
    directory.join("/")
for fits_file in os.listdir(dir_path):
    if fits_file.endswith(".fits"):
        print(fits_file)
        table = import_lightcurve(dir_path+fits_file)
        t,flux = clean_data(table)
        flux = normalise_flux(flux)
        flux = fourier_filter(flux,8)
        T = test_statistic_array(flux,60)
        params = double_gaussian_curve_fit(T)
        print(interpret(params))

