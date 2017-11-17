#!/usr/bin/env python3
from analysis_tools_cython import *
import os
import sys
import traceback

if len(sys.argv) > 0:
    directory = sys.argv[1]
else:
    print("Usage:",sys.argv[0],"[directory]",file=sys.stderr)

if not os.path.isdir(directory):
    print("Directory not found.",file=sys.stderr)
    sys.exit()

if not directory.endswith("/"):
    directory += "/"

for fits_file in os.listdir(directory):
    if fits_file.endswith(".fits"):
        try:
            table = import_lightcurve(directory+fits_file)
            t,flux = clean_data(table)
            flux = normalise_flux(flux)
            flux = fourier_filter(flux,8)
            T = test_statistic_array(flux,60)
            params = double_gaussian_curve_fit(T)
            ratio,separation = interpret(params)
            print(fits_file,ratio,separation,T.min())
        except (KeyboardInterrupt, SystemExit):
            raise
        except Exception as e:
            print(fits_file,file=sys.stderr)
            traceback.print_exc()
            print("/n",file=sys.stderr)

