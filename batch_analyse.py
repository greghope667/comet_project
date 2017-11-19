#!/usr/bin/env python3
from analysis_tools_cython import *
import os
import sys
import traceback
import argparse


parser = argparse.ArgumentParser(description='Analyse lightcurves in target directory.')
parser.add_argument('-d', help='target directory(s)',
                        default='.',nargs='+',dest='path')

# Still writing to stdout, no inbuilt logging yet.
#parser.add_argument('-o',default='output.txt',dest='of',help='output file')

args = parser.parse_args()

paths = []
for path in args.path:
    paths.append( os.path.expanduser(path) )


for path in paths:
    if not os.path.isdir(path):
        print(path,'not a directory, skipping.',file=sys.stderr)
        continue

    for fits_file in os.listdir(path):
        if not fits_file.endswith(".fits"):
            continue
        try:
            f = os.path.join(path,fits_file)
            table = import_lightcurve(f)
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
            print("\n",file=sys.stderr)

