#!/usr/bin/env python3
from analysis_tools_cython import *
from multiprocessing import Pool,Lock
from itertools import product
from functools import partial
import os
import sys
import traceback
import argparse


parser = argparse.ArgumentParser(description='Analyse lightcurves in target directory.')
parser.add_argument('-d', help='target directory(s)',
                        default='.',nargs='+',dest='path')

parser.add_argument('-t', help='number of threads to use',default=1,
                        dest='threads',type=int)

# Still writing to stdout, no inbuilt logging yet.
#parser.add_argument('-o',default='output.txt',dest='of',help='output file')


def process_file(path,f):
    f_path = os.path.join(path,f)
    table = import_lightcurve(f_path)
    t,flux = clean_data(table)
    flux = normalise_flux(flux)
    flux = fourier_filter(flux,8)
    T = test_statistic_array(flux,60)
    params = double_gaussian_curve_fit(T)
    ratio,separation = interpret(params)

#    lock.acquire()
    print(f,ratio,separation,T.min())
#    lock.release()


# Get directories from command line arguments.
args = parser.parse_args()

paths = []
for path in args.path:
    paths.append( os.path.expanduser(path) )


## Prepare multithreading.
#  Locking not working for now, probably fine...
#
#def init(l):
#    global lock
#    lock=l
#
#l = Lock()
pool = Pool(processes=args.threads)#, initializer=init, initargs=(l,))


for path in paths:
    if not os.path.isdir(path):
        print(path,'not a directory, skipping.',file=sys.stderr)
        continue

    fits_files = [f for f in os.listdir(path) if f.endswith('.fits')]

    process_in_dir = partial(process_file,path)
    try:
        pool.map(process_in_dir,fits_files)
    except (KeyboardInterrupt, SystemExit):
        raise
    except Exception as e:
        traceback.print_exc()
        print("\n",file=sys.stderr)

