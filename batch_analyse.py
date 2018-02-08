#!/usr/bin/env python3
from analysis_tools_cython import *
from functools import partial
#from scipy.stats import skew
import multiprocessing
import os
import sys
import traceback
import argparse


parser = argparse.ArgumentParser(description='Analyse lightcurves in target directory.')
parser.add_argument(help='target directory(s)',
                        default='.',nargs='+',dest='path')

parser.add_argument('-t', help='number of threads to use',default=1,
                        dest='threads',type=int)

parser.add_argument('-o',default='output.txt',dest='of',help='output file')


def process_file(lock,of,path,f):

    try:
        f_path = os.path.join(path,f)
        table = import_lightcurve(f_path)
        t,flux,real = clean_data(table)
        flux = normalise_flux(flux)
        lombscargle_filter(t,flux,real,0.05)
        flux = flux*real
        T = test_statistic_array(flux,60)

        Ts = nonzero(T).std()
        m,n = np.unravel_index(T.argmin(),T.shape)
        Tm = T[m,n]
        Tm_time = t[n]
        Tm_duration = m*calculate_timestep(table)
        Tm_start = n-math.floor((m-1)/2)
        Tm_end = Tm_start + m
        Tm_depth = flux[Tm_start:Tm_end].mean()

        s = classify(m,n,real)
        asym = calc_asymmetry(m,n,t,flux)

#        result_str = ' '.join([f, str(Tm), str(Tm/Ts),str(asym),
#                            str(Tm_time),str(Tm_duration),str(Tm_depth),s])

        result_str = f+' '+\
            ' '.join([str(round(a,8)) for a in 
                [Tm, Tm/Ts, asym,Tm_time,Tm_duration,Tm_depth]])+' '+s

        lock.acquire()
        with open(of,'a') as out_file:
            out_file.write(result_str+'\n')
        lock.release()
    except (KeyboardInterrupt, SystemExit):
        print("Process terminated early, exiting",file=sys.stderr)
        raise
    except Exception as e:
        traceback.print_exc()
        print("\n",file=sys.stderr)


# Get directories from command line arguments.
args = parser.parse_args()

paths = []
for path in args.path:
    paths.append( os.path.expanduser(path) )


## Prepare multithreading.
pool = multiprocessing.Pool(processes=args.threads)
m = multiprocessing.Manager()
l = m.Lock()


for path in paths:
    if not os.path.isdir(path):
        print(path,'not a directory, skipping.',file=sys.stderr)
        continue

    fits_files = [f for f in os.listdir(path) if f.endswith('.fits')]

    process_partial = partial(process_file,l,args.of,path)
    
    pool.map(process_partial,fits_files)

