#!/usr/bin/env python3
from distutils.core import setup
from Cython.Build import cythonize

setup(
        name = 'light_curve_analysis',
        ext_modules = cythonize("analysis_tools_cython.pyx"),
)

