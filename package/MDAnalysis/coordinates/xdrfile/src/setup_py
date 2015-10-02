# Python bindings for libxdrfile2
# Copyright (c) 2010 Oliver Beckstein <orbeckst@gmail.com>
# Copyright (c) 2013,2014 Manuel Melo <manuel.nuno.melo@gmail.com>
# Published under the GNU GENERAL PUBLIC LICENSE Version 2 (or higher)
"""Distutils setup for the stand-alone xdrfile package

Copy this file to the parent directory of xdrfile and rename it to
setup.py. See xdrfile/README for instructions.

"""

# For setuptools (easy install) uncomment
##from ez_setup import use_setuptools
##use_setuptools()
##from setuptools import setup, Extension

# For basic distutils use
from distutils.core import setup, Extension

import sys, os

# Obtain the numpy include directory.  This logic works across numpy versions.
import numpy
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

include_dirs = [numpy_include]


if __name__ == '__main__':
    setup(name              = 'xdrfile',
          version           = '0.1-libxdrfile2',
          description       = 'Python interface to libxdrfile2 library',
          author            = 'David van der Spoel, Erik Lindahl, Oliver Beckstein, Manuel Melo',
          author_email      = 'orbeckst@gmail.com',
          url               = 'http://www.MDAnalysis.org/',
          license           = 'GPL 2',
          packages          = ['xdrfile'],
          package_dir       = {'xdrfile': 'xdrfile'},
          ext_package       = 'xdrfile',
          ext_modules       = [Extension('_libxdrfile2',
                                         sources=['xdrfile/libxdrfile2_wrap.c',
                                                  'xdrfile/xdrfile.c', 	
                                                  'xdrfile/xdrfile_trr.c',
                                                  'xdrfile/xdrfile_xtc.c'],
                                         include_dirs = include_dirs),],
          )
