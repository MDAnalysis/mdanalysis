# $Id$
"""Setuptools based setup script for MDAnalysis (for developers)

This uses setuptools <http://pypi.python.org/pypi/setuptools> and
Cython <http://cython.org> to build all files from primary sources.

A working installation of NumPy <http://numpy.scipy.org> is required.

For the easiest installation just type the command:

  python setup.py install

The details of such an "EasyInstall" installation procedure are shown on

  http://peak.telecommunity.com/DevCenter/EasyInstall

For more in-depth instructions, see the installation section at the
MDAnalysis Wiki:

  http://code.google.com/p/mdanalysis/wiki/Install

Or, if all else fails, feel free to ask on the MDAnalysis mailing list
for help:

  http://groups.google.com/group/mdnalysis-discussion

(Note that the group really is called `mdnalysis-discussion' because
Google groups forbids any name that contains the string `anal'.)
"""
# EasyInstall installation:
from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup, Extension

# If you don't require EasyInstall features you can also comment out
# the previous three lines and use the following standard distutils-based
# installation:
###from distutils.core import setup, Extension

import sys, os
import glob

# Make sure I have the right Python version.
if sys.version_info[:2] < (2, 5):
    print "MDAnalysis requires Python 2.5 or better. Python %d.%d detected" % \
        sys.version_info[:2]
    print "Please upgrade your version of Python."
    sys.exit(-1)

try:
    # Obtain the numpy include directory.  This logic works across numpy versions.
    import numpy
except ImportError:
    # TODO: somehow fix this so that easy_install could get numpy if needed
    print "*** package 'numpy' not found ***"
    print "MDAnalysis requires a version of NumPy, even for setup."
    print "Please get it from http://numpy.scipy.org/ or install it through your package manager."
    sys.exit(-1)

try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

try:
    from Cython.Distutils import build_ext
except ImportError:
    print "*** package 'Cython' not found ***"
    print "MDAnalysis requires Cython at the setup and build stage."
    print "Please get it from http://cython.org/ or install it through your package manager."
    sys.exit(-1)

include_dirs = [numpy_include]


if __name__ == '__main__':
    RELEASE = "0.7.3-devel"
    LONG_DESCRIPTION = \
"""MDAnalysis is a tool for analyzing molecular dynamics trajectories.
"""
    CLASSIFIERS = ['Development Status :: 1 - Alpha',
                   'Environment :: Workstation',
                   'Intended Audience :: Scientists',
                   'License :: OSI Approved :: GPL License',
                   'Operating System :: OS Independent',
                   'Programming Language :: Python',
                   'Topic :: Scientific Software :: Biology',
                   'Topic :: Scientific Software :: Chemistry',]

    if 'DEBUG_CFLAGS' in os.environ:
        extra_compile_args = '\
            -std=c99 -pedantic -Wall -Wcast-align -Wcast-qual -Wpointer-arith \
            -Wchar-subscripts -Winline -Wnested-externs -Wbad-function-cast \
            -Wunreachable-code -Werror'
        define_macros = [('DEBUG', '1')]
    else:
        extra_compile_args = ''
        define_macros = []

    extensions = [Extension('coordinates._dcdmodule', ['src/dcd/dcd.c'],
                            include_dirs = include_dirs+['src/dcd/include'],
                            define_macros=define_macros,
                            extra_compile_args=extra_compile_args),
                  Extension('coordinates.dcdtimeseries', ['src/dcd/dcdtimeseries.pyx'],
                            include_dirs = include_dirs+['src/dcd/include'],
                            define_macros=define_macros,
                            extra_compile_args=extra_compile_args),
                  Extension('core.distances', ['src/numtools/distances.pyx'],
                            include_dirs = include_dirs+['src/numtools'],
                            libraries = ['m'],
                            define_macros=define_macros,
                            extra_compile_args=extra_compile_args),
                  Extension('core.qcprot', ['src/pyqcprot/pyqcprot.pyx'],
                            include_dirs=include_dirs,
                            extra_compile_args=["-O3","-ffast-math"]),
                  Extension('core._transformations', ['src/transformations/transformations.c'],
                            libraries = ['m'],
                            define_macros=define_macros,
                            include_dirs = include_dirs,
                            extra_compile_args=extra_compile_args),
                  #Extension('util.delaunay', ['src/delaunay/delaunay.pyx', 'src/delaunay/blas.c', 'src/delaunay/tess.c'],
                  #          libraries = ['m'],
                  #          define_macros=define_macros,
                  #          include_dirs = include_dirs+fast_numeric_include+['src/delaunay'],
                  #          extra_link_args=fast_numeric_link,
                  #          extra_compile_args=extra_compile_args),
                  Extension('KDTree._CKDTree',
                            ["src/KDTree/KDTree.cpp",
                             "src/KDTree/KDTree.swig.cpp"],
                            include_dirs = include_dirs,
                            libraries=["stdc++"],
                            language="c++"),
                  Extension('coordinates.xdrfile._libxdrfile',
                            sources=['src/xdrfile/libxdrfile_wrap.c',
                                     'src/xdrfile/xdrfile.c',
                                     'src/xdrfile/xdrfile_trr.c',
                                     'src/xdrfile/xdrfile_xtc.c'],
                            include_dirs = include_dirs),
                  ]

    setup(name              = 'MDAnalysis',
          version           = RELEASE,
          description       = 'Python tools to support analysis of trajectories',
          author            = 'Naveen Michaud-Agrawal',
          author_email      = 'naveen.michaudagrawal@gmail.com',
          url               = 'http://mdanalysis.googlecode.com/',
          license           = 'GPL 2',
          packages          = [ 'MDAnalysis', 'MDAnalysis.core', 'MDAnalysis.topology',
                                'MDAnalysis.selections',
                                'MDAnalysis.coordinates',
                                'MDAnalysis.coordinates.xdrfile',
                                'MDAnalysis.coordinates.pdb',
                                'MDAnalysis.util', 'MDAnalysis.KDTree',
                                'MDAnalysis.analysis',
                                'MDAnalysis.tests'],
          package_dir       = {'MDAnalysis': 'MDAnalysis'},
          package_data      = {'MDAnalysis':
                                   ['tests/data/*.psf','tests/data/*.dcd','tests/data/*.pdb',
                                    'tests/data/*.gro', 'tests/data/*.xtc','tests/data/*.trr',
                                    'tests/data/*.crd', 'tests/data/*.xyz',
                                    'tests/data/*.prmtop', 'tests/data/*.trj', 'tests/data/*.mdcrd',
                                    'tests/data/*.pqr', 'tests/data/*.bz2',
                                    ],
                               },
          ext_package       = 'MDAnalysis',
          ext_modules       = extensions,
          classifiers       = CLASSIFIERS,
          long_description  = LONG_DESCRIPTION,
          cmdclass = {'build_ext': build_ext},
          install_requires = ['numpy>=1.0.3',  # currently not useful because without numpy we don't get here
                              'biopython',   # required for standard PDB reader
                              ],
          extras_require = {
                'tests': ['nose>=0.10'],
                'analysis': ['networkx>=1.0',  # LeafletFinder
                             'scipy',          # sparse contact matrix
                             'GridDataFormats',# http://github.com/orbeckst/GridDataFormats
                             ],
                },
          zip_safe = False,     # as a zipped egg the *.so files are not found (at least in Ubuntu/Linux)
          )
