#! /usr/bin/python
"""Setuptools-based setup script for MDAnalysis.

A working installation of NumPy <http://numpy.scipy.org> is required.

For a basic installation just type the command::

  python setup.py install

For more in-depth instructions, see the installation section at the
MDAnalysis Wiki:

  http://code.google.com/p/mdanalysis/wiki/Install

Also free to ask on the MDAnalysis mailing list for help:

  http://groups.google.com/group/mdnalysis-discussion

(Note that the group really is called `mdnalysis-discussion' because
Google groups forbids any name that contains the string `anal'.)

By default we use setuptools <http://pypi.python.org/pypi/setuptools>.  The
details of such an "EasyInstall" installation procedure are shown on

  http://peak.telecommunity.com/DevCenter/EasyInstall

By changing the code below you can also switch to a standard distutils
installation.
"""
from __future__ import with_statement

#------------------------------------------------------------
# selection of the installation system
#------------------------------------------------------------
#
# Standard distutils-based installation:
#
##from distutils.core import setup, Extension

# setuptools ("EasyInstall") installation:
#
# If you want EasyInstall features then enable the next three lines and comment
# out the preceding line 'from distutils.core import ...'
#
from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup, Extension
#
#------------------------------------------------------------

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
    print "*** package 'numpy' not found ***"
    print "MDAnalysis requires a version of NumPy (>=1.0.3), even for setup."
    print "Please get it from http://numpy.scipy.org/ or install it through your package manager."
    sys.exit(-1)

try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

include_dirs = [numpy_include]

# Handle cython modules
try:
    from Cython.Distutils import build_ext
    use_cython = True
    cmdclass = {'build_ext': build_ext}
except ImportError:
    use_cython = False
    cmdclass = {}

if __name__ == '__main__':
    RELEASE = "0.7.5-devel"
    with open("SUMMARY.txt") as summary:
        LONG_DESCRIPTION = summary.read()
    CLASSIFIERS = ['Development Status :: 4 - Beta',
                   'Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: GNU General Public License (GPL)',
                   'Operating System :: POSIX',
                   'Operating System :: MacOS :: MacOS X',
                   'Programming Language :: Python',
                   'Programming Language :: C',
                   'Topic :: Scientific/Engineering :: Bio-Informatics',
                   'Topic :: Scientific/Engineering :: Chemistry',
                   'Topic :: Software Development :: Libraries :: Python Modules',
                   ]

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
                  Extension('coordinates.dcdtimeseries', ['src/dcd/dcdtimeseries.%s' % ("pyx" if use_cython else "c")],
                            include_dirs = include_dirs+['src/dcd/include'],
                            define_macros=define_macros,
                            extra_compile_args=extra_compile_args),
                  Extension('core.distances', ['src/numtools/distances.%s' % ("pyx" if use_cython else "c")],
                            include_dirs = include_dirs+['src/numtools'],
                            libraries = ['m'],
                            define_macros=define_macros,
                            extra_compile_args=extra_compile_args),
                  Extension('core.qcprot', ['src/pyqcprot/pyqcprot.%s' % ("pyx" if use_cython else "c")],
                            include_dirs=include_dirs,
                            extra_compile_args=["-O3","-ffast-math"]),
                  Extension('core._transformations', ['src/transformations/transformations.c'],
                            libraries = ['m'],
                            define_macros=define_macros,
                            include_dirs = include_dirs,
                            extra_compile_args=extra_compile_args),
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
          description       = 'An object-oriented toolkit to analyze molecular dynamics trajectories generated by CHARMM, Gromacs, NAMD, LAMMPS, or Amber.',
          author            = 'Naveen Michaud-Agrawal',
          author_email      = 'naveen.michaudagrawal@gmail.com',
          url               = 'http://mdanalysis.googlecode.com/',
          requires          = ['numpy (>=1.0.3)', 'biopython',
                               'networkx (>=1.0)', 'scipy',
                               'GridDataFormats'],
          provides          = ['MDAnalysis'],
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
          ext_package       = 'MDAnalysis',
          ext_modules       = extensions,
          classifiers       = CLASSIFIERS,
          long_description  = LONG_DESCRIPTION,
          cmdclass          = cmdclass,
          # all standard requirements are available through PyPi and
          # typically can be installed without difficulties through setuptools
          install_requires = ['numpy>=1.0.3',   # currently not useful because without numpy we don't get here
                              'biopython',      # required for standard PDB reader
                              'networkx>=1.0',  # LeafletFinder
                              'GridDataFormats>=0.2.2', # volumes and densities
                              ],
          # extras can be difficult to install through setuptools and/or
          # you might prefer to use the version available through your
          # packaging system
          extras_require = {
                'analysis': ['matplotlib',
                             'scipy',          # sparse contact matrix
                             ],
                },
          test_suite = "nose.collector",
          tests_require = ['nose>=0.10',
                           'MDAnalysisTests==0.7.5',
                           ],
          zip_safe = False,     # as a zipped egg the *.so files are not found (at least in Ubuntu/Linux)
          )
