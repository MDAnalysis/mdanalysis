# $Id$
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


import ConfigParser

include_dirs = [numpy_include]

if sys.platform == "darwin": # Mac OS X
    fast_numeric_include = ['/System/Library/Frameworks/vecLib.framework/Versions/A/Headers']
    fast_numeric_link = ["-framework","vecLib"]
elif sys.platform[:5] == "linux":
    parser = ConfigParser.ConfigParser()
    parser.read("setup.cfg")
    try:
        fast_numeric_include = parser.get("linux","fast_numeric_include").split()
        linkpath = ["-L"+path for path in parser.get("linux","fast_numeric_linkpath").split()]
        linklibs = ["-l"+lib for lib in parser.get("linux","fast_numeric_libs").split()]
        fast_numeric_link = linkpath + linklibs
    except ConfigParser.NoSectionError:
        fast_numeric_include = []
        fast_numeric_link = ["-llapack"]
else:
    fast_numeric_include = []
    fast_numeric_link = ["-llapack"]

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
                  Extension('coordinates.dcdtimeseries', ['src/dcd/dcdtimeseries.c'],
                            include_dirs = include_dirs+['src/dcd/include'],
                            define_macros=define_macros,
                            extra_compile_args=extra_compile_args),
                  Extension('core.distances', ['src/numtools/distances.c'],
                            include_dirs = include_dirs+['src/numtools'],
                            libraries = ['m'],
                            define_macros=define_macros,
                            extra_compile_args=extra_compile_args),
                  Extension('core.rms_fitting', ['src/numtools/rms_fitting.c'],
                            libraries = ['m'],
                            define_macros=define_macros,
                            include_dirs = include_dirs+fast_numeric_include,
                            extra_link_args=fast_numeric_link,
                            extra_compile_args=extra_compile_args),
                  Extension('core.qcprot', ['src/pyqcprot/pyqcprot.c'],
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
                                    'tests/data/*.pqr', 'tests/data/*.pdbqt', 'tests/data/*.bz2',
                                    ],
                               },
          ext_package       = 'MDAnalysis',
          ext_modules       = extensions,
          classifiers       = CLASSIFIERS,
          long_description  = LONG_DESCRIPTION,
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
