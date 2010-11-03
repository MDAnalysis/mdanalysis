# $Id$
"""Setuptools based setup script for MDAnalysis.

This uses setuptools <http://pypi.python.org/pypi/setuptools>.

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
if sys.version_info[:2] < (2, 4):
    print "MDAnalysis requires Python 2.4 or better.  Python %d.%d detected" % \
        sys.version_info[:2]
    print "Please upgrade your version of python."
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

from Pyrex.Distutils import build_ext
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
    RELEASE = "0.7.1-devel"
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
                  Extension('core.rms_fitting', ['src/numtools/rms_fitting.pyx'],
                            libraries = ['m'],
                            define_macros=define_macros,
                            include_dirs = include_dirs+fast_numeric_include,
                            extra_link_args=fast_numeric_link,
                            extra_compile_args=extra_compile_args),
                  Extension('core._transformations', ['src/transformations/transformations.c'],
                            libraries = ['m'],
                            define_macros=define_macros,
                            include_dirs = include_dirs+fast_numeric_include,
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
                                    'tests/data/*.crd', 'tests/data/*.xyz', 'tests/data/*bz2'],
                               },
          ext_package       = 'MDAnalysis',
          ext_modules       = extensions,
          classifiers       = CLASSIFIERS,
          long_description  = LONG_DESCRIPTION,
          cmdclass = {'build_ext': build_ext},
          install_requires = ['numpy>=1.0',  # currently not useful because without numpy we don't get here
                              'biopython',   # required for standard PDB reader
                              ],
          extras_require = {
                'tests': ['nose>=0.10'],
                'analysis': ['networkx>=1.0',  # LeafletFinder
                             'scipy',          # sparse contact matrix
                             ],
                },          
          zip_safe = False,     # as a zipped egg the *.so files are not found (at least in Ubuntu/Linux)
          )
