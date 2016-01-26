#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org

# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
# Beckstein and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""Setuptools-based setup script for MDAnalysis.

A working installation of NumPy <http://numpy.scipy.org> is required.

For a basic installation just type the command::

  python setup.py install

For more in-depth instructions, see the installation section at the
MDAnalysis Wiki:

  https://github.com/MDAnalysis/mdanalysis/wiki/INSTALL

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

from __future__ import print_function
from setuptools import setup, Extension, find_packages
from distutils.ccompiler import new_compiler
import os
import sys
import shutil
import tempfile

# Make sure I have the right Python version.
if sys.version_info[:2] < (2, 7):
    print('MDAnalysis requires Python 2.7 or better. Python %d.%d detected' %
          sys.version_info[:2])
    print('Please upgrade your version of Python.')
    sys.exit(-1)

if sys.version_info[0] < 3:
    import ConfigParser as configparser
    open_kwargs = {}
else:
    import configparser
    open_kwargs = {'encoding': 'utf-8'}

try:
    # Obtain the numpy include directory. This logic works across numpy
    # versions.
    import numpy as np
except ImportError:
    print('*** package "numpy" not found ***')
    print('MDAnalysis requires a version of NumPy (>=1.5.0), even for setup.')
    print('Please get it from http://numpy.scipy.org/ or install it through '
          'your package manager.')
    sys.exit(-1)

# Handle cython modules
try:
    from Cython.Distutils import build_ext
    cython_found = True
    cmdclass = {'build_ext': build_ext}
except ImportError:
    cython_found = False
    cmdclass = {}

if cython_found:
    # cython has to be >=0.16 to support cython.parallel
    import Cython
    from distutils.version import LooseVersion

    required_version = "0.16"

    if not LooseVersion(Cython.__version__) >= LooseVersion(required_version):
        raise ImportError(
            "Cython version {0} (found {1}) is required because it offers "
            "a handy parallelisation module".format(
                required_version, Cython.__version__))
    del Cython
    del LooseVersion


class Config(object):
    """Config wrapper class to get build options

    This class looks for options in the environment variables and the
    'setup.cfg' file. The order how we look for an option is.

    1. Environment Variable
    2. set in 'setup.cfg'
    3. given default

    Environment variables should start with 'MDA_' and be all uppercase.

    """

    def __init__(self, fname='setup.cfg'):
        if os.path.exists(fname):
            self.config = configparser.SafeConfigParser()
            self.config.read(fname)

    def get(self, option_name, default=None):
        environ_name = 'MDA_' + option_name.upper()
        if environ_name in os.environ:
            return os.environ[environ_name]
        try:
            option = self.config.get('options', option_name)
            return option
        except:
            return default


def get_numpy_include():
    try:
        numpy_include = np.get_include()
    except AttributeError:
        numpy_include = np.get_numpy_include()
    return numpy_include


def hasfunction(cc, funcname, include=None, extra_postargs=None):
    # From http://stackoverflow.com/questions/
    #            7018879/disabling-output-when-compiling-with-distutils
    tmpdir = tempfile.mkdtemp(prefix='hasfunction-')
    devnull = oldstderr = None
    try:
        try:
            fname = os.path.join(tmpdir, 'funcname.c')
            with open(fname, 'w') as f:
                if include is not None:
                    f.write('#include %s\n' % include)
                f.write('int main(void) {\n')
                f.write('    %s;\n' % funcname)
                f.write('}\n')
            # Redirect stderr to /dev/null to hide any error messages
            # from the compiler.
            # This will have to be changed if we ever have to check
            # for a function on Windows.
            devnull = open('/dev/null', 'w')
            oldstderr = os.dup(sys.stderr.fileno())
            os.dup2(devnull.fileno(), sys.stderr.fileno())
            objects = cc.compile([fname], output_dir=tmpdir,
                                 extra_postargs=extra_postargs)
            cc.link_executable(objects, os.path.join(tmpdir, "a.out"))
        except Exception:
            return False
        return True
    finally:
        if oldstderr is not None:
            os.dup2(oldstderr, sys.stderr.fileno())
        if devnull is not None:
            devnull.close()
        shutil.rmtree(tmpdir)


def detect_openmp():
    """Does this compiler support OpenMP parallelization?"""
    print("Attempting to autodetect OpenMP support... ", end="")
    compiler = new_compiler()
    compiler.add_library('gomp')
    include = '<omp.h>'
    extra_postargs = ['-fopenmp']
    hasopenmp = hasfunction(compiler, 'omp_get_num_threads()', include=include,
                            extra_postargs=extra_postargs)
    if hasopenmp:
        print("Compiler supports OpenMP")
    else:
        print("Did not detect OpenMP support.")
    return hasopenmp


def extensions(config):
    use_cython = config.get('use_cython', default=True)
    use_openmp = config.get('use_openmp', default=True)

    if config.get('debug_cflags', default=False):
        extra_compile_args = '\
            -std=c99 -pedantic -Wall -Wcast-align -Wcast-qual -Wpointer-arith \
            -Wchar-subscripts -Winline -Wnested-externs -Wbad-function-cast \
            -Wunreachable-code -Werror'
        define_macros = [('DEBUG', '1')]
    else:
        extra_compile_args = ''
        define_macros = []

    # Needed for large-file seeking under 32bit systems (for xtc/trr indexing
    # and access).
    largefile_macros = [
        ('_LARGEFILE_SOURCE', None),
        ('_LARGEFILE64_SOURCE', None),
        ('_FILE_OFFSET_BITS', '64')
    ]

    has_openmp = detect_openmp()

    if use_openmp and not has_openmp:
        print('No openmp compatible compiler found default to serial build.')

    parallel_args = ['-fopenmp'] if has_openmp and use_openmp else []
    parallel_libraries = ['gomp'] if has_openmp and use_openmp else []
    parallel_macros = [('PARALLEL', None)] if has_openmp and use_openmp else []

    if use_cython:
        if not cython_found:
            print("Couldn't find Cython installation. "
                  "Not recompiling cython extension")
            use_cython = False

    source_suffix = '.pyx' if use_cython else '.c'

    include_dirs = [get_numpy_include(), 'src/clustering']

    dcd = Extension('coordinates._dcdmodule',
                    ['MDAnalysis/coordinates/src/dcd.c'],
                    include_dirs=include_dirs + ['MDAnalysis/coordinates/include'],
                    define_macros=define_macros,
                    extra_compile_args=extra_compile_args)
    dcd_time = Extension('coordinates.dcdtimeseries',
                         ['MDAnalysis/coordinates/dcdtimeseries' + source_suffix],
                         include_dirs=include_dirs + ['MDAnalysis/coordinates/include'],
                         define_macros=define_macros,
                         extra_compile_args=extra_compile_args)
    distances = Extension('lib._distances',
                          ['MDAnalysis/lib/distances' + source_suffix],
                          include_dirs=include_dirs + ['MDAnalysis/lib/include'],
                          libraries=['m'],
                          define_macros=define_macros,
                          extra_compile_args=extra_compile_args)
    distances_omp = Extension('lib._distances_openmp',
                              ['MDAnalysis/lib/distances_openmp' + source_suffix],
                              include_dirs=include_dirs + ['MDAnalysis/lib/include'],
                              libraries=['m'] + parallel_libraries,
                              define_macros=define_macros + parallel_macros,
                              extra_compile_args=parallel_args,
                              extra_link_args=parallel_args)
    parallel_dist = Extension("lib.parallel.distances",
                              ['MDAnalysis/lib/distances_parallel' + source_suffix],
                              include_dirs=include_dirs,
                              libraries=['m'] + parallel_libraries,
                              extra_compile_args=parallel_args,
                              extra_link_args=parallel_args)
    qcprot = Extension('lib.qcprot',
                       ['MDAnalysis/lib/src/pyqcprot/pyqcprot' + source_suffix],
                       include_dirs=include_dirs,
                       extra_compile_args=["-O3", "-ffast-math"])
    transformation = Extension('lib._transformations',
                               ['MDAnalysis/lib/src/transformations/transformations.c'],
                               libraries=['m'],
                               define_macros=define_macros,
                               include_dirs=include_dirs,
                               extra_compile_args=extra_compile_args)
    xdr = Extension('coordinates.xdrfile._libxdrfile2',
                    sources=['MDAnalysis/coordinates/xdrfile/src/' + f
                           for f in ('libxdrfile2_wrap.c',
                                     'xdrfile.c',
                                     'xdrfile_trr.c',
                                     'xdrfile_xtc.c')])
    encore_utils = Extension('analysis.encore.cutils',
                            sources = ['MDAnalysis/lib/src/encore_cutils/cutils' + source_suffix],
                            include_dirs = include_dirs,
                            extra_compile_args = ["-O3", "-ffast-math"])
    ap_clustering = Extension('analysis.encore.clustering.affinityprop',
                            sources = ['MDAnalysis/lib/src/clustering/affinityprop' + source_suffix, "MDAnalysis/lib/src/clustering/ap.c"],
                            include_dirs = include_dirs,
                            libraries=["m"],
                            extra_compile_args=["-O3", "-ffast-math","-std=c99"])
    spe_dimred = Extension('analysis.encore.dimensionality_reduction.stochasticproxembed',
                            sources = ['MDAnalysis/lib/src/dimensionality_reduction/stochasticproxembed' + source_suffix, "MDAnalysis/lib/src/dimensionality_reduction/spe.c"],
                            include_dirs = include_dirs,
                            libraries=["m"],
                            extra_compile_args=["-O3", "-ffast-math","-std=c99"])

    return [dcd, dcd_time, distances, distances_omp, parallel_dist, qcprot,
            transformation, xdr, encore_utils, ap_clustering, spe_dimred]

if __name__ == '__main__':
    # NOTE: keep in sync with MDAnalysis.__version__ in version.py
    RELEASE = "0.12.1"
    with open("SUMMARY.txt") as summary:
        LONG_DESCRIPTION = summary.read()
    CLASSIFIERS = [
        'Development Status :: 4 - Beta',
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

    config = Config()

    setup(name='MDAnalysis',
          version=RELEASE,
          description='An object-oriented toolkit to analyze molecular dynamics '
          'trajectories generated by CHARMM, Gromacs, NAMD, LAMMPS, or Amber.',
          long_description=LONG_DESCRIPTION,
          author='Naveen Michaud-Agrawal',
          author_email='naveen.michaudagrawal@gmail.com',
          maintainer='Richard Gowers',
          maintainer_email='mdnalysis-discussion@googlegroups.com',
          url='http://www.mdanalysis.org',
          provides=['MDAnalysis'],
          license='GPL 2',
          packages=find_packages(),
          package_dir={'MDAnalysis': 'MDAnalysis'},
          ext_package='MDAnalysis',
          ext_modules=extensions(config),
          classifiers=CLASSIFIERS,
          cmdclass=cmdclass,
          requires=['numpy (>=1.5.0)', 'biopython',
                    'networkx (>=1.0)', 'GridDataFormats'],
          # all standard requirements are available through PyPi and
          # typically can be installed without difficulties through setuptools
          install_requires=[
              'biopython>=1.59',
              'networkx>=1.0',
              'GridDataFormats>=0.2.2',
          ],
          # extras can be difficult to install through setuptools and/or
          # you might prefer to use the version available through your
          # packaging system
          extras_require={
              'AMBER': ['netCDF4>=1.0'],  # for AMBER netcdf, also needs HDF5
                                          # and netcdf-4
              'analysis': [
                  'matplotlib',
                  'scipy',
                  'seaborn',  # for annotated heat map and nearest neighbor
                              # plotting in PSA
              ],
          },
          test_suite="MDAnalysisTests",
          tests_require=[
              'nose>=1.3.7',
              'MDAnalysisTests=={}'.format(RELEASE),  # same as this release!
          ],
          zip_safe=False,  # as a zipped egg the *.so files are not found (at
                           # least in Ubuntu/Linux)
    )
