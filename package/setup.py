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
"""
from __future__ import print_function
from setuptools import setup, Extension, find_packages
from distutils.ccompiler import new_compiler
import codecs
import os
import sys
import shutil
import tempfile
import warnings

# Make sure I have the right Python version.
if sys.version_info[:2] < (2, 7):
    print('MDAnalysis requires Python 2.7 or better. Python {0:d}.{1:d} detected'.format(*
          sys.version_info[:2]))
    print('Please upgrade your version of Python.')
    sys.exit(-1)

if sys.version_info[0] < 3:
    import ConfigParser as configparser
    open_kwargs = {}
else:
    import configparser
    open_kwargs = {'encoding': 'utf-8'}

# Handle cython modules
try:
    from Cython.Distutils import build_ext
    cython_found = True
    cmdclass = {'build_ext': build_ext}
except ImportError:
    cython_found = False
    cmdclass = {}

# NOTE: keep in sync with MDAnalysis.__version__ in version.py
RELEASE = "0.15.0"

is_release = not 'dev' in RELEASE

if cython_found:
    # cython has to be >=0.16 to support cython.parallel
    import Cython
    from Cython.Build import cythonize
    from distutils.version import LooseVersion

    required_version = "0.16"

    if not LooseVersion(Cython.__version__) >= LooseVersion(required_version):
        # We don't necessarily die here. Maybe we already have
        #  the cythonized '.c' files.
        print("Cython version {0} was found but won't be used: version {1} "
              "or greater is required because it offers a handy "
              "parallelization module".format(
               Cython.__version__, required_version))
        cython_found = False
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
    Values passed to environment variables are checked (case-insensitively)
    for specific strings with boolean meaning: 'True' or '1' will cause `True`
    to be returned. '0' or 'False' cause `False` to be returned.

    """

    def __init__(self, fname='setup.cfg'):
        if os.path.exists(fname):
            self.config = configparser.SafeConfigParser()
            self.config.read(fname)

    def get(self, option_name, default=None):
        environ_name = 'MDA_' + option_name.upper()
        if environ_name in os.environ:
            val = os.environ[environ_name]
            if val.upper() in ('1', 'TRUE'):
                return True
            elif val.upper() in ('0', 'FALSE'):
                return False
            return val
        try:
            option = self.config.get('options', option_name)
            return option
        except configparser.NoOptionError:
            return default

class MDAExtension(Extension, object):
    """Derived class to cleanly handle setup-time (numpy) dependencies.
    """
    # The only setup-time numpy dependency comes when setting up its
    #  include dir.
    # The actual numpy import and call can be delayed until after pip
    #  has figured it must install numpy.
    # This is accomplished by passing the get_numpy_include function
    #  as one of the include_dirs. This derived Extension class takes
    #  care of calling it when needed.
    def __init__(self, *args, **kwargs):
        self._mda_include_dirs = []
        super(MDAExtension, self).__init__(*args, **kwargs)

    @property
    def include_dirs(self):
        if not self._mda_include_dirs:
            for item in self._mda_include_dir_args:
                try:
                    self._mda_include_dirs.append(item()) #The numpy callable
                except TypeError:
                    self._mda_include_dirs.append(item)
        return self._mda_include_dirs

    @include_dirs.setter
    def include_dirs(self, val):
        self._mda_include_dir_args = val

def get_numpy_include():
    # Obtain the numpy include directory. This logic works across numpy
    # versions.
    # setuptools forgets to unset numpy's setup flag and we get a crippled
    # version of it unless we do it ourselves.
    try:
        # Python 3 renamed the ``__builin__`` module into ``builtins``.
        # Here we import the python 2 or the python 3 version of the module
        # with the python 3 name. This could be done with ``six`` but that
        # module may not be installed at that point.
        import __builtin__ as builtins
    except ImportError:
        import builtins
    builtins.__NUMPY_SETUP__ = False
    try:
        import numpy as np
    except ImportError:
        print('*** package "numpy" not found ***')
        print('MDAnalysis requires a version of NumPy (>=1.5.0), even for setup.')
        print('Please get it from http://numpy.scipy.org/ or install it through '
              'your package manager.')
        sys.exit(-1)
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
                    f.write('#include {0!s}\n'.format(include))
                f.write('int main(void) {\n')
                f.write('    {0!s};\n'.format(funcname))
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
    # dev installs must build their own cythonized files.
    use_cython = config.get('use_cython', default=not is_release)
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
        print('Will attempt to use Cython.')
        if not cython_found:
            print("Couldn't find a Cython installation. "
                  "Not recompiling cython extensions.")
            use_cython = False
    else:
        print('Will not attempt to use Cython.')

    source_suffix = '.pyx' if use_cython else '.c'

    # The callable is passed so that it is only evaluated at install time.
    include_dirs = [get_numpy_include]

    dcd = MDAExtension('coordinates._dcdmodule',
                       ['MDAnalysis/coordinates/src/dcd.c'],
                       include_dirs=include_dirs + ['MDAnalysis/coordinates/include'],
                       define_macros=define_macros,
                       extra_compile_args=extra_compile_args)
    dcd_time = MDAExtension('coordinates.dcdtimeseries',
                         ['MDAnalysis/coordinates/dcdtimeseries' + source_suffix],
                         include_dirs=include_dirs + ['MDAnalysis/coordinates/include'],
                         define_macros=define_macros,
                         extra_compile_args=extra_compile_args)
    distances = MDAExtension('lib.c_distances',
                             ['MDAnalysis/lib/c_distances' + source_suffix],
                             include_dirs=include_dirs + ['MDAnalysis/lib/include'],
                             libraries=['m'],
                             define_macros=define_macros,
                             extra_compile_args=extra_compile_args)
    distances_omp = MDAExtension('lib.c_distances_openmp',
                                 ['MDAnalysis/lib/c_distances_openmp' + source_suffix],
                                 include_dirs=include_dirs + ['MDAnalysis/lib/include'],
                                 libraries=['m'] + parallel_libraries,
                                 define_macros=define_macros + parallel_macros,
                                 extra_compile_args=parallel_args,
                                 extra_link_args=parallel_args)
    qcprot = MDAExtension('lib.qcprot',
                          ['MDAnalysis/lib/qcprot' + source_suffix],
                          include_dirs=include_dirs,
                          extra_compile_args=["-O3", "-ffast-math"])
    transformation = MDAExtension('lib._transformations',
                                  ['MDAnalysis/lib/src/transformations/transformations.c'],
                                  libraries=['m'],
                                  define_macros=define_macros,
                                  include_dirs=include_dirs,
                                  extra_compile_args=extra_compile_args)
    libmdaxdr = MDAExtension('lib.formats.libmdaxdr',
                          sources=['MDAnalysis/lib/formats/libmdaxdr' + source_suffix,
                                   'MDAnalysis/lib/formats/src/xdrfile.c',
                                   'MDAnalysis/lib/formats/src/xdrfile_xtc.c',
                                   'MDAnalysis/lib/formats/src/xdrfile_trr.c',
                                   'MDAnalysis/lib/formats/src/trr_seek.c',
                                   'MDAnalysis/lib/formats/src/xtc_seek.c',
                                   ],
                          include_dirs=include_dirs + ['MDAnalysis/lib/formats/include',
                                                       'MDAnalysis/lib/formats'],
                          define_macros=largefile_macros)
    util = MDAExtension('lib.formats.cython_util',
                        sources=['MDAnalysis/lib/formats/cython_util' + source_suffix],
                        include_dirs=include_dirs)

    pre_exts = [dcd, dcd_time, distances, distances_omp, qcprot,
                  transformation, libmdaxdr, util]
    cython_generated = []
    if use_cython:
        extensions = cythonize(pre_exts)
        for pre_ext, post_ext in zip(pre_exts, extensions):
            for source in post_ext.sources:
                if source not in pre_ext.sources:
                    cython_generated.append(source)
    else:
        #Let's check early for missing .c files
        extensions = pre_exts
        for ext in extensions:
            for source in ext.sources:
                if not (os.path.isfile(source) and
                        os.access(source, os.R_OK)):
                    raise IOError("Source file '{}' not found. This might be "
                                "caused by a missing Cython install, or a "
                                "failed/disabled Cython build.".format(source))
    return extensions, cython_generated


def dynamic_author_list():
    """Generate __authors__ from AUTHORS

    This function generates authors.py that contains the list of the
    authors from the AUTHORS file. This avoids having that list maintained in
    several places. Note that AUTHORS is sorted chronologically while we want
    __authors__ in authors.py to be sorted alphabetically.

    The authors are written in AUTHORS as bullet points under the
    "Chronological list of authors" title.
    """
    authors = []
    with codecs.open('AUTHORS', encoding='utf-8') as infile:
        # An author is a bullet point under the title "Chronological list of
        # authors". We first want move the cursor down to the title of
        # interest.
        for line_no, line in enumerate(infile, start=1):
            if line[:-1] == "Chronological list of authors":
                break
        else:
            # If we did not break, it means we did not find the authors.
            raise IOError('EOF before the list of authors')
        # Skip the next line as it is the title underlining
        line = next(infile)
        line_no += 1
        if line[:4] != '----':
            raise IOError('Unexpected content on line {0}, '
                          'should be a string of "-".'.format(line_no))
        # Add each bullet point as an author until the next title underlining
        for line in infile:
            if line[:4] in ('----', '====', '~~~~'):
                # The previous line was a title, hopefully it did not start as
                # a bullet point so it got ignored. Since we hit a title, we
                # are done reading the list of authors.
                break
            elif line.strip()[:2] == '- ':
                # This is a bullet point, so it should be an author name.
                name = line.strip()[2:].strip()
                authors.append(name)

    # So far, the list of authors is sorted chronologically. We want it
    # sorted alphabetically of the last name.
    authors.sort(key=lambda name: name.split()[-1])
    # Move Naveen and Elizabeth first, and Oliver last.
    authors.remove('Naveen Michaud-Agrawal')
    authors.remove('Elizabeth J. Denning')
    authors.remove('Oliver Beckstein')
    authors = (['Naveen Michaud-Agrawal', 'Elizabeth J. Denning']
               + authors + ['Oliver Beckstein'])

    # Write the authors.py file.
    out_path = 'MDAnalysis/authors.py'
    with codecs.open(out_path, 'w', encoding='utf-8') as outfile:
        # Write the header
        header = '''\
#-*- coding:utf-8 -*-

# This file is generated from the AUTHORS file during the installation process.
# Do not edit it as your changes will be overwritten.
'''
        print(header, file=outfile)

        # Write the list of authors as a python list
        template = u'__authors__ = [\n{}\n]'
        author_string = u',\n'.join(u'    u"{}"'.format(name)
                                    for name in authors)
        print(template.format(author_string), file=outfile)


if __name__ == '__main__':
    try:
        dynamic_author_list()
    except (OSError, IOError):
        warnings.warn('Cannot write the list of authors.')

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
    exts, cythonfiles = extensions(config)

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
          download_url='https://github.com/MDAnalysis/mdanalysis/releases',
          provides=['MDAnalysis'],
          license='GPL 2',
          packages=find_packages(),
          package_dir={'MDAnalysis': 'MDAnalysis'},
          ext_package='MDAnalysis',
          ext_modules=exts,
          classifiers=CLASSIFIERS,
          cmdclass=cmdclass,
          requires=['numpy (>=1.5.0)', 'biopython',
                    'networkx (>=1.0)', 'GridDataFormats (>=0.3.2)'],
          # all standard requirements are available through PyPi and
          # typically can be installed without difficulties through setuptools
          setup_requires=[
              'numpy>=1.5.0',
          ],
          install_requires=[
              'numpy>=1.5.0',
              'biopython>=1.59',
              'networkx>=1.0',
              'GridDataFormats>=0.3.2',
              'six>=1.4.0',
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
              'MDAnalysisTests=={0}'.format(RELEASE),  # same as this release!
          ],
          zip_safe=False,  # as a zipped egg the *.so files are not found (at
                           # least in Ubuntu/Linux)
    )

    # Releases keep their cythonized stuff for shipping.
    if not config.get('keep_cythonized', default=is_release):
        for cythonized in cythonfiles:
            try:
                os.unlink(cythonized)
            except OSError as err:
                print("Warning: failed to delete cythonized file {0}: {1}. "
                    "Moving on.".format(cythonized, err.strerror))

