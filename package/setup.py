#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""Setuptools-based setup script for MDAnalysis.

A working installation of NumPy <http://numpy.scipy.org> is required.

For a basic installation just type the command::

  pip install .

For more in-depth instructions, see the installation section at the
MDAnalysis User Guide:

  https://userguide.mdanalysis.org/stable/installation.html

Also free to ask on GitHub Discussions for help:

  https://github.com/MDAnalysis/mdanalysis/discussions/categories/installation

"""

from setuptools import setup, Extension, find_packages
from distutils.ccompiler import new_compiler
from distutils.sysconfig import customize_compiler
from packaging.version import Version
import codecs
import os
import sys
import re
import shutil
import tempfile
import warnings
import platform

import configparser
from subprocess import getoutput

# NOTE: keep in sync with MDAnalysis.__version__ in version.py
RELEASE = "2.8.0-dev0"

is_release = 'dev' not in RELEASE

# Handle cython modules
try:
    # cython has to be >=0.16 <0.28 to support cython.parallel
    # minimum cython version now set to 0.28 to match pyproject.toml
    import Cython
    from Cython.Build import cythonize
    cython_found = True

    required_version = "0.28"
    if not Version(Cython.__version__) >= Version(required_version):
        # We don't necessarily die here. Maybe we already have
        #  the cythonized '.c' files.
        print("Cython version {0} was found but won't be used: version {1} "
              "or greater is required because it offers a handy "
              "parallelization module".format(
               Cython.__version__, required_version))
        cython_found = False
    cython_linetrace = bool(os.environ.get('CYTHON_TRACE_NOGIL', False))
except ImportError:
    cython_found = False
    if not is_release:
        print("*** package: Cython not found ***")
        print("MDAnalysis requires cython for development builds")
        sys.exit(1)
    cython_linetrace = False

def abspath(file):
    return os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        file)

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
        fname = abspath(fname)
        if os.path.exists(fname):
            self.config = configparser.ConfigParser()
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
    def __init__(self, name, sources, *args, **kwargs):
        self._mda_include_dirs = []
        # don't abspath sources else packaging fails on Windows (issue #3129)
        super(MDAExtension, self).__init__(name, sources, *args, **kwargs)

    @property
    def include_dirs(self):
        if not self._mda_include_dirs:
            for item in self._mda_include_dir_args:
                try:
                    self._mda_include_dirs.append(item()) #The numpy callable
                except TypeError:
                    item = abspath(item)
                    self._mda_include_dirs.append((item))
        return self._mda_include_dirs

    @include_dirs.setter
    def include_dirs(self, val):
        self._mda_include_dir_args = val


def get_numpy_include():
    # Obtain the numpy include directory. This logic works across numpy
    # versions.
    # setuptools forgets to unset numpy's setup flag and we get a crippled
    # version of it unless we do it ourselves.
    import builtins

    builtins.__NUMPY_SETUP__ = False
    try:
        import numpy as np
    except ImportError:
        print('*** package "numpy" not found ***')
        print('MDAnalysis requires a version of NumPy (>=1.21.0), even for setup.')
        print('Please get it from http://numpy.scipy.org/ or install it through '
              'your package manager.')
        sys.exit(-1)
    return np.get_include()


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
    customize_compiler(compiler)
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


def using_clang():
    """Will we be using a clang compiler?"""
    compiler = new_compiler()
    customize_compiler(compiler)
    compiler_ver = getoutput("{0} -v".format(compiler.compiler[0]))
    if 'Spack GCC' in compiler_ver:
        # when gcc toolchain is built from source with spack
        # using clang, the 'clang' string may be present in
        # the compiler metadata, but it is not clang
        is_clang = False
    elif 'clang' in compiler_ver:
        # by default, Apple will typically alias gcc to
        # clang, with some mention of 'clang' in the
        # metadata
        is_clang = True
    else:
        is_clang = False
    return is_clang


def extensions(config):
    # usually (except coming from release tarball) cython files must be generated
    use_cython = config.get('use_cython', default=cython_found)
    use_openmp = config.get('use_openmp', default=True)
    annotate_cython = config.get('annotate_cython', default=False)

    extra_compile_args = ['-std=c99', '-O3', '-funroll-loops',
                          '-fsigned-zeros'] # see #2722
    define_macros = []
    if config.get('debug_cflags', default=False):
        extra_compile_args.extend(['-Wall', '-pedantic'])
        define_macros.extend([('DEBUG', '1')])

    # encore is sensitive to floating point accuracy, especially on non-x86
    # to avoid reducing optimisations on everything, we make a set of compile
    # args specific to encore see #2997 for an example of this.
    encore_compile_args = [a for a in extra_compile_args if 'O3' not in a]
    if platform.machine() == 'aarch64' or platform.machine() == 'ppc64le':
        encore_compile_args.append('-O1')
    else:
        encore_compile_args.append('-O3')

    # allow using custom c/c++ flags and architecture specific instructions.
    # This allows people to build optimized versions of MDAnalysis.
    # Do here so not included in encore
    extra_cflags = config.get('extra_cflags', default=False)
    if extra_cflags:
        flags = extra_cflags.split()
        extra_compile_args.extend(flags)

    cpp_extra_compile_args = [a for a in extra_compile_args if 'std' not in a]
    cpp_extra_compile_args.append('-std=c++11')
    cpp_extra_link_args=[]
    # needed to specify c++ runtime library on OSX
    if platform.system() == 'Darwin' and using_clang():
        cpp_extra_compile_args.append('-stdlib=libc++')
        cpp_extra_compile_args.append('-mmacosx-version-min=10.9')
        cpp_extra_link_args.append('-stdlib=libc++')
        cpp_extra_link_args.append('-mmacosx-version-min=10.9')

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
    cpp_source_suffix = '.pyx' if use_cython else '.cpp'

    # The callable is passed so that it is only evaluated at install time.

    include_dirs = [get_numpy_include]
    # Windows automatically handles math library linking
    # and will not build MDAnalysis if we try to specify one
    if os.name == 'nt':
        mathlib = []
    else:
        mathlib = ['m']

    if cython_linetrace:
        extra_compile_args.append("-DCYTHON_TRACE_NOGIL")
        cpp_extra_compile_args.append("-DCYTHON_TRACE_NOGIL")

    libdcd = MDAExtension('MDAnalysis.lib.formats.libdcd',
                          ['MDAnalysis/lib/formats/libdcd' + source_suffix],
                          include_dirs=include_dirs + ['MDAnalysis/lib/formats/include'],
                          define_macros=define_macros,
                          extra_compile_args=extra_compile_args)
    distances = MDAExtension('MDAnalysis.lib.c_distances',
                             ['MDAnalysis/lib/c_distances' + source_suffix],
                             include_dirs=include_dirs + ['MDAnalysis/lib/include'],
                             libraries=mathlib,
                             define_macros=define_macros,
                             extra_compile_args=extra_compile_args)
    distances_omp = MDAExtension('MDAnalysis.lib.c_distances_openmp',
                                 ['MDAnalysis/lib/c_distances_openmp' + source_suffix],
                                 include_dirs=include_dirs + ['MDAnalysis/lib/include'],
                                 libraries=mathlib + parallel_libraries,
                                 define_macros=define_macros + parallel_macros,
                                 extra_compile_args=parallel_args + extra_compile_args,
                                 extra_link_args=parallel_args)
    qcprot = MDAExtension('MDAnalysis.lib.qcprot',
                          ['MDAnalysis/lib/qcprot' + source_suffix],
                          include_dirs=include_dirs,
                          define_macros=define_macros,
                          extra_compile_args=extra_compile_args)
    transformation = MDAExtension('MDAnalysis.lib._transformations',
                                  ['MDAnalysis/lib/src/transformations/transformations.c'],
                                  libraries=mathlib,
                                  define_macros=define_macros,
                                  include_dirs=include_dirs,
                                  extra_compile_args=extra_compile_args)
    libmdaxdr = MDAExtension('MDAnalysis.lib.formats.libmdaxdr',
                             sources=['MDAnalysis/lib/formats/libmdaxdr' + source_suffix,
                                      'MDAnalysis/lib/formats/src/xdrfile.c',
                                      'MDAnalysis/lib/formats/src/xdrfile_xtc.c',
                                      'MDAnalysis/lib/formats/src/xdrfile_trr.c',
                                      'MDAnalysis/lib/formats/src/trr_seek.c',
                                      'MDAnalysis/lib/formats/src/xtc_seek.c',
                             ],
                             include_dirs=include_dirs + ['MDAnalysis/lib/formats/include',
                                                          'MDAnalysis/lib/formats'],
                             define_macros=largefile_macros + define_macros,
                             extra_compile_args=extra_compile_args)
    util = MDAExtension('MDAnalysis.lib.formats.cython_util',
                        sources=['MDAnalysis/lib/formats/cython_util' + source_suffix],
                        include_dirs=include_dirs,
                        define_macros=define_macros,
                        extra_compile_args=extra_compile_args)
    cutil = MDAExtension('MDAnalysis.lib._cutil',
                         sources=['MDAnalysis/lib/_cutil' + cpp_source_suffix],
                         language='c++',
                         libraries=mathlib,
                         include_dirs=include_dirs + ['MDAnalysis/lib/include'],
                         define_macros=define_macros,
                         extra_compile_args=cpp_extra_compile_args,
                         extra_link_args= cpp_extra_link_args)
    augment = MDAExtension('MDAnalysis.lib._augment',
                         sources=['MDAnalysis/lib/_augment' + cpp_source_suffix],
                         language='c++',
                         include_dirs=include_dirs,
                         define_macros=define_macros,
                         extra_compile_args=cpp_extra_compile_args,
                         extra_link_args= cpp_extra_link_args)
    timestep = MDAExtension('MDAnalysis.coordinates.timestep',
                         sources=['MDAnalysis/coordinates/timestep' + cpp_source_suffix],
                         language='c++',
                         include_dirs=include_dirs,
                         define_macros=define_macros,
                         extra_compile_args=cpp_extra_compile_args,
                         extra_link_args= cpp_extra_link_args)


    encore_utils = MDAExtension('MDAnalysis.analysis.encore.cutils',
                                sources=['MDAnalysis/analysis/encore/cutils' + source_suffix],
                                include_dirs=include_dirs,
                                define_macros=define_macros,
                                extra_compile_args=encore_compile_args)
    ap_clustering = MDAExtension('MDAnalysis.analysis.encore.clustering.affinityprop',
                                 sources=['MDAnalysis/analysis/encore/clustering/affinityprop' + source_suffix,
                                          'MDAnalysis/analysis/encore/clustering/src/ap.c'],
                                 include_dirs=include_dirs+['MDAnalysis/analysis/encore/clustering/include'],
                                 libraries=mathlib,
                                 define_macros=define_macros,
                                 extra_compile_args=encore_compile_args)
    spe_dimred = MDAExtension('MDAnalysis.analysis.encore.dimensionality_reduction.stochasticproxembed',
                              sources=['MDAnalysis/analysis/encore/dimensionality_reduction/stochasticproxembed' + source_suffix,
                                       'MDAnalysis/analysis/encore/dimensionality_reduction/src/spe.c'],
                              include_dirs=include_dirs+['MDAnalysis/analysis/encore/dimensionality_reduction/include'],
                              libraries=mathlib,
                              define_macros=define_macros,
                              extra_compile_args=encore_compile_args)
    nsgrid = MDAExtension('MDAnalysis.lib.nsgrid',
                             ['MDAnalysis/lib/nsgrid' + cpp_source_suffix],
                             include_dirs=include_dirs + ['MDAnalysis/lib/include'],
                             language='c++',
                             define_macros=define_macros,
                             extra_compile_args=cpp_extra_compile_args,
                             extra_link_args= cpp_extra_link_args)
    pre_exts = [libdcd, distances, distances_omp, qcprot,
                transformation, libmdaxdr, util, encore_utils,
                ap_clustering, spe_dimred, cutil, augment, nsgrid, timestep]


    cython_generated = []
    if use_cython:
        extensions = cythonize(
            pre_exts,
            annotate=annotate_cython,
            compiler_directives={'linetrace': cython_linetrace,
                                 'embedsignature': False,
                                 'language_level': '3'},
        )
        if cython_linetrace:
            print("Cython coverage will be enabled")
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
    with codecs.open(abspath('AUTHORS'), encoding='utf-8') as infile:
        # An author is a bullet point under the title "Chronological list of
        # authors". We first want move the cursor down to the title of
        # interest.
        for line_no, line in enumerate(infile, start=1):
            if line.rstrip() == "Chronological list of authors":
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
    out_path = abspath('MDAnalysis/authors.py')
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


def long_description(readme):
    """Create reST SUMMARY file for PyPi."""

    with open(abspath(readme)) as summary:
        buffer = summary.read()
    # remove top heading that messes up pypi display
    m = re.search('====*\n[^\n]*README[^\n]*\n=====*\n', buffer,
                  flags=re.DOTALL)
    assert m, "README.rst does not contain a level-1 heading"
    return buffer[m.end():]


if __name__ == '__main__':
    try:
        dynamic_author_list()
    except (OSError, IOError):
        warnings.warn('Cannot write the list of authors.')

    try:
        # when building from repository for creating the distribution
        LONG_DESCRIPTION = long_description("../README.rst")
    except OSError:
        # when building from a tar file for installation
        # (LONG_DESCRIPTION is not really needed)
        LONG_DESCRIPTION = "MDAnalysis -- https://www.mdanalysis.org/"

    config = Config()
    exts, cythonfiles = extensions(config)

    setup(name='MDAnalysis',
          version=RELEASE,
          long_description=LONG_DESCRIPTION,
          long_description_content_type='text/x-rst',
          license='GPL-3.0-or-later',
          provides=['MDAnalysis'],
          packages=find_packages(),
          package_data={'MDAnalysis':
                        ['analysis/data/*.npy',
                        ],
          },
          ext_modules=exts,
          test_suite="MDAnalysisTests",
          tests_require=[
              'MDAnalysisTests=={0!s}'.format(RELEASE),  # same as this release!
          ],
    )

    # Releases keep their cythonized stuff for shipping.
    if not config.get('keep_cythonized', default=is_release) and not cython_linetrace:
        for cythonized in cythonfiles:
            try:
                os.unlink(cythonized)
            except OSError as err:
                print("Warning: failed to delete cythonized file {0}: {1}. "
                    "Moving on.".format(cythonized, err.strerror))
