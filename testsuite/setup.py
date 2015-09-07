#!/usr/bin/env python
# coding=utf-8
"""Setuptools-based setup script for tests of MDAnalysis.

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

# ------------------------------------------------------------
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
from __future__ import print_function
from ez_setup import use_setuptools

use_setuptools()
from setuptools import setup, Extension
#
#------------------------------------------------------------

import sys
import os
import glob

# Make sure I have the right Python version.
if sys.version_info[:2] < (2, 6):
    print("MDAnalysis requires Python 2.6 or better. Python %d.%d detected" %
          sys.version_info[:2])
    print("Please upgrade your version of Python.")
    sys.exit(-1)


if __name__ == '__main__':
    RELEASE = "0.11.1-dev"  # this must be in-sync with MDAnalysis
    LONG_DESCRIPTION = \
        """MDAnalysis is a tool for analyzing molecular dynamics trajectories.

This package (MDAnalysisTests) contains the test code and the trajectory data
that are used for the test cases. In order to make downloads and binary package
maintenance more efficient, these tests were moved into this package.

For details see the report for `Issue 87`_.

.. _`Issue 87`: https://github.com/MDAnalysis/mdanalysis/issues/87
"""
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

    setup(name='MDAnalysisTests',
          version=RELEASE,
          description='Python tools to support analysis of trajectories (test cases)',
          author='Naveen Michaud-Agrawal',
          author_email='naveen.michaudagrawal@gmail.com',
          url='http://www.mdanalysis.org',
          license='GPL 2',
          packages=['MDAnalysisTests', 'MDAnalysisTests.plugins'],
          package_dir={'MDAnalysisTests': 'MDAnalysisTests',
                       'MDAnalysisTests.plugins': 'MDAnalysisTests/plugins'},
          package_data={'MDAnalysisTests':
              [
                  'data/*.psf', 'data/*.dcd', 'data/*.pdb',
                  'data/tprs/*.tpr', 'data/*.tpr',
                  'data/*.gro', 'data/*.xtc', 'data/*.trr', 'data/*npy',
                  'data/*.crd', 'data/*.xyz',
                  'data/*.prmtop', 'data/*.top', 'data/*.trj', 'data/*.mdcrd', 'data/*.ncdf',
                  'data/*.pqr', 'data/*.pdbqt', 'data/*.bz2',
                  'data/*.fasta',
                  'data/*.dms',
                  'data/merge/2zmm/*.pdb',
                  'data/*.trz',
                  'data/mol2/*.mol2',
                  'data/capping/*.gro', 'data/capping/*.pdb',
                  'data/*.data',
                  'data/gms/*.xyz', 'data/gms/*.gms', 'data/gms/*.gms.gz',
                  'data/*.inpcrd',
                  'data/dlpoly/CONFIG*',
                  'data/dlpoly/HISTORY*',
                  'data/*.xml',
              ],
          },
          classifiers=CLASSIFIERS,
          long_description=LONG_DESCRIPTION,
          install_requires=[
              'MDAnalysis==%s' % RELEASE,  # same as this release!
              'numpy>=1.5',
              'nose>=1.3.7',
          ],
          zip_safe=False,  # had 'KeyError' as zipped egg (2MB savings are not worth the trouble)
          )
