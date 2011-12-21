"""Setuptools-based setup script for tests of MDAnalysis.

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


if __name__ == '__main__':
    RELEASE = "0.7.5"
    LONG_DESCRIPTION = \
"""MDAnalysis is a tool for analyzing molecular dynamics trajectories.

This package contains the test code and the trajectory data that are
used for the test cases. In order to make downloads and binary package
maintenance more efficient, these tests were moved into this packe.
"""
    CLASSIFIERS = ['Development Status :: 4 - Beta',
                   'Environment :: Console',
                   'Intended Audience :: Scientists',
                   'License :: OSI Approved :: GNU General Public License (GPL)',
                   'Operating System :: POSIX',
                   'Operating System :: MacOS :: MacOS X',
                   'Programming Language :: Python',
                   'Topic :: Scientific/Engineering :: Bio-Informatics',
                   'Topic :: Scientific/Engineering :: Chemistry',
                   ]

    setup(name              = 'MDAnalysisTests',
          version           = RELEASE,
          description       = 'Python tools to support analysis of trajectories (test cases)',
          author            = 'Naveen Michaud-Agrawal',
          author_email      = 'naveen.michaudagrawal@gmail.com',
          url               = 'http://mdanalysis.googlecode.com/',
          license           = 'GPL 2',
          packages          = ['MDAnalysisTests'],
          package_dir       = {'MDAnalysisTests': 'MDAnalysisTests'},
          package_data      = {'MDAnalysisTests':
                                   ['data/*.psf','data/*.dcd','data/*.pdb',
                                    'data/*.gro', 'data/*.xtc','data/*.trr',
                                    'data/*.crd', 'data/*.xyz',
                                    'data/*.prmtop', 'data/*.trj', 'data/*.mdcrd',
                                    'data/*.pqr', 'data/*.pdbqt', 'data/*.bz2',
                                    ],
                               },
          classifiers       = CLASSIFIERS,
          long_description  = LONG_DESCRIPTION,
          zip_safe = False,   # had 'KeyError' as zipped egg (2MB savings are not worth the trouble)
          )
