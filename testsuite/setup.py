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
"""
from setuptools import setup, find_packages
from setuptools.command import sdist

import os
import shutil
import codecs
import sys
import warnings


class MDA_SDist(sdist.sdist):
    # To avoid having duplicate AUTHORS file...
    def run(self):
        here = os.path.dirname(os.path.abspath(__file__))
        has_authors = os.path.exists(os.path.join(here, 'AUTHORS'))

        if not has_authors:
            # If there is no AUTHORS file here, lets hope we're in
            # a repo checkout and grab from '../package'
            print("Grabbing AUTHORS file...")
            repo_root = os.path.split(here)[0]
            try:
                shutil.copyfile(
                    os.path.join(repo_root, 'package', 'AUTHORS'),
                    os.path.join(here, 'AUTHORS'))
            except:
                raise IOError("Couldn't grab AUTHORS file")
            else:
                copied_authors = True
        try:
            super(MDA_SDist, self).run()
        finally:
            if not has_authors and copied_authors:
                os.remove(os.path.join(here, 'AUTHORS'))


# Make sure I have the right Python version.
if sys.version_info[:2] < (3, 6):
    print("MDAnalysis requires Python 3.6 or better. "
          "Python {0:d}.{1:d} detected".format(*sys.version_info[:2]))
    print("Please upgrade your version of Python.")
    sys.exit(-1)


if __name__ == '__main__':
    # this must be in-sync with MDAnalysis
    RELEASE = "2.0.0-dev0"
    with open("README") as summary:
        LONG_DESCRIPTION = summary.read()

    CLASSIFIERS = [
        'Development Status :: 6 - Mature',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
        'Operating System :: POSIX',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows ',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: C',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ]

    setup(name='MDAnalysisTests',
          version=RELEASE,
          description='MDAnalysis testsuite',
          long_description=LONG_DESCRIPTION,
          long_description_content_type='text/x-rst',
          author='Naveen Michaud-Agrawal',
          author_email='naveen.michaudagrawal@gmail.com',
          maintainer='Richard Gowers',
          maintainer_email='mdnalysis-discussion@googlegroups.com',
          url='https://www.mdanalysis.org',
          download_url='https://github.com/MDAnalysis/mdanalysis/releases',
          project_urls={'Documentation': 'https://github.com/MDAnalysis/mdanalysis/wiki/UnitTests',
                        'CI Tests': 'https://travis-ci.org/MDAnalysis/mdanalysis',
                        'CI Coverage': 'https://codecov.io/gh/MDAnalysis/mdanalysis',
                        'Developer Group': 'https://groups.google.com/forum/#!forum/mdnalysis-devel',
                        'Issue Tracker': 'https://github.com/mdanalysis/mdanalysis/issues',
                        'Source': 'https://github.com/mdanalysis/mdanalysis',
                        },
          license='GPL 2',
          classifiers=CLASSIFIERS,
          packages=find_packages(),
          package_dir={'MDAnalysisTests': 'MDAnalysisTests',
                       'MDAnalysisTests.plugins': 'MDAnalysisTests/plugins'},
          package_data={'MDAnalysisTests':
                        ['data/*.psf', 'data/*.dcd', 'data/*.pdb',
                         'data/tprs/*.tpr',
                         'data/tprs/all_bonded/*.gro',
                         'data/tprs/all_bonded/*.top',
                         'data/tprs/all_bonded/*.mdp', 'data/*.tpr',
                         'data/tprs/*/*.tpr',
                         'data/*.gro', 'data/*.xtc', 'data/*.trr', 'data/*npy',
                         'data/*.crd', 'data/*.xyz',
                         'data/Amber/*.bz2',
                         'data/Amber/*.prmtop', 'data/Amber/*.top',
                         'data/Amber/*.parm7',
                         'data/Amber/*.rst7',
                         'data/Amber/*.trj', 'data/Amber/*.mdcrd',
                         'data/Amber/*.ncdf', 'data/Amber/*.nc',
                         'data/Amber/*.inpcrd',
                         'data/*.pqr', 'data/*.pdbqt', 'data/*.bz2', 'data/*.gz',
                         'data/*.ent',
                         'data/*.fasta',
                         'data/*.dat',
                         'data/*.dms',
                         'data/merge/2zmm/*.pdb',
                         'data/*.trz',
                         'data/mol2/*.mol2',
                         'data/contacts/*.gro.bz2', 'data/contacts/*.dat',
                         'data/capping/*.gro', 'data/capping/*.pdb',
                         'data/lammps/*',
                         'data/gms/*.xyz', 'data/gms/*.gms',
                         'data/gms/*.gms.gz',
                         'data/*.inpcrd',
                         'data/dlpoly/CONFIG*',
                         'data/dlpoly/HISTORY*',
                         'data/*.xml',
                         'data/coordinates/*',
                         'data/*xvg',
                         'data/*.mmtf', 'data/*.mmtf.gz',
                         'data/analysis/*',
                         'data/*.gsd',
                         'data/windows/*',
                         'data/*.itp', "data/gromacs/gromos54a7_edited.ff/*",
                         'data/*.coor',
                         'data/*.h5md',
                         'data/*.in',
                         'data/*.top',
                         'data/*.sdf',
                        ],
          },
          install_requires=[
              'MDAnalysis=={0!s}'.format(RELEASE),  # same as this release!
              'pytest>=3.3.0', # Raised to 3.3.0 due to Issue 2329
              'hypothesis',
              'psutil>=4.0.2',
          ],
          # had 'KeyError' as zipped egg (2MB savings are not worth the
          # trouble)
          zip_safe=False,
          cmdclass={'sdist': MDA_SDist},
          )
