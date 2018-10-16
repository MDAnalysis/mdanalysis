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
from __future__ import print_function
from setuptools import setup, find_packages

import codecs
import sys
import warnings


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
    authors = (['Naveen Michaud-Agrawal', 'Elizabeth J. Denning'] +
               authors + ['Oliver Beckstein'])

    # Write the authors.py file.
    out_path = 'MDAnalysisTests/authors.py'
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


# Make sure I have the right Python version.
if sys.version_info[:2] < (2, 7):
    print("MDAnalysis requires Python 2.7 or better. "
          "Python {0:d}.{1:d} detected".format(*sys.version_info[:2]))
    print("Please upgrade your version of Python.")
    sys.exit(-1)


if __name__ == '__main__':
    try:
        dynamic_author_list()
    except (OSError, IOError) as e:
        warnings.warn('Cannot write the list of authors. '
                      '{}'.format(e))

    # this must be in-sync with MDAnalysis
    RELEASE = "0.19.1-dev"
    with open("README") as summary:
        LONG_DESCRIPTION = summary.read()

    CLASSIFIERS = [
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
        'Operating System :: POSIX',
        'Operating System :: MacOS :: MacOS X',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
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
                         'data/tprs/*.tpr', 'data/tprs/all_bonded/*.tpr',
                         'data/tprs/all_bonded/*.gro',
                         'data/tprs/all_bonded/*.top',
                         'data/tprs/all_bonded/*.mdp', 'data/*.tpr',
                         'data/*.gro', 'data/*.xtc', 'data/*.trr', 'data/*npy',
                         'data/*.crd', 'data/*.xyz',
                         'data/Amber/*.bz2',
                         'data/Amber/*.prmtop', 'data/Amber/*.top',
                         'data/Amber/*.parm7',
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
                        ],
          },
          install_requires=[
              'MDAnalysis=={0!s}'.format(RELEASE),  # same as this release!
              'pytest>=3.1.2',
              'hypothesis',
              'psutil>=4.0.2',
              'mock>=2.0.0',  # replace with unittest.mock in python 3 only version
          ],
          # had 'KeyError' as zipped egg (2MB savings are not worth the
          # trouble)
          zip_safe=False,
          )
