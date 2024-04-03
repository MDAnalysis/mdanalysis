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

A working installation of MDAnalysis of the same version is required.

For a basic installation just type the command::

  pip install .

For more in-depth instructions, see the installation section in the
MDAnalysis User Guide:

  https://userguide.mdanalysis.org/stable/installation.html

Also free to ask on GitHub Discussions for help:

  https://github.com/MDAnalysis/mdanalysis/discussions

"""
from setuptools import setup
from setuptools.command import sdist

import os
import shutil


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


if __name__ == '__main__':
    # this must be in-sync with MDAnalysis
    RELEASE = "2.8.0-dev0"

    setup(
        version=RELEASE,
        install_requires=[
            'MDAnalysis=={0!s}'.format(RELEASE),  # same as this release!
            'pytest>=3.3.0', # Raised to 3.3.0 due to Issue 2329
            'hypothesis',
        ],
        cmdclass={'sdist': MDA_SDist},
    )
