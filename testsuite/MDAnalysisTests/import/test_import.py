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
from __future__ import absolute_import, print_function
import sys
import os
import subprocess

"""Test if importing MDAnalysis has unwanted side effects (PR #1794)."""

class TestMDAImport(object):
    # Tests concerning importing MDAnalysis.
    def test_os_dot_fork_not_called(self):
        # This test has to run in a separate Python instance and is therefore
        # offloaded to the script "fork_called.py".
        loc = __file__
        if loc.endswith('.pyc') and os.path.exists(loc[:-1]):
            loc = loc[:-1]
        loc = os.path.dirname(os.path.realpath(loc))
        script = os.path.join(loc, 'fork_called.py')
        encoding = sys.stdout.encoding
        if encoding is None:
            encoding = "utf-8"
        try:
            out = subprocess.check_output([sys.executable, script],
                                          stderr=subprocess.STDOUT)\
                                         .decode(encoding)
        except subprocess.CalledProcessError as err:
            print(err.output)
            raise(err)

    def test_os_dot_fork_not_none(self):
        # This test has to run in a separate Python instance and is therefore
        # offloaded to the script "fork_restored.py".
        loc = __file__
        if loc.endswith('.pyc') and os.path.exists(loc[:-1]):
            loc = loc[:-1]
        loc = os.path.dirname(os.path.realpath(loc))
        script = os.path.join(loc, 'fork_restored.py')
        encoding = sys.stdout.encoding
        if encoding is None:
            encoding = "utf-8"
        try:
            out = subprocess.check_output([sys.executable, script],
                                          stderr=subprocess.STDOUT)\
                                         .decode(encoding)
        except subprocess.CalledProcessError as err:
            print(err.output)
            raise(err)
