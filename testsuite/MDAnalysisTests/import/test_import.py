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
from __future__ import absolute_import, print_function
import sys
import os
import subprocess
import pytest

"""Test if importing MDAnalysis has unwanted side effects (PR #1794)."""

@pytest.mark.skipif(os.name == 'nt',
                    reason="fork-related import checks irrelevant on Windows")
class TestMDAImport(object):
    # Tests concerning importing MDAnalysis.
    def test_os_dot_fork_not_called(self):
        # Importing MDAnalysis shall not trigger calls to os.fork (see PR #1794)
        # This test has to run in a separate Python instance to ensure that
        # no previously imported modules interfere with it. It is therefore
        # offloaded to the script "fork_called.py".
        loc = os.path.dirname(os.path.realpath(__file__))
        script = os.path.join(loc, 'fork_called.py')
        encoding = sys.stdout.encoding
        if encoding is None:
            encoding = "utf-8"
        # If the script we are about to call fails, we want its stderr to show
        # up in the pytest log. This is accomplished by redirecting the stderr
        # of the subprocess to its stdout and catching the raised
        # CalledProcessError. That error's output member then contains the
        # failed script's stderr and we can print it:
        try:
            out = subprocess.check_output([sys.executable, script],
                                          stderr=subprocess.STDOUT)\
                                         .decode(encoding)
        except subprocess.CalledProcessError as err:
            print(err.output)
            raise(err)

    def test_os_dot_fork_not_none(self):
        # In MDAnalysis.core.universe, os.fork is set to None prior to importing
        # the uuid module and restored afterwards (see PR #1794 for details).
        # This tests asserts that os.fork has been restored.
        import MDAnalysis
        assert os.fork is not None
