# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import print_function

import MDAnalysis.analysis.helanal
from MDAnalysis import FinishTimeException

from numpy.testing import TestCase, assert_raises
import numpy as np

import os
import tempdir

from MDAnalysisTests.datafiles import GRO, XTC


class Test_Helanal(TestCase):
    def setUp(self):
        self.universe = MDAnalysis.Universe(GRO, XTC)
        self.selection = 'name CA'

    def tearDown(self):
        del self.universe

    def test_xtc_striding(self):
        """MDAnalysis.analysis.helanal: Check for resolution of Issue 188."""
        u = self.universe
        u.trajectory[1]

        with tempdir.in_tempdir():
            assert_raises(FinishTimeException,
                          MDAnalysis.analysis.helanal.helanal_trajectory,
                          u, selection=self.selection, finish=5)

        #with assert_raises(FinishTimeException):
        #    try:
        #        MDAnalysis.analysis.helanal.helanal_trajectory(u, selection=sel, finish=5)
         #   except IndexError:
         #       self.fail("IndexError consistent with Issue 188.")

