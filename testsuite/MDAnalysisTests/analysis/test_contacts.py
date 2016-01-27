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

import MDAnalysis
import MDAnalysis.analysis.contacts
from MDAnalysis import SelectionError

from numpy.testing import (TestCase, dec,
                           assert_almost_equal, assert_raises, assert_equal)
import numpy as np
import nose
from nose.plugins.attrib import attr

import os
import tempdir

from MDAnalysisTests.datafiles import PSF, DCD
from MDAnalysisTests import executable_not_found, parser_not_found

class TestContactAnalysis1(TestCase):
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def setUp(self):
        self.universe = MDAnalysis.Universe(PSF, DCD)
        self.trajectory = self.universe.trajectory

    def tearDown(self):
        del self.universe
        del self.trajectory

    def _run_ContactAnalysis1(self, **runkwargs):
        sel_basic = "(resname ARG or resname LYS) and (name NH* or name NZ)"
        sel_acidic = "(resname ASP or resname GLU) and (name OE* or name OD*)"
        acidic = self.universe.select_atoms(sel_acidic)
        basic = self.universe.select_atoms(sel_basic)
        outfile = 'qsalt.dat'
        CA1 = MDAnalysis.analysis.contacts.ContactAnalysis1(
            self.universe,
            selection=(sel_acidic, sel_basic), refgroup=(acidic, basic),
            radius=6.0, outfile=outfile)
        kwargs = runkwargs.copy()
        kwargs['force'] = True
        CA1.run(**kwargs)
        return CA1

    def test_startframe(self):
        """test_startframe: TestContactAnalysis1: start frame set to 0 (resolution of Issue #624)"""
        with tempdir.in_tempdir():
            CA1 = self._run_ContactAnalysis1()
            self.assertEqual(CA1.timeseries.shape[1], self.universe.trajectory.n_frames)

    def test_end_zero(self):
        """test_end_zero: TestContactAnalysis1: stop frame 0 is not ignored"""
        with tempdir.in_tempdir():
            CA1 = self._run_ContactAnalysis1(stop=0)
            self.assertEqual(len(CA1.timeseries), 0)

    def test_slicing(self):
        start, stop, step = 10, 30, 5
        with tempdir.in_tempdir():
            CA1 = self._run_ContactAnalysis1(start=start, stop=stop, step=step)
            frames = np.arange(self.universe.trajectory.n_frames)[start:stop:step]
            self.assertEqual(CA1.timeseries.shape[1], len(frames))
