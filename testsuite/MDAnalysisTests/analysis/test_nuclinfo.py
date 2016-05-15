# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
# Beckstein and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

from numpy.testing import (
    assert_,
    assert_almost_equal,
)

import MDAnalysis as mda
from MDAnalysis.analysis import nuclinfo

from MDAnalysisTests.datafiles import NUCL


class TestNuclinfo(object):
    def setUp(self):
        self.u = mda.Universe(NUCL)

    def tearDown(self):
        del self.u

    def test_wc_pair_1(self):
        val = nuclinfo.wc_pair(self.u, 1, 2, seg1='RNAA', seg2='RNAA')

        assert_almost_equal(val, 4.449, decimal=3)

    def test_wc_pair_2(self):
        val = nuclinfo.wc_pair(self.u, 22, 23, seg1='RNAA', seg2='RNAA')
        assert_almost_equal(val, 4.601, decimal=3)
