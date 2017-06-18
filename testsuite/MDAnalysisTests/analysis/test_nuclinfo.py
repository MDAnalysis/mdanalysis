# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.mdanalysis.org
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
from __future__ import absolute_import

from numpy.testing import (
    assert_,
    assert_almost_equal,
    assert_allclose,
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

    def test_minor_pair_1(self):
        val = nuclinfo.minor_pair(self.u, 3, 17, seg1='RNAA', seg2='RNAA')
        assert_almost_equal(val, 16.649, decimal=3)

    def test_minor_pair_2(self):
        val = nuclinfo.minor_pair(self.u, 20, 5, seg1='RNAA', seg2='RNAA')
        assert_almost_equal(val, 4.356, decimal=3)

    def test_major_pair_1(self):
        val = nuclinfo.major_pair(self.u, 2, 12, seg1='RNAA', seg2='RNAA')
        assert_almost_equal(val, 30.001, decimal=3)

    def test_major_pair_2(self):
        val = nuclinfo.major_pair(self.u, 5, 9, seg1='RNAA', seg2='RNAA')
        assert_almost_equal(val, 12.581, decimal=3)

    def test_phase_cp_1(self):
        val = nuclinfo.phase_cp(self.u, seg='RNAA', i=9)
        assert_almost_equal(val, 14.846, decimal=3)

    def test_phase_cp_2(self):
        val = nuclinfo.phase_cp(self.u, seg='RNAA', i=21)
        assert_almost_equal(val, 16.199, decimal=3)

    def test_phase_as_1(self):
        val = nuclinfo.phase_as(self.u, seg='RNAA', i=1)
        assert_almost_equal(val, 8.565, decimal=3)

    def test_phase_as_2(self):
        val = nuclinfo.phase_as(self.u, seg='RNAA', i=11)
        assert_almost_equal(val, 151.397, decimal=3)

    def test_tors_1(self):
        val = nuclinfo.tors(self.u, seg='RNAA', i=5)
        assert_allclose(val, [ 296.103016,  179.084427,
                   47.712818,   81.49588 ,  205.350479,
                  284.046562,  198.58165 ], rtol=1e-03)

    def test_tors_2(self):
        val = nuclinfo.tors(self.u, seg='RNAA', i=21)
        assert_allclose(val, [ 300.047634,  172.476593,
                   50.744877,   82.144211,  206.106445,
                  286.208565,  195.652176], rtol=1e-03)

    def test_tors_alpha_1(self):
        val = nuclinfo.tors_alpha(self.u, seg='RNAA', i=6)
        assert_almost_equal(val, 299.278, decimal=3)

    def test_tors_alpha_2(self):
        val = nuclinfo.tors_alpha(self.u, seg='RNAA', i=18)
        assert_almost_equal(val, 304.093, decimal=3)

    def test_tors_beta_1(self):
        val = nuclinfo.tors_beta(self.u, seg='RNAA', i=7)
        assert_almost_equal(val, 177.041, decimal=3)

    def test_tors_beta_2(self):
        val = nuclinfo.tors_beta(self.u, seg='RNAA', i=15)
        assert_almost_equal(val, 162.441, decimal=3)

    def test_tors_gamma_1(self):
        val = nuclinfo.tors_gamma(self.u, seg='RNAA', i=7)
        assert_almost_equal(val, 47.201, decimal=3)

    def test_tors_gamma_2(self):
        val = nuclinfo.tors_gamma(self.u, seg='RNAA', i=15)
        assert_almost_equal(val, 64.276, decimal=3)

    def test_tors_delta_1(self):
        val = nuclinfo.tors_delta(self.u, seg='RNAA', i=7)
        assert_almost_equal(val, 80.223, decimal=3)

    def test_tors_delta_2(self):
        val = nuclinfo.tors_delta(self.u, seg='RNAA', i=15)
        assert_almost_equal(val, 81.553, decimal=3)

    def test_tors_eps_1(self):
        val = nuclinfo.tors_eps(self.u, seg='RNAA', i=7)
        assert_almost_equal(val, 207.476, decimal=3)

    def test_tors_eps_2(self):
        val = nuclinfo.tors_eps(self.u, seg='RNAA', i=15)
        assert_almost_equal(val, 202.352, decimal=3)

    def test_tors_zeta_1(self):
        val = nuclinfo.tors_zeta(self.u, seg='RNAA', i=7)
        assert_almost_equal(val, 289.429, decimal=3)

    def test_tors_zeta_2(self):
        val = nuclinfo.tors_zeta(self.u, seg='RNAA', i=15)
        assert_almost_equal(val, 289.318, decimal=3)

    def test_tors_chi_1(self):
        val = nuclinfo.tors_chi(self.u, seg='RNAA', i=1)
        assert_almost_equal(val, 198.601, decimal=3)

    def test_tors_chi_2(self):
        val = nuclinfo.tors_chi(self.u, seg='RNAA', i=2)
        assert_almost_equal(val, 198.932, decimal=3)

    def test_hydroxyl_1(self):
        val = nuclinfo.hydroxyl(self.u, seg='RNAA', i=20)
        assert_almost_equal(val, 122.689, decimal=3)

    def test_hydroxyl_2(self):
        val = nuclinfo.hydroxyl(self.u, seg='RNAA', i=5)
        assert_almost_equal(val, 122.965, decimal=3)

    def test_pseudo_dihe_baseflip_1(self):
        val = nuclinfo.pseudo_dihe_baseflip(self.u, 16, 2, 3, seg1='RNAA',seg2='RNAA', seg3='RNAA' )
        assert_almost_equal(val, 308.244, decimal=3)

    def test_pseudo_dihe_baseflip_2(self):
        val = nuclinfo.pseudo_dihe_baseflip(self.u, 8, 9, 10, seg1='RNAA',seg2='RNAA', seg3='RNAA')
        assert_almost_equal(val, 25.224, decimal=3)
