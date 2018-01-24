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
from __future__ import absolute_import

import MDAnalysis as mda
import pytest
from MDAnalysis.analysis import nuclinfo
from MDAnalysisTests.datafiles import NUCL
from numpy.testing import (
    assert_almost_equal,
    assert_allclose,
)


@pytest.fixture(scope='module')
def u():
    return mda.Universe(NUCL)


@pytest.mark.parametrize('i, bp, seg1, seg2, expected_value', (
        (1, 2, 'RNAA', 'RNAA', 4.449),
        (22, 23, 'RNAA', 'RNAA', 4.601)
))
def test_wc_pair(u, i, bp, seg1, seg2, expected_value):
    val = nuclinfo.wc_pair(u, i, bp, seg1=seg1, seg2=seg2)
    assert_almost_equal(val, expected_value, decimal=3)


@pytest.mark.parametrize('i, bp, seg1, seg2, expected_value', (
        (3, 17, 'RNAA', 'RNAA', 16.649),
        (20, 5, 'RNAA', 'RNAA', 4.356)
))
def test_minor_pair(u, i, bp, seg1, seg2, expected_value):
    val = nuclinfo.minor_pair(u, i, bp, seg1=seg1, seg2=seg2)
    assert_almost_equal(val, expected_value, decimal=3)


@pytest.mark.parametrize('i, bp, seg1, seg2, expected_value', (
        (2, 12, 'RNAA', 'RNAA', 30.001),
        (5, 9, 'RNAA', 'RNAA', 12.581)
))
def test_major_pair(u, i, bp, seg1, seg2, expected_value):
    val = nuclinfo.major_pair(u, i, bp, seg1=seg1, seg2=seg2)
    assert_almost_equal(val, expected_value, decimal=3)


@pytest.mark.parametrize('seg, i, expected_value', (
        ('RNAA', 9, 14.846),
        ('RNAA', 21, 16.199)
))
def test_phase_cp(u, seg, i, expected_value):
    val = nuclinfo.phase_cp(u, seg=seg, i=i)
    assert_almost_equal(val, expected_value, decimal=3)


@pytest.mark.parametrize('seg, i, expected_value', (
        ('RNAA', 1, 8.565),
        ('RNAA', 11, 151.397)
))
def test_phase_as(u, seg, i, expected_value):
    val = nuclinfo.phase_as(u, seg=seg, i=i)
    assert_almost_equal(val, expected_value, decimal=3)


@pytest.mark.parametrize('seg, i, expected_value', (
        ('RNAA', 5, [296.103016, 179.084427, 47.712818, 81.49588, 205.350479,
                     284.046562, 198.58165]),
        ('RNAA', 21, [300.047634, 172.476593, 50.744877, 82.144211, 206.106445,
                      286.208565, 195.652176])
))
def test_tors(u, seg, i, expected_value):
    val = nuclinfo.tors(u, seg=seg, i=i)
    assert_allclose(val, expected_value, rtol=1e-03)


@pytest.mark.parametrize('seg, i, expected_value', (
        ('RNAA', 6, 299.278),
        ('RNAA', 18, 304.093)
))
def test_tors_alpha(u, seg, i, expected_value):
    val = nuclinfo.tors_alpha(u, seg=seg, i=i)
    assert_almost_equal(val, expected_value, decimal=3)


@pytest.mark.parametrize('seg, i, expected_value', (
        ('RNAA', 7, 177.041),
        ('RNAA', 15, 162.441)
))
def test_tors_beta(u, seg, i, expected_value):
    val = nuclinfo.tors_beta(u, seg=seg, i=i)
    assert_almost_equal(val, expected_value, decimal=3)


@pytest.mark.parametrize('seg, i, expected_value', (
        ('RNAA', 7, 47.201),
        ('RNAA', 15, 64.276)
))
def test_tors_gamma(u, seg, i, expected_value):
    val = nuclinfo.tors_gamma(u, seg=seg, i=i)
    assert_almost_equal(val, expected_value, decimal=3)


@pytest.mark.parametrize('seg, i, expected_value', (
        ('RNAA', 7, 80.223),
        ('RNAA', 15, 81.553)
))
def test_tors_delta(u, seg, i, expected_value):
    val = nuclinfo.tors_delta(u, seg=seg, i=i)
    assert_almost_equal(val, expected_value, decimal=3)


@pytest.mark.parametrize('seg, i, expected_value', (
        ('RNAA', 7, 207.476),
        ('RNAA', 15, 202.352)
))
def test_tors_eps(u, seg, i, expected_value):
    val = nuclinfo.tors_eps(u, seg=seg, i=i)
    assert_almost_equal(val, expected_value, decimal=3)


@pytest.mark.parametrize('seg, i, expected_value', (
        ('RNAA', 7, 289.429),
        ('RNAA', 15, 289.318)
))
def test_tors_zeta(u, seg, i, expected_value):
    val = nuclinfo.tors_zeta(u, seg=seg, i=i)
    assert_almost_equal(val, expected_value, decimal=3)


@pytest.mark.parametrize('seg, i, expected_value', (
        ('RNAA', 1, 198.601),
        ('RNAA', 2, 198.932)
))
def test_tors_chi(u, seg, i, expected_value):
    val = nuclinfo.tors_chi(u, seg=seg, i=i)
    assert_almost_equal(val, expected_value, decimal=3)


@pytest.mark.parametrize('seg, i, expected_value', (
        ('RNAA', 20, 122.689),
        ('RNAA', 5, 122.965)
))
def test_hydroxyl(u, seg, i, expected_value):
    val = nuclinfo.hydroxyl(u, seg=seg, i=i)
    assert_almost_equal(val, expected_value, decimal=3)


@pytest.mark.parametrize('bp1, bp2, i, seg1, seg2, seg3, expected_value', (
        (16, 2, 3, 'RNAA', 'RNAA', 'RNAA', 308.244),
        (8, 9, 10, 'RNAA', 'RNAA', 'RNAA', 25.224),
))
def test_pseudo_dihe_baseflip(u, bp1, bp2, i, seg1, seg2, seg3, expected_value):
    val = nuclinfo.pseudo_dihe_baseflip(u, bp1, bp2, i, seg1, seg2, seg3)
    assert_almost_equal(val, expected_value, decimal=3)
