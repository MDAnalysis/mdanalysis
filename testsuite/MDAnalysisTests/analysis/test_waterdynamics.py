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
from __future__ import print_function, absolute_import
import MDAnalysis
from MDAnalysis.analysis import waterdynamics

from MDAnalysisTests.datafiles import waterPSF, waterDCD
from MDAnalysisTests.datafiles import PDB, XTC

import pytest
import numpy as np
from mock import patch
from mock import Mock
from numpy.testing import assert_almost_equal

SELECTION1 = "byres name OH2"
SELECTION2 = "byres name P1"


@pytest.fixture(scope='module')
def universe():
    return MDAnalysis.Universe(waterPSF, waterDCD)


def test_HydrogenBondLifetimes(universe):
    hbl = waterdynamics.HydrogenBondLifetimes(
        universe, SELECTION1, SELECTION1, 0, 5, 3)
    hbl.run()
    assert_almost_equal(hbl.timeseries[2][1], 0.75, 5)


def test_WaterOrientationalRelaxation(universe):
    wor = waterdynamics.WaterOrientationalRelaxation(
        universe, SELECTION1, 0, 5, 2)
    wor.run()
    assert_almost_equal(wor.timeseries[1][2], 0.35887,
                        decimal=5)


def test_WaterOrientationalRelaxation_zeroMolecules(universe):
    wor_zero = waterdynamics.WaterOrientationalRelaxation(
        universe, SELECTION2, 0, 5, 2)
    wor_zero.run()
    assert_almost_equal(wor_zero.timeseries[1], (0.0, 0.0, 0.0))


def test_AngularDistribution(universe):
    ad = waterdynamics.AngularDistribution(universe, SELECTION1, 40)
    ad.run()
    # convert a string with two "floats" into a float array
    result = np.array(ad.graph[0][39].split(), dtype=np.float64)
    assert_almost_equal(result, (0.951172947884, 0.48313682125))


def test_MeanSquareDisplacement(universe):
    msd = waterdynamics.MeanSquareDisplacement(universe, SELECTION1, 0, 10, 2)
    msd.run()
    assert_almost_equal(msd.timeseries[1], 0.03984,
                        decimal=5)


def test_MeanSquareDisplacement_zeroMolecules(universe):
    msd_zero = waterdynamics.MeanSquareDisplacement(
        universe, SELECTION2, 0, 10, 2)
    msd_zero.run()
    assert_almost_equal(msd_zero.timeseries[1], 0.0)


def test_SurvivalProbability_t0tf(universe):
    with patch.object(universe, 'select_atoms') as select_atoms_mock:
        ids = [(0, ), (0, ), (7, 6, 5), (6, 5, 4), (5, 4, 3), (4, 3, 2), (3, 2, 1), (0, )]
        select_atoms_mock.side_effect = lambda selection: Mock(ids=ids.pop(2))   # atom IDs fed set by set
        sp = waterdynamics.SurvivalProbability(universe, "")
        sp.run(tau_max=3, start=2, stop=6)
        assert_almost_equal(sp.sp_timeseries, [2 / 3.0, 1 / 3.0, 0])


def test_SurvivalProbability_definedTaus(universe):
    with patch.object(universe, 'select_atoms') as select_atoms_mock:
        ids = [(9, 8, 7), (8, 7, 6), (7, 6, 5), (6, 5, 4), (5, 4, 3), (4, 3, 2), (3, 2, 1)]
        select_atoms_mock.side_effect = lambda selection: Mock(ids=ids.pop())   # atom IDs fed set by set
        sp = waterdynamics.SurvivalProbability(universe, "")
        sp.run(tau_max=3, start=0, stop=6)
        assert_almost_equal(sp.sp_timeseries, [2 / 3.0, 1 / 3.0, 0])


def test_SurvivalProbability_zeroMolecules(universe):
    with patch.object(universe, 'select_atoms') as select_atoms_mock:
        # no atom IDs found
        select_atoms_mock.return_value = Mock(ids=[])
        sp = waterdynamics.SurvivalProbability(universe, "")
        sp.run(tau_max=3, start=3, stop=6)
        assert all(np.isnan(sp.sp_timeseries))


def test_SurvivalProbability_alwaysPresent(universe):
    with patch.object(universe, 'select_atoms') as select_atoms_mock:
        # always the same atom IDs found
        select_atoms_mock.return_value = Mock(ids=[7, 8])
        sp = waterdynamics.SurvivalProbability(universe, "")
        sp.run(tau_max=3, start=0, stop=6)
        assert all(np.equal(sp.sp_timeseries, 1))