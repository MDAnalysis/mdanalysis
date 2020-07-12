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
from __future__ import print_function, absolute_import
import MDAnalysis
from MDAnalysis.analysis import waterdynamics
from MDAnalysis.lib.correlations import autocorrelation, correct_intermittency

from MDAnalysisTests.datafiles import waterPSF, waterDCD

import pytest
import numpy as np
from mock import patch
from mock import Mock
from numpy.testing import assert_almost_equal, assert_equal

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


def test_HydrogenBondLifetimes_deprecation(universe):
    with pytest.deprecated_call():
        hbl = waterdynamics.HydrogenBondLifetimes(
            universe, SELECTION1, SELECTION1, 0, 5, 3)


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


def test_SurvivalProbability_intermittency1and2(universe):
    """
    Intermittency of 2 means that we still count an atom if it is not present for up to 2 consecutive frames,
    but then returns at the following step.
    """
    with patch.object(universe, 'select_atoms') as select_atoms_mock:
        ids = [(9, 8), (), (8,), (9,), (8,), (), (9, 8), (), (8,), (9, 8)]
        select_atoms_mock.side_effect = lambda selection: Mock(ids=ids.pop())   # atom IDs fed set by set
        sp = waterdynamics.SurvivalProbability(universe, "")
        sp.run(tau_max=3, stop=10, verbose=True, intermittency=2)
        assert all(x == {9, 8} for x in sp._intermittent_selected_ids)
        assert_almost_equal(sp.sp_timeseries, [1, 1, 1, 1])


def test_SurvivalProbability_intermittency2lacking(universe):
    """
    If an atom is not present for more than 2 consecutive frames,
    it is considered to have left the region.
    """
    with patch.object(universe, 'select_atoms') as select_atoms_mock:
        ids = [(9,), (), (), (), (9,), (), (), (), (9,)]
        select_atoms_mock.side_effect = lambda selection: Mock(ids=ids.pop())   # atom IDs fed set by set
        sp = waterdynamics.SurvivalProbability(universe, "")
        sp.run(tau_max=3, stop=9, verbose=True, intermittency=2)
        assert_almost_equal(sp.sp_timeseries, [1, 0, 0, 0])


def test_SurvivalProbability_intermittency1_step5_noSkipping(universe):
    """
    Step leads to skipping frames if (tau_max + 1) + (intermittency * 2) < step.
    No frames should be skipped.
    """
    with patch.object(universe, 'select_atoms') as select_atoms_mock:
        ids = [(2, 3), (3,), (2, 3), (3,), (2,), (3,), (2, 3), (3,), (2, 3), (2, 3)]
        select_atoms_mock.side_effect = lambda selection: Mock(ids=ids.pop())   # atom IDs fed set by set
        sp = waterdynamics.SurvivalProbability(universe, "")
        sp.run(tau_max=2, stop=10, verbose=True, intermittency=1, step=5)
        assert all((x == {2, 3} for x in sp._intermittent_selected_ids))
        assert_almost_equal(sp.sp_timeseries, [1, 1, 1])


def test_SurvivalProbability_intermittency1_step5_Skipping(universe):
    """
    Step leads to skipping frames if (tau_max + 1) * (intermittency * 2) < step.
    In this case one frame will be skipped per window.
    """
    with patch.object(universe, 'select_atoms') as select_atoms_mock:
        ids = [(1,), (), (1,), (), (1,), (), (1,), (), (1,), (1,)]
        beforepopsing = len(ids) - 2
        select_atoms_mock.side_effect = lambda selection: Mock(ids=ids.pop())   # atom IDs fed set by set
        sp = waterdynamics.SurvivalProbability(universe, "")
        sp.run(tau_max=1, stop=10, verbose=True, intermittency=1, step=5)
        assert all((x == {1} for x in sp._intermittent_selected_ids))
        assert len(sp._selected_ids) == beforepopsing
        assert_almost_equal(sp.sp_timeseries, [1, 1])


def test_intermittency_none():
    # No changes asked - returns the same data
    input_ids = [{1}, {1}, {1}, set(), set(), {1}, {1}, {1}, set(), set(), {1}, set(), set(), {1}]
    corrected = correct_intermittency(input_ids, intermittency=0)
    assert all(x == y for x,y in zip(input_ids, corrected))


def test_intermittency_1and2():
    # The maximum gap in the dataset is 2, so the IDs are always present after correction
    input_ids = [{9, 8}, set(), {8, }, {9, }, {8, }, set(), {9, 8}, set(), {8, }, {9, 8, }]
    corrected = correct_intermittency(input_ids, intermittency=2)
    assert all((x == {9, 8} for x in corrected))


def test_intermittency_2tooShort():
    #The IDs are abscent for too long/
    input_ids = [{9,}, {}, {}, {}, {9,}, {}, {}, {}, {9,}]
    corrected = correct_intermittency(input_ids, intermittency=2)
    assert all(x == y for x, y in zip(input_ids, corrected))


def test_intermittency_setsOfSets():
    # Verificaiton for the case of hydrogen bonds (sets of sets)
    input_ids = [{frozenset({1,2}), frozenset({3, 4})},set(), set(),
                 {frozenset({1, 2}), frozenset({3, 4})}, set(), set(),
                 {frozenset({1, 2}), frozenset({3, 4})}, set(), set(),
                 {frozenset({1, 2}), frozenset({3, 4})}]
    corrected = correct_intermittency(input_ids, intermittency=2)
    assert all((x == {frozenset({1, 2}), frozenset({3, 4})} for x in corrected))


def test_autocorrelation_alwaysPresent():
    input = [{1, 2}, {1, 2}, {1, 2}, {1, 2}, {1, 2}, {1, 2}, {1, 2}]
    tau_timeseries, sp_timeseries, sp_timeseries_data = autocorrelation(input, tau_max=3)
    assert all(np.equal(sp_timeseries, 1))


def test_autocorrelation_definedTaus():
    input_ids = [{9, 8, 7}, {8, 7, 6}, {7, 6, 5}, {6, 5, 4}, {5, 4, 3}, {4, 3, 2}, {3, 2, 1}]
    tau_timeseries, sp_timeseries, sp_timeseries_data = autocorrelation(input_ids, tau_max=3)
    assert_almost_equal(sp_timeseries, [1, 2/3., 1/3., 0])


def test_autocorrelation_intermittency1_windowJump_intermittencyAll():
    """
    Step leads to skipping frames if (tau_max + 1) + (intermittency * 2) < step.
    No frames should be skipped so intermittency should be applied to all.
    """
    input_ids = [{2, 3}, {3,}, {2, 3}, {3,}, {2,}, {3,}, {2, 3}, {3,}, {2, 3}, {2, 3}]
    corrected = correct_intermittency(input_ids, intermittency=1)
    tau_timeseries, sp_timeseries, sp_timeseries_data = autocorrelation(corrected, tau_max=2,
                                                                        window_step=5)
    assert all((x == {2, 3} for x in corrected))
    assert_almost_equal(sp_timeseries, [1, 1, 1])


def test_autocorrelation_windowBigJump():
    #The empty sets are ignored (no intermittency)
    input_ids = [{1}, {1}, {1}, set(), set(), {1}, {1}, {1}, set(), set(), {1}, {1}, {1}]
    tau_timeseries, sp_timeseries, sp_timeseries_data = autocorrelation(input_ids, tau_max=2, window_step=5)
    assert_almost_equal(sp_timeseries, [1, 1, 1])


def test_autocorrelation_windowBigJump_absence():
    # In the last frame the molecules are absent
    input_ids = [{1}, {1}, {1}, set(), set(), {1}, {1}, {1}, set(), set(), {1}, set(), set()]
    tau_timeseries, sp_timeseries, sp_timeseries_data = autocorrelation(input_ids, tau_max=2, window_step=5)
    assert_almost_equal(sp_timeseries, [1, 2/3., 2/3.])


def test_autocorrelation_intermittency1_many():
    input_ids = [{1}, set(), {1}, set(), {1}, set(), {1}, set(), {1}, set(), {1}, set(), {1}, set(), {1}]
    corrected = correct_intermittency(input_ids, intermittency=1)
    tau_timeseries, sp_timeseries, sp_timeseries_data = autocorrelation(corrected, tau_max=14,
                                                                        window_step=5)
    assert_almost_equal(sp_timeseries, [1] * 15)


def test_autocorrelation_intermittency2_windowBigJump():
    # The intermittency corrects the last frame
    input_ids = [{1}, {1}, {1}, set(), set(), {1}, {1}, {1}, set(), set(), {1}, set(), set(), {1}]
    corrected = correct_intermittency(input_ids, intermittency=2)
    tau_timeseries, sp_timeseries, sp_timeseries_data = autocorrelation(corrected, tau_max=2,
                                                                        window_step=5)
    assert_almost_equal(sp_timeseries, [1, 1, 1])


def test_SurvivalProbability_t0tf(universe):
    with patch.object(universe, 'select_atoms') as select_atoms_mock:
        ids = [(0, ), (0, ), (7, 6, 5), (6, 5, 4), (5, 4, 3), (4, 3, 2), (3, 2, 1), (0, )]
        select_atoms_mock.side_effect = lambda selection: Mock(ids=ids.pop(2))   # atom IDs fed set by set
        sp = waterdynamics.SurvivalProbability(universe, "")
        sp.run(tau_max=3, start=2, stop=7)
        assert_almost_equal(sp.sp_timeseries, [1, 2 / 3.0, 1 / 3.0, 0])


def test_SurvivalProbability_definedTaus(universe):
    with patch.object(universe, 'select_atoms') as select_atoms_mock:
        ids = [(9, 8, 7), (8, 7, 6), (7, 6, 5), (6, 5, 4), (5, 4, 3), (4, 3, 2), (3, 2, 1)]
        select_atoms_mock.side_effect = lambda selection: Mock(ids=ids.pop())   # atom IDs fed set by set
        sp = waterdynamics.SurvivalProbability(universe, "")
        sp.run(tau_max=3, start=0, stop=7, verbose=True)
        assert_almost_equal(sp.sp_timeseries, [1, 2 / 3.0, 1 / 3.0, 0])


def test_SurvivalProbability_zeroMolecules(universe):
    # no atom IDs found
    with patch.object(universe, 'select_atoms', return_value=Mock(ids=[])) as select_atoms_mock:
        sp = waterdynamics.SurvivalProbability(universe, "")
        sp.run(tau_max=3, start=3, stop=7, verbose=True)
        assert all(np.isnan(sp.sp_timeseries[1:]))


def test_SurvivalProbability_alwaysPresent(universe):
    # always the same atom IDs found, 7 and 8
    with patch.object(universe, 'select_atoms', return_value=Mock(ids=[7, 8])) as select_atoms_mock:
        sp = waterdynamics.SurvivalProbability(universe, "")
        sp.run(tau_max=3, start=0, stop=7, verbose=True)
        assert all(np.equal(sp.sp_timeseries, 1))


def test_SurvivalProbability_stepLargerThanDtmax(universe):
    # Testing if the frames are skipped correctly
    with patch.object(universe, 'select_atoms', return_value=Mock(ids=(1,))) as select_atoms_mock:
        sp = waterdynamics.SurvivalProbability(universe, "")
        sp.run(tau_max=2, step=5, stop=10, verbose=True)
        assert_equal(sp.sp_timeseries, [1, 1, 1])
        # with tau_max=2 for all the frames we only read 6 of them
        # this is because the frames which are not used are skipped, and therefore 'select_atoms'
        assert universe.trajectory.n_frames > 6
        assert_equal(select_atoms_mock.call_count, 6)


def test_SurvivalProbability_stepEqualDtMax(universe):
    with patch.object(universe, 'select_atoms', return_value=Mock(ids=(1,))) as select_atoms_mock:
        sp = waterdynamics.SurvivalProbability(universe, "")
        sp.run(tau_max=4, step=5, stop=10, verbose=True)
        # all frames from 0, with 9 inclusive
        assert_equal(select_atoms_mock.call_count, 10)
