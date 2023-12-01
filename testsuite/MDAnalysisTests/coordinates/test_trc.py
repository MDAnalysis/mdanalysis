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
import pytest
import numpy as np
from numpy.testing import assert_allclose

import MDAnalysis as mda
from MDAnalysis.coordinates.TRC import TRCReader
from MDAnalysisTests.datafiles import TRC_PDB_VAC, TRC_TRAJ1_VAC, TRC_TRAJ2_VAC
from MDAnalysisTests.datafiles import TRC_CLUSTER_VAC
from MDAnalysisTests.datafiles import TRC_PDB_SOLV, TRC_TRAJ_SOLV


class TestTRCReaderVacuumBox:
    @pytest.fixture(scope='class')
    def TRC_U(self):
        return mda.Universe(TRC_PDB_VAC, [TRC_TRAJ1_VAC, TRC_TRAJ2_VAC],
                            continuous=True)

    def test_initial_frame_is_0(self, TRC_U):
        assert TRC_U.trajectory.ts.frame == 0

    def test_trc_positions(self, TRC_U):
        # first frame first particle
        TRC_U.trajectory[0]
        assert_allclose(TRC_U.atoms.positions[0],
                        [2.19782507, 24.65064345, 29.39783426])
        # fith frame first particle
        TRC_U.trajectory[4]
        assert_allclose(TRC_U.atoms.positions[0],
                        [0.37026654, 22.78805010, 3.69695262])

    def test_trc_dimensions(self, TRC_U):
        assert TRC_U.trajectory[0].dimensions is None

    def test_trc_n_frames(self, TRC_U):
        assert len(TRC_U.trajectory) == 6
        assert (TRC_U.trajectory.n_frames) == 6

    def test_trc_n_atoms(self, TRC_U):
        assert (TRC_U.trajectory.n_atoms) == 73

    def test_trc_frame(self, TRC_U):
        assert TRC_U.trajectory[0].frame == 0
        assert TRC_U.trajectory[4].frame == 4

    def test_trc_time(self, TRC_U):
        assert TRC_U.trajectory[0].time == 0
        assert TRC_U.trajectory[4].time == 80

    def test_trc_dt(self, TRC_U):
        time_array = np.array([ts.time for ts in TRC_U.trajectory])
        assert_allclose(time_array, np.arange(6)*20.0)

        dt_array = np.diff(time_array)
        assert_allclose(dt_array, np.full(5, 20.0))

    def test_trc_data_step(self, TRC_U):
        assert TRC_U.trajectory[0].data['step'] == 0
        assert TRC_U.trajectory[4].data['step'] == 10000

    def test_periodic(self, TRC_U):
        assert TRC_U.trajectory.periodic is False

    def test_rewind(self, TRC_U):
        TRC_U.trajectory[0]
        trc = TRC_U.trajectory
        trc.next()
        trc.next()
        trc.next()
        trc.next()
        assert trc.ts.frame == 4, "trajectory.next() did not forward to frameindex 4"
        trc.rewind()
        assert trc.ts.frame == 0, "trajectory.rewind() failed to rewind to first frame"

        assert np.any(TRC_U.atoms.positions != 0), "The atom positions " \
                                                   "are not populated"

    def test_random_access(self, TRC_U):
        TRC_U.trajectory[0]
        pos0 = TRC_U.atoms.positions
        TRC_U.trajectory.next()
        TRC_U.trajectory.next()
        pos2 = TRC_U.atoms.positions

        TRC_U.trajectory[0]
        assert_allclose(TRC_U.atoms.positions, pos0)

        TRC_U.trajectory[2]
        assert_allclose(TRC_U.atoms.positions, pos2)

    def test_read_frame_reopens(self, TRC_U):
        TRC_U.trajectory._reopen()
        TRC_U.trajectory[4]
        assert TRC_U.trajectory.ts.frame == 4


class TestTRCReaderSolvatedBox:
    @pytest.fixture(scope='class')
    def TRC_U(self):
        return mda.Universe(TRC_PDB_SOLV, TRC_TRAJ_SOLV)

    def test_trc_n_atoms(self, TRC_U):
        assert TRC_U.trajectory.n_atoms == 2797

    def test_periodic(self, TRC_U):
        assert TRC_U.trajectory.periodic is True

    def test_trc_dimensions(self, TRC_U):
        ts = TRC_U.trajectory[1]
        assert_allclose(
            ts.dimensions,
            [30.54416298, 30.54416298, 30.54416298, 90., 90., 90.]
        )

    def test_open_twice(self, TRC_U):
        TRC_U.trajectory._reopen()
        with pytest.raises(IOError):
            TRC_U.trajectory.open_trajectory()


class TestTRCReaderClusterTrajectory:
    @pytest.fixture(scope='class')
    def TRC_U(self):
        with pytest.warns(UserWarning) as record_of_warnings:
            TRC_U = mda.Universe(TRC_PDB_VAC, TRC_CLUSTER_VAC, format="TRC")

        warning_strings = ["The trajectory does not contain TIMESTEP blocks!",
                           "POSITION block is not supported!"]

        warning_strings_found = 0
        for w in record_of_warnings:
            if (str(w.message) in warning_strings):
                warning_strings_found += 1

        assert warning_strings_found == 2
        return TRC_U

    def test_trc_n_atoms(self, TRC_U):
        assert TRC_U.trajectory.n_atoms == 73

    def test_periodic(self, TRC_U):
        assert TRC_U.trajectory.periodic is False

    def test_trc_dimensions(self, TRC_U):
        ts = TRC_U.trajectory[1]
        assert ts.dimensions is None

    def test_trc_n_frames(self, TRC_U):
        assert len(TRC_U.trajectory) == 3
        assert TRC_U.trajectory.n_frames == 3

    def test_trc_frame(self, TRC_U):
        with pytest.warns(UserWarning, match="POSITION block is not supported!"):
            assert TRC_U.trajectory[0].frame == 0
            assert TRC_U.trajectory[2].frame == 2

    def test_trc_time(self, TRC_U):
        with pytest.warns(UserWarning, match="POSITION block is not supported!"):
            assert TRC_U.trajectory[0].time == 0
            assert TRC_U.trajectory[2].time == 0
