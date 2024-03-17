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
import MDAnalysis as mda
from MDAnalysis.lib.mdamath import triclinic_box
import numpy as np
import pytest

from numpy.testing import (assert_equal, assert_allclose)

from MDAnalysisTests.datafiles import (DLP_CONFIG, DLP_CONFIG_minimal,
                                       DLP_CONFIG_order, DLP_HISTORY,
                                       DLP_HISTORY_minimal, DLP_HISTORY_order,
                                       DLP_HISTORY_minimal_cell)


class _DLPConfig(object):
    @pytest.fixture()
    def r(self):
        return mda.coordinates.DLPoly.ConfigReader

    @pytest.fixture()
    def rd(self, r):
        return r(self.f)

    @pytest.fixture()
    def ts(self, rd):
        return rd.ts

    def test_read_unitcell(self, ts):
        ref = np.array([[18.6960000000, 0.0000000000, 0.0000000000],
                        [0.0000000000, 18.6960000000, 0.0000000000],
                        [0.0000000000, 0.0000000000, 18.6960000000]])
        assert_allclose(ts.dimensions, triclinic_box(*ref))

    def test_positions(self, ts):
        ref = np.array([-7.608595309, -7.897790000, -7.892053559])
        assert_allclose(ts._pos[0], ref)

    def test_velocities(self, ts):
        ref = np.array([1.056610291, -1.218664448, 3.345828610])
        assert_allclose(ts._velocities[0], ref)

    def test_forces(self, ts):
        ref = np.array([-1979.558687, 739.7961625, 1027.996603])
        assert_allclose(ts._forces[0], ref)


class TestConfigReader(_DLPConfig):
    f = DLP_CONFIG

    def test_read(self, rd):
        assert_equal(rd.title, "DL_POLY: Potassium Chloride Test Case")


class TestConfigOrder(_DLPConfig):
    f = DLP_CONFIG_order


class TestConfigMinimal(_DLPConfig):
    f = DLP_CONFIG_minimal

    def test_read_unitcell(self):
        pass

    # cythonised class can no longer raise AttributeError
    # so changed to test of has_... properties
    def test_velocities(self, ts):
        assert(ts.has_velocities == False)

    def test_forces(self, ts):
        assert(ts.has_forces == False)


class _DLPConfig2(object):
    @pytest.fixture()
    def u(self):
        return mda.Universe(self.f, format='CONFIG')

    def test_names(self, u):
        ref = ['C', 'B', 'A']
        assert_equal([a.name for a in u.atoms], ref)

    def test_pos(self, u):
        ref = np.array([-7.821414265, -4.635443539, -4.732164540])
        assert_allclose(u.atoms[2].position, ref)

    def test_vel(self, u):
        ref = np.array([2.637614561, 0.5778767520E-01, -1.704765568])
        assert_allclose(u.atoms[2].velocity, ref)

    def test_for(self, u):
        ref = np.array([150.3309776, -812.6932914, 1429.413120])
        assert_allclose(u.atoms[2].force, ref)

    def test_number(self, u):
        ref = [0, 1, 2]
        assert_equal([a.index for a in u.atoms], ref)


class TestConfigReader2(_DLPConfig2):

    f = DLP_CONFIG_order


class TestConfigReaderMinimal2(_DLPConfig2):

    f = DLP_CONFIG_minimal

    def test_vel(self):
        pass

    def test_for(self):
        pass


class _DLHistory(object):
    @pytest.fixture()
    def u(self):
        return mda.Universe(self.f, format='HISTORY')


    def test_len(self, u):
        assert_equal(len(u.trajectory), 3)
        assert_equal([ts.frame for ts in u.trajectory], [0, 1, 2])

    def test_getting(self, u):
        ts = u.trajectory[1]
        assert_equal(ts.frame, 1)

    def test_slicing(self, u):
        nums = [ts.frame for ts in u.trajectory[::2]]
        assert_equal(nums, [0, 2])

    def test_slicing_2(self, u):
        nums = [ts.frame for ts in u.trajectory[1::-2]]
        assert_equal(nums, [1])

    def test_position(self, u):
        ref = np.array([[-7.595541651, -7.898808509, -7.861763110
                         ], [-7.019565641, -7.264933320, -7.045213551],
                        [-6.787470785, -6.912685099, -6.922156843]])
        for ts, r in zip(u.trajectory, ref):
            assert_allclose(u.atoms[0].position, r)

    def test_velocity(self, u):
        ref = np.array([[1.109901682, -1.500264697, 4.752251711
                         ], [-1.398479696, 2.091141311, 1.957430003],
                        [0.2570827995, -0.7146878577, -3.547444215]])
        for ts, r in zip(u.trajectory, ref):
            assert_allclose(u.atoms[0].velocity, r)

    def test_force(self, u):
        ref = np.array([[-2621.386432, 1579.334443, 1041.103241
                         ], [-1472.262341, 2450.379615, -8149.916193],
                        [2471.802059, -3828.467296, 3596.679326]])
        for ts, r in zip(u.trajectory, ref):
            assert_allclose(u.atoms[0].force, r)

    def test_unitcell(self, u):
        ref1 = np.array([[18.6796195135, 0.0000058913, -0.0000139999
                          ], [0.0000058913, 18.6794658887, -0.0000016255],
                         [-0.0000139999, -0.0000016255, 18.6797229304]])
        ref2 = np.array([[17.2277221163, -0.0044216126, -0.0003229237
                          ], [-0.0044205826, 17.2124253987, 0.0019439244],
                         [-0.0003226531, 0.0019445826, 17.2416976104]])
        ref3 = np.array([[16.5435673205, -0.0108424742, 0.0014935464
                          ], [-0.0108333201, 16.5270298891, 0.0011094612],
                         [0.0014948739, 0.0011058349, 16.5725517831]])
        for ts, r in zip(u.trajectory, [ref1, ref2, ref3]):
            assert_allclose(ts.dimensions, triclinic_box(*r))


class TestDLPolyHistory(_DLHistory):
    f = DLP_HISTORY


class TestDLPolyHistoryOrder(_DLHistory):
    f = DLP_HISTORY_order


class TestDLPolyHistoryMinimal(_DLHistory):
    f = DLP_HISTORY_minimal

    def test_velocity(self, u):
        with pytest.raises(mda.NoDataError):
            getattr(u.atoms[0], 'velocity')

    def test_force(self, u):
        with pytest.raises(mda.NoDataError):
            getattr(u.atoms[0], 'force')

    def test_unitcell(self):
        pass


class TestDLPolyHistoryMinimalCell(_DLHistory):
    f = DLP_HISTORY_minimal_cell

    def test_velocity(self, u):
        with pytest.raises(mda.NoDataError):
            getattr(u.atoms[0], 'velocity')

    def test_force(self, u):
        with pytest.raises(mda.NoDataError):
            getattr(u.atoms[0], 'force')
