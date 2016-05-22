import MDAnalysis as mda
import numpy as np
from six.moves import zip

from numpy.testing import (assert_equal, assert_raises, assert_allclose)

from MDAnalysisTests.datafiles import (DLP_CONFIG, DLP_CONFIG_minimal,
                                       DLP_CONFIG_order, DLP_HISTORY,
                                       DLP_HISTORY_minimal, DLP_HISTORY_order)
from MDAnalysisTests.coordinates.base import BaseTimestepTest


class _DLPConfig(object):
    def setUp(self):
        self.r = mda.coordinates.DLPoly.ConfigReader
        rd = self.rd = self.r(self.f)
        self.ts = rd.ts

    def tearDown(self):
        del self.r
        del self.rd
        del self.ts

    def test_read_unitcell(self):
        ref = np.array([[18.6960000000, 0.0000000000, 0.0000000000
                         ], [0.0000000000, 18.6960000000, 0.0000000000],
                        [0.0000000000, 0.0000000000, 18.6960000000]])
        assert_allclose(self.ts._unitcell, ref)

    def test_positions(self):
        ref = np.array([-7.608595309, -7.897790000, -7.892053559])
        assert_allclose(self.ts._pos[0], ref)

    def test_velocities(self):
        ref = np.array([1.056610291, -1.218664448, 3.345828610])
        assert_allclose(self.ts._velocities[0], ref)

    def test_forces(self):
        ref = np.array([-1979.558687, 739.7961625, 1027.996603])
        assert_allclose(self.ts._forces[0], ref)


class TestConfigReader(_DLPConfig):
    f = DLP_CONFIG

    def test_read(self):
        assert_equal(self.rd.title, "DL_POLY: Potassium Chloride Test Case")


class TestConfigOrder(_DLPConfig):
    f = DLP_CONFIG_order


class TestConfigMinimal(_DLPConfig):
    f = DLP_CONFIG_minimal

    def test_read_unitcell(self):
        pass

    def test_velocities(self):
        assert_raises(AttributeError, getattr, self.ts, "_velocities")

    def test_forces(self):
        assert_raises(AttributeError, getattr, self.ts, "_forces")


class _DLPConfig2(object):
    def setUp(self):
        self.u = mda.Universe(self.f, format='CONFIG')

    def tearDown(self):
        del self.u

    def test_names(self):
        ref = ['C', 'B', 'A']
        assert_equal([a.name for a in self.u.atoms], ref)

    def test_pos(self):
        ref = np.array([-7.821414265, -4.635443539, -4.732164540])
        assert_allclose(self.u.atoms[2].pos, ref)

    def test_vel(self):
        ref = np.array([2.637614561, 0.5778767520E-01, -1.704765568])
        assert_allclose(self.u.atoms[2].velocity, ref)

    def test_for(self):
        ref = np.array([150.3309776, -812.6932914, 1429.413120])
        assert_allclose(self.u.atoms[2].force, ref)

    def test_number(self):
        ref = [0, 1, 2]
        assert_equal([a.index for a in self.u.atoms], ref)


class TestConfigReader2(_DLPConfig2):
    f = DLP_CONFIG_order


class TestConfigReaderMinimal2(_DLPConfig2):
    f = DLP_CONFIG_minimal

    def test_vel(self):
        pass

    def test_for(self):
        pass


class _DLHistory(object):
    def setUp(self):
        self.u = mda.Universe(self.f, format='HISTORY')

    def tearDown(self):
        self.u.trajectory.close()
        del self.u

    def test_len(self):
        assert_equal(len(self.u.trajectory), 3)
        assert_equal([ts.frame for ts in self.u.trajectory], [0, 1, 2])

    def test_getting(self):
        ts = self.u.trajectory[1]
        assert_equal(ts.frame, 1)

    def test_slicing(self):
        nums = [ts.frame for ts in self.u.trajectory[::2]]
        assert_equal(nums, [0, 2])

    def test_slicing_2(self):
        nums = [ts.frame for ts in self.u.trajectory[1::-2]]
        assert_equal(nums, [1])

    def test_position(self):
        ref = np.array([[-7.595541651, -7.898808509, -7.861763110
                         ], [-7.019565641, -7.264933320, -7.045213551],
                        [-6.787470785, -6.912685099, -6.922156843]])
        for ts, r in zip(self.u.trajectory, ref):
            assert_allclose(self.u.atoms[0].pos, r)

    def test_velocity(self):
        ref = np.array([[1.109901682, -1.500264697, 4.752251711
                         ], [-1.398479696, 2.091141311, 1.957430003],
                        [0.2570827995, -0.7146878577, -3.547444215]])
        for ts, r in zip(self.u.trajectory, ref):
            assert_allclose(self.u.atoms[0].velocity, r)

    def test_force(self):
        ref = np.array([[-2621.386432, 1579.334443, 1041.103241
                         ], [-1472.262341, 2450.379615, -8149.916193],
                        [2471.802059, -3828.467296, 3596.679326]])
        for ts, r in zip(self.u.trajectory, ref):
            assert_allclose(self.u.atoms[0].force, r)

    def test_unitcell(self):
        ref1 = np.array([[18.6796195135, 0.0000058913, -0.0000139999
                          ], [0.0000058913, 18.6794658887, -0.0000016255],
                         [-0.0000139999, -0.0000016255, 18.6797229304]])
        ref2 = np.array([[17.2277221163, -0.0044216126, -0.0003229237
                          ], [-0.0044205826, 17.2124253987, 0.0019439244],
                         [-0.0003226531, 0.0019445826, 17.2416976104]])
        ref3 = np.array([[16.5435673205, -0.0108424742, 0.0014935464
                          ], [-0.0108333201, 16.5270298891, 0.0011094612],
                         [0.0014948739, 0.0011058349, 16.5725517831]])
        for ts, r in zip(self.u.trajectory, [ref1, ref2, ref3]):
            assert_allclose(ts._unitcell, r)


class TestDLPolyHistory(_DLHistory):
    f = DLP_HISTORY


class TestDLPolyHistoryOrder(_DLHistory):
    f = DLP_HISTORY_order


class TestDLPolyHistoryMinimal(_DLHistory):
    f = DLP_HISTORY_minimal

    def test_velocity(self):
        assert_raises(mda.NoDataError, getattr, self.u.atoms[0], 'velocity')

    def test_force(self):
        assert_raises(mda.NoDataError, getattr, self.u.atoms[0], 'force')

    def test_unitcell(self):
        pass


class TestDLPolyTimestep(BaseTimestepTest):
    Timestep = mda.coordinates.DLPoly.Timestep
    name = "DLPoly"
    has_box = True
    set_box = True
    unitcell = np.array([[10., 0., 0.],
                         [0., 11., 0.],
                         [0., 0., 12.]])
