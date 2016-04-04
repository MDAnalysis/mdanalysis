from six.moves import range

import MDAnalysis as mda
import numpy as np

from numpy.testing import (assert_equal, assert_almost_equal)
from unittest import TestCase

from MDAnalysisTests.datafiles import (GMS_ASYMOPT, GMS_ASYMSURF, GMS_SYMOPT)


class TestGMSReader(TestCase):
    ''' Test cases for GAMESS output log-files '''

    def setUp(self):
        # optimize no-symmetry
        self.u_aso = mda.Universe(GMS_ASYMOPT,
                                  GMS_ASYMOPT,
                                  format='GMS',
                                  topology_format='GMS')
        self.u_so = mda.Universe(GMS_SYMOPT, GMS_SYMOPT)
        self.u_ass = mda.Universe(GMS_ASYMSURF, GMS_ASYMSURF)

    def test_n_frames(self):
        desired = [21, 8, 10]
        assert_equal(
            self.u_aso.trajectory.n_frames,
            desired[0],
            err_msg="Wrong number of frames read from GAMESS C1 optimization")
        assert_equal(
            self.u_so.trajectory.n_frames,
            desired[1],
            err_msg="Wrong number of frames read from GAMESS D4H optimization")
        assert_equal(
            self.u_ass.trajectory.n_frames,
            desired[2],
            err_msg="Wrong number of frames read from GAMESS C1 surface")

    def test_step5distances_asymopt(self):
        '''TestGMSReader: C1 optimization:
            distance between 1st and 4th atoms changes after 5 steps '''
        desired = -0.0484664
        assert_almost_equal(self.__calcFD(self.u_aso), desired, decimal=5,
                            err_msg="Wrong 1-4 atom distance change after "
                            "5 steps for GAMESS C1 optimization")

    def test_step5distances_symopt(self):
        '''TestGMSReader: Symmetry-input optimization:
            distance between 1st and 4th atoms changes after 5 steps '''
        desired = 0.227637
        assert_almost_equal(self.__calcFD(self.u_so), desired, decimal=5,
                            err_msg="Wrong 1-4 atom distance change after 5 "
                            "steps for GAMESS D4H optimization")

    def test_step5distances_asymsurf(self):
        '''TestGMSReader: Symmetry-input potential-energy surface:
            distance between 1st and 4th atoms changes after 5 steps '''
        desired = -0.499996
        assert_almost_equal(self.__calcFD(self.u_ass), desired, decimal=5,
                            err_msg="Wrong 1-4 atom distance change after 5 "
                            "steps for GAMESS C1 surface")

    def __calcFD(self, u):
        u.trajectory.rewind()
        pp = (u.trajectory.ts._pos[0] - u.trajectory.ts._pos[3])
        z1 = np.sqrt(sum(pp ** 2))
        for i in range(5):
            u.trajectory.next()
        pp = (u.trajectory.ts._pos[0] - u.trajectory.ts._pos[3])
        z2 = np.sqrt(sum(pp ** 2))
        return z1 - z2

    def test_rewind(self):
        self.u_aso.trajectory.rewind()
        assert_equal(self.u_aso.trajectory.ts.frame, 0, "rewinding to frame 0")

    def test_next(self):
        self.u_aso.trajectory.rewind()
        self.u_aso.trajectory.next()
        assert_equal(self.u_aso.trajectory.ts.frame, 1, "loading frame 1")

    def test_dt(self):
        assert_almost_equal(self.u_aso.trajectory.dt,
                            1.0,
                            4,
                            err_msg="wrong timestep dt")

    def tearDown(self):
        del self.u_aso
        del self.u_so
        del self.u_ass
