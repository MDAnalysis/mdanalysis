import MDAnalysis as mda
import numpy as np

from numpy.testing import (assert_equal, assert_array_equal, assert_raises)
from unittest import TestCase

from MDAnalysis.lib.mdamath import triclinic_vectors

from MDAnalysisTests.datafiles import (DMS)
from MDAnalysisTests.coordinates.base import BaseTimestepTest

class TestDMSReader(TestCase):
    def setUp(self):
        self.universe = mda.Universe(DMS)
        self.ts = self.universe.trajectory.ts

    def tearDown(self):
        del self.universe
        del self.ts

    def test_global_cell(self):
        assert_equal(self.ts.dimensions, [0., 0., 0., 0., 0., 0.])

    def test_velocities(self):
        assert_equal(hasattr(self.ts, "_velocities"), False)

    def test_number_of_coords(self):
        # Desired value taken from VMD
        #      Info)    Atoms: 3341
        assert_equal(len(self.universe.atoms), 3341)

    def test_coords_atom_0(self):
        # Desired coordinates taken directly from the SQLite file. Check unit
        # conversion
        coords_0 = np.array([-11.0530004501343,
                             26.6800003051758,
                             12.7419996261597, ],
                            dtype=np.float32)
        assert_array_equal(self.universe.atoms[0].pos, coords_0)

    def test_n_frames(self):
        assert_equal(self.universe.trajectory.n_frames, 1,
                     "wrong number of frames in pdb")

    def test_time(self):
        assert_equal(self.universe.trajectory.time, 0.0,
                     "wrong time of the frame")

    def test_frame(self):
        assert_equal(self.universe.trajectory.frame, 0, "wrong frame number "
                     "(0-based, should be 0 for single frame readers)")

    def test_frame_index_0(self):
        self.universe.trajectory[0]
        assert_equal(self.universe.trajectory.ts.frame, 0,
                     "frame number for frame index 0 should be 0")

    def test_frame_index_1_raises_IndexError(self):
        def go_to_2(traj=self.universe.trajectory):
            traj[1]

        assert_raises(IndexError, go_to_2)


class TestDMSTimestep(BaseTimestepTest):
    Timestep = mda.coordinates.DMS.Timestep
    name = "DMS"
    has_box = True
    unitcell = {'x':np.array([10., 0, 0]),
                'y':np.array([0, 11., 0]),
                'z':np.array([0, 0, 12.])}
    uni_args = (DMS,)

    def test_dimensions_set_box(self):
        self.ts.dimensions = self.newbox
        assert_equal(self.ts.dimensions, self.newbox)
        assert_equal(self.ts._unitcell, self.unitcell)

    def test_set_triclinic_vectors(self):
        ref_vec = triclinic_vectors(self.newbox)
        self.ts.triclinic_dimensions = ref_vec
        assert_equal(self.ts.dimensions, self.newbox)
        assert_equal(self.ts._unitcell, self.unitcell)
