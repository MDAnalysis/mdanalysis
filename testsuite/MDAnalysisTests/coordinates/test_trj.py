import MDAnalysis as mda
import numpy as np

from numpy.testing import (assert_equal, assert_,
                           assert_almost_equal, assert_raises)

from MDAnalysisTests.coordinates.reference import RefACHE, RefCappedAla
from MDAnalysisTests.datafiles import (PRM, TRJ, TRJ_bz2, PRMpbc, TRJpbc_bz2)
from MDAnalysisTests.coordinates.base import BaseTimestepTest


class _TRJReaderTest(object):
    # use as a base class (override setUp()) and mixin a reference
    def tearDown(self):
        del self.universe

    def test_load_prm(self):
        U = self.universe
        assert_equal(len(U.atoms), self.ref_n_atoms,
                     "load Universe from PRM and TRJ")

    def test_n_atoms(self):
        assert_equal(self.universe.trajectory.n_atoms, self.ref_n_atoms,
                     "wrong number of atoms")

    def test_n_frames(self):
        assert_equal(self.universe.trajectory.n_frames, self.ref_n_frames,
                     "wrong number of frames in xyz")

    def test_periodic(self):
        assert_equal(self.universe.trajectory.periodic, self.ref_periodic)

    def test_amber_proteinselection(self):
        protein = self.universe.select_atoms('protein')
        assert_equal(protein.n_atoms, self.ref_proteinatoms,
                     "error in protein selection (HIS or termini?)")

    def test_sum_centres_of_geometry(self):
        protein = self.universe.select_atoms('protein')
        total = np.sum([protein.center_of_geometry() for ts in
                        self.universe.trajectory])
        assert_almost_equal(total, self.ref_sum_centre_of_geometry, self.prec,
                            err_msg="sum of centers of geometry over the "
                            "trajectory do not match")

    def test_initial_frame_is_0(self):
        assert_equal(self.universe.trajectory.ts.frame, 0,
                     "initial frame is not 0 but {0}".format(
                         self.universe.trajectory.ts.frame))

    def test_starts_with_first_frame(self):
        """Test that coordinate arrays are filled as soon as the trajectory
        has been opened."""
        assert_(np.any(self.universe.atoms.positions > 0),
                "Reader does not populate positions right away.")

    def test_rewind(self):
        trj = self.universe.trajectory
        trj.next()
        trj.next()  # for readers that do not support indexing
        assert_equal(trj.ts.frame, 2,
                     "failed to forward to frame 2 (frameindex 2)")
        trj.rewind()
        assert_equal(trj.ts.frame, 0, "failed to rewind to first frame")
        assert_(np.any(self.universe.atoms.positions > 0),
                "Reader does not populate positions after rewinding.")

    def test_full_slice(self):
        trj_iter = self.universe.trajectory[:]
        frames = [ts.frame for ts in trj_iter]
        assert_equal(frames, np.arange(self.universe.trajectory.n_frames))

    def test_random_access(self):
        u = self.universe

        pos1 = u.atoms[0].position
        u.trajectory.next()
        u.trajectory.next()
        pos3 = u.atoms[0].position

        u.trajectory[0]

        assert_equal(u.atoms[0].position, pos1)

        u.trajectory[2]

        assert_equal(u.atoms[0].position, pos3)


class TestTRJReader(_TRJReaderTest, RefACHE):
    def setUp(self):
        self.universe = mda.Universe(PRM, TRJ)
        self.prec = 3

    def test_read_frame_reopens(self):
        # should automatically reopen
        u = self.universe

        u.trajectory.close()

        u.trajectory[2]

        assert_(u.trajectory.ts.frame == 2)


class TestBzippedTRJReader(TestTRJReader):
    def setUp(self):
        self.universe = mda.Universe(PRM, TRJ_bz2)
        self.prec = 3


class TestBzippedTRJReaderPBC(_TRJReaderTest, RefCappedAla):
    def setUp(self):
        self.universe = mda.Universe(PRMpbc, TRJpbc_bz2)
        self.prec = 3



class TestTRJTimestep(BaseTimestepTest):
    Timestep = mda.coordinates.TRJ.Timestep
    name = "TRJ"
    has_box = True
    set_box = True
    unitcell = np.array([10., 11., 12., 90., 90., 90.])
    uni_args = (PRM, TRJ)

def test_trj_no_natoms():
    assert_raises(ValueError, mda.coordinates.TRJ.TRJReader, 'somefile.txt')

