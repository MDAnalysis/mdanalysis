import MDAnalysis as mda
import numpy as np
import os
from six.moves import zip

from numpy.testing import (assert_equal, assert_array_almost_equal,
                           assert_almost_equal)
import tempdir
from unittest import TestCase

from MDAnalysisTests.coordinates.reference import Ref2r9r
from MDAnalysisTests.datafiles import (XYZ, XYZ_psf, XYZ_bz2, XYZ_mini)


class TestXYZReader(TestCase, Ref2r9r):
    def setUp(self):
        self.universe = mda.Universe(XYZ_psf, XYZ)
        self.prec = 3  # 4 decimals in xyz file

    def tearDown(self):
        del self.universe

    def test_load_xyz(self):
        U = self.universe
        assert_equal(len(U.atoms), self.ref_n_atoms,
                     "load Universe from PSF and XYZ")

    def test_n_atoms(self):
        assert_equal(self.universe.trajectory.n_atoms, self.ref_n_atoms,
                     "wrong number of atoms")

    def test_n_frames(self):
        assert_equal(self.universe.trajectory.n_frames, self.ref_n_frames,
                     "wrong number of frames in xyz")

    def test_sum_centres_of_geometry(self):
        centreOfGeometry = 0

        for i in self.universe.trajectory:
            sel = self.universe.select_atoms("all")
            centreOfGeometry += sum(sel.center_of_geometry())

        assert_almost_equal(centreOfGeometry, self.ref_sum_centre_of_geometry,
                            self.prec, err_msg="sum of centers of geometry "
                            "over the trajectory do not match")

    def test_full_slice(self):
        trj_iter = self.universe.trajectory[:]
        frames = [ts.frame for ts in trj_iter]
        assert_equal(frames, np.arange(self.universe.trajectory.n_frames))

    def test_slice(self):
        frames = [ts.frame for ts in self.universe.trajectory[::2]]
        assert_equal(frames, np.arange(len(self.universe.trajectory))[::2])


class TestXYZReaderAsTopology(object):
    """Test that an XYZ file can act as its own topology"""

    def setUp(self):
        self.universe = mda.Universe(XYZ_mini)

    def tearDown(self):
        del self.universe

    def test_coords(self):
        ref = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [2.0, 2.0, 2.0]],
                       dtype=np.float32)
        assert_array_almost_equal(self.universe.atoms.positions, ref)

    def test_rewind(self):
        self.universe.trajectory.rewind()
        assert_equal(self.universe.trajectory.ts.frame, 0,
                     "rewinding to frame 0")

    def test_dt(self):
        assert_almost_equal(self.universe.trajectory.dt,
                            1.0,
                            4,
                            err_msg="wrong timestep dt")


class TestCompressedXYZReader(TestCase, Ref2r9r):
    def setUp(self):
        self.universe = mda.Universe(XYZ_psf, XYZ_bz2)
        self.prec = 3  # 4 decimals in xyz file

    def tearDown(self):
        del self.universe

    def test_load_xyz(self):
        U = self.universe
        assert_equal(len(U.atoms), self.ref_n_atoms,
                     "load Universe from PSF and XYZ")

    def test_n_atoms(self):
        assert_equal(self.universe.trajectory.n_atoms, self.ref_n_atoms,
                     "wrong number of atoms")

    def test_n_frames(self):
        assert_equal(self.universe.trajectory.n_frames, self.ref_n_frames,
                     "wrong number of frames in xyz")

    def test_sum_centres_of_geometry(self):
        centreOfGeometry = 0

        for i in self.universe.trajectory:
            sel = self.universe.select_atoms("all")
            centreOfGeometry += sum(sel.center_of_geometry())

        assert_almost_equal(centreOfGeometry, self.ref_sum_centre_of_geometry,
                            self.prec, err_msg="sum of centers of geometry "
                            "over the trajectory do not match")

    def test_full_slice(self):
        trj_iter = self.universe.trajectory[:]
        frames = [ts.frame for ts in trj_iter]
        assert_equal(frames, np.arange(self.universe.trajectory.n_frames))

    def test_slice(self):
        frames = [ts.frame for ts in self.universe.trajectory[::2]]
        assert_equal(frames, np.arange(len(self.universe.trajectory))[::2])

    def test_rewind(self):
        self.universe.trajectory.rewind()
        assert_equal(self.universe.trajectory.ts.frame, 0,
                     "rewinding to frame 0")

    def test_next(self):
        self.universe.trajectory.rewind()
        self.universe.trajectory.next()
        assert_equal(self.universe.trajectory.ts.frame, 1, "loading frame 1")

    def test_dt(self):
        assert_almost_equal(self.universe.trajectory.dt,
                            1.0,
                            4,
                            err_msg="wrong timestep dt")


class TestXYZWriter(TestCase, Ref2r9r):
    def setUp(self):
        self.universe = mda.Universe(XYZ_psf, XYZ_bz2)
        self.prec = 3  # 4 decimals in xyz file
        ext = ".xyz"
        self.tmpdir = tempdir.TempDir()
        self.outfile = self.tmpdir.name + '/xyz-writer' + ext
        self.Writer = mda.coordinates.XYZ.XYZWriter

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except OSError:
            pass
        del self.universe
        del self.tmpdir

    def test_write_trajectory_timestep(self):
        W = self.Writer(self.outfile)
        self._copy_traj(W)

    def test_write_trajectory_atomgroup(self):
        W = self.Writer(self.outfile)
        for ts in self.universe.trajectory:
            W.write(self.universe.atoms)
        W.close()
        self._check_copy()

    def test_ReaderWriter(self):
        t = self.universe.trajectory
        W = t.Writer(self.outfile)
        self._copy_traj(W)

    def _copy_traj(self, writer):
        for ts in self.universe.trajectory:
            writer.write_next_timestep(ts)
        writer.close()
        self._check_copy()

    def _check_copy(self):
        uw = mda.Universe(XYZ_psf, self.outfile)
        assert_equal(self.universe.trajectory.n_frames, uw.trajectory.n_frames)
        # check that the trajectories are identical for each time step
        for orig_ts, written_ts in zip(self.universe.trajectory,
                                       uw.trajectory):
            assert_array_almost_equal(written_ts._pos, orig_ts._pos, self.prec,
                                      err_msg="coordinate mismatch between "
                                      "original and written trajectory at "
                                      "frame %d (orig) vs %d (written)" % (
                                           orig_ts.frame, written_ts.frame))
