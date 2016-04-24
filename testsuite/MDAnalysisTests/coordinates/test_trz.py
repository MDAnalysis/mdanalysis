import MDAnalysis as mda
import os
from six.moves import zip

from numpy.testing import (assert_equal, assert_array_almost_equal,
                           assert_almost_equal, assert_raises)
import numpy as np

from unittest import TestCase

from MDAnalysis import NoDataError

from MDAnalysisTests.coordinates.reference import RefTRZ
from MDAnalysisTests.coordinates.base import BaseTimestepTest
from MDAnalysisTests.datafiles import (TRZ_psf, TRZ, two_water_gro)
from MDAnalysisTests import tempdir


class TestTRZReader(TestCase, RefTRZ):
    def setUp(self):
        self.universe = mda.Universe(TRZ_psf, TRZ)
        self.trz = self.universe.trajectory
        self.ts = self.universe.trajectory.ts
        self.prec = 3

    def tearDown(self):
        del self.universe
        del self.trz
        del self.ts

    def test_load_trz(self):
        U = self.universe
        assert_equal(len(U.atoms), self.ref_n_atoms,
                     "load Universe from PSF and TRZ")

    def test_next_trz(self):
        assert_equal(self.ts.frame, 0, "starts at first frame")
        self.trz.next()
        assert_equal(self.ts.frame, 1, "next returns frame index 1")

    def test_rewind_trz(self):
        # move to different frame and rewind to get first frame back
        self.trz[2]
        self.trz.rewind()
        assert_equal(self.ts.frame, 0, "rewinding to frame 1")

    def test_n_frames(self):
        assert_equal(self.universe.trajectory.n_frames, self.ref_n_frames,
                     "wrong number of frames in trz")

    def test_seeking(self):
        self.universe.trajectory[3]
        assert_equal(self.ts.frame, 3, "loading frame 3")

        orig = self.universe.atoms[0:3].positions.copy()

        self.universe.trajectory[4]
        assert_equal(self.ts.frame, 4, "loading frame 4")
        self.universe.trajectory[3]

        assert_almost_equal(self.universe.atoms[0:3].positions, orig,
                            self.prec)

        self.universe.trajectory[0]
        assert_equal(self.ts.frame, 0, "loading frame 0")
        self.universe.trajectory[3]

        assert_almost_equal(self.universe.atoms[0:3].positions, orig,
                            self.prec)

    def test_volume(self):
        # Lower precision here because errors seem to accumulate and
        # throw this off (is rounded value**3)
        assert_almost_equal(self.ts.volume, self.ref_volume, 1,
                            "wrong volume for trz")

    def test_unitcell(self):
        assert_almost_equal(self.ts.dimensions, self.ref_dimensions, self.prec,
                            "wrong dimensions for trz")

    def test_coordinates(self):
        fortytwo = self.universe.atoms[41]  # 41 because is 0 based
        assert_almost_equal(fortytwo.pos, self.ref_coordinates, self.prec,
                            "wrong coordinates in trz")

    def test_velocities(self):
        fortytwo = self.universe.select_atoms('bynum 42')
        assert_almost_equal(fortytwo.velocities, self.ref_velocities,
                            self.prec, "wrong velocities in trz")

    def test_delta(self):
        assert_almost_equal(self.trz.delta, self.ref_delta, self.prec,
                            "wrong time delta in trz")

    def test_time(self):
        assert_almost_equal(self.trz.time, self.ref_time, self.prec,
                            "wrong time value in trz")

    def test_title(self):
        assert_equal(self.ref_title, self.trz.title, "wrong title in trz")

    def test_get_writer(self):
        with tempdir.in_tempdir():
            self.outfile = 'test-trz-writer.trz'
            W = self.trz.Writer(self.outfile)
            assert_equal(isinstance(W, mda.coordinates.TRZ.TRZWriter), True)
            assert_equal(W.n_atoms, self.trz.n_atoms)
            try:
                os.unlink(self.outfile)
            except OSError:
                pass

    def test_get_writer_2(self):
        with tempdir.in_tempdir():
            self.outfile = 'test-trz-writer-1.trz'
            W = self.trz.Writer(self.outfile, n_atoms=100)
            assert_equal(isinstance(W, mda.coordinates.TRZ.TRZWriter), True)
            assert_equal(W.n_atoms, 100)
            try:
                os.unlink(self.outfile)
            except OSError:
                pass


class TestTRZWriter(TestCase, RefTRZ):
    def setUp(self):
        self.universe = mda.Universe(TRZ_psf, TRZ)
        self.prec = 3
        self.tmpdir = tempdir.TempDir()
        self.outfile = self.tmpdir.name + '/test-trz-writer.trz'
        self.outfile_long = self.tmpdir.name + '/test-trz-writer-long.trz'
        self.Writer = mda.coordinates.TRZ.TRZWriter
        self.title_to_write = 'Test title TRZ'

    def tearDown(self):
        del self.universe
        del self.prec
        try:
            os.unlink(self.outfile)
            os.unlink(self.outfile_long)
        except OSError:
            pass
        del self.Writer
        del self.tmpdir

    def test_write_trajectory(self):
        t = self.universe.trajectory
        W = self.Writer(self.outfile, t.n_atoms, title=self.title_to_write)
        self._copy_traj(W)

    def _copy_traj(self, writer):
        for ts in self.universe.trajectory:
            writer.write_next_timestep(ts)
        writer.close()

        uw = mda.Universe(TRZ_psf, self.outfile)

        assert_equal(uw.trajectory.title, self.title_to_write,
                     "Title mismatch between original and written files.")

        for orig_ts, written_ts in zip(self.universe.trajectory,
                                       uw.trajectory):
            assert_array_almost_equal(orig_ts._pos, written_ts._pos, self.prec,
                                      err_msg="Coordinate mismatch between "
                                      "orig and written at frame %d" %
                                      orig_ts.frame)
            assert_array_almost_equal(orig_ts._velocities,
                                      written_ts._velocities, self.prec,
                                      err_msg="Coordinate mismatch between "
                                      "orig and written at frame %d" %
                                      orig_ts.frame)
            assert_array_almost_equal(orig_ts._unitcell, written_ts._unitcell,
                                      self.prec, err_msg="Unitcell mismatch "
                                      "between orig and written at frame %d" %
                                      orig_ts.frame)
            for att in orig_ts.data:
                assert_array_almost_equal(orig_ts.data[att],
                                          written_ts.data[att], self.prec,
                                          err_msg="TS equal failed for {0!s}".format(att))

    def test_long_title(self):
        title = '*' * 81
        assert_raises(ValueError,
                      self.Writer, self.outfile, self.ref_n_atoms, title=title)


class TestTRZWriter2(object):
    def setUp(self):
        self.u = mda.Universe(two_water_gro)

    def tearDown(self):
        del self.u

    def test_writer_trz_from_other(self):
        with tempdir.in_tempdir():
            outfile = 'trz-writer-2.trz'
            W = mda.coordinates.TRZ.TRZWriter(outfile,
                                              n_atoms=len(self.u.atoms))

            W.write(self.u.trajectory.ts)
            W.close()

            u2 = mda.Universe(two_water_gro, outfile)

            assert_array_almost_equal(self.u.atoms.positions,
                                      u2.atoms.positions, 3)


class TestWrite_Partial_Timestep(TestCase):
    """Test writing a partial timestep made by passing only an atomgroup to
    Writer. (Issue 163)

    The contents of the AtomGroup.ts are checked in test_atomgroup, this test
    just checks that Writer is receiving this information properly.

    """

    def setUp(self):
        self.universe = mda.Universe(TRZ_psf, TRZ)
        self.ag = self.universe.select_atoms('name N')
        self.prec = 3
        self.tmpdir = tempdir.TempDir()
        self.outfile = self.tmpdir.name + '/partial-write-test.pdb'
        self.Writer = mda.Writer(self.outfile, n_atoms=len(self.ag))

    def tearDown(self):
        del self.universe
        del self.ag
        del self.prec
        try:
            os.unlink(self.outfile)
        except OSError:
            pass
        del self.Writer
        del self.tmpdir

    def test_write_trajectory(self):
        self.Writer.write(self.ag)
        self.Writer.close()

        u_ag = mda.Universe(self.outfile)

        assert_array_almost_equal(self.ag.positions,
                                  u_ag.atoms.positions,
                                  self.prec,
                                  err_msg="Writing AtomGroup timestep failed.")


class TestTRZTimestep(BaseTimestepTest):
    Timestep = mda.coordinates.TRZ.Timestep
    name = "TRZ"
    has_box = True
    set_box = True
    unitcell = np.array([10., 0., 0.,
                         0., 11., 0.,
                         0., 0., 12.])
    uni_args = (TRZ_psf, TRZ)
