import errno
import MDAnalysis as mda
import numpy as np
import os
from six.moves import zip

from nose.plugins.attrib import attr
from numpy.testing import (assert_equal, assert_array_almost_equal, dec,
                           assert_almost_equal, assert_raises,
                           assert_array_equal)
import tempdir
from unittest import TestCase

from MDAnalysisTests.datafiles import (PDB_sub_dry, PDB_sub_sol, TRR_sub_sol,
                                       TRR, XTC, GRO, PDB, CRD, PRMncdf, NCDF)


class TestTRRReader_Sub(TestCase):
    def setUp(self):
        """
        grab values from selected atoms from full solvated traj,
        later compare to using 'sub'
        """
        usol = mda.Universe(PDB_sub_sol, TRR_sub_sol)
        atoms = usol.select_atoms("not resname SOL")
        self.pos = atoms.positions
        self.vel = atoms.velocities
        self.force = atoms.forces
        self.sub = atoms.indices
        # universe from un-solvated protein
        self.udry = mda.Universe(PDB_sub_dry)

    def test_load_new_raises_ValueError(self):
        # should fail if we load universe with a trajectory with different
        # number of atoms when NOT using sub, same as before.
        def load_new_without_sub():
            self.udry.load_new(TRR_sub_sol)

        assert_raises(ValueError, load_new_without_sub)

    def test_sub_coordinates(self):
        """
        load solvated trajectory into universe with unsolvated protein.
        """
        self.udry.load_new(TRR_sub_sol, sub=self.sub)
        assert_array_almost_equal(self.pos,
                                  self.udry.atoms.positions,
                                  err_msg="positions differ")
        assert_array_almost_equal(self.vel,
                                  self.udry.atoms.velocities,
                                  err_msg="positions differ")
        assert_array_almost_equal(self.force,
                                  self.udry.atoms.forces,
                                  err_msg="positions differ")


class _GromacsReader(TestCase):
    # This base class assumes same lengths and dt for XTC and TRR test cases!
    filename = None
    ref_unitcell = np.array([80.017, 80.017, 80.017, 60., 60., 90.],
                            dtype=np.float32)
    # computed with Gromacs: 362.26999999999998 nm**3 * 1000 A**3/nm**3
    ref_volume = 362270.0

    def setUp(self):
        # loading from GRO is 4x faster than the PDB reader
        self.universe = mda.Universe(GRO, self.filename, convert_units=True)
        self.trajectory = self.universe.trajectory
        self.prec = 3
        self.ts = self.universe.coord
        # dummy output file
        ext = os.path.splitext(self.filename)[1]
        self.tmpdir = tempdir.TempDir()
        self.outfile = self.tmpdir.name + '/xdr-reader-test' + ext

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except:
            pass
        del self.universe
        del self.tmpdir

    @dec.slow
    def test_flag_convert_lengths(self):
        assert_equal(mda.core.flags['convert_lengths'], True,
                     "MDAnalysis.core.flags['convert_lengths'] should be "
                     "True by default")

    @dec.slow
    def test_rewind_xdrtrj(self):
        self.trajectory.rewind()
        assert_equal(self.ts.frame, 0, "rewinding to frame 1")

    @dec.slow
    def test_next_xdrtrj(self):
        self.trajectory.rewind()
        self.trajectory.next()
        assert_equal(self.ts.frame, 1, "loading frame 1")

    @dec.slow
    def test_jump_xdrtrj(self):
        self.trajectory[4]  # index is 0-based and frames are 0-based
        assert_equal(self.ts.frame, 4, "jumping to frame 4")

    @dec.slow
    def test_jump_lastframe_xdrtrj(self):
        self.trajectory[-1]
        assert_equal(self.ts.frame, 9,
                     "indexing last frame with trajectory[-1]")

    @dec.slow
    def test_slice_xdrtrj(self):
        frames = [ts.frame for ts in self.trajectory[2:9:3]]
        assert_equal(frames, [2, 5, 8], "slicing xdrtrj [2:9:3]")

    @dec.slow
    def test_reverse_xdrtrj(self):
        frames = [ts.frame for ts in self.trajectory[::-1]]
        assert_equal(frames, range(9, -1, -1), "slicing xdrtrj [::-1]")

    @dec.slow
    def test_coordinates(self):
        ca_nm = np.array([[6.043369675, 7.385184479, 1.381425762]],
                         dtype=np.float32)
        # coordinates in the base unit (needed for True)
        ca_Angstrom = ca_nm * 10.0
        U = self.universe
        T = U.trajectory
        T.rewind()
        T.next()
        T.next()
        assert_equal(self.ts.frame, 2, "failed to step to frame 3")
        ca = U.select_atoms('name CA and resid 122')
        # low precision match (2 decimals in A, 3 in nm) because the above are
        # the trr coords
        assert_array_almost_equal(ca.coordinates(), ca_Angstrom, 2,
                                  err_msg="coords of Ca of resid 122 do not "
                                  "match for frame 3")

    @dec.slow
    @attr('issue')
    def test_unitcell(self):
        """Test that xtc/trr unitcell is read correctly (Issue 34)"""
        self.universe.trajectory.rewind()
        uc = self.ts.dimensions
        assert_array_almost_equal(
            uc,
            self.ref_unitcell,
            self.prec,
            err_msg="unit cell dimensions (rhombic dodecahedron)")

    @dec.slow
    def test_volume(self):
        # need to reduce precision for test (nm**3 <--> A**3)
        self.universe.trajectory.rewind()
        vol = self.ts.volume
        assert_array_almost_equal(
            vol,
            self.ref_volume,
            0,
            err_msg="unit cell volume (rhombic dodecahedron)")

    @dec.slow
    def test_dt(self):
        assert_almost_equal(self.universe.trajectory.dt,
                            100.0,
                            4,
                            err_msg="wrong timestep dt")

    @dec.slow
    def test_totaltime(self):
        # test_totaltime(): need to reduce precision because dt is only precise
        # to ~4 decimals and accumulating the inaccuracy leads to even lower
        # precision in the totaltime (consequence of fixing Issue 64)
        assert_almost_equal(self.universe.trajectory.totaltime,
                            1000.0,
                            3,
                            err_msg="wrong total length of trajectory")

    @dec.slow
    def test_frame(self):
        self.trajectory[4]  # index is 0-based and frames are 0-based
        assert_equal(self.universe.trajectory.frame, 4, "wrong frame number")

    @dec.slow
    def test_time(self):
        self.trajectory[4]
        assert_almost_equal(self.universe.trajectory.time,
                            400.0,
                            3,
                            err_msg="wrong time of frame")

    @dec.slow
    def test_get_Writer(self):
        W = self.universe.trajectory.Writer(self.outfile)
        assert_equal(self.universe.trajectory.format, W.format)
        assert_equal(self.universe.atoms.n_atoms, W.n_atoms)
        W.close()

    @dec.slow
    def test_Writer(self):
        W = self.universe.trajectory.Writer(self.outfile)
        W.write(self.universe.atoms)
        self.universe.trajectory.next()
        W.write(self.universe.atoms)
        W.close()
        self.universe.trajectory.rewind()
        u = mda.Universe(GRO, self.outfile)
        assert_equal(u.trajectory.n_frames, 2)
        # prec = 6: TRR test fails; here I am generous and take self.prec =
        # 3...
        assert_almost_equal(u.atoms.coordinates(),
                            self.universe.atoms.coordinates(), self.prec)

    @dec.slow
    def test_EOFraisesIOErrorEIO(self):
        def go_beyond_EOF():
            self.universe.trajectory[-1]
            self.universe.trajectory.next()

        assert_raises(IOError, go_beyond_EOF)
        try:
            go_beyond_EOF()
        except IOError as err:
            assert_equal(err.errno, errno.EIO,
                         "IOError produces wrong error code")


class TestXTCReader(_GromacsReader):
    filename = XTC


class TestXTCReaderClass(TestCase):
    def test_with_statement(self):
        from MDAnalysis.coordinates.XTC import XTCReader

        try:
            with XTCReader(XTC) as trj:
                N = trj.n_frames
                frames = [ts.frame for ts in trj]
        except:
            raise AssertionError("with_statement not working for XTCReader")
        assert_equal(
            N,
            10,
            err_msg="with_statement: XTCReader reads wrong number of frames")
        assert_array_equal(
            frames,
            np.arange(0, N),
            err_msg="with_statement: XTCReader does not read all frames")


class TestTRRReader(_GromacsReader):
    filename = TRR

    @dec.slow
    def test_velocities(self):
        # frame 0, v in nm/ps
        # from gmxdump -f MDAnalysisTests/data/adk_oplsaa.trr
        #      v[47675]={-7.86469e-01,  1.57479e+00,  2.79722e-01}
        #      v[47676]={ 2.70593e-08,  1.08052e-06,  6.97028e-07}
        v_native = np.array(
            [
                [-7.86469e-01, 1.57479e+00, 2.79722e-01
                 ], [2.70593e-08, 1.08052e-06, 6.97028e-07]
            ],
            dtype=np.float32)

        # velocities in the MDA base unit A/ps (needed for True)
        v_base = v_native * 10.0
        self.universe.trajectory.rewind()
        assert_equal(self.ts.frame, 0, "failed to read frame 1")

        assert_array_almost_equal(
            self.universe.trajectory.ts._velocities[[47675, 47676]], v_base,
            self.prec, err_msg="ts._velocities for indices 47675,47676 do not "
            "match known values")

        assert_array_almost_equal(
            self.universe.atoms.velocities[[47675, 47676]], v_base,
            self.prec, err_msg="velocities for indices 47675,47676 do not "
            "match known values")

        for index, v_known in zip([47675, 47676], v_base):
            assert_array_almost_equal(
                self.universe.atoms[index].velocity,
                v_known,
                self.prec,
                err_msg="atom[%d].velocity does not match known values" %
                index)


class _XDRNoConversion(TestCase):
    filename = None

    def setUp(self):
        self.universe = mda.Universe(PDB, self.filename, convert_units=False)
        self.ts = self.universe.trajectory.ts

    def tearDown(self):
        del self.universe
        del self.ts

    @dec.slow
    def test_coordinates(self):
        # note: these are the native coordinates in nm
        ca_nm = np.array([[6.043369675, 7.385184479, 1.381425762]],
                         dtype=np.float32)
        U = self.universe
        T = U.trajectory
        T.rewind()
        T.next()
        T.next()
        assert_equal(self.ts.frame, 2, "failed to step to frame 3")
        ca = U.select_atoms('name CA and resid 122')
        # low precision match because we also look at the trr: only 3 decimals
        # in nm in xtc!
        assert_array_almost_equal(ca.coordinates(), ca_nm, 3,
                                  err_msg="native coords of Ca of resid 122 "
                                  "do not match for frame 3 with "
                                  "convert_units=False")


class TestXTCNoConversion(_XDRNoConversion):
    filename = XTC


class TestTRRNoConversion(_XDRNoConversion):
    filename = TRR


class _GromacsWriter(TestCase):
    infilename = None  # XTC or TRR
    Writers = {
        '.trr': mda.coordinates.TRR.TRRWriter,
        '.xtc': mda.coordinates.XTC.XTCWriter,
    }

    def setUp(self):
        self.universe = mda.Universe(GRO, self.infilename)
        ext = os.path.splitext(self.infilename)[1]
        self.tmpdir = tempdir.TempDir()
        self.outfile = self.tmpdir.name + '/xdr-writer-test' + ext
        self.Writer = self.Writers[ext]

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except:
            pass
        del self.universe
        del self.Writer
        del self.tmpdir

    @dec.slow
    @attr('issue')
    def test_write_trajectory(self):
        """Test writing Gromacs trajectories (Issue 38)"""
        t = self.universe.trajectory
        W = self.Writer(self.outfile, t.n_atoms, dt=t.dt)
        for ts in self.universe.trajectory:
            W.write_next_timestep(ts)
        W.close()

        uw = mda.Universe(GRO, self.outfile)

        # check that the coordinates are identical for each time step
        for orig_ts, written_ts in zip(self.universe.trajectory,
                                       uw.trajectory):
            assert_array_almost_equal(written_ts._pos, orig_ts._pos, 3,
                                      err_msg="coordinate mismatch between "
                                      "original and written trajectory at "
                                      "frame %d (orig) vs %d (written)" % (
                                           orig_ts.frame, written_ts.frame))

    @dec.slow
    def test_timestep_not_modified_by_writer(self):
        trj = self.universe.trajectory
        ts = trj.ts

        trj[-1]  # last timestep (so that time != 0)
        x = ts._pos.copy()
        time = ts.time

        W = self.Writer(self.outfile, trj.n_atoms, dt=trj.dt)
        # last timestep (so that time != 0) (say it again, just in case...)
        trj[-1]
        W.write_next_timestep(ts)
        W.close()

        assert_equal(ts._pos,
                     x,
                     err_msg="Positions in Timestep were modified by writer.")
        assert_equal(ts.time,
                     time,
                     err_msg="Time in Timestep was modified by writer.")


class TestXTCWriter(_GromacsWriter):
    infilename = XTC


class TestTRRWriter(_GromacsWriter):
    infilename = TRR

    def test_velocities(self):
        t = self.universe.trajectory
        W = self.Writer(self.outfile, t.n_atoms, dt=t.dt)
        for ts in self.universe.trajectory:
            W.write_next_timestep(ts)
        W.close()

        uw = mda.Universe(GRO, self.outfile)

        # check that the velocities are identical for each time step
        for orig_ts, written_ts in zip(self.universe.trajectory,
                                       uw.trajectory):
            assert_array_almost_equal(written_ts._velocities,
                                      orig_ts._velocities, 3,
                                      err_msg="velocities mismatch between "
                                      "original and written trajectory at "
                                      "frame %d (orig) vs %d (written)" % (
                                          orig_ts.frame, written_ts.frame))

    def test_gaps(self):
        """Tests the writing and reading back of TRRs with gaps in any of
        the coordinates/velocities properties."""
        t = self.universe.trajectory
        W = self.Writer(self.outfile, t.n_atoms, dt=t.dt)
        for ts in self.universe.trajectory:
            # Inset some gaps in the properties: coords every 4 steps, vels
            # every 2.
            if not ts.frame % 4:
                ts.has_positions = False
            if not ts.frame % 2:
                ts.has_velocities = False
            W.write_next_timestep(ts)
        W.close()

        uw = mda.Universe(GRO, self.outfile)

        # check that the velocities are identical for each time step, except
        # for the gaps (that we must make sure to raise exceptions on).
        for orig_ts, written_ts in zip(self.universe.trajectory,
                                       uw.trajectory):
            if ts.frame % 4:
                assert_array_almost_equal(written_ts.positions,
                                          orig_ts.positions, 3,
                                          err_msg="coordinates mismatch "
                                          "between original and written "
                                          "trajectory at frame {} (orig) "
                                          "vs {} (written)".format(
                                              orig_ts.frame, written_ts.frame))
            else:
                assert_raises(mda.NoDataError, getattr, written_ts,
                              'positions')

            if ts.frame % 2:
                assert_array_almost_equal(written_ts.velocities,
                                          orig_ts.velocities, 3,
                                          err_msg="velocities mismatch "
                                          "between original and written "
                                          "trajectory at frame {} (orig) "
                                          "vs {} (written)".format(
                                              orig_ts.frame, written_ts.frame))
            else:
                assert_raises(mda.NoDataError, getattr, written_ts,
                              'velocities')


class _GromacsWriterIssue101(TestCase):
    Writers = {
        '.trr': mda.coordinates.TRR.TRRWriter,
        '.xtc': mda.coordinates.XTC.XTCWriter,
    }
    ext = None  # set to '.xtc' or '.trr'
    prec = 3

    def setUp(self):
        self.tmpdir = tempdir.TempDir()
        self.outfile = self.tmpdir.name + '/xdr-writer-issue101' + self.ext
        self.Writer = self.Writers[self.ext]

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except:
            pass
        del self.Writer
        del self.tmpdir

    @dec.slow
    @attr('issue')
    def test_single_frame_GRO(self):
        self._single_frame(GRO)

    @dec.slow
    @attr('issue')
    def test_single_frame_PDB(self):
        self._single_frame(PDB)

    @attr('issue')
    def test_single_frame_CRD(self):
        self._single_frame(CRD)

    def _single_frame(self, filename):
        u = mda.Universe(filename)
        with self.Writer(self.outfile, u.atoms.n_atoms) as W:
            W.write(u.atoms)
        w = mda.Universe(filename, self.outfile)
        assert_equal(w.trajectory.n_frames, 1,
                     "single frame trajectory has wrong number of frames")
        assert_almost_equal(
            w.atoms.coordinates(),
            u.atoms.coordinates(),
            self.prec,
            err_msg="coordinates do not match for %r" % filename)


class TestXTCWriterSingleFrame(_GromacsWriterIssue101):
    ext = ".xtc"
    prec = 2


class TestTRRWriterSingleFrame(_GromacsWriterIssue101):
    ext = ".trr"


class _GromacsWriterIssue117(TestCase):
    """Issue 117: Cannot write XTC or TRR from AMBER NCDF"""
    ext = None
    prec = 5

    def setUp(self):
        self.universe = mda.Universe(PRMncdf, NCDF)
        self.tmpdir = tempdir.TempDir()
        self.outfile = self.tmpdir.name + '/xdr-writer-issue117' + self.ext
        self.Writer = mda.Writer(self.outfile,
                                 n_atoms=self.universe.atoms.n_atoms)

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except:
            pass
        del self.universe
        del self.Writer

    @attr('issue')
    def test_write_trajectory(self):
        """Test writing Gromacs trajectories from AMBER NCDF (Issue 117)"""
        self.universe.trajectory
        for ts in self.universe.trajectory:
            self.Writer.write_next_timestep(ts)
        self.Writer.close()

        uw = mda.Universe(PRMncdf, self.outfile)

        # check that the coordinates are identical for each time step
        for orig_ts, written_ts in zip(self.universe.trajectory,
                                       uw.trajectory):
            assert_array_almost_equal(written_ts._pos, orig_ts._pos,
                                      self.prec, err_msg="coordinate mismatch "
                                      "between original and written "
                                      "trajectory at frame %d (orig) vs %d "
                                      "(written)" % (
                                           orig_ts.frame, written_ts.frame))


class TestXTCWriterIssue117(_GromacsWriterIssue117):
    ext = ".xtc"
    prec = 2


class TestTRRWriterIssue117(_GromacsWriterIssue117):
    ext = ".trr"


@attr('issue')
def test_triclinic_box():
    """Test coordinates.core.triclinic_box() (Issue 61)"""
    unitcell = np.array([80.017, 55, 100.11, 60.00, 30.50, 90.00])
    box = mda.coordinates.core.triclinic_vectors(unitcell)
    new_unitcell = mda.coordinates.core.triclinic_box(box[0], box[1],
                                                      box[2])
    assert_array_almost_equal(
        new_unitcell,
        unitcell,
        3,
        err_msg="unitcell round-trip connversion failed (Issue 61)")
