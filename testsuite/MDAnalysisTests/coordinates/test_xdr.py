from six.moves import zip, range

import errno
import MDAnalysis as mda
from MDAnalysis.coordinates.base import Timestep
import numpy as np
import os
import shutil
import warnings

from nose.plugins.attrib import attr
from numpy.testing import (assert_equal, assert_array_almost_equal, dec,
                           assert_almost_equal, assert_raises,
                           assert_array_equal)
from unittest import TestCase


from MDAnalysisTests import module_not_found

from MDAnalysisTests.datafiles import (PDB_sub_dry, PDB_sub_sol, TRR_sub_sol,
                                       TRR, XTC, GRO, PDB, CRD, PRMncdf, NCDF,
                                       XTC_sub_sol)

from MDAnalysisTests.datafiles import (COORDINATES_XTC, COORDINATES_TOPOLOGY,
                                       COORDINATES_TRR)
from MDAnalysisTests.coordinates.base import (BaseReaderTest, BaseReference,
                                              BaseWriterTest,
                                              assert_timestep_almost_equal)
from MDAnalysisTests import tempdir

import MDAnalysis.core.AtomGroup
from MDAnalysis.coordinates import XDR

# I want to catch all warnings in the tests. If this is not set at the start it
# could cause test that check for warnings to fail.
warnings.simplefilter('always')


class _XDRReader_Sub(TestCase):

    def setUp(self):
        """
        grab values from selected atoms from full solvated traj,
        later compare to using 'sub'
        """
        usol = mda.Universe(PDB_sub_sol, self.XDR_SUB_SOL)
        atoms = usol.select_atoms("not resname SOL")
        self.ts = atoms.ts
        self.sub = atoms.indices
        # universe from un-solvated protein
        self.udry = mda.Universe(PDB_sub_dry)

    def test_load_new_raises_ValueError(self):
        # should fail if we load universe with a trajectory with different
        # number of atoms when NOT using sub, same as before.
        def load_new_without_sub():
            self.udry.load_new(self.XDR_SUB_SOL)

        assert_raises(ValueError, load_new_without_sub)

    def test_sub_coordinates(self):
        """
        load solvated trajectory into universe with unsolvated protein.
        """
        self.udry.load_new(self.XDR_SUB_SOL, sub=self.sub)
        ts = self.udry.atoms.ts
        assert_timestep_almost_equal(ts, self.ts)


class TestTRRReader_Sub(_XDRReader_Sub):
    XDR_SUB_SOL = TRR_sub_sol


class TestXTCReader_Sub(_XDRReader_Sub):
    XDR_SUB_SOL = XTC_sub_sol


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
        assert_equal(frames, list(range(9, -1, -1)), "slicing xdrtrj [::-1]")

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
        assert_array_almost_equal(ca.positions, ca_Angstrom, 2,
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
        assert_almost_equal(u.atoms.positions,
                            self.universe.atoms.positions, self.prec)

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
                err_msg="atom[{0:d}].velocity does not match known values".format(
                index))


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
        assert_array_almost_equal(ca.positions, ca_nm, 3,
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
        except OSError:
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
            if ts.frame % 4 == 0:
                ts.has_positions = False
            if ts.frame % 2 == 0:
                ts.has_velocities = False
            W.write_next_timestep(ts)
        W.close()

        uw = mda.Universe(GRO, self.outfile)
        # check that the velocities are identical for each time step, except
        # for the gaps (that we must make sure to raise exceptions on).
        for orig_ts, written_ts in zip(self.universe.trajectory,
                                       uw.trajectory):
            if ts.frame % 4 != 0:
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

            if ts.frame % 2 != 0:
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
            w.atoms.positions,
            u.atoms.positions,
            self.prec,
            err_msg="coordinates do not match for {0!r}".format(filename))


class TestXTCWriterSingleFrame(_GromacsWriterIssue101):
    ext = ".xtc"
    prec = 2


class TestTRRWriterSingleFrame(_GromacsWriterIssue101):
    ext = ".trr"


class _GromacsWriterIssue117(TestCase):
    """Issue 117: Cannot write XTC or TRR from AMBER NCDF"""
    ext = None
    prec = 5

    @dec.skipif(module_not_found("netCDF4"), "Test skipped because netCDF is not available.")
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


class XTCReference(BaseReference):
    def __init__(self):
        super(XTCReference, self).__init__()
        self.trajectory = COORDINATES_XTC
        self.topology = COORDINATES_TOPOLOGY
        self.reader = mda.coordinates.XTC.XTCReader
        self.writer = mda.coordinates.XTC.XTCWriter
        self.ext = 'xtc'
        self.prec = 3
        self.changing_dimensions = True


class TestXTCReader_2(BaseReaderTest):
    def __init__(self, reference=None):
        if reference is None:
            reference = XTCReference()
        super(TestXTCReader_2, self).__init__(reference)


class TestXTCWriter_2(BaseWriterTest):
    def __init__(self, reference=None):
        if reference is None:
            reference = XTCReference()
        super(TestXTCWriter_2, self).__init__(reference)

    def test_different_precision(self):
        out = self.tmp_file('precision-test')
        # store more then 9 atoms to enable compression
        n_atoms = 40
        with self.ref.writer(out, n_atoms, precision=5) as w:
            ts = Timestep(n_atoms=n_atoms)
            ts.positions = np.random.random(size=(n_atoms, 3))
            w.write(ts)
        xtc = mda.lib.formats.libmdaxdr.XTCFile(out)
        frame = xtc.read()
        assert_equal(len(xtc), 1)
        assert_equal(xtc.n_atoms, n_atoms)
        assert_equal(frame.prec, 10.0 ** 5)


class TRRReference(BaseReference):
    def __init__(self):
        super(TRRReference, self).__init__()
        self.trajectory = COORDINATES_TRR
        self.topology = COORDINATES_TOPOLOGY
        self.changing_dimensions = True
        self.reader = mda.coordinates.TRR.TRRReader
        self.writer = mda.coordinates.TRR.TRRWriter
        self.ext = 'xtc'
        self.prec = 3
        self.first_frame.velocities = self.first_frame.positions / 10
        self.first_frame.forces = self.first_frame.positions / 100

        self.second_frame.velocities = self.second_frame.positions / 10
        self.second_frame.forces = self.second_frame.positions / 100

        self.last_frame.velocities = self.last_frame.positions / 10
        self.last_frame.forces = self.last_frame.positions / 100

        self.jump_to_frame.velocities = self.jump_to_frame.positions / 10
        self.jump_to_frame.forces = self.jump_to_frame.positions / 100

    def iter_ts(self, i):
        ts = self.first_frame.copy()
        ts.positions = 2**i * self.first_frame.positions
        ts.velocities = ts.positions / 10
        ts.forces = ts.positions / 100
        ts.time = i
        ts.frame = i
        return ts


class TestTRRReader_2(BaseReaderTest):
    def __init__(self, reference=None):
        if reference is None:
            reference = TRRReference()
        super(TestTRRReader_2, self).__init__(reference)


class TestTRRWriter_2(BaseWriterTest):
    def __init__(self, reference=None):
        if reference is None:
            reference = TRRReference()
        super(TestTRRWriter_2, self).__init__(reference)


class _GromacsReader_offsets(TestCase):
    # This base class assumes same lengths and dt for XTC and TRR test cases!
    filename = None
    ref_unitcell = np.array([80.017, 80.017, 80.017, 60., 60., 90.],
                            dtype=np.float32)
    # computed with Gromacs: 362.26999999999998 nm**3 * 1000 A**3/nm**3
    ref_volume = 362270.0
    ref_offsets = None
    _reader = None

    def setUp(self):
        # since offsets are automatically generated in the same directory
        # as the trajectory, we do everything from a temporary directory
        self.tmpdir = tempdir.TempDir()
        shutil.copy(self.filename, self.tmpdir.name)

        self.traj = os.path.join(self.tmpdir.name,
                                 os.path.basename(self.filename))

        self.trajectory = self._reader(self.traj)
        self.prec = 3
        self.ts = self.trajectory.ts

    def tearDown(self):
        del self.tmpdir
        del self.trajectory

    @dec.slow
    def test_offsets(self):
        self.trajectory._read_offsets(store=True)
        assert_array_almost_equal(self.trajectory._xdr.offsets,
                                  self.ref_offsets,
                                  err_msg="wrong frame offsets")

        outfile_offsets = XDR.offsets_filename(self.traj)
        saved_offsets = XDR.read_numpy_offsets(outfile_offsets)

        assert_array_almost_equal(self.trajectory._xdr.offsets,
                                  saved_offsets['offsets'],
                                  err_msg="error saving frame offsets")
        assert_array_almost_equal(self.ref_offsets, saved_offsets['offsets'],
                                  err_msg="saved frame offsets don't match "
                                  "the known ones")

        self.trajectory._load_offsets()
        assert_array_almost_equal(self.trajectory._xdr.offsets,
                                  self.ref_offsets,
                                  err_msg="error loading frame offsets")
        assert_equal(saved_offsets['ctime'], os.path.getctime(self.traj))
        assert_equal(saved_offsets['size'], os.path.getsize(self.traj))

    def test_reload_offsets(self):
        self._reader(self.traj, refresh_offsets=True)

    @dec.slow
    def test_persistent_offsets_size_mismatch(self):
        # check that stored offsets are not loaded when trajectory
        # size differs from stored size
        fname = XDR.offsets_filename(self.traj)
        saved_offsets = XDR.read_numpy_offsets(fname)
        saved_offsets['size'] += 1
        with open(fname, 'wb') as f:
            np.savez(f, **saved_offsets)

        with warnings.catch_warnings(record=True) as warn:
            warnings.simplefilter('always')
            self._reader(self.traj)
        assert_equal(warn[0].message.args,
                     ('Reload offsets from trajectory\n ctime or size or n_atoms did not match', ))

    @dec.slow
    def test_persistent_offsets_ctime_mismatch(self):
        # check that stored offsets are not loaded when trajectory
        # ctime differs from stored ctime
        fname = XDR.offsets_filename(self.traj)
        saved_offsets = XDR.read_numpy_offsets(fname)
        saved_offsets['ctime'] += 1
        with open(fname, 'wb') as f:
            np.savez(f, **saved_offsets)

        with warnings.catch_warnings(record=True) as warn:
            warnings.simplefilter('always')
            self._reader(self.traj)
        assert_equal(warn[0].message.args,
                     ('Reload offsets from trajectory\n ctime or size or n_atoms did not match', ))

    @dec.slow
    def test_persistent_offsets_natoms_mismatch(self):
        # check that stored offsets are not loaded when trajectory
        # ctime differs from stored ctime
        fname = XDR.offsets_filename(self.traj)
        saved_offsets = XDR.read_numpy_offsets(fname)
        saved_offsets['n_atoms'] += 1
        np.savez(fname, **saved_offsets)

        with warnings.catch_warnings(record=True) as warn:
            warnings.simplefilter('always')
            self._reader(self.traj)
        assert_equal(warn[0].message.args,
                     ('Reload offsets from trajectory\n ctime or size or n_atoms did not match', ))

    @dec.slow
    def test_persistent_offsets_last_frame_wrong(self):
        fname = XDR.offsets_filename(self.traj)
        saved_offsets = XDR.read_numpy_offsets(fname)

        idx_frame = 3
        saved_offsets['offsets'][idx_frame] += 42
        np.savez(fname, **saved_offsets)

        with warnings.catch_warnings(record=True) as warn:
            warnings.simplefilter('always')
            reader = self._reader(self.traj)
            reader[idx_frame]

        assert_equal(warn[0].message.args[0],
                     'seek failed, recalculating offsets and retrying')

    @dec.slow
    def test_unsupported_format(self):
        fname = XDR.offsets_filename(self.traj)
        saved_offsets = XDR.read_numpy_offsets(fname)

        idx_frame = 3
        saved_offsets.pop('n_atoms')
        np.savez(fname, **saved_offsets)

        # ok as long as this doesn't throw
        reader = self._reader(self.traj)
        reader[idx_frame]

    @dec.slow
    def test_persistent_offsets_readonly(self):
        os.remove(XDR.offsets_filename(self.traj))
        assert_equal(os.path.exists(
            XDR.offsets_filename(self.trajectory.filename)), False)

        os.chmod(self.tmpdir.name, 0o555)
        self.trajectory._read_offsets(store=True)
        assert_equal(os.path.exists(
            XDR.offsets_filename(self.trajectory.filename)), False)


class TestXTCReader_offsets(_GromacsReader_offsets):
    filename = XTC
    ref_offsets = np.array([0, 165188, 330364, 495520, 660708, 825872, 991044,
                            1156212, 1321384, 1486544])
    _reader = MDAnalysis.coordinates.XTC.XTCReader


class TestTRRReader_offsets(_GromacsReader_offsets):
    filename = TRR
    ref_offsets = np.array([0, 1144464, 2288928, 3433392, 4577856, 5722320,
                            6866784, 8011248, 9155712, 10300176])
    _reader = MDAnalysis.coordinates.TRR.TRRReader
