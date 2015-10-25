import MDAnalysis as mda
import numpy as np
from six.moves import zip
from MDAnalysis.coordinates.base import Timestep

from numpy.testing import (assert_equal, assert_raises, assert_almost_equal,
                           assert_array_almost_equal, raises)
from unittest import TestCase
import tempdir

from MDAnalysisTests.coordinates.reference import RefAdKSmall


class _SingleFrameReader(TestCase, RefAdKSmall):
    # see TestPDBReader how to set up!

    def tearDown(self):
        del self.universe

    def test_flag_permissive_pdb_reader(self):
        """test_flag_permissive_pdb_reader: permissive_pdb_reader==True enables
        primitive PDB reader"""
        assert_equal(mda.core.flags['permissive_pdb_reader'], True,
                     "'permissive_pdb_reader' flag should be True as "
                     "MDAnalysis default")

    def test_load_file(self):
        U = self.universe
        assert_equal(len(U.atoms), self.ref_n_atoms,
                     "load Universe from file %s" % U.trajectory.filename)
        assert_equal(U.atoms.select_atoms('resid 150 and name HA2').atoms[0],
                     U.atoms[self.ref_E151HA2_index], "Atom selections")

    def test_n_atoms(self):
        assert_equal(self.universe.trajectory.n_atoms, self.ref_n_atoms,
                     "wrong number of atoms")

    def test_numres(self):
        assert_equal(self.universe.atoms.n_residues, 214,
                     "wrong number of residues")

    def test_n_frames(self):
        assert_equal(self.universe.trajectory.n_frames, 1,
                     "wrong number of frames in pdb")

    def test_time(self):
        assert_equal(self.universe.trajectory.time, 0.0,
                     "wrong time of the frame")

    def test_frame(self):
        assert_equal(
            self.universe.trajectory.frame, 0,
            "wrong frame number (0-based, should be 0 for single frame "
            "readers)")

    def test_frame_index_0(self):
        self.universe.trajectory[0]
        assert_equal(self.universe.trajectory.ts.frame, 0,
                     "frame number for frame index 0 should be 0")

    def test_frame_index_1_raises_IndexError(self):
        def go_to_2(traj=self.universe.trajectory):
            traj[1]

        assert_raises(IndexError, go_to_2)

    def test_dt(self):
        """testing that accessing universe.trajectory.dt gives 1.0
        (the default)"""
        assert_equal(1.0, self.universe.trajectory.dt)

    def test_coordinates(self):
        A10CA = self.universe.atoms.CA[10]
        # restrict accuracy to maximum in PDB files (3 decimals)
        assert_almost_equal(A10CA.pos,
                            self.ref_coordinates['A10CA'],
                            3,
                            err_msg="wrong coordinates for A10:CA")

    def test_distances(self):
        NTERM = self.universe.atoms.N[0]
        CTERM = self.universe.atoms.C[-1]
        d = mda.lib.mdamath.norm(NTERM.position - CTERM.position)
        assert_almost_equal(d,
                            self.ref_distances['endtoend'],
                            self.prec,
                            err_msg="distance between M1:N and G214:C")

    def test_full_slice(self):
        trj_iter = self.universe.trajectory[:]
        frames = [ts.frame for ts in trj_iter]
        assert_equal(frames, np.arange(self.universe.trajectory.n_frames))

    def test_last_slice(self):
        # should be same as above: only 1 frame!
        trj_iter = self.universe.trajectory[-1:]
        frames = [ts.frame for ts in trj_iter]
        assert_equal(frames, np.arange(self.universe.trajectory.n_frames))


class BaseReference(object):
    def __init__(self):
        self.trajectory = None
        self.n_atoms = 5
        self.n_frames = 5
        # default for the numpy test functions
        self.prec = 6

        self.first_frame = Timestep(self.n_atoms)
        self.first_frame.positions = np.arange(
            3 * self.n_atoms).reshape(self.n_atoms, 3)
        self.first_frame.frame = 0

        self.second_frame = self.first_frame.copy()
        self.second_frame.positions = 2 ** 1 * self.first_frame.positions
        self.second_frame.frame = 1

        self.last_frame = self.first_frame.copy()
        self.last_frame.positions = 2 ** 4 * self.first_frame.positions
        self.last_frame.frame = self.n_frames - 1

        # remember frames are 0 indexed
        self.jump_to_frame = self.first_frame.copy()
        self.jump_to_frame.positions = 2 ** 3 * self.first_frame.positions
        self.jump_to_frame.frame = 3

        self.dimensions = np.array([80, 80, 80, 60, 60, 90], dtype=np.float32)
        self.volume = mda.lib.mdamath.box_volume(self.dimensions)
        self.time = 0
        self.dt = 1
        self.totaltime = 5

    def iter_ts(self, i):
        ts = self.first_frame.copy()
        ts.positions = 2**i * self.first_frame.positions
        ts.time = i
        ts.frame = i
        return ts


class BaseReaderTest(object):
    def __init__(self, reference):
        self.ref = reference
        self.reader = self.ref.reader(self.ref.trajectory)

    def test_n_atoms(self):
        assert_equal(self.reader.n_atoms, self.ref.n_atoms)

    def test_n_frames(self):
        assert_equal(len(self.reader), self.ref.n_frames)

    def test_first_frame(self):
        self.reader.rewind()
        assert_timestep_almost_equal(self.reader.ts, self.ref.first_frame,
                                     decimal=self.ref.prec)

    def test_double_close(self):
        self.reader.close()
        self.reader.close()
        self.reader._reopen()

    def test_reopen(self):
        self.reader.close()
        self.reader._reopen()
        ts = self.reader.next()
        assert_timestep_almost_equal(ts, self.ref.first_frame,
                                     decimal=self.ref.prec)

    def test_last_frame(self):
        ts = self.reader[-1]
        assert_timestep_almost_equal(ts, self.ref.last_frame,
                                     decimal=self.ref.prec)

    def test_next_gives_second_frame(self):
        reader = self.ref.reader(self.ref.trajectory)
        ts = reader.next()
        assert_timestep_almost_equal(ts, self.ref.second_frame,
                                     decimal=self.ref.prec)

    @raises(IndexError)
    def test_go_over_last_frame(self):
        self.reader[self.ref.n_frames + 1]

    def test_frame_jump(self):
        ts = self.reader[self.ref.jump_to_frame.frame]
        assert_timestep_almost_equal(ts, self.ref.jump_to_frame,
                                     decimal=self.ref.prec)

    def test_get_writer(self):
        with tempdir.in_tempdir():
            self.outfile = 'test-writer' + self.ref.ext
            with self.reader.Writer(self.outfile) as W:
                assert_equal(isinstance(W, self.ref.writer), True)
                assert_equal(W.n_atoms, self.reader.n_atoms)

    def test_get_writer_2(self):
        with tempdir.in_tempdir():
            self.outfile = 'test-writer' + self.ref.ext
            with self.reader.Writer(self.outfile, n_atoms=100) as W:
                assert_equal(isinstance(W, self.ref.writer), True)
                assert_equal(W.n_atoms, 100)

    def test_dt(self):
        assert_equal(self.reader.dt, self.ref.dt)

    def test_total_time(self):
        assert_equal(self.reader.totaltime, self.ref.totaltime)

    def test_dimensions(self):
        assert_array_almost_equal(self.reader.ts.dimensions,
                                  self.ref.dimensions,
                                  decimal=self.ref.prec)

    def test_volume(self):
        self.reader.rewind()
        vol = self.reader.ts.volume
        assert_array_almost_equal(vol, self.ref.volume)

    def test_iter(self):
        for i, ts in enumerate(self.reader):
            assert_timestep_almost_equal(ts, self.ref.iter_ts(i),
                                         decimal=self.ref.prec)



class BaseWriterTest(object):
    def __init__(self, reference):
        self.ref = reference
        self.tmpdir = tempdir.TempDir()
        self.reader = self.ref.reader(self.ref.trajectory)

    def tmp_file(self, name):
        return self.tmpdir.name + name + '.' + self.ref.ext

    def test_write_trajectory_timestep(self):
        outfile = self.tmp_file('write-timestep-test')
        with self.ref.writer(outfile) as W:
            for ts in self.reader:
                W.write(ts)
        self._check_copy(outfile)

    def test_write_trajectory_atomgroup(self):
        uni = mda.Universe(self.ref.topology, self.ref.trajectory)
        outfile = self.tmp_file('write-atoms-test')
        with self.ref.writer(outfile) as w:
            for ts in uni.trajectory:
                w.write(uni.atoms)
        self._check_copy(outfile)

    def test_write_trajectory_universe(self):
        uni = mda.Universe(self.ref.topology, self.ref.trajectory)
        outfile = self.tmp_file('write-uni-test')
        with self.ref.writer(outfile) as w:
            for ts in uni.trajectory:
                w.write(uni)
        self._check_copy(outfile)

    def test_write_selection(self):
        uni = mda.Universe(self.ref.topology, self.ref.trajectory)
        sel_str = 'resid 1'
        sel = uni.select_atoms(sel_str)
        outfile = self.tmp_file('write-selection-test')

        with self.ref.writer(outfile) as W:
            for ts in uni.trajectory:
                W.write(sel.atoms)

        copy = self.ref.reader(outfile)
        for orig_ts, copy_ts in zip(uni.trajectory, copy):
            assert_array_almost_equal(
                copy_ts._pos, sel.atoms.positions, self.ref.prec,
                err_msg="coordinate mismatch between original and written "
                "trajectory at frame {} (orig) vs {} (copy)".format(
                    orig_ts.frame, copy_ts.frame))

    def _check_copy(self, fname):
        copy = self.ref.reader(fname)
        assert_equal(self.reader.n_frames, copy.n_frames)
        for orig_ts, copy_ts in zip(self.reader, copy):
            assert_timestep_almost_equal(
                copy_ts, orig_ts, decimal=self.ref.prec)

    @raises(mda.NoDataError)
    def test_write_none(self):
        outfile = self.tmp_file('write-none')
        with self.ref.writer(outfile) as w:
            w.write(None)


def assert_timestep_almost_equal(A, B, decimal=6, verbose=True):
    if not isinstance(A, Timestep):
        raise AssertionError('A is not of type Timestep')
    if not isinstance(B, Timestep):
        raise AssertionError('B is not of type Timestep')

    if A.frame != B.frame:
        raise AssertionError('A and B refer to different frames: '
                             'A.frame = {}, B.frame={}'.format(
                                 A.frame, B.frame))

    if A.time != B.time:
        raise AssertionError('A and B refer to different times: '
                             'A.time = {}, B.time={}'.format(
                                 A.time, B.time))

    if A.n_atoms != B.n_atoms:
        raise AssertionError('A and B have a differnent number of atoms: '
                             'A.n_atoms = {}, B.n_atoms = {}'.format(
                                 A.n_atoms, B.n_atoms))

    assert_array_almost_equal(A.positions, B.positions, decimal=decimal,
                              err_msg='Timestep positions', verbose=verbose)

    if A.has_velocities != B.has_velocities:
        raise AssertionError('Only one Timestep has velocities:'
                             'A.has_velocities = {}, B.has_velocities'.format(
                                 A.has_velocities, B.has_velocities))
    if A.has_velocities:
        assert_array_almost_equal(A.velocities, B.velocities, decimal=decimal,
                                  err_msg='Timestep velocities',
                                  verbose=verbose)

    if A.has_forces != B.has_forces:
        raise AssertionError('Only one Timestep has forces:'
                             'A.has_forces = {}, B.has_forces'.format(
                                 A.has_forces, B.has_forces))
    if A.has_forces:
        assert_array_almost_equal(A.forces, B.forces, decimal=decimal,
                                  err_msg='Timestep forces', verbose=verbose)
