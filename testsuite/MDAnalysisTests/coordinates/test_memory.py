from numpy.testing import raises

import MDAnalysis as mda
from MDAnalysisTests.datafiles import DCD, PDB_small
from MDAnalysisTests.coordinates.base import (BaseReference,
                                              assert_timestep_almost_equal)
from MDAnalysis.coordinates.memory import MemoryReader
from numpy.testing import assert_equal

from unittest import TestCase


class MemoryReference(BaseReference):
    def __init__(self):
        super(MemoryReference, self).__init__()
        self.universe = mda.Universe(PDB_small, DCD)
        self.trajectory = \
            self.universe.trajectory.timeseries(self.universe.atoms)
        self.n_atoms = self.universe.trajectory.n_atoms
        self.n_frames = self.universe.trajectory.n_frames
        self.topology = PDB_small
        self.reader = mda.coordinates.memory.MemoryReader

        self.first_frame = MemoryReader.MemoryTimestep(self.n_atoms)
        self.first_frame.positions = self.trajectory[:,0,:]
        self.first_frame.frame = 0

        self.second_frame = MemoryReader.MemoryTimestep(self.n_atoms)
        self.second_frame.positions = self.trajectory[:,1,:]
        self.second_frame.frame = 1

        self.last_frame = MemoryReader.MemoryTimestep(self.n_atoms)
        self.last_frame.positions = self.trajectory[:,self.n_frames - 1,:]
        self.last_frame.frame = self.n_frames - 1

        self.jump_to_frame = self.first_frame.copy()
        self.jump_to_frame.positions = self.trajectory[:,3,:]
        self.jump_to_frame.frame = 3


class TestArrayReader(TestCase):
    def setUp(self):
        reference = MemoryReference()
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

    def test_iteration(self):
        frames = 0
        for i, frame in enumerate(self.reader):
            frames += 1
        assert_equal(frames, self.ref.n_frames)

    def test_extract_array_afc(self):
        assert_equal(self.reader.timeseries(format='afc').shape, (3341, 98, 3))

    def test_extract_array_fac(self):
        assert_equal(self.reader.timeseries(format='fac').shape, (98, 3341, 3))

    def test_extract_array_cfa(self):
        assert_equal(self.reader.timeseries(format='cfa').shape, (3, 98, 3341))

    def test_extract_array_acf(self):
        assert_equal(self.reader.timeseries(format='acf').shape, (3341, 3, 98))

    def test_extract_array_fca(self):
        assert_equal(self.reader.timeseries(format='fca').shape, (98, 3, 3341))

    def test_extract_array_caf(self):
        assert_equal(self.reader.timeseries(format='caf').shape, (3, 3341, 98))

    def test_timeseries_skip1(self):
        assert_equal(self.reader.timeseries(self.ref.universe.atoms).shape,
                     (3341, 98, 3))

    def test_timeseries_skip10(self):
        assert_equal(self.reader.timeseries(skip=10).shape,
                     (3341, 9, 3))
        assert_equal(self.ref.universe.trajectory.timeseries(skip=10)[0,:,0],
                     self.reader.timeseries(skip=10)[0,:,0])

    def test_timeseries_view(self):
        assert_equal(self.reader.timeseries().base is self.reader.get_array(),
                     True)

    def test_timeseries_view2(self):
        assert_equal(
            self.reader.timeseries(start=5,
                                   stop=15,
                                   skip=2,
                                   format='fac').base is self.reader.get_array(),
                     True)

    def test_timeseries_noview2(self):
        assert_equal(
            self.reader.timeseries(
                asel=self.ref.universe.
                    select_atoms("name CA")).base is self.reader.get_array(),
                    False)

    def test_repr(self):
        str_rep = str(self.reader)
        expected = "<MemoryReader with 98 frames of 3341 atoms>"
        assert_equal(str_rep, expected)
