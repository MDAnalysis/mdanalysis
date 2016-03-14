from numpy.testing import raises

import MDAnalysis as mda
from MDAnalysisTests.datafiles import DCD, PDB_small
from MDAnalysisTests.coordinates.base import (BaseReference,
                                              assert_timestep_almost_equal)
from MDAnalysis.coordinates.array import ArrayReader
from numpy.testing import assert_equal

from unittest import TestCase


class ArrayReference(BaseReference):
    def __init__(self):
        super(ArrayReference, self).__init__()
        universe = mda.Universe(PDB_small, DCD)
        self.trajectory = universe.trajectory.timeseries(universe.atoms)
        self.n_atoms = universe.trajectory.n_atoms
        self.n_frames = universe.trajectory.n_frames
        self.topology = PDB_small
        self.reader = mda.coordinates.array.ArrayReader

        self.first_frame = ArrayReader.ArrayTimestep(self.n_atoms)
        self.first_frame.positions = self.trajectory[:,0,:]
        self.first_frame.frame = 0

        self.second_frame = ArrayReader.ArrayTimestep(self.n_atoms)
        self.second_frame.positions = self.trajectory[:,1,:]
        self.second_frame.frame = 1

        self.last_frame = ArrayReader.ArrayTimestep(self.n_atoms)
        self.last_frame.positions = self.trajectory[:,self.n_frames - 1,:]
        self.last_frame.frame = self.n_frames - 1

        self.jump_to_frame = self.first_frame.copy()
        self.jump_to_frame.positions = self.trajectory[:,3,:]
        self.jump_to_frame.frame = 3


class TestArrayReader(TestCase):
    def setUp(self):
        reference = ArrayReference()
        self.ref = reference
        self.reader = self.ref.reader(self.ref.trajectory)
        print self.reader.ts

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