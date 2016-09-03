import numpy as np
import logging
from numpy.testing import raises

import MDAnalysis as mda
from MDAnalysisTests.datafiles import DCD, PSF
from MDAnalysisTests.coordinates.base import (BaseReference,
                                              BaseReaderTest)
from MDAnalysis.coordinates.memory import MemoryReader
from numpy.testing import assert_equal


class MemoryReference(BaseReference):
    def __init__(self):
        super(MemoryReference, self).__init__()
        
        self.topology = PSF
        self.trajectory = DCD
        self.universe = mda.Universe(PSF, DCD)

        self.n_atoms = self.universe.trajectory.n_atoms
        self.n_frames = self.universe.trajectory.n_frames

        self.dt = self.universe.trajectory.ts.dt
        self.dimensions = self.universe.trajectory.ts.dimensions
        self.totaltime = self.universe.trajectory.totaltime
        self.volume = self.universe.trajectory.ts.volume

        self.first_frame = MemoryReader.MemoryTimestep(self.n_atoms)
        self.first_frame.positions = np.array(self.universe.trajectory[0])
        self.first_frame.frame = 0
        self.first_frame.time = self.first_frame.frame*self.dt

        self.second_frame = MemoryReader.MemoryTimestep(self.n_atoms)
        self.second_frame.positions = np.array(self.universe.trajectory[1])
        self.second_frame.frame = 1
        self.second_frame.time = self.second_frame.frame*self.dt

        self.last_frame = MemoryReader.MemoryTimestep(self.n_atoms)
        self.last_frame.positions = \
            np.array(self.universe.trajectory[self.n_frames - 1])
        self.last_frame.frame = self.n_frames - 1
        self.last_frame.time = self.last_frame.frame*self.dt

        self.jump_to_frame = self.first_frame.copy()
        self.jump_to_frame.positions = np.array(self.universe.trajectory[3])
        self.jump_to_frame.frame = 3
        self.jump_to_frame.time = self.jump_to_frame.frame*self.dt

    def reader(self, trajectory):
        return mda.Universe(self.topology,
                            trajectory, in_memory=True).trajectory

    def iter_ts(self, i):
        ts = self.universe.trajectory[i]
        return ts



class TestMemoryReader(BaseReaderTest):
    def __init__(self):

        reference = MemoryReference()
        super(TestMemoryReader, self).__init__(reference)

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
        # Check that timeseries skip works similar to numpy slicing
        array1 = self.reader.timeseries(step=10)
        array2 = self.reader.timeseries()[:,::10,:]
        assert_equal(array1, array2)

    def test_timeseries_view(self):
        assert_equal(self.reader.timeseries().base is self.reader.get_array(),
                     True)

    def test_timeseries_view2(self):
        assert_equal(
            self.reader.timeseries(start=5,
                                   stop=15,
                                   step=2,
                                   format='fac').base is self.reader.get_array(),
                     True)

    def test_timeseries_view3(self):
        selection = self.ref.universe.atoms
        assert_equal(self.reader.timeseries(
            asel=selection).base is self.reader.get_array(),
            True)

    def test_timeseries_view4(self):
        selection = self.ref.universe.select_atoms("all")
        assert_equal(self.reader.timeseries(
            asel=selection).base is self.reader.get_array(),
            True)

    def test_timeseries_noview(self):
        selection = self.ref.universe.select_atoms("name CA")
        assert_equal(self.reader.timeseries(
            asel=selection).base is self.reader.get_array(),
            False)

    def test_repr(self):
        str_rep = str(self.reader)
        expected = "<MemoryReader with 98 frames of 3341 atoms>"
        assert_equal(str_rep, expected)

    def test_get_writer_1(self):
        pass

    def test_get_writer_2(self):
        pass

