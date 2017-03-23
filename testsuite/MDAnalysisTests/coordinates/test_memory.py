import numpy as np

import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysisTests.datafiles import DCD, PSF
from MDAnalysisTests.coordinates.base import (BaseReference,
                                              MultiframeReaderTest)
from MDAnalysis.coordinates.memory import Timestep
from numpy.testing import assert_equal, dec
from MDAnalysisTests import parser_not_found


class MemoryReference(BaseReference):
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
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

        self.first_frame = Timestep(self.n_atoms)
        self.first_frame.positions = np.array(self.universe.trajectory[0])
        self.first_frame.frame = 0
        self.first_frame.time = self.first_frame.frame*self.dt

        self.second_frame = Timestep(self.n_atoms)
        self.second_frame.positions = np.array(self.universe.trajectory[1])
        self.second_frame.frame = 1
        self.second_frame.time = self.second_frame.frame*self.dt

        self.last_frame = Timestep(self.n_atoms)
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


class TestMemoryReader(MultiframeReaderTest):
    def __init__(self):
        reference = MemoryReference()
        super(TestMemoryReader, self).__init__(reference)

    def test_filename_transefer_to_memory(self):
        # MemoryReader should have a filename attribute set to the trajaectory filename
        universe = mda.Universe(PSF, DCD)
        universe.transfer_to_memory()
        assert_equal(universe.trajectory.filename, DCD)

    def test_filename_array(self):
        # filename attribute of MemoryReader should be None when generated from an array
        universe = mda.Universe(PSF, DCD)
        coordinates = universe.trajectory.timeseries(universe.atoms)
        universe2 = mda.Universe(PSF, coordinates, format=MemoryReader, order='afc')
        assert_equal(universe2.trajectory.filename, None)

    def test_default_memory_layout(self):
        universe1 = mda.Universe(PSF, DCD, in_memory=True)
        universe2 = mda.Universe(PSF, DCD, in_memory=True, order='fac')
        assert_equal(universe1.trajectory.get_array().shape,
                     universe2.trajectory.get_array().shape)
        
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
        # timeseries() is expected to provide a view of the underlying array
        assert_equal(self.reader.timeseries().base is self.reader.get_array(),
                     True)

    def test_timeseries_subarray_view(self):
        # timeseries() is expected to provide a view of the underlying array
        # also in the case where we slice the array using the start, stop and
        # step options.
        assert_equal(
            self.reader.timeseries(start=5,
                                   stop=15,
                                   step=2,
                                   format='fac').base is self.reader.get_array(),
                     True)

    def test_timeseries_view_from_universe_atoms(self):
        # timeseries() is expected to provide a view of the underlying array
        # also in the special case when asel=universe.atoms.
        selection = self.ref.universe.atoms
        assert_equal(self.reader.timeseries(
            asel=selection).base is self.reader.get_array(),
            True)

    def test_timeseries_view_from_select_all(self):
        # timeseries() is expected to provide a view of the underlying array
        # also in the special case when using "all" in selections.
        selection = self.ref.universe.select_atoms("all")
        assert_equal(self.reader.timeseries(
            asel=selection).base is self.reader.get_array(),
            True)

    def test_timeseries_noview(self):
        # timeseries() is expected NOT to provide a view of the underlying array
        # for any other selection than "all".
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

    def test_float32(self):
        # Check that we get float32 positions even when initializing with float64
        coordinates = np.random.uniform(size=(100, self.ref.universe.atoms.n_atoms, 3)).cumsum(0)
        universe = mda.Universe(self.ref.universe.filename, coordinates, format=MemoryReader)
        assert_equal(universe.trajectory.get_array().dtype, np.dtype('float32'))
