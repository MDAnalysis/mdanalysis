# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
import numpy as np

import MDAnalysis as mda
import pytest
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysisTests.datafiles import DCD, PSF
from MDAnalysisTests.coordinates.base import (BaseReference,
                                              MultiframeReaderTest)
from MDAnalysis.coordinates.memory import Timestep
from numpy.testing import assert_equal, assert_almost_equal


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
        # correct time because memory reader doesn't read the correct time
        ts.time = ts.frame * self.dt
        return ts


class TestMemoryReader(MultiframeReaderTest):
    @staticmethod
    @pytest.fixture(scope='class')
    def ref():
        return MemoryReference()

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
        assert universe2.trajectory.filename is None

    def test_default_memory_layout(self):
        universe1 = mda.Universe(PSF, DCD, in_memory=True)
        universe2 = mda.Universe(PSF, DCD, in_memory=True, order='fac')
        assert_equal(universe1.trajectory.get_array().shape,
                     universe2.trajectory.get_array().shape)

    def test_iteration(self, ref, reader):
        frames = 0
        for i, frame in enumerate(reader):
            frames += 1
        assert frames == ref.n_frames

    def test_extract_array_afc(self,reader):
        assert_equal(reader.timeseries(order='afc').shape, (3341, 98, 3))

    def test_extract_array_afc(self, reader):
        assert_equal(reader.timeseries(order='afc').shape, (3341, 98, 3))

    def test_extract_array_fac(self, reader):
        assert_equal(reader.timeseries(order='fac').shape, (98, 3341, 3))

    def test_extract_array_cfa(self, reader):
        assert_equal(reader.timeseries(order='cfa').shape, (3, 98, 3341))

    def test_extract_array_acf(self, reader):
        assert_equal(reader.timeseries(order='acf').shape, (3341, 3, 98))

    def test_extract_array_fca(self, reader):
        assert_equal(reader.timeseries(order='fca').shape, (98, 3, 3341))

    def test_extract_array_caf(self, reader):
        assert_equal(reader.timeseries(order='caf').shape, (3, 3341, 98))

    def test_timeseries_skip1(self, ref, reader):
        assert_equal(reader.timeseries(ref.universe.atoms).shape,
                     (3341, 98, 3))

    def test_timeseries_skip10(self, reader):
        # Check that timeseries skip works similar to numpy slicing
        array1 = reader.timeseries(step=10)
        array2 = reader.timeseries()[:,::10,:]
        assert_equal(array1, array2)

    def test_timeseries_view(self, reader):
        # timeseries() is expected to provide a view of the underlying array
        assert reader.timeseries().base is reader.get_array()

    def test_timeseries_subarray_view(self, reader):
        # timeseries() is expected to provide a view of the underlying array
        # also in the case where we slice the array using the start, stop and
        # step options.
        assert reader.timeseries(start=5,stop=15,step=2,order='fac').base is\
               reader.get_array()

    def test_timeseries_view_from_universe_atoms(self, ref, reader):
        # timeseries() is expected to provide a view of the underlying array
        # also in the special case when asel=universe.atoms.
        selection = ref.universe.atoms
        assert reader.timeseries(asel=selection).base is reader.get_array()

    def test_timeseries_view_from_select_all(self, ref, reader):
        # timeseries() is expected to provide a view of the underlying array
        # also in the special case when using "all" in selections.
        selection = ref.universe.select_atoms("all")
        assert_equal(reader.timeseries(
            asel=selection).base is reader.get_array(),
            True)

    def test_timeseries_noview(self, ref, reader):
        # timeseries() is expected NOT to provide a view of the underlying array
        # for any other selection than "all".
        selection = ref.universe.select_atoms("name CA")
        assert reader.timeseries(asel=selection).base is not reader.get_array()

    def test_repr(self, reader):
        str_rep = str(reader)
        expected = "<MemoryReader with 98 frames of 3341 atoms>"
        assert_equal(str_rep, expected)

    def test_get_writer_1(self):
        pass

    def test_get_writer_2(self):
        pass

    def test_float32(self, ref):
        # Check that we get float32 positions even when initializing with float64
        coordinates = np.random.uniform(size=(100, ref.universe.atoms.n_atoms, 3)).cumsum(0)
        universe = mda.Universe(ref.universe.filename, coordinates, format=MemoryReader)
        assert_equal(universe.trajectory.get_array().dtype, np.dtype('float32'))

    def test_position_assignation(self, reader):
        # When coordinates are assigned to a timestep, is the change persistent?
        new_positions = np.ones_like(reader.ts.positions, dtype=np.float32)
        reader.ts.positions = new_positions
        reader[0]
        assert_almost_equal(reader.ts.positions, new_positions)

    def test_timeseries_warns_deprecation(self, reader):
        with pytest.warns(DeprecationWarning, match="MemoryReader.timeseries "
                          "inclusive"):
            reader.timeseries(start=0, stop=3, step=1)

    def test_timeseries_asel_warns_deprecation(self, ref, reader):
        selection = ref.universe.atoms
        with pytest.warns(DeprecationWarning, match="asel argument to"):
            reader.timeseries(asel=selection)
    
    def test_timeseries_atomgroup(self, ref, reader):
        selection = ref.universe.atoms
        reader.timeseries(atomgroup=selection)
    
    def test_timeseries_atomgroup_asel_mutex(self, ref, reader):
        selection = ref.universe.atoms
        with pytest.raises(ValueError, match="Cannot provide both"):
            reader.timeseries(atomgroup=selection, asel=selection)


class TestMemoryReaderVelsForces(object):
    @staticmethod
    @pytest.fixture(params=['2d', '3d'])
    def ref_pos(request):
        if request.param == '2d':
            return np.arange(30).reshape(10, 3)
        elif request.param == '3d':
            return np.arange(30).reshape(1, 10, 3)

    @staticmethod
    @pytest.fixture(params=['2d', '3d'])
    def ref_vels(request):
        if request.param == '2d':
            return np.arange(30).reshape(10, 3) + 100
        elif request.param == '3d':
            return np.arange(30).reshape(1, 10, 3) + 100

    @staticmethod
    @pytest.fixture(params=['2d', '3d'])
    def ref_forces(request):
        if request.param == '2d':
            return np.arange(30).reshape(10, 3) + 1000
        elif request.param == '3d':
            return np.arange(30).reshape(1, 10, 3) + 1000

    @staticmethod
    def assert_equal_dims(arr1, arr2):
        if arr2.ndim == 3:
            assert_equal(arr1, arr2[0])
        elif arr2.ndim == 2:
            assert_equal(arr1, arr2)

    def test_velocities(self, ref_pos, ref_vels):
        mr = MemoryReader(ref_pos,
                          velocities=ref_vels)

        assert mr.ts.has_velocities
        self.assert_equal_dims(mr.ts.velocities, ref_vels)
        assert not mr.ts.has_forces

    def test_forces(self, ref_pos, ref_forces):
        mr = MemoryReader(ref_pos,
                          forces=ref_forces)

        assert not mr.ts.has_velocities
        assert mr.ts.has_forces
        self.assert_equal_dims(mr.ts.forces, ref_forces)

    def test_both(self, ref_pos, ref_vels, ref_forces):
        mr = MemoryReader(ref_pos,
                          velocities=ref_vels,
                          forces=ref_forces)
        assert mr.ts.has_velocities
        self.assert_equal_dims(mr.ts.velocities, ref_vels)
        assert mr.ts.has_forces
        self.assert_equal_dims(mr.ts.forces, ref_forces)

    @pytest.mark.parametrize('param', ['velocities', 'forces'])
    def test_wrongshape(self, ref_pos, param):
        with pytest.raises(ValueError):
            mr = MemoryReader(ref_pos, **{param: np.zeros((3, 2, 1))})


class TestDimensions(object):
    @staticmethod
    @pytest.fixture
    def ref_pos():
        return np.arange(270).reshape(3, 30, 3)

    @staticmethod
    @pytest.fixture
    def ref_box():
        return np.arange(18).reshape(3, 6)

    def test_single_box(self, ref_pos):
        box = np.array([3, 4, 5, 90, 90, 90])

        mr = MemoryReader(ref_pos, dimensions=box)

        for ts in mr:
            assert_equal(ts.dimensions, box)

    def test_varying_box(self, ref_pos, ref_box):
        mr = MemoryReader(ref_pos, dimensions=ref_box)

        for i, ts in enumerate(mr):
            assert_equal(ts.dimensions, ref_box[i])

    def test_wrong_length(self, ref_pos):
        bad_box = np.arange(12).reshape(2, 6)

        with pytest.raises(ValueError):
            mr = MemoryReader(ref_pos, dimensions=bad_box)

    def test_wrong_shape(self, ref_pos):
        bad_box = np.arange(15).reshape(3, 5)

        with pytest.raises(ValueError):
            mr = MemoryReader(ref_pos, dimensions=bad_box)


class TestMemoryReaderModifications(object):
    # check that modifying MR things behaves as expected
    # in general, modifying the Timestep should be *permanent*
    # this is unlike other Readers!
    n_atoms = 10
    n_frames = 4

    @pytest.fixture()
    def mr_reader(self):
        pos = np.arange(self.n_frames * self.n_atoms * 3).reshape(
            self.n_frames, self.n_atoms, 3)
        vel = np.arange(self.n_frames * self.n_atoms * 3).reshape(
            self.n_frames, self.n_atoms, 3) + 200
        frc = np.arange(self.n_frames * self.n_atoms * 3).reshape(
            self.n_frames, self.n_atoms, 3) + 400
        box = np.arange(self.n_frames * 6).reshape(self.n_frames, 6) + 600

        return MemoryReader(pos,
                            velocities=vel,
                            forces=frc,
                            dimensions=box)

    @pytest.fixture()
    def mr_universe(self, mr_reader):
        u = mda.Universe.empty(self.n_atoms)
        u.trajectory = mr_reader

        return u

    @pytest.mark.parametrize('attr', ['positions', 'velocities', 'forces', 'dimensions'])
    def test_copying(self, mr_reader, attr):
        mr2 = mr_reader.copy()
        # update the attribute
        ts = mr2.ts
        setattr(ts, attr, 7)
        # check the change worked
        assert_almost_equal(getattr(ts, attr), 7)
        assert ts.positions.shape == (self.n_atoms, 3)
        assert ts.velocities.shape == (self.n_atoms, 3)
        assert ts.forces.shape == (self.n_atoms, 3)
        assert ts.dimensions.shape == (6,)
        # move the Reader around, forcing updates of ts
        ts = mr2[2]
        ts = mr2[0]
        # check our old change is still there
        assert_almost_equal(getattr(ts, attr), 7)

    @pytest.mark.parametrize('attr', ['positions', 'velocities', 'forces', 'dimensions'])
    def test_attr_set(self, mr_universe, attr):
        # same as above, but via a Universe/AtomGroup
        u = mr_universe
        ts = u.trajectory[0]

        setattr(ts, attr, 7)

        assert_almost_equal(getattr(ts, attr), 7)

        ts = u.trajectory[2]
        ts = u.trajectory[0]

        assert_almost_equal(getattr(ts, attr), 7)
        assert u.atoms.positions.shape == (self.n_atoms, 3)
        assert u.atoms.velocities.shape == (self.n_atoms, 3)
        assert u.atoms.forces.shape == (self.n_atoms, 3)
        assert u.atoms.dimensions.shape == (6,)

    @pytest.mark.parametrize('attr', ['velocities', 'forces', 'dimensions'])
    def test_non_numpy_arr(self, attr):
        with pytest.raises(TypeError):
            mr = MemoryReader(np.zeros((10, 30, 3)),
                              **{attr: 'not an array'})
