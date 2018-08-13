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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import absolute_import
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

    def test_timeseries_deprecation(self, reader):
        with pytest.warns(DeprecationWarning):
            reader.timeseries(format='fac')

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
