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
import itertools
import pickle

import numpy as np
import pytest
from unittest import TestCase
from numpy.testing import (assert_equal, assert_almost_equal,
                           assert_array_almost_equal, assert_allclose)

import MDAnalysis as mda
from MDAnalysis.coordinates.timestep import Timestep
from MDAnalysis.transformations import translate


from MDAnalysisTests.coordinates.reference import RefAdKSmall
from MDAnalysisTests.datafiles import AUX_XVG_HIGHF, AUX_XVG_LOWF
from MDAnalysisTests import make_Universe


class _SingleFrameReader(TestCase, RefAdKSmall):
    # see TestPDBReader how to set up!
    __test__ = False

    def tearDown(self):
        del self.universe

    def test_load_file(self):
        U = self.universe
        assert_equal(len(U.atoms), self.ref_n_atoms,
                     "load Universe from file {0!s}".format(U.trajectory.filename))
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
        with pytest.raises(IndexError):
            go_to_2()

    def test_dt(self):
        """testing that accessing universe.trajectory.dt gives 1.0
        (the default)"""
        assert_equal(1.0, self.universe.trajectory.dt)

    def test_coordinates(self):
        A10CA = self.universe.select_atoms('name CA')[10]
        # restrict accuracy to maximum in PDB files (3 decimals)
        assert_almost_equal(A10CA.position,
                            self.ref_coordinates['A10CA'],
                            3,
                            err_msg="wrong coordinates for A10:CA")

    def test_distances(self):
        NTERM = self.universe.select_atoms('name N')[0]
        CTERM = self.universe.select_atoms('name C')[-1]
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

    def test_pickle_singleframe_reader(self):
        reader = self.universe.trajectory
        reader_p = pickle.loads(pickle.dumps(reader))
        reader_p_p = pickle.loads(pickle.dumps(reader_p))
        assert_equal(len(reader), len(reader_p))
        assert_equal(reader.ts, reader_p.ts,
                     "Single-frame timestep is changed after pickling")
        assert_equal(len(reader), len(reader_p_p))
        assert_equal(reader.ts, reader_p_p.ts,
                     "Single-frame timestep is changed after double pickling")
        reader.ts.positions[0] = np.array([1, 2, 3])
        reader_p = pickle.loads(pickle.dumps(reader))
        assert_equal(reader.ts, reader_p.ts,
                     "Modification of ts not preserved after serialization")


class BaseReference(object):
    def __init__(self):
        self.trajectory = None
        self.n_atoms = 5
        self.n_frames = 5
        # default for the numpy test functions
        self.prec = 6
        self.container_format = False
        self.changing_dimensions = False

        # for testing auxiliary addition
        self.aux_lowf = AUX_XVG_LOWF  # test auxiliary with lower frequency
        self.aux_lowf_dt = 2  # has steps at 0ps, 2ps, 4ps
        # representative data for each trajectory frame, assuming 'closest' option
        self.aux_lowf_data = [[2 ** 0],  # frame 0 = 0ps = step 0
                              [np.nan],  # frame 1 = 1ps = no step
                              [2 ** 1],  # frame 2 = 2ps = step 1
                              [np.nan],  # frame 3 = 3ps = no step
                              [2 ** 2],  # frame 4 = 4ps = step 2
                              ]
        self.aux_lowf_frames_with_steps = [0, 2, 4]  # trajectory frames with
        # corresponding auxiliary steps

        self.aux_highf = AUX_XVG_HIGHF  # test auxiliary with higher frequency
        self.aux_highf_dt = 0.5  # has steps at 0, 0.5, 1, ... 3.5, 4ps
        self.aux_highf_data = [[2 ** 0],  # frame 0 = 0ps = step 0
                               [2 ** 2],  # frame 1 = 1ps = step 2
                               [2 ** 4],  # frame 2 = 2ps = step 4
                               [2 ** 6],  # frame 3 = 3ps = step 6
                               [2 ** 8],  # frame 4 = 4ps = step 8
                               ]
        self.aux_highf_n_steps = 10
        self.aux_highf_all_data = [[2 ** i] for i in range(self.aux_highf_n_steps)]

        self.aux_offset_by = 0.25

        self.first_frame = Timestep(self.n_atoms)
        self.first_frame.positions = np.arange(
            3 * self.n_atoms).reshape(self.n_atoms, 3)
        self.first_frame.frame = 0
        self.first_frame.aux.lowf = self.aux_lowf_data[0]
        self.first_frame.aux.highf = self.aux_highf_data[0]

        self.second_frame = self.first_frame.copy()
        self.second_frame.positions = 2 ** 1 * self.first_frame.positions
        self.second_frame.frame = 1
        self.second_frame.aux.lowf = self.aux_lowf_data[1]
        self.second_frame.aux.highf = self.aux_highf_data[1]

        self.last_frame = self.first_frame.copy()
        self.last_frame.positions = 2 ** 4 * self.first_frame.positions
        self.last_frame.frame = self.n_frames - 1
        self.last_frame.aux.lowf = self.aux_lowf_data[-1]
        self.last_frame.aux.highf = self.aux_highf_data[-1]

        # remember frames are 0 indexed
        self.jump_to_frame = self.first_frame.copy()
        self.jump_to_frame.positions = 2 ** 3 * self.first_frame.positions
        self.jump_to_frame.frame = 3
        self.jump_to_frame.aux.lowf = self.aux_lowf_data[3]
        self.jump_to_frame.aux.highf = self.aux_highf_data[3]

        self.dimensions = np.array([81.1, 82.2, 83.3, 75, 80, 85],
                                   dtype=np.float32)
        self.dimensions_second_frame = np.array([82.1, 83.2, 84.3, 75.1, 80.1,
                                                 85.1], dtype=np.float32)
        self.volume = mda.lib.mdamath.box_volume(self.dimensions)
        self.time = 0
        self.dt = 1
        self.totaltime = 4

    def iter_ts(self, i):
        ts = self.first_frame.copy()
        ts.positions = 2 ** i * self.first_frame.positions
        ts.time = i
        ts.frame = i
        ts.aux.lowf = np.array(self.aux_lowf_data[i])
        ts.aux.highf = np.array(self.aux_highf_data[i])
        return ts


class BaseReaderTest(object):
    @staticmethod
    @pytest.fixture()
    def reader(ref):
        reader = ref.reader(ref.trajectory)
        reader.add_auxiliary('lowf', ref.aux_lowf, dt=ref.aux_lowf_dt, initial_time=0, time_selector=None)
        reader.add_auxiliary('highf', ref.aux_highf, dt=ref.aux_highf_dt, initial_time=0, time_selector=None)
        return reader

    @staticmethod
    @pytest.fixture()
    def transformed(ref):
        transformed = ref.reader(ref.trajectory)
        transformed.add_transformations(translate([1,1,1]), translate([0,0,0.33]))
        return transformed

    def test_n_atoms(self, ref, reader):
        assert_equal(reader.n_atoms, ref.n_atoms)

    def test_n_frames(self, ref, reader):
        assert_equal(len(reader), ref.n_frames)

    def test_first_frame(self, ref, reader):
        reader.rewind()
        assert_timestep_almost_equal(reader.ts, ref.first_frame,
                                     decimal=ref.prec)

    def test_double_close(self, reader):
        reader.close()
        reader.close()
        reader._reopen()

    def test_get_writer_1(self, ref, reader, tmpdir):
        with tmpdir.as_cwd():
            outfile = 'test-writer.' + ref.ext
            with reader.Writer(outfile) as W:
                assert_equal(isinstance(W, ref.writer), True)
                assert_equal(W.n_atoms, reader.n_atoms)

    def test_get_writer_2(self, ref, reader, tmpdir):
        with tmpdir.as_cwd():
            outfile = 'test-writer.' + ref.ext
            with reader.Writer(outfile, n_atoms=100) as W:
                assert_equal(isinstance(W, ref.writer), True)
                assert_equal(W.n_atoms, 100)

    def test_dt(self, ref, reader):
        assert_almost_equal(reader.dt, ref.dt, decimal=ref.prec)

    def test_ts_dt_matches_reader(self, reader):
        assert_equal(reader.ts.dt, reader.dt)

    def test_total_time(self, ref, reader):
        assert_almost_equal(reader.totaltime, ref.totaltime, decimal=ref.prec)

    def test_first_dimensions(self, ref, reader):
        reader.rewind()
        if ref.dimensions is None:
            assert reader.ts.dimensions is None
        else:
            assert_array_almost_equal(reader.ts.dimensions,
                                      ref.dimensions,
                                      decimal=ref.prec)

    def test_changing_dimensions(self, ref, reader):
        if ref.changing_dimensions:
            reader.rewind()
            if ref.dimensions is None:
                assert reader.ts.dimensions is None
            else:
                assert_array_almost_equal(reader.ts.dimensions,
                                          ref.dimensions,
                                          decimal=ref.prec)
            reader[1]
            if ref.dimensions_second_frame is None:
                assert reader.ts.dimensions is None
            else:
                assert_array_almost_equal(reader.ts.dimensions,
                                          ref.dimensions_second_frame,
                                          decimal=ref.prec)

    def test_volume(self, ref, reader):
        reader.rewind()
        vol = reader.ts.volume
        # Here we can only be sure about the numbers upto the decimal point due
        # to floating point impressions.
        assert_almost_equal(vol, ref.volume, 0)

    def test_iter(self, ref, reader):
        for i, ts in enumerate(reader):
            assert_timestep_almost_equal(ts, ref.iter_ts(i),
                                         decimal=ref.prec)

    def test_add_same_auxname_raises_ValueError(self, ref, reader):
        with pytest.raises(ValueError):
            reader.add_auxiliary('lowf', ref.aux_lowf)

    def test_remove_auxiliary(self, reader):
        reader.remove_auxiliary('lowf')
        with pytest.raises(AttributeError):
            getattr(reader._auxs, 'lowf')
        with pytest.raises(AttributeError):
            getattr(reader.ts.aux, 'lowf')

    def test_remove_nonexistant_auxiliary_raises_ValueError(self, reader):
        with pytest.raises(ValueError):
            reader.remove_auxiliary('nonexistant')

    def test_iter_auxiliary(self, ref, reader):
        # should go through all steps in 'highf'
        for i, auxstep in enumerate(reader.iter_auxiliary('highf')):
            assert_almost_equal(auxstep.data, ref.aux_highf_all_data[i],
                                err_msg="Auxiliary data does not match for "
                                        "step {}".format(i))

    def test_get_aux_attribute(self, ref, reader):
        assert_equal(reader.get_aux_attribute('lowf', 'dt'),
                     ref.aux_lowf_dt)

    def test_iter_as_aux_cutoff(self, ref, reader):
        # load an auxiliary with the same dt but offset from trajectory, and a
        # cutoff of 0
        reader.add_auxiliary('offset', ref.aux_lowf,
                                  dt=ref.dt, time_selector=None,
                                  initial_time=ref.aux_offset_by,
                                  cutoff=0)
        # no auxiliary steps will fall within the cutoff for any frame, so
        # iterating using iter_as_aux should give us nothing
        num_frames = len([i for i in reader.iter_as_aux('offset')])
        assert_equal(num_frames, 0, "iter_as_aux should iterate over 0 frames,"
                                    " not {}".format(num_frames))

    def test_reload_auxiliaries_from_description(self, ref, reader):
        # get auxiliary desscriptions form existing reader
        descriptions = reader.get_aux_descriptions()
        # load a new reader, without auxiliaries
        reader = ref.reader(ref.trajectory)
        # load auxiliaries into new reader, using description...
        for aux in descriptions:
            reader.add_auxiliary(**aux)
        # should have the same number of auxiliaries
        assert_equal(reader.aux_list, reader.aux_list,
                     'Number of auxiliaries does not match')
        # each auxiliary should be the same
        for auxname in reader.aux_list:
            assert_equal(reader._auxs[auxname], reader._auxs[auxname],
                         'AuxReaders do not match')

    def test_stop_iter(self, reader):
        # reset to 0
        reader.rewind()
        for ts in reader[:-1]:
            pass
        assert_equal(reader.frame, 0)

    def test_transformations_iter(self, ref, transformed):
        # test one iteration and see if the transformations are applied
        v1 = np.float32((1,1,1))
        v2 = np.float32((0,0,0.33))
        for i, ts in enumerate(transformed):
            idealcoords = ref.iter_ts(i).positions + v1 + v2
            assert_array_almost_equal(ts.positions, idealcoords, decimal=ref.prec)

    def test_transformations_2iter(self, ref, transformed):
        # Are the transformations applied and
        # are the coordinates "overtransformed"?
        v1 = np.float32((1,1,1))
        v2 = np.float32((0,0,0.33))
        idealcoords=[]
        for i, ts in enumerate(transformed):
            idealcoords.append(ref.iter_ts(i).positions + v1 + v2)
            assert_array_almost_equal(ts.positions, idealcoords[i], decimal=ref.prec)

        for i, ts in enumerate(transformed):
            assert_almost_equal(ts.positions, idealcoords[i], decimal=ref.prec)

    def test_transformations_slice(self, ref, transformed):
        # Are the transformations applied when iterating over a slice of the trajectory?
        v1 = np.float32((1,1,1))
        v2 = np.float32((0,0,0.33))
        for i,ts in enumerate(transformed[2:3:1]):
            idealcoords = ref.iter_ts(ts.frame).positions + v1 + v2
            assert_array_almost_equal(ts.positions, idealcoords, decimal = ref.prec)

    def test_transformations_switch_frame(self, ref, transformed):
        # This test checks if the transformations are applied and if the coordinates
        # "overtransformed" on different situations
        # Are the transformations applied when we switch to a different frame?
        v1 = np.float32((1,1,1))
        v2 = np.float32((0,0,0.33))
        first_ideal = ref.iter_ts(0).positions + v1 + v2
        if len(transformed)>1:
            assert_array_almost_equal(transformed[0].positions, first_ideal, decimal = ref.prec)
            second_ideal = ref.iter_ts(1).positions + v1 + v2
            assert_array_almost_equal(transformed[1].positions, second_ideal, decimal = ref.prec)

            # What if we comeback to the previous frame?
            assert_array_almost_equal(transformed[0].positions, first_ideal, decimal = ref.prec)

            # How about we switch the frame to itself?
            assert_array_almost_equal(transformed[0].positions, first_ideal, decimal = ref.prec)
        else:
            assert_array_almost_equal(transformed[0].positions, first_ideal, decimal = ref.prec)

    def test_transformation_rewind(self,ref, transformed):
        # this test checks if the transformations are applied after rewinding the
        # trajectory
        v1 = np.float32((1,1,1))
        v2 = np.float32((0,0,0.33))
        ideal_coords = ref.iter_ts(0).positions + v1 + v2
        transformed.rewind()
        assert_array_almost_equal(transformed[0].positions, ideal_coords, decimal = ref.prec)

    def test_transformations_copy(self,ref,transformed):
        # this test checks if transformations are carried over a copy and if the
        # coordinates of the copy are equal to the original's
        v1 = np.float32((1,1,1))
        v2 = np.float32((0,0,0.33))
        new = transformed.copy()
        assert_equal(transformed.transformations, new.transformations,
                     "transformations are not equal")
        for i, ts in enumerate(new):
            ideal_coords = ref.iter_ts(i).positions + v1 + v2
            assert_array_almost_equal(ts.positions, ideal_coords, decimal = ref.prec)

    def test_add_another_transformations_raises_ValueError(self, transformed):
        # After defining the transformations, the workflow cannot be changed
        with pytest.raises(ValueError):
            transformed.add_transformations(translate([2,2,2]))

    def test_pickle_reader(self, reader):
        reader_p = pickle.loads(pickle.dumps(reader))
        assert_equal(len(reader), len(reader_p))
        assert_equal(reader.ts, reader_p.ts,
                     "Timestep is changed after pickling")
        reader_p_p = pickle.loads(pickle.dumps(reader_p))
        assert_equal(len(reader), len(reader_p_p))
        assert_equal(reader.ts, reader_p_p.ts,
                     "Timestep is changed after double pickling")
        reader.ts.positions[0] = np.array([1, 2, 3])
        reader_p = pickle.loads(pickle.dumps(reader))
        assert_equal(reader.ts, reader_p.ts,
                     "Modification of ts not preserved after serialization")

    def test_frame_collect_all_same(self, reader):
        # check that the timestep resets so that the base reference is the same
        # for all timesteps in a collection with the exception of memoryreader
        # and DCDReader
        if isinstance(reader, mda.coordinates.memory.MemoryReader):
            pytest.xfail("memoryreader allows independent coordinates")
        if isinstance(reader, mda.coordinates.DCD.DCDReader):
            pytest.xfail("DCDReader allows independent coordinates."
                          "This behaviour is deprecated and will be changed"
                          "in 3.0")
        collected_ts = []
        for i, ts in enumerate(reader):
            collected_ts.append(ts.positions)
        for array  in collected_ts:
            assert_allclose(array, collected_ts[0])

    @pytest.mark.parametrize('order', ('fac', 'fca', 'afc', 'acf', 'caf', 'cfa'))
    def test_timeseries_shape(self, reader, order):
        timeseries = reader.timeseries(order=order)
        a_index = order.index('a')
        f_index = order.index('f')
        c_index = order.index('c')
        assert(timeseries.shape[a_index] == reader.n_atoms)
        assert(timeseries.shape[f_index] == len(reader))
        assert(timeseries.shape[c_index] == 3)

    @pytest.mark.parametrize('slice', ([0,2,1], [0,10,2], [0,10,3]))
    def test_timeseries_values(self, reader, slice):
        ts_positions = []
        if isinstance(reader, mda.coordinates.memory.MemoryReader):
            pytest.xfail("MemoryReader uses deprecated stop inclusive"
                         " indexing, see Issue #3893")
        if slice[1] > len(reader):
            pytest.skip("too few frames in reader")
        for i in range(slice[0], slice[1], slice[2]):
            ts = reader[i]
            ts_positions.append(ts.positions.copy())
        positions = np.asarray(ts_positions)
        timeseries = reader.timeseries(start=slice[0], stop=slice[1],
                                       step=slice[2], order='fac')
        assert_allclose(timeseries, positions)

    @pytest.mark.parametrize('asel', ("index 1", "index 2", "index 1 to 3"))
    def test_timeseries_asel_shape(self, reader, asel):
        atoms = mda.Universe(reader.filename).select_atoms(asel)
        timeseries = reader.timeseries(atoms, order='fac')
        assert(timeseries.shape[0] == len(reader))
        assert(timeseries.shape[1] == len(atoms))
        assert(timeseries.shape[2] == 3)

    def test_timeseries_empty_asel(self, reader):
        with pytest.warns(UserWarning,
                          match="Empty string to select atoms, empty group returned."):
            atoms = mda.Universe(reader.filename).select_atoms(None)
        with pytest.raises(ValueError, match="Timeseries requires at least"):
            reader.timeseries(asel=atoms)

    def test_timeseries_empty_atomgroup(self, reader):
        with pytest.warns(UserWarning,
                          match="Empty string to select atoms, empty group returned."):
            atoms = mda.Universe(reader.filename).select_atoms(None)
        with pytest.raises(ValueError, match="Timeseries requires at least"):
            reader.timeseries(atomgroup=atoms)

    def test_timeseries_asel_warns_deprecation(self, reader):
        atoms = mda.Universe(reader.filename).select_atoms("index 1")
        with pytest.warns(DeprecationWarning, match="asel argument to"):
            timeseries = reader.timeseries(asel=atoms, order='fac')

    def test_timeseries_atomgroup(self, reader):
        atoms = mda.Universe(reader.filename).select_atoms("index 1")
        timeseries = reader.timeseries(atomgroup=atoms, order='fac')

    def test_timeseries_atomgroup_asel_mutex(self, reader):
        atoms = mda.Universe(reader.filename).select_atoms("index 1")
        with pytest.raises(ValueError, match="Cannot provide both"):
            timeseries = reader.timeseries(atomgroup=atoms, asel=atoms, order='fac')


class MultiframeReaderTest(BaseReaderTest):
    def test_last_frame(self, ref, reader):
        ts = reader[-1]
        assert_timestep_almost_equal(ts, ref.last_frame,
                                     decimal=ref.prec)

    def test_go_over_last_frame(self, ref, reader):
        with pytest.raises(IndexError):
            reader[ref.n_frames + 1]

    def test_frame_jump(self, ref, reader):
        ts = reader[ref.jump_to_frame.frame]
        assert_timestep_almost_equal(ts, ref.jump_to_frame,
                                     decimal=ref.prec)

    def test_frame_jump_issue1942(self, ref, reader):
        """Test for issue 1942 (especially XDR on macOS)"""
        reader.rewind()
        try:
            for ii in range(ref.n_frames + 2):
                reader[0]
        except StopIteration:
            pytest.fail("Frame-seeking wrongly iterated (#1942)")

    def test_next_gives_second_frame(self, ref, reader):
        reader = ref.reader(ref.trajectory)
        ts = reader.next()
        assert_timestep_almost_equal(ts, ref.second_frame,
                                     decimal=ref.prec)

    def test_reopen(self, ref, reader):
        reader.close()
        reader._reopen()
        ts = reader.next()
        assert_timestep_almost_equal(ts, ref.first_frame,
                                     decimal=ref.prec)

    def test_rename_aux(self, ref, reader):
        reader.rename_aux('lowf', 'lowf_renamed')
        # data should now be in aux namespace under new name
        assert_equal(reader.ts.aux.lowf_renamed,
                    ref.aux_lowf_data[0])
        # old name should be removed
        with pytest.raises(AttributeError):
            getattr(reader.ts.aux, 'lowf')
        # new name should be retained
        next(reader)
        assert_equal(reader.ts.aux.lowf_renamed,
                     ref.aux_lowf_data[1])

    def test_iter_as_aux_highf(self, ref, reader):
        # auxiliary has a higher frequency, so iter_as_aux should behave the
        # same as regular iteration over the trjectory
        for i, ts in enumerate(reader.iter_as_aux('highf')):
            assert_timestep_almost_equal(ts, ref.iter_ts(i),
                                         decimal=ref.prec)

    def test_iter_as_aux_lowf(self, ref, reader):
        # auxiliary has a lower frequency, so iter_as_aux should iterate over
        # only frames where there is a corresponding auxiliary value
        for i, ts in enumerate(reader.iter_as_aux('lowf')):
            assert_timestep_almost_equal(ts,
                                         ref.iter_ts(ref.aux_lowf_frames_with_steps[i]),
                                         decimal=ref.prec)

    @pytest.mark.parametrize("accessor", [
              lambda traj: traj[[0, 1, 2]],
              lambda traj: traj[:3],
              lambda traj: traj],
              ids=["indexed", "sliced", "all"])
    def test_iter_rewinds(self, reader, accessor):
        for ts_indices in accessor(reader):
            pass
        assert_equal(ts_indices.frame, 0)

    #  To make sure we not only save the current timestep information,
    #  but also maintain its relative position.
    def test_pickle_next_ts_reader(self, reader):
        reader_p = pickle.loads(pickle.dumps(reader))
        assert_equal(next(reader), next(reader_p),
                     "Next timestep is changed after pickling")
        reader_p_p = pickle.loads(pickle.dumps(reader_p))
        assert_equal(next(reader), next(reader_p_p),
                     "Next timestep is changed after double pickling")

    #  To make sure pickle works for last frame.
    def test_pickle_last_ts_reader(self, reader):
        #  move current ts to last frame.
        reader[-1]
        reader_p = pickle.loads(pickle.dumps(reader))
        assert_equal(len(reader), len(reader_p),
                     "Last timestep is changed after pickling")
        assert_equal(reader.ts, reader_p.ts,
                     "Last timestep is changed after pickling")


class BaseWriterTest(object):
    @staticmethod
    @pytest.fixture()
    def reader(ref):
        return ref.reader(ref.trajectory)

    @staticmethod
    @pytest.fixture()
    def u_no_resnames():
        return make_Universe(['names', 'resids'], trajectory=True)

    @staticmethod
    @pytest.fixture()
    def u_no_resids():
        return make_Universe(['names', 'resnames'], trajectory=True)

    @staticmethod
    @pytest.fixture()
    def u_no_names():
        return make_Universe(['resids', 'resnames'],
                             trajectory=True)

    @staticmethod
    @pytest.fixture()
    def universe(ref):
        return mda.Universe(ref.topology, ref.trajectory)

    def test_write_different_box(self, ref, universe, tmpdir):
        if ref.changing_dimensions:
            outfile = 'write-dimensions-test' + ref.ext
            with tmpdir.as_cwd():
                with ref.writer(outfile, universe.atoms.n_atoms) as W:
                    for ts in universe.trajectory:
                        universe.dimensions[:3] += 1
                        W.write(universe)

                written = ref.reader(outfile)

                for ts_ref, ts_w in zip(universe.trajectory, written):
                    universe.dimensions[:3] += 1
                    assert_array_almost_equal(universe.dimensions,
                                              ts_w.dimensions,
                                              decimal=ref.prec)

    def test_write_trajectory_atomgroup(self, ref,reader, universe, tmpdir):
        outfile = 'write-atoms-test.' + ref.ext
        with tmpdir.as_cwd():
            with ref.writer(outfile, universe.atoms.n_atoms) as w:
                for ts in universe.trajectory:
                    w.write(universe.atoms)
            self._check_copy(outfile, ref, reader)

    def test_write_trajectory_universe(self, ref, reader, universe, tmpdir):
        outfile = 'write-uni-test.' + ref.ext
        with tmpdir.as_cwd():
            with ref.writer(outfile, universe.atoms.n_atoms) as w:
                for ts in universe.trajectory:
                    w.write(universe)
            self._check_copy(outfile, ref, reader)

    def test_write_selection(self, ref, reader, universe, u_no_resnames,
                             u_no_resids, u_no_names, tmpdir):
        sel_str = 'resid 1'
        sel = universe.select_atoms(sel_str)
        outfile = 'write-selection-test.' + ref.ext


        with tmpdir.as_cwd():
            with ref.writer(outfile, sel.n_atoms) as W:
                for ts in universe.trajectory:
                    W.write(sel.atoms)

            copy = ref.reader(outfile)
            for orig_ts, copy_ts in zip(universe.trajectory, copy):
                assert_array_almost_equal(
                    copy_ts._pos, sel.atoms.positions, ref.prec,
                    err_msg="coordinate mismatch between original and written "
                            "trajectory at frame {} (orig) vs {} (copy)".format(
                        orig_ts.frame, copy_ts.frame))

    def _check_copy(self, fname, ref, reader):
        copy = ref.reader(fname)
        assert_equal(reader.n_frames, copy.n_frames)
        for orig_ts, copy_ts in zip(reader, copy):
            assert_timestep_almost_equal(
                copy_ts, orig_ts, decimal=ref.prec)

    def test_write_none(self, ref, tmpdir):
        outfile = 'write-none.' + ref.ext
        with tmpdir.as_cwd():
            with pytest.raises(TypeError):
                with ref.writer(outfile, 42) as w:
                    w.write(None)

    def test_no_container(self, ref, tmpdir):
        with tmpdir.as_cwd():
            if ref.container_format:
                ref.writer('foo')
            else:
                with pytest.raises(TypeError):
                    ref.writer('foo')

    def test_write_not_changing_ts(self, ref, universe, tmpdir):
        outfile = 'write-not-changing-ts.' + ref.ext

        copy_ts = universe.trajectory.ts.copy()
        with tmpdir.as_cwd():
            with ref.writer(outfile, n_atoms=5) as W:
                W.write(universe)
                assert_timestep_almost_equal(copy_ts, universe.trajectory.ts)


def assert_timestep_equal(A, B, msg=''):
    """ assert that two timesteps are exactly equal and commutative
    """
    assert A == B, msg
    assert B == A, msg


def assert_timestep_almost_equal(A, B, decimal=6, verbose=True):
    if not isinstance(A, mda.coordinates.timestep.Timestep):
        raise AssertionError('A is not of type Timestep')
    if not isinstance(B, mda.coordinates.timestep.Timestep):
        raise AssertionError('B is not of type Timestep')

    if A.frame != B.frame:
        raise AssertionError('A and B refer to different frames: '
                             'A.frame = {}, B.frame={}'.format(
                                 A.frame, B.frame))

    if not np.allclose(A.time, B.time):
        raise AssertionError('A and B refer to different times: '
                             'A.time = {}, B.time={}'.format(
                                 A.time, B.time))

    if A.n_atoms != B.n_atoms:
        raise AssertionError('A and B have a differnent number of atoms: '
                             'A.n_atoms = {}, B.n_atoms = {}'.format(
                                 A.n_atoms, B.n_atoms))

    if A.has_positions != B.has_positions:
        raise AssertionError('Only one Timestep has positions:'
                             'A.has_positions = {}, B.has_positions = {}'.format(
                                 A.has_positions, B.has_positions))

    if A.has_positions:
        assert_array_almost_equal(A.positions, B.positions, decimal=decimal,
                                  err_msg='Timestep positions',
                                  verbose=verbose)

    if A.has_velocities != B.has_velocities:
        raise AssertionError('Only one Timestep has velocities:'
                             'A.has_velocities = {}, B.has_velocities = {}'.format(
                                 A.has_velocities, B.has_velocities))
    if A.has_velocities:
        assert_array_almost_equal(A.velocities, B.velocities, decimal=decimal,
                                  err_msg='Timestep velocities',
                                  verbose=verbose)

    if A.has_forces != B.has_forces:
        raise AssertionError('Only one Timestep has forces:'
                             'A.has_forces = {}, B.has_forces = {}'.format(
                                 A.has_forces, B.has_forces))
    if A.has_forces:
        assert_array_almost_equal(A.forces, B.forces, decimal=decimal,
                                  err_msg='Timestep forces', verbose=verbose)

    # Check we've got auxiliaries before comparing values (auxiliaries aren't written
    # so we won't have aux values to compare when testing writer)
    if len(A.aux) > 0 and len(B.aux) > 0:
        assert_equal(A.aux, B.aux, err_msg='Auxiliary values do not match: '
                                  'A.aux = {}, B.aux = {}'.format(A.aux, B.aux))
