# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.mdanalysis.org
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
import itertools
import numpy as np
from six.moves import zip, range
from nose.plugins.attrib import attr
from unittest import TestCase
from numpy.testing import (assert_equal, assert_raises, assert_almost_equal,
                           assert_array_almost_equal, raises, assert_allclose,
                           assert_, dec)

import MDAnalysis as mda
from MDAnalysis.coordinates.base import Timestep
from MDAnalysis import NoDataError
from MDAnalysis.lib.mdamath import triclinic_vectors

from MDAnalysisTests.coordinates.reference import RefAdKSmall
from MDAnalysisTests.datafiles import AUX_XVG_HIGHF, AUX_XVG_LOWF
from MDAnalysisTests import tempdir


class _SingleFrameReader(TestCase, RefAdKSmall):
    # see TestPDBReader how to set up!

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

        assert_raises(IndexError, go_to_2)

    def test_dt(self):
        """testing that accessing universe.trajectory.dt gives 1.0
        (the default)"""
        assert_equal(1.0, self.universe.trajectory.dt)

    def test_coordinates(self):
        A10CA = self.universe.atoms.CA[10]
        # restrict accuracy to maximum in PDB files (3 decimals)
        assert_almost_equal(A10CA.position,
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
    def __init__(self, reference):
        self.ref = reference
        self.reader = self.ref.reader(self.ref.trajectory)
        self.reader.add_auxiliary('lowf', self.ref.aux_lowf,
                                  dt=self.ref.aux_lowf_dt,
                                  initial_time=0, time_selector=None)
        self.reader.add_auxiliary('highf', self.ref.aux_highf,
                                  dt=self.ref.aux_highf_dt,
                                  initial_time=0, time_selector=None)

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

    def test_get_writer_1(self):
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

    def test_ts_dt_matches_reader(self):
        assert_equal(self.reader.ts.dt, self.reader.dt)

    def test_total_time(self):
        assert_equal(self.reader.totaltime, self.ref.totaltime)

    def test_first_dimensions(self):
        self.reader.rewind()
        assert_array_almost_equal(self.reader.ts.dimensions,
                                  self.ref.dimensions,
                                  decimal=self.ref.prec)

    def test_changing_dimensions(self):
        if self.ref.changing_dimensions:
            self.reader.rewind()
            assert_array_almost_equal(self.reader.ts.dimensions,
                                      self.ref.dimensions,
                                      decimal=self.ref.prec)
            self.reader[1]
            assert_array_almost_equal(self.reader.ts.dimensions,
                                      self.ref.dimensions_second_frame,
                                      decimal=self.ref.prec)

    def test_volume(self):
        self.reader.rewind()
        vol = self.reader.ts.volume
        # Here we can only be sure about the numbers upto the decimal point due
        # to floating point impressions.
        assert_almost_equal(vol, self.ref.volume, 0)

    def test_iter(self):
        for i, ts in enumerate(self.reader):
            assert_timestep_almost_equal(ts, self.ref.iter_ts(i),
                                         decimal=self.ref.prec)

    @raises(ValueError)
    def test_add_same_auxname_raises_ValueError(self):
        self.reader.add_auxiliary('lowf', self.ref.aux_lowf)

    def test_remove_auxiliary(self):
        self.reader.remove_auxiliary('lowf')
        assert_raises(AttributeError, getattr, self.reader._auxs, 'lowf')
        assert_raises(AttributeError, getattr, self.reader.ts.aux, 'lowf')

    @raises(ValueError)
    def test_remove_nonexistant_auxiliary_raises_ValueError(self):
        self.reader.remove_auxiliary('nonexistant')

    def test_iter_auxiliary(self):
        # should go through all steps in 'highf'
        for i, auxstep in enumerate(self.reader.iter_auxiliary('highf')):
            assert_almost_equal(auxstep.data, self.ref.aux_highf_all_data[i],
                                err_msg="Auxiliary data does not match for "
                                        "step {}".format(i))

    def test_get_aux_attribute(self):
        assert_equal(self.reader.get_aux_attribute('lowf', 'dt'),
                     self.ref.aux_lowf_dt)

    def test_iter_as_aux_cutoff(self):
        # load an auxiliary with the same dt but offset from trajectory, and a
        # cutoff of 0
        self.reader.add_auxiliary('offset', self.ref.aux_lowf,
                                  dt=self.ref.dt, time_selector=None,
                                  initial_time=self.ref.aux_offset_by,
                                  cutoff=0)
        # no auxiliary steps will fall within the cutoff for any frame, so
        # iterating using iter_as_aux should give us nothing
        num_frames = len([i for i in self.reader.iter_as_aux('offset')])
        assert_equal(num_frames, 0, "iter_as_aux should iterate over 0 frames,"
                                    " not {}".format(num_frames))

    def test_reload_auxiliaries_from_description(self):
        # get auxiliary desscriptions form existing reader
        descriptions = self.reader.get_aux_descriptions()
        # load a new reader, without auxiliaries
        reader = self.ref.reader(self.ref.trajectory)
        # load auxiliaries into new reader, using description...
        for aux in descriptions:
            reader.add_auxiliary(**aux)
        # should have the same number of auxiliaries
        assert_equal(reader.aux_list, self.reader.aux_list,
                     'Number of auxiliaries does not match')
        # each auxiliary should be the same
        for auxname in reader.aux_list:
            assert_equal(reader._auxs[auxname], self.reader._auxs[auxname],
                         'AuxReaders do not match')

    def test_stop_iter(self):
        # reset to 0
        self.reader.rewind()
        for ts in self.reader[:-1]:
            pass
        assert_equal(self.reader.frame, 0)


class MultiframeReaderTest(BaseReaderTest):
    def test_last_frame(self):
        ts = self.reader[-1]
        assert_timestep_almost_equal(ts, self.ref.last_frame,
                                     decimal=self.ref.prec)

    @raises(IndexError)
    def test_go_over_last_frame(self):
        self.reader[self.ref.n_frames + 1]

    def test_frame_jump(self):
        ts = self.reader[self.ref.jump_to_frame.frame]
        assert_timestep_almost_equal(ts, self.ref.jump_to_frame,
                                     decimal=self.ref.prec)

    def test_next_gives_second_frame(self):
        reader = self.ref.reader(self.ref.trajectory)
        ts = reader.next()
        assert_timestep_almost_equal(ts, self.ref.second_frame,
                                     decimal=self.ref.prec)

    def test_reopen(self):
        self.reader.close()
        self.reader._reopen()
        ts = self.reader.next()
        assert_timestep_almost_equal(ts, self.ref.first_frame,
                                     decimal=self.ref.prec)

    def test_rename_aux(self):
        self.reader.rename_aux('lowf', 'lowf_renamed')
        # data should now be in aux namespace under new name
        assert_equal(self.reader.ts.aux.lowf_renamed,
                     self.ref.aux_lowf_data[0])
        # old name should be removed
        assert_raises(AttributeError, getattr, self.reader.ts.aux, 'lowf')
        # new name should be retained
        next(self.reader)
        assert_equal(self.reader.ts.aux.lowf_renamed,
                     self.ref.aux_lowf_data[1])

    def test_iter_as_aux_highf(self):
        # auxiliary has a higher frequency, so iter_as_aux should behave the
        # same as regular iteration over the trjectory
        for i, ts in enumerate(self.reader.iter_as_aux('highf')):
            assert_timestep_almost_equal(ts, self.ref.iter_ts(i),
                                         decimal=self.ref.prec)

    def test_iter_as_aux_lowf(self):
        # auxiliary has a lower frequency, so iter_as_aux should iterate over
        # only frames where there is a corresponding auxiliary value
        for i, ts in enumerate(self.reader.iter_as_aux('lowf')):
            assert_timestep_almost_equal(ts,
                                         self.ref.iter_ts(self.ref.aux_lowf_frames_with_steps[i]),
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
        with self.ref.writer(outfile, self.reader.n_atoms) as W:
            for ts in self.reader:
                W.write(ts)
        self._check_copy(outfile)

    def test_write_different_box(self):
        if self.ref.changing_dimensions:
            outfile = self.tmp_file('write-dimensions-test')
            with self.ref.writer(outfile, self.reader.n_atoms) as W:
                for ts in self.reader:
                    ts.dimensions[:3] += 1
                    W.write(ts)

            written = self.ref.reader(outfile)

            for ts_ref, ts_w in zip(self.reader, written):
                ts_ref.dimensions[:3] += 1
                assert_array_almost_equal(ts_ref.dimensions,
                                          ts_w.dimensions,
                                          decimal=self.ref.prec)

    def test_write_trajectory_atomgroup(self):
        uni = mda.Universe(self.ref.topology, self.ref.trajectory)
        outfile = self.tmp_file('write-atoms-test')
        with self.ref.writer(outfile, uni.atoms.n_atoms) as w:
            for ts in uni.trajectory:
                w.write(uni.atoms)
        self._check_copy(outfile)

    def test_write_trajectory_universe(self):
        uni = mda.Universe(self.ref.topology, self.ref.trajectory)
        outfile = self.tmp_file('write-uni-test')
        with self.ref.writer(outfile, uni.atoms.n_atoms) as w:
            for ts in uni.trajectory:
                w.write(uni)
        self._check_copy(outfile)

    def test_write_selection(self):
        uni = mda.Universe(self.ref.topology, self.ref.trajectory)
        sel_str = 'resid 1'
        sel = uni.select_atoms(sel_str)
        outfile = self.tmp_file('write-selection-test')

        with self.ref.writer(outfile, sel.n_atoms) as W:
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

    @raises(TypeError)
    def test_write_none(self):
        outfile = self.tmp_file('write-none')
        with self.ref.writer(outfile, 42) as w:
            w.write(None)

    def test_no_container(self):
        with tempdir.in_tempdir():
            if self.ref.container_format:
                self.ref.writer('foo')
            else:
                assert_raises(TypeError, self.ref.writer, 'foo')

    def test_write_not_changing_ts(self):
        outfile = self.tmp_file('write-not-changing-ts')
        ts = self.reader.ts.copy()
        copy_ts = ts.copy()
        with self.ref.writer(outfile, n_atoms=5) as W:
            W.write(ts)
            assert_timestep_almost_equal(copy_ts, ts)


class BaseTimestepTest(object):
    """Test all the base functionality of a Timestep

    All Timesteps must pass these tests!

    These test the Timestep independent of the Reader which it
    comes into contact with.  Failures here are the Timesteps fault.
    """
    # define the class made in test
    Timestep = mda.coordinates.base.Timestep
    name = "base"  # for error messages only
    size = 10  # size of arrays, 10 is enough to allow slicing etc
    # each coord is unique
    refpos = np.arange(size * 3, dtype=np.float32).reshape(size, 3) * 1.234
    refvel = np.arange(size * 3, dtype=np.float32).reshape(size, 3) * 2.345
    reffor = np.arange(size * 3, dtype=np.float32).reshape(size, 3) * 3.456

    has_box = True
    set_box = True  # whether you can set dimensions info.
    # If you can set box, what the underlying unitcell should be
    # if dimensions are:
    newbox = np.array([10., 11., 12., 90., 90., 90.])
    unitcell = np.array([10., 11., 12., 90., 90., 90.])
    ref_volume = 1320.  # what the volume is after setting newbox
    uni_args = None

    def setUp(self):
        self.ts = self.Timestep(self.size)
        self.ts.frame += 1
        self.ts.positions = self.refpos

    def tearDown(self):
        del self.ts

    def test_getitem(self):
        assert_equal(self.ts[1], self.refpos[1])

    def test_getitem_neg(self):
        assert_equal(self.ts[-1], self.refpos[-1])

    def test_getitem_neg_IE(self):
        assert_raises(IndexError, self.ts.__getitem__, -(self.size + 1))

    def test_getitem_pos_IE(self):
        assert_raises(IndexError, self.ts.__getitem__, (self.size + 1))

    def test_getitem_slice(self):
        assert_equal(len(self.ts[:2]), len(self.refpos[:2]))
        assert_allclose(self.ts[:2], self.refpos[:2])

    def test_getitem_slice2(self):
        assert_equal(len(self.ts[1::2]), len(self.refpos[1::2]))
        assert_allclose(self.ts[1::2], self.refpos[1::2])

    def test_getitem_ndarray(self):
        sel = np.array([0, 1, 4])
        assert_equal(len(self.ts[sel]), len(self.refpos[sel]))
        assert_allclose(self.ts[sel], self.refpos[sel])

    def test_getitem_TE(self):
        assert_raises(TypeError, self.ts.__getitem__, 'string')

    def test_len(self):
        assert_equal(len(self.ts), self.size)

    def test_iter(self):
        for a, b in zip(self.ts, self.refpos):
            assert_allclose(a, b)
        assert_equal(len(list(self.ts)), self.size)

    def test_repr(self):
        assert_equal(type(repr(self.ts)), str)

    # Dimensions has 2 possible cases
    # Timestep doesn't do dimensions,
    # should raise NotImplementedError for .dimension and .volume
    # Timestep does do them, should return values properly
    def test_dimensions(self):
        if self.has_box:
            assert_allclose(self.ts.dimensions, np.zeros(6, dtype=np.float32))
        else:
            assert_raises(NotImplementedError, getattr, self.ts, "dimensions")

    def test_dimensions_set_box(self):
        if self.set_box:
            self.ts.dimensions = self.newbox
            assert_allclose(self.ts._unitcell, self.unitcell)
            assert_allclose(self.ts.dimensions, self.newbox)
        else:
            pass

    def test_volume(self):
        if self.has_box and self.set_box:
            self.ts.dimensions = self.newbox
            assert_equal(self.ts.volume, self.ref_volume)
        elif self.has_box and not self.set_box:
            pass  # How to test volume of box when I don't set unitcell first?
        else:
            assert_raises(NotImplementedError, getattr, self.ts, "volume")

    def test_triclinic_vectors(self):
        assert_allclose(self.ts.triclinic_dimensions,
                        triclinic_vectors(self.ts.dimensions))

    def test_set_triclinic_vectors(self):
        ref_vec = triclinic_vectors(self.newbox)
        self.ts.triclinic_dimensions = ref_vec
        assert_equal(self.ts.dimensions, self.newbox)
        assert_allclose(self.ts._unitcell, self.unitcell)

    @attr('issue')
    def test_coordinate_getter_shortcuts(self):
        """testing that reading _x, _y, and _z works as expected
        # (Issue 224) (TestTimestep)"""
        assert_allclose(self.ts._x, self.ts._pos[:, 0])
        assert_allclose(self.ts._y, self.ts._pos[:, 1])
        assert_allclose(self.ts._z, self.ts._pos[:, 2])

    def test_coordinate_setter_shortcuts(self):
        # Check that _x _y and _z are read only
        for coordinate in ('_x', '_y', '_z'):
            random_positions = np.arange(self.size).astype(np.float32)
            assert_raises(AttributeError, setattr,
                          self.ts, coordinate, random_positions)

    # n_atoms should be a read only property
    # all Timesteps require this attribute
    def test_n_atoms(self):
        assert_equal(self.ts.n_atoms, self.ts._n_atoms)

    def test_n_atoms_readonly(self):
        assert_raises(AttributeError, self.ts.__setattr__, 'n_atoms', 20)

    def test_n_atoms_presence(self):
        assert_equal(hasattr(self.ts, '_n_atoms'), True)

    def test_unitcell_presence(self):
        assert_equal(hasattr(self.ts, '_unitcell'), True)

    def test_data_presence(self):
        assert_equal(hasattr(self.ts, 'data'), True)
        assert_equal(isinstance(self.ts.data, dict), True)

    def test_allocate_velocities(self):
        assert_equal(self.ts.has_velocities, False)
        assert_raises(NoDataError, getattr, self.ts, 'velocities')

        self.ts.has_velocities = True
        assert_equal(self.ts.has_velocities, True)
        assert_equal(self.ts.velocities.shape, (self.ts.n_atoms, 3))

    def test_allocate_forces(self):
        assert_equal(self.ts.has_forces, False)
        assert_raises(NoDataError, getattr, self.ts, 'forces')

        self.ts.has_forces = True
        assert_equal(self.ts.has_forces, True)
        assert_equal(self.ts.forces.shape, (self.ts.n_atoms, 3))

    def test_velocities_remove(self):
        ts = self.Timestep(10, velocities=True)
        ts.frame += 1
        assert_equal(ts.has_velocities, True)

        ts.has_velocities = False
        assert_equal(ts.has_velocities, False)
        assert_raises(NoDataError, getattr, ts, 'velocities')

    def test_forces_remove(self):
        ts = self.Timestep(10, forces=True)
        ts.frame += 1
        assert_equal(ts.has_forces, True)

        ts.has_forces = False
        assert_equal(ts.has_forces, False)
        assert_raises(NoDataError, getattr, ts, 'forces')

    def _empty_ts(self):
        assert_raises(ValueError, self.Timestep.from_coordinates,
                      None, None, None)

    def _from_coords(self, p, v, f):
        posi = self.refpos if p else None
        velo = self.refvel if v else None
        forc = self.reffor if f else None

        ts = self.Timestep.from_coordinates(posi, velo, forc)

        return ts

    def _check_from_coordinates(self, p, v, f):
        ts = self._from_coords(p, v, f)

        if p:
            assert_array_almost_equal(ts.positions, self.refpos)
        else:
            assert_raises(NoDataError, getattr, ts, 'positions')
        if v:
            assert_array_almost_equal(ts.velocities, self.refvel)
        else:
            assert_raises(NoDataError, getattr, ts, 'velocities')
        if f:
            assert_array_almost_equal(ts.forces, self.reffor)
        else:
            assert_raises(NoDataError, getattr, ts, 'forces')

    def test_from_coordinates(self):
        # Check all combinations of creating a Timestep from data
        # 8 possibilites of with and without 3 data categories
        for p, v, f in itertools.product([True, False], repeat=3):
            if not any([p, v, f]):
                yield self._empty_ts
            else:
                yield self._check_from_coordinates, p, v, f

    def test_from_coordinates_mismatch(self):
        velo = self.refvel[:2]

        assert_raises(ValueError, self.Timestep.from_coordinates,
                      self.refpos, velo)

    def test_from_coordinates_nodata(self):
        assert_raises(ValueError, self.Timestep.from_coordinates)

    def _check_from_timestep(self, p, v, f):
        ts = self._from_coords(p, v, f)
        ts2 = self.Timestep.from_timestep(ts)

        assert_timestep_almost_equal(ts, ts2)

    def test_from_timestep(self):
        for p, v, f in itertools.product([True, False], repeat=3):
            if not any([p, v, f]):
                continue
            yield self._check_from_timestep, p, v, f

    # Time related tests
    def test_supply_dt(self):
        # Check that this gets stored in data properly
        ts = self.Timestep(20, dt=0.04)

        assert_equal(ts.data['dt'], 0.04)
        assert_equal(ts.dt, 0.04)

    def test_redefine_dt(self):
        ts = self.Timestep(20, dt=0.04)
        assert_equal(ts.data['dt'], 0.04)
        assert_equal(ts.dt, 0.04)
        ts.dt = refdt = 0.46
        assert_equal(ts.data['dt'], refdt)
        assert_equal(ts.dt, refdt)

    def test_delete_dt(self):
        ts = self.Timestep(20, dt=0.04)
        assert_equal(ts.data['dt'], 0.04)
        assert_equal(ts.dt, 0.04)
        del ts.dt
        assert_equal('dt' in ts.data, False)
        assert_equal(ts.dt, 1.0)  # default value

    def test_supply_time_offset(self):
        ts = self.Timestep(20, time_offset=100.0)

        assert_equal(ts.data['time_offset'], 100.0)

    def test_time(self):
        ts = self.Timestep(20)
        ts.frame = 0
        ts.time = reftime = 1234.0

        assert_equal(ts.time, reftime)

    def test_delete_time(self):
        ts = self.Timestep(20)
        ts.frame = 0
        ts.time = reftime = 1234.0

        assert_equal(ts.time, reftime)
        del ts.time
        assert_equal(ts.time, 0.0)  # default to 1.0 (dt) * 0 (frame)

    def test_time_with_offset(self):
        ref_offset = 2345.0
        ts = self.Timestep(20, time_offset=ref_offset)
        ts.frame = 0
        ts.time = reftime = 1234.0

        assert_equal(ts.time, reftime + ref_offset)

    def test_dt(self):
        ref_dt = 45.0
        ts = self.Timestep(20, dt=ref_dt)

        for i in range(10):
            ts.frame = i
            assert_equal(ts.time, i * ref_dt)

    def test_change_dt(self):
        ref_dt = 45.0
        ts = self.Timestep(20, dt=ref_dt)

        for i in range(10):
            ts.frame = i
            assert_equal(ts.time, i * ref_dt)

        ts.dt = ref_dt = 77.0

        for i in range(10):
            ts.frame = i
            assert_equal(ts.time, i * ref_dt)

    def test_dt_with_offset(self):
        ref_dt = 45.0
        ref_offset = 2345.0
        ts = self.Timestep(20, dt=ref_dt, time_offset=ref_offset)

        for i in range(10):
            ts.frame = i
            assert_equal(ts.time, i * ref_dt + ref_offset)

    def test_time_overrides_dt_with_offset(self):
        ref_dt = 45.0
        ref_offset = 2345.0
        ts = self.Timestep(20, dt=ref_dt, time_offset=ref_offset)

        ts.frame = 0
        ts.time = reftime = 456.7

        assert_equal(ts.time, reftime + ref_offset)

    def _check_copy(self, name, ref_ts):
        """Check basic copy"""
        ts2 = ref_ts.copy()

        err_msg = ("Timestep copy failed for format {form}"
                   " on attribute {att}")

        # eq method checks:
        # - frame
        # - n_atoms
        # - positions, vels and forces
        assert_(ref_ts == ts2)

        assert_array_almost_equal(ref_ts.dimensions, ts2.dimensions,
                                  decimal=4)

        # Check things not covered by eq
        for d in ref_ts.data:
            assert_(d in ts2.data)
            if isinstance(ref_ts.data[d], np.ndarray):
                assert_array_almost_equal(
                    ref_ts.data[d], ts2.data[d])
            else:
                assert_(ref_ts.data[d] == ts2.data[d])

    def _check_independent(self, name, ts):
        """Check that copies made are independent"""
        ts2 = ts.copy()

        if ts.has_positions:
            self._check_array(ts.positions, ts2.positions)
        if ts.has_velocities:
            self._check_array(ts.velocities, ts2.velocities)
        if ts.has_forces:
            self._check_array(ts.forces, ts2.forces)
        self._check_array(ts.dimensions, ts2.dimensions)

    def _check_array(self, arr1, arr2):
        """Check modifying one array doesn't change other"""
        ref = arr1.copy()
        arr2 += 1.0
        assert_array_almost_equal(ref, arr1)

    def _check_copy_slice_indices(self, name, ts):
        sl = slice(0, len(ts), 3)
        ts2 = ts.copy_slice(sl)
        self._check_slice(ts, ts2, sl)

    def _check_copy_slice_slice(self, name, ts):
        sl = [0, 1, 3, 5, 6, 7]
        ts2 = ts.copy_slice(sl)
        self._check_slice(ts, ts2, sl)

    def _check_npint_slice(self, name, ts):
        for integers in [np.byte, np.short, np.intc, np.int_, np.longlong,
                         np.intp, np.int8, np.int16, np.int32, np.int64, 
                         np.ubyte, np.ushort, np.uintc, np.ulonglong,
                         np.uintp, np.uint8, np.uint16, np.uint32, np.uint64]:
            sl = slice(1, 2, 1)
            ts2 = ts.copy_slice(slice(integers(1), integers(2), integers(1)))
            self._check_slice(ts, ts2, sl)
        
    def _check_slice(self, ts1, ts2, sl):
        if ts1.has_positions:
            assert_array_almost_equal(ts1.positions[sl], ts2.positions)
        if ts1.has_velocities:
            assert_array_almost_equal(ts1.velocities[sl], ts2.velocities)
        if ts1.has_forces:
            assert_array_almost_equal(ts1.forces[sl], ts2.forces)

    def test_copy(self):
        if self.uni_args is None:
            return
        u = mda.Universe(*self.uni_args)
        ts = u.trajectory.ts
        
        yield self._check_copy, self.name, ts
        yield self._check_independent, self.name, ts
        yield self._check_copy_slice_indices, self.name, ts
        yield self._check_copy_slice_slice, self.name, ts
        yield self._check_npint_slice, self.name, ts

    def test_copy_slice(self):
        for p, v, f in itertools.product([True, False], repeat=3):
            if not any([p, v, f]):
                continue
            ts = self._from_coords(p, v, f)
            yield self._check_copy, self.name, ts
            yield self._check_independent, self.name, ts
            yield self._check_copy_slice_indices, self.name, ts
            yield self._check_copy_slice_slice, self.name, ts

    def _check_bad_slice(self, p, v, f):
        ts = self._from_coords(p, v, f)
        sl = ['this', 'is', 'silly']
        assert_raises(TypeError, ts.copy_slice, sl)

    def test_bad_copy_slice(self):
        for p, v, f in itertools.product([True, False], repeat=3):
            if not any([p, v, f]):
                continue
            yield self._check_bad_slice, p, v, f

    def _get_pos(self):
        # Get generic reference positions
        return np.arange(30).reshape(10, 3) * 1.234

    def _check_ts_equal(self, a, b, err_msg):
        assert_(a == b, err_msg)
        assert_(b == a, err_msg)

    def test_check_equal(self):
        for p, v, f in itertools.product([True, False], repeat=3):
            if not any([p, v, f]):
                continue

            ts1 = self.Timestep(self.size,
                                positions=p,
                                velocities=v,
                                forces=f)
            ts2 = self.Timestep(self.size,
                                positions=p,
                                velocities=v,
                                forces=f)
            if p:
                ts1.positions = self.refpos.copy()
                ts2.positions = self.refpos.copy()
            if v:
                ts1.velocities = self.refvel.copy()
                ts2.velocities = self.refvel.copy()
            if f:
                ts1.forces = self.reffor.copy()
                ts2.forces = self.reffor.copy()

            yield (self._check_ts_equal, ts1, ts2,
                   'Failed on {0}'.format(self.name))

    def test_wrong_class_equality(self):
        ts1 = self.Timestep(self.size)
        ts1.positions = self._get_pos()

        b = tuple([0, 1, 2, 3])

        assert_(ts1 != b)
        assert_(b != ts1)

    def test_wrong_frame_equality(self):
        ts1 = self.Timestep(self.size)
        ts1.positions = self._get_pos()

        ts2 = self.Timestep(self.size)
        ts2.positions = self._get_pos()
        ts2.frame = 987

        assert_(ts1 != ts2)
        assert_(ts2 != ts1)

    def test_wrong_n_atoms_equality(self):
        ts1 = self.Timestep(self.size)
        ts1.positions = self._get_pos()

        ts3 = self.Timestep(self.size * 2)

        assert_(ts1 != ts3)
        assert_(ts3 != ts1)

    def test_wrong_pos_equality(self):
        ts1 = self.Timestep(self.size)
        ts1.positions = self._get_pos()

        ts2 = self.Timestep(self.size)
        ts2.positions = self._get_pos() + 1.0

        assert_(ts1 != ts2)
        assert_(ts2 != ts1)

    def test_check_vels_equality(self):
        ts1 = self.Timestep(self.size, velocities=True)
        ts2 = self.Timestep(self.size, velocities=True)

        ts1.velocities = self._get_pos()
        ts2.velocities = self._get_pos()

        assert_(ts1 == ts2)
        assert_(ts2 == ts1)

    def test_check_mismatched_vels_equality(self):
        ts1 = self.Timestep(self.size, velocities=True)
        ts2 = self.Timestep(self.size, velocities=False)

        ts1.velocities = self._get_pos()

        assert_(ts1 != ts2)
        assert_(ts2 != ts1)

    def test_check_wrong_vels_equality(self):
        ts1 = self.Timestep(self.size, velocities=True)
        ts2 = self.Timestep(self.size, velocities=True)

        ts1.velocities = self._get_pos()
        ts2.velocities = self._get_pos() + 1.0

        assert_(ts1 != ts2)
        assert_(ts2 != ts1)

    def test_check_forces_equality(self):
        ts1 = self.Timestep(self.size, forces=True)
        ts2 = self.Timestep(self.size, forces=True)

        ts1.forces = self._get_pos()
        ts2.forces = self._get_pos()

        assert_(ts1 == ts2)
        assert_(ts2 == ts1)

    def test_check_mismatched_forces_equality(self):
        ts1 = self.Timestep(self.size, forces=True)
        ts2 = self.Timestep(self.size, forces=False)

        ts1.forces = self._get_pos()

        assert_(ts1 != ts2)
        assert_(ts2 != ts1)

    def test_check_wrong_forces_equality(self):
        ts1 = self.Timestep(self.size, forces=True)
        ts2 = self.Timestep(self.size, forces=True)

        ts1.forces = self._get_pos()
        ts2.forces = self._get_pos() + 1.0

        assert_(ts1 != ts2)
        assert_(ts2 != ts1)


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
