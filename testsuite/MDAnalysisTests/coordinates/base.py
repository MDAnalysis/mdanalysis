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
import itertools
import numpy as np
import pytest
from six.moves import zip, range
from unittest import TestCase
from numpy.testing import (assert_equal, assert_almost_equal,
                           assert_array_almost_equal, assert_allclose)

import MDAnalysis as mda
from MDAnalysis.coordinates.base import Timestep
from MDAnalysis import NoDataError


from MDAnalysisTests.coordinates.reference import RefAdKSmall
from MDAnalysisTests.datafiles import AUX_XVG_HIGHF, AUX_XVG_LOWF
from MDAnalysisTests import tempdir, make_Universe


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
    @staticmethod
    @pytest.fixture()
    def reader(ref):
        reader = ref.reader(ref.trajectory)
        reader.add_auxiliary('lowf', ref.aux_lowf, dt=ref.aux_lowf_dt, initial_time=0, time_selector=None)
        reader.add_auxiliary('highf', ref.aux_highf, dt=ref.aux_highf_dt, initial_time=0, time_selector=None)
        return reader

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

    def test_get_writer_1(self, ref, reader):
        with tempdir.in_tempdir():
            outfile = 'test-writer' + ref.ext
            with reader.Writer(outfile) as W:
                assert_equal(isinstance(W, ref.writer), True)
                assert_equal(W.n_atoms, reader.n_atoms)

    def test_get_writer_2(self, ref, reader):
        with tempdir.in_tempdir():
            outfile = 'test-writer' + ref.ext
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
        assert_array_almost_equal(reader.ts.dimensions,
                                  ref.dimensions,
                                  decimal=ref.prec)

    def test_changing_dimensions(self, ref, reader):
        if ref.changing_dimensions:
            reader.rewind()
            assert_array_almost_equal(reader.ts.dimensions,
                                      ref.dimensions,
                                      decimal=ref.prec)
            reader[1]
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
    def tempdir():
        return tempdir.TempDir()

    def tmp_file(self, name, ref, tempdir):
        return tempdir.name + name + '.' + ref.ext

    def test_write_trajectory_timestep(self,ref, reader, tempdir):
        outfile = self.tmp_file('write-timestep-test', ref, tempdir)
        with ref.writer(outfile, reader.n_atoms) as W:
            for ts in reader:
                W.write(ts)
        self._check_copy(outfile, ref, reader)

    def test_write_different_box(self, ref, reader, tempdir):
        if ref.changing_dimensions:
            outfile = self.tmp_file('write-dimensions-test', ref, tempdir)
            with ref.writer(outfile, reader.n_atoms) as W:
                for ts in reader:
                    ts.dimensions[:3] += 1
                    W.write(ts)

            written = ref.reader(outfile)

            for ts_ref, ts_w in zip(reader, written):
                ts_ref.dimensions[:3] += 1
                assert_array_almost_equal(ts_ref.dimensions,
                                          ts_w.dimensions,
                                          decimal=ref.prec)

    def test_write_trajectory_atomgroup(self, ref,reader, tempdir):
        uni = mda.Universe(ref.topology, ref.trajectory)
        outfile = self.tmp_file('write-atoms-test', ref, tempdir)
        with ref.writer(outfile, uni.atoms.n_atoms) as w:
            for ts in uni.trajectory:
                w.write(uni.atoms)
        self._check_copy(outfile, ref, reader)

    def test_write_trajectory_universe(self, ref, reader, tempdir):
        uni = mda.Universe(ref.topology, ref.trajectory)
        outfile = self.tmp_file('write-uni-test', ref, tempdir)
        with ref.writer(outfile, uni.atoms.n_atoms) as w:
            for ts in uni.trajectory:
                w.write(uni)
        self._check_copy(outfile, ref, reader)

    def test_write_selection(self, ref, reader, u_no_resnames, u_no_resids, u_no_names, tempdir):
        uni = mda.Universe(ref.topology, ref.trajectory)
        sel_str = 'resid 1'
        sel = uni.select_atoms(sel_str)
        outfile = self.tmp_file('write-selection-test', ref, tempdir)

        with ref.writer(outfile, sel.n_atoms) as W:
            for ts in uni.trajectory:
                W.write(sel.atoms)

        copy = ref.reader(outfile)
        for orig_ts, copy_ts in zip(uni.trajectory, copy):
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

    def test_write_none(self, ref, tempdir):
        outfile = self.tmp_file('write-none', ref, tempdir)
        with pytest.raises(TypeError):
            with ref.writer(outfile, 42) as w:
                w.write(None)

    def test_no_container(self, ref):
        with tempdir.in_tempdir():
            if ref.container_format:
                ref.writer('foo')
            else:
                with pytest.raises(TypeError):
                    ref.writer('foo')

    def test_write_not_changing_ts(self, ref, reader, tempdir):
        outfile = self.tmp_file('write-not-changing-ts', ref, tempdir)
        ts = reader.ts.copy()
        copy_ts = ts.copy()
        with ref.writer(outfile, n_atoms=5) as W:
            W.write(ts)
            assert_timestep_almost_equal(copy_ts, ts)


def assert_timestep_equal(A, B, msg=''):
    """ assert that two timesteps are exactly equal and commutative
    """
    assert A == B, msg
    assert B == A, msg


def assert_timestep_almost_equal(A, B, decimal=6, verbose=True):
    if not isinstance(A, Timestep):
        raise AssertionError('A is not of type Timestep')
    if not isinstance(B, Timestep):
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
