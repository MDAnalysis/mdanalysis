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
import sys
import warnings
import numpy as np
import os

import pytest

from numpy.testing import (assert_equal, assert_almost_equal)

import MDAnalysis as mda
from MDAnalysis.transformations import translate
from MDAnalysisTests.datafiles import (PDB, PSF, CRD, DCD,
                                       GRO, XTC, TRR, PDB_small,
                                       LAMMPS_chain, LAMMPSDUMP_chain1,
                                       LAMMPSDUMP_chain2,)
from MDAnalysisTests.util import no_warning


class TestChainReader(object):
    prec = 3

    @pytest.fixture()
    def universe(self):
        return mda.Universe(PSF,
                            [DCD, CRD, DCD, CRD, DCD, CRD, CRD])

    @pytest.fixture()
    def transformed(ref):
        return mda.Universe(PSF,
                            [DCD, CRD, DCD, CRD, DCD, CRD, CRD],
                            transformations=[translate([10,10,10])])
                            
    def test_regular_repr(self):
        u = mda.Universe(PSF, [DCD, CRD, DCD])
        assert_equal("<ChainReader containing adk_dims.dcd, adk_open.crd, adk_dims.dcd with 197 frames of 3341 atoms>", u.trajectory.__repr__())
        
                                
    def test_truncated_repr(self, universe):
        assert_equal("<ChainReader containing adk_dims.dcd and 6 more with 298 frames of 3341 atoms>", universe.trajectory.__repr__())

    def test_next_trajectory(self, universe):
        universe.trajectory.rewind()
        universe.trajectory.next()
        assert_equal(universe.trajectory.ts.frame, 1, "loading frame 2")

    def test_n_atoms(self, universe):
        assert_equal(universe.trajectory.n_atoms, 3341,
                     "wrong number of atoms")

    def test_n_frames(self, universe):
        assert_equal(universe.trajectory.n_frames, 3 * 98 + 4,
                     "wrong number of frames in chained dcd")

    def test_iteration(self, universe):
        for ts in universe.trajectory:
            pass  # just forward to last frame
        assert_equal(universe.trajectory.n_frames - 1, ts.frame,
                     "iteration yielded wrong number of frames ({0:d}), "
                     "should be {1:d}".format(ts.frame,
                                              universe.trajectory.n_frames))

    def test_jump_lastframe_trajectory(self, universe):
        universe.trajectory[-1]
        assert_equal(universe.trajectory.ts.frame + 1,
                     universe.trajectory.n_frames,
                     "indexing last frame with trajectory[-1]")

    def test_slice_trajectory(self, universe):
        frames = [ts.frame for ts in universe.trajectory[5:17:3]]
        assert_equal(frames, [5, 8, 11, 14], "slicing dcd [5:17:3]")

    def test_full_slice(self, universe):
        trj_iter = universe.trajectory[:]
        frames = [ts.frame for ts in trj_iter]
        assert_equal(frames, np.arange(universe.trajectory.n_frames))

    def test_frame_numbering(self, universe):
        universe.trajectory[98]  # index is 0-based and frames are 0-based
        assert_equal(universe.trajectory.frame, 98, "wrong frame number")

    def test_next_after_frame_numbering(self, universe):
        universe.trajectory[98]  # index is 0-based and frames are 0-based
        universe.trajectory.next()
        assert_equal(universe.trajectory.frame, 99, "wrong frame number")

    def test_frame(self, universe):
        universe.trajectory[0]
        coord0 = universe.atoms.positions.copy()
        # forward to frame where we repeat original dcd again:
        # dcd:0..97 crd:98 dcd:99..196
        universe.trajectory[99]
        assert_equal(universe.atoms.positions, coord0,
                     "coordinates at frame 1 and 100 should be the same!")

    def test_time(self, universe):
        universe.trajectory[30]  # index and frames 0-based
        assert_almost_equal(
            universe.trajectory.time, 30.0, 5, err_msg="Wrong time of frame")

    def test_write_dcd(self, universe, tmpdir):
        """test that ChainReader written dcd (containing crds) is correct
        (Issue 81)"""
        outfile = str(tmpdir) + "chain-reader.dcd"
        with mda.Writer(outfile, universe.atoms.n_atoms) as W:
            for ts in universe.trajectory:
                W.write(universe)
        universe.trajectory.rewind()
        u = mda.Universe(PSF, outfile)
        for (ts_orig, ts_new) in zip(universe.trajectory, u.trajectory):
            assert_almost_equal(
                ts_orig._pos,
                ts_new._pos,
                self.prec,
                err_msg="Coordinates disagree at frame {0:d}".format(
                    ts_orig.frame))       
    
    def test_transform_iteration(self, universe, transformed):
        vector = np.float32([10,10,10])
        # # Are the transformations applied and
        # are the coordinates "overtransformed"?
        # iterate once:
        for ts in transformed.trajectory:
            frame = ts.frame
            ref = universe.trajectory[frame].positions + vector
            assert_almost_equal(ts.positions, ref, decimal = 6)
        # iterate again:
        for ts in transformed.trajectory:
            frame = ts.frame
            ref = universe.trajectory[frame].positions + vector
            assert_almost_equal(ts.positions, ref, decimal = 6)
    
    def test_transform_slice(self, universe, transformed):
        vector = np.float32([10,10,10])
        # what happens when we slice the trajectory?
        for ts in transformed.trajectory[5:17:3]:
            frame = ts.frame
            ref = universe.trajectory[frame].positions + vector
            assert_almost_equal(ts.positions, ref, decimal = 6)
    
    def test_transform_switch(self, universe, transformed):
        vector = np.float32([10,10,10])
        # grab a frame:
        ref = universe.trajectory[2].positions + vector
        assert_almost_equal(transformed.trajectory[2].positions, ref, decimal = 6)
        # now switch to another frame
        newref = universe.trajectory[10].positions + vector
        assert_almost_equal(transformed.trajectory[10].positions, newref, decimal = 6)
        # what happens when we comeback to the previous frame?
        assert_almost_equal(transformed.trajectory[2].positions, ref, decimal = 6)
    
    def test_transfrom_rewind(self, universe, transformed):
        vector = np.float32([10,10,10])
        ref = universe.trajectory[0].positions + vector
        transformed.trajectory.rewind()
        assert_almost_equal(transformed.trajectory.ts.positions, ref, decimal = 6)


class TestChainReaderCommonDt(object):
    common_dt = 100.0
    prec = 3

    @pytest.fixture()
    def trajectory(self):
        universe = mda.Universe(
            PSF, [DCD, CRD, DCD, CRD, DCD, CRD, CRD], dt=self.common_dt)
        return universe.trajectory

    def test_time(self, trajectory):
        # We test this for the beginning, middle and end of the trajectory.
        for frame_n in (0, trajectory.n_frames // 2, -1):
            trajectory[frame_n]
            assert_almost_equal(
                trajectory.time,
                trajectory.frame * self.common_dt,
                5,
                err_msg="Wrong time for frame {0:d}".format(frame_n))


class TestChainReaderFormats(object):
    """Test of ChainReader with explicit formats (Issue 76)."""

    def test_set_all_format_tuples(self):
        universe = mda.Universe(GRO, [(PDB, 'pdb'), (XTC, 'xtc'), (TRR,
                                                                   'trr')])
        assert universe.trajectory.n_frames == 21
        assert_equal(universe.trajectory.filenames, [PDB, XTC, TRR])

    def test_set_all_format_lammps(self):
        universe = mda.Universe(
            LAMMPS_chain, [LAMMPSDUMP_chain1, LAMMPSDUMP_chain2], 
            format="LAMMPSDUMP", continuous=True)
        assert universe.trajectory.n_frames == 11

        # Test whether the amount of time between frames is consistently = 1
        time_values = np.array([ts.time for ts in universe.trajectory])
        dt_array = time_values[1:]-time_values[:-1]
        assert np.unique(dt_array) == 1
        assert_equal(time_values, np.arange(11))

    def test_set_format_tuples_and_format(self):
        universe = mda.Universe(GRO, [(PDB, 'pdb'), GRO, GRO, (XTC, 'xtc'), 
                                      (TRR, 'trr')], format='gro')
        assert universe.trajectory.n_frames == 23
        assert_equal(universe.trajectory.filenames, [PDB, GRO, GRO, XTC, TRR])
        
        with pytest.raises(TypeError) as errinfo:
            mda.Universe(GRO, [(PDB, 'pdb'), GRO, GRO, (XTC, 'xtc'), 
                                      (TRR, 'trr')], format='pdb')
        assert 'Unable to read' in str(errinfo.value)


    def test_set_one_format_tuple(self):
        universe = mda.Universe(PSF, [(PDB_small, 'pdb'), DCD])
        assert universe.trajectory.n_frames == 99

    def test_set_all_formats(self):
        with pytest.raises(TypeError) as errinfo:
            mda.Universe(PDB, [PDB, GRO], format='gro')
        assert 'Unable to read' in str(errinfo.value)

        universe = mda.Universe(GRO, [PDB, PDB, PDB], format='pdb')
        assert_equal(universe.trajectory.filenames, [PDB, PDB, PDB])


def build_trajectories(folder, sequences, fmt='xtc'):
    """
    A scenario is given as a series of time sequences. The result is
    returned as a list of time and origin of the frame. Each element of that
    result list correspond to a frame of the conatenated trajectory, The
    first element of each tuple is the time it corresponds to, the second
    element is the index of the input the frame comes from.

    See original gist and PR 1728 for comparison with gromacs trjcat.

    https://gist.github.com/jbarnoud/cacd0957d3df01d1577f640b20e86039
    """
    template = 'trjcat_test_{{}}.{}'.format(fmt)
    template = os.path.join(folder, template)

    # Use an empty universe to have a topology
    utop = mda.Universe.empty(1, trajectory=True)

    # Create synthetic trajectories. The times come from the user input,
    # the coordinates indicate the index of the input.
    fnames = []
    for index, subseq in enumerate(sequences):
        coords = np.zeros((len(subseq), 1, 3), dtype=np.float32) + index
        u = mda.Universe(utop._topology, coords)
        out_traj = mda.Writer(template.format(index), n_atoms=len(u.atoms))
        fnames.append(out_traj.filename)
        with out_traj:
            for ts, time in zip(u.trajectory, subseq):
                # The time step needs a box to avoid a warning
                ts.dimensions = [10, 10, 10, 90, 90, 90]
                # The time is set from the user input
                ts.time = time
                out_traj.write(u)
    return utop, fnames


class TestChainReaderContinuous(object):
    class SequenceInfo(object):
        def __init__(self, seq, n_frames, order):
            self.seq = seq
            self.n_frames = n_frames
            self.order = order

    @pytest.mark.parametrize('fmt', ('xtc', 'trr'))
    @pytest.mark.parametrize('seq_info',[
        SequenceInfo(seq=([0, 1, 2, 3], [2, 3, 4, 5], [4, 5, 6, 7]),
                     n_frames=8,
                     order=[0, 0, 1, 1, 2, 2, 2, 2]),
        SequenceInfo(seq=([0, 1, 2, 4], [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
                     n_frames=10,
                     order=np.ones(10)),
        SequenceInfo(seq=([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], [0, 1, 2, 3]),
                     n_frames=4,
                     order=np.ones(4)),
        SequenceInfo(seq=([5, 6, 7, 8, 9], [2, 3, 4, 5, 6], [0, 1, 2, 3]),
                     n_frames=10,
                     order=[2, 2, 1, 1, 1, 0, 0, 0, 0, 0]),
        SequenceInfo(seq=([0, 1, 2],) * 3,
                     n_frames=3,
                     order=[2, 2, 2]),
        SequenceInfo(seq=([0, 1, 2, 3,], [3, 4], [4, 5, 6, 7]),
                     n_frames=8,
                     order=[0, 0, 0, 1, 2, 2, 2, 2]),
        SequenceInfo(seq=([5, 6, 7, 8, 9], [2, 3, 4, 5, 6], [0, 1, 2, 3]),
                     n_frames=10,
                     order=[2, 2, 1, 1, 1, 0, 0, 0, 0, 0]),
        SequenceInfo(seq=[list(range(0, 6)), list(range(2, 5)),
                          list(range(2, 5)), list(range(2, 5)),
                          list(range(3, 8))],
                     n_frames=8,
                     order=[0, 0, 3, 4, 4, 4, 4, 4]),
    ])
    def test_order(self, seq_info, tmpdir, fmt):
        folder = str(tmpdir)
        utop, fnames = build_trajectories(folder, sequences=seq_info.seq, fmt=fmt)
        u = mda.Universe(utop._topology, fnames, continuous=True)
        assert u.trajectory.n_frames == seq_info.n_frames
        for i, ts in enumerate(u.trajectory):
            assert_almost_equal(i, ts.time, decimal=4)
            # check we have used the right trajectory
            assert seq_info.order[i] == int(ts.positions[0, 0])

    def test_start_frames(self, tmpdir):
        folder = str(tmpdir)
        sequences = ([0, 1, 2, 3], [2, 3, 4, 5], [4, 5, 6, 7])
        utop, fnames = build_trajectories(folder, sequences=sequences,)
        u = mda.Universe(utop._topology, fnames, continuous=True)
        assert_equal(u.trajectory._start_frames, [0, 2, 4])

    def test_missing(self, tmpdir):
        folder = str(tmpdir)
        sequences = ([0, 1, 2, 3], [5, 6, 7, 8, 9])
        utop, fnames = build_trajectories(folder, sequences=sequences,)
        u = mda.Universe(utop._topology, fnames, continuous=True)
        assert u.trajectory.n_frames == 9

    def test_warning(self, tmpdir):
        folder = str(tmpdir)
        # this sequence *should* trigger a warning
        sequences = ([0, 1, 2, 3], [5, 6, 7])
        utop, fnames = build_trajectories(folder, sequences=sequences,)
        with pytest.warns(UserWarning):
            mda.Universe(utop._topology, fnames, continuous=True)

    def test_interleaving_error(self, tmpdir):
        folder = str(tmpdir)
        # interleaving is not supported by chainreader
        sequences = ([0, 2, 4, 6], [1, 3, 5, 7])
        utop, fnames = build_trajectories(folder, sequences=sequences,)
        with pytest.raises(RuntimeError):
            mda.Universe(utop._topology, fnames, continuous=True)

    def test_easy_trigger_warning(self, tmpdir):
        folder = str(tmpdir)
        # this sequence shouldn't trigger a warning
        sequences = ([0, 1, 2, 3], [2, 3, 4, 5], [4, 5, 6, 7])
        utop, fnames = build_trajectories(folder, sequences=sequences,)
        with no_warning(UserWarning):
            with warnings.catch_warnings():
                # for windows Python 3.10 ignore:
                # ImportWarning('_SixMetaPathImporter.find_spec() not found
                # TODO: remove when we no longer have a dependency
                # that still imports six
                if sys.version_info >= (3, 10):
                    warnings.filterwarnings(
                            action='ignore',
                            category=ImportWarning)
                mda.Universe(utop._topology, fnames, continuous=True)

    def test_single_frames(self, tmpdir):
        folder = str(tmpdir)
        sequences = ([0, 1, 2, 3], [5, ])
        utop, fnames = build_trajectories(folder, sequences=sequences,)
        with pytest.raises(RuntimeError):
            mda.Universe(utop._topology, fnames, continuous=True)

    def test_mixed_filetypes(self):
        with pytest.raises(ValueError):
            mda.Universe(PDB, [XTC, TRR], continuous=True)

    def test_unsupported_filetypes(self):
        with pytest.raises(NotImplementedError):
            mda.Universe(PSF, [DCD, DCD], continuous=True)
        # see issue 2353. The PDB reader has multiple format endings. To ensure 
        # the not implemented error is thrown  we  do a check here. A more  
        # careful test in the future would be a dummy reader with multiple 
        # formats, just in case PDB will allow continuous reading in the future.
        with pytest.raises(ValueError):
            mda.Universe(PDB, [PDB, XTC], continuous=True)


@pytest.mark.parametrize('l, ref', ([((0, 3), (3, 3), (4, 7)), (0, 1, 2)],
                                    [((0, 9), (0, 4)), (0, 1)],
                                    [((0, 3), (2, 2), (3, 3), (2, 6), (5, 9)), (0, 1, 3, 2, 4)],
                                    [((0, 2), (4, 9), (0, 4), (7, 9)), (0, 2, 1, 3)],
                                    [((0, 5), (2, 4), (2, 4), (2, 4), (3, 7)), (0, 1, 2, 3, 4)]
                                    ))
def test_multilevel_arg_sort(l, ref):
    indices = mda.coordinates.chain.multi_level_argsort(l)
    assert_equal(indices, ref)


@pytest.mark.parametrize('l, ref', ([((0, 4), (3, 6), (6, 9)), (0, 1, 2)],
                                    [((0, 3), (3, 4), (4, 7)), (0, 1, 2)],
                                    [((0, 3), (3, 5), (4, 7)), (0, 1, 2)],
                                    [((0, 3), (0, 4), (3, 5), (4, 7)), (1, 2, 3)],
                                    [((0, 3), (0, 4)), (1,)],
                                    [((0, 3), (0, 3)), (1,)],
                                    [((1, 3), (0, 4)), (0,)],
                                    [((0, 3), ), (0, )],
                                    [((0, 3), (5, 9)), (0, 1)],
                                    [((0, 3), (0, 3), (5, 9)), (1, 2)],
                                    [((0, 3), (0, 3), (0, 3), (5, 9)), (2, 3)],
                                    [((0, 5), (2, 4), (2, 4), (2, 4), (3, 7)), (0, 3, 4)],
                                    [((0, 3), (2, 4), (2, 4), (2, 4), (3, 7)), (0, 3, 4)],
                                    [((0, 3), (2, 4), (4, 7), (4, 7), (4, 7), (6, 9)), (0, 1, 4, 5)],
                                    [((0, 6), (2, 5), (2, 5), (2, 5), (3, 8)), (0, 3, 4)],
                                    ))
def test_filter_times(l, ref):
    indices = mda.coordinates.chain.filter_times(l, dt=1)
    assert_equal(indices, ref)


def test_issue_4008():
    # test for transformation on a chainreader trajectory
    # in issue #4008 it was observed that the transformation did not work on
    # the frames which were "borders" of the chained trajectories
    def shift_by_group(ag):
        # this just shifts all atoms in the box by their center of mass,
        # i.e. afterwards the center of mass of all atoms should be [0,0,0]
        def wrapped(ts):
            shift = ag.center_of_mass()
            ts.positions[:] -= shift
            return ts

        return wrapped

    u = mda.Universe(GRO, [XTC, XTC])

    trafos = (shift_by_group(u.atoms),)
    u.trajectory.add_transformations(*trafos)
    com = np.zeros((len(u.trajectory), 3))

    for i, ts in enumerate(u.trajectory):
        com[i] = u.atoms.center_of_mass()

    # see issue for details, but in general the com shouldn't be far from zero
    # in the issue report values of ~25 were seen before the fix
    # with the fix all values should be <1e-6
    assert np.abs(com).max() < 1e-6


def test_issue_3349():
    # test for repr of a chainreader when the subreaders are memoryreaders
    cdx = np.random.random(3341 * 3 * 10)
    cdx = cdx.reshape(-1, 3341, 3).astype(np.float32) * 10
    u = mda.Universe(PSF, (cdx, DCD))
    u2 = mda.Universe(PSF, (cdx, cdx, cdx, DCD))
    u_expected_filenames = np.array([None, str(DCD)])
    u2_expected_filenames = np.array([None, None, None, str(DCD)])

    assert_equal("<ChainReader containing numpy.ndarray, adk_dims.dcd with " +
                 "108 frames of 3341 atoms>",
                 u.trajectory.__repr__())
    assert_equal(u_expected_filenames, u.trajectory.filenames)
    assert_equal("<ChainReader containing numpy.ndarray and 3 more with " +
                 "128 frames of 3341 atoms>",
                 u2.trajectory.__repr__())
    assert_equal(u2_expected_filenames, u2.trajectory.filenames)
