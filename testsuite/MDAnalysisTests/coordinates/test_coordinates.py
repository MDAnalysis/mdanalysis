# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
# Beckstein and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

import MDAnalysis as mda
import numpy as np
import os
from six.moves import zip

from nose.plugins.attrib import attr
from numpy.testing import (assert_allclose, assert_equal, assert_array_equal,
                           assert_almost_equal, dec)
from unittest import TestCase

from MDAnalysisTests.datafiles import (PDB, INPCRD, XYZ_five, PSF, CRD, DCD,
                                       GRO, XTC, TRR, PDB_small, PDB_closed)
from MDAnalysisTests.plugins.knownfailure import knownfailure
from MDAnalysisTests import parser_not_found, tempdir


class TestINPCRDReader(TestCase):
    """Test reading Amber restart coordinate files"""

    def _check_ts(self, ts):
        # Check a ts has the right values in
        ref_pos = np.array([[6.6528795, 6.6711416, -8.5963255],
                            [7.3133773, 5.8359736, -8.8294175],
                            [8.3254058, 6.2227613, -8.7098593],
                            [7.0833200, 5.5038197, -9.8417650],
                            [7.1129439, 4.6170351, -7.9729560]])
        for ref, val in zip(ref_pos, ts._pos):
            assert_allclose(ref, val)

    def test_reader(self):
        from MDAnalysis.coordinates.INPCRD import INPReader

        r = INPReader(INPCRD)

        assert_equal(r.n_atoms, 5)
        self._check_ts(r.ts)

    def test_universe_inpcrd(self):
        u = mda.Universe(XYZ_five, INPCRD)

        self._check_ts(u.trajectory.ts)

    def test_universe_restrt(self):
        u = mda.Universe(XYZ_five, INPCRD, format='RESTRT')
        self._check_ts(u.trajectory.ts)


class TestChainReader(TestCase):

    @dec.skipif(parser_not_found('DCD'),
                'DCD parset not available. Are you using python 3?')
    def setUp(self):
        self.universe = mda.Universe(PSF,
                                     [DCD, CRD, DCD, CRD, DCD, CRD, CRD])
        self.trajectory = self.universe.trajectory
        self.prec = 3
        # dummy output DCD file
        self.tmpdir = tempdir.TempDir()
        self.outfile = self.tmpdir.name + '/chain-reader.dcd'

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except OSError:
            pass
        del self.universe
        del self.tmpdir

    def test_next_trajectory(self):
        self.trajectory.rewind()
        self.trajectory.next()
        assert_equal(self.trajectory.ts.frame, 1, "loading frame 2")

    def test_n_atoms(self):
        assert_equal(self.universe.trajectory.n_atoms, 3341,
                     "wrong number of atoms")

    def test_n_frames(self):
        assert_equal(self.universe.trajectory.n_frames, 3 * 98 + 4,
                     "wrong number of frames in chained dcd")

    def test_iteration(self):
        for ts in self.trajectory:
            pass  # just forward to last frame
        assert_equal(
            self.trajectory.n_frames - 1, ts.frame,
            "iteration yielded wrong number of frames ({0:d}), should be {1:d}".format(ts.frame, self.trajectory.n_frames))

    def test_jump_lastframe_trajectory(self):
        self.trajectory[-1]
        assert_equal(self.trajectory.ts.frame + 1, self.trajectory.n_frames,
                     "indexing last frame with trajectory[-1]")

    def test_slice_trajectory(self):
        frames = [ts.frame for ts in self.trajectory[5:17:3]]
        assert_equal(frames, [5, 8, 11, 14], "slicing dcd [5:17:3]")

    def test_full_slice(self):
        trj_iter = self.universe.trajectory[:]
        frames = [ts.frame for ts in trj_iter]
        assert_equal(frames, np.arange(self.universe.trajectory.n_frames))

    def test_frame_numbering(self):
        self.trajectory[98]  # index is 0-based and frames are 0-based
        assert_equal(self.universe.trajectory.frame, 98, "wrong frame number")

    def test_frame(self):
        self.trajectory[0]
        coord0 = self.universe.atoms.positions.copy()
        # forward to frame where we repeat original dcd again:
        # dcd:0..97 crd:98 dcd:99..196
        self.trajectory[99]
        assert_array_equal(
            self.universe.atoms.positions, coord0,
            "coordinates at frame 1 and 100 should be the same!")

    def test_time(self):
        self.trajectory[30]  # index and frames 0-based
        assert_almost_equal(self.universe.trajectory.time,
                            30.0,
                            5,
                            err_msg="Wrong time of frame")

    @dec.slow
    def test_write_dcd(self):
        """test that ChainReader written dcd (containing crds) is correct
        (Issue 81)"""
        W = mda.Writer(self.outfile, self.universe.atoms.n_atoms)
        for ts in self.universe.trajectory:
            W.write(self.universe)
        W.close()
        self.universe.trajectory.rewind()
        u = mda.Universe(PSF, self.outfile)
        for (ts_orig, ts_new) in zip(self.universe.trajectory,
                                     u.trajectory):
            assert_almost_equal(
                ts_orig._pos,
                ts_new._pos,
                self.prec,
                err_msg="Coordinates disagree at frame {0:d}".format(ts_orig.frame))

class TestChainReaderCommonDt(TestCase):

    @dec.skipif(parser_not_found('DCD'),
                'DCD parset not available. Are you using python 3?')
    def setUp(self):
        self.common_dt = 100.0
        self.universe = mda.Universe(PSF,
                                     [DCD, CRD, DCD, CRD, DCD, CRD, CRD],
                                     dt=self.common_dt)
        self.trajectory = self.universe.trajectory
        self.prec = 3

    def test_time(self):
        # We test this for the beginning, middle and end of the trajectory.
        for frame_n in (0, self.trajectory.n_frames/2, -1):
            self.trajectory[frame_n]
            assert_almost_equal(self.trajectory.time,
                            self.trajectory.frame*self.common_dt,
                            5,
                            err_msg="Wrong time for frame {0:d}".format(frame_n) )


class TestChainReaderFormats(TestCase):
    """Test of ChainReader with explicit formats (Issue 76)."""

    @attr('issue')
    def test_set_all_format_tuples(self):
        universe = mda.Universe(GRO, [(PDB, 'pdb'), (XTC, 'xtc'),
                                      (TRR, 'trr')])
        assert_equal(universe.trajectory.n_frames, 21)

    @attr('issue')
    @dec.skipif(parser_not_found('DCD'),
                'DCD parset not available. Are you using python 3?')
    def test_set_one_format_tuple(self):
        universe = mda.Universe(PSF, [(PDB_small, 'pdb'), DCD])
        assert_equal(universe.trajectory.n_frames, 99)

    @attr('issue')
    def test_set_all_formats(self):
        universe = mda.Universe(PSF, [PDB_small, PDB_closed], format='pdb')
        assert_equal(universe.trajectory.n_frames, 2)
