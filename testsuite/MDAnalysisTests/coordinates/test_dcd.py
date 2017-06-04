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
from __future__ import absolute_import, print_function
import MDAnalysis as mda
from MDAnalysis.coordinates.DCD import DCDReader
import numpy as np
import os
from six.moves import zip, range

from nose.plugins.attrib import attr
from numpy.testing import (assert_equal, assert_array_equal, assert_raises,
                           assert_almost_equal, assert_array_almost_equal,
                           assert_allclose, dec)
from unittest import TestCase

from MDAnalysisTests.datafiles import (DCD, PSF, DCD_empty, CRD, PRMncdf, NCDF,
                                       COORDINATES_TOPOLOGY, COORDINATES_DCD)
from MDAnalysisTests.coordinates.reference import (RefCHARMMtriclinicDCD,
                                                   RefNAMDtriclinicDCD)
# from MDAnalysisTests.coordinates.base import BaseTimestepTest

from MDAnalysisTests.coordinates.base import (MultiframeReaderTest, BaseReference,
                                              BaseWriterTest,
                                              assert_timestep_almost_equal)

from MDAnalysisTests import module_not_found, tempdir
# from MDAnalysisTests.plugins.knownfailure import knownfailure


class DCDReference(BaseReference):
    def __init__(self):
        super(DCDReference, self).__init__()
        self.trajectory = COORDINATES_DCD
        self.topology = COORDINATES_TOPOLOGY
        self.reader = mda.coordinates.DCD.DCDReader
        self.writer = mda.coordinates.DCD.DCDWriter
        self.ext = 'xtc'
        self.prec = 3
        self.changing_dimensions = True


class TestDCDReader(MultiframeReaderTest):
    def __init__(self, reference=None):
        if reference is None:
            reference = DCDReference()
        super(TestDCDReader, self).__init__(reference)

    @staticmethod
    def test_empty_dcd():
        assert_raises(IOError, mda.Universe, PSF, DCD_empty)

    @staticmethod
    def test_with_statement():

        try:
            with DCDReader(DCD) as trj:
                N = trj.n_frames
                frames = [ts.frame for ts in trj]
        except:
            raise AssertionError("with_statement not working for DCDReader")
        assert_equal(
            N,
            98,
            err_msg="with_statement: DCDReader reads wrong number of frames")
        assert_array_equal(
            frames,
            np.arange(0, N),
            err_msg="with_statement: DCDReader does not read all frames")


class TestDCDWriter(BaseWriterTest):
    def __init__(self, reference=None):
        if reference is None:
            reference = DCDReference()
        super(TestDCDWriter, self).__init__(reference)


################
# Legacy tests #
################


class TestDCDReaderClass(TestCase):
    pass


class TestDCDReader_old(TestCase):
    def setUp(self):
        self.universe = mda.Universe(PSF, DCD)
        self.dcd = self.universe.trajectory

    def tearDown(self):
        del self.universe
        del self.dcd

    def test_rewind_dcd(self):
        self.dcd.rewind()
        assert_equal(self.dcd.ts.frame, 0, "rewinding to frame 0")

    def test_next_dcd(self):
        self.dcd.rewind()
        self.dcd.next()
        assert_equal(self.dcd.ts.frame, 1, "loading frame 1")

    def test_jump_lastframe_dcd(self):
        self.dcd[-1]
        assert_equal(self.dcd.ts.frame, 97, "indexing last frame with dcd[-1]")

    def test_slice_dcd(self):
        frames = [ts.frame for ts in self.dcd[5:17:3]]
        assert_equal(frames, [5, 8, 11, 14], "slicing dcd [5:17:3]")

    def test_list_trajectory(self):
        frames = [ts.frame for ts in self.dcd[[0, 3, 4, 5]]]
        assert_equal(frames, [0, 3, 4, 5])

    def test_array_trajectory(self):
        frames = [ts.frame for ts in self.dcd[np.array([0, 3, 4, 5])]]
        assert_equal(frames, [0, 3, 4, 5])

    def test_list_reverse_trajectory(self):
        frames = [ts.frame for ts in self.dcd[[0, 4, 2, 3, 0, 1]]]
        assert_equal(frames, [0, 4, 2, 3, 0, 1])

    def test_list_repeated_trajectory(self):
        frames = [ts.frame for ts in self.dcd[[0, 0, 1, 1, 2, 1, 1]]]
        assert_equal(frames, [0, 0, 1, 1, 2, 1, 1])

    def test_reverse_dcd(self):
        frames = [ts.frame for ts in self.dcd[20:5:-1]]
        assert_equal(frames, list(range(20, 5, -1)),
                     "reversing dcd [20:5:-1]")

    def test_timeseries_slicing(self):
        # check that slicing behaves correctly
        # should  before issue #914 resolved
        x = [(0, 1, 1), (1,1,1), (1, 2, 1), (1, 4, 2), (1, 4, 4),
             (0, 5, 5), (3, 5, 1), (None, None, None)]
        for start, stop, step in x:
            yield self._slice_generation_test, start, stop, step

    def test_backwards_stepping(self):
        x = [(4, 0, -1), (5, 0, -2), (5, 0, -4)]
        for start, stop, step in x:
            yield self._failed_slices_test, start, stop, step

    def _slice_generation_test(self, start, stop, step):
        self.u = mda.Universe(PSF, DCD)
        ts = self.u.trajectory.timeseries(self.u.atoms)
        ts_skip = self.u.trajectory.timeseries(self.u.atoms, start, stop, step)
        assert_array_almost_equal(ts[:,start:stop:step], ts_skip, 5)

    def _failed_slices_test(self, start, stop, step):
        self.u = mda.Universe(PSF, DCD)
        ts = self.u.trajectory.timeseries(self.u.atoms)
        ts_skip = self.u.trajectory.timeseries(self.u.atoms, start, stop, step)
        assert_array_almost_equal(ts[:, start:stop:step,:], ts_skip, 5)


def test_DCDReader_set_dt(dt=100., frame=3):
    u = mda.Universe(PSF, DCD, dt=dt)
    assert_almost_equal(u.trajectory[frame].time, frame*dt,
                        err_msg="setting time step dt={0} failed: "
                        "actually used dt={1}".format(
            dt, u.trajectory._ts_kwargs['dt']))
    assert_almost_equal(u.trajectory.dt, dt,
                        err_msg="trajectory.dt does not match set dt")

class TestDCDWriter_old(object):
    def setUp(self):
        self.universe = mda.Universe(PSF, DCD)
        ext = ".dcd"
        self.tmpdir = tempdir.TempDir()
        self.outfile = self.tmpdir.name + '/dcd-writer' + ext
        self.Writer = mda.coordinates.DCD.DCDWriter

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except OSError:
            pass
        del self.universe
        del self.Writer
        del self.tmpdir

    def test_dt(self):
        DT = 5.0
        t = self.universe.trajectory
        with self.Writer(self.outfile,
                         t.n_atoms,
                         dt=DT) as W:  # set time step to 5 ps
            for ts in self.universe.trajectory:
                W.write_next_timestep(ts)

        uw = mda.Universe(PSF, self.outfile)
        assert_almost_equal(uw.trajectory.totaltime,
                            (uw.trajectory.n_frames - 1) * DT, 5)
        times = np.array([uw.trajectory.time for ts in uw.trajectory])
        frames = np.array([ts.frame for ts in uw.trajectory])
        assert_array_almost_equal(times, frames * DT, 5)

    def test_OtherWriter(self):
        t = self.universe.trajectory
        W = t.OtherWriter(self.outfile)
        for ts in self.universe.trajectory:
            W.write_next_timestep(ts)
        W.close()

        uw = mda.Universe(PSF, self.outfile)

        # check that the coordinates are identical for each time step
        for orig_ts, written_ts in zip(self.universe.trajectory,
                                       uw.trajectory):
            assert_array_almost_equal(written_ts._pos, orig_ts._pos, 3,
                                      err_msg="coordinate mismatch between "
                                      "original and written trajectory at "
                                      "frame %d (orig) vs %d (written)" % (
                                          orig_ts.frame, written_ts.frame))

    def test_single_frame(self):
        u = mda.Universe(PSF, CRD)
        W = mda.Writer(self.outfile, u.atoms.n_atoms)
        W.write(u.atoms)
        W.close()
        w = mda.Universe(PSF, self.outfile)
        assert_equal(w.trajectory.n_frames, 1,
                     "single frame trajectory has wrong number of frames")
        assert_almost_equal(w.atoms.positions,
                            u.atoms.positions,
                            3,
                            err_msg="coordinates do not match")


class TestDCDWriter_Issue59(object):
    def setUp(self):
        """Generate input xtc."""
        self.u = mda.Universe(PSF, DCD)
        self.tmpdir = tempdir.TempDir()
        self.xtc = self.tmpdir.name + '/dcd-writer-issue59-test.xtc'
        wXTC = mda.Writer(self.xtc, self.u.atoms.n_atoms)
        for ts in self.u.trajectory:
            wXTC.write(ts)
        wXTC.close()

    def tearDown(self):
        try:
            os.unlink(self.xtc)
        except OSError:
            pass
        try:
            os.unlink(self.dcd)
        except (AttributeError, OSError):
            pass
        del self.u
        del self.tmpdir

    @attr('issue')
    def test_issue59(self):
        """Test writing of XTC to DCD (Issue 59)"""
        xtc = mda.Universe(PSF, self.xtc)
        self.dcd = self.tmpdir.name + '/dcd-writer-issue59-test.dcd'
        wDCD = mda.Writer(self.dcd, xtc.atoms.n_atoms)
        for ts in xtc.trajectory:
            wDCD.write(ts)
        wDCD.close()

        dcd = mda.Universe(PSF, self.dcd)

        xtc.trajectory.rewind()
        dcd.trajectory.rewind()

        assert_array_almost_equal(
            xtc.atoms.positions,
            dcd.atoms.positions,
            3,
            err_msg="XTC -> DCD: DCD coordinates are messed up (Issue 59)")

    def test_OtherWriter(self):
        dcd = self.u
        wXTC = dcd.trajectory.OtherWriter(self.xtc)
        for ts in dcd.trajectory:
            wXTC.write(ts)
        wXTC.close()

        xtc = mda.Universe(PSF, self.xtc)
        xtc.trajectory.rewind()
        dcd.trajectory.rewind()

        assert_array_almost_equal(
            dcd.atoms.positions,
            xtc.atoms.positions,
            2,
            err_msg="DCD -> XTC: coordinates are messed up (frame {0:d})".format(
            dcd.trajectory.frame))
        xtc.trajectory[3]
        dcd.trajectory[3]
        assert_array_almost_equal(
            dcd.atoms.positions,
            xtc.atoms.positions,
            2,
            err_msg="DCD -> XTC: coordinates are messed up (frame {0:d})".format(
            dcd.trajectory.frame))


class _TestDCDReader_TriclinicUnitcell(TestCase):
    __test__ = False
    def setUp(self):
        self.u = mda.Universe(self.topology, self.trajectory)
        self.tempdir = tempdir.TempDir()
        self.dcd = self.tempdir.name + '/dcd-reader-triclinic.dcd'

    def tearDown(self):
        try:
            os.unlink(self.dcd)
        except (AttributeError, OSError):
            pass
        del self.u
        del self.tempdir

    @attr('issue')
    def test_read_triclinic(self):
        """test reading of triclinic unitcell (Issue 187) for NAMD or new
        CHARMM format (at least since c36b2)"""
        for ts, box in zip(self.u.trajectory,
                           self.ref_dimensions[:, 1:]):
            assert_array_almost_equal(ts.dimensions, box, 4,
                                      err_msg="box dimensions A,B,C,alpha,"
                                      "beta,gamma not identical at frame "
                                      "{}".format(ts.frame))

    @attr('issue')
    def test_write_triclinic(self):
        """test writing of triclinic unitcell (Issue 187) for NAMD or new
        CHARMM format (at least since c36b2)"""
        print("writing")
        with self.u.trajectory.OtherWriter(self.dcd) as w:
            for ts in self.u.trajectory:
                w.write(ts)
        w = mda.Universe(self.topology, self.dcd)
        print("reading\n")
        for ts_orig, ts_copy in zip(self.u.trajectory,
                                    w.trajectory):
            assert_almost_equal(ts_orig.dimensions, ts_copy.dimensions, 4,
                                err_msg="DCD->DCD: unit cell dimensions wrong "
                                "at frame {0}".format(ts_orig.frame))
        del w


class TestDCDReader_CHARMM_Unitcell(_TestDCDReader_TriclinicUnitcell,
                                    RefCHARMMtriclinicDCD):
    __test__ = True


class TestDCDReader_NAMD_Unitcell(_TestDCDReader_TriclinicUnitcell,
                                  RefNAMDtriclinicDCD):
    __test__ = True


class TestNCDF2DCD(object):
    @dec.skipif(module_not_found("netCDF4"),
                "Test skipped because netCDF is not available.")
    def setUp(self):
        self.u = mda.Universe(PRMncdf, NCDF)
        # create the DCD
        self.tmpdir = tempdir.TempDir()
        self.dcd = self.tmpdir.name + '/ncdf-2-dcd.dcd'
        DCD = mda.Writer(self.dcd, n_atoms=self.u.atoms.n_atoms)
        for ts in self.u.trajectory:
            DCD.write(ts)
        DCD.close()
        self.w = mda.Universe(PRMncdf, self.dcd)

    def tearDown(self):
        try:
            os.unlink(self.dcd)
        except (AttributeError, OSError):
            pass
        del self.u
        del self.w
        del self.tmpdir

    @attr('issue')
    def test_unitcell(self):
        """NCDFReader: Test that DCDWriter correctly writes the CHARMM
        unit cell"""
        for ts_orig, ts_copy in zip(self.u.trajectory,
                                    self.w.trajectory):
            assert_almost_equal(
                ts_orig.dimensions,
                ts_copy.dimensions,
                3,
                err_msg="NCDF->DCD: unit cell dimensions wrong at frame {0:d}".format(
                ts_orig.frame))

    def test_coordinates(self):
        for ts_orig, ts_copy in zip(self.u.trajectory,
                                    self.w.trajectory):
            assert_almost_equal(
                self.u.atoms.positions,
                self.w.atoms.positions,
                3,
                err_msg="NCDF->DCD: coordinates wrong at frame {0:d}".format(
                ts_orig.frame))


