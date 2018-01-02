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
from __future__ import absolute_import, print_function, division
from six.moves import zip, range
import numpy as np

import MDAnalysis as mda
from MDAnalysis.coordinates.DCD import DCDReader
from MDAnalysis.exceptions import NoDataError

from numpy.testing import (assert_equal, assert_array_equal,
                           assert_almost_equal, assert_array_almost_equal)

from MDAnalysisTests.datafiles import (DCD, PSF, DCD_empty, PRMncdf, NCDF,
                                       COORDINATES_TOPOLOGY, COORDINATES_DCD,
                                       PSF_TRICLINIC, DCD_TRICLINIC,
                                       PSF_NAMD_TRICLINIC, DCD_NAMD_TRICLINIC)
from MDAnalysisTests.coordinates.base import (MultiframeReaderTest,
                                              BaseReference,
                                              BaseWriterTest)

import pytest


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
    @staticmethod
    @pytest.fixture()
    def ref():
        return DCDReference()

    def test_empty_dcd(self):
        with pytest.raises(IOError):
            mda.Universe(PSF, DCD_empty)

    def test_with_statement(self):

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
    @staticmethod
    @pytest.fixture()
    def ref():
        return DCDReference()


def test_write_random_unitcell(tmpdir):
    testname = str(tmpdir.join('test.dcd'))
    rstate = np.random.RandomState(1178083)
    random_unitcells = rstate.uniform(
        high=80, size=(98, 6)).astype(np.float64)

    u = mda.Universe(PSF, DCD)
    with mda.Writer(testname, n_atoms=u.atoms.n_atoms) as w:
        for index, ts in enumerate(u.trajectory):
            u.atoms.dimensions = random_unitcells[index]
            w.write(u.atoms)

    u2 = mda.Universe(PSF, testname)
    for index, ts in enumerate(u2.trajectory):
        assert_array_almost_equal(u2.trajectory.dimensions,
                                  random_unitcells[index],
                                  decimal=5)


################
# Legacy tests #
################

@pytest.fixture(scope='module')
def universe_dcd():
    return mda.Universe(PSF, DCD)


def test_rewind(universe_dcd):
    universe_dcd.trajectory.rewind()
    assert universe_dcd.trajectory.ts.frame == 0


def test_next(universe_dcd):
    universe_dcd.trajectory.rewind()
    universe_dcd.trajectory.next()
    assert universe_dcd.trajectory.ts.frame == 1


def test_jump_last_frame(universe_dcd):
    universe_dcd.trajectory[-1]
    assert universe_dcd.trajectory.ts.frame == 97


@pytest.mark.parametrize("start, stop, step", ((5, 17, 3),
                                               (20, 5, -1)))
def test_slice(universe_dcd, start, stop, step):
    frames = [ts.frame for ts in universe_dcd.trajectory[start:stop:step]]
    assert_array_equal(frames, np.arange(start, stop, step))


@pytest.mark.parametrize("array_like", [list, np.array])
def test_array_like(universe_dcd, array_like):
    ar = array_like([0, 3, 4, 5])
    frames = [ts.frame for ts in universe_dcd.trajectory[ar]]
    assert_array_equal(frames, ar)


@pytest.mark.parametrize("indices", ([0, 4, 2, 3, 0, 1],
                                     [0, 0, 1, 1, 2, 1, 1]))
def test_list_indices(universe_dcd, indices):
    frames = [ts.frame for ts in universe_dcd.trajectory[indices]]
    assert_array_equal(frames, indices)


@pytest.mark.parametrize(
    "slice, length",
    [([None, None, None], 98), ([0, None, None], 98), ([None, 98, None], 98),
     ([None, None, 1], 98), ([None, None, -1], 98), ([2, 6, 2], 2),
     ([0, 10, None], 10), ([2, 10, None], 8), ([0, 1, 1], 1), ([1, 1, 1], 0),
     ([1, 2, 1], 1), ([1, 2, 2], 1), ([1, 4, 2], 2), ([1, 4, 4], 1),
     ([0, 5, 5], 1), ([3, 5, 1], 2), ([4, 0, -1], 4), ([5, 0, -2], 3),
     ([5, 0, -4], 2)])
def test_timeseries_slices(slice, length, universe_dcd):
    start, stop, step = slice
    allframes = universe_dcd.trajectory.timeseries(order='fac')
    xyz = universe_dcd.trajectory.timeseries(start=start, stop=stop, step=step,
                                             order='fac')
    assert len(xyz) == length
    assert_array_almost_equal(xyz, allframes[start:stop:step])


def test_timeseries_deprecation(universe_dcd):
    with pytest.warns(DeprecationWarning):
        universe_dcd.trajectory.timeseries(format='fac')


@pytest.mark.parametrize("order, shape", (
    ('fac', (98, 3341, 3)),
    ('fca', (98, 3, 3341)),
    ('afc', (3341, 98, 3)),
    ('acf', (3341, 3, 98)),
    ('caf', (3, 3341, 98)),
    ('cfa', (3, 98, 3341)), ))
def test_timeseries_order(order, shape, universe_dcd):
    x = universe_dcd.trajectory.timeseries(order=order)
    assert x.shape == shape


@pytest.mark.parametrize("indices", [[1, 2, 3, 4], [5, 10, 15, 19],
                                     [9, 4, 2, 0, 50]])
def test_timeseries_atomindices(indices, universe_dcd):
        allframes = universe_dcd.trajectory.timeseries(order='afc')
        asel = universe_dcd.atoms[indices]
        xyz = universe_dcd.trajectory.timeseries(asel=asel, order='afc')
        assert len(xyz) == len(indices)
        assert_array_almost_equal(xyz, allframes[indices])


def test_timeseries_empty_selection(universe_dcd):
    with pytest.raises(NoDataError):
        asel = universe_dcd.select_atoms('name FOO')
        universe_dcd.trajectory.timeseries(asel=asel)


def test_timeseries_skip(universe_dcd):
    with pytest.warns(DeprecationWarning):
        xyz = universe_dcd.trajectory.timeseries(skip=2, order='fac')
    assert len(xyz) == universe_dcd.trajectory.n_frames / 2


def test_reader_set_dt():
    dt = 100
    frame = 3
    u = mda.Universe(PSF, DCD, dt=dt)
    assert_almost_equal(u.trajectory[frame].time, frame*dt,
                        err_msg="setting time step dt={0} failed: "
                        "actually used dt={1}".format(
                            dt, u.trajectory._ts_kwargs['dt']))
    assert_almost_equal(u.trajectory.dt, dt,
                        err_msg="trajectory.dt does not match set dt")


@pytest.mark.parametrize("ext, decimal", (("dcd", 5),
                                          ("xtc", 3)))
def test_writer_dt(tmpdir, ext, decimal):
    dt = 5.0  # set time step to 5 ps
    universe_dcd = mda.Universe(PSF, DCD, dt=dt)
    t = universe_dcd.trajectory
    outfile = str(tmpdir.join("test.{}".format(ext)))
    with mda.Writer(outfile, n_atoms=t.n_atoms, dt=dt) as W:
        for ts in universe_dcd.trajectory:
            W.write(universe_dcd.atoms)

    uw = mda.Universe(PSF, outfile)
    assert_almost_equal(uw.trajectory.totaltime,
                        (uw.trajectory.n_frames - 1) * dt, decimal)
    times = np.array([uw.trajectory.time for ts in uw.trajectory])
    frames = np.array([ts.frame for ts in uw.trajectory])
    assert_array_almost_equal(times, frames * dt, decimal)


@pytest.mark.parametrize("ext, decimal", (("dcd", 5),
                                          ("xtc", 2)))
def test_other_writer(universe_dcd, tmpdir, ext, decimal):
    t = universe_dcd.trajectory
    outfile = str(tmpdir.join("test.{}".format(ext)))
    with t.OtherWriter(outfile) as W:
        for ts in universe_dcd.trajectory:
            W.write_next_timestep(ts)

    uw = mda.Universe(PSF, outfile)
    # check that the coordinates are identical for each time step
    for orig_ts, written_ts in zip(universe_dcd.trajectory,
                                   uw.trajectory):
        assert_array_almost_equal(written_ts.positions, orig_ts.positions,
                                  decimal,
                                  err_msg="coordinate mismatch between "
                                  "original and written trajectory at "
                                  "frame {} (orig) vs {} (written)".format(
                                      orig_ts.frame, written_ts.frame))


def test_single_frame(universe_dcd, tmpdir):
    u = universe_dcd
    outfile = str(tmpdir.join("test.dcd"))
    with mda.Writer(outfile, u.atoms.n_atoms) as W:
        W.write(u.atoms)
    w = mda.Universe(PSF, outfile)
    assert w.trajectory.n_frames == 1
    assert_almost_equal(w.atoms.positions,
                        u.atoms.positions,
                        3,
                        err_msg="coordinates do not match")


def test_write_no_natoms():
    with pytest.raises(ValueError):
        mda.Writer('foobar.dcd')


def test_writer_trajectory_no_natoms(tmpdir, universe_dcd):
    with tmpdir.as_cwd():
        universe_dcd.trajectory.Writer("foo.dcd")


class RefCHARMMtriclinicDCD(object):
    topology = PSF_TRICLINIC
    trajectory = DCD_TRICLINIC
    # time(ps) A B C alpha beta gamma (length in Angstrome, angles in degrees)
    # dcd starts at t = 1ps
    ref_dimensions = np.array([
        [1., 35.44604, 35.06156, 34.1585, 91.32802, 61.73521, 44.40703],
        [2., 34.65957, 34.22689, 33.09897, 90.56206, 61.79192, 44.14549],
        [3., 34.52772, 34.66422, 33.53881, 90.55859, 63.11228, 40.14044],
        [4., 34.43749, 33.38432, 34.02133, 88.82457, 64.98057, 36.77397],
        [5., 33.73129, 32.47752, 34.18961, 89.88102, 65.89032, 36.10921],
        [6., 33.78703, 31.90317, 34.98833, 90.03092, 66.12877, 35.07141],
        [7., 33.24708, 31.18271, 34.9654, 93.11122, 68.17743, 35.73643],
        [8., 32.92599, 30.31393, 34.99197, 93.89051, 69.3799, 33.48945],
        [9., 32.15295, 30.43056, 34.96157, 96.01416, 71.50115, 32.56111],
        [10., 31.99748, 30.21518, 35.24292, 95.85821, 71.08429, 31.85939]
    ])


class RefNAMDtriclinicDCD(object):
    topology = PSF_NAMD_TRICLINIC
    trajectory = DCD_NAMD_TRICLINIC
    # vmd topology trajectory
    # molinfo 0 get {a b c alpha beta gamma}
    # time(ps) A B C alpha beta gamma (length in Angstrome, angles in degrees)
    ref_dimensions = np.array([
        [1., 38.426594, 38.393101, 44.759800, 90.000000, 90.000000, 60.028915],
    ])


@pytest.mark.parametrize("ref", (RefCHARMMtriclinicDCD, RefNAMDtriclinicDCD))
def test_read_unitcell_triclinic(ref):
    u = mda.Universe(ref.topology, ref.trajectory)
    for ts, box in zip(u.trajectory, ref.ref_dimensions[:, 1:]):
        assert_array_almost_equal(ts.dimensions, box, 4,
                                  err_msg="box dimensions A,B,C,alpha,"
                                  "beta,gamma not identical at frame "
                                  "{}".format(ts.frame))


@pytest.mark.parametrize("ref", (RefCHARMMtriclinicDCD, RefNAMDtriclinicDCD))
def test_write_unitcell_triclinic(ref, tmpdir):
    u = mda.Universe(ref.topology, ref.trajectory)
    outfile = 'triclinic.dcd'
    with tmpdir.as_cwd():
        with u.trajectory.OtherWriter(outfile) as w:
            for ts in u.trajectory:
                w.write(ts)

        w = mda.Universe(ref.topology, outfile)
        for ts_orig, ts_copy in zip(u.trajectory, w.trajectory):
            assert_almost_equal(ts_orig.dimensions, ts_copy.dimensions, 4,
                                err_msg="DCD->DCD: unit cell dimensions wrong "
                                "at frame {0}".format(ts_orig.frame))


@pytest.fixture(scope='module')
def ncdf2dcd(tmpdir_factory):
    pytest.importorskip("netCDF4")
    testfile = tmpdir_factory.mktemp('dcd').join('ncdf2dcd.dcd')
    testfile = str(testfile)
    ncdf = mda.Universe(PRMncdf, NCDF)
    with mda.Writer(testfile, n_atoms=ncdf.atoms.n_atoms) as w:
        for ts in ncdf.trajectory:
            w.write(ts)
    return ncdf, mda.Universe(PRMncdf, testfile)


def test_ncdf2dcd_unitcell(ncdf2dcd):
    ncdf, dcd = ncdf2dcd
    for ts_ncdf, ts_dcd in zip(ncdf.trajectory, dcd.trajectory):
        assert_almost_equal(ts_ncdf.dimensions,
                            ts_dcd.dimensions,
                            3)


def test_ncdf2dcd_coords(ncdf2dcd):
    ncdf, dcd = ncdf2dcd
    for ts_ncdf, ts_dcd in zip(ncdf.trajectory, dcd.trajectory):
        assert_almost_equal(ts_ncdf.positions,
                            ts_dcd.positions,
                            3)
