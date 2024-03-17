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
import pytest
from unittest.mock import patch

import re
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
from numpy.testing import (assert_equal,
                           assert_almost_equal,
                           assert_allclose)

from MDAnalysisTests import make_Universe
from MDAnalysisTests.datafiles import (
    PDB_sub_dry, PDB_sub_sol, TRR_sub_sol, TRR, XTC, GRO, PDB, CRD, PRMncdf,
    NCDF, XTC_sub_sol, COORDINATES_XTC, COORDINATES_TOPOLOGY, COORDINATES_TRR)
from MDAnalysisTests.coordinates.base import (MultiframeReaderTest,
                                              BaseReference, BaseWriterTest,
                                              assert_timestep_almost_equal)

import MDAnalysis as mda
from MDAnalysis.coordinates.base import Timestep
from MDAnalysis.lib.mdamath import triclinic_box, triclinic_vectors
from MDAnalysis.coordinates import XDR
from MDAnalysisTests.util import get_userid


@pytest.mark.parametrize("filename,kwargs,reference", [
    ("foo.xtc", {}, ".foo.xtc_offsets.npz"),
    ("foo.xtc", {"ending": "npz"}, ".foo.xtc_offsets.npz"),
    ("bar.0001.trr", {"ending": "npzzzz"}, ".bar.0001.trr_offsets.npzzzz"),
])
def test_offsets_filename(filename, kwargs, reference):
    fn = XDR.offsets_filename(filename, **kwargs)
    assert fn == reference


class _XDRReader_Sub(object):
    @pytest.fixture()
    def atoms(self):
        usol = mda.Universe(PDB_sub_sol, self.XDR_SUB_SOL)
        return usol.select_atoms("not resname SOL")

    def test_load_new_raises_ValueError(self):
        # should fail if we load universe with a trajectory with different
        # number of atoms when NOT using sub, same as before.
        udry = mda.Universe(PDB_sub_dry)
        with pytest.raises(ValueError):
            udry.load_new(self.XDR_SUB_SOL)

    def test_sub_coordinates(self, atoms):
        """
        load solvated trajectory into universe with unsolvated protein.
        """
        udry = mda.Universe(PDB_sub_dry)
        udry.load_new(self.XDR_SUB_SOL, sub=atoms.indices)
        ts = udry.atoms.ts
        assert_timestep_almost_equal(ts, atoms.ts)


class TestTRRReader_Sub(_XDRReader_Sub):
    XDR_SUB_SOL = TRR_sub_sol


class TestXTCReader_Sub(_XDRReader_Sub):
    XDR_SUB_SOL = XTC_sub_sol


class _GromacsReader(object):
    # This base class assumes same lengths and dt for XTC and TRR test cases!
    filename = None
    ref_unitcell = np.array(
        [80.017, 80.017, 80.017, 60., 60., 90.], dtype=np.float32)
    # computed with Gromacs: 362.26999999999998 nm**3 * 1000 A**3/nm**3
    ref_volume = 362270.0

    prec = 3

    @pytest.fixture(scope='class')
    def universe(self):
        return mda.Universe(GRO, self.filename, convert_units=True)

    def test_rewind_xdrtrj(self, universe):
        universe.trajectory.rewind()
        assert_equal(universe.coord.frame, 0, "rewinding to frame 1")
        assert universe.trajectory._xdr._has_offsets == 1

    def test_next_xdrtrj(self, universe):
        universe.trajectory.rewind()
        universe.trajectory.next()
        assert_equal(universe.coord.frame, 1, "loading frame 1")

    def test_jump_xdrtrj(self, universe):
        universe.trajectory[4]  # index is 0-based and frames are 0-based
        assert_equal(universe.coord.frame, 4, "jumping to frame 4")

    def test_jump_lastframe_xdrtrj(self, universe):
        universe.trajectory[-1]
        assert_equal(universe.coord.frame, 9,
                     "indexing last frame with trajectory[-1]")

    def test_slice_xdrtrj(self, universe):
        frames = [ts.frame for ts in universe.trajectory[2:9:3]]
        assert_equal(frames, [2, 5, 8], "slicing xdrtrj [2:9:3]")

    def test_reverse_xdrtrj(self, universe):
        frames = [ts.frame for ts in universe.trajectory[::-1]]
        assert_equal(frames, list(range(9, -1, -1)), "slicing xdrtrj [::-1]")

    def test_coordinates(self, universe):
        ca_nm = np.array(
            [[6.043369675, 7.385184479, 1.381425762]], dtype=np.float32)
        # coordinates in the base unit (needed for True)
        ca_Angstrom = ca_nm * 10.0
        universe.trajectory.rewind()
        universe.trajectory.next()
        universe.trajectory.next()
        assert_equal(universe.coord.frame, 2, "failed to step to frame 3")
        ca = universe.select_atoms('name CA and resid 122')
        # low precision match (2 decimals in A, 3 in nm) because the above are
        # the trr coords
        assert_almost_equal(
            ca.positions,
            ca_Angstrom,
            2,
            err_msg="coords of Ca of resid 122 do not "
            "match for frame 3")

    def test_unitcell(self, universe):
        """Test that xtc/trr unitcell is read correctly (Issue 34)"""
        universe.trajectory.rewind()
        uc = universe.coord.dimensions
        assert_almost_equal(
            uc,
            self.ref_unitcell,
            self.prec,
            err_msg="unit cell dimensions (rhombic dodecahedron)")

    def test_volume(self, universe):
        # need to reduce precision for test (nm**3 <--> A**3)
        universe.trajectory.rewind()
        vol = universe.coord.volume
        assert_almost_equal(
            vol,
            self.ref_volume,
            0,
            err_msg="unit cell volume (rhombic dodecahedron)")

    def test_dt(self, universe):
        assert_almost_equal(
            universe.trajectory.dt, 100.0, 4, err_msg="wrong timestep dt")

    def test_totaltime(self, universe):
        # test_totaltime(): need to reduce precision because dt is only precise
        # to ~4 decimals and accumulating the inaccuracy leads to even lower
        # precision in the totaltime (consequence of fixing Issue 64)
        assert_almost_equal(
            universe.trajectory.totaltime,
            900.0,
            3,
            err_msg="wrong total length of trajectory")

    def test_frame(self, universe):
        universe.trajectory[4]  # index is 0-based and frames are 0-based
        assert_equal(universe.trajectory.frame, 4, "wrong frame number")

    def test_time(self, universe):
        universe.trajectory[4]
        assert_almost_equal(
            universe.trajectory.time, 400.0, 3, err_msg="wrong time of frame")

    def test_get_Writer(self, universe, tmpdir):
        ext = os.path.splitext(self.filename)[1]
        outfile = str(tmpdir.join('xdr-reader-test' + ext))
        with universe.trajectory.Writer(outfile) as W:
            assert_equal(universe.trajectory.format, W.format)
            assert_equal(universe.atoms.n_atoms, W.n_atoms)

    def test_Writer(self, tmpdir):
        universe = mda.Universe(GRO, self.filename, convert_units=True)
        ext = os.path.splitext(self.filename)[1]
        outfile = str(tmpdir.join('/xdr-reader-test' + ext))
        with universe.trajectory.Writer(outfile) as W:
            W.write(universe.atoms)
            universe.trajectory.next()
            W.write(universe.atoms)

        universe.trajectory.rewind()
        u = mda.Universe(GRO, outfile)
        assert_equal(u.trajectory.n_frames, 2)
        # prec = 6: TRR test fails; here I am generous and take self.prec =
        # 3...
        assert_almost_equal(u.atoms.positions, universe.atoms.positions,
                            self.prec)

    def test_EOFraisesStopIteration(self, universe):
        def go_beyond_EOF():
            universe.trajectory[-1]
            universe.trajectory.next()

        with pytest.raises(StopIteration):
            go_beyond_EOF()

    def test_read_next_timestep_ts_no_positions(self, universe):
        # primarily tests branching on ts.has_positions in _read_next_timestep
        ts = universe.trajectory[0]
        ts.has_positions=False
        ts_passed_in = universe.trajectory._read_next_timestep(ts=ts).copy()
        universe.trajectory.rewind()
        ts_returned = universe.trajectory._read_next_timestep(ts=None).copy()
        assert(ts_passed_in == ts_returned)

class TestXTCReader(_GromacsReader):
    filename = XTC


class TestXTCReaderClass(object):
    def test_with_statement(self):
        from MDAnalysis.coordinates.XTC import XTCReader

        try:
            with XTCReader(XTC) as trj:
                N = trj.n_frames
                frames = [ts.frame for ts in trj]
        except:
            raise AssertionError("with_statement not working for XTCReader")
        assert_equal(
            N,
            10,
            err_msg="with_statement: XTCReader reads wrong number of frames")
        assert_equal(
            frames,
            np.arange(0, N),
            err_msg="with_statement: XTCReader does not read all frames")


class TestTRRReader(_GromacsReader):
    filename = TRR

    def test_velocities(self, universe):
        # frame 0, v in nm/ps
        # from gmxdump -f MDAnalysisTests/data/adk_oplsaa.trr
        #      v[47675]={-7.86469e-01,  1.57479e+00,  2.79722e-01}
        #      v[47676]={ 2.70593e-08,  1.08052e-06,  6.97028e-07}
        v_native = np.array(
            [[-7.86469e-01, 1.57479e+00, 2.79722e-01],
             [2.70593e-08, 1.08052e-06, 6.97028e-07]],
            dtype=np.float32)

        # velocities in the MDA base unit A/ps (needed for True)
        v_base = v_native * 10.0
        universe.trajectory.rewind()
        assert_equal(universe.coord.frame, 0, "failed to read frame 1")

        assert_almost_equal(
            universe.trajectory.ts._velocities[[47675, 47676]],
            v_base,
            self.prec,
            err_msg="ts._velocities for indices 47675,47676 do not "
            "match known values")

        assert_almost_equal(
            universe.atoms.velocities[[47675, 47676]],
            v_base,
            self.prec,
            err_msg="velocities for indices 47675,47676 do not "
            "match known values")

        for index, v_known in zip([47675, 47676], v_base):
            assert_almost_equal(
                universe.atoms[index].velocity,
                v_known,
                self.prec,
                err_msg="atom[{0:d}].velocity does not match known values".
                format(index))


class _XDRNoConversion(object):
    filename = None

    @pytest.fixture()
    def universe(self):
        return mda.Universe(PDB, self.filename, convert_units=False)

    def test_coordinates(self, universe):
        # note: these are the native coordinates in nm
        ca_nm = np.array(
            [[6.043369675, 7.385184479, 1.381425762]], dtype=np.float32)
        universe.trajectory.rewind()
        universe.trajectory.next()
        universe.trajectory.next()
        assert_equal(universe.trajectory.ts.frame, 2,
                     "failed to step to frame 3")
        ca = universe.select_atoms('name CA and resid 122')
        # low precision match because we also look at the trr: only 3 decimals
        # in nm in xtc!
        assert_almost_equal(
            ca.positions,
            ca_nm,
            3,
            err_msg="native coords of Ca of resid 122 "
            "do not match for frame 3 with "
            "convert_units=False")


class TestXTCNoConversion(_XDRNoConversion):
    filename = XTC


class TestTRRNoConversion(_XDRNoConversion):
    filename = TRR


class _GromacsWriter(object):
    infilename = None  # XTC or TRR
    Writers = {
        '.trr': mda.coordinates.TRR.TRRWriter,
        '.xtc': mda.coordinates.XTC.XTCWriter,
    }

    @pytest.fixture(scope='class')
    def universe(self):
        return mda.Universe(GRO, self.infilename)

    @pytest.fixture()
    def Writer(self):
        ext = os.path.splitext(self.infilename)[1]
        return self.Writers[ext]

    @pytest.fixture()
    def outfile(self, tmpdir):
        ext = os.path.splitext(self.infilename)[1]
        return str(tmpdir.join('xdr-writer-test' + ext))

    def test_write_trajectory(self, universe, Writer, outfile):
        """Test writing Gromacs trajectories (Issue 38)"""
        with Writer(outfile, universe.atoms.n_atoms, dt=universe.trajectory.dt) as W:
            for ts in universe.trajectory:
                W.write(universe)

        uw = mda.Universe(GRO, outfile)

        # check that the coordinates are identical for each time step
        for orig_ts, written_ts in zip(universe.trajectory, uw.trajectory):
            assert_almost_equal(
                written_ts._pos,
                orig_ts._pos,
                3,
                err_msg="coordinate mismatch between "
                "original and written trajectory at "
                "frame %d (orig) vs %d (written)" % (orig_ts.frame,
                                                     written_ts.frame))

    def test_timestep_not_modified_by_writer(self, universe, Writer, outfile):
        trj = universe.trajectory
        ts = trj.ts

        trj[-1]  # last timestep (so that time != 0)
        x = ts._pos.copy()
        time = ts.time

        with Writer(outfile, trj.n_atoms, dt=trj.dt) as W:
            # last timestep (so that time != 0) (say it again, just in case...)
            trj[-1]
            W.write(universe)

        assert_equal(
            ts._pos,
            x,
            err_msg="Positions in Timestep were modified by writer.")
        assert_equal(
            ts.time, time, err_msg="Time in Timestep was modified by writer.")


class TestXTCWriter(_GromacsWriter):
    __test__ = True
    infilename = XTC


class TestTRRWriter(_GromacsWriter):
    __test__ = True
    infilename = TRR

    def test_velocities(self, universe, Writer, outfile):
        with Writer(outfile, universe.atoms.n_atoms, dt=universe.trajectory.dt) as W:
            for ts in universe.trajectory:
                W.write(universe)

        uw = mda.Universe(GRO, outfile)

        # check that the velocities are identical for each time step
        for orig_ts, written_ts in zip(universe.trajectory, uw.trajectory):
            assert_almost_equal(
                written_ts._velocities,
                orig_ts._velocities,
                3,
                err_msg="velocities mismatch between "
                "original and written trajectory at "
                "frame %d (orig) vs %d (written)" % (orig_ts.frame,
                                                     written_ts.frame))

    def test_gaps(self, universe, Writer, outfile):
        """Tests the writing and reading back of TRRs with gaps in any of
        the coordinates/velocities properties."""
        with Writer(outfile, universe.atoms.n_atoms, dt=universe.trajectory.dt) as W:
            for ts in universe.trajectory:
                # Inset some gaps in the properties: coords every 4 steps, vels
                # every 2.
                if ts.frame % 4 == 0:
                    ts.has_positions = False
                if ts.frame % 2 == 0:
                    ts.has_velocities = False
                W.write(universe)

        uw = mda.Universe(GRO, outfile)
        # check that the velocities are identical for each time step, except
        # for the gaps (that we must make sure to raise exceptions on).
        for orig_ts, written_ts in zip(universe.trajectory, uw.trajectory):
            if ts.frame % 4 != 0:
                assert_almost_equal(
                    written_ts.positions,
                    orig_ts.positions,
                    3,
                    err_msg="coordinates mismatch "
                    "between original and written "
                    "trajectory at frame {} (orig) "
                    "vs {} (written)".format(orig_ts.frame, written_ts.frame))
            else:
                with pytest.raises(mda.NoDataError):
                    getattr(written_ts, 'positions')

            if ts.frame % 2 != 0:
                assert_almost_equal(
                    written_ts.velocities,
                    orig_ts.velocities,
                    3,
                    err_msg="velocities mismatch "
                    "between original and written "
                    "trajectory at frame {} (orig) "
                    "vs {} (written)".format(orig_ts.frame, written_ts.frame))
            else:
                with pytest.raises(mda.NoDataError):
                    getattr(written_ts, 'velocities')

    def test_data_preservation(self, universe, Writer, outfile):

        with Writer(outfile, universe.atoms.n_atoms, dt=universe.trajectory.dt) as W:
            for ts in universe.trajectory:
                W.write(universe)

        uw = mda.Universe(GRO, outfile)

        assert np.isclose(ts.data['time'], 0.0)
        assert ts.data['step'] == 0
        assert np.isclose(ts.data['lambda'], 0.0)
        assert np.isclose(ts.data['dt'], 100.0)

        # check that the data are identical for each time step
        for orig_ts, written_ts in zip(universe.trajectory, uw.trajectory):
            # data lengths must be the same
            assert len(written_ts.data) == len(orig_ts.data)

            # check that the keys exist in both dictionaries
            for k in orig_ts.data:
                assert k in written_ts.data

            err_msg = ('mismatch between '
                       'original and written trajectory at '
                       f'frame {orig_ts.frame} vs {written_ts.frame}')

            # check that each value is the same
            for k in orig_ts.data:
                assert_allclose(orig_ts.data[k],
                                written_ts.data[k],
                                err_msg=err_msg)


class _GromacsWriterIssue101(object):
    Writers = {
        '.trr': mda.coordinates.TRR.TRRWriter,
        '.xtc': mda.coordinates.XTC.XTCWriter,
    }
    ext = None  # set to '.xtc' or '.trr'
    prec = 3

    @pytest.fixture()
    def Writer(self):
        return self.Writers[self.ext]

    @pytest.fixture()
    def outfile(self, tmpdir):
        return str(tmpdir.join('/xdr-writer-issue101' + self.ext))

    def test_single_frame_GRO(self, Writer, outfile):
        self._single_frame(GRO, Writer, outfile)

    def test_single_frame_PDB(self, Writer, outfile):
        self._single_frame(PDB, Writer, outfile)

    def test_single_frame_CRD(self, Writer, outfile):
        self._single_frame(CRD, Writer, outfile)

    def _single_frame(self, filename, Writer, outfile):
        u = mda.Universe(filename)
        with Writer(outfile, u.atoms.n_atoms) as W:
            W.write(u.atoms)
        w = mda.Universe(filename, outfile)
        assert_equal(w.trajectory.n_frames, 1,
                     "single frame trajectory has wrong number of frames")
        assert_almost_equal(
            w.atoms.positions,
            u.atoms.positions,
            self.prec,
            err_msg="coordinates do not match for {0!r}".format(filename))


class TestXTCWriterSingleFrame(_GromacsWriterIssue101):
    ext = ".xtc"
    prec = 2


class TestTRRWriterSingleFrame(_GromacsWriterIssue101):
    ext = ".trr"


class _GromacsWriterIssue117(object):
    """Issue 117: Cannot write XTC or TRR from AMBER NCDF"""
    ext = None
    prec = 5

    @pytest.fixture()
    def universe(self):
        return mda.Universe(PRMncdf, NCDF)

    @pytest.mark.filterwarnings("ignore: ATOMIC_NUMBER record not found")
    def test_write_trajectory(self, universe, tmpdir):
        """Test writing Gromacs trajectories from AMBER NCDF (Issue 117)"""
        outfile = str(tmpdir.join('xdr-writer-issue117' + self.ext))
        with mda.Writer(outfile, n_atoms=universe.atoms.n_atoms) as W:
            for ts in universe.trajectory:
                W.write(universe)

        uw = mda.Universe(PRMncdf, outfile)

        # check that the coordinates are identical for each time step
        for orig_ts, written_ts in zip(universe.trajectory, uw.trajectory):
            assert_almost_equal(
                written_ts._pos,
                orig_ts._pos,
                self.prec,
                err_msg=("coordinate mismatch between original and written "
                         f"trajectory at frame {orig_ts.frame:d} (orig) vs "
                         f"{orig_ts.frame:d} (written)"))


class TestXTCWriterIssue117(_GromacsWriterIssue117):
    __test__ = True
    ext = ".xtc"
    prec = 2


class TestTRRWriterIssue117(_GromacsWriterIssue117):
    __test__ = True
    ext = ".trr"


def test_triclinic_box():
    """Test coordinates.core.triclinic_box() (Issue 61)"""
    unitcell = np.array([80.017, 55, 100.11, 60.00, 30.50, 90.00])
    box = triclinic_vectors(unitcell)
    new_unitcell = triclinic_box(box[0], box[1], box[2])
    assert_almost_equal(
        new_unitcell,
        unitcell,
        3,
        err_msg="unitcell round-trip connversion failed (Issue 61)")


class XTCReference(BaseReference):
    def __init__(self):
        super(XTCReference, self).__init__()
        self.trajectory = COORDINATES_XTC
        self.topology = COORDINATES_TOPOLOGY
        self.reader = mda.coordinates.XTC.XTCReader
        self.writer = mda.coordinates.XTC.XTCWriter
        self.ext = 'xtc'
        self.prec = 3
        self.changing_dimensions = True


class TestXTCReader_2(MultiframeReaderTest):
    @staticmethod
    @pytest.fixture()
    def ref():
        return XTCReference()


class TestXTCWriter_2(BaseWriterTest):
    @staticmethod
    @pytest.fixture()
    def ref():
        return XTCReference()

    def test_different_precision(self, ref, tmpdir):
        out = 'precision-test' + ref.ext
        # store more then 9 atoms to enable compression
        n_atoms = 40
        with tmpdir.as_cwd():
            with ref.writer(out, n_atoms, precision=5) as w:
                u = make_Universe(size=(n_atoms, 1, 1), trajectory=True)
                u.trajectory.ts.positions = np.random.random(size=(n_atoms, 3))
                w.write(u)
            xtc = mda.lib.formats.libmdaxdr.XTCFile(out)
            frame = xtc.read()
            assert_equal(len(xtc), 1)
            assert_equal(xtc.n_atoms, n_atoms)
            assert_equal(frame.prec, 10.0**5)


class TRRReference(BaseReference):
    def __init__(self):
        super(TRRReference, self).__init__()
        self.trajectory = COORDINATES_TRR
        self.topology = COORDINATES_TOPOLOGY
        self.changing_dimensions = True
        self.reader = mda.coordinates.TRR.TRRReader
        self.writer = mda.coordinates.TRR.TRRWriter
        self.ext = 'trr'
        self.prec = 3
        self.first_frame.velocities = self.first_frame.positions / 10
        self.first_frame.forces = self.first_frame.positions / 100

        self.second_frame.velocities = self.second_frame.positions / 10
        self.second_frame.forces = self.second_frame.positions / 100

        self.last_frame.velocities = self.last_frame.positions / 10
        self.last_frame.forces = self.last_frame.positions / 100

        self.jump_to_frame.velocities = self.jump_to_frame.positions / 10
        self.jump_to_frame.forces = self.jump_to_frame.positions / 100

    def iter_ts(self, i):
        ts = self.first_frame.copy()
        ts.positions = 2**i * self.first_frame.positions
        ts.velocities = ts.positions / 10
        ts.forces = ts.positions / 100
        ts.time = i
        ts.frame = i
        return ts


class TestTRRReader_2(MultiframeReaderTest):
    @staticmethod
    @pytest.fixture()
    def ref():
        return TRRReference()


class TestTRRWriter_2(BaseWriterTest):
    @staticmethod
    @pytest.fixture()
    def ref():
        return TRRReference()

    # tests writing and reading in one!
    def test_lambda(self, ref, universe, tmpdir):
        outfile = 'write-lambda-test' + ref.ext

        with tmpdir.as_cwd():
            with ref.writer(outfile, universe.trajectory.n_atoms) as W:
                for i, ts in enumerate(universe.trajectory):
                    ts.data['lambda'] = i / float(universe.trajectory.n_frames)
                    W.write(universe)

            reader = ref.reader(outfile)
            for i, ts in enumerate(reader):
                assert_almost_equal(ts.data['lambda'], i / float(reader.n_frames))


class _GromacsReader_offsets(object):

    # This base class assumes same lengths and dt for XTC and TRR test cases!
    filename = None
    ref_unitcell = np.array(
        [80.017, 80.017, 80.017, 60., 60., 90.], dtype=np.float32)
    # computed with Gromacs: 362.26999999999998 nm**3 * 1000 A**3/nm**3
    ref_volume = 362270.0
    ref_offsets = None
    _reader = None
    prec = 3

    @pytest.fixture(scope='class')
    def traj(self, tmpdir_factory):
        # copy of original test trajectory in a temporary folder. This is
        # needed since offsets are automatically generated in the same
        # directory. Here we also clean up nicely all files we generate
        tmpdir = tmpdir_factory.mktemp('xtc')
        shutil.copy(self.filename, str(tmpdir))
        traj = str(tmpdir.join(os.path.basename(self.filename)))
        # ensure initialization of offsets
        self._reader(traj)
        return traj

    @pytest.fixture()
    def trajectory(self, traj):
        return self._reader(traj)

    def test_offsets(self, trajectory, traj):
        trajectory._read_offsets(store=True)
        assert_almost_equal(
            trajectory._xdr.offsets,
            self.ref_offsets,
            err_msg="wrong frame offsets")

        outfile_offsets = XDR.offsets_filename(traj)
        saved_offsets = XDR.read_numpy_offsets(outfile_offsets)

        assert isinstance(saved_offsets, dict), \
            "read_numpy_offsets did not return a dict"
        assert_almost_equal(
            trajectory._xdr.offsets,
            saved_offsets['offsets'],
            err_msg="error saving frame offsets")
        assert_almost_equal(
            self.ref_offsets,
            saved_offsets['offsets'],
            err_msg="saved frame offsets don't match "
            "the known ones")

        trajectory._load_offsets()
        assert_almost_equal(
            trajectory._xdr.offsets,
            self.ref_offsets,
            err_msg="error loading frame offsets")
        assert_equal(saved_offsets['ctime'], os.path.getctime(traj))
        assert_equal(saved_offsets['size'], os.path.getsize(traj))

    def test_reload_offsets(self, traj):
        self._reader(traj, refresh_offsets=True)

    def test_nonexistent_offsets_file(self, traj):
        # assert that a nonexistent file returns False during read-in
        outfile_offsets = XDR.offsets_filename(traj)
        with patch.object(np, "load") as np_load_mock:
            np_load_mock.side_effect = IOError
            with pytest.warns(UserWarning, match=re.escape(
                    f"Failed to load offsets file {outfile_offsets}")):
                saved_offsets = XDR.read_numpy_offsets(outfile_offsets)
        assert saved_offsets == False

    def test_corrupted_offsets_file(self, traj):
        # assert that a corrupted file returns False during read-in
        # Issue #3230
        outfile_offsets = XDR.offsets_filename(traj)
        with patch.object(np, "load") as np_load_mock:
            np_load_mock.side_effect = ValueError
            with pytest.warns(UserWarning, match=re.escape(
                    f"Failed to load offsets file {outfile_offsets}")):
                saved_offsets = XDR.read_numpy_offsets(outfile_offsets)
        assert saved_offsets == False

    def test_reload_offsets_if_offsets_readin_io_fails(self, trajectory):
        # force the np.load call that is called in read_numpy_offsets
        # during _load_offsets to give an IOError
        # ensure that offsets are then read-in from the trajectory
        with patch.object(np, "load") as np_load_mock:
            np_load_mock.side_effect = IOError
            with (pytest.warns(UserWarning,
                               match="Failed to load offsets file") and
                  pytest.warns(UserWarning,
                               match="reading offsets from trajectory instead")):
                trajectory._load_offsets()

            assert_almost_equal(
                trajectory._xdr.offsets,
                self.ref_offsets,
                err_msg="error loading frame offsets")

    def test_reload_offsets_if_offsets_readin_value_fails(self, trajectory):
        # force the np.load call that is called in read_numpy_offsets
        # during _load_offsets to give an ValueError (Issue #3230)
        # ensure that offsets are then read-in from the trajectory
        with patch.object(np, "load") as np_load_mock:
            np_load_mock.side_effect = ValueError
            with pytest.warns(UserWarning, match="Failed to load offsets"):
                trajectory._load_offsets()
            assert_almost_equal(
                trajectory._xdr.offsets,
                self.ref_offsets,
                err_msg="error loading frame offsets")

    def test_persistent_offsets_size_mismatch(self, traj):
        # check that stored offsets are not loaded when trajectory
        # size differs from stored size
        fname = XDR.offsets_filename(traj)
        saved_offsets = XDR.read_numpy_offsets(fname)

        assert isinstance(saved_offsets, dict), \
            "read_numpy_offsets did not return a dict"

        saved_offsets['size'] += 1
        with open(fname, 'wb') as f:
            np.savez(f, **saved_offsets)

        with pytest.warns(UserWarning, match="Reload offsets"):
            self._reader(traj)

    def test_persistent_offsets_ctime_mismatch(self, traj):
        # check that stored offsets are not loaded when trajectory
        # ctime differs from stored ctime
        fname = XDR.offsets_filename(traj)
        saved_offsets = XDR.read_numpy_offsets(fname)

        assert isinstance(saved_offsets, dict), \
            "read_numpy_offsets did not return a dict"

        saved_offsets['ctime'] += 1
        with open(fname, 'wb') as f:
            np.savez(f, **saved_offsets)

        with pytest.warns(UserWarning, match="Reload offsets"):
            self._reader(traj)

    def test_persistent_offsets_natoms_mismatch(self, traj):
        # check that stored offsets are not loaded when trajectory
        # ctime differs from stored ctime
        fname = XDR.offsets_filename(traj)
        saved_offsets = XDR.read_numpy_offsets(fname)

        assert isinstance(saved_offsets, dict), \
            "read_numpy_offsets did not return a dict"

        saved_offsets['n_atoms'] += 1
        np.savez(fname, **saved_offsets)

        with pytest.warns(UserWarning, match="Reload offsets"):
            self._reader(traj)

    def test_persistent_offsets_last_frame_wrong(self, traj):
        fname = XDR.offsets_filename(traj)
        saved_offsets = XDR.read_numpy_offsets(fname)

        assert isinstance(saved_offsets, dict), \
            "read_numpy_offsets did not return a dict"

        idx_frame = 3
        saved_offsets['offsets'][idx_frame] += 42
        np.savez(fname, **saved_offsets)

        with pytest.warns(UserWarning, match="seek failed"):
            reader = self._reader(traj)
            reader[idx_frame]

    def test_unsupported_format(self, traj):
        fname = XDR.offsets_filename(traj)
        saved_offsets = XDR.read_numpy_offsets(fname)

        assert isinstance(saved_offsets, dict), \
            "read_numpy_offsets did not return a dict"

        idx_frame = 3
        saved_offsets.pop('n_atoms')
        np.savez(fname, **saved_offsets)

        # ok as long as this doesn't throw
        with pytest.warns(UserWarning, match="Reload offsets from trajectory"):
            reader = self._reader(traj)
        reader[idx_frame]

    @pytest.mark.skipif(get_userid() == 0, reason="cannot readonly as root")
    def test_persistent_offsets_readonly(self, tmpdir):
        shutil.copy(self.filename, str(tmpdir))

        if os.name == 'nt':
            # Windows platform has a unique way to deny write access
            subprocess.call("icacls {fname} /deny Users:W".format(fname=tmpdir),
                            shell=True)
        else:
            os.chmod(str(tmpdir), 0o555)

        filename = str(tmpdir.join(os.path.basename(self.filename)))
        # try to write a offsets file
        with (pytest.warns(UserWarning, match="Couldn't save offsets") and
              pytest.warns(UserWarning, match="Cannot write")):
            self._reader(filename)
        assert_equal(os.path.exists(XDR.offsets_filename(filename)), False)
        # check the lock file is not created as well.
        assert_equal(os.path.exists(XDR.offsets_filename(filename,
                                                    ending='.lock')), False)

        # pre-teardown permission fix - leaving permission blocked dir
        # is problematic on py3.9 + Windows it seems. See issue
        # [4123](https://github.com/MDAnalysis/mdanalysis/issues/4123)
        # for more details.
        if os.name == 'nt':
            subprocess.call(f"icacls {tmpdir} /grant Users:W", shell=True)
        else:
            os.chmod(str(tmpdir), 0o777)

        shutil.rmtree(tmpdir)

    def test_offset_lock_created(self):
        assert os.path.exists(XDR.offsets_filename(self.filename,
                                                   ending='lock'))


class TestXTCReader_offsets(_GromacsReader_offsets):
    __test__ = True
    filename = XTC
    ref_offsets = np.array([
        0, 165188, 330364, 495520, 660708, 825872, 991044, 1156212, 1321384,
        1486544
    ])
    _reader = mda.coordinates.XTC.XTCReader


class TestTRRReader_offsets(_GromacsReader_offsets):
    __test__ = True
    filename = TRR
    ref_offsets = np.array([
        0, 1144464, 2288928, 3433392, 4577856, 5722320, 6866784, 8011248,
        9155712, 10300176
    ])
    _reader = mda.coordinates.TRR.TRRReader


def test_pathlib():
    # regression test for XDR path of
    # gh-2497
    top = Path(GRO)
    traj = Path(XTC)
    u = mda.Universe(top, traj)
    # we really only care that pathlib
    # object handling worked
    assert u.atoms.n_atoms == 47681
