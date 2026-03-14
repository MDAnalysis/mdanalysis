# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
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
import sqlite3

import MDAnalysis as mda
import numpy as np
import pytest

from numpy.testing import assert_allclose, assert_equal

from MDAnalysis.lib.mdamath import triclinic_vectors

from MDAnalysisTests.datafiles import DMS


def _make_dms_with_box(path, box_vectors):
    """Create a minimal DMS SQLite file with one atom and real box vectors.

    Parameters
    ----------
    path : str or Path
        Destination file path.
    box_vectors : list of list
        Three row vectors [[x1,y1,z1],[x2,y2,z2],[x3,y3,z3]] for global_cell.
    """
    with sqlite3.connect(str(path)) as con:
        cur = con.cursor()
        cur.execute(
            "CREATE TABLE particle ("
            "id integer primary key, anum integer, "
            "x float, y float, z float, "
            "vx float, vy float, vz float, "
            "mass float, charge float, "
            "name text, resname text, resid integer, chain text, segid text)"
        )
        cur.execute(
            "INSERT INTO particle VALUES (0, 7, 1.0, 2.0, 3.0, "
            "0.0, 0.0, 0.0, 14.007, -0.3, 'N', 'MET', 1, 'A', 'SYSTEM')"
        )
        cur.execute(
            "CREATE TABLE global_cell "
            "(id integer primary key, x float, y float, z float)"
        )
        for i, row in enumerate(box_vectors, start=1):
            cur.execute(
                "INSERT INTO global_cell VALUES (?, ?, ?, ?)",
                (i, row[0], row[1], row[2]),
            )
        cur.execute(
            "CREATE TABLE bond (p0 integer, p1 integer, 'order' integer)"
        )


class TestDMSReader(object):
    @pytest.fixture()
    def universe(self):
        return mda.Universe(DMS)

    @pytest.fixture()
    def ts(self, universe):
        return universe.trajectory.ts

    def test_global_cell(self, ts):
        assert ts.dimensions is None

    # cythonised class can no longer raise AttributeError
    # so changed to test of has_velocities
    def test_velocities(self, ts):
        assert_equal(ts.has_velocities, False)

    def test_number_of_coords(self, universe):
        # Desired value taken from VMD
        #      Info)    Atoms: 3341
        assert_equal(len(universe.atoms), 3341)

    def test_coords_atom_0(self, universe):
        # Desired coordinates taken directly from the SQLite file. Check unit
        # conversion
        coords_0 = np.array(
            [
                -11.0530004501343,
                26.6800003051758,
                12.7419996261597,
            ],
            dtype=np.float32,
        )
        assert_equal(universe.atoms[0].position, coords_0)

    def test_n_frames(self, universe):
        assert_equal(
            universe.trajectory.n_frames, 1, "wrong number of frames in pdb"
        )

    def test_time(self, universe):
        assert_equal(universe.trajectory.time, 0.0, "wrong time of the frame")

    def test_frame(self, universe):
        assert_equal(
            universe.trajectory.frame,
            0,
            "wrong frame number "
            "(0-based, should be 0 for single frame readers)",
        )

    def test_frame_index_0(self, universe):
        universe.trajectory[0]
        assert_equal(
            universe.trajectory.ts.frame,
            0,
            "frame number for frame index 0 should be 0",
        )

    def test_frame_index_1_raises_IndexError(self, universe):
        with pytest.raises(IndexError):
            universe.trajectory[1]


class TestDMSReaderWithBox(object):
    """Tests for the convert_pos_from_native code path that runs when a DMS
    file contains a non-zero unit cell (ts.dimensions is not None)."""

    @pytest.fixture()
    def universe_with_box(self, tmp_path):
        dms_path = tmp_path / "box.dms"
        # Orthorhombic 30 x 40 x 50 Angstrom box.
        box_vectors = [
            [30.0, 0.0, 0.0],
            [0.0, 40.0, 0.0],
            [0.0, 0.0, 50.0],
        ]
        _make_dms_with_box(dms_path, box_vectors)
        return mda.Universe(str(dms_path))

    def test_dimensions_not_none(self, universe_with_box):
        assert universe_with_box.trajectory.ts.dimensions is not None

    def test_dimensions_values(self, universe_with_box):
        expected = np.array([30.0, 40.0, 50.0, 90.0, 90.0, 90.0], dtype=np.float32)
        assert_allclose(
            universe_with_box.trajectory.ts.dimensions,
            expected,
            atol=1e-4,
        )

    def test_position_read_with_box(self, universe_with_box):
        # Verify that coordinates are read correctly when a box is present.
        expected = np.array([1.0, 2.0, 3.0], dtype=np.float32)
        assert_allclose(
            universe_with_box.atoms[0].position,
            expected,
            atol=1e-4,
        )
