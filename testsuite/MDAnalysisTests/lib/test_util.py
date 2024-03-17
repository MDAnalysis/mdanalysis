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
from io import StringIO
import pytest
import os
import warnings
import re
import textwrap
from unittest.mock import Mock, patch
import sys
import copy
import shutil

import numpy as np
from numpy.testing import (assert_equal, assert_almost_equal,
                           assert_array_almost_equal, assert_array_equal,
                           assert_allclose)
from itertools import combinations_with_replacement as comb_wr

import MDAnalysis as mda
from MDAnalysis.core._get_readers import get_parser_for
import MDAnalysis.lib.util as util
import MDAnalysis.lib.mdamath as mdamath
from MDAnalysis.lib.util import (cached, static_variables, warn_if_not_unique,
                                 check_coords, store_init_arguments, 
                                 check_atomgroup_not_empty,)
from MDAnalysis.core.topologyattrs import Bonds
from MDAnalysis.exceptions import NoDataError, DuplicateWarning
from MDAnalysis.core.groups import AtomGroup
from MDAnalysisTests.datafiles import (PSF, DCD,
    Make_Whole, TPR, GRO, fullerene, two_water_gro,
)

def test_absence_cutil():
    with patch.dict('sys.modules', {'MDAnalysis.lib._cutil':None}):
        import importlib
        with pytest.raises(ImportError):
            importlib.reload(sys.modules['MDAnalysis.lib.util'])

def test_presence_cutil():
    mock = Mock()
    with patch.dict('sys.modules', {'MDAnalysis.lib._cutil':mock}):
        try:
            import MDAnalysis.lib._cutil
        except ImportError:
            pytest.fail(msg='''MDAnalysis.lib._cutil should not raise
                         an ImportError if cutil is available.''')

def convert_aa_code_long_data():
    aa = [
        ('H',
         ('HIS', 'HISA', 'HISB', 'HSE', 'HSD', 'HIS1', 'HIS2', 'HIE', 'HID')),
        ('K', ('LYS', 'LYSH', 'LYN')),
        ('A', ('ALA',)),
        ('D', ('ASP', 'ASPH', 'ASH')),
        ('E', ('GLU', 'GLUH', 'GLH')),
        ('N', ('ASN',)),
        ('Q', ('GLN',)),
        ('C', ('CYS', 'CYSH', 'CYS1', 'CYS2')),
    ]
    for resname1, strings in aa:
        for resname3 in strings:
            yield (resname3, resname1)


class TestStringFunctions(object):
    # (1-letter, (canonical 3 letter, other 3/4 letter, ....))
    aa = [
        ('H',
         ('HIS', 'HISA', 'HISB', 'HSE', 'HSD', 'HIS1', 'HIS2', 'HIE', 'HID')),
        ('K', ('LYS', 'LYSH', 'LYN')),
        ('A', ('ALA',)),
        ('D', ('ASP', 'ASPH', 'ASH')),
        ('E', ('GLU', 'GLUH', 'GLH')),
        ('N', ('ASN',)),
        ('Q', ('GLN',)),
        ('C', ('CYS', 'CYSH', 'CYS1', 'CYS2')),
    ]

    residues = [
        ("LYS300:HZ1", ("LYS", 300, "HZ1")),
        ("K300:HZ1", ("LYS", 300, "HZ1")),
        ("K300", ("LYS", 300, None)),
        ("LYS 300:HZ1", ("LYS", 300, "HZ1")),
        ("M1:CA", ("MET", 1, "CA")),
    ]

    @pytest.mark.parametrize('rstring, residue', residues)
    def test_parse_residue(self, rstring, residue):
        assert util.parse_residue(rstring) == residue

    def test_parse_residue_ValueError(self):
        with pytest.raises(ValueError):
            util.parse_residue('ZZZ')

    @pytest.mark.parametrize('resname3, resname1', convert_aa_code_long_data())
    def test_convert_aa_3to1(self, resname3, resname1):
        assert util.convert_aa_code(resname3) == resname1

    @pytest.mark.parametrize('resname1, strings', aa)
    def test_convert_aa_1to3(self, resname1, strings):
        assert util.convert_aa_code(resname1) == strings[0]

    @pytest.mark.parametrize('x', (
        'XYZXYZ',
        '£'
    ))
    def test_ValueError(self, x):
        with pytest.raises(ValueError):
            util.convert_aa_code(x)


def test_greedy_splitext(inp="foo/bar/boing.2.pdb.bz2",
                         ref=["foo/bar/boing", ".2.pdb.bz2"]):
    inp = os.path.normpath(inp)
    ref[0] = os.path.normpath(ref[0])
    ref[1] = os.path.normpath(ref[1])
    root, ext = util.greedy_splitext(inp)
    assert root == ref[0], "root incorrect"
    assert ext == ref[1], "extension incorrect"


@pytest.mark.parametrize('iterable, value', [
    ([1, 2, 3], True),
    ([], True),
    ((1, 2, 3), True),
    ((), True),
    (range(3), True),
    (np.array([1, 2, 3]), True),
    (123, False),
    ("byte string", False),
    (u"unicode string", False)
])
def test_iterable(iterable, value):
    assert util.iterable(iterable) == value


class TestFilename(object):
    root = "foo"
    filename = "foo.psf"
    ext = "pdb"
    filename2 = "foo.pdb"

    @pytest.mark.parametrize('name, ext, keep, actual_name', [
        (filename, None, False, filename),
        (filename, ext, False, filename2),
        (filename, ext, True, filename),
        (root, ext, False, filename2),
        (root, ext, True, filename2)
    ])
    def test_string(self, name, ext, keep, actual_name):
        file_name = util.filename(name, ext, keep)
        assert file_name == actual_name

    def test_named_stream(self):
        ns = util.NamedStream(StringIO(), self.filename)
        fn = util.filename(ns, ext=self.ext)
        # assert_equal replace by this if loop to avoid segfault on some systems
        if fn != ns:
            pytest.fail("fn and ns are different")
        assert str(fn) == self.filename2
        assert ns.name == self.filename2


class TestGeometryFunctions(object):
    e1, e2, e3 = np.eye(3)
    a = np.array([np.cos(np.pi / 3), np.sin(np.pi / 3), 0])
    null = np.zeros(3)

    @pytest.mark.parametrize('x_axis, y_axis, value', [
        # Unit vectors
        (e1, e2, np.pi / 2),
        (e1, a, np.pi / 3),
        # Angle vectors
        (2 * e1, e2, np.pi / 2),
        (-2 * e1, e2, np.pi - np.pi / 2),
        (23.3 * e1, a, np.pi / 3),
        # Null vector
        (e1, null, np.nan),
        # Coleniar
        (a, a, 0.0)
    ])
    def test_vectors(self, x_axis, y_axis, value):
        assert_allclose(mdamath.angle(x_axis, y_axis), value)

    @pytest.mark.parametrize('x_axis, y_axis, value', [
        (-2.3456e7 * e1, 3.4567e-6 * e1, np.pi),
        (2.3456e7 * e1, 3.4567e-6 * e1, 0.0)
    ])
    def test_angle_pi(self, x_axis, y_axis, value):
        assert_almost_equal(mdamath.angle(x_axis, y_axis), value)

    @pytest.mark.parametrize('x', np.linspace(0, np.pi, 20))
    def test_angle_range(self, x):
        r = 1000.
        v = r * np.array([np.cos(x), np.sin(x), 0])
        assert_almost_equal(mdamath.angle(self.e1, v), x, 6)

    @pytest.mark.parametrize('vector, value', [
        (e3, 1),
        (a, np.linalg.norm(a)),
        (null, 0.0)
    ])
    def test_norm(self, vector, value):
        assert mdamath.norm(vector) == value

    @pytest.mark.parametrize('x', np.linspace(0, np.pi, 20))
    def test_norm_range(self, x):
        r = 1000.
        v = r * np.array([np.cos(x), np.sin(x), 0])
        assert_almost_equal(mdamath.norm(v), r, 6)

    @pytest.mark.parametrize('vec1, vec2, value', [
        (e1, e2, e3),
        (e1, null, 0.0)
    ])
    def test_normal(self, vec1, vec2, value):
        assert_allclose(mdamath.normal(vec1, vec2), value)
        # add more non-trivial tests

    def test_angle_lower_clip(self):
        a = np.array([0.1, 0, 0.2])
        x = np.dot(a**0.5, -(a**0.5)) / \
            (mdamath.norm(a**0.5) * mdamath.norm(-(a**0.5)))
        assert x < -1.0
        assert mdamath.angle(a, -(a)) == np.pi
        assert mdamath.angle(a**0.5, -(a**0.5)) == np.pi

    def test_stp(self):
        assert mdamath.stp(self.e1, self.e2, self.e3) == 1.0
        # add more non-trivial tests

    def test_dihedral(self):
        ab = self.e1
        bc = ab + self.e2
        cd = bc + self.e3
        assert_almost_equal(mdamath.dihedral(ab, bc, cd), -np.pi / 2)

    def test_pdot(self):
        arr = np.random.rand(4, 3)
        matrix_dot = mdamath.pdot(arr, arr)
        list_dot = [np.dot(a, a) for a in arr]
        assert_almost_equal(matrix_dot, list_dot)

    def test_pnorm(self):
        arr = np.random.rand(4, 3)
        matrix_norm = mdamath.pnorm(arr)
        list_norm = [np.linalg.norm(a) for a in arr]
        assert_almost_equal(matrix_norm, list_norm)


class TestMatrixOperations(object):

    def ref_trivecs(self, box):
        box = np.asarray(box, dtype=np.float64)
        x, y, z, a, b, c = box
        # Only positive edge lengths and angles in (0, 180) are allowed:
        if np.any(box <= 0) or a >= 180 or b >= 180 or c >= 180:
            ref = np.zeros((3, 3), dtype=np.float32)
        # detect orthogonal boxes:
        elif a == 90 and b == 90 and c == 90:
            ref = np.diag(box[:3].astype(np.float32))
        else:
            ref = np.zeros((3, 3), dtype=np.float64)
            cos_a = 0.0 if a == 90 else np.cos(np.deg2rad(a))
            cos_b = 0.0 if b == 90 else np.cos(np.deg2rad(b))
            cos_c = 0.0 if c == 90 else np.cos(np.deg2rad(c))
            sin_c = 1.0 if c == 90 else np.sin(np.deg2rad(c))
            ref[0, 0] = x
            ref[1, 0] = y * cos_c
            ref[1, 1] = y * sin_c
            ref[2, 0] = z * cos_b
            ref[2, 1] = z * (cos_a - cos_b * cos_c) / sin_c
            ref[2, 2] = np.sqrt(z * z - ref[2, 0] ** 2 - ref[2, 1] ** 2)
            if ref[2, 2] == 0 or np.isnan(ref[2, 2]):
                ref[:, :] = 0.0
            ref = ref.astype(np.float32)
        return ref

    def ref_trivecs_unsafe(self, box):
        box = np.asarray(box, dtype=np.float64)
        x, y, z, a, b, c = box
        # detect orthogonal boxes:
        if a == 90 and b == 90 and c == 90:
            ref = np.diag(box[:3].astype(np.float32))
        else:
            ref = np.zeros((3, 3), dtype=np.float64)
            cos_a = 0.0 if a == 90 else np.cos(np.deg2rad(a))
            cos_b = 0.0 if b == 90 else np.cos(np.deg2rad(b))
            cos_c = 0.0 if c == 90 else np.cos(np.deg2rad(c))
            sin_c = 1.0 if c == 90 else np.sin(np.deg2rad(c))
            ref[0, 0] = x
            ref[1, 0] = y * cos_c
            ref[1, 1] = y * sin_c
            ref[2, 0] = z * cos_b
            ref[2, 1] = z * (cos_a - cos_b * cos_c) / sin_c
            with np.errstate(invalid="ignore"):
                ref[2, 2] = np.sqrt(z * z - ref[2, 0] ** 2 - ref[2, 1] ** 2)
            ref = ref.astype(np.float32)
        return ref

    def ref_tribox(self, tri_vecs):
        tri_vecs = tri_vecs.astype(np.float64)
        x, y, z = np.linalg.norm(tri_vecs, axis=1)
        a = np.rad2deg(np.arccos(np.dot(tri_vecs[1], tri_vecs[2]) / (y * z)))
        b = np.rad2deg(np.arccos(np.dot(tri_vecs[0], tri_vecs[2]) / (x * z)))
        c = np.rad2deg(np.arccos(np.dot(tri_vecs[0], tri_vecs[1]) / (x * y)))
        box = np.array([x, y, z, a, b, c], dtype=np.float32)
        if not (np.all(box > 0) and a < 180 and b < 180 and c < 180):
            box = np.zeros(6, dtype=np.float32)
        return box

    @pytest.mark.parametrize('lengths', comb_wr([-1, 0, 1, 2], 3))
    @pytest.mark.parametrize('angles',
                             comb_wr([-10, 0, 20, 70, 90, 120, 180], 3))
    def test_triclinic_vectors(self, lengths, angles):
        box = lengths + angles
        ref = self.ref_trivecs(box)
        res = mdamath.triclinic_vectors(box)
        assert_array_equal(res, ref)
        # check for default dtype:
        assert res.dtype == np.float32
        # belts and braces, make sure upper triangle is always zero:
        assert not(res[0, 1] or res[0, 2] or res[1, 2])

    @pytest.mark.parametrize('alpha', (60, 90))
    @pytest.mark.parametrize('beta', (60, 90))
    @pytest.mark.parametrize('gamma', (60, 90))
    def test_triclinic_vectors_right_angle_zeros(self, alpha, beta, gamma):
        angles = [alpha, beta, gamma]
        box = [10, 20, 30] + angles
        mat = mdamath.triclinic_vectors(box)
        if 90 in angles:
            if gamma == 90:
                assert not mat[1, 0]
                if alpha == 90:
                    assert not mat[2, 1]
                    if beta == 90:
                        assert not mat[2, 0]
                    else:
                        assert mat[2, 0]
                else:
                    assert mat[2, 1]
            else:
                assert mat[1, 0]
                if beta == 90:
                    assert not mat[2, 0]
                    if alpha == 90:
                        assert not mat[2, 1]
                    else:
                        assert mat[2, 1]
                else:
                    assert mat[2, 0]
                    # 2, 1 cannot be zero here regardless of alpha
                    assert mat[2, 1]
        else:
            assert mat[1, 0] and mat[2, 0] and mat[2, 1]

    @pytest.mark.parametrize('dtype', (int, float, np.float32, np.float64))
    def test_triclinic_vectors_retval(self, dtype):
        # valid box
        box = [1, 1, 1, 70, 80, 90]
        res = mdamath.triclinic_vectors(box, dtype=dtype)
        assert res.shape == (3, 3)
        assert res.dtype == dtype
        # zero box
        box = [0, 0, 0, 0, 0, 0]
        res = mdamath.triclinic_vectors(box, dtype=dtype)
        assert res.shape == (3, 3)
        assert res.dtype == dtype
        assert np.all(res == 0)
        # invalid box angles
        box = [1, 1, 1, 40, 40, 90]
        res = mdamath.triclinic_vectors(box, dtype=dtype)
        assert res.shape == (3, 3)
        assert res.dtype == dtype
        assert np.all(res == 0)
        # invalid box lengths:
        box = [-1, 1, 1, 70, 80, 90]
        res = mdamath.triclinic_vectors(box, dtype=dtype)
        assert res.shape == (3, 3)
        assert res.dtype == dtype
        assert np.all(res == 0)

    def test_triclinic_vectors_box_cycle(self):
        max_error = 0.0
        for a in range(10, 91, 10):
            for b in range(10, 91, 10):
                for g in range(10, 91, 10):
                    ref = np.array([1, 1, 1, a, b, g], dtype=np.float32)
                    res = mdamath.triclinic_box(
                        *mdamath.triclinic_vectors(ref))
                    if not np.all(res == 0.0):
                        assert_almost_equal(res, ref, 5)

    @pytest.mark.parametrize('angles', ([70, 70, 70],
                                        [70, 70, 90],
                                        [70, 90, 70],
                                        [90, 70, 70],
                                        [70, 90, 90],
                                        [90, 70, 90],
                                        [90, 90, 70]))
    def test_triclinic_vectors_box_cycle_exact(self, angles):
        # These cycles were inexact prior to PR #2201
        ref = np.array([10.1, 10.1, 10.1] + angles, dtype=np.float32)
        res = mdamath.triclinic_box(*mdamath.triclinic_vectors(ref))
        assert_allclose(res, ref)

    @pytest.mark.parametrize('lengths', comb_wr([-1, 0, 1, 2], 3))
    @pytest.mark.parametrize('angles',
                             comb_wr([-10, 0, 20, 70, 90, 120, 180], 3))
    def test_triclinic_box(self, lengths, angles):
        tri_vecs = self.ref_trivecs_unsafe(lengths + angles)
        ref = self.ref_tribox(tri_vecs)
        res = mdamath.triclinic_box(*tri_vecs)
        assert_array_equal(res, ref)
        assert res.dtype == ref.dtype

    @pytest.mark.parametrize('lengths', comb_wr([-1, 0, 1, 2], 3))
    @pytest.mark.parametrize('angles',
                             comb_wr([-10, 0, 20, 70, 90, 120, 180], 3))
    def test_box_volume(self, lengths, angles):
        box = np.array(lengths + angles, dtype=np.float32)
        assert_almost_equal(mdamath.box_volume(box),
                            np.linalg.det(self.ref_trivecs(box)),
                            decimal=5)

    def test_sarrus_det(self):
        comb = comb_wr(np.linspace(-133.7, 133.7, num=5), 9)
        # test array of matrices:
        matrix = np.array(tuple(comb)).reshape((-1, 5, 3, 3))
        ref = np.linalg.det(matrix)
        res = mdamath.sarrus_det(matrix)
        assert_almost_equal(res, ref, 7)
        assert ref.dtype == res.dtype == np.float64
        # test single matrices:
        matrix = matrix.reshape(-1, 3, 3)
        ref = ref.ravel()
        res = np.array([mdamath.sarrus_det(m) for m in matrix])
        assert_almost_equal(res, ref, 7)
        assert ref.dtype == res.dtype == np.float64

    @pytest.mark.parametrize('shape', ((0,), (3, 2), (2, 3), (1, 1, 3, 1)))
    def test_sarrus_det_wrong_shape(self, shape):
        matrix = np.zeros(shape)
        with pytest.raises(ValueError):
            mdamath.sarrus_det(matrix)


class TestMakeWhole(object):
    """Set up a simple system:

    +-----------+
    |           |
    | 6       3 | 6
    | !       ! | !
    |-5-8   1-2-|-5-8
    | !       ! | !
    | 7       4 | 7
    |           |
    +-----------+
    """

    prec = 5

    @pytest.fixture()
    def universe(self):
        universe = mda.Universe(Make_Whole)
        bondlist = [(0, 1), (1, 2), (1, 3), (1, 4), (4, 5), (4, 6), (4, 7)]
        universe.add_TopologyAttr(Bonds(bondlist))
        return universe

    def test_return_value(self, universe):
        ag = universe.residues[0].atoms
        orig_pos = ag.positions.copy()
        retval = mdamath.make_whole(ag)
        assert retval.dtype == np.float32
        assert_array_equal(ag.positions, retval)
        assert np.any(ag.positions != orig_pos)

    def test_single_atom_no_bonds(self):
        # Call make_whole on single atom with no bonds, shouldn't move
        u = mda.Universe(Make_Whole)
        # Atom0 is isolated
        bondlist = [(1, 2), (1, 3), (1, 4), (4, 5), (4, 6), (4, 7)]
        u.add_TopologyAttr(Bonds(bondlist))

        ag = u.atoms[[0]]
        refpos = ag.positions.copy()
        mdamath.make_whole(ag)

        assert_array_equal(ag.positions, refpos)  # must be untouched

    def test_empty_ag(self, universe):
        ag = mda.AtomGroup([], universe)
        retval = mdamath.make_whole(ag)
        assert retval.dtype == np.float32
        assert_array_equal(retval, np.empty((0, 3), dtype=np.float32))

    def test_scrambled_ag(self, universe):
        # if order of atomgroup is mixed
        ag = universe.atoms[[1, 3, 2, 4, 0, 6, 5, 7]]

        mdamath.make_whole(ag)

        # artificial system which uses 1nm bonds, so
        # largest bond should be 20A
        assert ag.bonds.values().max() < 20.1

    def test_out_of_place(self, universe):
        ag = universe.residues[0].atoms
        orig_pos = ag.positions.copy()
        mdamath.make_whole(ag, inplace=False)
        # positions must be untouched:
        assert_array_equal(ag.positions, orig_pos)

    def test_double_precision_box(self):
        # This test could in principle be removed since PR #2213
        # universe with double precision box containing a 2-atom molecule
        # broken accross a corner:
        u = mda.Universe.empty(
            n_atoms=2,
            n_residues=1,
            n_segments=1,
            atom_resindex=[0, 0],
            residue_segindex=[0],
            trajectory=True,
            velocities=False,
            forces=False)
        ts = u.trajectory.ts
        ts.frame = 0
        ts.dimensions = [10, 10, 10, 90, 90, 90]
        # assert ts.dimensions.dtype == np.float64
        # not applicable since #2213
        ts.positions = np.array([[1, 1, 1, ], [9, 9, 9]], dtype=np.float32)
        u.add_TopologyAttr(Bonds([(0, 1)]))
        mdamath.make_whole(u.atoms)
        assert_array_almost_equal(u.atoms.positions,
                                  np.array([[1, 1, 1, ], [-1, -1, -1]],
                                           dtype=np.float32))

    @staticmethod
    @pytest.fixture()
    def ag(universe):
        return universe.residues[0].atoms

    def test_no_bonds(self):
        # NoData caused by no bonds
        universe = mda.Universe(Make_Whole)
        ag = universe.residues[0].atoms
        with pytest.raises(NoDataError):
            mdamath.make_whole(ag)

    def test_zero_box_size(self, universe, ag):
        universe.dimensions = [0., 0., 0., 90., 90., 90.]
        with pytest.raises(ValueError):
            mdamath.make_whole(ag)

    def test_wrong_reference_atom(self, universe, ag):
        # Reference atom not in atomgroup
        with pytest.raises(ValueError):
            mdamath.make_whole(ag, reference_atom=universe.atoms[-1])

    def test_impossible_solve(self, universe):
        # check that the algorithm sees the bad walk
        with pytest.raises(ValueError):
            mdamath.make_whole(universe.atoms)

    def test_solve_1(self, universe, ag):
        # regular usage of function

        refpos = universe.atoms[:4].positions.copy()

        mdamath.make_whole(ag)

        assert_array_almost_equal(universe.atoms[:4].positions, refpos)
        assert_array_almost_equal(universe.atoms[4].position,
                                  np.array([110.0, 50.0, 0.0]), decimal=self.prec)
        assert_array_almost_equal(universe.atoms[5].position,
                                  np.array([110.0, 60.0, 0.0]), decimal=self.prec)
        assert_array_almost_equal(universe.atoms[6].position,
                                  np.array([110.0, 40.0, 0.0]), decimal=self.prec)
        assert_array_almost_equal(universe.atoms[7].position,
                                  np.array([120.0, 50.0, 0.0]), decimal=self.prec)

    def test_solve_2(self, universe, ag):
        # use but specify the center atom

        refpos = universe.atoms[4:8].positions.copy()

        mdamath.make_whole(ag, reference_atom=universe.residues[0].atoms[4])

        assert_array_almost_equal(universe.atoms[4:8].positions, refpos)
        assert_array_almost_equal(universe.atoms[0].position,
                                  np.array([-20.0, 50.0, 0.0]), decimal=self.prec)
        assert_array_almost_equal(universe.atoms[1].position,
                                  np.array([-10.0, 50.0, 0.0]), decimal=self.prec)
        assert_array_almost_equal(universe.atoms[2].position,
                                  np.array([-10.0, 60.0, 0.0]), decimal=self.prec)
        assert_array_almost_equal(universe.atoms[3].position,
                                  np.array([-10.0, 40.0, 0.0]), decimal=self.prec)

    def test_solve_3(self, universe):
        # put in a chunk that doesn't need any work

        refpos = universe.atoms[:1].positions.copy()

        mdamath.make_whole(universe.atoms[:1])

        assert_array_almost_equal(universe.atoms[:1].positions, refpos)

    def test_solve_4(self, universe):
        # Put in only some of a fragment,
        # check that not everything gets moved

        chunk = universe.atoms[:7]
        refpos = universe.atoms[7].position.copy()

        mdamath.make_whole(chunk)

        assert_array_almost_equal(universe.atoms[7].position, refpos)
        assert_array_almost_equal(universe.atoms[4].position,
                                  np.array([110.0, 50.0, 0.0]))
        assert_array_almost_equal(universe.atoms[5].position,
                                  np.array([110.0, 60.0, 0.0]))
        assert_array_almost_equal(universe.atoms[6].position,
                                  np.array([110.0, 40.0, 0.0]))

    def test_double_frag_short_bonds(self, universe, ag):
        # previous bug where if two fragments are given
        # but all bonds were short, the algorithm didn't
        # complain
        mdamath.make_whole(ag)
        with pytest.raises(ValueError):
            mdamath.make_whole(universe.atoms)

    def test_make_whole_triclinic(self):
        u = mda.Universe(TPR, GRO)
        thing = u.select_atoms('not resname SOL NA+')
        mdamath.make_whole(thing)

        blengths = thing.bonds.values()

        assert blengths.max() < 2.0

    def test_make_whole_fullerene(self):
        # lots of circular bonds as a nice pathological case
        u = mda.Universe(fullerene)

        bbox = u.atoms.bbox()
        u.dimensions = np.r_[bbox[1] - bbox[0], [90]*3]

        blengths = u.atoms.bonds.values()
        # kaboom
        u.atoms[::2].translate([u.dimensions[0], -2 * u.dimensions[1], 0.0])
        u.atoms[1::2].translate(
            [0.0, 7 * u.dimensions[1], -5 * u.dimensions[2]])

        mdamath.make_whole(u.atoms)

        assert_array_almost_equal(
            u.atoms.bonds.values(), blengths, decimal=self.prec)

    def test_make_whole_multiple_molecules(self):
        u = mda.Universe(two_water_gro, guess_bonds=True)

        for f in u.atoms.fragments:
            mdamath.make_whole(f)

        assert u.atoms.bonds.values().max() < 2.0


class Class_with_Caches(object):
    def __init__(self):
        self._cache = dict()
        self.ref1 = 1.0
        self.ref2 = 2.0
        self.ref3 = 3.0
        self.ref4 = 4.0
        self.ref5 = 5.0
        self.ref6 = 6.0
        # For universe-validated caches
        # One-line lambda-like class
        self.universe = type('Universe', (), dict())()
        self.universe._cache = {'_valid': {}}

    @cached('val1')
    def val1(self):
        return self.ref1

    # Do one with property decorator as these are used together often
    @property
    @cached('val2')
    def val2(self):
        return self.ref2

    # Check use of property setters
    @property
    @cached('val3')
    def val3(self):
        return self.ref3

    @val3.setter
    def val3(self, new):
        self._clear_caches('val3')
        self._fill_cache('val3', new)

    @val3.deleter
    def val3(self):
        self._clear_caches('val3')

    # Check that args are passed through to underlying functions
    @cached('val4')
    def val4(self, n1, n2):
        return self._init_val_4(n1, n2)

    def _init_val_4(self, m1, m2):
        return self.ref4 + m1 + m2

    # Args and Kwargs
    @cached('val5')
    def val5(self, n, s=None):
        return self._init_val_5(n, s=s)

    def _init_val_5(self, n, s=None):
        return n * s

    # Property decorator and universally-validated cache
    @property
    @cached('val6', universe_validation=True)
    def val6(self):
        return self.ref5 + 1.0

    # These are designed to mimic the AG and Universe cache methods
    def _clear_caches(self, *args):
        if len(args) == 0:
            self._cache = dict()
        else:
            for name in args:
                try:
                    del self._cache[name]
                except KeyError:
                    pass

    def _fill_cache(self, name, value):
        self._cache[name] = value


class TestCachedDecorator(object):
    @pytest.fixture()
    def obj(self):
        return Class_with_Caches()

    def test_val1_lookup(self, obj):
        obj._clear_caches()
        assert 'val1' not in obj._cache
        assert obj.val1() == obj.ref1
        ret = obj.val1()
        assert 'val1' in obj._cache
        assert obj._cache['val1'] == ret
        assert obj.val1() is obj._cache['val1']

    def test_val1_inject(self, obj):
        # Put something else into the cache and check it gets returned
        # this tests that the cache is blindly being used
        obj._clear_caches()
        ret = obj.val1()
        assert 'val1' in obj._cache
        assert ret == obj.ref1
        new = 77.0
        obj._fill_cache('val1', new)
        assert obj.val1() == new

    # Managed property
    def test_val2_lookup(self, obj):
        obj._clear_caches()
        assert 'val2' not in obj._cache
        assert obj.val2 == obj.ref2
        ret = obj.val2
        assert 'val2' in obj._cache
        assert obj._cache['val2'] == ret

    def test_val2_inject(self, obj):
        obj._clear_caches()
        ret = obj.val2
        assert 'val2' in obj._cache
        assert ret == obj.ref2
        new = 77.0
        obj._fill_cache('val2', new)
        assert obj.val2 == new

        # Setter on cached attribute

    def test_val3_set(self, obj):
        obj._clear_caches()
        assert obj.val3 == obj.ref3
        new = 99.0
        obj.val3 = new
        assert obj.val3 == new
        assert obj._cache['val3'] == new

    def test_val3_del(self, obj):
        # Check that deleting the property removes it from cache,
        obj._clear_caches()
        assert obj.val3 == obj.ref3
        assert 'val3' in obj._cache
        del obj.val3
        assert 'val3' not in obj._cache
        # But allows it to work as usual afterwards
        assert obj.val3 == obj.ref3
        assert 'val3' in obj._cache

    # Pass args
    def test_val4_args(self, obj):
        obj._clear_caches()
        assert obj.val4(1, 2) == 1 + 2 + obj.ref4
        # Further calls should yield the old result
        # this arguably shouldn't be cached...
        assert obj.val4(3, 4) == 1 + 2 + obj.ref4

    # Pass args and kwargs
    def test_val5_kwargs(self, obj):
        obj._clear_caches()
        assert obj.val5(5, s='abc') == 5 * 'abc'

        assert obj.val5(5, s='!!!') == 5 * 'abc'

    # property decorator, with universe validation
    def test_val6_universe_validation(self, obj):
        obj._clear_caches()
        assert not hasattr(obj, '_cache_key')
        assert 'val6' not in obj._cache
        assert 'val6' not in obj.universe._cache['_valid']

        ret = obj.val6  # Trigger caching
        assert obj.val6 == obj.ref6
        assert ret is obj.val6
        assert 'val6' in obj._cache
        assert 'val6' in obj.universe._cache['_valid']
        assert obj._cache_key in obj.universe._cache['_valid']['val6']
        assert obj._cache['val6'] is ret

        # Invalidate cache at universe level
        obj.universe._cache['_valid']['val6'].clear()
        ret2 = obj.val6
        assert ret2 is obj.val6
        assert ret2 is not ret

        # Clear obj cache and access again
        obj._clear_caches()
        ret3 = obj.val6
        assert ret3 is obj.val6
        assert ret3 is not ret2
        assert ret3 is not ret


class TestConvFloat(object):
    @pytest.mark.parametrize('s, output', [
        ('0.45', 0.45),
        ('.45', 0.45),
        ('a.b', 'a.b')
    ])
    def test_float(self, s, output):
        assert util.conv_float(s) == output

    @pytest.mark.parametrize('input, output', [
        (('0.45', '0.56', '6.7'), [0.45, 0.56, 6.7]),
        (('0.45', 'a.b', '!!'), [0.45, 'a.b', '!!'])
    ])
    def test_map(self, input, output):
        ret = [util.conv_float(el) for el in input]
        assert ret == output


class TestFixedwidthBins(object):
    def test_keys(self):
        ret = util.fixedwidth_bins(0.5, 1.0, 2.0)
        for k in ['Nbins', 'delta', 'min', 'max']:
            assert k in ret

    def test_ValueError(self):
        with pytest.raises(ValueError):
            util.fixedwidth_bins(0.1, 5.0, 4.0)

    @pytest.mark.parametrize(
        'delta, xmin, xmax, output_Nbins, output_delta, output_min, output_max',
        [
            (0.1, 4.0, 5.0, 10, 0.1, 4.0, 5.0),
            (0.4, 4.0, 5.0, 3, 0.4, 3.9, 5.1)
        ])
    def test_usage(self, delta, xmin, xmax, output_Nbins, output_delta,
                   output_min, output_max):
        ret = util.fixedwidth_bins(delta, xmin, xmax)
        assert ret['Nbins'] == output_Nbins
        assert ret['delta'] == output_delta
        assert ret['min'], output_min
        assert ret['max'], output_max


@pytest.fixture
def atoms():
    from MDAnalysisTests import make_Universe
    u = make_Universe(extras=("masses",), size=(3, 1, 1))
    return u.atoms


@pytest.mark.parametrize('weights,result',
                         [
                             (None, None),
                             ("mass", np.array([5.1, 4.2, 3.3])),
                             (np.array([12.0, 1.0, 12.0]),
                              np.array([12.0, 1.0, 12.0])),
                             ([12.0, 1.0, 12.0], np.array([12.0, 1.0, 12.0])),
                             (range(3), np.arange(3, dtype=int)),
                         ])
def test_check_weights_ok(atoms, weights, result):
    assert_array_equal(util.get_weights(atoms, weights), result)


@pytest.mark.parametrize('weights',
                         [42,
                          "geometry",
                          np.array(1.0),
                          np.array([12.0, 1.0, 12.0, 1.0]),
                          [12.0, 1.0],
                          np.array([[12.0, 1.0, 12.0]]),
                          np.array([[12.0, 1.0, 12.0], [12.0, 1.0, 12.0]]),
                          ])
def test_check_weights_raises_ValueError(atoms, weights):
    with pytest.raises(ValueError):
        util.get_weights(atoms, weights)


class TestGuessFormat(object):
    """Test guessing of format from filenames

    Tests also getting the appropriate Parser and Reader from a
    given filename
    """
    # list of known formats, followed by the desired Parser and Reader
    # None indicates that there isn't a Reader for this format
    # All formats call fallback to the MinimalParser
    formats = [
        ('CHAIN', mda.topology.MinimalParser.MinimalParser,
         mda.coordinates.chain.ChainReader),
        ('CONFIG', mda.topology.DLPolyParser.ConfigParser,
         mda.coordinates.DLPoly.ConfigReader),
        ('CRD', mda.topology.CRDParser.CRDParser, mda.coordinates.CRD.CRDReader),
        ('DATA', mda.topology.LAMMPSParser.DATAParser,
         mda.coordinates.LAMMPS.DATAReader),
        ('DCD', mda.topology.MinimalParser.MinimalParser,
         mda.coordinates.DCD.DCDReader),
        ('DMS', mda.topology.DMSParser.DMSParser, mda.coordinates.DMS.DMSReader),
        ('GMS', mda.topology.GMSParser.GMSParser, mda.coordinates.GMS.GMSReader),
        ('GRO', mda.topology.GROParser.GROParser, mda.coordinates.GRO.GROReader),
        ('HISTORY', mda.topology.DLPolyParser.HistoryParser,
         mda.coordinates.DLPoly.HistoryReader),
        ('INPCRD', mda.topology.MinimalParser.MinimalParser,
         mda.coordinates.INPCRD.INPReader),
        ('LAMMPS', mda.topology.MinimalParser.MinimalParser,
         mda.coordinates.LAMMPS.DCDReader),
        ('MDCRD', mda.topology.MinimalParser.MinimalParser,
         mda.coordinates.TRJ.TRJReader),
        ('MMTF', mda.topology.MMTFParser.MMTFParser,
         mda.coordinates.MMTF.MMTFReader),
        ('MOL2', mda.topology.MOL2Parser.MOL2Parser,
         mda.coordinates.MOL2.MOL2Reader),
        ('NC', mda.topology.MinimalParser.MinimalParser,
         mda.coordinates.TRJ.NCDFReader),
        ('NCDF', mda.topology.MinimalParser.MinimalParser,
         mda.coordinates.TRJ.NCDFReader),
        ('PDB', mda.topology.PDBParser.PDBParser, mda.coordinates.PDB.PDBReader),
        ('PDBQT', mda.topology.PDBQTParser.PDBQTParser,
         mda.coordinates.PDBQT.PDBQTReader),
        ('PRMTOP', mda.topology.TOPParser.TOPParser, None),
        ('PQR', mda.topology.PQRParser.PQRParser, mda.coordinates.PQR.PQRReader),
        ('PSF', mda.topology.PSFParser.PSFParser, None),
        ('RESTRT', mda.topology.MinimalParser.MinimalParser,
         mda.coordinates.INPCRD.INPReader),
        ('TOP', mda.topology.TOPParser.TOPParser, None),
        ('TPR', mda.topology.TPRParser.TPRParser, None),
        ('TRJ', mda.topology.MinimalParser.MinimalParser,
         mda.coordinates.TRJ.TRJReader),
        ('TRR', mda.topology.MinimalParser.MinimalParser,
         mda.coordinates.TRR.TRRReader),
        ('XML', mda.topology.HoomdXMLParser.HoomdXMLParser, None),
        ('XPDB', mda.topology.ExtendedPDBParser.ExtendedPDBParser,
         mda.coordinates.PDB.ExtendedPDBReader),
        ('XTC', mda.topology.MinimalParser.MinimalParser,
         mda.coordinates.XTC.XTCReader),
        ('XYZ', mda.topology.XYZParser.XYZParser, mda.coordinates.XYZ.XYZReader),
        ('TRZ', mda.topology.MinimalParser.MinimalParser,
         mda.coordinates.TRZ.TRZReader),
    ]
    # list of possible compressed extensions
    # include no extension too!
    compressed_extensions = ['.bz2', '.gz']

    @pytest.mark.parametrize('extention',
                             [format_tuple[0].upper() for format_tuple in
                              formats] +
                             [format_tuple[0].lower() for format_tuple in
                              formats])
    def test_get_extention(self, extention):
        """Check that get_ext works"""
        file_name = 'file.{0}'.format(extention)
        a, b = util.get_ext(file_name)

        assert a == 'file'
        assert b == extention.lower()

    @pytest.mark.parametrize('extention',
                             [format_tuple[0].upper() for format_tuple in
                              formats] +
                             [format_tuple[0].lower() for format_tuple in
                              formats])
    def test_compressed_without_compression_extention(self, extention):
        """Check that format suffixed by compressed extension works"""
        file_name = 'file.{0}'.format(extention)
        a = util.format_from_filename_extension(file_name)
        # expect answer to always be uppercase
        assert a == extention.upper()

    @pytest.mark.parametrize('extention',
                             [format_tuple[0].upper() for format_tuple in
                              formats] +
                             [format_tuple[0].lower() for format_tuple in
                              formats])
    @pytest.mark.parametrize('compression_extention', compressed_extensions)
    def test_compressed(self, extention, compression_extention):
        """Check that format suffixed by compressed extension works"""
        file_name = 'file.{0}{1}'.format(extention, compression_extention)
        a = util.format_from_filename_extension(file_name)
        # expect answer to always be uppercase
        assert a == extention.upper()

    @pytest.mark.parametrize('extention',
                             [format_tuple[0].upper() for format_tuple in
                              formats] + [format_tuple[0].lower() for
                                          format_tuple in formats])
    def test_guess_format(self, extention):
        file_name = 'file.{0}'.format(extention)
        a = util.guess_format(file_name)
        # expect answer to always be uppercase
        assert a == extention.upper()

    @pytest.mark.parametrize('extention',
                             [format_tuple[0].upper() for format_tuple in
                              formats] + [format_tuple[0].lower() for
                                          format_tuple in formats])
    @pytest.mark.parametrize('compression_extention', compressed_extensions)
    def test_guess_format_compressed(self, extention, compression_extention):
        file_name = 'file.{0}{1}'.format(extention, compression_extention)
        a = util.guess_format(file_name)
        # expect answer to always be uppercase
        assert a == extention.upper()

    @pytest.mark.parametrize('extention, parser',
                             [(format_tuple[0], format_tuple[1]) for
                              format_tuple in formats if
                              format_tuple[1] is not None]
                             )
    def test_get_parser(self, extention, parser):
        file_name = 'file.{0}'.format(extention)
        a = get_parser_for(file_name)

        assert a == parser

    @pytest.mark.parametrize('extention, parser',
                             [(format_tuple[0], format_tuple[1]) for
                              format_tuple in formats if
                              format_tuple[1] is not None]
                             )
    @pytest.mark.parametrize('compression_extention', compressed_extensions)
    def test_get_parser_compressed(self, extention, parser,
                                   compression_extention):
        file_name = 'file.{0}{1}'.format(extention, compression_extention)
        a = get_parser_for(file_name)

        assert a == parser

    @pytest.mark.parametrize('extention',
                             [(format_tuple[0], format_tuple[1]) for
                              format_tuple in formats if
                              format_tuple[1] is None]
                             )
    def test_get_parser_invalid(self, extention):
        file_name = 'file.{0}'.format(extention)
        with pytest.raises(ValueError):
            get_parser_for(file_name)

    @pytest.mark.parametrize('extention, reader',
                             [(format_tuple[0], format_tuple[2]) for
                              format_tuple in formats if
                              format_tuple[2] is not None]
                             )
    def test_get_reader(self, extention, reader):
        file_name = 'file.{0}'.format(extention)
        a = mda.coordinates.core.get_reader_for(file_name)

        assert a == reader

    @pytest.mark.parametrize('extention, reader',
                             [(format_tuple[0], format_tuple[2]) for
                              format_tuple in formats if
                              format_tuple[2] is not None]
                             )
    @pytest.mark.parametrize('compression_extention', compressed_extensions)
    def test_get_reader_compressed(self, extention, reader,
                                   compression_extention):
        file_name = 'file.{0}{1}'.format(extention, compression_extention)
        a = mda.coordinates.core.get_reader_for(file_name)

        assert a == reader

    @pytest.mark.parametrize('extention',
                             [(format_tuple[0], format_tuple[2]) for
                              format_tuple in formats if
                              format_tuple[2] is None]
                             )
    def test_get_reader_invalid(self, extention):
        file_name = 'file.{0}'.format(extention)
        with pytest.raises(ValueError):
            mda.coordinates.core.get_reader_for(file_name)

    def test_check_compressed_format_TypeError(self):
        with pytest.raises(TypeError):
            util.check_compressed_format(1234, 'bz2')

    def test_format_from_filename_TypeError(self):
        with pytest.raises(TypeError):
            util.format_from_filename_extension(1234)

    def test_guess_format_stream_ValueError(self):
        # This stream has no name, so can't guess format
        s = StringIO('this is a very fun file')
        with pytest.raises(ValueError):
            util.guess_format(s)

    def test_from_ndarray(self):
        fn = np.zeros((3, 3))
        rd = mda.coordinates.core.get_reader_for(fn)
        assert rd == mda.coordinates.memory.MemoryReader


class TestUniqueRows(object):
    def test_unique_rows_2(self):
        a = np.array([[0, 1], [1, 2], [2, 1], [0, 1], [0, 1], [2, 1]])

        assert_array_equal(util.unique_rows(a),
                           np.array([[0, 1], [1, 2], [2, 1]]))

    def test_unique_rows_3(self):
        a = np.array([[0, 1, 2], [0, 1, 2], [2, 3, 4], [0, 1, 2]])

        assert_array_equal(util.unique_rows(a),
                           np.array([[0, 1, 2], [2, 3, 4]]))

    def test_unique_rows_with_view(self):
        # unique_rows doesn't work when flags['OWNDATA'] is False,
        # happens when second dimension is created through broadcast
        a = np.array([1, 2])

        assert_array_equal(util.unique_rows(a[None, :]),
                           np.array([[1, 2]]))


class TestGetWriterFor(object):
    def test_no_filename_argument(self):
        # Does ``get_writer_for`` fails as expected when provided no
        # filename arguments
        with pytest.raises(TypeError):
            mda.coordinates.core.get_writer_for()

    def test_precedence(self):
        writer = mda.coordinates.core.get_writer_for('test.pdb', 'GRO')
        assert writer == mda.coordinates.GRO.GROWriter
        # Make sure ``get_writer_for`` uses *format* if provided

    def test_missing_extension(self):
        # Make sure ``get_writer_for`` behave as expected if *filename*
        # has no extension
        with pytest.raises(ValueError):
            mda.coordinates.core.get_writer_for(filename='test', format=None)

    def test_extension_empty_string(self):
        """
        Test format=''.

        Raises TypeError because format can be only None or
        valid formats.
        """
        with pytest.raises(ValueError):
            mda.coordinates.core.get_writer_for(filename='test', format='')

    def test_file_no_extension(self):
        """No format given"""
        with pytest.raises(ValueError):
            mda.coordinates.core.get_writer_for('outtraj')

    def test_wrong_format(self):
        # Make sure ``get_writer_for`` fails if the format is unknown
        with pytest.raises(TypeError):
            mda.coordinates.core.get_writer_for(filename="fail_me",
                                                format='UNK')

    def test_compressed_extension(self):
        for ext in ('.gz', '.bz2'):
            fn = 'test.gro' + ext
            writer = mda.coordinates.core.get_writer_for(filename=fn)
            assert writer == mda.coordinates.GRO.GROWriter
            # Make sure ``get_writer_for`` works with compressed file file names

    def test_compressed_extension_fail(self):
        for ext in ('.gz', '.bz2'):
            fn = 'test.unk' + ext
            # Make sure ``get_writer_for`` fails if an unknown format is compressed
            with pytest.raises(TypeError):
                mda.coordinates.core.get_writer_for(filename=fn)

    def test_non_string_filename(self):
        # Does ``get_writer_for`` fails with non string filename, no format
        with pytest.raises(ValueError):
            mda.coordinates.core.get_writer_for(filename=StringIO(),
                                                format=None)

    def test_multiframe_failure(self):
        # does ``get_writer_for`` fail with invalid format and multiframe not None
        with pytest.raises(TypeError):
            mda.coordinates.core.get_writer_for(filename="fail_me",
                                                format='UNK', multiframe=True)
            mda.coordinates.core.get_writer_for(filename="fail_me",
                                                format='UNK', multiframe=False)

    def test_multiframe_nonsense(self):
        with pytest.raises(ValueError):
            mda.coordinates.core.get_writer_for(filename='this.gro',
                                                multiframe='sandwich')

    formats = [
        # format name, related class, singleframe, multiframe
        ('CRD', mda.coordinates.CRD.CRDWriter, True, False),
        ('DATA', mda.coordinates.LAMMPS.DATAWriter, True, False),
        ('DCD', mda.coordinates.DCD.DCDWriter, True, True),
        # ('ENT', mda.coordinates.PDB.PDBWriter, True, False),
        ('GRO', mda.coordinates.GRO.GROWriter, True, False),
        ('LAMMPS', mda.coordinates.LAMMPS.DCDWriter, True, True),
        ('MOL2', mda.coordinates.MOL2.MOL2Writer, True, True),
        ('NCDF', mda.coordinates.TRJ.NCDFWriter, True, True),
        ('NULL', mda.coordinates.null.NullWriter, True, True),
        # ('PDB', mda.coordinates.PDB.PDBWriter, True, True), special case, done separately
        ('PDBQT', mda.coordinates.PDBQT.PDBQTWriter, True, False),
        ('PQR', mda.coordinates.PQR.PQRWriter, True, False),
        ('TRR', mda.coordinates.TRR.TRRWriter, True, True),
        ('XTC', mda.coordinates.XTC.XTCWriter, True, True),
        ('XYZ', mda.coordinates.XYZ.XYZWriter, True, True),
        ('TRZ', mda.coordinates.TRZ.TRZWriter, True, True),
    ]

    @pytest.mark.parametrize('format, writer',
                             [(format_tuple[0], format_tuple[1]) for
                              format_tuple in formats if
                              format_tuple[2] is True])
    def test_singleframe(self, format, writer):
        assert mda.coordinates.core.get_writer_for('this', format=format,
                                                   multiframe=False) == writer

    @pytest.mark.parametrize('format', [(format_tuple[0], format_tuple[1]) for
                                        format_tuple in formats if
                                        format_tuple[2] is False])
    def test_singleframe_fails(self, format):
        with pytest.raises(TypeError):
            mda.coordinates.core.get_writer_for('this', format=format,
                                                multiframe=False)

    @pytest.mark.parametrize('format, writer',
                             [(format_tuple[0], format_tuple[1]) for
                              format_tuple in formats if
                              format_tuple[3] is True])
    def test_multiframe(self, format, writer):
        assert mda.coordinates.core.get_writer_for('this', format=format,
                                                   multiframe=True) == writer

    @pytest.mark.parametrize('format',
                             [format_tuple[0] for format_tuple in formats if
                              format_tuple[3] is False])
    def test_multiframe_fails(self, format):
        with pytest.raises(TypeError):
            mda.coordinates.core.get_writer_for('this', format=format,
                                                multiframe=True)

    def test_get_writer_for_pdb(self):
        assert mda.coordinates.core.get_writer_for('this', format='PDB',
                                                   multiframe=False) == mda.coordinates.PDB.PDBWriter
        assert mda.coordinates.core.get_writer_for('this', format='PDB',
                                                   multiframe=True) == mda.coordinates.PDB.MultiPDBWriter
        assert mda.coordinates.core.get_writer_for('this', format='ENT',
                                                   multiframe=False) == mda.coordinates.PDB.PDBWriter
        assert mda.coordinates.core.get_writer_for('this', format='ENT',
                                                   multiframe=True) == mda.coordinates.PDB.MultiPDBWriter


class TestBlocksOf(object):
    def test_blocks_of_1(self):
        arr = np.arange(16).reshape(4, 4)

        view = util.blocks_of(arr, 1, 1)

        assert view.shape == (4, 1, 1)
        assert_array_almost_equal(view,
                                  np.array([[[0]], [[5]], [[10]], [[15]]]))

        # Change my view, check changes are reflected in arr
        view[:] = 1001

        assert_array_almost_equal(arr,
                                  np.array([[1001, 1, 2, 3],
                                            [4, 1001, 6, 7],
                                            [8, 9, 1001, 11],
                                            [12, 13, 14, 1001]]))

    def test_blocks_of_2(self):
        arr = np.arange(16).reshape(4, 4)

        view = util.blocks_of(arr, 2, 2)

        assert view.shape == (2, 2, 2)
        assert_array_almost_equal(view, np.array([[[0, 1], [4, 5]],
                                                  [[10, 11], [14, 15]]]))

        view[0] = 100
        view[1] = 200

        assert_array_almost_equal(arr,
                                  np.array([[100, 100, 2, 3],
                                            [100, 100, 6, 7],
                                            [8, 9, 200, 200],
                                            [12, 13, 200, 200]]))

    def test_blocks_of_3(self):
        # testing non square array
        arr = np.arange(32).reshape(8, 4)

        view = util.blocks_of(arr, 2, 1)

        assert view.shape == (4, 2, 1)

    def test_blocks_of_4(self):
        # testing block exceeding array size results in empty view
        arr = np.arange(4).reshape(2, 2)
        view = util.blocks_of(arr, 3, 3)
        assert view.shape == (0, 3, 3)
        view[:] = 100
        assert_array_equal(arr, np.arange(4).reshape(2, 2))

    def test_blocks_of_ValueError(self):
        arr = np.arange(16).reshape(4, 4)
        with pytest.raises(ValueError):
            util.blocks_of(arr, 2, 1)  # blocks don't fit
        with pytest.raises(ValueError):
            util.blocks_of(arr[:, ::2], 2, 1)  # non-contiguous input


@pytest.mark.parametrize('arr,answer', [
    ([2, 3, 4, 7, 8, 9, 10, 15, 16], [[2, 3, 4], [7, 8, 9, 10], [15, 16]]),
    ([11, 12, 13, 14, 15, 16], [[11, 12, 13, 14, 15, 16]]),
    ([1, 2, 2, 2, 3, 6], [[1, 2, 2, 2, 3], [6]])
])
def test_group_same_or_consecutive_integers(arr, answer):
    assert_equal(util.group_same_or_consecutive_integers(arr), answer)


class TestNamespace(object):
    @staticmethod
    @pytest.fixture()
    def ns():
        return util.Namespace()

    def test_getitem(self, ns):
        ns.this = 42
        assert ns['this'] == 42

    def test_getitem_KeyError(self, ns):
        with pytest.raises(KeyError):
            dict.__getitem__(ns, 'this')

    def test_setitem(self, ns):
        ns['this'] = 42

        assert ns['this'] == 42

    def test_delitem(self, ns):
        ns['this'] = 42
        assert 'this' in ns
        del ns['this']
        assert 'this' not in ns

    def test_delitem_AttributeError(self, ns):
        with pytest.raises(AttributeError):
            del ns.this

    def test_setattr(self, ns):
        ns.this = 42

        assert ns.this == 42

    def test_getattr(self, ns):
        ns['this'] = 42

        assert ns.this == 42

    def test_getattr_AttributeError(self, ns):
        with pytest.raises(AttributeError):
            getattr(ns, 'this')

    def test_delattr(self, ns):
        ns['this'] = 42

        assert 'this' in ns
        del ns.this
        assert 'this' not in ns

    def test_eq(self, ns):
        ns['this'] = 42

        ns2 = util.Namespace()
        ns2['this'] = 42

        assert ns == ns2

    def test_len(self, ns):
        assert len(ns) == 0
        ns['this'] = 1
        ns['that'] = 2
        assert len(ns) == 2

    def test_iter(self, ns):
        ns['this'] = 12
        ns['that'] = 24
        ns['other'] = 48

        seen = []
        for val in ns:
            seen.append(val)
        for val in ['this', 'that', 'other']:
            assert val in seen


class TestTruncateInteger(object):
    @pytest.mark.parametrize('a, b', [
        ((1234, 1), 4),
        ((1234, 2), 34),
        ((1234, 3), 234),
        ((1234, 4), 1234),
        ((1234, 5), 1234),
    ])
    def test_ltruncate_int(self, a, b):
        assert util.ltruncate_int(*a) == b


class TestFlattenDict(object):
    def test_flatten_dict(self):
        d = {
            'A': {1: ('a', 'b', 'c')},
            'B': {2: ('c', 'd', 'e')},
            'C': {3: ('f', 'g', 'h')}
        }
        result = util.flatten_dict(d)

        for k in result:
            assert type(k) == tuple
            assert len(k) == 2
            assert k[0] in d
            assert k[1] in d[k[0]]
            assert result[k] in d[k[0]].values()


class TestStaticVariables(object):
    """Tests concerning the decorator @static_variables
    """

    def test_static_variables(self):
        x = [0]

        @static_variables(foo=0, bar={'test': x})
        def myfunc():
            assert myfunc.foo == 0
            assert type(myfunc.bar) is type(dict())
            if 'test2' not in myfunc.bar:
                myfunc.bar['test2'] = "a"
            else:
                myfunc.bar['test2'] += "a"
            myfunc.bar['test'][0] += 1
            return myfunc.bar['test']

        assert hasattr(myfunc, 'foo')
        assert hasattr(myfunc, 'bar')

        y = myfunc()
        assert y is x
        assert x[0] == 1
        assert myfunc.bar['test'][0] == 1
        assert myfunc.bar['test2'] == "a"

        x = [0]
        y = myfunc()
        assert y is not x
        assert myfunc.bar['test'][0] == 2
        assert myfunc.bar['test2'] == "aa"


class TestWarnIfNotUnique(object):
    """Tests concerning the decorator @warn_if_not_unique
    """

    def warn_msg(self, func, group, group_name):
        msg = ("{}.{}(): {} {} contains duplicates. Results might be "
               "biased!".format(group.__class__.__name__, func.__name__,
                                group_name, group.__repr__()))
        return msg

    def test_warn_if_not_unique(self, atoms):
        # Check that the warn_if_not_unique decorator has a "static variable"
        # warn_if_not_unique.warned:
        assert hasattr(warn_if_not_unique, 'warned')
        assert warn_if_not_unique.warned is False

    def test_warn_if_not_unique_once_outer(self, atoms):

        # Construct a scenario with two nested functions, each one decorated
        # with @warn_if_not_unique:

        @warn_if_not_unique
        def inner(group):
            if not group.isunique:
                # The inner function should not trigger a warning, and the state
                # of warn_if_not_unique.warned should reflect that:
                assert warn_if_not_unique.warned is True
            return 0

        @warn_if_not_unique
        def outer(group):
            return inner(group)

        # Check that no warning is raised for a unique group:
        assert atoms.isunique

        with warnings.catch_warnings():
            warnings.simplefilter("error")
            x = outer(atoms)
            assert x == 0

        # Check that a warning is raised for a group with duplicates:
        ag = atoms + atoms[0]
        msg = self.warn_msg(outer, ag, "'ag'")
        with pytest.warns(DuplicateWarning) as w:
            assert warn_if_not_unique.warned is False
            x = outer(ag)
            # Assert that the "warned" state is restored:
            assert warn_if_not_unique.warned is False
            # Check correct function execution:
            assert x == 0
            # Only one warning must have been raised:
            assert len(w) == 1
            # For whatever reason pytest.warns(DuplicateWarning, match=msg)
            # doesn't work, so we compare the recorded warning message instead:
            assert w[0].message.args[0] == msg
            # Make sure the warning uses the correct stacklevel and references
            # this file instead of MDAnalysis/lib/util.py:
            assert w[0].filename == __file__

    def test_warned_state_restored_on_failure(self, atoms):

        # A decorated function raising an exception:
        @warn_if_not_unique
        def thisfails(group):
            raise ValueError()

        ag = atoms + atoms[0]
        msg = self.warn_msg(thisfails, ag, "'ag'")
        with pytest.warns(DuplicateWarning) as w:
            assert warn_if_not_unique.warned is False
            with pytest.raises(ValueError):
                thisfails(ag)
            # Assert that the "warned" state is restored despite `thisfails`
            # raising an exception:
            assert warn_if_not_unique.warned is False
            assert len(w) == 1
            assert w[0].message.args[0] == msg
            assert w[0].filename == __file__

    def test_warn_if_not_unique_once_inner(self, atoms):

        # Construct a scenario with two nested functions, each one decorated
        # with @warn_if_not_unique, but the outer function adds a duplicate
        # to the group:

        @warn_if_not_unique
        def inner(group):
            return 0

        @warn_if_not_unique
        def outer(group):
            dupgroup = group + group[0]
            return inner(dupgroup)

        # Check that even though outer() is called the warning is raised for
        # inner():
        msg = self.warn_msg(inner, atoms + atoms[0], "'dupgroup'")
        with pytest.warns(DuplicateWarning) as w:
            assert warn_if_not_unique.warned is False
            x = outer(atoms)
            # Assert that the "warned" state is restored:
            assert warn_if_not_unique.warned is False
            # Check correct function execution:
            assert x == 0
            # Only one warning must have been raised:
            assert len(w) == 1
            assert w[0].message.args[0] == msg
            assert w[0].filename == __file__

    def test_warn_if_not_unique_multiple_references(self, atoms):
        ag = atoms + atoms[0]
        aag = ag
        aaag = aag

        @warn_if_not_unique
        def func(group):
            return group.isunique

        # Check that the warning message contains the names of all references to
        # the group in alphabetic order:
        msg = self.warn_msg(func, ag, "'aaag' a.k.a. 'aag' a.k.a. 'ag'")
        with pytest.warns(DuplicateWarning) as w:
            x = func(ag)
            # Assert that the "warned" state is restored:
            assert warn_if_not_unique.warned is False
            # Check correct function execution:
            assert x is False
            # Check warning message:
            assert w[0].message.args[0] == msg
            # Check correct file referenced:
            assert w[0].filename == __file__

    def test_warn_if_not_unique_unnamed(self, atoms):

        @warn_if_not_unique
        def func(group):
            pass

        msg = self.warn_msg(func, atoms + atoms[0],
                            "'unnamed {}'".format(atoms.__class__.__name__))
        with pytest.warns(DuplicateWarning) as w:
            func(atoms + atoms[0])
            # Check warning message:
            assert w[0].message.args[0] == msg

    def test_warn_if_not_unique_fails_for_non_groupmethods(self):

        @warn_if_not_unique
        def func(group):
            pass

        class dummy(object):
            pass

        with pytest.raises(AttributeError):
            func(dummy())

    def test_filter_duplicate_with_userwarning(self, atoms):

        @warn_if_not_unique
        def func(group):
            pass

        with warnings.catch_warnings(record=True) as record:
            warnings.resetwarnings()
            warnings.filterwarnings("ignore", category=UserWarning)
            with warnings.catch_warnings():
                warnings.simplefilter("error")
                func(atoms)
            assert len(record) == 0


class TestCheckCoords(object):
    """Tests concerning the decorator @check_coords
    """

    prec = 6

    def test_default_options(self):
        a_in = np.zeros(3, dtype=np.float32)
        b_in = np.ones(3, dtype=np.float32)
        b_in2 = np.ones((2, 3), dtype=np.float32)

        @check_coords('a', 'b')
        def func(a, b):
            # check that enforce_copy is True by default:
            assert a is not a_in
            assert b is not b_in
            # check that convert_single is True by default:
            assert a.shape == (1, 3)
            assert b.shape == (1, 3)
            return a + b

        # check that allow_single is True by default:
        res = func(a_in, b_in)
        # check that reduce_result_if_single is True by default:
        assert res.shape == (3,)
        # check correct function execution:
        assert_array_equal(res, b_in)

        # check that check_lenghts_match is True by default:
        with pytest.raises(ValueError):
            res = func(a_in, b_in2)

    @pytest.fixture()
    def atomgroup(self):
        u = mda.Universe(PSF, DCD)
        return u.atoms

    # check atomgroup handling with every option except allow_atomgroup
    @pytest.mark.parametrize('enforce_copy', [True, False])
    @pytest.mark.parametrize('enforce_dtype', [True, False])
    @pytest.mark.parametrize('allow_single', [True, False])
    @pytest.mark.parametrize('convert_single', [True, False])
    @pytest.mark.parametrize('reduce_result_if_single', [True, False])
    @pytest.mark.parametrize('check_lengths_match', [True, False])
    def test_atomgroup(self, atomgroup, enforce_copy, enforce_dtype,
                       allow_single, convert_single, reduce_result_if_single,
                       check_lengths_match):
        ag1 = atomgroup
        ag2 = atomgroup

        @check_coords('ag1', 'ag2', enforce_copy=enforce_copy,
                      enforce_dtype=enforce_dtype, allow_single=allow_single,
                      convert_single=convert_single,
                      reduce_result_if_single=reduce_result_if_single,
                      check_lengths_match=check_lengths_match,
                      allow_atomgroup=True)
        def func(ag1, ag2):
            assert_allclose(ag1, ag2)
            assert isinstance(ag1, np.ndarray)
            assert isinstance(ag2, np.ndarray)
            assert ag1.dtype == ag2.dtype == np.float32
            return ag1 + ag2

        res = func(ag1, ag2)

        assert_allclose(res, atomgroup.positions*2)

    def test_atomgroup_not_allowed(self, atomgroup):

        @check_coords('ag1', allow_atomgroup=False)
        def func(ag1):
            return ag1

        with pytest.raises(TypeError, match="allow_atomgroup is False"):
            _ = func(atomgroup)

    def test_atomgroup_not_allowed_default(self, atomgroup):

        @check_coords('ag1')
        def func_default(ag1):
            return ag1

        with pytest.raises(TypeError, match="allow_atomgroup is False"):
            _ = func_default(atomgroup)

    def test_enforce_copy(self):

        a_2d = np.ones((1, 3), dtype=np.float32)
        b_1d = np.zeros(3, dtype=np.float32)
        c_2d = np.zeros((1, 6), dtype=np.float32)[:, ::2]
        d_2d = np.zeros((1, 3), dtype=np.int64)

        @check_coords('a', 'b', 'c', 'd', enforce_copy=False)
        def func(a, b, c, d):
            # Assert that if enforce_copy is False:
            # no copy is made if input shape, order, and dtype are correct:
            assert a is a_2d
            # a copy is made if input shape has to be changed:
            assert b is not b_1d
            # a copy is made if input order has to be changed:
            assert c is not c_2d
            # a copy is made if input dtype has to be changed:
            assert d is not d_2d
            # Assert correct dtype conversion:
            assert d.dtype == np.float32
            assert_almost_equal(d, d_2d, self.prec)
            # Assert all shapes are converted to (1, 3):
            assert a.shape == b.shape == c.shape == d.shape == (1, 3)
            return a + b + c + d

        # Call func() to:
        # - test the above assertions
        # - ensure that input of single coordinates is simultaneously possible
        #   with different shapes (3,) and (1, 3)
        res = func(a_2d, b_1d, c_2d, d_2d)
        # Since some inputs are not 1d, even though reduce_result_if_single is
        # True, the result must have shape (1, 3):
        assert res.shape == (1, 3)
        # check correct function execution:
        assert_array_equal(res, a_2d)

    def test_no_allow_single(self):

        @check_coords('a', allow_single=False)
        def func(a):
            pass

        with pytest.raises(ValueError) as err:
            func(np.zeros(3, dtype=np.float32))
            assert err.msg == ("func(): a.shape must be (n, 3), got (3,).")

    def test_no_convert_single(self):

        a_1d = np.arange(-3, 0, dtype=np.float32)

        @check_coords('a', enforce_copy=False, convert_single=False)
        def func(a):
            # assert no conversion and no copy were performed:
            assert a is a_1d
            return a

        res = func(a_1d)
        # Assert result has been reduced:
        assert res == a_1d[0]
        assert type(res) is np.float32

    def test_no_reduce_result_if_single(self):

        a_1d = np.zeros(3, dtype=np.float32)

        # Test without shape conversion:
        @check_coords('a', enforce_copy=False, convert_single=False,
                      reduce_result_if_single=False)
        def func(a):
            return a

        res = func(a_1d)
        # make sure the input array is just passed through:
        assert res is a_1d

        # Test with shape conversion:
        @check_coords('a', enforce_copy=False, reduce_result_if_single=False)
        def func(a):
            return a

        res = func(a_1d)
        assert res.shape == (1, 3)
        assert_array_equal(res[0], a_1d)

    def test_no_check_lengths_match(self):

        a_2d = np.zeros((1, 3), dtype=np.float32)
        b_2d = np.zeros((3, 3), dtype=np.float32)

        @check_coords('a', 'b', enforce_copy=False, check_lengths_match=False)
        def func(a, b):
            return a, b

        res_a, res_b = func(a_2d, b_2d)
        # Assert arrays are just passed through:
        assert res_a is a_2d
        assert res_b is b_2d

    def test_atomgroup_mismatched_lengths(self):
        u = mda.Universe(PSF, DCD)
        ag1 = u.select_atoms("index 0 to 10")
        ag2 = u.atoms

        @check_coords('ag1', 'ag2', check_lengths_match=True,
                      allow_atomgroup=True)
        def func(ag1, ag2):

            return ag1, ag2

        with pytest.raises(ValueError, match="must contain the same number of "
                           "coordinates"):
            _, _ = func(ag1, ag2)

    def test_invalid_input(self):

        a_inv_dtype = np.array([['hello', 'world', '!']])
        a_inv_type = [[0., 0., 0.]]
        a_inv_shape_1d = np.zeros(6, dtype=np.float32)
        a_inv_shape_2d = np.zeros((3, 2), dtype=np.float32)

        @check_coords('a')
        def func(a):
            pass

        with pytest.raises(TypeError) as err:
            func(a_inv_dtype)
            assert err.msg.startswith("func(): a.dtype must be convertible to "
                                      "float32, got ")

        with pytest.raises(TypeError) as err:
            func(a_inv_type)
            assert err.msg == ("func(): Parameter 'a' must be a numpy.ndarray, "
                               "got <class 'list'>.")

        with pytest.raises(ValueError) as err:
            func(a_inv_shape_1d)
            assert err.msg == ("func(): a.shape must be (3,) or (n, 3), got "
                               "(6,).")

        with pytest.raises(ValueError) as err:
            func(a_inv_shape_2d)
            assert err.msg == ("func(): a.shape must be (3,) or (n, 3), got "
                               "(3, 2).")

    def test_usage_with_kwargs(self):

        a_2d = np.zeros((1, 3), dtype=np.float32)

        @check_coords('a', enforce_copy=False)
        def func(a, b, c=0):
            return a, b, c

        # check correct functionality if passed as keyword argument:
        a, b, c = func(a=a_2d, b=0, c=1)
        assert a is a_2d
        assert b == 0
        assert c == 1

    def test_wrong_func_call(self):

        @check_coords('a', enforce_copy=False)
        def func(a, b, c=0):
            pass

        # Make sure invalid call marker is present:
        func._invalid_call = False

        # usage with posarg doubly defined:
        assert not func._invalid_call
        with pytest.raises(TypeError):
            func(0, a=0)  # pylint: disable=redundant-keyword-arg
        assert func._invalid_call
        func._invalid_call = False

        # usage with missing posargs:
        assert not func._invalid_call
        with pytest.raises(TypeError):
            func(0)
        assert func._invalid_call
        func._invalid_call = False

        # usage with missing posargs (supplied as kwargs):
        assert not func._invalid_call
        with pytest.raises(TypeError):
            func(a=0, c=1)
        assert func._invalid_call
        func._invalid_call = False

        # usage with too many posargs:
        assert not func._invalid_call
        with pytest.raises(TypeError):
            func(0, 0, 0, 0)
        assert func._invalid_call
        func._invalid_call = False

        # usage with unexpected kwarg:
        assert not func._invalid_call
        with pytest.raises(TypeError):
            func(a=0, b=0, c=1, d=1)  # pylint: disable=unexpected-keyword-arg
        assert func._invalid_call
        func._invalid_call = False

    def test_wrong_decorator_usage(self):

        # usage without parantheses:
        @check_coords
        def func():
            pass

        with pytest.raises(TypeError):
            func()

        # usage without arguments:
        with pytest.raises(ValueError) as err:
            @check_coords()
            def func():
                pass

            assert err.msg == ("Decorator check_coords() cannot be used "
                               "without positional arguments.")

        # usage with defaultarg:
        with pytest.raises(ValueError) as err:
            @check_coords('a')
            def func(a=1):
                pass

            assert err.msg == ("In decorator check_coords(): Name 'a' doesn't "
                               "correspond to any positional argument of the "
                               "decorated function func().")

        # usage with invalid parameter name:
        with pytest.raises(ValueError) as err:
            @check_coords('b')
            def func(a):
                pass

            assert err.msg == ("In decorator check_coords(): Name 'b' doesn't "
                               "correspond to any positional argument of the "
                               "decorated function func().")


@pytest.mark.parametrize("old_name", (None, "MDAnalysis.Universe"))
@pytest.mark.parametrize("new_name", (None, "Multiverse"))
@pytest.mark.parametrize("remove", (None, "99.0.0", 2099))
@pytest.mark.parametrize("message", (None, "use the new stuff"))
def test_deprecate(old_name, new_name, remove, message, release="2.7.1"):
    def AlternateUniverse(anything):
        # important: first line needs to be """\ so that textwrap.dedent()
        # works
        """\
        AlternateUniverse provides a true view of the Universe.

        Parameters
        ----------
        anything : object

        Returns
        -------
        truth

        """
        return True

    oldfunc = util.deprecate(AlternateUniverse, old_name=old_name,
                             new_name=new_name,
                             release=release, remove=remove,
                             message=message)
    # match_expr changed to match (Issue 2329)
    with pytest.warns(DeprecationWarning, match="`.+` is deprecated"):
        oldfunc(42)

    doc = oldfunc.__doc__
    name = old_name if old_name else AlternateUniverse.__name__

    deprecation_line_1 = ".. deprecated:: {0}".format(release)
    assert re.search(deprecation_line_1, doc)

    if message:
        deprecation_line_2 = message
    else:
        if new_name is None:
            default_message = "`{0}` is deprecated!".format(name)
        else:
            default_message = "`{0}` is deprecated, use `{1}` instead!".format(
                name, new_name)
        deprecation_line_2 = default_message
    assert re.search(deprecation_line_2, doc)

    if remove:
        deprecation_line_3 = "`{0}` will be removed in release {1}".format(
            name,  remove)
        assert re.search(deprecation_line_3, doc)

    # check that the old docs are still present
    assert re.search(textwrap.dedent(AlternateUniverse.__doc__), doc)


def test_deprecate_missing_release_ValueError():
    with pytest.raises(ValueError):
        util.deprecate(mda.Universe)


def test_set_function_name(name="bar"):
    def foo():
        pass
    util._set_function_name(foo, name)
    assert foo.__name__ == name


@pytest.mark.parametrize("text",
                         ("",
                          "one line text",
                          "  one line with leading space",
                          "multiline\n\n   with some\n   leading space",
                          "   multiline\n\n   with all\n   leading space"))
def test_dedent_docstring(text):
    doc = util.dedent_docstring(text)
    for line in doc.splitlines():
        assert line == line.lstrip()


class TestCheckBox(object):

    prec = 6
    ref_ortho = np.ones(3, dtype=np.float32)
    ref_tri_vecs = np.array([[1, 0, 0], [0, 1, 0], [0, 2 ** 0.5, 2 ** 0.5]],
                            dtype=np.float32)

    @pytest.mark.parametrize('box',
                             ([1, 1, 1, 90, 90, 90],
                              (1, 1, 1, 90, 90, 90),
                                 ['1', '1', 1, 90, '90', '90'],
                                 ('1', '1', 1, 90, '90', '90'),
                                 np.array(['1', '1', 1, 90, '90', '90']),
                                 np.array([1, 1, 1, 90, 90, 90],
                                          dtype=np.float32),
                                 np.array([1, 1, 1, 90, 90, 90],
                                          dtype=np.float64),
                                 np.array([1, 1, 1, 1, 1, 1,
                                           90, 90, 90, 90, 90, 90],
                                          dtype=np.float32)[::2]))
    def test_check_box_ortho(self, box):
        boxtype, checked_box = util.check_box(box)
        assert boxtype == 'ortho'
        assert_allclose(checked_box, self.ref_ortho)
        assert checked_box.dtype == np.float32
        assert checked_box.flags['C_CONTIGUOUS']

    def test_check_box_None(self):
        with pytest.raises(ValueError, match="Box is None"):
            util.check_box(None)

    @pytest.mark.parametrize('box',
                             ([1, 1, 2, 45, 90, 90],
                              (1, 1, 2, 45, 90, 90),
                                 ['1', '1', 2, 45, '90', '90'],
                                 ('1', '1', 2, 45, '90', '90'),
                                 np.array(['1', '1', 2, 45, '90', '90']),
                                 np.array([1, 1, 2, 45, 90, 90],
                                          dtype=np.float32),
                                 np.array([1, 1, 2, 45, 90, 90],
                                          dtype=np.float64),
                                 np.array([1, 1, 1, 1, 2, 2,
                                           45, 45, 90, 90, 90, 90],
                                          dtype=np.float32)[::2]))
    def test_check_box_tri_vecs(self, box):
        boxtype, checked_box = util.check_box(box)
        assert boxtype == 'tri_vecs'
        assert_almost_equal(checked_box, self.ref_tri_vecs, self.prec)
        assert checked_box.dtype == np.float32
        assert checked_box.flags['C_CONTIGUOUS']

    def test_check_box_wrong_data(self):
        with pytest.raises(ValueError):
            wrongbox = ['invalid', 1, 1, 90, 90, 90]
            boxtype, checked_box = util.check_box(wrongbox)

    def test_check_box_wrong_shape(self):
        with pytest.raises(ValueError):
            wrongbox = np.ones((3, 3), dtype=np.float32)
            boxtype, checked_box = util.check_box(wrongbox)


class StoredClass:
    """
    A simple class that takes positional and keyword arguments of various types
    """
    @store_init_arguments
    def __init__(self, a, b, /, *args, c="foo", d="bar", e="foobar", **kwargs):
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.e = e
        self.args = args
        self.kwargs = kwargs

    def copy(self):
        kwargs = copy.deepcopy(self._kwargs)
        args = kwargs.pop('args', tuple())
        new = self.__class__(kwargs.pop('a'), kwargs.pop('b'),
                             *args, **kwargs)
        return new


class TestStoreInitArguments:
    def test_store_arguments_default(self):
        store = StoredClass('parsnips', ['roast'])
        assert store.a == store._kwargs['a'] == 'parsnips'
        assert store.b is store._kwargs['b'] == ['roast']
        assert store._kwargs['c'] == 'foo'
        assert store._kwargs['d'] == 'bar'
        assert store._kwargs['e'] == 'foobar'
        assert 'args' not in store._kwargs.keys()
        assert 'kwargs' not in store._kwargs.keys()
        assert store.args is ()

        store2 = store.copy()
        assert store2.__dict__ == store.__dict__
        assert store2.__dict__["b"] is not store.__dict__["b"]

    def test_store_arguments_withkwargs(self):
        store = StoredClass('parsnips', 'roast', 'honey', 'glaze', c='richard',
                            d='has', e='a', f='recipe', g='allegedly')
        assert store.a == store._kwargs['a'] == "parsnips"
        assert store.b == store._kwargs['b'] == "roast"
        assert store.c == store._kwargs['c'] == "richard"
        assert store.d == store._kwargs['d'] == "has"
        assert store.e == store._kwargs['e'] == "a"
        assert store.kwargs['f'] == store._kwargs['f'] == "recipe"
        assert store.kwargs['g'] == store._kwargs['g'] == "allegedly"
        assert store.args[0] == store._kwargs['args'][0] == "honey"
        assert store.args[1] == store._kwargs['args'][1] == "glaze"

        store2 = store.copy()
        assert store2.__dict__ == store.__dict__


@pytest.mark.xfail(os.name == 'nt',
                   reason="util.which does not get right binary on Windows")
def test_which():
    wmsg = "This method is deprecated"

    with pytest.warns(DeprecationWarning, match=wmsg):
        assert util.which('python') == shutil.which('python')
