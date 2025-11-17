import numpy as np
import pytest

import MDAnalysis as mda
from MDAnalysis.core.topologyattrs import Bonds
from MDAnalysis.exceptions import NoDataError
from MDAnalysisTests.datafiles import (
    Make_Whole,
    PSF,
    DCD,
    TPR,
    GRO,
    fullerene,
    two_water_gro,
)

from MDAnalysis.lib import mdamath
from numpy.testing import (
    assert_array_equal,
    assert_almost_equal,
    assert_allclose,
    assert_array_almost_equal,
)
from itertools import combinations_with_replacement as comb_wr


class TestGeometryFunctions:
    e1, e2, e3 = np.eye(3)
    a = np.array([np.cos(np.pi / 3), np.sin(np.pi / 3), 0])
    null = np.zeros(3)

    @pytest.mark.parametrize(
        "x_axis, y_axis, value",
        [
            (e1, e2, np.pi / 2),
            (e1, a, np.pi / 3),
            (2 * e1, e2, np.pi / 2),
            (-2 * e1, e2, np.pi - np.pi / 2),
            (23.3 * e1, a, np.pi / 3),
            (e1, null, np.nan),
            (a, a, 0.0),
        ],
    )
    def test_angle(self, x_axis, y_axis, value):
        result = mdamath.angle(x_axis, y_axis)
        if np.isnan(value):
            assert np.isnan(result)
        else:
            np.testing.assert_allclose(result, value)

    @pytest.mark.parametrize(
        "x_axis, y_axis, value",
        [
            (-2.3456e7 * e1, 3.4567e-6 * e1, np.pi),
            (2.3456e7 * e1, 3.4567e-6 * e1, 0.0),
        ],
    )
    def test_angle_collinear(self, x_axis, y_axis, value):
        result = mdamath.angle(x_axis, y_axis)
        np.testing.assert_allclose(result, value)

    @pytest.mark.parametrize("x", np.linspace(0, np.pi, 20))
    def test_angle_symmetry(self, x):
        v1 = np.array([np.cos(x), np.sin(x), 0])
        v2 = np.array([1, 0, 0])
        assert np.allclose(mdamath.angle(v1, v2), x)

    @pytest.mark.parametrize(
        "vector, value", [(e3, 1), (a, np.linalg.norm(a)), (null, 0.0)]
    )
    def test_norm(self, vector, value):
        assert np.allclose(mdamath.norm(vector), value)

    @pytest.mark.parametrize("x", np.linspace(0, np.pi, 20))
    def test_normal_unit(self, x):
        v1 = np.array([np.cos(x), np.sin(x), 0])
        v2 = np.array([1, 0, 0])
        n = mdamath.normal(v1, v2)
        if np.allclose(v1, v2):
            assert np.allclose(n, np.zeros(3))
        else:
            assert np.allclose(np.linalg.norm(n), 1)

    @pytest.mark.parametrize(
        "vec1, vec2, value", [(e1, e2, e3), (e1, null, 0.0)]
    )
    def test_normal(self, vec1, vec2, value):
        n = mdamath.normal(vec1, vec2)
        if isinstance(value, float):
            assert np.allclose(np.linalg.norm(n), value)
        else:
            assert np.allclose(n, value)

    def test_angle_lower_clip(self):
        # Test for values slightly less than -1.0
        a = np.array([1.0, 0.0, 0.0])
        b = np.array([-1.0, 0.0, 0.0])
        result = mdamath.angle(a, b)
        assert np.allclose(result, np.pi)

    def test_stp(self):
        v1 = np.array([1, 0, 0])
        v2 = np.array([0, 1, 0])
        v3 = np.array([0, 0, 1])
        assert np.allclose(mdamath.stp(v1, v2, v3), 1.0)

    def test_dihedral(self):
        ab = np.array([1, 0, 0])
        bc = np.array([0, 1, 0])
        cd = np.array([0, 0, 1])
        result = mdamath.dihedral(ab, bc, cd)
        assert np.allclose(result, -np.pi / 2)

    def test_pdot(self):
        a = np.array([[1, 2, 3], [4, 5, 6]])
        b = np.array([[1, 0, 0], [0, 1, 0]])
        result = mdamath.pdot(a, b)
        np.testing.assert_array_equal(result, [1, 5])

    def test_pnorm(self):
        a = np.array([[3, 4, 0], [0, 0, 5]])
        result = mdamath.pnorm(a)
        np.testing.assert_allclose(result, [5, 5])


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

    @pytest.mark.parametrize("lengths", comb_wr([-1, 0, 1, 2], 3))
    @pytest.mark.parametrize(
        "angles", comb_wr([-10, 0, 20, 70, 90, 120, 180], 3)
    )
    def test_triclinic_vectors(self, lengths, angles):
        box = lengths + angles
        ref = self.ref_trivecs(box)
        res = mdamath.triclinic_vectors(box)
        assert_array_equal(res, ref)
        # check for default dtype:
        assert res.dtype == np.float32
        # belts and braces, make sure upper triangle is always zero:
        assert not (res[0, 1] or res[0, 2] or res[1, 2])

    @pytest.mark.parametrize("alpha", (60, 90))
    @pytest.mark.parametrize("beta", (60, 90))
    @pytest.mark.parametrize("gamma", (60, 90))
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

    @pytest.mark.parametrize("dtype", (int, float, np.float32, np.float64))
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
                        *mdamath.triclinic_vectors(ref)
                    )
                    if not np.all(res == 0.0):
                        assert_almost_equal(res, ref, 5)

    @pytest.mark.parametrize(
        "angles",
        (
            [70, 70, 70],
            [70, 70, 90],
            [70, 90, 70],
            [90, 70, 70],
            [70, 90, 90],
            [90, 70, 90],
            [90, 90, 70],
        ),
    )
    def test_triclinic_vectors_box_cycle_exact(self, angles):
        # These cycles were inexact prior to PR #2201
        ref = np.array([10.1, 10.1, 10.1] + angles, dtype=np.float32)
        res = mdamath.triclinic_box(*mdamath.triclinic_vectors(ref))
        assert_allclose(res, ref)

    @pytest.mark.parametrize("lengths", comb_wr([-1, 0, 1, 2], 3))
    @pytest.mark.parametrize(
        "angles", comb_wr([-10, 0, 20, 70, 90, 120, 180], 3)
    )
    def test_triclinic_box(self, lengths, angles):
        tri_vecs = self.ref_trivecs_unsafe(lengths + angles)
        ref = self.ref_tribox(tri_vecs)
        res = mdamath.triclinic_box(*tri_vecs)
        assert_array_equal(res, ref)
        assert res.dtype == ref.dtype

    @pytest.mark.parametrize("lengths", comb_wr([-1, 0, 1, 2], 3))
    @pytest.mark.parametrize(
        "angles", comb_wr([-10, 0, 20, 70, 90, 120, 180], 3)
    )
    def test_box_volume(self, lengths, angles):
        box = np.array(lengths + angles, dtype=np.float32)
        assert_almost_equal(
            mdamath.box_volume(box),
            np.linalg.det(self.ref_trivecs(box)),
            decimal=5,
        )

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

    @pytest.mark.parametrize("shape", ((0,), (3, 2), (2, 3), (1, 1, 3, 1)))
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
            forces=False,
        )
        ts = u.trajectory.ts
        ts.frame = 0
        ts.dimensions = [10, 10, 10, 90, 90, 90]
        # assert ts.dimensions.dtype == np.float64
        # not applicable since #2213
        ts.positions = np.array(
            [
                [1, 1, 1],
                [9, 9, 9],
            ],
            dtype=np.float32,
        )
        u.add_TopologyAttr(Bonds([(0, 1)]))
        mdamath.make_whole(u.atoms)
        assert_array_almost_equal(
            u.atoms.positions,
            np.array(
                [
                    [1, 1, 1],
                    [-1, -1, -1],
                ],
                dtype=np.float32,
            ),
        )

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
        universe.dimensions = [0.0, 0.0, 0.0, 90.0, 90.0, 90.0]
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
        assert_array_almost_equal(
            universe.atoms[4].position,
            np.array([110.0, 50.0, 0.0]),
            decimal=self.prec,
        )
        assert_array_almost_equal(
            universe.atoms[5].position,
            np.array([110.0, 60.0, 0.0]),
            decimal=self.prec,
        )
        assert_array_almost_equal(
            universe.atoms[6].position,
            np.array([110.0, 40.0, 0.0]),
            decimal=self.prec,
        )
        assert_array_almost_equal(
            universe.atoms[7].position,
            np.array([120.0, 50.0, 0.0]),
            decimal=self.prec,
        )

    def test_solve_2(self, universe, ag):
        # use but specify the center atom

        refpos = universe.atoms[4:8].positions.copy()

        mdamath.make_whole(ag, reference_atom=universe.residues[0].atoms[4])

        assert_array_almost_equal(universe.atoms[4:8].positions, refpos)
        assert_array_almost_equal(
            universe.atoms[0].position,
            np.array([-20.0, 50.0, 0.0]),
            decimal=self.prec,
        )
        assert_array_almost_equal(
            universe.atoms[1].position,
            np.array([-10.0, 50.0, 0.0]),
            decimal=self.prec,
        )
        assert_array_almost_equal(
            universe.atoms[2].position,
            np.array([-10.0, 60.0, 0.0]),
            decimal=self.prec,
        )
        assert_array_almost_equal(
            universe.atoms[3].position,
            np.array([-10.0, 40.0, 0.0]),
            decimal=self.prec,
        )

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
        assert_array_almost_equal(
            universe.atoms[4].position, np.array([110.0, 50.0, 0.0])
        )
        assert_array_almost_equal(
            universe.atoms[5].position, np.array([110.0, 60.0, 0.0])
        )
        assert_array_almost_equal(
            universe.atoms[6].position, np.array([110.0, 40.0, 0.0])
        )

    def test_double_frag_short_bonds(self, universe, ag):
        # previous bug where if two fragments are given
        # but all
        # were short, the algorithm didn't
        # complain
        mdamath.make_whole(ag)
        with pytest.raises(ValueError):
            mdamath.make_whole(universe.atoms)

    def test_make_whole_triclinic(self):
        u = mda.Universe(TPR, GRO)
        thing = u.select_atoms("not resname SOL NA+")
        mdamath.make_whole(thing)

        blengths = thing.bonds.values()

        assert blengths.max() < 2.0

    def test_make_whole_fullerene(self):
        # lots of circular bonds as a nice pathological case
        u = mda.Universe(fullerene)

        bbox = u.atoms.bbox()
        u.dimensions = np.r_[bbox[1] - bbox[0], [90] * 3]

        blengths = u.atoms.bonds.values()
        # kaboom
        u.atoms[::2].translate([u.dimensions[0], -2 * u.dimensions[1], 0.0])
        u.atoms[1::2].translate(
            [0.0, 7 * u.dimensions[1], -5 * u.dimensions[2]]
        )

        mdamath.make_whole(u.atoms)

        assert_array_almost_equal(
            u.atoms.bonds.values(), blengths, decimal=self.prec
        )

    def test_make_whole_multiple_molecules(self):
        u = mda.Universe(two_water_gro, guess_bonds=True)

        for f in u.atoms.fragments:
            mdamath.make_whole(f)

        assert u.atoms.bonds.values().max() < 2.0
