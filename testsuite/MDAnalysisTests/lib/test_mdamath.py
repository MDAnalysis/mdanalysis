import numpy as np
import pytest

from MDAnalysis.lib import mdamath
from numpy.testing import assert_array_equal, assert_almost_equal, assert_allclose
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
