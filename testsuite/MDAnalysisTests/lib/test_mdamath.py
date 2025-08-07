import numpy as np
import pytest

from MDAnalysis.lib import mdamath


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


class TestMatrixOperations:
    def test_triclinic_vectors_box_cycle(self):
        box = np.array([10.0, 20.0, 30.0, 90.0, 90.0, 90.0])
        tri_vecs = mdamath.triclinic_vectors(box)
        box2 = mdamath.triclinic_box(*tri_vecs)
        np.testing.assert_allclose(box, box2, rtol=1e-5)

    def test_triclinic_vectors_invalid(self):
        box = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        tri_vecs = mdamath.triclinic_vectors(box)
        assert np.allclose(tri_vecs, 0)

    def test_box_volume(self):
        box = np.array([10.0, 20.0, 30.0, 90.0, 90.0, 90.0])
        vol = mdamath.box_volume(box)
        assert np.allclose(vol, 6000.0)

    def test_box_volume_invalid(self):
        box = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        vol = mdamath.box_volume(box)
        assert np.allclose(vol, 0.0)

    def test_sarrus_det(self):
        m = np.eye(3)
        assert np.allclose(mdamath.sarrus_det(m), 1.0)
        m = np.zeros((3, 3))
        assert np.allclose(mdamath.sarrus_det(m), 0.0)
        m = np.array(
            [
                [[1, 2, 3], [4, 5, 6], [7, 8, 9]],
                [[2, 0, 1], [3, 0, 0], [5, 1, 1]],
            ]
        )
        dets = mdamath.sarrus_det(m)
        assert dets.shape == (2,)
