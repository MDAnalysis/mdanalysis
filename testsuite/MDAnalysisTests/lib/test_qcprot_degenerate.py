"""Regression checks for collinear QCP inputs (issues #945 and #1747)."""
import numpy as np
import pytest
from numpy.testing import assert_allclose

from MDAnalysis.lib import qcprot
from MDAnalysis.analysis import align, rms


def _issue_coordinates(case):
    if case == "945":
        x = np.ones((100, 3)) + np.arange(100)[:, None]
        y = np.zeros((100, 3)) + np.arange(100)[:, None]
    elif case == "1747":
        x = np.array([[.13662023, -.39526145, .75835566],
                      [-.13662023, .39526145, -.75835566]])
        y = np.array([[-.30452963, .78899221, -.18642158],
                      [.30452963, -.78899221, .18642158]])
    elif case == "perpendicular":
        x = np.array([[-1., 0., 0.], [1., 0., 0.]])
        y = np.array([[0., -1., 0.], [0., 1., 0.]])
    else:
        x = np.zeros((4, 3))
        y = np.ones((4, 3))
    return x, y


def _svd_reference(x, y, weights):
    u, _, vt = np.linalg.svd(x.T @ (weights[:, None] * y))
    rotation = u @ np.diag([1., 1., np.linalg.det(u @ vt)]) @ vt
    return np.sqrt(np.average(
        np.sum((x @ rotation - y)**2, axis=1), weights=weights))


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("case", ["945", "1747", "perpendicular", "collapsed"])
@pytest.mark.parametrize("scale", [1e-3, 1., 1e3])
@pytest.mark.parametrize("weighted", [False, True])
@pytest.mark.parametrize("with_rotation", [False, True])
def test_collinear_qcp(dtype, case, scale, weighted, with_rotation):
    x, y = _issue_coordinates(case)
    weights = np.linspace(1., 2., len(x)) if weighted else np.ones(len(x))
    weights /= weights.mean()
    x -= np.average(x, axis=0, weights=weights)
    y -= np.average(y, axis=0, weights=weights)
    x = np.ascontiguousarray(x * scale, dtype=dtype)
    y = np.ascontiguousarray(y * scale, dtype=dtype)
    low_weights = weights.astype(dtype) if weighted else None
    w = np.ones(len(x)) if low_weights is None else low_weights.astype(float)
    xd, yd = x.astype(float), y.astype(float)
    expected = _svd_reference(xd, yd, w)
    coordinate_scale = np.sqrt(np.average(
        np.sum(xd*xd + yd*yd, axis=1), weights=w) / 2)
    rot = np.zeros(9, dtype=dtype) if with_rotation else None
    result = qcprot.CalcRMSDRotationalMatrix(x, y, len(x), rot, low_weights)
    assert np.isfinite(result)
    # Accumulation and the E0-lambda subtraction limit RMSD accuracy near zero.
    squared_atol = 32 * len(x) * np.finfo(dtype).eps * coordinate_scale**2
    assert abs(result**2 - expected**2) <= squared_atol + 1e-24
    if with_rotation:
        matrix = rot.astype(float).reshape(3, 3)
        atol = 1e-6 if dtype is np.float32 else 1e-10
        assert_allclose(matrix @ matrix.T, np.eye(3), atol=atol, rtol=0)
        assert np.linalg.det(matrix) == pytest.approx(1., abs=atol)
        achieved = np.sqrt(np.average(
            np.sum((xd @ matrix.T - yd)**2, axis=1), weights=w))
        alignment_atol = (1e-6 if dtype is np.float32 else 1e-7)
        assert abs(achieved - expected) <= alignment_atol * coordinate_scale + 1e-12


@pytest.mark.parametrize("dtype", [np.float32, np.float64])
@pytest.mark.parametrize("case", ["945", "1747", "perpendicular"])
def test_collinear_public_api(dtype, case):
    x, y = _issue_coordinates(case)
    x, y = x.astype(dtype), y.astype(dtype)
    xc = x.astype(float) - x.astype(float).mean(axis=0)
    yc = y.astype(float) - y.astype(float).mean(axis=0)
    expected = _svd_reference(xc, yc, np.ones(len(x)))
    result = rms.rmsd(x, y, superposition=True)
    scale = np.sqrt(np.mean(np.sum(xc*xc + yc*yc, axis=1)) / 2)
    assert result == pytest.approx(expected, abs=5e-8 * scale)
    rotation, value = align.rotation_matrix(xc, yc)
    achieved = np.sqrt(np.mean(np.sum((xc @ rotation.T - yc)**2, axis=1)))
    assert achieved == pytest.approx(expected, abs=5e-8 * scale)
    assert value == pytest.approx(expected, abs=5e-8 * scale)


def test_scalar_fallback_skips_eigenvectors(monkeypatch):
    x, y = _issue_coordinates("perpendicular")
    original = np.linalg.eigvalsh
    calls = []

    def eigenvalues(matrix):
        calls.append(matrix.shape)
        return original(matrix)

    def unexpected_eigenvectors(matrix):
        pytest.fail("RMSD-only fallback should not compute eigenvectors")

    monkeypatch.setattr(np.linalg, "eigvalsh", eigenvalues)
    monkeypatch.setattr(np.linalg, "eigh", unexpected_eigenvectors)
    result = qcprot.CalcRMSDRotationalMatrix(x, y, len(x), None, None)
    assert result == pytest.approx(0., abs=1e-12)
    assert calls == [(4, 4)]
