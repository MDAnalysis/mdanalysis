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
from itertools import permutations

import numpy as np
import pytest
from numpy.testing import (assert_allclose, assert_equal, assert_almost_equal,
                           assert_array_equal)

from MDAnalysis.lib import transformations as t


"""
Testing transformations is weird because there are 2 versions of many of
these functions.  This is because both python and Cython versions of
these functions exist.  To test therefore, each test has to be done twice,
once for each backend.  This is done through parametrizing and passing
in both versions of the function as an argument.

This should ensure that both versions work and are covered!

.. versionchanged:: 1.0.0
   test_transformations_old_module was removed as core/transformations.py is
   gone 
"""

# tolerance for tests
_ATOL = 1e-06


@pytest.mark.parametrize('f', [
    t._py_identity_matrix,
    t.identity_matrix
])
def test_identity_matrix(f):
    I = f()
    assert_allclose(I, np.dot(I, I))
    assert_equal(np.sum(I), np.trace(I))
    assert_allclose(I, np.identity(4, dtype=np.float64))


@pytest.mark.parametrize('f', [
    t._py_translation_matrix,
    t.translation_matrix,
])
def test_translation_matrix(f):
    v = np.array([0.2, 0.2, 0.2])
    assert_allclose(v, f(v)[:3, 3])


def test_translation_from_matrix():
    # doesn't seem to have a Cython backend
    v0 = np.array([0.2, 0.2, 0.2])
    v1 = t.translation_from_matrix(t.translation_matrix(v0))
    assert_allclose(v0, v1)


@pytest.mark.parametrize('f', [
    t._py_reflection_matrix,
    t.reflection_matrix,
])
def test_reflection_matrix(f):
    v0 = np.array([0.2, 0.2, 0.2, 1.0])  # arbitrary values
    v1 = np.array([0.4, 0.4, 0.4])
    R = f(v0, v1)
    assert_allclose(2., np.trace(R))
    assert_allclose(v0, np.dot(R, v0))
    v2 = v0.copy()
    v2[:3] += v1
    v3 = v0.copy()
    v2[:3] -= v1
    assert_allclose(v2, np.dot(R, v3))


def test_reflection_from_matrix():
    v0 = np.array([0.2, 0.2, 0.2])  # arbitrary values
    v1 = np.array([0.4, 0.4, 0.4])
    M0 = t.reflection_matrix(v0, v1)
    point, normal = t.reflection_from_matrix(M0)
    M1 = t.reflection_matrix(point, normal)
    assert_equal(t.is_same_transform(M0, M1), True)


@pytest.mark.parametrize('f', [
    t._py_rotation_matrix,
    t.rotation_matrix,
])
def test_rotation_matrix(f):
    R = f(np.pi / 2.0, [0, 0, 1], [1, 0, 0])
    assert_allclose(np.dot(R, [0, 0, 0, 1]), [1., -1., 0., 1.])
    angle = 0.2 * 2 * np.pi  # arbitrary value
    direc = np.array([0.2, 0.2, 0.2])
    point = np.array([0.4, 0.4, 0.4])
    R0 = f(angle, direc, point)
    R1 = f(angle - 2 * np.pi, direc, point)
    assert_equal(t.is_same_transform(R0, R1), True)
    R0 = f(angle, direc, point)
    R1 = f(-angle, -direc, point)
    assert_equal(t.is_same_transform(R0, R1), True)
    I = np.identity(4, np.float64)
    assert_allclose(I, f(np.pi * 2, direc), atol=_ATOL)
    assert_allclose(2., np.trace(f(np.pi / 2, direc, point)))


def test_rotation_from_matrix():
    angle = 0.2 * 2 * np.pi  # arbitrary values
    direc = np.array([0.2, 0.2, 0.2])
    point = np.array([0.4, 0.4, 0.4])
    R0 = t.rotation_matrix(angle, direc, point)
    angle, direc, point = t.rotation_from_matrix(R0)
    R1 = t.rotation_matrix(angle, direc, point)
    assert_equal(t.is_same_transform(R0, R1), True)


@pytest.mark.parametrize('f', [
    t._py_scale_matrix,
    t.scale_matrix,
])
def test_scale_matrix(f):
    v = np.array([14.1, 15.1, 16.1, 1])
    S = f(-1.234)
    assert_allclose(np.dot(S, v)[:3], -1.234 * v[:3])


def test_scale_from_matrix():
    factor = 7
    origin = np.array([0.2, 0.2, 0.2])  # arbitrary values
    direct = np.array([0.4, 0.4, 0.4])
    S0 = t.scale_matrix(factor, origin)
    factor, origin, direction = t.scale_from_matrix(S0)
    S1 = t.scale_matrix(factor, origin, direction)
    assert_equal(t.is_same_transform(S0, S1), True)
    S0 = t.scale_matrix(factor, origin, direct)
    factor, origin, direction = t.scale_from_matrix(S0)
    S1 = t.scale_matrix(factor, origin, direction)
    assert_equal(t.is_same_transform(S0, S1), True)

@pytest.mark.parametrize('f', [
    t._py_projection_matrix,
    t.projection_matrix,
])
class TestProjectionMatrix(object):
    def test_projection_matrix_1(self, f):
        P = f((0, 0, 0), (1, 0, 0))
        assert_allclose(P[1:, 1:], np.identity(4)[1:, 1:], atol=_ATOL)

    def test_projection_matrix_2(self, f):
        point = np.array([0.2, 0.2, 0.2])  # arbitrary values
        normal = np.array([0.4, 0.4, 0.4])
        direct = np.array([0.6, 0.6, 0.6])
        persp = np.array([0.8, 0.8, 0.8])

        P0 = f(point, normal)
        # TODO: why isn't this used anymore?
        P1 = f(point, normal, direction=direct)
        P2 = f(point, normal, perspective=persp)
        P3 = f(point, normal, perspective=persp, pseudo=True)
        assert_equal(t.is_same_transform(P2, np.dot(P0, P3)), True)

    def test_projection_matrix_3(self, f):
        P = f((3, 0, 0), (1, 1, 0), (1, 0, 0))
        v0 = np.array([14.1, 15.1, 16.1, 1])  # arbitrary values
        v1 = np.dot(P, v0)
        assert_allclose(v1[1], v0[1], atol=_ATOL)
        assert_allclose(v1[0], 3.0 - v1[1], atol=_ATOL)



class TestProjectionFromMatrix(object):
    @staticmethod
    @pytest.fixture()
    def data():
        point = np.array([0.2, 0.2, 0.2])  # arbitrary values
        normal = np.array([0.4, 0.4, 0.4])
        direct = np.array([0.6, 0.6, 0.6])
        persp = np.array([0.8, 0.8, 0.8])
        return point, normal, direct, persp

    def test_projection_from_matrix_1(self, data):
        point, normal, direct, persp = data
        P0 = t.projection_matrix(point, normal)
        result = t.projection_from_matrix(P0)
        P1 = t.projection_matrix(*result)
        assert_equal(t.is_same_transform(P0, P1), True)

    def test_projection_from_matrix_2(self, data):
        point, normal, direct, persp = data
        P0 = t.projection_matrix(point, normal, direct)
        result = t.projection_from_matrix(P0)
        P1 = t.projection_matrix(*result)
        assert_equal(t.is_same_transform(P0, P1), True)

    def test_projection_from_matrix_3(self, data):
        point, normal, direct, persp = data
        P0 = t.projection_matrix(
            point, normal, perspective=persp, pseudo=False)
        result = t.projection_from_matrix(P0, pseudo=False)
        P1 = t.projection_matrix(*result)
        assert_equal(t.is_same_transform(P0, P1), True)

    def test_projection_from_matrix_4(self, data):
        point, normal, direct, persp = data
        P0 = t.projection_matrix(
            point, normal, perspective=persp, pseudo=True)
        result = t.projection_from_matrix(P0, pseudo=True)
        P1 = t.projection_matrix(*result)
        assert_equal(t.is_same_transform(P0, P1), True)



@pytest.mark.parametrize('f', [
    t._py_clip_matrix,
    t.clip_matrix,
])
class TestClipMatrix(object):
    def test_clip_matrix_1(self, f):
        frustrum = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])  # arbitrary values
        frustrum[1] += frustrum[0]
        frustrum[3] += frustrum[2]
        frustrum[5] += frustrum[4]
        M = f(perspective=False, *frustrum)
        assert_allclose(
            np.dot(M, [frustrum[0], frustrum[2], frustrum[4], 1.0]),
            np.array([-1., -1., -1., 1.]))
        assert_allclose(
            np.dot(M, [frustrum[1], frustrum[3], frustrum[5], 1.0]),
            np.array([1., 1., 1., 1.]))

    def test_clip_matrix_2(self, f):
        frustrum = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])  # arbitrary values
        frustrum[1] += frustrum[0]
        frustrum[3] += frustrum[2]
        frustrum[5] += frustrum[4]
        M = f(perspective=True, *frustrum)
        v = np.dot(M, [frustrum[0], frustrum[2], frustrum[4], 1.0])
        assert_allclose(v / v[3], np.array([-1., -1., -1., 1.]))
        v = np.dot(M, [frustrum[1], frustrum[3], frustrum[4], 1.0])
        assert_allclose(v / v[3], np.array([1., 1., -1., 1.]))

    def test_clip_matrix_frustrum_left_right_bounds(self, f):
        '''ValueError should be raised if left > right.'''
        frustrum = np.array([0.4, 0.3, 0.3, 0.7, 0.5, 1.1])
        with pytest.raises(ValueError):
            f(*frustrum)

    def test_clip_matrix_frustrum_bottom_top_bounds(self, f):
        '''ValueError should be raised if bottom > top.'''
        frustrum = np.array([0.1, 0.3, 0.71, 0.7, 0.5, 1.1])
        with pytest.raises(ValueError):
            f(*frustrum)

    def test_clip_matrix_frustrum_near_far_bounds(self, f):
        '''ValueError should be raised if near > far.'''
        frustrum = np.array([0.1, 0.3, 0.3, 0.7, 1.5, 1.1])
        with pytest.raises(ValueError):
            f(*frustrum)


@pytest.mark.parametrize('f', [
    t._py_shear_matrix,
    t.shear_matrix,
])
def test_shear_matrix(f):
    angle = 0.2 * 4 * np.pi  # arbitrary values
    direct = np.array([0.2, 0.2, 0.2])
    point = np.array([0.3, 0.4, 0.5])
    normal = np.cross(direct, np.array([0.8, 0.6, 0.4]))
    S = f(angle, direct, point, normal)
    assert_allclose(1.0, np.linalg.det(S), atol=_ATOL)


def test_shear_from_matrix():
    # This seems to fail sometimes if the random numbers
    # roll certain values....
    # angle = (random.random() - 0.5) * 4*np.pi
    # direct = np.random.random(3) - 0.5
    # point = np.random.random(3) - 0.5
    # normal = np.cross(direct, np.random.random(3))
    # In this random configuration the test will fail about 0.05% of all times.
    # Then we hit some edge-cases of the algorithm. The edge cases for these
    # values are slightly different for the linalg library used (MKL/LAPACK).
    # So here are some of my random numbers
    angle = 2.8969075413405783  # arbitrary values
    direct = np.array([-0.31117458, -0.41769518, -0.01188556])
    point = np.array([-0.0035982, -0.40997482, 0.42241425])
    normal = np.cross(direct, np.array([0.08122421, 0.4747914, 0.19851859]))

    S0 = t.shear_matrix(angle, direct, point, normal)
    angle, direct, point, normal = t.shear_from_matrix(S0)
    S1 = t.shear_matrix(angle, direct, point, normal)
    assert_equal(t.is_same_transform(S0, S1), True)


class TestDecomposeMatrix(object):
    def test_decompose_matrix_1(self):
        T0 = t.translation_matrix((1, 2, 3))
        scale, shear, angles, trans, persp = t.decompose_matrix(T0)
        T1 = t.translation_matrix(trans)
        assert_allclose(T0, T1)

    def test_decompose_matrix_2(self):
        S = t.scale_matrix(0.123)
        scale, shear, angles, trans, persp = t.decompose_matrix(S)
        assert_equal(scale[0], 0.123)

    def test_decompose_matrix_3(self):
        R0 = t.euler_matrix(1, 2, 3)
        scale, shear, angles, trans, persp = t.decompose_matrix(R0)
        R1 = t.euler_matrix(*angles)
        assert_allclose(R0, R1)


def test_compose_matrix():
    scale = np.array([0.2, 0.2, 0.2])  # arbitrary values
    shear = np.array([0.4, 0.4, 0.4])
    angles = np.array([0.6, 0.6, 0.6]) * 2 * np.pi
    trans = np.array([0.8, 0.8, 0.8])
    persp = np.array([0.9, 0.9, 0.9, 0.9])

    M0 = t.compose_matrix(scale, shear, angles, trans, persp)
    result = t.decompose_matrix(M0)
    M1 = t.compose_matrix(*result)
    assert_equal(t.is_same_transform(M0, M1), True)


@pytest.mark.parametrize('f', [
    t._py_orthogonalization_matrix,
    t.orthogonalization_matrix,
])
class TestOrthogonalizationMatrix(object):
    def test_orthogonalization_matrix_1(self, f):
        O = f((10., 10., 10.), (90., 90., 90.))
        assert_allclose(O[:3, :3], np.identity(3, float) * 10, atol=_ATOL)

    def test_orthogonalization_matrix_2(self, f):
        O = f([9.8, 12.0, 15.5], [87.2, 80.7, 69.7])
        assert_allclose(np.sum(O), 43.063229, atol=_ATOL)


@pytest.mark.parametrize('f', [
    t._py_superimposition_matrix,
    t.superimposition_matrix,
])
def test_superimposition_matrix(f):
    v0 = np.sin(np.linspace(0, 0.99, 30)).reshape(3, 10)  # arbitrary values
    M = f(v0, v0)
    assert_allclose(M, np.identity(4), atol=_ATOL)

    R = t.random_rotation_matrix(np.array([0.3, 0.4, 0.5]))
    v0 = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1))
    v1 = np.dot(R, v0)
    M = f(v0, v1)
    assert_allclose(v1, np.dot(M, v0), atol=_ATOL)

    v0 = np.sin(np.linspace(-1, 1, 400)).reshape(4, 100)
    v0[3] = 1.0
    v1 = np.dot(R, v0)
    M = f(v0, v1)
    assert_allclose(v1, np.dot(M, v0), atol=_ATOL)

    S = t.scale_matrix(0.45)
    T = t.translation_matrix(np.array([0.2, 0.2, 0.2]) - 0.5)
    M = t.concatenate_matrices(T, R, S)
    v1 = np.dot(M, v0)
    v0[:3] += np.sin(np.linspace(0.0, 1e-9, 300)).reshape(3, -1)
    M = f(v0, v1, scaling=True)
    assert_allclose(v1, np.dot(M, v0), atol=_ATOL)

    M = f(v0, v1, scaling=True, usesvd=False)
    assert_allclose(v1, np.dot(M, v0), atol=_ATOL)

    v = np.empty((4, 100, 3), dtype=np.float64)
    v[:, :, 0] = v0
    M = f(v0, v1, scaling=True, usesvd=False)
    assert_allclose(v1, np.dot(M, v[:, :, 0]), atol=_ATOL)


@pytest.mark.parametrize('f', [
    t._py_euler_matrix,
    t.euler_matrix,
])
class TestEulerMatrix(object):
    def test_euler_matrix_1(self, f):
        R = f(1, 2, 3, 'syxz')
        assert_allclose(np.sum(R[0]), -1.34786452)

    def test_euler_matrix_2(self, f):
        R = f(1, 2, 3, (0, 1, 0, 1))
        assert_allclose(np.sum(R[0]), -0.383436184)


@pytest.mark.parametrize('f', [
    t._py_euler_from_matrix,
    t.euler_from_matrix,
])
class TestEulerFromMatrix(object):
    def test_euler_from_matrix_1(self, f):
        R0 = t.euler_matrix(1, 2, 3, 'syxz')
        al, be, ga = f(R0, 'syxz')
        R1 = t.euler_matrix(al, be, ga, 'syxz')
        assert_allclose(R0, R1)

    def test_euler_from_matrix_2(self, f):
        angles = 4.0 * np.pi * np.array([-0.3, -0.3, -0.3])  # arbitrary values
        for axes in t._AXES2TUPLE.keys():
            R0 = t.euler_matrix(axes=axes, *angles)
            R1 = t.euler_matrix(axes=axes, *f(R0, axes))
            assert_allclose(R0, R1, err_msg=("{0} failed".format(axes)))


def test_euler_from_quaternion():
    angles = t.euler_from_quaternion([0.99810947, 0.06146124, 0, 0])
    assert_allclose(angles, [0.123, 0, 0], atol=_ATOL)


@pytest.mark.parametrize('f', [
    t._py_quaternion_from_euler,
    t.quaternion_from_euler,
])
def test_quaternion_from_euler(f):
    q = f(1, 2, 3, 'ryxz')
    assert_allclose(q, [0.435953, 0.310622, -0.718287, 0.444435], atol=_ATOL)


@pytest.mark.parametrize('f', [
    t._py_quaternion_about_axis,
    t.quaternion_about_axis,
])
def test_quaternion_about_axis(f):
    q = f(0.123, (1, 0, 0))
    assert_allclose(q, [0.99810947, 0.06146124, 0, 0], atol=_ATOL)


@pytest.mark.parametrize('f', [
    t._py_quaternion_matrix,
    t.quaternion_matrix,
])
class TestQuaternionMatrix(object):
    def test_quaternion_matrix_1(self, f):
        M = f([0.99810947, 0.06146124, 0, 0])
        assert_allclose(M, t.rotation_matrix(0.123, (1, 0, 0)), atol=_ATOL)

    def test_quaternion_matrix_2(self, f):
        M = f([1, 0, 0, 0])
        assert_allclose(M, t.identity_matrix(), atol=_ATOL)

    def test_quaternion_matrix_3(self, f):
        M = f([0, 1, 0, 0])
        assert_allclose(M, np.diag([1, -1, -1, 1]), atol=_ATOL)


@pytest.mark.parametrize('f', [
    t._py_quaternion_from_matrix,
    t.quaternion_from_matrix,
])
class TestQuaternionFromMatrix(object):
    def test_quaternion_from_matrix_1(self, f):
        q = f(t.identity_matrix(), True)
        assert_allclose(q, [1., 0., 0., 0.], atol=_ATOL)

    def test_quaternion_from_matrix_2(self, f):
        q = f(np.diag([1., -1., -1., 1.]))
        check = (np.allclose(
            q, [0, 1, 0, 0], atol=_ATOL) or np.allclose(
                q, [0, -1, 0, 0], atol=_ATOL))
        assert_equal(check, True)

    def test_quaternion_from_matrix_3(self, f):
        R = t.rotation_matrix(0.123, (1, 2, 3))
        q = f(R, True)
        assert_allclose(
            q, [0.9981095, 0.0164262, 0.0328524, 0.0492786], atol=_ATOL)

    def test_quaternion_from_matrix_4(self, f):
        R = [[-0.545, 0.797, 0.260, 0], [0.733, 0.603, -0.313, 0],
             [-0.407, 0.021, -0.913, 0], [0, 0, 0, 1]]
        q = f(R)
        assert_allclose(q, [0.19069, 0.43736, 0.87485, -0.083611], atol=_ATOL)

    def test_quaternion_from_matrix_5(self, f):
        R = [[0.395, 0.362, 0.843, 0], [-0.626, 0.796, -0.056, 0],
             [-0.677, -0.498, 0.529, 0], [0, 0, 0, 1]]
        q = f(R)
        assert_allclose(
            q, [0.82336615, -0.13610694, 0.46344705, -0.29792603], atol=_ATOL)

    def test_quaternion_from_matrix_6(self, f):
        R = t.random_rotation_matrix()
        q = f(R)
        assert_equal(t.is_same_transform(R, t.quaternion_matrix(q)), True)


@pytest.mark.parametrize('f', [
    t._py_quaternion_multiply,
    t.quaternion_multiply,
])
def test_quaternion_multiply(f):
    q = f([4, 1, -2, 3], [8, -5, 6, 7])
    assert_allclose(q, [28, -44, -14, 48])


@pytest.mark.parametrize('f', [
    t._py_quaternion_conjugate,
    t.quaternion_conjugate,
])
def test_quaternion_conjugate(f):
    q0 = t.random_quaternion()
    q1 = f(q0)
    check = q1[0] == q0[0] and all(q1[1:] == -q0[1:])
    assert_equal(check, True)


@pytest.mark.parametrize('f', [
    t._py_quaternion_inverse,
    t.quaternion_inverse,
])
def test_quaternion_inverse(f):
    q0 = t.random_quaternion()
    q1 = f(q0)
    assert_allclose(t.quaternion_multiply(q0, q1), [1, 0, 0, 0], atol=_ATOL)


def test_quaternion_real():
    assert_allclose(t.quaternion_real([3.0, 0.0, 1.0, 2.0]), 3.0)


def test_quaternion_imag():
    assert_allclose(t.quaternion_imag([3.0, 0.0, 1.0, 2.0]), [0.0, 1.0, 2.0])


@pytest.mark.parametrize('f', [
    t._py_quaternion_slerp,
    t.quaternion_slerp,
])
def test_quaternion_slerp(f):
    q0 = t.random_quaternion()
    q1 = t.random_quaternion()
    q = f(q0, q1, 0.0)
    assert_allclose(q, q0, atol=_ATOL)

    q = f(q0, q1, 1.0, 1)
    assert_allclose(q, q1, atol=_ATOL)

    q = f(q0, q1, 0.5)
    angle = np.arccos(np.dot(q0, q))

    check = (np.allclose(2.0, np.arccos(np.dot(q0, q1)) / angle) or
             np.allclose(2.0, np.arccos(-np.dot(q0, q1)) / angle))

    assert_equal(check, True)


@pytest.mark.parametrize('f', [
    t._py_random_quaternion,
    t.random_quaternion,
])
class TestRandomQuaternion(object):
    def test_random_quaternion_1(self, f):
        q = f()
        assert_allclose(1.0, t.vector_norm(q))

    def test_random_quaternion_2(self, f):
        q = f(np.array([0.2, 0.2, 0.2]))
        assert_equal(len(q.shape), 1)
        assert_equal(q.shape[0] == 4, True)


@pytest.mark.parametrize('f', [
    t._py_random_rotation_matrix,
    t.random_rotation_matrix,
])
def test_random_rotation_matrix(f):
    R = f()
    assert_allclose(np.dot(R.T, R), np.identity(4), atol=_ATOL)


@pytest.mark.parametrize('f', [
    t._py_inverse_matrix,
    t.inverse_matrix,
])
class TestInverseMatrix(object):
    @pytest.mark.parametrize('size', list(range(1, 7)))
    def test_inverse(self, size, f):
        # Create a known random state to generate numbers from
        # these numbers will then be uncorrelated but deterministic
        rs = np.random.RandomState(1234)
        M0 = rs.randn(size, size)
        M1 = f(M0)
        assert_allclose(M1, np.linalg.inv(M0), err_msg=str(size), atol=_ATOL)

    def test_inverse_matrix(self, f):
        M0 = t.random_rotation_matrix()
        M1 = f(M0.T)
        assert_allclose(M1, np.linalg.inv(M0.T))


@pytest.mark.parametrize('f', [
    t._py_is_same_transform,
    t.is_same_transform,
])
class TestIsSameTransform(object):
    def test_is_same_transform_1(self, f):
        assert_equal(f(np.identity(4), np.identity(4)), True)

    def test_is_same_transform_2(self, f):
        assert_equal(f(t.random_rotation_matrix(), np.identity(4)), False)


@pytest.mark.parametrize('f', [
    t._py_random_vector,
    t.random_vector,
])
class TestRandomVector(object):
    def test_random_vector_1(self, f):
        v = f(1000)
        check = np.all(v >= 0.0) and np.all(v < 1.0)
        assert_equal(check, True)

    def test_random_vector_2(self, f):
        v0 = f(10)
        v1 = f(10)
        assert_equal(np.any(v0 == v1), False)


@pytest.mark.parametrize('f', [
    t._py_unit_vector,
    t.unit_vector,
])
class TestUnitVector(object):
    def test_unit_vector_1(self, f):
        v0 = np.array([0.2, 0.2, 0.2])
        v1 = f(v0)
        assert_allclose(v1, v0 / np.linalg.norm(v0), atol=_ATOL)

    def test_unit_vector_2(self, f):
        v0 = np.sin(np.linspace(0, 10, 5 * 4 * 3)).reshape(5, 4, 3)
        v1 = f(v0, axis=-1)
        v2 = v0 / np.expand_dims(np.sqrt(np.sum(v0 * v0, axis=2)), 2)
        assert_allclose(v1, v2, atol=_ATOL)

    def test_unit_vector_3(self, f):
        v0 = np.sin(np.linspace(0, 10, 5 * 4 * 3)).reshape(5, 4, 3)
        v1 = f(v0, axis=1)
        v2 = v0 / np.expand_dims(np.sqrt(np.sum(v0 * v0, axis=1)), 1)
        assert_allclose(v1, v2, atol=_ATOL)

    def test_unit_vector_4(self, f):
        v0 = np.sin(np.linspace(0, 10, 5 * 4 * 3)).reshape(5, 4, 3)
        v1 = np.empty((5, 4, 3), dtype=np.float64)
        v2 = v0 / np.expand_dims(np.sqrt(np.sum(v0 * v0, axis=1)), 1)
        f(v0, axis=1, out=v1)
        assert_allclose(v1, v2, atol=_ATOL)

    def test_unit_vector_5(self, f):
        assert_equal(list(f([])), [])

    def test_unit_vector_6(self, f):
        assert_equal(list(f([1.0])), [1.0])


@pytest.mark.parametrize('f', [
    t._py_vector_norm,
    t.vector_norm,
])
class TestVectorNorm(object):
    def test_vector_norm_1(self, f):
        v = np.array([0.2, 0.2, 0.2])
        n = f(v)
        assert_allclose(n, np.linalg.norm(v), atol=_ATOL)

    def test_vector_norm_2(self, f):
        v = np.sin(np.linspace(0, 10, 6 * 5 * 3)).reshape(6, 5, 3)
        n = f(v, axis=-1)
        assert_allclose(n, np.sqrt(np.sum(v * v, axis=2)), atol=_ATOL)

    def test_vector_norm_3(self, f):
        v = np.sin(np.linspace(0, 10, 6 * 5 * 3)).reshape(6, 5, 3)
        n = f(v, axis=1)
        assert_allclose(n, np.sqrt(np.sum(v * v, axis=1)), atol=_ATOL)

    def test_vector_norm_4(self, f):
        v = np.sin(np.linspace(0, 10, 5 * 4 * 3)).reshape(5, 4, 3)
        n = np.empty((5, 3), dtype=np.float64)
        f(v, axis=1, out=n)
        assert_allclose(n, np.sqrt(np.sum(v * v, axis=1)), atol=_ATOL)

    def test_vector_norm_5(self, f):
        assert_equal(f([]), 0.0)

    def test_vector_norm_6(self, f):
        assert_equal(f([1.0]), 1.0)


class TestArcBall(object):
    def test_arcball_1(self):
        ball = t.Arcball()
        ball = t.Arcball(initial=np.identity(4))
        ball.place([320, 320], 320)
        ball.down([500, 250])
        ball.drag([475, 275])
        R = ball.matrix()
        assert_allclose(np.sum(R), 3.90583455, atol=_ATOL)

    def test_arcball_2(self):
        ball = t.Arcball(initial=[1, 0, 0, 0])
        ball.place([320, 320], 320)
        ball.setaxes([1, 1, 0], [-1, 1, 0])
        ball.setconstrain(True)
        ball.down([400, 200])
        ball.drag([200, 400])
        R = ball.matrix()
        assert_allclose(np.sum(R), 0.2055924)


def test_rotaxis_equal_vectors():
    a = np.arange(3)
    x = t.rotaxis(a, a)
    assert_array_equal(x, [1, 0, 0])


def test_rotaxis_different_vectors():
    # use random coordinate system
    e = np.eye(3)
    r = np.array([[0.69884766, 0.59804425, -0.39237102],
                  [0.18784672, 0.37585347, 0.90744023],
                  [0.69016342, -0.7078681, 0.15032367]])
    re = np.dot(r, e)

    for i, j, l in permutations(range(3)):
        x = t.rotaxis(re[i], re[j])
        # use abs since direction doesn't matter
        assert_almost_equal(np.abs(np.dot(x, re[l])), 1)
