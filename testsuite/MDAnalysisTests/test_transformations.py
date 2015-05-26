# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

import numpy as np
from numpy.testing import *
import math
import random

from MDAnalysis.core import transformations as t

# tolerance for tests
_ATOL = 1e-06

def test_identity_matrix():
    I = t.identity_matrix()
    assert_allclose(I, np.dot(I,I))
    assert_equal(np.sum(I), np.trace(I))
    assert_allclose(I, np.identity(4, dtype=np.float64))

def test_translation_matrix():
    v = np.random.random(3) - 0.5
    assert_allclose(v, t.translation_matrix(v)[:3, 3])

def test_translation_from_matrix():
    v0 = np.random.random(3) - 0.5
    v1 = t.translation_from_matrix(t.translation_matrix(v0))
    assert_allclose(v0, v1)

def test_reflection_matrix():
    v0 = np.random.random(4) - 0.5
    v0[3] = 1.0
    v1 = np.random.random(3) - 0.5
    R = t.reflection_matrix(v0, v1)
    assert_allclose(2., np.trace(R))
    assert_allclose(v0, np.dot(R, v0))
    v2 = v0.copy()
    v2[:3] += v1
    v3 = v0.copy()
    v2[:3] -= v1
    assert_allclose(v2, np.dot(R, v3))

def test_reflection_from_matrix():
    v0 = np.random.random(3) - 0.5
    v1 = np.random.random(3) - 0.5
    M0 = t.reflection_matrix(v0, v1)
    point, normal = t.reflection_from_matrix(M0)
    M1 = t.reflection_matrix(point, normal)
    assert_equal(t.is_same_transform(M0, M1), True)

def test_rotation_matrix():
    R = t.rotation_matrix(np.pi/2.0, [0, 0, 1], [1, 0, 0])
    assert_allclose(np.dot(R, [0, 0, 0, 1]), [ 1., -1.,  0.,  1.])
    angle = (np.random.random(1) - 0.5) * (2*np.pi)
    direc = np.random.random(3) - 0.5
    point = np.random.random(3) - 0.5
    R0 = t.rotation_matrix(angle, direc, point)
    R1 = t.rotation_matrix(angle-2*np.pi, direc, point)
    assert_equal(t.is_same_transform(R0, R1), True)
    R0 = t.rotation_matrix(angle, direc, point)
    R1 = t.rotation_matrix(-angle, -direc, point)
    assert_equal(t.is_same_transform(R0, R1), True)
    I = np.identity(4, np.float64)
    assert_allclose(I, t.rotation_matrix(np.pi*2, direc), atol=_ATOL)
    assert_allclose(2., np.trace(t.rotation_matrix(np.pi/2, direc, point)))

def test_rotation_from_matrix():
    angle = (np.random.random(1) - 0.5) * (2*np.pi)
    direc = np.random.random(3) - 0.5
    point = np.random.random(3) - 0.5
    R0 = t.rotation_matrix(angle, direc, point)
    angle, direc, point = t.rotation_from_matrix(R0)
    R1 = t.rotation_matrix(angle, direc, point)
    assert_equal(t.is_same_transform(R0, R1), True)

def test_scale_matrix():
    v = (np.random.rand(4, 5) - 0.5) * 20.0
    v[3] = 1.0
    S = t.scale_matrix(-1.234)
    assert_allclose(np.dot(S, v)[:3], -1.234*v[:3])

def test_scale_from_matrix():
    factor = np.random.random(1) * 10 - 5
    origin = np.random.random(3) - 0.5
    direct = np.random.random(3) - 0.5
    S0 = t.scale_matrix(factor, origin)
    factor, origin, direction = t.scale_from_matrix(S0)
    S1 = t.scale_matrix(factor, origin, direction)
    assert_equal(t.is_same_transform(S0, S1), True)
    S0 = t.scale_matrix(factor, origin, direct)
    factor, origin, direction = t.scale_from_matrix(S0)
    S1 = t.scale_matrix(factor, origin, direction)
    assert_equal(t.is_same_transform(S0, S1), True)

class TestProjectionMatrix(TestCase):
    def test_projection_matrix_1(self):
        P = t.projection_matrix((0, 0, 0), (1, 0, 0))
        assert_allclose(P[1:, 1:], np.identity(4)[1:, 1:])

    def test_projection_matrix_2(self):
        point = np.random.random(3) - 0.5
        normal = np.random.random(3) - 0.5
        direct = np.random.random(3) - 0.5
        persp = np.random.random(3) - 0.5
        P0 = t.projection_matrix(point, normal)
        P1 = t.projection_matrix(point, normal, direction=direct)
        P2 = t.projection_matrix(point, normal, perspective=persp)
        P3 = t.projection_matrix(point, normal, perspective=persp, pseudo=True)
        assert_equal(t.is_same_transform(P2, np.dot(P0, P3)), True)

    def test_projection_matrix_3(self):
        P = t.projection_matrix((3, 0, 0), (1, 1, 0), (1, 0, 0))
        v0 = (np.random.rand(4, 5) - 0.5) * 20.0
        v0[3] = 1.0
        v1 = np.dot(P, v0)
        assert_allclose(v1[1], v0[1])
        assert_allclose(v1[0], 3.0-v1[1])


class TestProjectionFromMatrix(TestCase):
    def setUp(self):
        self.point = np.random.random(3) - 0.5
        self.normal = np.random.random(3) - 0.5
        self.direct = np.random.random(3) - 0.5
        self.persp = np.random.random(3) - 0.5

    def test_projection_from_matrix_1(self):
        P0 = t.projection_matrix(self.point, self.normal)
        result = t.projection_from_matrix(P0)
        P1 = t.projection_matrix(*result)
        assert_equal(t.is_same_transform(P0, P1), True)

    def test_projection_from_matrix_2(self):
        P0 = t.projection_matrix(self.point, self.normal, self.direct)
        result = t.projection_from_matrix(P0)
        P1 = t.projection_matrix(*result)
        assert_equal(t.is_same_transform(P0, P1), True)

    def test_projection_from_matrix_3(self):
        P0 = t.projection_matrix(self.point, self.normal,
                                 perspective=self.persp, pseudo=False)
        result = t.projection_from_matrix(P0, pseudo=False)
        P1 = t.projection_matrix(*result)
        assert_equal(t.is_same_transform(P0, P1), True)

    def test_projection_from_matrix_4(self):
        P0 = t.projection_matrix(self.point, self.normal,
                                 perspective=self.persp, pseudo=True)
        result = t.projection_from_matrix(P0, pseudo=True)
        P1 = t.projection_matrix(*result)
        assert_equal(t.is_same_transform(P0, P1), True)


class TestClipMatrix(TestCase):
    def test_clip_matrix_1(self):
        frustrum = np.random.rand(6)
        frustrum[1] += frustrum[0]
        frustrum[3] += frustrum[2]
        frustrum[5] += frustrum[4]
        M = t.clip_matrix(perspective=False, *frustrum)
        assert_allclose(np.dot(M, [frustrum[0], frustrum[2], frustrum[4], 1.0]),
                        np.array([-1., -1., -1.,  1.]))
        assert_allclose(np.dot(M, [frustrum[1], frustrum[3], frustrum[5], 1.0]),
                         np.array([ 1.,  1.,  1.,  1.]))

    def test_clip_matrix_2(self):
        frustrum = np.random.rand(6)
        frustrum[1] += frustrum[0]
        frustrum[3] += frustrum[2]
        frustrum[5] += frustrum[4]
        M = t.clip_matrix(perspective=True, *frustrum)
        v = np.dot(M, [frustrum[0], frustrum[2], frustrum[4], 1.0])
        assert_allclose(v / v[3],
                        np.array([-1., -1., -1.,  1.]))
        v = np.dot(M, [frustrum[1], frustrum[3], frustrum[4], 1.0])
        assert_allclose(v / v[3],
                        np.array([ 1.,  1., -1.,  1.]))

def test_shear_matrix():
    angle = (np.random.random(1) - 0.5) * 4*np.pi
    direct = np.random.random(3) - 0.5
    point = np.random.random(3) - 0.5
    normal = np.cross(direct, np.random.random(3))
    S = t.shear_matrix(angle, direct, point, normal)
    assert_allclose(1.0, np.linalg.det(S))

def test_shear_from_matrix():
    angle = (np.random.random(1) - 0.5) * 4*np.pi
    direct = np.random.random(3) - 0.5
    point = np.random.random(3) - 0.5
    normal = np.cross(direct, np.random.random(3))
    S0 = t.shear_matrix(angle, direct, point, normal)
    angle, direct, point, normal = t.shear_from_matrix(S0)
    S1 = t.shear_matrix(angle, direct, point, normal)
    assert_equal(t.is_same_transform(S0, S1), True)

class TestDecomposeMatrix(TestCase):
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
    scale = np.random.random(3) - 0.5
    shear = np.random.random(3) - 0.5
    angles = (np.random.random(3) - 0.5) * (2*math.pi)
    trans = np.random.random(3) - 0.5
    persp = np.random.random(4) - 0.5
    M0 = t.compose_matrix(scale, shear, angles, trans, persp)
    result = t.decompose_matrix(M0)
    M1 = t.compose_matrix(*result)
    assert_equal(t.is_same_transform(M0, M1), True)

def test_orthogonalization_matrix_1():
    O = t.orthogonalization_matrix((10., 10., 10.), (90., 90., 90.))
    assert_allclose(O[:3, :3], np.identity(3, float) * 10, atol=_ATOL)

def test_orthogonalization_matrix_2():
    O = t.orthogonalization_matrix([9.8, 12.0, 15.5], [87.2, 80.7, 69.7])
    assert_allclose(np.sum(O), 43.063229)

def test_superimposition_matrix():
    v0 = np.random.rand(3, 10)
    M = t.superimposition_matrix(v0, v0)
    assert_allclose(M, np.identity(4))

    R = t.random_rotation_matrix(np.random.random(3))
    v0 = ((1,0,0), (0,1,0), (0,0,1), (1,1,1))
    v1 = np.dot(R, v0)
    M = t.superimposition_matrix(v0, v1)
    assert_allclose(v1, np.dot(M, v0))

    v0 = (np.random.rand(4, 100) - 0.5) * 20.0
    v0[3] = 1.0
    v1 = np.dot(R, v0)
    M = t.superimposition_matrix(v0, v1)
    assert_allclose(v1, np.dot(M, v0))

    S = t.scale_matrix(random.random())
    T = t.translation_matrix(np.random.random(3)-0.5)
    M = t.concatenate_matrices(T, R, S)
    v1 = np.dot(M, v0)
    v0[:3] += np.random.normal(0.0, 1e-9, 300).reshape(3, -1)
    M = t.superimposition_matrix(v0, v1, scaling=True)
    assert_allclose(v1, np.dot(M, v0))

    M = t.superimposition_matrix(v0, v1, scaling=True, usesvd=False)
    assert_allclose(v1, np.dot(M, v0))

    v = np.empty((4, 100, 3), dtype=np.float64)
    v[:, :, 0] = v0
    M = t.superimposition_matrix(v0, v1, scaling=True, usesvd=False)
    assert_allclose(v1, np.dot(M, v[:, :, 0]))


class TestEulerMatrix(TestCase):
    def test_euler_matrix_1(self):
        R = t.euler_matrix(1, 2, 3, 'syxz')
        assert_allclose(np.sum(R[0]), -1.34786452)

    def test_euler_matrix_2(self):
        R = t.euler_matrix(1, 2, 3, (0, 1, 0, 1))
        assert_allclose(np.sum(R[0]), -0.383436184)


def test_euler_from_matrix_1():
    R0 = t.euler_matrix(1, 2, 3, 'syxz')
    al, be, ga = t.euler_from_matrix(R0, 'syxz')
    R1 = t.euler_matrix(al, be, ga, 'syxz')
    assert_allclose(R0, R1)

def test_euler_from_matrix_2():
    angles = (4.0*math.pi) * (np.random.random(3) - 0.5)
    for axes in t._AXES2TUPLE.keys():
        R0 = t.euler_matrix(axes=axes, *angles)
        R1 = t.euler_matrix(axes=axes, *t.euler_from_matrix(R0, axes))
        assert_allclose(R0, R1, err_msg=("{0} failed".format(axes)))

def test_euler_from_quaternion():
    angles = t.euler_from_quaternion([0.99810947, 0.06146124, 0, 0])
    assert_allclose(angles, [0.123, 0, 0], atol=_ATOL)

def test_quaternion_from_euler():
    q = t.quaternion_from_euler(1, 2, 3, 'ryxz')
    assert_allclose(q, [0.435953, 0.310622, -0.718287, 0.444435], atol=_ATOL)

def test_quaternion_about_axis():
    q = t.quaternion_about_axis(0.123, (1, 0, 0))
    assert_allclose(q, [0.99810947, 0.06146124, 0, 0], atol=_ATOL)

class TestQuaternionMatrix(TestCase):
    def test_quaternion_matrix_1(self):
        M = t.quaternion_matrix([0.99810947, 0.06146124, 0, 0])
        assert_allclose(M, t.rotation_matrix(0.123, (1, 0, 0)), atol=_ATOL)

    def test_quaternion_matrix_2(self):
        M = t.quaternion_matrix([1, 0, 0, 0])
        assert_allclose(M, t.identity_matrix(), atol=_ATOL)

    def test_quaternion_matrix_3(self):
        M = t.quaternion_matrix([0, 1, 0, 0])
        assert_allclose(M, np.diag([1, -1, -1, 1]), atol=_ATOL)


class TestQuaternionFromMatrix(TestCase):
    def test_quaternion_from_matrix_1(self):
        q = t.quaternion_from_matrix(t.identity_matrix(), True)
        assert_allclose(q, [1., 0., 0., 0.], atol=_ATOL)

    def test_quaternion_from_matrix_2(self):
        q = t.quaternion_from_matrix(np.diag([1., -1., -1., 1.]))
        check = (np.allclose(q, [0, 1, 0, 0], atol=_ATOL) or
                 np.allclose(q, [0, -1, 0, 0], atol=_ATOL))
        assert_equal(check, True)

    def test_quaternion_from_matrix_3(self):
        R = t.rotation_matrix(0.123, (1, 2, 3))
        q = t.quaternion_from_matrix(R, True)
        assert_allclose(q, [0.9981095, 0.0164262, 0.0328524, 0.0492786], atol=_ATOL)

    def test_quaternion_from_matrix_4(self):
        R = [[-0.545, 0.797, 0.260, 0], [0.733, 0.603, -0.313, 0],
             [-0.407, 0.021, -0.913, 0], [0, 0, 0, 1]]
        q = t.quaternion_from_matrix(R)
        assert_allclose(q, [0.19069, 0.43736, 0.87485, -0.083611], atol=_ATOL)

    def test_quaternion_from_matrix_5(self):
        R = [[0.395, 0.362, 0.843, 0], [-0.626, 0.796, -0.056, 0],
             [-0.677, -0.498, 0.529, 0], [0, 0, 0, 1]]
        q = t.quaternion_from_matrix(R)
        assert_allclose(q, [0.82336615, -0.13610694, 0.46344705, -0.29792603],
                        atol=_ATOL)

    def test_quaternion_from_matrix_6(self):
        R = t.random_rotation_matrix()
        q = t.quaternion_from_matrix(R)
        assert_equal(t.is_same_transform(R, t.quaternion_matrix(q)), True)


def test_quaternion_multiply():
    q = t.quaternion_multiply([4, 1, -2, 3], [8, -5, 6, 7])
    assert_allclose(q, [28, -44, -14, 48])

def test_quaternion_conjugate():
    q0 = t.random_quaternion()
    q1 = t.quaternion_conjugate(q0)
    check = q1[0] == q0[0] and all(q1[1:] == -q0[1:])
    assert_equal(check, True)

def test_quaternion_inverse():
    q0 = t.random_quaternion()
    q1 = t.quaternion_inverse(q0)
    assert_allclose(t.quaternion_multiply(q0, q1), [1, 0, 0, 0], atol=_ATOL)

def test_quaternion_real():
    assert_allclose(t.quaternion_real([3.0, 0.0, 1.0, 2.0]),
                    3.0)

def test_quaternion_imag():
    assert_allclose(t.quaternion_imag([3.0, 0.0, 1.0, 2.0]),
                    [0.0, 1.0, 2.0])

def test_quaternion_slerp():
    q0 = t.random_quaternion()
    q1 = t.random_quaternion()
    q = t.quaternion_slerp(q0, q1, 0.0)
    assert_allclose(q, q0)

    q = t.quaternion_slerp(q0, q1, 1.0, 1)
    assert_allclose(q, q1)

    q = t.quaternion_slerp(q0, q1, 0.5)
    angle = math.acos(np.dot(q0, q))

    check = (np.allclose(2.0, math.acos(np.dot(q0, q1)) / angle) or
             np.allclose(2.0, math.acos(-np.dot(q0, q1)) / angle))

    assert_equal(check, True)

def test_random_quaternion():
    q = t.random_quaternion()
    assert_allclose(1.0, t.vector_norm(q))

def test_random_quaternion_2():
    q = t.random_quaternion(np.random.random(3))
    assert_equal(len(q.shape), 1)
    assert_equal(q.shape[0] == 4, True)

def test_random_rotation_matrix():
    R = t.random_rotation_matrix()
    assert_allclose(np.dot(R.T, R), np.identity(4), atol=_ATOL)

