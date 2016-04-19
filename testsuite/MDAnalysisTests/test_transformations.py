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

from six.moves import range

import numpy as np
from numpy.testing import assert_allclose, assert_equal, TestCase
import math
import random

from MDAnalysis.lib import transformations as t

"""
Testing transformations is weird because there are 2 versions of many of
these functions.  This is because both python and Cython versions of
these functions exist.  To test therefore, each test has to be done twice,
once for each backend.

The general pattern for this is,

1) Create tests which call self.f (the function)

2) Create mixins which define self.f (one of the two backends)

Eg:

class _ClipMatrix(object):
    def test_this(self):
        result = self.f(stuff)
        assert_awesome(me)

class TestClipMatrixNP(_ClipMatrix):
    f = staticmethod(MDAnalysis.lib.transformations._py_clip_matrix)

class TestClipMatrixCY(_ClipMatrix):
    f = staticmethod(MDAnalysis.lib.transformations.clip_matrix)

Note that the function to be tested needs to be defined as a static method!


This should ensure that both versions work and are covered!
"""

# tolerance for tests
_ATOL = 1e-06

class _IdentityMatrix(object):
    def test_identity_matrix(self):
        I = self.f()
        assert_allclose(I, np.dot(I,I))
        assert_equal(np.sum(I), np.trace(I))
        assert_allclose(I, np.identity(4, dtype=np.float64))

class TestIdentityMatrixNP(_IdentityMatrix):
    f = staticmethod(t._py_identity_matrix)

class TestIdentityMatrixCy(_IdentityMatrix):
    f = staticmethod(t.identity_matrix)


class _TranslationMatrix(object):
    def test_translation_matrix(self):
        v = np.random.random(3) - 0.5
        assert_allclose(v, self.f(v)[:3, 3])

class TestTranslationMatrixNP(_TranslationMatrix):
    f = staticmethod(t._py_translation_matrix)

class TestTranslationMatrixCy(_TranslationMatrix):
    f = staticmethod(t.translation_matrix)


def test_translation_from_matrix():
    # doesn't seem to have a Cython backend
    v0 = np.random.random(3) - 0.5
    v1 = t.translation_from_matrix(t.translation_matrix(v0))
    assert_allclose(v0, v1)

class _ReflectionMatrix(object):
    def test_reflection_matrix(self):
        v0 = np.random.random(4) - 0.5
        v0[3] = 1.0
        v1 = np.random.random(3) - 0.5
        R = self.f(v0, v1)
        assert_allclose(2., np.trace(R))
        assert_allclose(v0, np.dot(R, v0))
        v2 = v0.copy()
        v2[:3] += v1
        v3 = v0.copy()
        v2[:3] -= v1
        assert_allclose(v2, np.dot(R, v3))

class TestReflectionMatrixNP(_ReflectionMatrix):
    f = staticmethod(t._py_reflection_matrix)

class TestReflectionMatrixCy(_ReflectionMatrix):
    f = staticmethod(t.reflection_matrix)


def test_reflection_from_matrix():
    v0 = np.random.random(3) - 0.5
    v1 = np.random.random(3) - 0.5
    M0 = t.reflection_matrix(v0, v1)
    point, normal = t.reflection_from_matrix(M0)
    M1 = t.reflection_matrix(point, normal)
    assert_equal(t.is_same_transform(M0, M1), True)

class _RotationMatrix(object):
    def test_rotation_matrix(self):
        R = self.f(np.pi/2.0, [0, 0, 1], [1, 0, 0])
        assert_allclose(np.dot(R, [0, 0, 0, 1]), [ 1., -1.,  0.,  1.])
        angle = (random.random() - 0.5) * (2*np.pi)
        direc = np.random.random(3) - 0.5
        point = np.random.random(3) - 0.5
        R0 = self.f(angle, direc, point)
        R1 = self.f(angle-2*np.pi, direc, point)
        assert_equal(t.is_same_transform(R0, R1), True)
        R0 = self.f(angle, direc, point)
        R1 = self.f(-angle, -direc, point)
        assert_equal(t.is_same_transform(R0, R1), True)
        I = np.identity(4, np.float64)
        assert_allclose(I, self.f(np.pi*2, direc), atol=_ATOL)
        assert_allclose(2., np.trace(self.f(np.pi/2, direc, point)))

class TestRotationMatrixNP(_RotationMatrix):
    f = staticmethod(t._py_rotation_matrix)

class TestRotationMatrixCy(_RotationMatrix):
    f = staticmethod(t.rotation_matrix)


def test_rotation_from_matrix():
    angle = (random.random() - 0.5) * (2*np.pi)
    direc = np.random.random(3) - 0.5
    point = np.random.random(3) - 0.5
    R0 = t.rotation_matrix(angle, direc, point)
    angle, direc, point = t.rotation_from_matrix(R0)
    R1 = t.rotation_matrix(angle, direc, point)
    assert_equal(t.is_same_transform(R0, R1), True)

class _ScaleMatrix(object):
    def test_scale_matrix(self):
        v = (np.random.rand(4, 5) - 0.5) * 20.0
        v[3] = 1.0
        S = self.f(-1.234)
        assert_allclose(np.dot(S, v)[:3], -1.234*v[:3])

class TestScaleMatrixNP(_ScaleMatrix):
    f = staticmethod(t._py_scale_matrix)

class TestScaleMatrixCy(_ScaleMatrix):
    f = staticmethod(t.scale_matrix)

def test_scale_from_matrix():
    factor = random.random() * 10 - 5
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

class _ProjectionMatrix(object):
    def test_projection_matrix_1(self):
        P = self.f((0, 0, 0), (1, 0, 0))
        assert_allclose(P[1:, 1:], np.identity(4)[1:, 1:], atol=_ATOL)

    def test_projection_matrix_2(self):
        point = np.random.random(3) - 0.5
        normal = np.random.random(3) - 0.5
        direct = np.random.random(3) - 0.5
        persp = np.random.random(3) - 0.5
        P0 = self.f(point, normal)
        P1 = self.f(point, normal, direction=direct)
        P2 = self.f(point, normal, perspective=persp)
        P3 = self.f(point, normal, perspective=persp, pseudo=True)
        assert_equal(t.is_same_transform(P2, np.dot(P0, P3)), True)

    def test_projection_matrix_3(self):
        P = self.f((3, 0, 0), (1, 1, 0), (1, 0, 0))
        v0 = (np.random.rand(4, 5) - 0.5) * 20.0
        v0[3] = 1.0
        v1 = np.dot(P, v0)
        assert_allclose(v1[1], v0[1], atol=_ATOL)
        assert_allclose(v1[0], 3.0-v1[1], atol=_ATOL)

class TestProjectionMatrixNP(_ProjectionMatrix):
    f = staticmethod(t._py_projection_matrix)

class TestProjectionMatrixCy(_ProjectionMatrix):
    f = staticmethod(t.projection_matrix)


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


class _ClipMatrix(object):
    def test_clip_matrix_1(self):
        frustrum = np.random.rand(6)
        frustrum[1] += frustrum[0]
        frustrum[3] += frustrum[2]
        frustrum[5] += frustrum[4]
        M = self.f(perspective=False, *frustrum)
        assert_allclose(np.dot(M, [frustrum[0], frustrum[2], frustrum[4], 1.0]),
                        np.array([-1., -1., -1.,  1.]))
        assert_allclose(np.dot(M, [frustrum[1], frustrum[3], frustrum[5], 1.0]),
                         np.array([ 1.,  1.,  1.,  1.]))

    def test_clip_matrix_2(self):
        frustrum = np.random.rand(6)
        frustrum[1] += frustrum[0]
        frustrum[3] += frustrum[2]
        frustrum[5] += frustrum[4]
        M = self.f(perspective=True, *frustrum)
        v = np.dot(M, [frustrum[0], frustrum[2], frustrum[4], 1.0])
        assert_allclose(v / v[3],
                        np.array([-1., -1., -1.,  1.]))
        v = np.dot(M, [frustrum[1], frustrum[3], frustrum[4], 1.0])
        assert_allclose(v / v[3],
                        np.array([ 1.,  1., -1.,  1.]))

class TestClipMatrixNP(_ClipMatrix):
    f = staticmethod(t._py_clip_matrix)

class TestClipMatrixCy(_ClipMatrix):
    f = staticmethod(t.clip_matrix)


class _ShearMatrix(object):
    def test_shear_matrix(self):
        angle = (random.random() - 0.5) * 4*np.pi
        direct = np.random.random(3) - 0.5
        point = np.random.random(3) - 0.5
        normal = np.cross(direct, np.random.random(3))
        S = self.f(angle, direct, point, normal)
        assert_allclose(1.0, np.linalg.det(S), atol=_ATOL)

class TestShearMatrixNP(_ShearMatrix):
    f = staticmethod(t._py_shear_matrix)

class TestShearMatrixCy(_ShearMatrix):
    f = staticmethod(t.shear_matrix)


def test_shear_from_matrix():
    # This seems to fail sometimes if the random numbers
    # roll certain values....
    # angle = (random.random() - 0.5) * 4*math.pi
    # direct = np.random.random(3) - 0.5
    # point = np.random.random(3) - 0.5
    # normal = np.cross(direct, np.random.random(3))
    # In this random configuration the test will fail about 0.05% of all times.
    # Then we hit some edge-cases of the algorithm. The edge cases for these
    # values are slightly different for the linalg library used (MKL/LAPACK).
    # So here are some of my random numbers
    angle = 2.8969075413405783
    direct = np.array([-0.31117458, -0.41769518, -0.01188556])
    point = np.array([-0.0035982, -0.40997482,  0.42241425])
    normal = np.cross(direct, np.array([ 0.08122421,  0.4747914 ,  0.19851859]))

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

class _OrthogonalizationMatrix(object):
    def test_orthogonalization_matrix_1(self):
        O = self.f((10., 10., 10.), (90., 90., 90.))
        assert_allclose(O[:3, :3], np.identity(3, float) * 10, atol=_ATOL)

    def test_orthogonalization_matrix_2(self):
        O = self.f([9.8, 12.0, 15.5], [87.2, 80.7, 69.7])
        assert_allclose(np.sum(O), 43.063229, atol=_ATOL)

class TestOrthogonalizationMatrixNP(_OrthogonalizationMatrix):
    f = staticmethod(t._py_orthogonalization_matrix)

class TestOrthogonalizationMatrixCy(_OrthogonalizationMatrix):
    f = staticmethod(t.orthogonalization_matrix)


class _SuperimpositionMatrix(object):
    def test_superimposition_matrix(self):
        v0 = np.random.rand(3, 10)
        M = self.f(v0, v0)
        assert_allclose(M, np.identity(4), atol=_ATOL)

        R = t.random_rotation_matrix(np.random.random(3))
        v0 = ((1,0,0), (0,1,0), (0,0,1), (1,1,1))
        v1 = np.dot(R, v0)
        M = self.f(v0, v1)
        assert_allclose(v1, np.dot(M, v0), atol=_ATOL)

        v0 = (np.random.rand(4, 100) - 0.5) * 20.0
        v0[3] = 1.0
        v1 = np.dot(R, v0)
        M = self.f(v0, v1)
        assert_allclose(v1, np.dot(M, v0), atol=_ATOL)

        S = t.scale_matrix(random.random())
        T = t.translation_matrix(np.random.random(3)-0.5)
        M = t.concatenate_matrices(T, R, S)
        v1 = np.dot(M, v0)
        v0[:3] += np.random.normal(0.0, 1e-9, 300).reshape(3, -1)
        M = self.f(v0, v1, scaling=True)
        assert_allclose(v1, np.dot(M, v0), atol=_ATOL)

        M = self.f(v0, v1, scaling=True, usesvd=False)
        assert_allclose(v1, np.dot(M, v0), atol=_ATOL)

        v = np.empty((4, 100, 3), dtype=np.float64)
        v[:, :, 0] = v0
        M = self.f(v0, v1, scaling=True, usesvd=False)
        assert_allclose(v1, np.dot(M, v[:, :, 0]), atol=_ATOL)

class TestSuperimpositionMatrixNP(_SuperimpositionMatrix):
    f = staticmethod(t._py_superimposition_matrix)

class TestSuperimpositionMatrixCy(_SuperimpositionMatrix):
    f = staticmethod(t.superimposition_matrix)


class _EulerMatrix(object):
    def test_euler_matrix_1(self):
        R = self.f(1, 2, 3, 'syxz')
        assert_allclose(np.sum(R[0]), -1.34786452)

    def test_euler_matrix_2(self):
        R = self.f(1, 2, 3, (0, 1, 0, 1))
        assert_allclose(np.sum(R[0]), -0.383436184)

class TestEulerMatrixNP(_EulerMatrix):
    f = staticmethod(t._py_euler_matrix)

class TestEulerMatrixCy(_EulerMatrix):
    f = staticmethod(t.euler_matrix)


class _EulerFromMatrix(object):
    def test_euler_from_matrix_1(self):
        R0 = t.euler_matrix(1, 2, 3, 'syxz')
        al, be, ga = self.f(R0, 'syxz')
        R1 = t.euler_matrix(al, be, ga, 'syxz')
        assert_allclose(R0, R1)

    def test_euler_from_matrix_2(self):
        angles = (4.0*math.pi) * (np.random.random(3) - 0.5)
        for axes in t._AXES2TUPLE.keys():
            R0 = t.euler_matrix(axes=axes, *angles)
            R1 = t.euler_matrix(axes=axes, *self.f(R0, axes))
            assert_allclose(R0, R1, err_msg=("{0} failed".format(axes)))

class TestEulerFromMatrixNP(_EulerFromMatrix):
    f = staticmethod(t._py_euler_from_matrix)

class TestEulerFromMatrixCy(_EulerFromMatrix):
    f = staticmethod(t.euler_from_matrix)


def test_euler_from_quaternion():
    angles = t.euler_from_quaternion([0.99810947, 0.06146124, 0, 0])
    assert_allclose(angles, [0.123, 0, 0], atol=_ATOL)

class _QuaternionFromEuler(object):
    def test_quaternion_from_euler(self):
        q = self.f(1, 2, 3, 'ryxz')
        assert_allclose(q, [0.435953, 0.310622, -0.718287, 0.444435], atol=_ATOL)

class TestQuaternionFromEulerNP(_QuaternionFromEuler):
    f = staticmethod(t._py_quaternion_from_euler)

class TestQuaternionFromEulerCy(_QuaternionFromEuler):
    f = staticmethod(t.quaternion_from_euler)


class _QuaternionAboutAxis(object):
    def test_quaternion_about_axis(self):
        q = self.f(0.123, (1, 0, 0))
        assert_allclose(q, [0.99810947, 0.06146124, 0, 0], atol=_ATOL)

class TestQuaternionAboutAxisNP(_QuaternionAboutAxis):
    f = staticmethod(t._py_quaternion_about_axis)

class TestQuaternionAboutAxisCy(_QuaternionAboutAxis):
    f = staticmethod(t.quaternion_about_axis)


class _QuaternionMatrix(object):
    def test_quaternion_matrix_1(self):
        M = self.f([0.99810947, 0.06146124, 0, 0])
        assert_allclose(M, t.rotation_matrix(0.123, (1, 0, 0)), atol=_ATOL)

    def test_quaternion_matrix_2(self):
        M = self.f([1, 0, 0, 0])
        assert_allclose(M, t.identity_matrix(), atol=_ATOL)

    def test_quaternion_matrix_3(self):
        M = self.f([0, 1, 0, 0])
        assert_allclose(M, np.diag([1, -1, -1, 1]), atol=_ATOL)

class TestQuaternionMatrixNP(_QuaternionMatrix):
    f = staticmethod(t._py_quaternion_matrix)

class TestQuaternionMatrixCy(_QuaternionMatrix):
    f = staticmethod(t.quaternion_matrix)


class _QuaternionFromMatrix(object):
    def test_quaternion_from_matrix_1(self):
        q = self.f(t.identity_matrix(), True)
        assert_allclose(q, [1., 0., 0., 0.], atol=_ATOL)

    def test_quaternion_from_matrix_2(self):
        q = self.f(np.diag([1., -1., -1., 1.]))
        check = (np.allclose(q, [0, 1, 0, 0], atol=_ATOL) or
                 np.allclose(q, [0, -1, 0, 0], atol=_ATOL))
        assert_equal(check, True)

    def test_quaternion_from_matrix_3(self):
        R = t.rotation_matrix(0.123, (1, 2, 3))
        q = self.f(R, True)
        assert_allclose(q, [0.9981095, 0.0164262, 0.0328524, 0.0492786], atol=_ATOL)

    def test_quaternion_from_matrix_4(self):
        R = [[-0.545, 0.797, 0.260, 0], [0.733, 0.603, -0.313, 0],
             [-0.407, 0.021, -0.913, 0], [0, 0, 0, 1]]
        q = self.f(R)
        assert_allclose(q, [0.19069, 0.43736, 0.87485, -0.083611], atol=_ATOL)

    def test_quaternion_from_matrix_5(self):
        R = [[0.395, 0.362, 0.843, 0], [-0.626, 0.796, -0.056, 0],
             [-0.677, -0.498, 0.529, 0], [0, 0, 0, 1]]
        q = self.f(R)
        assert_allclose(q, [0.82336615, -0.13610694, 0.46344705, -0.29792603],
                        atol=_ATOL)

    def test_quaternion_from_matrix_6(self):
        R = t.random_rotation_matrix()
        q = self.f(R)
        assert_equal(t.is_same_transform(R, t.quaternion_matrix(q)), True)

class TestQuaternionFromMatrixNP(_QuaternionFromMatrix):
    f = staticmethod(t._py_quaternion_from_matrix)

class TestQuaternionFromMatrixCy(_QuaternionFromMatrix):
    f = staticmethod(t.quaternion_from_matrix)


class _QuaternionMultiply(object):
    def test_quaternion_multiply(self):
        q = self.f([4, 1, -2, 3], [8, -5, 6, 7])
        assert_allclose(q, [28, -44, -14, 48])

class TestQuaternionMultiplyNP(_QuaternionMultiply):
    f = staticmethod(t._py_quaternion_multiply)

class TestQuaternionMultiplyCy(_QuaternionMultiply):
    f = staticmethod(t.quaternion_multiply)


class _QuaternionConjugate(object):
    def test_quaternion_conjugate(self):
        q0 = t.random_quaternion()
        q1 = self.f(q0)
        check = q1[0] == q0[0] and all(q1[1:] == -q0[1:])
        assert_equal(check, True)

class TestQuaternionConjugateNP(_QuaternionConjugate):
    f = staticmethod(t._py_quaternion_conjugate)

class TestQuaternionConjugateCy(_QuaternionConjugate):
    f = staticmethod(t.quaternion_conjugate)


class _QuaternionInverse(object):
    def test_quaternion_inverse(self):
        q0 = t.random_quaternion()
        q1 = self.f(q0)
        assert_allclose(t.quaternion_multiply(q0, q1), [1, 0, 0, 0], atol=_ATOL)

class TestQuaternionInverseNP(_QuaternionInverse):
    f = staticmethod(t._py_quaternion_inverse)

class TestQuaternionInverseCy(_QuaternionInverse):
    f = staticmethod(t.quaternion_inverse)

def test_quaternion_real():
    assert_allclose(t.quaternion_real([3.0, 0.0, 1.0, 2.0]),
                    3.0)

def test_quaternion_imag():
    assert_allclose(t.quaternion_imag([3.0, 0.0, 1.0, 2.0]),
                    [0.0, 1.0, 2.0])

class _QuaternionSlerp(object):
    def test_quaternion_slerp(self):
        q0 = t.random_quaternion()
        q1 = t.random_quaternion()
        q = self.f(q0, q1, 0.0)
        assert_allclose(q, q0, atol=_ATOL)

        q = self.f(q0, q1, 1.0, 1)
        assert_allclose(q, q1, atol=_ATOL)

        q = self.f(q0, q1, 0.5)
        angle = math.acos(np.dot(q0, q))

        check = (np.allclose(2.0, math.acos(np.dot(q0, q1)) / angle) or
                 np.allclose(2.0, math.acos(-np.dot(q0, q1)) / angle))

        assert_equal(check, True)

class TestQuaternionSlerpNP(_QuaternionSlerp):
    f = staticmethod(t._py_quaternion_slerp)

class TestQuaternionSlerpCy(_QuaternionSlerp):
    f = staticmethod(t.quaternion_slerp)


class _RandomQuaternion(object):
    def test_random_quaternion_1(self):
        q = self.f()
        assert_allclose(1.0, t.vector_norm(q))

    def test_random_quaternion_2(self):
        q = self.f(np.random.random(3))
        assert_equal(len(q.shape), 1)
        assert_equal(q.shape[0] == 4, True)

class TestRandomQuaternionNP(_RandomQuaternion):
    f = staticmethod(t._py_random_quaternion)

class TestRandomQuaternionCy(_RandomQuaternion):
    f = staticmethod(t.random_quaternion)


class _RandomRotationMatrix(object):
    def test_random_rotation_matrix(self):
        R = self.f()
        assert_allclose(np.dot(R.T, R), np.identity(4), atol=_ATOL)

class TestRandomRotationMatrixNP(_RandomRotationMatrix):
    f = staticmethod(t._py_random_rotation_matrix)

class TestRandomRotationMatrixCy(_RandomRotationMatrix):
    f = staticmethod(t.random_rotation_matrix)


class _InverseMatrix(object):
    def test_inverse_matrix(self):
        M0 = t.random_rotation_matrix()
        M1 = self.f(M0.T)
        assert_allclose(M1, np.linalg.inv(M0.T))

        for size in range(1, 7):
            M0 = np.random.rand(size, size)
            M1 = self.f(M0)
            assert_allclose(M1, np.linalg.inv(M0), err_msg=str(size))

class TestInverseMatrixNP(_InverseMatrix):
    f = staticmethod(t._py_inverse_matrix)

class TestInverseMatrixCy(_InverseMatrix):
    f = staticmethod(t.inverse_matrix)


class _IsSameTransform(object):
    def test_is_same_transform_1(self):
        assert_equal(self.f(np.identity(4), np.identity(4)), True)

    def test_is_same_transform_2(self):
        assert_equal(self.f(t.random_rotation_matrix(), np.identity(4)), False)

class TestIsSameTransformNP(_IsSameTransform):
    f = staticmethod(t._py_is_same_transform)

class TestIsSameTransformCy(_IsSameTransform):
    f = staticmethod(t.is_same_transform)

class _RandomVector(object):
    def test_random_vector_1(self):
        v = self.f(1000)
        check = np.all(v >= 0.0) and np.all(v < 1.0)
        assert_equal(check, True)

    def test_random_vector_2(self):
        v0 = self.f(10)
        v1 = self.f(10)
        assert_equal(np.any(v0 == v1), False)

class TestRandomVectorNP(_RandomVector):
     f = staticmethod(t._py_random_vector)

class TestRandomVectorCy(_RandomVector):
     f = staticmethod(t.random_vector)


class _UnitVector(object):
    def test_unit_vector_1(self):
        v0 = np.random.random(3)
        v1 = self.f(v0)
        assert_allclose(v1, v0 / np.linalg.norm(v0), atol=_ATOL)

    def test_unit_vector_2(self):
        v0 = np.random.rand(5, 4, 3)
        v1 = self.f(v0, axis=-1)
        v2 = v0 / np.expand_dims(np.sqrt(np.sum(v0*v0, axis=2)), 2)
        assert_allclose(v1, v2, atol=_ATOL)

    def test_unit_vector_3(self):
        v0 = np.random.rand(5, 4, 3)
        v1 = self.f(v0, axis=1)
        v2 = v0 / np.expand_dims(np.sqrt(np.sum(v0*v0, axis=1)), 1)
        assert_allclose(v1, v2, atol=_ATOL)

    def test_unit_vector_4(self):
        v0 = np.random.rand(5, 4, 3)
        v1 = np.empty((5, 4, 3), dtype=np.float64)
        v2 = v0 / np.expand_dims(np.sqrt(np.sum(v0*v0, axis=1)), 1)
        self.f(v0, axis=1, out=v1)
        assert_allclose(v1, v2, atol=_ATOL)

    def test_unit_vector_5(self):
        assert_equal(list(self.f([])), [])

    def test_unit_vector_6(self):
        assert_equal(list(self.f([1.0])), [1.0])

class TestUnitVectorNP(_UnitVector):
    f = staticmethod(t._py_unit_vector)

class TestUnitVectorCy(_UnitVector):
    f = staticmethod(t.unit_vector)


class _VectorNorm(object):
    def test_vector_norm_1(self):
        v = np.random.random(3)
        n = self.f(v)
        assert_allclose(n, np.linalg.norm(v), atol=_ATOL)

    def test_vector_norm_2(self):
        v = np.random.rand(6, 5, 3)
        n = self.f(v, axis=-1)
        assert_allclose(n, np.sqrt(np.sum(v*v, axis=2)), atol=_ATOL)

    def test_vector_norm_3(self):
        v = np.random.rand(6, 5, 3)
        n = self.f(v, axis=1)
        assert_allclose(n, np.sqrt(np.sum(v*v, axis=1)), atol=_ATOL)

    def test_vector_norm_4(self):
        v = np.random.rand(5, 4, 3)
        n = np.empty((5, 3), dtype=np.float64)
        self.f(v, axis=1, out=n)
        assert_allclose(n, np.sqrt(np.sum(v*v, axis=1)), atol=_ATOL)

    def test_vector_norm_5(self):
        assert_equal(self.f([]), 0.0)

    def test_vector_norm_6(self):
        assert_equal(self.f([1.0]), 1.0)

class TestVectorNormNP(_VectorNorm):
    f = staticmethod(t._py_vector_norm)

class TestVectorNormCy(_VectorNorm):
    f = staticmethod(t.vector_norm)

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
        ball.setaxes([1,1,0], [-1, 1, 0])
        ball.setconstrain(True)
        ball.down([400, 200])
        ball.drag([200, 400])
        R = ball.matrix()
        assert_allclose(np.sum(R), 0.2055924)


def test_transformations_old_module():
    """test that MDAnalysis.core.transformations is still importable (deprecated for 1.0)"""
    try:
        import MDAnalysis.core.transformations
    except (ImportError, NameError):
        raise AssertionError("MDAnalysis.core.transformations not importable. Only remove for 1.0")

    # NOTE: removed this test with release 1.0 when we remove the stub
