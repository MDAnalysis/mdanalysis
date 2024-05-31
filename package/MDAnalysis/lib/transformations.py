# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

# transformations.py

# Copyright (c) 2006, Christoph Gohlke
# Copyright (c) 2006-2010, The Regents of the University of California
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# * Neither the name of the copyright holders nor the names of any
#   contributors may be used to endorse or promote products derived
#   from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

"""
Homogeneous Transformation Matrices and Quaternions --- :mod:`MDAnalysis.lib.transformations`
==============================================================================================

A library for calculating 4x4 matrices for translating, rotating, reflecting,
scaling, shearing, projecting, orthogonalizing, and superimposing arrays of
3D homogeneous coordinates as well as for converting between rotation matrices,
Euler angles, and quaternions. Also includes an Arcball control object and
functions to decompose transformation matrices.

:Authors:
  `Christoph Gohlke <http://www.lfd.uci.edu/~gohlke/>`__,
  Laboratory for Fluorescence Dynamics, University of California, Irvine
:Version: 2010.05.10
:Licence: BSD 3-clause

Requirements
------------

* `Python 2.6 or 3.1 <http://www.python.org>`__
* `Numpy 1.4 <http://numpy.scipy.org>`__
* `transformations.c 2010.04.10 <http://www.lfd.uci.edu/~gohlke/>`__
  (optional implementation of some functions in C)

Notes
-----

The API is not stable yet and is expected to change between revisions.

This Python code is not optimized for speed. Refer to the transformations.c
module for a faster implementation of some functions.

Documentation in HTML format can be generated with epydoc.

Matrices (M) can be inverted using ``numpy.linalg.inv(M)``, concatenated using
``numpy.dot(M0, M1)``, or used to transform homogeneous coordinates (v) using
``numpy.dot(M, v)`` for shape ``(4, *)`` "point of arrays", respectively
``numpy.dot(v, M.T)`` for shape ``(*, 4)`` "array of points".

Use the transpose of transformation matrices for OpenGL ``glMultMatrixd()``.

Calculations are carried out with ``numpy.float64`` precision.

Vector, point, quaternion, and matrix function arguments are expected to be
"array like", i.e. tuple, list, or numpy arrays.

Return types are numpy arrays unless specified otherwise.

Angles are in radians unless specified otherwise.

Quaternions w+ix+jy+kz are represented as ``[w, x, y, z]``.

A triple of Euler angles can be applied/interpreted in 24 ways, which can
be specified using a 4 character string or encoded 4-tuple:

  - *Axes 4-string*: e.g. 'sxyz' or 'ryxy'

    - first character : rotations are applied to 's'tatic or 'r'otating frame
    - remaining characters : successive rotation axis 'x', 'y', or 'z'

  - *Axes 4-tuple*: e.g. (0, 0, 0, 0) or (1, 1, 1, 1)

    - inner axis: code of axis ('x':0, 'y':1, 'z':2) of rightmost matrix.
    - parity : even (0) if inner axis 'x' is followed by 'y', 'y' is followed
      by 'z', or 'z' is followed by 'x'. Otherwise odd (1).
    - repetition : first and last axis are same (1) or different (0).
    - frame : rotations are applied to static (0) or rotating (1) frame.

.. rubric:: References

.. footbibliography::

Examples
--------

>>> from MDAnalysis.lib.transformations import *
>>> import numpy as np
>>> alpha, beta, gamma = 0.123, -1.234, 2.345
>>> origin, xaxis, yaxis, zaxis = (0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)
>>> I = identity_matrix()
>>> Rx = rotation_matrix(alpha, xaxis)
>>> Ry = rotation_matrix(beta, yaxis)
>>> Rz = rotation_matrix(gamma, zaxis)
>>> R = concatenate_matrices(Rx, Ry, Rz)
>>> euler = euler_from_matrix(R, 'rxyz')
>>> np.allclose([alpha, beta, gamma], euler)
True
>>> Re = euler_matrix(alpha, beta, gamma, 'rxyz')
>>> is_same_transform(R, Re)
True
>>> al, be, ga = euler_from_matrix(Re, 'rxyz')
>>> is_same_transform(Re, euler_matrix(al, be, ga, 'rxyz'))
True
>>> qx = quaternion_about_axis(alpha, xaxis)
>>> qy = quaternion_about_axis(beta, yaxis)
>>> qz = quaternion_about_axis(gamma, zaxis)
>>> q = quaternion_multiply(qx, qy)
>>> q = quaternion_multiply(q, qz)
>>> Rq = quaternion_matrix(q)
>>> is_same_transform(R, Rq)
True
>>> S = scale_matrix(1.23, origin)
>>> T = translation_matrix((1, 2, 3))
>>> Z = shear_matrix(beta, xaxis, origin, zaxis)
>>> R = random_rotation_matrix(np.random.rand(3))
>>> M = concatenate_matrices(T, R, Z, S)
>>> scale, shear, angles, trans, persp = decompose_matrix(M)
>>> np.allclose(scale, 1.23)
True
>>> np.allclose(trans, (1, 2, 3))
True
>>> np.allclose(shear, (0, math.tan(beta), 0))
True
>>> is_same_transform(R, euler_matrix(axes='sxyz', *angles))
True
>>> M1 = compose_matrix(scale, shear, angles, trans, persp)
>>> is_same_transform(M, M1)
True

Functions
---------

.. See `help(MDAnalysis.lib.transformations)` for a listing of functions or
.. the online help.


.. versionchanged:: 0.11.0
   Transformations library moved from MDAnalysis.core.transformations to
   MDAnalysis.lib.transformations
"""

import sys
import os
import warnings
import math
import numpy as np

def identity_matrix():
    """Return 4x4 identity/unit matrix.

    >>> from MDAnalysis.lib.transformations import identity_matrix
    >>> import numpy as np
    >>> I = identity_matrix()
    >>> np.allclose(I, np.dot(I, I))
    True
    >>> np.sum(I), np.trace(I)
    (4.0, 4.0)
    >>> np.allclose(I, np.identity(4, dtype=np.float64))
    True

    """
    return np.identity(4, dtype=np.float64)


def translation_matrix(direction):
    """Return matrix to translate by direction vector.

    >>> from MDAnalysis.lib.transformations import translation_matrix
    >>> import numpy as np
    >>> v = np.random.random(3) - 0.5
    >>> np.allclose(v, translation_matrix(v)[:3, 3])
    True

    """
    M = np.identity(4)
    M[:3, 3] = direction[:3]
    return M


def translation_from_matrix(matrix):
    """Return translation vector from translation matrix.

    >>> from MDAnalysis.lib.transformations import (translation_matrix,
    ...     translation_from_matrix)
    >>> import numpy as np
    >>> v0 = np.random.random(3) - 0.5
    >>> v1 = translation_from_matrix(translation_matrix(v0))
    >>> np.allclose(v0, v1)
    True

    """
    return np.array(matrix, copy=False)[:3, 3].copy()


def reflection_matrix(point, normal):
    """Return matrix to mirror at plane defined by point and normal vector.

    >>> from MDAnalysis.lib.transformations import reflection_matrix
    >>> import numpy as np
    >>> v0 = np.random.random(4) - 0.5
    >>> v0[3] = 1.0
    >>> v1 = np.random.random(3) - 0.5
    >>> R = reflection_matrix(v0, v1)
    >>> np.allclose(2., np.trace(R))
    True
    >>> np.allclose(v0, np.dot(R, v0))
    True
    >>> v2 = v0.copy()
    >>> v2[:3] += v1
    >>> v3 = v0.copy()
    >>> v2[:3] -= v1
    >>> np.allclose(v2, np.dot(R, v3))
    True

    """
    normal = unit_vector(normal[:3])
    M = np.identity(4)
    M[:3, :3] -= 2.0 * np.outer(normal, normal)
    M[:3, 3] = (2.0 * np.dot(point[:3], normal)) * normal
    return M


def reflection_from_matrix(matrix):
    """Return mirror plane point and normal vector from reflection matrix.

    >>> from MDAnalysis.lib.transformations import (reflection_matrix,
    ...     reflection_from_matrix, is_same_transform)
    >>> import numpy as np
    >>> v0 = np.random.random(3) - 0.5
    >>> v1 = np.random.random(3) - 0.5
    >>> M0 = reflection_matrix(v0, v1)
    >>> point, normal = reflection_from_matrix(M0)
    >>> M1 = reflection_matrix(point, normal)
    >>> is_same_transform(M0, M1)
    True

    """
    M = np.array(matrix, dtype=np.float64, copy=False)
    # normal: unit eigenvector corresponding to eigenvalue -1
    l, V = np.linalg.eig(M[:3, :3])
    i = np.where(abs(np.real(l) + 1.0) < 1e-8)[0]
    if not len(i):
        raise ValueError("no unit eigenvector corresponding to eigenvalue -1")
    normal = np.real(V[:, i[0]]).squeeze()
    # point: any unit eigenvector corresponding to eigenvalue 1
    l, V = np.linalg.eig(M)
    i = np.where(abs(np.real(l) - 1.0) < 1e-8)[0]
    if not len(i):
        raise ValueError("no unit eigenvector corresponding to eigenvalue 1")
    point = np.real(V[:, i[-1]]).squeeze()
    point /= point[3]
    return point, normal


def rotation_matrix(angle, direction, point=None):
    """Return matrix to rotate about axis defined by point and direction.

    >>> from MDAnalysis.lib.transformations import (rotation_matrix,
    ...     is_same_transform)
    >>> import random, math
    >>> import numpy as np
    >>> R = rotation_matrix(math.pi/2.0, [0, 0, 1], [1, 0, 0])
    >>> np.allclose(np.dot(R, [0, 0, 0, 1]), [ 1., -1.,  0.,  1.])
    True
    >>> angle = (random.random() - 0.5) * (2*math.pi)
    >>> direc = np.random.random(3) - 0.5
    >>> point = np.random.random(3) - 0.5
    >>> R0 = rotation_matrix(angle, direc, point)
    >>> R1 = rotation_matrix(angle-2*math.pi, direc, point)
    >>> is_same_transform(R0, R1)
    True
    >>> R0 = rotation_matrix(angle, direc, point)
    >>> R1 = rotation_matrix(-angle, -direc, point)
    >>> is_same_transform(R0, R1)
    True
    >>> I = np.identity(4, np.float64)
    >>> np.allclose(I, rotation_matrix(math.pi*2, direc))
    True
    >>> np.allclose(2., np.trace(rotation_matrix(math.pi/2,
    ...                                                direc, point)))
    True

    """
    sina = math.sin(angle)
    cosa = math.cos(angle)
    direction = unit_vector(direction[:3])
    # rotation matrix around unit vector
    R = np.array(
        (
            (cosa, 0.0, 0.0),
            (0.0, cosa, 0.0),
            (0.0, 0.0, cosa)), dtype=np.float64)
    R += np.outer(direction, direction) * (1.0 - cosa)
    direction *= sina
    R += np.array(
        ((0.0, -direction[2], direction[1]),
        (direction[2], 0.0, -direction[0]),
        (-direction[1], direction[0], 0.0)), dtype=np.float64)
    M = np.identity(4)
    M[:3, :3] = R
    if point is not None:
        # rotation not around origin
        point = np.array(point[:3], dtype=np.float64, copy=no_copy_shim)
        M[:3, 3] = point - np.dot(R, point)
    return M


def rotation_from_matrix(matrix):
    """Return rotation angle and axis from rotation matrix.

    >>> from MDAnalysis.lib.transformations import (rotation_matrix,
    ...     is_same_transform, rotation_from_matrix)
    >>> import random, math
    >>> import numpy as np
    >>> angle = (random.random() - 0.5) * (2*math.pi)
    >>> direc = np.random.random(3) - 0.5
    >>> point = np.random.random(3) - 0.5
    >>> R0 = rotation_matrix(angle, direc, point)
    >>> angle, direc, point = rotation_from_matrix(R0)
    >>> R1 = rotation_matrix(angle, direc, point)
    >>> is_same_transform(R0, R1)
    True

    """
    R = np.array(matrix, dtype=np.float64, copy=False)
    R33 = R[:3, :3]
    # direction: unit eigenvector of R33 corresponding to eigenvalue of 1
    l, W = np.linalg.eig(R33.T)
    i = np.where(abs(np.real(l) - 1.0) < 1e-8)[0]
    if not len(i):
        raise ValueError("no unit eigenvector corresponding to eigenvalue 1")
    direction = np.real(W[:, i[-1]]).squeeze()
    # point: unit eigenvector of R33 corresponding to eigenvalue of 1
    l, Q = np.linalg.eig(R)
    i = np.where(abs(np.real(l) - 1.0) < 1e-8)[0]
    if not len(i):
        raise ValueError("no unit eigenvector corresponding to eigenvalue 1")
    point = np.real(Q[:, i[-1]]).squeeze()
    point /= point[3]
    # rotation angle depending on direction
    cosa = (np.trace(R33) - 1.0) / 2.0
    if abs(direction[2]) > 1e-8:
        sina = (R[1, 0] + (cosa - 1.0) * direction[0] * direction[1]) / direction[2]
    elif abs(direction[1]) > 1e-8:
        sina = (R[0, 2] + (cosa - 1.0) * direction[0] * direction[2]) / direction[1]
    else:
        sina = (R[2, 1] + (cosa - 1.0) * direction[1] * direction[2]) / direction[0]
    angle = math.atan2(sina, cosa)
    return angle, direction, point


def scale_matrix(factor, origin=None, direction=None):
    """Return matrix to scale by factor around origin in direction.

    Use factor -1 for point symmetry.

    >>> from MDAnalysis.lib.transformations import scale_matrix
    >>> import random
    >>> import numpy as np
    >>> v = (np.random.rand(4, 5) - 0.5) * 20.0
    >>> v[3] = 1.0
    >>> S = scale_matrix(-1.234)
    >>> np.allclose(np.dot(S, v)[:3], -1.234*v[:3])
    True
    >>> factor = random.random() * 10 - 5
    >>> origin = np.random.random(3) - 0.5
    >>> direct = np.random.random(3) - 0.5
    >>> S = scale_matrix(factor, origin)
    >>> S = scale_matrix(factor, origin, direct)

    """
    if direction is None:
        # uniform scaling
        M = np.array(
            ((factor, 0.0, 0.0, 0.0),
            (0.0, factor, 0.0, 0.0),
            (0.0, 0.0, factor, 0.0),
            (0.0, 0.0, 0.0, 1.0)), dtype=np.float64)
        if origin is not None:
            M[:3, 3] = origin[:3]
            M[:3, 3] *= 1.0 - factor
    else:
        # nonuniform scaling
        direction = unit_vector(direction[:3])
        factor = 1.0 - factor
        M = np.identity(4)
        M[:3, :3] -= factor * np.outer(direction, direction)
        if origin is not None:
            M[:3, 3] = (factor * np.dot(origin[:3], direction)) * direction
    return M


def scale_from_matrix(matrix):
    """Return scaling factor, origin and direction from scaling matrix.

    >>> from MDAnalysis.lib.transformations import (scale_matrix,
    ...      scale_from_matrix, is_same_transform)
    >>> import random
    >>> import numpy as np
    >>> factor = random.random() * 10 - 5
    >>> origin = np.random.random(3) - 0.5
    >>> direct = np.random.random(3) - 0.5
    >>> S0 = scale_matrix(factor, origin)
    >>> factor, origin, direction = scale_from_matrix(S0)
    >>> S1 = scale_matrix(factor, origin, direction)
    >>> is_same_transform(S0, S1)
    True
    >>> S0 = scale_matrix(factor, origin, direct)
    >>> factor, origin, direction = scale_from_matrix(S0)
    >>> S1 = scale_matrix(factor, origin, direction)
    >>> is_same_transform(S0, S1)
    True

    """
    M = np.array(matrix, dtype=np.float64, copy=False)
    M33 = M[:3, :3]
    factor = np.trace(M33) - 2.0
    try:
        # direction: unit eigenvector corresponding to eigenvalue factor
        l, V = np.linalg.eig(M33)
        i = np.where(abs(np.real(l) - factor) < 1e-8)[0][0]
        direction = np.real(V[:, i]).squeeze()
        direction /= vector_norm(direction)
    except IndexError:
        # uniform scaling
        factor = (factor + 2.0) / 3.0
        direction = None
    # origin: any eigenvector corresponding to eigenvalue 1
    l, V = np.linalg.eig(M)
    i = np.where(abs(np.real(l) - 1.0) < 1e-8)[0]
    if not len(i):
        raise ValueError("no eigenvector corresponding to eigenvalue 1")
    origin = np.real(V[:, i[-1]]).squeeze()
    origin /= origin[3]
    return factor, origin, direction


def projection_matrix(point, normal, direction=None,
                      perspective=None, pseudo=False):
    """Return matrix to project onto plane defined by point and normal.

    Using either perspective point, projection direction, or none of both.

    If pseudo is True, perspective projections will preserve relative depth
    such that Perspective = dot(Orthogonal, PseudoPerspective).

    >>> from MDAnalysis.lib.transformations import (projection_matrix,
    ...     is_same_transform)
    >>> import numpy as np
    >>> P = projection_matrix((0, 0, 0), (1, 0, 0))
    >>> np.allclose(P[1:, 1:], np.identity(4)[1:, 1:])
    True
    >>> point = np.random.random(3) - 0.5
    >>> normal = np.random.random(3) - 0.5
    >>> direct = np.random.random(3) - 0.5
    >>> persp = np.random.random(3) - 0.5
    >>> P0 = projection_matrix(point, normal)
    >>> P1 = projection_matrix(point, normal, direction=direct)
    >>> P2 = projection_matrix(point, normal, perspective=persp)
    >>> P3 = projection_matrix(point, normal, perspective=persp, pseudo=True)
    >>> is_same_transform(P2, np.dot(P0, P3))
    True
    >>> P = projection_matrix((3, 0, 0), (1, 1, 0), (1, 0, 0))
    >>> v0 = (np.random.rand(4, 5) - 0.5) * 20.0
    >>> v0[3] = 1.0
    >>> v1 = np.dot(P, v0)
    >>> np.allclose(v1[1], v0[1])
    True
    >>> np.allclose(v1[0], 3.0-v1[1])
    True

    """
    M = np.identity(4)
    point = np.array(point[:3], dtype=np.float64, copy=no_copy_shim)
    normal = unit_vector(normal[:3])
    if perspective is not None:
        # perspective projection
        perspective = np.array(perspective[:3], dtype=np.float64,
                                  copy=False)
        M[0, 0] = M[1, 1] = M[2, 2] = np.dot(perspective - point, normal)
        M[:3, :3] -= np.outer(perspective, normal)
        if pseudo:
            # preserve relative depth
            M[:3, :3] -= np.outer(normal, normal)
            M[:3, 3] = np.dot(point, normal) * (perspective + normal)
        else:
            M[:3, 3] = np.dot(point, normal) * perspective
        M[3, :3] = -normal
        M[3, 3] = np.dot(perspective, normal)
    elif direction is not None:
        # parallel projection
        direction = np.array(direction[:3], dtype=np.float64, copy=no_copy_shim)
        scale = np.dot(direction, normal)
        M[:3, :3] -= np.outer(direction, normal) / scale
        M[:3, 3] = direction * (np.dot(point, normal) / scale)
    else:
        # orthogonal projection
        M[:3, :3] -= np.outer(normal, normal)
        M[:3, 3] = np.dot(point, normal) * normal
    return M


def projection_from_matrix(matrix, pseudo=False):
    """Return projection plane and perspective point from projection matrix.

    Return values are same as arguments for projection_matrix function:
    point, normal, direction, perspective, and pseudo.

    >>> from MDAnalysis.lib.transformations import (projection_matrix,
    ...     projection_from_matrix, is_same_transform)
    >>> import numpy as np
    >>> point = np.random.random(3) - 0.5
    >>> normal = np.random.random(3) - 0.5
    >>> direct = np.random.random(3) - 0.5
    >>> persp = np.random.random(3) - 0.5
    >>> P0 = projection_matrix(point, normal)
    >>> result = projection_from_matrix(P0)
    >>> P1 = projection_matrix(*result)
    >>> is_same_transform(P0, P1)
    True
    >>> P0 = projection_matrix(point, normal, direct)
    >>> result = projection_from_matrix(P0)
    >>> P1 = projection_matrix(*result)
    >>> is_same_transform(P0, P1)
    True
    >>> P0 = projection_matrix(point, normal, perspective=persp, pseudo=False)
    >>> result = projection_from_matrix(P0, pseudo=False)
    >>> P1 = projection_matrix(*result)
    >>> is_same_transform(P0, P1)
    True
    >>> P0 = projection_matrix(point, normal, perspective=persp, pseudo=True)
    >>> result = projection_from_matrix(P0, pseudo=True)
    >>> P1 = projection_matrix(*result)
    >>> is_same_transform(P0, P1)
    True

    """
    M = np.array(matrix, dtype=np.float64, copy=False)
    M33 = M[:3, :3]
    l, V = np.linalg.eig(M)
    i = np.where(abs(np.real(l) - 1.0) < 1e-8)[0]
    if not pseudo and len(i):
        # point: any eigenvector corresponding to eigenvalue 1
        point = np.real(V[:, i[-1]]).squeeze()
        point /= point[3]
        # direction: unit eigenvector corresponding to eigenvalue 0
        l, V = np.linalg.eig(M33)
        i = np.where(abs(np.real(l)) < 1e-8)[0]
        if not len(i):
            raise ValueError("no eigenvector corresponding to eigenvalue 0")
        direction = np.real(V[:, i[0]]).squeeze()
        direction /= vector_norm(direction)
        # normal: unit eigenvector of M33.T corresponding to eigenvalue 0
        l, V = np.linalg.eig(M33.T)
        i = np.where(abs(np.real(l)) < 1e-8)[0]
        if len(i):
            # parallel projection
            normal = np.real(V[:, i[0]]).squeeze()
            normal /= vector_norm(normal)
            return point, normal, direction, None, False
        else:
            # orthogonal projection, where normal equals direction vector
            return point, direction, None, None, False
    else:
        # perspective projection
        i = np.where(abs(np.real(l)) > 1e-8)[0]
        if not len(i):
            raise ValueError(
                "no eigenvector not corresponding to eigenvalue 0")
        point = np.real(V[:, i[-1]]).squeeze()
        point /= point[3]
        normal = - M[3, :3]
        perspective = M[:3, 3] / np.dot(point[:3], normal)
        if pseudo:
            perspective -= normal
        return point, normal, None, perspective, pseudo


def clip_matrix(left, right, bottom, top, near, far, perspective=False):
    """Return matrix to obtain normalized device coordinates from frustrum.

    The frustrum bounds are axis-aligned along x (left, right),
    y (bottom, top) and z (near, far).

    Normalized device coordinates are in range [-1, 1] if coordinates are
    inside the frustrum.

    If perspective is True the frustrum is a truncated pyramid with the
    perspective point at origin and direction along z axis, otherwise an
    orthographic canonical view volume (a box).

    Homogeneous coordinates transformed by the perspective clip matrix
    need to be dehomogenized (devided by w coordinate).

    >>> from MDAnalysis.lib.transformations import clip_matrix
    >>> import numpy as np
    >>> frustrum = np.random.rand(6)
    >>> frustrum[1] += frustrum[0]
    >>> frustrum[3] += frustrum[2]
    >>> frustrum[5] += frustrum[4]
    >>> M = clip_matrix(perspective=False, *frustrum)
    >>> np.dot(M, [frustrum[0], frustrum[2], frustrum[4], 1.0])
    array([-1., -1., -1.,  1.])
    >>> np.dot(M, [frustrum[1], frustrum[3], frustrum[5], 1.0])
    array([1., 1., 1., 1.])
    >>> M = clip_matrix(perspective=True, *frustrum)
    >>> v = np.dot(M, [frustrum[0], frustrum[2], frustrum[4], 1.0])
    >>> v / v[3]
    array([-1., -1., -1.,  1.])
    >>> v = np.dot(M, [frustrum[1], frustrum[3], frustrum[4], 1.0])
    >>> v / v[3]
    array([ 1.,  1., -1.,  1.])

    """
    if left >= right or bottom >= top or near >= far:
        raise ValueError("invalid frustrum")
    if perspective:
        if near <= _EPS:
            raise ValueError("invalid frustrum: near <= 0")
        t = 2.0 * near
        M = (
            (-t / (right - left), 0.0, (right + left) / (right - left), 0.0),
            (0.0, -t / (top - bottom), (top + bottom) / (top - bottom), 0.0),
            (0.0, 0.0, -(far + near) / (far - near), t * far / (far - near)),
            (0.0, 0.0, -1.0, 0.0))
    else:
        M = (
            (2.0 / (right - left), 0.0, 0.0, (right + left) / (left - right)),
            (0.0, 2.0 / (top - bottom), 0.0, (top + bottom) / (bottom - top)),
            (0.0, 0.0, 2.0 / (far - near), (far + near) / (near - far)),
            (0.0, 0.0, 0.0, 1.0))
    return np.array(M, dtype=np.float64)


def shear_matrix(angle, direction, point, normal):
    """Return matrix to shear by angle along direction vector on shear plane.

    The shear plane is defined by a point and normal vector. The direction
    vector must be orthogonal to the plane's normal vector.

    A point P is transformed by the shear matrix into P" such that
    the vector P-P" is parallel to the direction vector and its extent is
    given by the angle of P-P'-P", where P' is the orthogonal projection
    of P onto the shear plane.

    >>> from MDAnalysis.lib.transformations import shear_matrix
    >>> import random, math
    >>> import numpy as np
    >>> angle = (random.random() - 0.5) * 4*math.pi
    >>> direct = np.random.random(3) - 0.5
    >>> point = np.random.random(3) - 0.5
    >>> normal = np.cross(direct, np.random.random(3))
    >>> S = shear_matrix(angle, direct, point, normal)
    >>> np.allclose(1.0, np.linalg.det(S))
    True

    """
    normal = unit_vector(normal[:3])
    direction = unit_vector(direction[:3])
    if abs(np.dot(normal, direction)) > 1e-6:
        raise ValueError("direction and normal vectors are not orthogonal")
    angle = math.tan(angle)
    M = np.identity(4)
    M[:3, :3] += angle * np.outer(direction, normal)
    M[:3, 3] = -angle * np.dot(point[:3], normal) * direction
    return M


def shear_from_matrix(matrix):
    """Return shear angle, direction and plane from shear matrix.

    >>> from MDAnalysis.lib.transformations import (shear_matrix,
    ...     shear_from_matrix, is_same_transform)
    >>> import random, math
    >>> import numpy as np
    >>> angle = (random.random() - 0.5) * 4*math.pi
    >>> direct = np.random.random(3) - 0.5
    >>> point = np.random.random(3) - 0.5
    >>> normal = np.cross(direct, np.random.random(3))
    >>> S0 = shear_matrix(angle, direct, point, normal)
    >>> angle, direct, point, normal = shear_from_matrix(S0)
    >>> S1 = shear_matrix(angle, direct, point, normal)
    >>> is_same_transform(S0, S1)
    True

    """
    M = np.array(matrix, dtype=np.float64, copy=False)
    M33 = M[:3, :3]
    # normal: cross independent eigenvectors corresponding to the eigenvalue 1
    l, V = np.linalg.eig(M33)
    i = np.where(abs(np.real(l) - 1.0) < 1e-4)[0]
    if len(i) < 2:
        raise ValueError("no two linear independent eigenvectors found {0!s}".format(l))
    V = np.real(V[:, i]).squeeze().T
    lenorm = -1.0
    for i0, i1 in ((0, 1), (0, 2), (1, 2)):
        n = np.cross(V[i0], V[i1])
        l = vector_norm(n)
        if l > lenorm:
            lenorm = l
            normal = n
    normal /= lenorm
    # direction and angle
    direction = np.dot(M33 - np.identity(3), normal)
    angle = vector_norm(direction)
    direction /= angle
    angle = math.atan(angle)
    # point: eigenvector corresponding to eigenvalue 1
    l, V = np.linalg.eig(M)
    i = np.where(abs(np.real(l) - 1.0) < 1e-8)[0]
    if not len(i):
        raise ValueError("no eigenvector corresponding to eigenvalue 1")
    point = np.real(V[:, i[-1]]).squeeze()
    point /= point[3]
    return angle, direction, point, normal


def decompose_matrix(matrix):
    """Return sequence of transformations from transformation matrix.

    matrix : array_like
        Non-degenerative homogeneous transformation matrix

    Return tuple of:
        scale : vector of 3 scaling factors
        shear : list of shear factors for x-y, x-z, y-z axes
        angles : list of Euler angles about static x, y, z axes
        translate : translation vector along x, y, z axes
        perspective : perspective partition of matrix

    Raise ValueError if matrix is of wrong type or degenerative.

    >>> from MDAnalysis.lib.transformations import (translation_matrix,
    ...     decompose_matrix, scale_matrix, euler_matrix)
    >>> import numpy as np
    >>> T0 = translation_matrix((1, 2, 3))
    >>> scale, shear, angles, trans, persp = decompose_matrix(T0)
    >>> T1 = translation_matrix(trans)
    >>> np.allclose(T0, T1)
    True
    >>> S = scale_matrix(0.123)
    >>> scale, shear, angles, trans, persp = decompose_matrix(S)
    >>> scale[0]
    0.123
    >>> R0 = euler_matrix(1, 2, 3)
    >>> scale, shear, angles, trans, persp = decompose_matrix(R0)
    >>> R1 = euler_matrix(*angles)
    >>> np.allclose(R0, R1)
    True

    """
    # pylint: disable=unsubscriptable-object
    M = np.array(matrix, dtype=np.float64, copy=True).T
    if abs(M[3, 3]) < _EPS:
        raise ValueError("M[3, 3] is zero")
    M /= M[3, 3]
    P = M.copy()
    P[:, 3] = 0, 0, 0, 1
    if not np.linalg.det(P):
        raise ValueError("matrix is singular")

    scale = np.zeros((3, ), dtype=np.float64)
    shear = [0, 0, 0]
    angles = [0, 0, 0]

    if any(abs(M[:3, 3]) > _EPS):
        perspective = np.dot(M[:, 3], np.linalg.inv(P.T))
        M[:, 3] = 0, 0, 0, 1
    else:
        perspective = np.array((0, 0, 0, 1), dtype=np.float64)

    translate = M[3, :3].copy()
    M[3, :3] = 0

    row = M[:3, :3].copy()
    scale[0] = vector_norm(row[0])
    row[0] /= scale[0]
    shear[0] = np.dot(row[0], row[1])
    row[1] -= row[0] * shear[0]
    scale[1] = vector_norm(row[1])
    row[1] /= scale[1]
    shear[0] /= scale[1]
    shear[1] = np.dot(row[0], row[2])
    row[2] -= row[0] * shear[1]
    shear[2] = np.dot(row[1], row[2])
    row[2] -= row[1] * shear[2]
    scale[2] = vector_norm(row[2])
    row[2] /= scale[2]
    shear[1:] /= scale[2]

    if np.dot(row[0], np.cross(row[1], row[2])) < 0:
        scale *= -1
        row *= -1

    angles[1] = math.asin(-row[0, 2])
    if math.cos(angles[1]):
        angles[0] = math.atan2(row[1, 2], row[2, 2])
        angles[2] = math.atan2(row[0, 1], row[0, 0])
    else:
        #angles[0] = math.atan2(row[1, 0], row[1, 1])
        angles[0] = math.atan2(-row[2, 1], row[1, 1])
        angles[2] = 0.0

    return scale, shear, angles, translate, perspective


def compose_matrix(scale=None, shear=None, angles=None, translate=None,
                   perspective=None):
    """Return transformation matrix from sequence of transformations.

    This is the inverse of the decompose_matrix function.

    Sequence of transformations:
        scale : vector of 3 scaling factors
        shear : list of shear factors for x-y, x-z, y-z axes
        angles : list of Euler angles about static x, y, z axes
        translate : translation vector along x, y, z axes
        perspective : perspective partition of matrix

    >>> from MDAnalysis.lib.transformations import (compose_matrix,
    ...     decompose_matrix, is_same_transform)
    >>> import math
    >>> import numpy as np
    >>> scale = np.random.random(3) - 0.5
    >>> shear = np.random.random(3) - 0.5
    >>> angles = (np.random.random(3) - 0.5) * (2*math.pi)
    >>> trans = np.random.random(3) - 0.5
    >>> persp = np.random.random(4) - 0.5
    >>> M0 = compose_matrix(scale, shear, angles, trans, persp)
    >>> result = decompose_matrix(M0)
    >>> M1 = compose_matrix(*result)
    >>> is_same_transform(M0, M1)
    True

    """
    M = np.identity(4)
    if perspective is not None:
        P = np.identity(4)
        P[3, :] = perspective[:4]
        M = np.dot(M, P)
    if translate is not None:
        T = np.identity(4)
        T[:3, 3] = translate[:3]
        M = np.dot(M, T)
    if angles is not None:
        R = euler_matrix(angles[0], angles[1], angles[2], 'sxyz')
        M = np.dot(M, R)
    if shear is not None:
        Z = np.identity(4)
        Z[1, 2] = shear[2]
        Z[0, 2] = shear[1]
        Z[0, 1] = shear[0]
        M = np.dot(M, Z)
    if scale is not None:
        S = np.identity(4)
        S[0, 0] = scale[0]
        S[1, 1] = scale[1]
        S[2, 2] = scale[2]
        M = np.dot(M, S)
    M /= M[3, 3]
    return M


def orthogonalization_matrix(lengths, angles):
    """Return orthogonalization matrix for crystallographic cell coordinates.

    Angles are expected in degrees.

    The de-orthogonalization matrix is the inverse.

    >>> from MDAnalysis.lib.transformations import orthogonalization_matrix
    >>> import numpy as np
    >>> O = orthogonalization_matrix((10., 10., 10.), (90., 90., 90.))
    >>> np.allclose(O[:3, :3], np.identity(3, float) * 10)
    True
    >>> O = orthogonalization_matrix([9.8, 12.0, 15.5], [87.2, 80.7, 69.7])
    >>> np.allclose(np.sum(O), 43.063229)
    True

    """
    a, b, c = lengths
    angles = np.radians(angles)
    sina, sinb, _ = np.sin(angles)
    cosa, cosb, cosg = np.cos(angles)
    co = (cosa * cosb - cosg) / (sina * sinb)
    return np.array(
        (
            (a * sinb * math.sqrt(1.0 - co * co), 0.0, 0.0, 0.0),
            (-a * sinb * co, b * sina, 0.0, 0.0),
            (a * cosb, b * cosa, c, 0.0),
            (0.0, 0.0, 0.0, 1.0)),
        dtype=np.float64)


def superimposition_matrix(v0, v1, scaling=False, usesvd=True):
    """Return matrix to transform given vector set into second vector set.

    `v0` and `v1` are shape `(3, *)` or `(4, *)` arrays of at least 3 vectors.

    If `usesvd` is ``True``, the weighted sum of squared deviations (RMSD) is
    minimized according to the algorithm by W. Kabsch [8]. Otherwise the
    quaternion based algorithm by B. Horn [9] is used (slower when using
    this Python implementation).

    The returned matrix performs rotation, translation and uniform scaling
    (if specified).

    >>> from MDAnalysis.lib.transformations import (superimposition_matrix,
    ... random_rotation_matrix, scale_matrix, translation_matrix,
    ... concatenate_matrices)
    >>> import random
    >>> import numpy as np
    >>> v0 = np.random.rand(3, 10)
    >>> M = superimposition_matrix(v0, v0)
    >>> np.allclose(M, np.identity(4))
    True
    >>> R = random_rotation_matrix(np.random.random(3))
    >>> v0 = ((1,0,0), (0,1,0), (0,0,1), (1,1,1))
    >>> v1 = np.dot(R, v0)
    >>> M = superimposition_matrix(v0, v1)
    >>> np.allclose(v1, np.dot(M, v0))
    True
    >>> v0 = (np.random.rand(4, 100) - 0.5) * 20.0
    >>> v0[3] = 1.0
    >>> v1 = np.dot(R, v0)
    >>> M = superimposition_matrix(v0, v1)
    >>> np.allclose(v1, np.dot(M, v0))
    True
    >>> S = scale_matrix(random.random())
    >>> T = translation_matrix(np.random.random(3)-0.5)
    >>> M = concatenate_matrices(T, R, S)
    >>> v1 = np.dot(M, v0)
    >>> v0[:3] += np.random.normal(0.0, 1e-9, 300).reshape(3, -1)
    >>> M = superimposition_matrix(v0, v1, scaling=True)
    >>> np.allclose(v1, np.dot(M, v0))
    True
    >>> M = superimposition_matrix(v0, v1, scaling=True, usesvd=False)
    >>> np.allclose(v1, np.dot(M, v0))
    True
    >>> v = np.empty((4, 100, 3), dtype=np.float64)
    >>> v[:, :, 0] = v0
    >>> M = superimposition_matrix(v0, v1, scaling=True, usesvd=False)
    >>> np.allclose(v1, np.dot(M, v[:, :, 0]))
    True

    """
    v0 = np.array(v0, dtype=np.float64, copy=no_copy_shim)[:3]
    v1 = np.array(v1, dtype=np.float64, copy=no_copy_shim)[:3]

    if v0.shape != v1.shape or v0.shape[1] < 3:
        raise ValueError("vector sets are of wrong shape or type")

    # move centroids to origin
    t0 = np.mean(v0, axis=1)
    t1 = np.mean(v1, axis=1)
    v0 = v0 - t0.reshape(3, 1)
    v1 = v1 - t1.reshape(3, 1)

    if usesvd:
        # Singular Value Decomposition of covariance matrix
        u, s, vh = np.linalg.svd(np.dot(v1, v0.T))
        # rotation matrix from SVD orthonormal bases
        R = np.dot(u, vh)
        if np.linalg.det(R) < 0.0:
            # R does not constitute right handed system
            R -= np.outer(u[:, 2], vh[2, :] * 2.0)
            s[-1] *= -1.0
        # homogeneous transformation matrix
        M = np.identity(4)
        M[:3, :3] = R
    else:
        # compute symmetric matrix N
        xx, yy, zz = np.einsum('ij,ij->i', v0 , v1)
        xy, yz, zx = np.einsum('ij,ij->i', v0, np.roll(v1, -1, axis=0))
        xz, yx, zy = np.einsum('ij,ij->i', v0, np.roll(v1, -2, axis=0))
        N = (
            (xx + yy + zz, 0.0, 0.0, 0.0),
            (yz - zy, xx - yy - zz, 0.0, 0.0),
            (zx - xz, xy + yx, -xx + yy - zz, 0.0),
            (xy - yx, zx + xz, yz + zy, -xx - yy + zz))
        # quaternion: eigenvector corresponding to most positive eigenvalue
        l, V = np.linalg.eigh(N)
        q = V[:, np.argmax(l)]
        q /= vector_norm(q)  # unit quaternion
        # homogeneous transformation matrix
        M = quaternion_matrix(q)

    # scale: ratio of rms deviations from centroid
    if scaling:
        M[:3, :3] *= math.sqrt(np.einsum('ij,ij->',v1,v1) / 
                                np.einsum('ij,ij->',v0,v0))


    # translation
    M[:3, 3] = t1
    T = np.identity(4)
    T[:3, 3] = -t0
    M = np.dot(M, T)
    return M


def euler_matrix(ai, aj, ak, axes='sxyz'):
    """Return homogeneous rotation matrix from Euler angles and axis sequence.

    ai, aj, ak : Euler's roll, pitch and yaw angles
    axes : One of 24 axis sequences as string or encoded tuple

    >>> from MDAnalysis.lib.transformations import (euler_matrix,
    ... _AXES2TUPLE, _TUPLE2AXES)
    >>> import math
    >>> import numpy as np
    >>> R = euler_matrix(1, 2, 3, 'syxz')
    >>> np.allclose(np.sum(R[0]), -1.34786452)
    True
    >>> R = euler_matrix(1, 2, 3, (0, 1, 0, 1))
    >>> np.allclose(np.sum(R[0]), -0.383436184)
    True
    >>> ai, aj, ak = (4.0*math.pi) * (np.random.random(3) - 0.5)
    >>> for axes in _AXES2TUPLE.keys():
    ...    R = euler_matrix(ai, aj, ak, axes)
    >>> for axes in _TUPLE2AXES.keys():
    ...    R = euler_matrix(ai, aj, ak, axes)

    """
    try:
        firstaxis, parity, repetition, frame = _AXES2TUPLE[axes]
    except (AttributeError, KeyError):
        _ = _TUPLE2AXES[axes]
        firstaxis, parity, repetition, frame = axes

    i = firstaxis
    j = _NEXT_AXIS[i + parity]
    k = _NEXT_AXIS[i - parity + 1]

    if frame:
        ai, ak = ak, ai
    if parity:
        ai, aj, ak = -ai, -aj, -ak

    si, sj, sk = math.sin(ai), math.sin(aj), math.sin(ak)
    ci, cj, ck = math.cos(ai), math.cos(aj), math.cos(ak)
    cc, cs = ci * ck, ci * sk
    sc, ss = si * ck, si * sk

    M = np.identity(4)
    if repetition:
        M[i, i] = cj
        M[i, j] = sj * si
        M[i, k] = sj * ci
        M[j, i] = sj * sk
        M[j, j] = -cj * ss + cc
        M[j, k] = -cj * cs - sc
        M[k, i] = -sj * ck
        M[k, j] = cj * sc + cs
        M[k, k] = cj * cc - ss
    else:
        M[i, i] = cj * ck
        M[i, j] = sj * sc - cs
        M[i, k] = sj * cc + ss
        M[j, i] = cj * sk
        M[j, j] = sj * ss + cc
        M[j, k] = sj * cs - sc
        M[k, i] = -sj
        M[k, j] = cj * si
        M[k, k] = cj * ci
    return M


def euler_from_matrix(matrix, axes='sxyz'):
    """Return Euler angles from rotation matrix for specified axis sequence.

    axes : One of 24 axis sequences as string or encoded tuple

    Note that many Euler angle triplets can describe one matrix.

    >>> from MDAnalysis.lib.transformations import (euler_matrix,
    ... euler_from_matrix, _AXES2TUPLE)
    >>> import math
    >>> import numpy as np
    >>> R0 = euler_matrix(1, 2, 3, 'syxz')
    >>> al, be, ga = euler_from_matrix(R0, 'syxz')
    >>> R1 = euler_matrix(al, be, ga, 'syxz')
    >>> np.allclose(R0, R1)
    True
    >>> angles = (4.0*math.pi) * (np.random.random(3) - 0.5)
    >>> for axes in _AXES2TUPLE.keys():
    ...    R0 = euler_matrix(axes=axes, *angles)
    ...    R1 = euler_matrix(axes=axes, *euler_from_matrix(R0, axes))
    ...    if not np.allclose(R0, R1): print(axes, "failed")

    """
    try:
        firstaxis, parity, repetition, frame = _AXES2TUPLE[axes.lower()]
    except (AttributeError, KeyError):
        _ = _TUPLE2AXES[axes]
        firstaxis, parity, repetition, frame = axes

    i = firstaxis
    j = _NEXT_AXIS[i + parity]
    k = _NEXT_AXIS[i - parity + 1]

    M = np.array(matrix, dtype=np.float64, copy=False)[:3, :3]
    if repetition:
        sy = math.sqrt(M[i, j] * M[i, j] + M[i, k] * M[i, k])
        if sy > _EPS:
            ax = math.atan2(M[i, j], M[i, k])
            ay = math.atan2(sy, M[i, i])
            az = math.atan2(M[j, i], -M[k, i])
        else:
            ax = math.atan2(-M[j, k], M[j, j])
            ay = math.atan2(sy, M[i, i])
            az = 0.0
    else:
        cy = math.sqrt(M[i, i] * M[i, i] + M[j, i] * M[j, i])
        if cy > _EPS:
            ax = math.atan2(M[k, j], M[k, k])
            ay = math.atan2(-M[k, i], cy)
            az = math.atan2(M[j, i], M[i, i])
        else:
            ax = math.atan2(-M[j, k], M[j, j])
            ay = math.atan2(-M[k, i], cy)
            az = 0.0

    if parity:
        ax, ay, az = -ax, -ay, -az
    if frame:
        ax, az = az, ax
    return ax, ay, az


def euler_from_quaternion(quaternion, axes='sxyz'):
    """Return Euler angles from quaternion for specified axis sequence.

    >>> from MDAnalysis.lib.transformations import euler_from_quaternion
    >>> import numpy as np
    >>> angles = euler_from_quaternion([0.99810947, 0.06146124, 0, 0])
    >>> np.allclose(angles, [0.123, 0, 0])
    True

    """
    return euler_from_matrix(quaternion_matrix(quaternion), axes)


def quaternion_from_euler(ai, aj, ak, axes='sxyz'):
    """Return quaternion from Euler angles and axis sequence.

    ai, aj, ak : Euler's roll, pitch and yaw angles
    axes : One of 24 axis sequences as string or encoded tuple

    >>> from MDAnalysis.lib.transformations import quaternion_from_euler
    >>> q = quaternion_from_euler(1, 2, 3, 'ryxz')
    >>> np.allclose(q, [0.435953, 0.310622, -0.718287, 0.444435])
    True

    """
    try:
        firstaxis, parity, repetition, frame = _AXES2TUPLE[axes.lower()]
    except (AttributeError, KeyError):
        _ = _TUPLE2AXES[axes]
        firstaxis, parity, repetition, frame = axes

    i = firstaxis + 1
    j = _NEXT_AXIS[i + parity - 1] + 1
    k = _NEXT_AXIS[i - parity] + 1

    if frame:
        ai, ak = ak, ai
    if parity:
        aj = -aj

    ai /= 2.0
    aj /= 2.0
    ak /= 2.0
    ci = math.cos(ai)
    si = math.sin(ai)
    cj = math.cos(aj)
    sj = math.sin(aj)
    ck = math.cos(ak)
    sk = math.sin(ak)
    cc = ci * ck
    cs = ci * sk
    sc = si * ck
    ss = si * sk

    quaternion = np.empty((4, ), dtype=np.float64)
    if repetition:
        quaternion[0] = cj * (cc - ss)
        quaternion[i] = cj * (cs + sc)
        quaternion[j] = sj * (cc + ss)
        quaternion[k] = sj * (cs - sc)
    else:
        quaternion[0] = cj * cc + sj * ss
        quaternion[i] = cj * sc - sj * cs
        quaternion[j] = cj * ss + sj * cc
        quaternion[k] = cj * cs - sj * sc
    if parity:
        quaternion[j] *= -1

    return quaternion


def quaternion_about_axis(angle, axis):
    """Return quaternion for rotation about axis.

    >>> from MDAnalysis.lib.transformations import quaternion_about_axis
    >>> import numpy as np
    >>> q = quaternion_about_axis(0.123, (1, 0, 0))
    >>> np.allclose(q, [0.99810947, 0.06146124, 0, 0])
    True

    """
    quaternion = np.zeros((4, ), dtype=np.float64)
    quaternion[1] = axis[0]
    quaternion[2] = axis[1]
    quaternion[3] = axis[2]
    qlen = vector_norm(quaternion)
    if qlen > _EPS:
        quaternion *= math.sin(angle / 2.0) / qlen
    quaternion[0] = math.cos(angle / 2.0)
    return quaternion


def quaternion_matrix(quaternion):
    """Return homogeneous rotation matrix from quaternion.

    >>> from MDAnalysis.lib.transformations import (identity_matrix,
    ...     quaternion_matrix, rotation_matrix)
    >>> import numpy as np
    >>> M = quaternion_matrix([0.99810947, 0.06146124, 0, 0])
    >>> np.allclose(M, rotation_matrix(0.123, (1, 0, 0)))
    True
    >>> M = quaternion_matrix([1, 0, 0, 0])
    >>> np.allclose(M, identity_matrix())
    True
    >>> M = quaternion_matrix([0, 1, 0, 0])
    >>> np.allclose(M, np.diag([1, -1, -1, 1]))
    True

    """
    q = np.array(quaternion[:4], dtype=np.float64, copy=True)
    nq = np.dot(q, q)
    if nq < _EPS:
        return np.identity(4)
    q *= math.sqrt(2.0 / nq)
    q = np.outer(q, q)
    return np.array(
        (
            (1.0 - q[2, 2] - q[3, 3], q[1, 2] - q[3, 0], q[1, 3] + q[2, 0], 0.0),
            (q[1, 2] + q[3, 0], 1.0 - q[1, 1] - q[3, 3], q[2, 3] - q[1, 0], 0.0),
            (q[1, 3] - q[2, 0], q[2, 3] + q[1, 0], 1.0 - q[1, 1] - q[2, 2], 0.0),
            (0.0, 0.0, 0.0, 1.0)
        ), dtype=np.float64)


def quaternion_from_matrix(matrix, isprecise=False):
    """Return quaternion from rotation matrix.

    If isprecise=True, the input matrix is assumed to be a precise rotation
    matrix and a faster algorithm is used.

    >>> from MDAnalysis.lib.transformations import (identity_matrix,
    ...     quaternion_from_matrix, rotation_matrix, random_rotation_matrix,
    ...     is_same_transform, quaternion_matrix)
    >>> import numpy as np
    >>> q = quaternion_from_matrix(identity_matrix(), True)
    >>> np.allclose(q, [1., 0., 0., 0.])
    True
    >>> q = quaternion_from_matrix(np.diag([1., -1., -1., 1.]))
    >>> np.allclose(q, [0, 1, 0, 0]) or np.allclose(q, [0, -1, 0, 0])
    True
    >>> R = rotation_matrix(0.123, (1, 2, 3))
    >>> q = quaternion_from_matrix(R, True)
    >>> np.allclose(q, [0.9981095, 0.0164262, 0.0328524, 0.0492786])
    True
    >>> R = [[-0.545, 0.797, 0.260, 0], [0.733, 0.603, -0.313, 0],
    ...      [-0.407, 0.021, -0.913, 0], [0, 0, 0, 1]]
    >>> q = quaternion_from_matrix(R)
    >>> np.allclose(q, [0.19069, 0.43736, 0.87485, -0.083611])
    True
    >>> R = [[0.395, 0.362, 0.843, 0], [-0.626, 0.796, -0.056, 0],
    ...      [-0.677, -0.498, 0.529, 0], [0, 0, 0, 1]]
    >>> q = quaternion_from_matrix(R)
    >>> np.allclose(q, [0.82336615, -0.13610694, 0.46344705, -0.29792603])
    True
    >>> R = random_rotation_matrix()
    >>> q = quaternion_from_matrix(R)
    >>> is_same_transform(R, quaternion_matrix(q))
    True

    """
    M = np.array(matrix, dtype=np.float64, copy=no_copy_shim)[:4, :4]
    if isprecise:
        q = np.empty((4, ), dtype=np.float64)
        t = np.trace(M)
        if t > M[3, 3]:
            q[0] = t
            q[3] = M[1, 0] - M[0, 1]
            q[2] = M[0, 2] - M[2, 0]
            q[1] = M[2, 1] - M[1, 2]
        else:
            i, j, k = 1, 2, 3
            if M[1, 1] > M[0, 0]:
                i, j, k = 2, 3, 1
            if M[2, 2] > M[i, i]:
                i, j, k = 3, 1, 2
            t = M[i, i] - (M[j, j] + M[k, k]) + M[3, 3]
            q[i] = t
            q[j] = M[i, j] + M[j, i]
            q[k] = M[k, i] + M[i, k]
            q[3] = M[k, j] - M[j, k]
        q *= 0.5 / math.sqrt(t * M[3, 3])
    else:
        m00 = M[0, 0]
        m01 = M[0, 1]
        m02 = M[0, 2]
        m10 = M[1, 0]
        m11 = M[1, 1]
        m12 = M[1, 2]
        m20 = M[2, 0]
        m21 = M[2, 1]
        m22 = M[2, 2]
        # symmetric matrix K
        K = np.array((
            (m00 - m11 - m22, 0.0, 0.0, 0.0),
            (m01 + m10, m11 - m00 - m22, 0.0, 0.0),
            (m02 + m20, m12 + m21, m22 - m00 - m11, 0.0),
            (m21 - m12, m02 - m20, m10 - m01, m00 + m11 + m22)))
        K /= 3.0
        # quaternion is eigenvector of K that corresponds to largest eigenvalue
        l, V = np.linalg.eigh(K)
        q = V[[3, 0, 1, 2], np.argmax(l)]

    if q[0] < 0.0:
        q *= -1.0
    return q


def quaternion_multiply(quaternion1, quaternion0):
    """Return multiplication of two quaternions.

    >>> from MDAnalysis.lib.transformations import quaternion_multiply
    >>> import numpy as np
    >>> q = quaternion_multiply([4, 1, -2, 3], [8, -5, 6, 7])
    >>> np.allclose(q, [28, -44, -14, 48])
    True

    """
    w0, x0, y0, z0 = quaternion0
    w1, x1, y1, z1 = quaternion1
    return np.array(
        (
            -x1 * x0 - y1 * y0 - z1 * z0 + w1 * w0,
            x1 * w0 + y1 * z0 - z1 * y0 + w1 * x0,
            -x1 * z0 + y1 * w0 + z1 * x0 + w1 * y0,
            x1 * y0 - y1 * x0 + z1 * w0 + w1 * z0), dtype=np.float64)


def quaternion_conjugate(quaternion):
    """Return conjugate of quaternion.

    >>> from MDAnalysis.lib.transformations import (random_quaternion,
    ...     quaternion_conjugate)
    >>> import numpy as np
    >>> q0 = random_quaternion()
    >>> q1 = quaternion_conjugate(q0)
    >>> q1[0] == q0[0] and all(q1[1:] == -q0[1:])
    True

    """
    return np.array(
        (
            quaternion[0], -quaternion[1],
            -quaternion[2], -quaternion[3]), dtype=np.float64)


def quaternion_inverse(quaternion):
    """Return inverse of quaternion.

    >>> from MDAnalysis.lib.transformations import (random_quaternion,
    ...     quaternion_inverse)
    >>> import numpy as np
    >>> q0 = random_quaternion()
    >>> q1 = quaternion_inverse(q0)
    >>> np.allclose(quaternion_multiply(q0, q1), [1, 0, 0, 0])
    True

    """
    return quaternion_conjugate(quaternion) / np.dot(quaternion, quaternion)


def quaternion_real(quaternion):
    """Return real part of quaternion.

    >>> from MDAnalysis.lib.transformations import quaternion_real
    >>> quaternion_real([3.0, 0.0, 1.0, 2.0])
    3.0

    """
    return quaternion[0]


def quaternion_imag(quaternion):
    """Return imaginary part of quaternion.

    >>> from MDAnalysis.lib.transformations import quaternion_imag
    >>> quaternion_imag([3.0, 0.0, 1.0, 2.0])
    [0.0, 1.0, 2.0]

    """
    return quaternion[1:4]


def quaternion_slerp(quat0, quat1, fraction, spin=0, shortestpath=True):
    r"""Return spherical linear interpolation between two quaternions.

    >>> from MDAnalysis.lib.transformations import (random_quaternion,
    ...     quaternion_slerp)
    >>> import math
    >>> import numpy as np
    >>> q0 = random_quaternion()
    >>> q1 = random_quaternion()
    >>> q = quaternion_slerp(q0, q1, 0.0)
    >>> np.allclose(q, q0)
    True
    >>> q = quaternion_slerp(q0, q1, 1.0, 1)
    >>> np.allclose(q, q1)
    True
    >>> q = quaternion_slerp(q0, q1, 0.5)
    >>> angle = math.acos(np.dot(q0, q))
    >>> np.allclose(2.0, math.acos(np.dot(q0, q1)) / angle) or \
    ... np.allclose(2.0, math.acos(-np.dot(q0, q1)) / angle)
    True

    """
    q0 = unit_vector(quat0[:4])
    q1 = unit_vector(quat1[:4])
    if fraction == 0.0:
        return q0
    elif fraction == 1.0:
        return q1
    d = np.dot(q0, q1)
    if abs(abs(d) - 1.0) < _EPS:
        return q0
    if shortestpath and d < 0.0:
        # invert rotation
        d = -d
        q1 *= -1.0
    angle = math.acos(d) + spin * math.pi
    if abs(angle) < _EPS:
        return q0
    isin = 1.0 / math.sin(angle)
    q0 *= math.sin((1.0 - fraction) * angle) * isin
    q1 *= math.sin(fraction * angle) * isin
    q0 += q1
    return q0


def random_quaternion(rand=None):
    """Return uniform random unit quaternion.

    rand: array like or None
        Three independent random variables that are uniformly distributed
        between 0 and 1.

    >>> from MDAnalysis.lib.transformations import (random_quaternion,
    ...     vector_norm)
    >>> import numpy as np
    >>> q = random_quaternion()
    >>> np.allclose(1.0, vector_norm(q))
    True
    >>> q = random_quaternion(np.random.random(3))
    >>> len(q.shape), q.shape[0]==4
    (1, True)

    """
    if rand is None:
        rand = np.random.rand(3)
    else:
        assert len(rand) == 3
    r1 = np.sqrt(1.0 - rand[0])
    r2 = np.sqrt(rand[0])
    pi2 = math.pi * 2.0
    t1 = pi2 * rand[1]
    t2 = pi2 * rand[2]
    return np.array(
        (
            np.cos(t2) * r2,
            np.sin(t1) * r1,
            np.cos(t1) * r1,
            np.sin(t2) * r2), dtype=np.float64)


def random_rotation_matrix(rand=None):
    """Return uniform random rotation matrix.

    rnd: array like
        Three independent random variables that are uniformly distributed
        between 0 and 1 for each returned quaternion.

    >>> from MDAnalysis.lib.transformations import random_rotation_matrix
    >>> import numpy as np
    >>> R = random_rotation_matrix()
    >>> np.allclose(np.dot(R.T, R), np.identity(4))
    True

    """
    return quaternion_matrix(random_quaternion(rand))


class Arcball(object):
    """Virtual Trackball Control.

    >>> from MDAnalysis.lib.transformations import Arcball
    >>> import numpy as np
    >>> ball = Arcball()
    >>> ball = Arcball(initial=np.identity(4))
    >>> ball.place([320, 320], 320)
    >>> ball.down([500, 250])
    >>> ball.drag([475, 275])
    >>> R = ball.matrix()
    >>> np.allclose(np.sum(R), 3.90583455)
    True
    >>> ball = Arcball(initial=[1, 0, 0, 0])
    >>> ball.place([320, 320], 320)
    >>> ball.setaxes([1,1,0], [-1, 1, 0])
    >>> ball.setconstrain(True)
    >>> ball.down([400, 200])
    >>> ball.drag([200, 400])
    >>> R = ball.matrix()
    >>> np.allclose(np.sum(R), 0.2055924)
    True
    >>> ball.next()

    """

    def __init__(self, initial=None):
        """Initialize virtual trackball control.

        initial : quaternion or rotation matrix

        """
        self._axis = None
        self._axes = None
        self._radius = 1.0
        self._center = [0.0, 0.0]
        self._vdown = np.array([0, 0, 1], dtype=np.float64)
        self._constrain = False

        if initial is None:
            self._qdown = np.array([1, 0, 0, 0], dtype=np.float64)
        else:
            initial = np.array(initial, dtype=np.float64)
            if initial.shape == (4, 4):
                self._qdown = quaternion_from_matrix(initial)
            elif initial.shape == (4, ):
                initial /= vector_norm(initial)
                self._qdown = initial
            else:
                raise ValueError("initial not a quaternion or matrix")

        self._qnow = self._qpre = self._qdown

    def place(self, center, radius):
        """Place Arcball, e.g. when window size changes.

        center : sequence[2]
            Window coordinates of trackball center.
        radius : float
            Radius of trackball in window coordinates.

        """
        self._radius = float(radius)
        self._center[0] = center[0]
        self._center[1] = center[1]

    def setaxes(self, *axes):
        """Set axes to constrain rotations."""
        if axes is None:
            self._axes = None
        else:
            self._axes = [unit_vector(axis) for axis in axes]

    def setconstrain(self, constrain):
        """Set state of constrain to axis mode."""
        self._constrain = constrain is True

    def getconstrain(self):
        """Return state of constrain to axis mode."""
        return self._constrain

    def down(self, point):
        """Set initial cursor window coordinates and pick constrain-axis."""
        self._vdown = arcball_map_to_sphere(point, self._center, self._radius)
        self._qdown = self._qpre = self._qnow

        if self._constrain and self._axes is not None:
            self._axis = arcball_nearest_axis(self._vdown, self._axes)
            self._vdown = arcball_constrain_to_axis(self._vdown, self._axis)
        else:
            self._axis = None

    def drag(self, point):
        """Update current cursor window coordinates."""
        vnow = arcball_map_to_sphere(point, self._center, self._radius)

        if self._axis is not None:
            vnow = arcball_constrain_to_axis(vnow, self._axis)

        self._qpre = self._qnow

        t = np.cross(self._vdown, vnow)
        if np.dot(t, t) < _EPS:
            self._qnow = self._qdown
        else:
            q = [np.dot(self._vdown, vnow), t[0], t[1], t[2]]
            self._qnow = quaternion_multiply(q, self._qdown)

    def next(self, acceleration=0.0):
        """Continue rotation in direction of last drag."""
        q = quaternion_slerp(self._qpre, self._qnow, 2.0 + acceleration, False)
        self._qpre, self._qnow = self._qnow, q

    def matrix(self):
        """Return homogeneous rotation matrix."""
        return quaternion_matrix(self._qnow)


def arcball_map_to_sphere(point, center, radius):
    """Return unit sphere coordinates from window coordinates."""
    v = np.array(
        (
            (point[0] - center[0]) / radius,
            (center[1] - point[1]) / radius,
            0.0
        ), dtype=np.float64
    )
    n = v[0] * v[0] + v[1] * v[1]
    if n > 1.0:
        v /= math.sqrt(n)  # position outside of sphere
    else:
        v[2] = math.sqrt(1.0 - n)
    return v


def arcball_constrain_to_axis(point, axis):
    """Return sphere point perpendicular to axis."""
    v = np.array(point, dtype=np.float64, copy=True)
    a = np.array(axis, dtype=np.float64, copy=True)
    v -= a * np.dot(a, v)  # on plane
    n = vector_norm(v)
    if n > _EPS:
        if v[2] < 0.0:
            v *= -1.0
        v /= n
        return v
    if a[2] == 1.0:
        return np.array([1, 0, 0], dtype=np.float64)
    return unit_vector([-a[1], a[0], 0])


def arcball_nearest_axis(point, axes):
    """Return axis, which arc is nearest to point."""
    point = np.array(point, dtype=np.float64, copy=False)
    nearest = None
    mx = -1.0
    for axis in axes:
        t = np.dot(arcball_constrain_to_axis(point, axis), point)
        if t > mx:
            nearest = axis
            mx = t
    return nearest


# epsilon for testing whether a number is close to zero
_EPS = np.finfo(float).eps * 4.0

# axis sequences for Euler angles
_NEXT_AXIS = [1, 2, 0, 1]

# map axes strings to/from tuples of inner axis, parity, repetition, frame
_AXES2TUPLE = {
    'sxyz': (0, 0, 0, 0), 'sxyx': (0, 0, 1, 0), 'sxzy': (0, 1, 0, 0),
    'sxzx': (0, 1, 1, 0), 'syzx': (1, 0, 0, 0), 'syzy': (1, 0, 1, 0),
    'syxz': (1, 1, 0, 0), 'syxy': (1, 1, 1, 0), 'szxy': (2, 0, 0, 0),
    'szxz': (2, 0, 1, 0), 'szyx': (2, 1, 0, 0), 'szyz': (2, 1, 1, 0),
    'rzyx': (0, 0, 0, 1), 'rxyx': (0, 0, 1, 1), 'ryzx': (0, 1, 0, 1),
    'rxzx': (0, 1, 1, 1), 'rxzy': (1, 0, 0, 1), 'ryzy': (1, 0, 1, 1),
    'rzxy': (1, 1, 0, 1), 'ryxy': (1, 1, 1, 1), 'ryxz': (2, 0, 0, 1),
    'rzxz': (2, 0, 1, 1), 'rxyz': (2, 1, 0, 1), 'rzyz': (2, 1, 1, 1)}

_TUPLE2AXES = dict((v, k) for k, v in _AXES2TUPLE.items())


def vector_norm(data, axis=None, out=None):
    """Return length, i.e. eucledian norm, of ndarray along axis.

    >>> from MDAnalysis.lib.transformations import vector_norm
    >>> import numpy as np
    >>> v = np.random.random(3)
    >>> n = vector_norm(v)
    >>> np.allclose(n, np.linalg.norm(v))
    True
    >>> v = np.random.rand(6, 5, 3)
    >>> n = vector_norm(v, axis=-1)
    >>> np.allclose(n, np.sqrt(np.sum(v*v, axis=2)))
    True
    >>> n = vector_norm(v, axis=1)
    >>> np.allclose(n, np.sqrt(np.sum(v*v, axis=1)))
    True
    >>> v = np.random.rand(5, 4, 3)
    >>> n = np.empty((5, 3), dtype=np.float64)
    >>> vector_norm(v, axis=1, out=n)
    >>> np.allclose(n, np.sqrt(np.sum(v*v, axis=1)))
    True
    >>> vector_norm([])
    0.0
    >>> vector_norm([1.0])
    1.0

    """
    data = np.array(data, dtype=np.float64, copy=True)
    if out is None:
        if data.ndim == 1:
            return math.sqrt(np.dot(data, data))
        data *= data
        out = np.atleast_1d(np.sum(data, axis=axis))
        np.sqrt(out, out)
        return out
    else:
        data *= data
        np.sum(data, axis=axis, out=out)
        np.sqrt(out, out)


def unit_vector(data, axis=None, out=None):
    """Return ndarray normalized by length, i.e. eucledian norm, along axis.

    >>> from MDAnalysis.lib.transformations import unit_vector
    >>> import numpy as np
    >>> v0 = np.random.random(3)
    >>> v1 = unit_vector(v0)
    >>> np.allclose(v1, v0 / np.linalg.norm(v0))
    True
    >>> v0 = np.random.rand(5, 4, 3)
    >>> v1 = unit_vector(v0, axis=-1)
    >>> v2 = v0 / np.expand_dims(np.sqrt(np.sum(v0*v0, axis=2)), 2)
    >>> np.allclose(v1, v2)
    True
    >>> v1 = unit_vector(v0, axis=1)
    >>> v2 = v0 / np.expand_dims(np.sqrt(np.sum(v0*v0, axis=1)), 1)
    >>> np.allclose(v1, v2)
    True
    >>> v1 = np.empty((5, 4, 3), dtype=np.float64)
    >>> unit_vector(v0, axis=1, out=v1)
    >>> np.allclose(v1, v2)
    True
    >>> list(unit_vector([]))
    []
    >>> list(unit_vector([1.0]))
    [1.0]

    """
    if out is None:
        data = np.array(data, dtype=np.float64, copy=True)
        if data.ndim == 1:
            data /= math.sqrt(np.dot(data, data))
            return data
    else:
        if out is not data:
            out[:] = np.array(data, copy=False)
        data = out
    length = np.atleast_1d(np.sum(data * data, axis))
    np.sqrt(length, length)
    if axis is not None:
        length = np.expand_dims(length, axis)
    data /= length
    if out is None:
        return data


def random_vector(size):
    """Return array of random doubles in the half-open interval [0.0, 1.0).

    >>> from MDAnalysis.lib.transformations import random_vector
    >>> import numpy as np
    >>> v = random_vector(10000)
    >>> np.all(v >= 0.0) and np.all(v < 1.0)
    True
    >>> v0 = random_vector(10)
    >>> v1 = random_vector(10)
    >>> np.any(v0 == v1)
    False

    """
    return np.random.random(size)


def inverse_matrix(matrix):
    """Return inverse of square transformation matrix.

    >>> from MDAnalysis.lib.transformations import (random_rotation_matrix,
    ...     inverse_matrix)
    >>> import numpy as np
    >>> M0 = random_rotation_matrix()
    >>> M1 = inverse_matrix(M0.T)
    >>> np.allclose(M1, np.linalg.inv(M0.T))
    True
    >>> for size in range(1, 7):
    ...     M0 = np.random.rand(size, size)
    ...     M1 = inverse_matrix(M0)
    ...     if not np.allclose(M1, np.linalg.inv(M0)): print(size)

    """
    return np.linalg.inv(matrix)


def concatenate_matrices(*matrices):
    """Return concatenation of series of transformation matrices.

    >>> from MDAnalysis.lib.transformations import concatenate_matrices
    >>> import numpy as np
    >>> M = np.random.rand(16).reshape((4, 4)) - 0.5
    >>> np.allclose(M, concatenate_matrices(M))
    True
    >>> np.allclose(np.dot(M, M.T), concatenate_matrices(M, M.T))
    True

    """
    M = np.identity(4)
    for i in matrices:
        M = np.dot(M, i)
    return M


def is_same_transform(matrix0, matrix1):
    """Return True if two matrices perform same transformation.

    >>> from MDAnalysis.lib.transformations import (is_same_transform,
    ...     random_rotation_matrix)
    >>> import numpy as np
    >>> is_same_transform(np.identity(4), np.identity(4))
    True
    >>> is_same_transform(np.identity(4), random_rotation_matrix())
    False

    """
    matrix0 = np.array(matrix0, dtype=np.float64, copy=True)
    matrix0 /= matrix0[3, 3]
    matrix1 = np.array(matrix1, dtype=np.float64, copy=True)
    matrix1 /= matrix1[3, 3]
    return np.allclose(matrix0, matrix1)


def _import_module(module_name, warn=True, prefix='_py_', ignore='_'):
    """Try import all public attributes from module into global namespace.

    Existing attributes with name clashes are renamed with prefix.
    Attributes starting with underscore are ignored by default.

    Return True on successful import.

    """
    sys.path.append(os.path.dirname(__file__))
    try:
        module = __import__(module_name)
    except ImportError:
        sys.path.pop()
        if warn:
            warnings.warn("failed to import module " + module_name)
    else:
        sys.path.pop()
        for attr in dir(module):
            if ignore and attr.startswith(ignore):
                continue
            if prefix:
                if attr in globals():
                    globals()[prefix + attr] = globals()[attr]
                elif warn:
                    warnings.warn("no Python implementation of " + attr)
            globals()[attr] = getattr(module, attr)
        return True


# orbeckst --- some simple geometry

def rotaxis(a, b):
    """Return the rotation axis to rotate vector a into b.

    Parameters
    ----------
    a, b : array_like
        two vectors

    Returns
    -------
    c : np.ndarray
        vector to rotate a into b


    Note
    ----
    If a == b this will always return [1, 0, 0]

    """
    if np.allclose(a, b):
        return np.array([1, 0, 0])
    c = np.cross(a, b)
    return c / np.linalg.norm(c)


_import_module('_transformations')

# Documentation in HTML format can be generated with Epydoc
__docformat__ = "restructuredtext en"
