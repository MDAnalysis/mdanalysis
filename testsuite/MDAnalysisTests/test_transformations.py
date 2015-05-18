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

from MDAnalysis.core import transformations as t

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
    assert_allclose(I, t.rotation_matrix(np.pi*2, direc), atol=1e-06)
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
