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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

from __future__ import absolute_import
import pytest
import numpy as np
from numpy.testing import assert_equal

import MDAnalysis as mda

@pytest.mark.parametrize('coord_dtype', (np.float32, np.float64))
def test_transform_StoR_pass(coord_dtype):
    box = np.array([10, 7, 3, 45, 60, 90], dtype=np.float32)
    s = np.array([[0.5, -0.1, 0.5]], dtype=coord_dtype)

    original_r = np.array([[ 5.75,  0.36066014, 0.75000012]], dtype=np.float32)

    test_r = mda.lib.distances.transform_StoR(s, box)

    assert_equal(original_r, test_r)


# different boxlengths to shift a coordinate
shifts = [
    ((0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0)), # no shifting
    ((1, 0, 0), (0, 1, 1), (0, 0, 1), (1, 1, 0)), # single box lengths
    ((-1, 0, 1), (0, -1, 0), (1, 0, 1), (-1, -1, -1)), # negative single
    ((4, 3, -2), (-2, 2, 2), (-5, 2, 2), (0, 2, 2)),  # multiple boxlengths
]


@pytest.mark.parametrize('shift', shifts)
@pytest.mark.parametrize('periodic', [True, False])
def test_calc_distance(shift, periodic):
    box = np.array([10, 10, 10, 90., 90., 90.], dtype=np.float32)

    coords = np.array([[1, 1, 1], [3, 1, 1]], dtype=np.float32)

    shift1, shift2, _, _ = shift

    coords[0] += shift1 * box[:3]
    coords[1] += shift2 * box[:3]

    box = box if periodic else None
    result = mda.lib.distances.calc_distance(coords[0], coords[1], box)

    reference = 2.0 if periodic else np.linalg.norm(coords[0] - coords[1])

    assert result == pytest.approx(reference)


@pytest.mark.parametrize('case', [
    # 90 degree angle
    (np.array([[1, 1, 1], [1, 2, 1], [2, 2, 1]], dtype=np.float32), 90.),
    # straight line / 180.
    (np.array([[1, 1, 1], [1, 2, 1], [1, 3, 1]], dtype=np.float32), 180.),
])
@pytest.mark.parametrize('shift', shifts)
@pytest.mark.parametrize('periodic', [True, False])
def test_calc_angle(case, shift, periodic):
    def manual_angle(x, y, z):
        return np.rad2deg(mda.lib.mdamath.angle(y - x, y - z))

    box = np.array([10, 20, 30, 90., 90., 90.], dtype=np.float32)
    (a, b, c), ref = case

    shift1, shift2, shift3, _ = shift

    a += shift1 * box[:3]
    b += shift2 * box[:3]
    c += shift3 * box[:3]

    box = box if periodic else None
    result = mda.lib.distances.calc_angle(a, b, c, box)

    reference = ref if periodic else manual_angle(a, b, c)

    assert result == pytest.approx(reference, abs=1e-3)


@pytest.mark.parametrize('case', [
    # 0 degree angle (cis)
    (np.array([[1, 2, 1], [1, 1, 1], [2, 1, 1], [2, 2, 1]], dtype=np.float32), 0.),
    # 180 degree (trans)
    (np.array([[1, 2, 1], [1, 1, 1], [2, 1, 1], [2, 0, 1]], dtype=np.float32), 180.),
])
@pytest.mark.parametrize('shift', shifts)
@pytest.mark.parametrize('periodic', [True, False])
def test_calc_dihedral(case, shift, periodic):
    def manual_dihedral(a, b, c, d):
        return np.rad2deg(mda.lib.mdamath.dihedral(b - a, c - b, d - c))

    box = np.array([10., 10., 10., 90., 90., 90.], dtype=np.float32)

    (a, b, c, d), ref = case

    shift1, shift2, shift3, shift4 = shift

    a += shift1 * box[:3]
    b += shift2 * box[:3]
    c += shift3 * box[:3]
    d += shift4 * box[:3]

    box = box if periodic else None
    result = mda.lib.distances.calc_dihedral(a, b, c, d, box)

    reference = ref if periodic else manual_dihedral(a, b, c, d)

    assert abs(result) == pytest.approx(abs(reference), abs=1e-3)
