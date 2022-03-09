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

import os
import pytest
import numpy as np
from numpy.testing import assert_almost_equal, assert_equal

from MDAnalysis.lib._augment import augment_coordinates, undo_augment
from MDAnalysis.lib.distances import apply_PBC, transform_StoR

# Find images for several query points,
# here in fractional coordinates
# Every element of qres tuple is (query, images)
qres = (
       ([0.1, 0.5, 0.5], [[1.1, 0.5, 0.5]]),     # box face
       ([0.5, 0.5, 0.5], []),                    # box center
       ([0.5, -0.1, 0.5], [[0.5, -0.1, 0.5]]),   # box face
       ([0.1, 0.1, 0.5], [[1.1, 0.1, 0.5],
                          [0.1, 1.1, 0.5],
                          [1.1, 1.1, 0.5]]),     # box edge
       ([0.5, -0.1, 1.1], [[0.5, -0.1, 0.1],
                           [0.5, 0.9, 1.1],
                           [0.5, -0.1, 1.1]]),   # box edge
       ([0.1, 0.1, 0.1], [[1.1, 0.1, 0.1],
                          [0.1, 1.1, 0.1],
                          [0.1, 0.1, 1.1],
                          [0.1, 1.1, 1.1],
                          [1.1, 1.1, 0.1],
                          [1.1, 0.1, 1.1],
                          [1.1, 1.1, 1.1]]),     # box vertex
       ([0.1, -0.1, 1.1], [[1.1, 0.9, 0.1],
                           [0.1, -0.1, 0.1],
                           [0.1, 0.9, 1.1],
                           [0.1, -0.1, 1.1],
                           [1.1, -0.1, 0.1],
                           [1.1, 0.9, 1.1],
                           [1.1, -0.1, 1.1]]),   # box vertex
       ([2.1, -3.1, 0.1], [[1.1, 0.9, 0.1],
                           [0.1, -0.1, 0.1],
                           [0.1, 0.9, 1.1],
                           [0.1, -0.1, 1.1],
                           [1.1, -0.1, 0.1],
                           [1.1, 0.9, 1.1],
                           [1.1, -0.1, 1.1]]),   # box vertex
       ([[0.1, 0.5, 0.5],
         [0.5, -0.1, 0.5]], [[1.1, 0.5, 0.5],
                             [0.5, -0.1, 0.5]])  # multiple queries
       )


@pytest.mark.xfail(os.name == "nt",
                   reason="see gh-3248")
@pytest.mark.parametrize('b', (
                         np.array([10, 10, 10, 90, 90, 90], dtype=np.float32),
                         np.array([10, 10, 10, 45, 60, 90], dtype=np.float32)
                         ))
@pytest.mark.parametrize('q, res', qres)
def test_augment(b, q, res):
    radius = 1.5
    q = transform_StoR(np.array(q, dtype=np.float32), b)
    if q.shape == (3, ):
        q = q.reshape((1, 3))
    q = apply_PBC(q, b)
    aug, mapping = augment_coordinates(q, b, radius)
    if aug.size > 0:
        aug = np.sort(aug, axis=0)
    else:
        aug = list()
    if len(res) > 0:
        cs = transform_StoR(np.array(res, dtype=np.float32), b)
        cs = np.sort(cs, axis=0)
    else:
        cs = list()
    assert_almost_equal(aug, cs, decimal=5)


@pytest.mark.parametrize('b', (
                         np.array([10, 10, 10, 90, 90, 90], dtype=np.float32),
                         np.array([10, 10, 10, 45, 60, 90], dtype=np.float32)
                         ))
@pytest.mark.parametrize('qres', qres)
def test_undoaugment(b, qres):
    radius = 1.5
    q = transform_StoR(np.array(qres[0], dtype=np.float32), b)
    if q.shape == (3, ):
        q = q.reshape((1, 3))
    q = apply_PBC(q, b)
    aug, mapping = augment_coordinates(q, b, radius)
    for idx, val in enumerate(aug):
        imageid = np.asarray([len(q) + idx], dtype=np.intp)
        assert_equal(mapping[idx], undo_augment(imageid, mapping, len(q))[0])
