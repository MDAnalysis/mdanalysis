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

from __future__ import absolute_import

import pytest
import numpy as np
from numpy.testing import assert_almost_equal, assert_equal
from itertools import product

from MDAnalysis.lib._augment import augment, undo_augment
from MDAnalysis.lib.distances import apply_PBC, transform_StoR
from MDAnalysis.lib.mdamath import triclinic_vectors

boxes = ([10, 10, 10, 90, 90, 90],  # ortho
         [10, 10, 10, 45, 60, 90])  # tri_box

# Find images for several query points,
# here in fractional coordinates using augment
queries = ([0.1, 0.5, 0.5],  # box face
           [0.5, 0.5, 0.5],  # box center
           [0.5, -0.1, 0.5],  # box face
           [0.1, 0.1, 0.5],  # box edge
           [0.5, -0.1, 1.1],  # box edge
           [0.1, 0.1, 0.1],  # box vertex
           [0.1, -0.1, 1.1],  # box vertex
           [2.1, -3.1, 0.1]  # box vertex
           )

# Images for the previous query vectors, here in fractional coordinates.
images = (([1.1, 0.5, 0.5],),
          (),
          ([0.5, -0.1, 0.5],),
          ([1.1, 0.1, 0.5], [0.1, 1.1, 0.5], [1.1, 1.1, 0.5]),
          ([0.5, -0.1, 0.1], [0.5, 0.9, 1.1], [0.5, -0.1, 1.1]),
          ([1.1, 0.1, 0.1], [0.1, 1.1, 0.1], [0.1, 0.1, 1.1],
           [0.1, 1.1, 1.1], [1.1, 1.1, 0.1], [1.1, 0.1, 1.1],
           [1.1, 1.1, 1.1]),
          ([1.1, 0.9, 0.1], [0.1, -0.1, 0.1], [0.1, 0.9, 1.1],
           [0.1, -0.1, 1.1], [1.1, -0.1, 0.1], [1.1, 0.9, 1.1],
           [1.1, -0.1, 1.1]),
          ([1.1, 0.9, 0.1], [0.1, -0.1, 0.1], [0.1, 0.9, 1.1],
           [0.1, -0.1, 1.1], [1.1, -0.1, 0.1], [1.1, 0.9, 1.1],
           [1.1, -0.1, 1.1]))

radius = 1.5


@pytest.mark.parametrize('b, qres', product(boxes, zip(queries, images)))
def test_augment(b, qres):
    b = np.array(b, dtype=np.float32)
    q = transform_StoR(np.array(qres[0], dtype=np.float32), b)
    if q.shape == (3, ):
        q = q.reshape((1, 3))
    q = apply_PBC(q, b)
    dm = triclinic_vectors(b)
    aug, mapping = augment(q, dm, radius)
    if aug.size > 0:
        aug = np.sort(aug, axis=0)
    else:
        aug = list()
    cs = transform_StoR(np.array(qres[1], dtype=np.float32), b)
    if cs.size > 0:
        cs = np.sort(cs, axis=0)
    else:
        cs = list()
    assert_almost_equal(aug, cs, decimal=5)


@pytest.mark.parametrize('b, qres', product(boxes, zip(queries, images)))
def test_undoaugment(b, qres):
    b = np.array(b, dtype=np.float32)
    q = transform_StoR(np.array(qres[0], dtype=np.float32), b)
    if q.shape == (3, ):
        q = q.reshape((1, 3))
    q = apply_PBC(q, b)
    dm = triclinic_vectors(b)
    aug, mapping = augment(q, dm, radius)
    for idx, val in enumerate(aug):
        imageid = np.asarray([len(q) + idx], dtype=np.int32)
        assert_equal(mapping[idx], undo_augment(imageid, mapping, len(q))[0])
