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

import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal

import MDAnalysis as mda
from MDAnalysis.transformations import rotateby
from MDAnalysis.lib.transformations import rotation_matrix
from MDAnalysisTests import make_Universe

@pytest.fixture()
def rotateby_universes():
    # create the Universe objects for the tests
    reference = make_Universe(trajectory=True)
    transformed = make_Universe(trajectory=True)
    return reference, transformed

def test_rotateby_custom_position():
    # what happens when we use a custom position for the axis of rotation?
    ref_u, trans_u = rotateby_universes()
    trans = trans_u.trajectory.ts
    ref = ref_u.trajectory.ts
    vector = [1,0,0]
    pos = [0,0,0]
    matrix = rotation_matrix(np.pi / 2, vector, pos)[:3, :3]
    ref.positions = np.dot(ref.positions, matrix)
    transformed = rotateby(np.pi / 2, vector, position=pos)(trans)
    assert_array_almost_equal(transformed.positions, ref.positions, decimal=6)
    
def test_rotateby_atomgroup():
    # what happens when we rotate arround the center of geometry of a residue?
    ref_u, trans_u = rotateby_universes()
    trans = trans_u.trajectory.ts
    ref = ref_u.trajectory.ts
    center_pos = [6,7,8]
    vector = [1,0,0]
    matrix = rotation_matrix(np.pi, vector, center_pos)[:3, :3]
    ref.positions = np.dot(ref.positions, matrix)
    selection = trans_u.residues[0].atoms
    transformed = rotateby(np.pi, vector, ag=selection, center='geometry')(trans) 
    assert_array_almost_equal(transformed.positions, ref.positions, decimal=6)
    
    