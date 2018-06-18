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
import math
from numpy.testing import assert_array_almost_equal

import MDAnalysis as mda
from MDAnalysis.transformations import rotateby
from MDAnalysis.lib.transformations import rotation_matrix
from MDAnalysisTests.datafiles import GRO

def test_rotate_custom_position():
    u = mda.Universe(GRO)
    vector = np.float32([1,0,0])
    pos = np.float32([0,0,0])
    matrix = rotation_matrix(math.pi/2, vector, pos)[:3, :3]
    ref = u.trajectory.ts
    ref.positions = np.dot(ref.positions, matrix)
    transformed = rotateby(math.pi/2, vector, position = pos)(ref) 
    assert_array_almost_equal(transformed, ref.positions, decimal=6)
    
def test_rotate_atomgroup():
    u = mda.Universe(GRO)
    vector = np.float32([1,0,0])
    ref_pos = np.float32([52.953686, 44.179996, 29.397898])
    matrix = rotation_matrix(math.pi/2, vector, ref_pos)[:3, :3]
    ref = u.trajectory.ts
    ref.positions = np.dot(ref.positions, matrix)
    selection = u.atoms.select_atoms('resid 1')
    transformed = rotateby(math.pi/2, vector, ag = selection, center='geometry')(ref) 
    assert_array_almost_equal(transformed, ref.positions, decimal=6)
    
    