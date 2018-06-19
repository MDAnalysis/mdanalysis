#-*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
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
from MDAnalysis.transformations import translate, center_in_box
from MDAnalysisTests import make_Universe

@pytest.fixture()
def translate_universes():
    # create the Universe objects for the tests
    reference = make_Universe(trajectory=True)
    transformed = make_Universe(trajectory=True)
    transformed.trajectory.ts.dimensions = np.array([372.,373.,374.,90,90,90])
    return reference, transformed

def test_translate():
    ref_u, trans_u = translate_universes()
    ref = ref_u.trajectory.ts
    vector = np.float32([1,2,3])
    ref.positions += vector
    trans = translate(vector)(trans_u.trajectory.ts)
    assert_array_almost_equal(trans.positions, ref.positions, decimal=6)
    

def test_center():
    # what happens when we center the coordinates arround the center of geometry of a residue?
    ref_u, trans_u = translate_universes()
    ref = ref_u.trajectory.ts
    ref_center = np.float32([6,7,8])
    box_center = np.float32([186.,186.5,187.])
    ref.positions += box_center - ref_center
    ag = trans_u.residues[0].atoms
    trans = center_in_box(ag)(trans_u.trajectory.ts)
    assert_array_almost_equal(trans.positions, ref.positions, decimal=6)
    
    