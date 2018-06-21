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
    transformed = make_Universe(['masses'], trajectory=True)
    transformed.trajectory.ts.dimensions = np.array([372., 373., 374., 90, 90, 90])
    return reference, transformed

def test_translate_coords(translate_universes):
    ref_u = translate_universes[0]
    trans_u = translate_universes[1]
    ref = ref_u.trajectory.ts
    vector = np.float32([1, 2, 3])
    ref.positions += vector
    trans = translate(vector)(trans_u.trajectory.ts)
    assert_array_almost_equal(trans.positions, ref.positions, decimal=6)

def test_translate_vector(translate_universes):
    # what happens if the vector argument is smaller than 3?
    ts = translate_universes[0].trajectory.ts
    vector = [0, 1]
    with pytest.raises(ValueError):
        translate(vector)(ts)
    # what happens if the vector argument is bigger than 3?
    vector = [0, 1, 2, 3]
    with pytest.raises(ValueError):
        translate(vector)(ts)

def test_center_in_box_bad_ag(translate_universes):
    # this universe as a box size zero
    ts = translate_universes[0].trajectory.ts
    ag = translate_universes[0].residues[0].atoms
    # what happens if something other than an AtomGroup is given?
    bad_ag = 1
    with pytest.raises(ValueError): 
        center_in_box(bad_ag)(ts)

def test_center_in_box_bad_box(translate_universes):
    # this universe as a box size zero
    ts = translate_universes[0].trajectory.ts
    ag = translate_universes[0].residues[0].atoms
    # what if the box is in the wrong format?
    bad_box = [1, 2, 3]
    with pytest.raises(ValueError): 
        center_in_box(ag, box=bad_box)(ts)
    
def test_center_in_box_bad_box(translate_universes):    
    # this universe as a box size zero
    ts = translate_universes[0].trajectory.ts
    ag = translate_universes[0].residues[0].atoms
    # is pbc passed to the center methods?
    # if yes it should raise an exception for boxes that are zero in size
    with pytest.raises(ValueError): 
        center_in_box(ag, pbc=True)(ts)

def test_center_in_box_bad_box(translate_universes):
    # this universe as a box size zero
    ts = translate_universes[0].trajectory.ts
    ag = translate_universes[0].residues[0].atoms
    # what if a wrong center type name is passed?
    bad_center = " "
    with pytest.raises(ValueError): 
        center_in_box(ag, center=bad_center)(ts)
    
def test_center_in_box_bad_box(translate_universes):   
    # this universe as a box size zero
    ts = translate_universes[0].trajectory.ts
    ag = translate_universes[0].residues[0].atoms
    # if the universe has no masses and `mass` is passed as the center arg
    bad_center = "mass"
    with pytest.raises(AttributeError): 
        center_in_box(ag, center=bad_center)(ts)

def test_center_in_box_coords_no_options(translate_universes):
    # what happens when we center the coordinates arround the center of geometry of a residue?
    ref_u = translate_universes[0]
    trans_u = translate_universes[1]
    ref = ref_u.trajectory.ts
    ref_center = np.float32([6, 7, 8])
    box_center = np.float32([186., 186.5, 187.])
    ref.positions += box_center - ref_center
    ag = trans_u.residues[0].atoms
    trans = center_in_box(ag)(trans_u.trajectory.ts)
    assert_array_almost_equal(trans.positions, ref.positions, decimal=6)

def test_center_in_box_coords_with_pbc(translate_universes):
    # what happens when we center the coordinates arround the center of geometry of a residue?
    # using pbc into account for center of geometry calculation
    ref_u = translate_universes[0]
    ref = ref_u.trajectory.ts
    trans_u = translate_universes[1]
    trans_u.dimensions = [363., 364., 365., 90., 90., 90.]
    ag = trans_u.residues[24].atoms
    box_center = np.float32([181.5, 182., 182.5])
    ref_center = np.float32([75.6, 75.8, 76.])
    ref.positions += box_center - ref_center
    trans = center_in_box(ag, pbc=True)(trans_u.trajectory.ts)
    assert_array_almost_equal(trans.positions, ref.positions, decimal=6)

def test_center_in_box_coords_with_mass(translate_universes):   
    # using masses for calculating the center of the atomgroup
    ref_u = translate_universes[0]
    ref = ref_u.trajectory.ts
    trans_u = translate_universes[1]
    ag = trans_u.residues[24].atoms
    box_center = np.float32([186., 186.5, 187.])
    ref_center = ag.center_of_mass()
    ref.positions += box_center - ref_center
    trans = center_in_box(ag, center="mass")(trans_u.trajectory.ts)
    assert_array_almost_equal(trans.positions, ref.positions, decimal=6)

def test_center_in_box_coords_with_box(translate_universes):   
    # using masses for calculating the center of the atomgroup
    ref_u = translate_universes[0]
    ref = ref_u.trajectory.ts
    trans_u = translate_universes[1]
    ag = trans_u.residues[0].atoms
    newbox = [1000, 1000, 1000, 90, 90, 90]
    box_center = np.float32([500, 500, 500])
    ref_center = np.float32([6, 7, 8])
    ref.positions += box_center - ref_center
    trans = center_in_box(ag, box=newbox)(trans_u.trajectory.ts)
    assert_array_almost_equal(trans.positions, ref.positions, decimal=6)

def test_center_in_box_coords_all_options(translate_universes):
    # what happens when we center the coordinates arround the center of geometry of a residue?
    # using pbc into account for center of geometry calculation
    ref_u = translate_universes[0]
    ref = ref_u.trajectory.ts
    trans_u = translate_universes[1]
    ag = trans_u.residues[24].atoms
    newbox = [1000, 1000, 1000, 90, 90, 90]
    box_center = np.float32([500, 500, 500])
    ref_center = ag.center_of_mass(pbc=True)
    ref.positions += box_center - ref_center
    trans = center_in_box(ag, center='mass', pbc=True, box=newbox)(trans_u.trajectory.ts)
    assert_array_almost_equal(trans.positions, ref.positions, decimal=6)
