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
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal

import MDAnalysis as mda
from MDAnalysis.transformations import translate, center_in_box
from MDAnalysisTests import make_Universe


@pytest.fixture()
def translate_universes():
    # create the Universe objects for the tests
    # this universe has no masses and some tests need it as such
    reference = make_Universe(trajectory=True)
    transformed = make_Universe(['masses'], trajectory=True)
    transformed.trajectory.ts.dimensions = np.array([372., 373., 374., 90, 90, 90])
    
    return reference, transformed


def test_translate_coords(translate_universes):
    ref_u, trans_u = translate_universes
    ref = ref_u.trajectory.ts
    vector = np.float32([1, 2, 3])
    ref.positions += vector
    trans = translate(vector)(trans_u.trajectory.ts)
    assert_array_almost_equal(trans.positions, ref.positions, decimal=6)


@pytest.mark.parametrize('vector', (
    [0, 1],
    [0, 1, 2, 3, 4],
    np.array([0, 1]),
    np.array([0, 1, 2, 3, 4]),
    np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]]),
    np.array([[0], [1], [2]]))
)
def test_translate_vector(translate_universes, vector):
    # what happens if the vector argument is of wrong size?
    ts = translate_universes[0].trajectory.ts
    with pytest.raises(ValueError):
        translate(vector)(ts)

      
def test_translate_transformations_api(translate_universes):
    # test if the translate transformation works when using the 
    # on-the-fly transformations API
    ref_u, trans_u = translate_universes
    ref = ref_u.trajectory.ts
    vector = np.float32([1, 2, 3])
    ref.positions += vector
    trans_u.trajectory.add_transformations(translate(vector))
    assert_array_almost_equal(trans_u.trajectory.ts.positions, ref.positions, decimal=6)


def test_center_in_box_bad_ag(translate_universes):
    # this universe has a box size zero
    ts = translate_universes[0].trajectory.ts
    # what happens if something other than an AtomGroup is given?
    bad_ag = 1
    with pytest.raises(ValueError): 
        center_in_box(bad_ag)(ts)


@pytest.mark.parametrize('point', (
    [0, 1],
    [0, 1, 2, 3, 4],
    np.array([0, 1]),
    np.array([0, 1, 2, 3, 4]),
    np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]]),
    np.array([[0], [1], [2]]))
)
def test_center_in_box_bad_point(translate_universes, point):
    ts = translate_universes[0].trajectory.ts
    ag = translate_universes[0].residues[0].atoms
    # what if the box is in the wrong format?
    with pytest.raises(ValueError): 
        center_in_box(ag, point=point)(ts)

   
def test_center_in_box_bad_pbc(translate_universes):    
    # this universe has a box size zero
    ts = translate_universes[0].trajectory.ts
    ag = translate_universes[0].residues[0].atoms
    # is pbc passed to the center methods?
    # if yes it should raise an exception for boxes that are zero in size
    with pytest.raises(ValueError): 
        center_in_box(ag, wrap=True)(ts)


def test_center_in_box_bad_center(translate_universes):
    # this universe has a box size zero
    ts = translate_universes[0].trajectory.ts
    ag = translate_universes[0].residues[0].atoms
    # what if a wrong center type name is passed?
    bad_center = " "
    with pytest.raises(ValueError): 
        center_in_box(ag, center=bad_center)(ts)


def test_center_in_box_no_masses(translate_universes):   
    # this universe has no masses
    ts = translate_universes[0].trajectory.ts
    ag = translate_universes[0].residues[0].atoms
    # if the universe has no masses and `mass` is passed as the center arg
    bad_center = "mass"
    with pytest.raises(mda.NoDataError):
        center_in_box(ag, center=bad_center)(ts)


def test_center_in_box_coords_no_options(translate_universes):
    # what happens when we center the coordinates arround the center of geometry of a residue?
    ref_u, trans_u = translate_universes
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
    ref_u, trans_u = translate_universes
    ref = ref_u.trajectory.ts
    trans_u.dimensions = [363., 364., 365., 90., 90., 90.]
    ag = trans_u.residues[24].atoms
    box_center = np.float32([181.5, 182., 182.5])
    ref_center = np.float32([75.6, 75.8, 76.])
    ref.positions += box_center - ref_center
    trans = center_in_box(ag, wrap=True)(trans_u.trajectory.ts)
    assert_array_almost_equal(trans.positions, ref.positions, decimal=6)


def test_center_in_box_coords_with_mass(translate_universes):   
    # using masses for calculating the center of the atomgroup
    ref_u, trans_u = translate_universes
    ref = ref_u.trajectory.ts
    ag = trans_u.residues[24].atoms
    box_center = np.float32([186., 186.5, 187.])
    ref_center = ag.center_of_mass()
    ref.positions += box_center - ref_center
    trans = center_in_box(ag, center="mass")(trans_u.trajectory.ts)
    assert_array_almost_equal(trans.positions, ref.positions, decimal=6)


def test_center_in_box_coords_with_box(translate_universes):   
    # using masses for calculating the center of the atomgroup
    ref_u, trans_u = translate_universes
    ref = ref_u.trajectory.ts
    ag = trans_u.residues[0].atoms
    newpoint = [1000, 1000, 1000]
    box_center = np.float32(newpoint)
    ref_center = np.float32([6, 7, 8])
    ref.positions += box_center - ref_center
    trans = center_in_box(ag, point=newpoint)(trans_u.trajectory.ts)
    assert_array_almost_equal(trans.positions, ref.positions, decimal=6)


def test_center_in_box_coords_all_options(translate_universes):
    # what happens when we center the coordinates arround the center of geometry of a residue?
    # using pbc into account for center of geometry calculation
    ref_u, trans_u = translate_universes
    ref = ref_u.trajectory.ts
    ag = trans_u.residues[24].atoms
    newpoint = [1000, 1000, 1000]
    box_center = np.float32(newpoint)
    ref_center = ag.center_of_mass(pbc=True)
    ref.positions += box_center - ref_center
    trans = center_in_box(ag, center='mass', wrap=True, point=newpoint)(trans_u.trajectory.ts)
    assert_array_almost_equal(trans.positions, ref.positions, decimal=6)


def test_center_transformations_api(translate_universes):
    # test if the translate transformation works when using the 
    # on-the-fly transformations API
    ref_u, trans_u = translate_universes
    ref = ref_u.trajectory.ts
    ref_center = np.float32([6, 7, 8])
    box_center = np.float32([186., 186.5, 187.])
    ref.positions += box_center - ref_center
    ag = trans_u.residues[0].atoms
    trans_u.trajectory.add_transformations(center_in_box(ag))
    assert_array_almost_equal(trans_u.trajectory.ts.positions, ref.positions, decimal=6)
