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
from MDAnalysis.transformations import translate, center_in_box, center_in_plane, center_in_axis
from MDAnalysisTests import make_Universe
from MDAnalysisTests.datafiles import fullerene


@pytest.fixture()
def translate_universes():
    # create the Universe objects for the tests
    # this universe has no masses and some tests need it as such
    reference = make_Universe(trajectory=True)
    transformed = make_Universe(['masses'], trajectory=True)
    transformed.trajectory.ts.dimensions = np.array([372., 373., 374., 90, 90, 90])
    
    return reference, transformed


@pytest.fixture()
def translate_unwrap_universes():
    # create the Universe objects for the tests
    # this universe is used for the unwrap testing cases
    reference = mda.Universe(fullerene)
    transformed = mda.Universe(fullerene)
    transformed.dimensions = np.asarray([10, 10, 10, 90, 90, 90], np.float32)
    transformed.atoms.wrap()
    
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


def test_translate_coords_dtype(translate_universes):
    # is the dtype of the coordinates correct? 
    trans_u = translate_universes[1]
    vector = np.float32([1, 2, 3])
    trans = translate(vector)(trans_u.trajectory.ts)
    assert(trans.positions.dtype == np.float32)


def test_translate_transformations_api(translate_universes):
    # test if the translate transformation works when using the 
    # on-the-fly transformations API
    ref_u, trans_u = translate_universes
    ref = ref_u.trajectory.ts
    vector = np.float32([1, 2, 3])
    ref.positions += vector
    trans_u.trajectory.add_transformations(translate(vector))
    assert_array_almost_equal(trans_u.trajectory.ts.positions, ref.positions, decimal=6)


@pytest.mark.parametrize('ag', (
    [0, 1],
    [0, 1, 2, 3, 4],
    np.array([0, 1]),
    np.array([0, 1, 2, 3, 4]),
    np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]]),
    np.array([[0], [1], [2]]),
    'thisisnotanag',
    1)
)
def test_center_in_box_bad_ag(translate_universes, ag):
    # this universe has a box size zero
    ts = translate_universes[0].trajectory.ts
    # what happens if something other than an AtomGroup is given?
    bad_ag = ag
    with pytest.raises(ValueError): 
        center_in_box(bad_ag)(ts)


@pytest.mark.parametrize('center_to', (
    [0, 1],
    [0, 1, 2, 3, 4],
    np.array([0, 1]),
    np.array([0, 1, 2, 3, 4]),
    np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]]),
    np.array([[0], [1], [2]]),
    'thisisnotacenter',
    1)
)
def test_center_in_box_bad_center_to(translate_universes, center_to):
    ts = translate_universes[0].trajectory.ts
    ag = translate_universes[0].residues[0].atoms
    bad_center_to = center_to
    # what if the box is in the wrong format?
    with pytest.raises(ValueError): 
        center_in_box(ag, center_to=bad_center_to)(ts)

   
def test_center_in_box_bad_wrap(translate_universes):    
    # this universe has a box size zero
    ts = translate_universes[0].trajectory.ts
    ag = translate_universes[0].residues[0].atoms
    # is wrap passed to the center methods?
    # if yes it should raise an exception for boxes that are zero in size
    with pytest.raises(ValueError): 
        center_in_box(ag, wrap=True)(ts)


def test_center_in_box_bad_unwrap(translate_universes):
    # this universe has a box size zero
    ts = translate_universes[0].trajectory.ts
    ag = translate_universes[0].residues[0].atoms
    # is wrap passed to the center methods?
    # if yes it should raise an exception for boxes that are zero in size
    with pytest.raises(ValueError): 
        center_in_box(ag, unwrap=True)(ts)
 

def test_center_in_box_bad_wrap_unwrap(translate_universes):
    # this universe has a box size zero
    ts = translate_universes[1].trajectory.ts
    ag = translate_universes[1].residues[0].atoms
    # the two options are not compatible so it should throw an error
    with pytest.raises(ValueError): 
        center_in_box(ag, wrap=True, unwrap=True)(ts)
       

@pytest.mark.parametrize('weights', (
    " ",
    "totallynotmasses",
    123456789,
    [0, 1, 2, 3, 4],
    np.array([0, 1]),
    np.array([0, 1, 2, 3, 4]),
    np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]]))
)
def test_center_in_box_bad_weights(translate_universes, weights):
    # this universe has a box size zero
    ts = translate_universes[0].trajectory.ts
    ag = translate_universes[0].residues[0].atoms
    # what if a wrong center type name is passed?
    bad_weights = " "
    with pytest.raises(TypeError): 
        center_in_box(ag, weights=bad_weights)(ts)


def test_center_in_box_no_masses(translate_universes):   
    # this universe has no masses
    ts = translate_universes[0].trajectory.ts
    ag = translate_universes[0].residues[0].atoms
    # if the universe has no masses and `mass` is passed as the center arg
    bad_center = "mass"
    with pytest.raises(TypeError): 
        center_in_box(ag, center_of=bad_center)(ts)


def test_center_in_box_coords_dtype(translate_universes):
    # is the dtype of the coordinates correct? 
    trans_u = translate_universes[1]
    ag = trans_u.residues[0].atoms
    trans = center_in_box(ag)(trans_u.trajectory.ts)
    assert(trans.positions.dtype == np.float32)
    

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


def test_center_in_box_coords_with_wrap(translate_universes):
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


def test_center_in_box_coords_with_unwrap(translate_unwrap_universes):
    # what happens when we center the coordinates arround the center of geometry of a residue?
    # using pbc into account for center of geometry calculation
    ref_u, trans_u = translate_unwrap_universes
    ref = ref_u.trajectory.ts
    ag = trans_u.atoms
    box_center = np.float32([5, 5, 5])
    ref_center = ref_u.atoms.center(weights=None)
    ref.positions += box_center - ref_center
    trans = center_in_box(ag, unwrap=True)(trans_u.trajectory.ts)
    assert_array_almost_equal(trans.positions, ref.positions, decimal=6)


def test_center_in_box_coords_with_mass(translate_universes):   
    # using masses for calculating the center of the atomgroup
    ref_u, trans_u = translate_universes
    ref = ref_u.trajectory.ts
    ag = trans_u.residues[24].atoms
    box_center = np.float32([186., 186.5, 187.])
    ref_center = ag.center(weights=ag.masses)
    ref.positions += box_center - ref_center
    trans = center_in_box(ag, weights="mass")(trans_u.trajectory.ts)
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
    trans = center_in_box(ag, center_to=newpoint)(trans_u.trajectory.ts)
    assert_array_almost_equal(trans.positions, ref.positions, decimal=6)


def test_center_in_box_coords_all_options_wrap(translate_universes):
    # what happens when we center the coordinates arround the center of geometry of a residue?
    # using pbc into account for center of geometry calculation
    ref_u, trans_u = translate_universes
    ref = ref_u.trajectory.ts
    ag = trans_u.residues[24].atoms
    newpoint = [1000, 1000, 1000]
    box_center = np.float32(newpoint)
    ref_center = ag.center(weights=ag.masses, pbc=True)
    ref.positions += box_center - ref_center
    trans = center_in_box(ag, weights='mass', wrap=True, center_to=newpoint)(trans_u.trajectory.ts)
    assert_array_almost_equal(trans.positions, ref.positions, decimal=6)


def test_center_in_box_coords_all_options_unwrap(translate_unwrap_universes):
    # what happens when we center the coordinates arround the center of geometry of a residue?
    # using pbc into account for center of geometry calculation
    ref_u, trans_u = translate_unwrap_universes
    ref = ref_u.trajectory.ts
    ag = trans_u.atoms
    newpoint = [1000, 1000, 1000]
    box_center = np.float32(newpoint)
    ref_center = ref_u.atoms.center(weights=ref_u.atoms.masses)
    ref.positions += box_center - ref_center
    trans = center_in_box(ag, weights='mass', unwrap=True, center_to=newpoint)(trans_u.trajectory.ts)
    assert_array_almost_equal(trans.positions, ref.positions, decimal=6)


@pytest.mark.parametrize('ag', (
    [0, 1],
    [0, 1, 2, 3, 4],
    np.array([0, 1]),
    np.array([0, 1, 2, 3, 4]),
    np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]]),
    np.array([[0], [1], [2]]),
    'thisisnotanag',
    1)
)
def test_center_in_axis_bad_ag(translate_universes, ag):
    # this universe has a box size zero
    ts = translate_universes[0].trajectory.ts
    # what happens if something other than an AtomGroup is given?
    bad_ag = ag
    axis = 'x'
    with pytest.raises(ValueError): 
        center_in_axis(bad_ag, axis)(ts)


@pytest.mark.parametrize('center_to', (
    [0, 1],
    [0, 1, 2, 3, 4],
    np.array([0, 1]),
    np.array([0, 1, 2, 3, 4]),
    np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]]),
    np.array([[0], [1], [2]]),
    'thisisnotacenter',
    1)
)
def test_center_in_axis_bad_center_to(translate_universes, center_to):
    ts = translate_universes[0].trajectory.ts
    ag = translate_universes[0].residues[0].atoms
    # what if the box is in the wrong format?
    axis = 'x'
    with pytest.raises(ValueError): 
        center_in_axis(ag, axis, center_to=center_to)(ts)


@pytest.mark.parametrize('axis', (
    1,
    [0, 1],
    [0, 1, 2, 3, 4],
    np.array([0, 1]),
    "xyz",
    "notanaxis")
)
def test_center_in_axis_bad_axis(translate_universes, axis):
    ts = translate_universes[0].trajectory.ts
    ag = translate_universes[0].residues[0].atoms
    # what if the box is in the wrong format?
    bad_axis = axis
    with pytest.raises(ValueError): 
        center_in_axis(ag, bad_axis)(ts)

   
def test_center_in_axis_bad_wrap(translate_universes):    
    # this universe has a box size zero
    ts = translate_universes[0].trajectory.ts
    ag = translate_universes[0].residues[0].atoms
    # is wrap passed to the center methods?
    # if yes it should raise an exception for boxes that are zero in size
    axis = 'x'
    with pytest.raises(ValueError): 
        center_in_axis(ag,axis, wrap=True)(ts)


def test_center_in_axis_bad_unwrap(translate_universes):    
    # this universe has a box size zero
    ts = translate_universes[0].trajectory.ts
    ag = translate_universes[0].residues[0].atoms
    # is unwrap passed to the center methods?
    # if yes it should raise an exception for boxes that are zero in size
    axis = 'x'
    with pytest.raises(ValueError): 
        center_in_axis(ag,axis, unwrap=True)(ts)


def test_center_in_axis_bad_wrap_unwrap(translate_universes):    
    # this universe has a box size zero
    ts = translate_universes[1].trajectory.ts
    ag = translate_universes[1].residues[0].atoms
    # the two options are not compatible so it should throw an error
    axis = 'x'
    with pytest.raises(ValueError): 
        center_in_axis(ag,axis, wrap=True, unwrap=True)(ts)


@pytest.mark.parametrize('weights', (
    " ",
    "totallynotmasses",
    123456789,
    [0, 1, 2, 3, 4],
    np.array([0, 1]),
    np.array([0, 1, 2, 3, 4]),
    np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]]))
)
def test_center_in_axis_bad_weights(translate_universes, weights):
    # this universe has a box size zero
    ts = translate_universes[0].trajectory.ts
    ag = translate_universes[0].residues[0].atoms
    # what if a wrong center type name is passed?
    bad_weights = " "
    axis = 'x'
    with pytest.raises(TypeError):
        center_in_axis(ag, axis, weights=bad_weights)(ts)


def test_center_in_axis_no_masses(translate_universes):   
    # this universe has no masses
    ts = translate_universes[0].trajectory.ts
    ag = translate_universes[0].residues[0].atoms
    # if the universe has no masses and `mass` is passed as the center arg
    bad_weights = "mass"
    axis = 'x'
    with pytest.raises(TypeError): 
        center_in_axis(ag, axis, weights=bad_weights)(ts)

def test_center_in_axis_coords_dtype(translate_universes):
    # is the dtype of the coordinates correct? 
    trans_u = translate_universes[1]
    axis = 'x'
    ag = trans_u.residues[0].atoms
    trans = center_in_axis(ag, axis)(trans_u.trajectory.ts)
    assert(trans.positions.dtype == np.float32)


def test_center_in_axis_coords_no_options(translate_universes):
    # what happens when we center the coordinates on an axis without any other options?
    ref_u, trans_u = translate_universes
    ref = ref_u.trajectory.ts
    ag = trans_u.residues[0].atoms
    box_center = np.float32([186., 186.5, 187.])
    center_to = box_center
    ref_center = np.float32([6, 7, 8])
    axis = 'x'
    axis_point = np.float32([ref_center[0], center_to[1], center_to[2]])
    ref.positions += axis_point - ref_center
    trans = center_in_axis(ag, axis)(trans_u.trajectory.ts)
    assert_array_almost_equal(trans.positions, ref.positions, decimal=6)


def test_center_in_axis_coords_with_wrap(translate_universes):
    # what happens when we center the coordinates arround the center of geometry of a residue
    # taking pbc into account for center of geometry calculation?
    ref_u, trans_u = translate_universes
    ref = ref_u.trajectory.ts
    trans_u.dimensions = [363., 364., 365., 90., 90., 90.]
    ag = trans_u.residues[24].atoms
    box_center = np.float32([181.5, 182., 182.5])
    center_to = box_center
    ref_center = ag.center_of_geometry(pbc=True)
    axis = 'x'
    axis_point = np.float32([ref_center[0], center_to[1], center_to[2]])
    ref.positions += axis_point - ref_center
    trans = center_in_axis(ag, axis, wrap=True)(trans_u.trajectory.ts)
    assert_array_almost_equal(trans.positions, ref.positions, decimal=6)


def test_center_in_axis_coords_with_unwrap(translate_unwrap_universes):
    # what happens when we center the coordinates arround the center of geometry of a residue
    # taking pbc into account for center of geometry calculation?
    ref_u, trans_u = translate_unwrap_universes
    ref = ref_u.trajectory.ts
    ag = trans_u.atoms
    box_center = np.float32([5, 5, 5])
    center_to = box_center
    ref_center = ref_u.atoms.center(weights=None)
    axis = 'x'
    axis_point = np.float32([ref_center[0], center_to[1], center_to[2]])
    ref.positions += axis_point - ref_center
    trans = center_in_axis(ag, axis, unwrap=True)(trans_u.trajectory.ts)
    assert_array_almost_equal(trans.positions, ref.positions, decimal=6)


def test_center_in_axis_coords_with_mass(translate_universes):   
    # using masses for calculating the center of the atomgroup
    ref_u, trans_u = translate_universes
    ref = ref_u.trajectory.ts
    ag = trans_u.residues[24].atoms
    ref_center = ag.center(weights=ag.masses)
    box_center = np.float32([186., 186.5, 187.])
    center_to = box_center
    axis = 'x'
    axis_point = np.float32([ref_center[0], center_to[1], center_to[2]])
    ref.positions += axis_point - ref_center
    trans = center_in_axis(ag, axis, weights="mass")(trans_u.trajectory.ts)
    assert_array_almost_equal(trans.positions, ref.positions, decimal=6)


def test_center_in_axis_coords_with_coord_center_to(translate_universes):   
    # what happens when we use a custom coordinate as center_to point?
    ref_u, trans_u = translate_universes
    ref = ref_u.trajectory.ts
    ag = trans_u.residues[0].atoms
    center_to = [1000, 1000, 1000]
    ref_center = np.float32([6, 7, 8])
    axis = 'x'
    axis_point = np.float32([ref_center[0], center_to[1], center_to[2]])
    ref.positions += axis_point - ref_center
    trans = center_in_axis(ag, axis, center_to=center_to)(trans_u.trajectory.ts)
    assert_array_almost_equal(trans.positions, ref.positions, decimal=6)


def test_center_in_axis_coords_all_options_wrap(translate_universes):
    # what happens when we center the coordinates arround the center of geometry of a residue?
    # using pbc into account for center of geometry calculation
    ref_u, trans_u = translate_universes
    ref = ref_u.trajectory.ts
    ag = trans_u.residues[24].atoms
    neworigin = [1000, 1000, 1000]
    ref_center = ag.center(weights=ag.masses, pbc=True)
    axis = 'x'
    axis_point = np.float32([ref_center[0], neworigin[1], neworigin[2]])
    ref.positions += axis_point - ref_center
    trans = center_in_axis(ag, axis, weights='mass', wrap=True, center_to=neworigin)(trans_u.trajectory.ts)
    assert_array_almost_equal(trans.positions, ref.positions, decimal=6)


def test_center_in_axis_coords_all_options_unwrap(translate_unwrap_universes):
    # what happens when we center the coordinates arround the center of geometry of a residue?
    # using pbc into account for center of geometry calculation
    ref_u, trans_u = translate_unwrap_universes
    ref = ref_u.trajectory.ts
    ag = trans_u.atoms
    neworigin = [1000, 1000, 1000]
    ref_center = ref_u.atoms.center(weights=ref_u.atoms.masses)
    axis = 'x'
    axis_point = np.float32([ref_center[0], neworigin[1], neworigin[2]])
    ref.positions += axis_point - ref_center
    trans = center_in_axis(ag, axis, weights='mass', unwrap=True, center_to=neworigin)(trans_u.trajectory.ts)
    assert_array_almost_equal(trans.positions, ref.positions, decimal=6)


@pytest.mark.parametrize('ag', (
    [0, 1],
    [0, 1, 2, 3, 4],
    np.array([0, 1]),
    np.array([0, 1, 2, 3, 4]),
    np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]]),
    np.array([[0], [1], [2]]),
    'thisisnotanag',
    1)
)
def test_center_in_plane_bad_ag(translate_universes, ag):
    # this universe has a box size zero
    ts = translate_universes[0].trajectory.ts
    # what happens if something other than an AtomGroup is given?
    bad_ag = ag
    plane = 'yz'
    with pytest.raises(ValueError):
        center_in_plane(bad_ag, plane)(ts)


@pytest.mark.parametrize('center_to', (
    [0, 1],
    [0, 1, 2, 3, 4],
    np.array([0, 1]),
    np.array([0, 1, 2, 3, 4]),
    np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]]),
    np.array([[0], [1], [2]]),
    'thisisnotacenter',
    1)
)
def test_center_in_plane_bad_center_to(translate_universes, center_to):
    ts = translate_universes[0].trajectory.ts
    ag = translate_universes[0].residues[0].atoms
    # what if the box is in the wrong format?
    plane = 'yz'
    with pytest.raises(ValueError): 
        center_in_plane(ag, plane, center_to=center_to)(ts)


@pytest.mark.parametrize('plane', (
    1,
    [0, 1],
    [0, 1, 2, 3, 4],
    np.array([0, 1]),
    "xyz",
    "notaplane")
)
def test_center_in_plane_bad_plane(translate_universes, plane):
    ts = translate_universes[0].trajectory.ts
    ag = translate_universes[0].residues[0].atoms
    # what if the box is in the wrong format?
    bad_plane = plane
    with pytest.raises(ValueError): 
        center_in_plane(ag, bad_plane)(ts)

   
def test_center_in_plane_bad_wrap(translate_universes):    
    # this universe has a box size zero
    ts = translate_universes[0].trajectory.ts
    ag = translate_universes[0].residues[0].atoms
    # is wrap passed to the center methods?
    # if yes it should raise an exception for boxes that are zero in size
    plane = 'yz'
    with pytest.raises(ValueError): 
        center_in_plane(ag, plane, wrap=True)(ts)


def test_center_in_plane_bad_unwrap(translate_universes):    
    # this universe has a box size zero
    ts = translate_universes[0].trajectory.ts
    ag = translate_universes[0].residues[0].atoms
    # is unwrap passed to the center methods?
    # if yes it should raise an exception for boxes that are zero in size
    plane = 'yz'
    with pytest.raises(ValueError): 
        center_in_plane(ag, plane, unwrap=True)(ts)


def test_center_in_plane_bad_weights(translate_universes):
    # this universe has a box size zero
    ts = translate_universes[0].trajectory.ts
    ag = translate_universes[0].residues[0].atoms
    # what if a wrong center type name is passed?
    bad_weights = " "
    plane = 'yz'
    with pytest.raises(TypeError):
        center_in_plane(ag, plane, weights=bad_weights)(ts)


def test_center_in_planes_no_masses(translate_universes):   
    # this universe has no masses
    ts = translate_universes[0].trajectory.ts
    ag = translate_universes[0].residues[0].atoms
    # if the universe has no masses and `mass` is passed as the center arg
    weights = "mass"
    plane = 'yz'
    with pytest.raises(TypeError): 
        center_in_plane(ag, plane, weights=weights)(ts)

def test_center_in_plane_coords_dtype(translate_universes):
    # is the dtype of the coordinates correct? 
    trans_u = translate_universes[1]
    plane = 'yz'
    ag = trans_u.residues[0].atoms
    trans = center_in_plane(ag, plane)(trans_u.trajectory.ts)
    assert(trans.positions.dtype == np.float32)
    
    
def test_center_in_plane_coords_no_options(translate_universes):
    # what happens when we center the coordinates on a plane without special options?
    ref_u, trans_u = translate_universes
    ref = ref_u.trajectory.ts
    ag = trans_u.residues[0].atoms
    box_center = np.float32([186., 186.5, 187.])
    ref_center = np.float32([6, 7, 8])
    plane = 'yz'
    plane_point = np.asarray([box_center[0], ref_center[1], ref_center[2]], np.float32)
    ref.positions += plane_point - ref_center
    trans = center_in_plane(ag, plane)(trans_u.trajectory.ts)
    assert_array_almost_equal(trans.positions, ref.positions, decimal=6)


def test_center_in_plane_coords_with_wrap(translate_universes):
    # what happens when we center the coordinates arround the center of geometry of a residue?
    # using pbc into account for center of geometry calculation
    ref_u, trans_u = translate_universes
    ref = ref_u.trajectory.ts
    trans_u.dimensions = [363., 364., 365., 90., 90., 90.]
    box_center = np.float32([181.5, 182., 182.5])
    ag = trans_u.residues[24].atoms
    ref_center = ag.center(weights=None, pbc=True)
    plane = 'yz'
    plane_point = np.asarray([box_center[0], ref_center[1], ref_center[2]], np.float32)
    ref.positions += plane_point - ref_center
    trans = center_in_plane(ag, plane, wrap=True)(trans_u.trajectory.ts)
    assert_array_almost_equal(trans.positions, ref.positions, decimal=6)


def test_center_in_plane_coords_with_unwrap(translate_unwrap_universes):
    # what happens when we center the coordinates arround the center of geometry of a residue?
    # using pbc into account for center of geometry calculation
    ref_u, trans_u = translate_unwrap_universes
    ref = ref_u.trajectory.ts
    box_center = np.float32([5, 5, 5])
    ag = trans_u.atoms
    ref_center = ref_u.atoms.center(weights=None)
    plane = 'yz'
    plane_point = np.asarray([box_center[0], ref_center[1], ref_center[2]], np.float32)
    ref.positions += plane_point - ref_center
    trans = center_in_plane(ag, plane, unwrap=True)(trans_u.trajectory.ts)
    assert_array_almost_equal(trans.positions, ref.positions, decimal=6)


def test_center_in_plane_coords_with_mass(translate_universes):   
    # using masses for calculating the center of the atomgroup
    ref_u, trans_u = translate_universes
    ref = ref_u.trajectory.ts
    ag = trans_u.residues[24].atoms
    box_center = np.float32([186., 186.5, 187.])
    ref_center = ag.center(weights=ag.masses)
    plane = 'yz'
    plane_point = np.asarray([box_center[0], ref_center[1], ref_center[2]], np.float32)
    ref.positions += plane_point - ref_center
    trans = center_in_plane(ag, plane, weights="mass")(trans_u.trajectory.ts)
    assert_array_almost_equal(trans.positions, ref.positions, decimal=6)


def test_center_in_plane_coords_with_coord_center_to(translate_universes):   
    # what happens when we use a custom coordinate as center_to point?
    ref_u, trans_u = translate_universes
    ref = ref_u.trajectory.ts
    ag = trans_u.residues[0].atoms
    center_to = [1000, 1000, 1000]
    box_center = np.float32([186., 186.5, 187.])
    ref_center = np.float32([6, 7, 8])
    plane = 'yz'
    plane_point = np.float32([center_to[0], ref_center[1], ref_center[2]])
    ref.positions += plane_point - ref_center
    trans = center_in_plane(ag, plane, center_to=center_to)(trans_u.trajectory.ts)
    assert_array_almost_equal(trans.positions, ref.positions, decimal=6)


def test_center_in_plane_coords_all_options_wrap(translate_universes):
    # what happens when we center the coordinates arround the center of geometry of a residue?
    # using pbc into account for center of geometry calculation
    ref_u, trans_u = translate_universes
    ref = ref_u.trajectory.ts
    ag = trans_u.residues[24].atoms
    trans_u.dimensions = [363., 364., 365., 90., 90., 90.]
    box_center = np.float32([181.5, 182., 182.5])
    center_to = [1000, 1000, 1000]
    ref_center = ag.center(weights=ag.masses, pbc=True)
    plane = 'yz'
    plane_point = np.float32([center_to[0], ref_center[1], ref_center[2]])
    ref.positions += plane_point - ref_center
    trans = center_in_plane(ag, plane, center_to=center_to, weights='mass', wrap=True)(trans_u.trajectory.ts)
    assert_array_almost_equal(trans.positions, ref.positions, decimal=6)


def test_center_in_plane_coords_all_options_unwrap(translate_unwrap_universes):
    # what happens when we center the coordinates arround the center of geometry of a residue?
    # using pbc into account for center of geometry calculation
    ref_u, trans_u = translate_unwrap_universes
    ref = ref_u.trajectory.ts
    ag = trans_u.atoms
    box_center = np.float32([5, 5, 5])
    center_to = [1000, 1000, 1000]
    ref_center = ref_u.atoms.center(weights=ref_u.atoms.masses)
    plane = 'yz'
    plane_point = np.float32([center_to[0], ref_center[1], ref_center[2]])
    ref.positions += plane_point - ref_center
    trans = center_in_plane(ag, plane, center_to=center_to, weights='mass', unwrap=True)(trans_u.trajectory.ts)
    assert_array_almost_equal(trans.positions, ref.positions, decimal=6)


def test_center_in_box_transformations_api(translate_universes):
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


def test_center_in_axis_transformations_api(translate_universes):
    # what happens when we center the coordinates arround the center of geometry of a residue?
    ref_u, trans_u = translate_universes
    ref = ref_u.trajectory.ts
    ref_center = np.float32([6, 7, 8])
    box_center = np.float32([186., 186.5, 187.])
    center_to = box_center
    axis = 'x'
    axis_point = np.float32([ref_center[0], center_to[1], center_to[2]])
    ref.positions += axis_point - ref_center
    ag = trans_u.residues[0].atoms
    trans_u.trajectory.add_transformations(center_in_axis(ag, axis))
    assert_array_almost_equal(trans_u.trajectory.ts.positions, ref.positions, decimal=6)


def test_center_in_plane_transformations_api(translate_universes):
    # what happens when we center the coordinates arround the center of geometry of a residue?
    ref_u, trans_u = translate_universes
    ref = ref_u.trajectory.ts
    ref_center = np.float32([6, 7, 8])
    box_center = np.float32([186., 186.5, 187.])
    plane = 'yz'
    plane_point = np.asarray([box_center[0], ref_center[1], ref_center[2]], np.float32)
    ref.positions += plane_point - ref_center
    ag = trans_u.residues[0].atoms
    trans_u.trajectory.add_transformations(center_in_plane(ag, plane))
    assert_array_almost_equal(trans_u.trajectory.ts.positions, ref.positions, decimal=6)
