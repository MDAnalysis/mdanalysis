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
import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal

from MDAnalysisTests import make_Universe
from MDAnalysis.transformations.fit import fit_translation, fit_rot_trans
from MDAnalysis.lib.transformations import rotation_matrix


@pytest.fixture()
def fit_universe():
    # make a test universe
    test = make_Universe(('masses', ), trajectory=True)
    ref = make_Universe(('masses', ), trajectory=True)
    ref.atoms.positions += np.asarray([10, 10, 10], np.float32)
    return test, ref


@pytest.mark.parametrize('universe', (
    [0, 1],
    [0, 1, 2, 3, 4],
    np.array([0, 1]),
    np.array([0, 1, 2, 3, 4]),
    np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]]),
    np.array([[0], [1], [2]]),
    'thisisnotanag',
    1)
)
def test_fit_translation_bad_ag(fit_universe, universe):
    ts = fit_universe[0].trajectory.ts
    test_u = fit_universe[0]
    ref_u  = fit_universe[1]
    # what happens if something other than an AtomGroup or Universe is given?
    with pytest.raises(AttributeError):
        fit_translation(universe, ref_u)(ts)


@pytest.mark.parametrize('weights', (
    " ",
    "totallynotmasses",
    123456789,
    [0, 1, 2, 3, 4],
    np.array([0, 1]),
    np.array([0, 1, 2, 3, 4]),
    np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]]))
)
def test_fit_translation_bad_weights(fit_universe, weights):
    ts = fit_universe[0].trajectory.ts
    test_u = fit_universe[0]
    ref_u  = fit_universe[1]
    # what happens if a bad string for center is given?
    with pytest.raises(ValueError):
        fit_translation(test_u, ref_u, weights=weights)(ts)


@pytest.mark.parametrize('plane', (
    1,
    [0, 1],
    [0, 1, 2, 3, 4],
    np.array([0, 1]),
    "xyz",
    "notaplane")
)
def test_fit_translation_bad_plane(fit_universe, plane):
    ts = fit_universe[0].trajectory.ts
    test_u = fit_universe[0]
    ref_u  = fit_universe[1]
    # what happens if a bad string for center is given?
    with pytest.raises(ValueError):
        fit_translation(test_u, ref_u, plane=plane)(ts)


def test_fit_translation_no_masses(fit_universe):
    ts = fit_universe[0].trajectory.ts
    test_u = fit_universe[0]
    # create a universe without masses
    ref_u  = make_Universe()
    # what happens Universe without masses is given?
    with pytest.raises(TypeError) as exc:
        fit_translation(test_u, ref_u, weights="mass")(ts)
    assert 'atoms.masses is missing' in str(exc.value)


def test_fit_translation_no_options(fit_universe):
    test_u = fit_universe[0]
    ref_u = fit_universe[1]
    fit_translation(test_u, ref_u)(test_u.trajectory.ts)
    # what happens when no options are passed?
    assert_array_almost_equal(test_u.trajectory.ts.positions, ref_u.trajectory.ts.positions, decimal=6)

def test_fit_translation_residue_mismatch(fit_universe):
    test_u = fit_universe[0]
    ref_u = fit_universe[1].residues[:-1].atoms
    with pytest.raises(ValueError, match='number of residues'):
        fit_translation(test_u, ref_u)(test_u.trajectory.ts)

def test_fit_translation_com(fit_universe):
    test_u = fit_universe[0]
    ref_u = fit_universe[1]
    fit_translation(test_u, ref_u, weights="mass")(test_u.trajectory.ts)
    # what happens when the center o mass is used?
    assert_array_almost_equal(test_u.trajectory.ts.positions, ref_u.trajectory.ts.positions, decimal=6)


@pytest.mark.parametrize('plane', (
    "yz",
    "xz",
    "xy")
)
def test_fit_translation_plane(fit_universe, plane):
    test_u = fit_universe[0]
    ref_u = fit_universe[1]
    axes = {'yz' : 0, 'xz' : 1, 'xy' : 2}
    idx = axes[plane]
    # translate the test universe on the plane coordinates only
    fit_translation(test_u, ref_u, plane=plane)(test_u.trajectory.ts)
    # the reference is 10 angstrom in all coordinates above the test universe
    shiftz = np.asanyarray([0, 0, 0], np.float32)
    shiftz[idx] = -10
    ref_coordinates = ref_u.trajectory.ts.positions + shiftz
    assert_array_almost_equal(test_u.trajectory.ts.positions, ref_coordinates, decimal=6)


def test_fit_translation_all_options(fit_universe):
    test_u = fit_universe[0]
    ref_u = fit_universe[1]
    # translate the test universe on the x and y coordinates only
    fit_translation(test_u, ref_u, plane="xy", weights="mass")(test_u.trajectory.ts)
    # the reference is 10 angstrom in the z coordinate above the test universe
    shiftz = np.asanyarray([0, 0, -10], np.float32)
    ref_coordinates = ref_u.trajectory.ts.positions + shiftz
    assert_array_almost_equal(test_u.trajectory.ts.positions, ref_coordinates, decimal=6)


def test_fit_translation_transformations_api(fit_universe):
    test_u = fit_universe[0]
    ref_u = fit_universe[1]
    transform = fit_translation(test_u, ref_u)
    test_u.trajectory.add_transformations(transform)
    assert_array_almost_equal(test_u.trajectory.ts.positions, ref_u.trajectory.ts.positions, decimal=6)


@pytest.mark.parametrize('universe', (
    [0, 1],
    [0, 1, 2, 3, 4],
    np.array([0, 1]),
    np.array([0, 1, 2, 3, 4]),
    np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]]),
    np.array([[0], [1], [2]]),
    'thisisnotanag',
    1)
)
def test_fit_rot_trans_bad_universe(fit_universe, universe):
    test_u = fit_universe[0]
    ref_u= universe
    with pytest.raises(AttributeError):
        fit_rot_trans(test_u, ref_u)(test_u.trajectory.ts)


def test_fit_rot_trans_shorter_universe(fit_universe):
    ref_u = fit_universe[1]
    bad_u =fit_universe[0].atoms[0:5]
    test_u= bad_u
    with pytest.raises(ValueError):
        fit_rot_trans(test_u, ref_u)(test_u.trajectory.ts)


@pytest.mark.parametrize('weights', (
    " ",
    "totallynotmasses",
    123456789,
    [0, 1, 2, 3, 4],
    np.array([0, 1]),
    np.array([0, 1, 2, 3, 4]),
    np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]]))
)
def test_fit_rot_trans_bad_weights(fit_universe, weights):
    test_u = fit_universe[0]
    ref_u = fit_universe[1]
    bad_weights = weights
    with pytest.raises(ValueError):
        fit_rot_trans(test_u, ref_u, weights=bad_weights)(test_u.trajectory.ts)


@pytest.mark.parametrize('plane', (
    " ",
    "totallynotaplane",
    "xyz",
    123456789,
    [0, 1, 2, 3, 4],
    np.array([0, 1]),
    np.array([0, 1, 2, 3, 4]),
    np.array([[0, 1, 2], [3, 4, 5], [6, 7, 8]]))
)
def test_fit_rot_trans_bad_plane(fit_universe, plane):
    test_u = fit_universe[0]
    ref_u = fit_universe[1]
    with pytest.raises(ValueError):
        fit_rot_trans(test_u, ref_u, plane=plane)(test_u.trajectory.ts)


def test_fit_rot_trans_no_options(fit_universe):
    test_u = fit_universe[0]
    ref_u = fit_universe[1]
    ref_com = ref_u.atoms.center(None)
    ref_u.trajectory.ts.positions -= ref_com
    R = rotation_matrix(np.pi/3, [1,0,0])[:3,:3]
    ref_u.trajectory.ts.positions = np.dot(ref_u.trajectory.ts.positions, R)
    ref_u.trajectory.ts.positions += ref_com
    fit_rot_trans(test_u, ref_u)(test_u.trajectory.ts)
    assert_array_almost_equal(test_u.trajectory.ts.positions, ref_u.trajectory.ts.positions, decimal=3)


@pytest.mark.parametrize('plane', (
    "yz",
    "xz",
    "xy")
)
def test_fit_rot_trans_plane(fit_universe, plane):
    # the reference is rotated in the x axis so removing the translations and rotations
    # in the yz plane should return the same as the fitting without specifying a plane
    test_u = fit_universe[0]
    ref_u = fit_universe[1]
    ref_com = ref_u.atoms.center(None)
    mobile_com = test_u.atoms.center(None)
    axes = {'yz' : 0, 'xz' : 1, 'xy' : 2}
    idx = axes[plane]
    rotaxis = np.asarray([0,0,0])
    rotaxis[idx]=1
    ref_u.trajectory.ts.positions -= ref_com
    R = rotation_matrix(np.pi/3, rotaxis)[:3,:3]
    ref_u.trajectory.ts.positions = np.dot(ref_u.trajectory.ts.positions, R)
    ref_com[idx] = mobile_com[idx]
    ref_u.trajectory.ts.positions += ref_com
    fit_rot_trans(test_u, ref_u, plane=plane)(test_u.trajectory.ts)
    assert_array_almost_equal(test_u.trajectory.ts.positions[:,idx], ref_u.trajectory.ts.positions[:,idx], decimal=3)


def test_fit_rot_trans_transformations_api(fit_universe):
    test_u = fit_universe[0]
    ref_u = fit_universe[1]
    ref_com = ref_u.atoms.center(None)
    ref_u.trajectory.ts.positions -= ref_com
    R = rotation_matrix(np.pi/3, [1,0,0])[:3,:3]
    ref_u.trajectory.ts.positions = np.dot(ref_u.trajectory.ts.positions, R)
    ref_u.trajectory.ts.positions += ref_com
    transform = fit_rot_trans(test_u, ref_u)
    test_u.trajectory.add_transformations(transform)
    assert_array_almost_equal(test_u.trajectory.ts.positions, ref_u.trajectory.ts.positions, decimal=3)
