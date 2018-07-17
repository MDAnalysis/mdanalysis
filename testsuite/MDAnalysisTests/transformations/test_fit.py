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

from MDAnalysisTests import make_Universe
from MDAnalysis.transformations.fit import alignto, fit_translation


@pytest.fixture()
def test_universe():
    # make a test universe
    test = make_Universe(('masses', ), trajectory=True)
    ref = make_Universe(('masses', ), trajectory=True)
    ref.trajectory.ts.positions += np.asarray([10, 10, 10], np.float32)
    return test, ref


def test_alignto_coords(test_universe):
    # when aligning the two universes do their coordinates become similar?
    test_u = test_universe[0]
    ref_u = test_universe[1]
    alignto(test_u, ref_u, weights="mass")(test_u.trajectory.ts)
    assert_array_almost_equal(test_u.trajectory.ts.positions, ref_u.trajectory.ts.positions, decimal=6)
    

def test_alignto_transformations_api(test_universe):
    test_u = test_universe[0]
    ref_u = test_universe[1]
    transform = alignto(test_u, ref_u, weights="mass")
    test_u.trajectory.add_transformations(transform)
    assert_array_almost_equal(test_u.trajectory.ts.positions, ref_u.trajectory.ts.positions, decimal=6)


def test_fit_translation_bad_ag(test_universe):
    ts = test_universe[0].trajectory.ts
    test_u = test_universe[0]
    ref_u  = test_universe[1]
    # what happens if something other than an AtomGroup or Universe is given?
    bad_ag = 1
    with pytest.raises(ValueError): 
        fit_translation(bad_ag, ref_u)(ts)


def test_fit_translation_bad_center(test_universe):
    ts = test_universe[0].trajectory.ts
    test_u = test_universe[0]
    ref_u  = test_universe[1]
    # what happens if a bad string for center is given?
    bad_center = "123"
    with pytest.raises(ValueError): 
        fit_translation(test_u, ref_u, center_of = bad_center)(ts)


def test_fit_translation_bad_plane(test_universe):
    ts = test_universe[0].trajectory.ts
    test_u = test_universe[0]
    ref_u  = test_universe[1]
    # what happens if a bad string for center is given?
    bad_plane = "123"
    with pytest.raises(ValueError): 
        fit_translation(test_u, ref_u, plane=bad_plane)(ts)


def test_fit_translation_no_masses(test_universe):
    ts = test_universe[0].trajectory.ts
    test_u = test_universe[0]
    # create a universe without masses
    ref_u  = make_Universe()
    # what happens Universe without masses is given?
    with pytest.raises(AttributeError): 
        fit_translation(test_u, ref_u, center_of="mass")(ts)


def test_fit_translation_no_options(test_universe):
    test_u = test_universe[0]
    ref_u = test_universe[1]
    fit_translation(test_u, ref_u)(test_u.trajectory.ts)
    # what happens when no options are passed?
    assert_array_almost_equal(test_u.trajectory.ts.positions, ref_u.trajectory.ts.positions, decimal=6)


def test_fit_translation_com(test_universe):
    test_u = test_universe[0]
    ref_u = test_universe[1]
    fit_translation(test_u, ref_u, center_of="mass")(test_u.trajectory.ts)
    # what happens when the center o mass is used?
    assert_array_almost_equal(test_u.trajectory.ts.positions, ref_u.trajectory.ts.positions, decimal=6)


def test_fit_translation_plane(test_universe):
    test_u = test_universe[0]
    ref_u = test_universe[1]
    # translate the test universe on the x and y coordinates only
    fit_translation(test_u, ref_u, plane="xy")(test_u.trajectory.ts)
    # the reference is 10 angstrom in the z coordinate above the test universe
    shiftz = np.asanyarray([0, 0, -10], np.float32)
    ref_coordinates = ref_u.trajectory.ts.positions + shiftz 
    assert_array_almost_equal(test_u.trajectory.ts.positions, ref_coordinates, decimal=6)


def test_fit_translation_all_options(test_universe):
    test_u = test_universe[0]
    ref_u = test_universe[1]
    # translate the test universe on the x and y coordinates only
    fit_translation(test_u, ref_u, plane="xy", center_of="mass")(test_u.trajectory.ts)
    # the reference is 10 angstrom in the z coordinate above the test universe
    shiftz = np.asanyarray([0, 0, -10], np.float32)
    ref_coordinates = ref_u.trajectory.ts.positions + shiftz 
    assert_array_almost_equal(test_u.trajectory.ts.positions, ref_coordinates, decimal=6)


def test_fit_translation_transformations_api(test_universe):
    test_u = test_universe[0]
    ref_u = test_universe[1]
    transform = fit_translation(test_u, ref_u)
    test_u.trajectory.add_transformations(transform)
    assert_array_almost_equal(test_u.trajectory.ts.positions, ref_u.trajectory.ts.positions, decimal=6)
