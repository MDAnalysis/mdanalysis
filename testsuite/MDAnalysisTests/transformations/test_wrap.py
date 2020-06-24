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

import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal

import MDAnalysis as mda
from MDAnalysis.transformations import wrap, unwrap
from MDAnalysis.lib._cutil import make_whole
from MDAnalysisTests.datafiles import fullerene
from MDAnalysisTests.datafiles import GRO, TPR


@pytest.fixture()
def wrap_universes():
    # create the Universe objects for the tests
    # this universe is used for the unwrap testing cases
    reference = mda.Universe(fullerene)
    transformed = mda.Universe(fullerene)
    transformed.dimensions = np.asarray([10, 10, 10, 90, 90, 90], np.float32)
    transformed.atoms.wrap()
    
    return reference, transformed


@pytest.fixture()
def compound_wrap_universes():
    # create universe for the compound tests
    # this universe has multiple residues, segments and
    # fragments
    transformed = mda.Universe(TPR, GRO)
    reference = mda.Universe(TPR, GRO)
    # wrap the atoms back into the unit cell
    # in this coordinate file only the protein
    # is broken across PBC however the system 
    # shape is not the same as the unit cell
    make_whole(reference.select_atoms("protein"))
    make_whole(transformed.select_atoms("protein"))
    
    return transformed, reference


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
def test_wrap_bad_ag(wrap_universes, ag):
    # this universe has a box size zero
    ts = wrap_universes[0].trajectory.ts
    # what happens if something other than an AtomGroup is given?
    bad_ag = ag
    with pytest.raises(AttributeError): 
        wrap(bad_ag)(ts)


def test_wrap_no_options(wrap_universes):
    # since were testing if the wrapping works
    # the reference and the transformed are switched
    trans, ref = wrap_universes
    trans.dimensions = ref.dimensions
    wrap(trans.atoms)(trans.trajectory.ts)
    assert_array_almost_equal(trans.trajectory.ts.positions, ref.trajectory.ts.positions, decimal=6)


@pytest.mark.parametrize('compound', (
    "group",
    "residues",
    "segments",
    "fragments")
)
def test_wrap_with_compounds(compound_wrap_universes, compound):
    trans, ref= compound_wrap_universes
    ref.select_atoms("not resname SOL").wrap(compound=compound)
    wrap(trans.select_atoms("not resname SOL"), compound=compound)(trans.trajectory.ts)
    assert_array_almost_equal(trans.trajectory.ts.positions, ref.trajectory.ts.positions, decimal=6)


def test_wrap_api(wrap_universes):
    trans, ref = wrap_universes
    trans.dimensions = ref.dimensions
    trans.trajectory.add_transformations(wrap(trans.atoms))
    assert_array_almost_equal(trans.trajectory.ts.positions, ref.trajectory.ts.positions, decimal=6)


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
def test_unwrap_bad_ag(wrap_universes, ag):
    # this universe has a box size zero
    ts = wrap_universes[0].trajectory.ts
    # what happens if something other than an AtomGroup is given?
    bad_ag = ag
    with pytest.raises(AttributeError): 
        unwrap(bad_ag)(ts)


def test_unwrap(wrap_universes):
    ref, trans = wrap_universes
    # after rebuild the trans molecule it should match the reference
    unwrap(trans.atoms)(trans.trajectory.ts)
    assert_array_almost_equal(trans.trajectory.ts.positions, ref.trajectory.ts.positions, decimal=6)


def test_unwrap_api(wrap_universes):
    ref, trans = wrap_universes
    # after rebuild the trans molecule it should match the reference
    trans.trajectory.add_transformations(unwrap(trans.atoms))
    assert_array_almost_equal(trans.trajectory.ts.positions, ref.trajectory.ts.positions, decimal=6)
