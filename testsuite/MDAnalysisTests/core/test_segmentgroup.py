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
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import absolute_import

from numpy.testing import (
    assert_equal,
)
import pytest

import MDAnalysis as mda

from MDAnalysisTests.datafiles import PSF, DCD


# Legacy tests from before 363
@pytest.fixture()
def universe():
    return mda.Universe(PSF, DCD)


@pytest.fixture()
def g(universe):
    return universe.atoms.segments

def test_newSegmentGroup(universe):
    """test that slicing a SegmentGroup returns a new SegmentGroup (Issue 135)"""
    g = universe.atoms.segments
    newg = g[:]
    assert isinstance(newg, mda.core.groups.SegmentGroup)
    assert_equal(len(newg), len(g))


def test_n_atoms(g):
    assert_equal(g.n_atoms, 3341)


def test_n_residues(g):
    assert_equal(g.n_residues, 214)


def test_resids_dim(g):
    assert_equal(len(g.resids), len(g))
    for seg, resids in zip(g, g.resids):
        assert len(resids) == len(seg.residues)
        assert_equal(seg.residues.resids, resids)


def test_resnums_dim(g):
    assert_equal(len(g.resnums), len(g))
    for seg, resnums in zip(g, g.resnums):
        assert len(resnums) == len(seg.residues)
        assert_equal(seg.residues.resnums, resnums)


def test_segids_dim(g):
    assert_equal(len(g.segids), len(g))


def test_set_segids(universe):
    s = universe.select_atoms('all').segments
    s.segids = 'ADK'
    assert_equal(universe.segments.segids, ['ADK'],
                 err_msg="failed to set_segid on segments")


def test_set_segid_updates_(universe):
    g = universe.select_atoms("resid 10:18").segments
    g.segids = 'ADK'
    assert_equal(g.segids, ['ADK'],
                 err_msg="old selection was not changed in place after set_segid")


def test_atom_order(universe):
    assert_equal(universe.segments.atoms.indices,
                 sorted(universe.segments.atoms.indices))
