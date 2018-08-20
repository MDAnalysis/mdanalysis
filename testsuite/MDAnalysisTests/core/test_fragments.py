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
from six.moves import range

import numpy as np
from numpy.testing import (
    assert_equal,
)
import pytest

import MDAnalysis as mda
from MDAnalysis.core.topologyattrs import Bonds
from MDAnalysis.core import groups
from MDAnalysis import NoDataError

from MDAnalysisTests import make_Universe
from MDAnalysisTests.datafiles import TPR, XTC


# Also used in topology/test_guessers
def make_starshape():
    u = make_Universe()
    bonds = []
    for seg in range(5):
        segbase = seg * 25
        for res in range(5):
            # offset for atoms in this res
            base = segbase + 5 * res
            bonds.append((0 + base, 1 + base))
            bonds.append((1 + base, 2 + base))
            bonds.append((1 + base, 3 + base))
            bonds.append((1 + base, 4 + base))
            if not res == 4:  # last res doesn't link onwards
                bonds.append((4 + base, 5 + base))
    u.add_TopologyAttr(Bonds(bonds))
    return u


def case1():
    return make_starshape()


def case2():
    u = make_Universe()
    bonds = []
    for seg in range(5):
        segbase = seg * 25
        for res in range(5):
            # offset for atoms in this res
            base = segbase + 5 * res
            bonds.append((0 + base, 1 + base))
            bonds.append((1 + base, 2 + base))
            bonds.append((2 + base, 3 + base))
            bonds.append((3 + base, 4 + base))
            bonds.append((1 + base, 4 + base))
            if not res == 4:  # last res doesn't link onwards
                bonds.append((0 + base, 5 + base))
    u.add_TopologyAttr(Bonds(bonds))
    return u


class TestFragments(object):
    """Use 125 atom test Universe

    5 segments of 5 residues of 5 atoms


    Case1
    -----

    Star shapes to try and test the branching prediction

      o   |   o   |   o
      |   |   |   |   |
    o-o-o-|-o-o-o-|-o-o-o
      |   |   |   |   |
      o   |   o   |x3 o


    Case2
    -----

    4-ring pendants to test cyclic conditions

      o------o------o
      |      |      |
      o      o      o
     / \    / \    / \
    o   o  o   o  o   o
     \ /    \ /    \ /
      o      o      o
    Test ring molecules?
    """

    @pytest.mark.parametrize('u', (
            case1(),
            case2()
    ))
    def test_total_frags(self, u):
        # should be 5 fragments of 25 atoms
        assert len(u.atoms.fragments) == 5
        for i, frag in enumerate(u.atoms.fragments):
            assert len(frag) == 25

    @pytest.mark.parametrize('u', (
            case1(),
            case2()
    ))
    def test_frag_external_ordering(self, u):
        # check fragments are sorted
        for i, frag in enumerate(u.atoms.fragments):
            assert frag[0].index == i * 25

    @pytest.mark.parametrize('u', (
            case1(),
            case2()
    ))
    def test_frag_internal_ordering(self, u):
        # check atoms are sorted within fragments
        for i, frag in enumerate(u.atoms.fragments):
            assert_equal(frag.ix, np.arange(25) + i * 25)

    @pytest.mark.parametrize('u', (
            case1(),
            case2()
    ))
    def test_atom_access(self, u):
        # check atom can access fragment
        for at in (u.atoms[0], u.atoms[76], u.atoms[111]):
            frag = at.fragment
            assert isinstance(frag, groups.AtomGroup)
            assert len(frag) == 25
            assert at in frag

    @pytest.mark.parametrize('u', (
            case1(),
            case2()
    ))
    def test_atomgroup_access(self, u):
        # check atomgroup can access fragments
        # first 60 atoms have 3 fragments, given as tuple
        # each fragment should still be 25 atoms
        frags = u.atoms[:60].fragments
        assert len(frags) == 3
        assert isinstance(frags, tuple)
        for frag in frags:
            assert len(frag) == 25

    def test_atomgroup_fragments_nobonds_NDE(self):
        # should raise NDE
        u = make_Universe()
        with pytest.raises(NoDataError):
            getattr(u.atoms[:10], 'fragments')

    def test_atom_fragment_nobonds_NDE(self):
        # should raise NDE
        u = make_Universe()
        with pytest.raises(NoDataError):
            getattr(u.atoms[10], 'fragment')


def test_tpr_fragments():
    frags = mda.Universe(TPR, XTC).atoms.fragments

    assert len(frags[0]) == 3341
