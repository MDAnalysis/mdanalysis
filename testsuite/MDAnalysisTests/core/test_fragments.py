# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.mdanalysis.org
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
    assert_,
    assert_raises,
    assert_array_equal,
)

from MDAnalysis.core.topologyattrs import Bonds
from MDAnalysis.core import groups
from MDAnalysis import NoDataError

from MDAnalysisTests import make_Universe

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
    @staticmethod
    def make_case1():
        return make_starshape()

    @staticmethod
    def make_case2():
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

    @staticmethod
    def _check_total_frags(u):
        # should be 5 fragments of 25 atoms
        assert_(len(u.atoms.fragments) == 5)
        for i, frag in enumerate(u.atoms.fragments):
            assert_(len(frag) == 25)

    @staticmethod
    def _check_frag_external_ordering(u):
        # check fragments are sorted
        for i, frag in enumerate(u.atoms.fragments):
            assert_(frag[0].index == i * 25)

    @staticmethod
    def _check_frag_internal_ordering(u):
        # check atoms are sorted within fragments
        for i, frag in enumerate(u.atoms.fragments):
            assert_array_equal(frag.ix, np.arange(25) + i * 25)

    @staticmethod
    def _check_atom_access(u):
        # check atom can access fragment
        for at in (u.atoms[0], u.atoms[76], u.atoms[111]):
            frag = at.fragment
            assert_(isinstance(frag, groups.AtomGroup))
            assert_(len(frag) == 25)
            assert_(at in frag)

    @staticmethod
    def _check_atomgroup_access(u):
        # check atomgroup can access fragments
        # first 60 atoms have 3 fragments, given as tuple
        # each fragment should still be 25 atoms
        frags = u.atoms[:60].fragments
        assert_(len(frags) == 3)
        assert_(isinstance(frags, tuple))
        for frag in frags:
            assert_(len(frag) == 25)

    def test_fragments(self):
        for case in (self.make_case1,self.make_case2):
            u = case()
            yield self._check_total_frags, u
            yield self._check_frag_internal_ordering, u
            yield self._check_frag_external_ordering, u
            yield self._check_atom_access, u
            yield self._check_atomgroup_access, u

    def test_atomgroup_fragments_nobonds_NDE(self):
        # should raise NDE
        u = make_Universe()

        assert_raises(NoDataError, getattr, u.atoms[:10], 'fragments')

    def test_atom_fragment_nobonds_NDE(self):
        # should raise NDE
        u = make_Universe()

        assert_raises(NoDataError, getattr, u.atoms[10], 'fragment')
