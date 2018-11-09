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

import pytest
from six.moves import cPickle

import MDAnalysis as mda
from numpy.testing import (
    TestCase,
    assert_equal
)

import gc

from MDAnalysisTests.datafiles import PDB_small


class TestAtomGroupPickle(object):
    """Set up hopefully unique universes."""

    # _n marks named universes/atomgroups/pickled strings

    @staticmethod
    @pytest.fixture()
    def universe():
        return mda.Universe(PDB_small, PDB_small, PDB_small)

    @staticmethod
    @pytest.fixture()
    def universe_n():
        return mda.Universe(PDB_small, PDB_small, PDB_small, anchor_name="test1")

    @staticmethod
    @pytest.fixture()
    def ag(universe):
        return universe.atoms[:20]

    @staticmethod
    @pytest.fixture()
    def ag_n(universe_n):
        return universe_n.atoms[:10]

    @staticmethod
    @pytest.fixture()
    def pickle_str(ag):
        return cPickle.dumps(ag, protocol=cPickle.HIGHEST_PROTOCOL)

    @staticmethod
    @pytest.fixture()
    def pickle_str_n(ag_n):
        return cPickle.dumps(ag_n, protocol=cPickle.HIGHEST_PROTOCOL)

    def test_unpickle(self, pickle_str, ag, universe):
        """Test that an AtomGroup can be unpickled (Issue 293)"""
        newag = cPickle.loads(pickle_str)
        # Can unpickle
        assert_equal(ag.indices, newag.indices)
        assert newag.universe is universe, "Unpickled AtomGroup on wrong Universe."

    def test_unpickle_named(self, pickle_str_n, ag_n, universe_n):
        """Test that an AtomGroup can be unpickled (Issue 293)"""
        newag = cPickle.loads(pickle_str_n)
        # Can unpickle
        assert_equal(ag_n.indices, newag.indices)
        assert newag.universe is universe_n, "Unpickled AtomGroup on wrong Universe."

    def test_unpickle_missing(self):
        universe = mda.Universe(PDB_small, PDB_small, PDB_small)
        universe_n = mda.Universe(PDB_small, PDB_small, PDB_small,
                                  anchor_name="test1")
        ag = universe.atoms[:20]  # prototypical AtomGroup
        ag_n = universe_n.atoms[:10]
        pickle_str_n = cPickle.dumps(ag_n, protocol=cPickle.HIGHEST_PROTOCOL)
        # Kill AtomGroup and Universe
        del ag_n
        del universe_n
        # and make sure they're very dead
        gc.collect()
        # we shouldn't be able to unpickle
        # assert_raises(RuntimeError, cPickle.loads, pickle_str_n)
        with pytest.raises(RuntimeError):
            cPickle.loads(pickle_str_n)

    def test_unpickle_noanchor(self, universe, pickle_str):
        # Shouldn't unpickle if the universe is removed from the anchors
        universe.remove_anchor()
        # In the complex (parallel) testing environment there's the risk of
        # other compatible Universes being available for anchoring even after
        # this one is expressly removed.
        # assert_raises(RuntimeError, cPickle.loads, pickle_str)
        with pytest.raises(RuntimeError):
            cPickle.loads(pickle_str)
        # If this fails to raise an exception either:
        # 1-the anchoring Universe failed to remove_anchor or 2-another
        # Universe with the same characteristics was created for testing and is
        # being used as anchor."

    def test_unpickle_reanchor(self, universe, pickle_str, ag):
        # universe is removed from the anchors
        universe.remove_anchor()
        # now it goes back into the anchor list again
        universe.make_anchor()
        newag = cPickle.loads(pickle_str)
        assert_equal(ag.indices, newag.indices)
        assert newag.universe is universe, "Unpickled AtomGroup on wrong Universe."

    def test_unpickle_wrongname(self, universe_n, pickle_str_n):
        # we change the universe's anchor_name
        universe_n.anchor_name = "test2"
        # shouldn't unpickle if no name matches, even if there's a compatible
        # universe in the unnamed anchor list.
        with pytest.raises(RuntimeError):
            cPickle.loads(pickle_str_n)

    def test_unpickle_rename(self, universe_n, universe, pickle_str_n, ag_n):
        # we change universe_n's anchor_name
        universe_n.anchor_name = "test2"
        # and make universe a named anchor
        universe.anchor_name = "test1"
        newag = cPickle.loads(pickle_str_n)
        assert_equal(ag_n.indices, newag.indices)
        assert newag.universe is universe, "Unpickled AtomGroup on wrong Universe."

    def test_pickle_unpickle_empty(self, universe):
        """Test that an empty AtomGroup can be pickled/unpickled (Issue 293)"""
        ag = universe.atoms[[]]
        pickle_str = cPickle.dumps(ag, protocol=cPickle.HIGHEST_PROTOCOL)
        newag = cPickle.loads(pickle_str)
        assert len(newag) == 0


class TestPicklingUpdatingAtomGroups(object):

    @staticmethod
    @pytest.fixture()
    def u():
        return mda.Universe(PDB_small)

    def test_pickling_uag(self, u):
        ag = u.atoms[:100]
        uag = ag.select_atoms('name C', updating=True)
        pickle_str = cPickle.dumps(uag, protocol=cPickle.HIGHEST_PROTOCOL)
        new_uag = cPickle.loads(pickle_str)

        assert_equal(uag.indices, new_uag.indices)

    def test_pickling_uag_of_uag(self, u):
        uag1 = u.select_atoms('name C or name H', updating=True)
        uag2 = uag1.select_atoms('name C', updating=True)
        pickle_str = cPickle.dumps(uag2, protocol=cPickle.HIGHEST_PROTOCOL)
        new_uag2 = cPickle.loads(pickle_str)

        assert_equal(uag2.indices, new_uag2.indices)
