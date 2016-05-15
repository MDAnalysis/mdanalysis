# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
# Beckstein and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from six.moves import cPickle

import MDAnalysis
from MDAnalysis.tests.datafiles import PDB_small, GRO, XTC, TRR
import MDAnalysis.core.AtomGroup
from MDAnalysis.core.AtomGroup import AtomGroup
from MDAnalysis.coordinates import XDR

import numpy as np
from numpy.testing import (assert_array_almost_equal, TestCase,
                           assert_array_equal, assert_, assert_raises,
                           assert_equal, dec)

import os
import gc
import shutil
import warnings


class TestAtomGroupPickle(TestCase):
    def setUp(self):
        """Set up hopefully unique universes."""
        # _n marks named universes/atomgroups/pickled strings
        self.universe = MDAnalysis.Universe(PDB_small, PDB_small, PDB_small)
        self.universe_n = MDAnalysis.Universe(PDB_small, PDB_small, PDB_small,
                                              is_anchor=False,
                                              anchor_name="test1")
        self.ag = self.universe.atoms[:20]  # prototypical AtomGroup
        self.ag_n = self.universe_n.atoms[:10]
        self.pickle_str = cPickle.dumps(self.ag,
                                        protocol=cPickle.HIGHEST_PROTOCOL)
        self.pickle_str_n = cPickle.dumps(self.ag_n,
                                          protocol=cPickle.HIGHEST_PROTOCOL)

    def tearDown(self):
        del self.universe
        del self.universe_n
        del self.ag
        del self.ag_n

    def test_unpickle(self):
        """Test that an AtomGroup can be unpickled (Issue 293)"""
        newag = cPickle.loads(self.pickle_str)
        # Can unpickle
        assert_array_equal(self.ag.indices, newag.indices)
        assert_(newag.universe is self.universe,
                "Unpickled AtomGroup on wrong Universe.")

    def test_unpickle_named(self):
        """Test that an AtomGroup can be unpickled (Issue 293)"""
        newag = cPickle.loads(self.pickle_str_n)
        # Can unpickle
        assert_array_equal(self.ag_n.indices, newag.indices)
        assert_(newag.universe is self.universe_n,
                "Unpickled AtomGroup on wrong Universe.")

    def test_unpickle_missing(self):
        # we kill the universes and ag's
        self.tearDown()
        # and make sure they're very dead
        gc.collect()
        # we shouldn't be able to unpickle
        assert_raises(RuntimeError, cPickle.loads, self.pickle_str)
        assert_raises(RuntimeError, cPickle.loads, self.pickle_str_n)
        # must reset these, otherwise tearDown fails
        self.setUp()

    def test_unpickle_noanchor(self):
        # Shouldn't unpickle if the universe is removed from the anchors
        self.universe.remove_anchor()
        # In the complex (parallel) testing environment there's the risk of
        # other compatible Universes being available for anchoring even after
        # this one is expressly removed.
        assert_raises(RuntimeError, cPickle.loads, self.pickle_str)
        # If this fails to raise an exception either:
        # 1-the anchoring Universe failed to remove_anchor or 2-another
        # Universe with the same characteristics was created for testing and is
        # being used as anchor."

    def test_unpickle_reanchor(self):
        # universe is removed from the anchors
        self.universe.remove_anchor()
        # now it goes back into the anchor list again
        self.universe.make_anchor()
        newag = cPickle.loads(self.pickle_str)
        assert_array_equal(self.ag.indices, newag.indices)
        assert_(newag.universe is self.universe,
                "Unpickled AtomGroup on wrong Universe.")

    def test_unpickle_reanchor_other(self):
        # universe is removed from the anchors
        self.universe.remove_anchor()
        # and universe_n goes into the anchor list
        self.universe_n.make_anchor()
        newag = cPickle.loads(self.pickle_str)
        assert_array_equal(self.ag.indices, newag.indices)
        assert_(newag.universe is self.universe_n,
                "Unpickled AtomGroup on wrong Universe.")

    def test_unpickle_wrongname(self):
        # we change the universe's anchor_name
        self.universe_n.anchor_name = "test2"
        # shouldn't unpickle if no name matches, even if there's a compatible
        # universe in the unnamed anchor list.
        assert_raises(RuntimeError, cPickle.loads, self.pickle_str_n)

    def test_unpickle_rename(self):
        # we change universe_n's anchor_name
        self.universe_n.anchor_name = "test2"
        # and make universe a named anchor
        self.universe.anchor_name = "test1"
        newag = cPickle.loads(self.pickle_str_n)
        assert_array_equal(self.ag_n.indices, newag.indices)
        assert_(newag.universe is self.universe,
                "Unpickled AtomGroup on wrong Universe.")


class TestEmptyAtomGroupPickle(TestCase):
    # This comes in a class just to get memleak testing
    def test_pickle_unpickle_empty(self):
        """Test that an empty AtomGroup can be pickled/unpickled (Issue 293)"""
        ag = AtomGroup([])
        pickle_str = cPickle.dumps(ag, protocol=cPickle.HIGHEST_PROTOCOL)
        newag = cPickle.loads(pickle_str)
        assert_equal(len(newag), 0)
