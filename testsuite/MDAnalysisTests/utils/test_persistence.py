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
from six.moves import cPickle

import MDAnalysis as mda
from MDAnalysis.core.groups import AtomGroup
from MDAnalysis.coordinates import XDR

import numpy as np
from numpy.testing import (assert_array_almost_equal, TestCase,
                           assert_array_equal, assert_, assert_raises,
                           assert_equal, dec)

import os
import gc
import shutil
import warnings

from MDAnalysisTests.datafiles import PDB_small, GRO, XTC, TRR
from MDAnalysisTests import make_Universe

class TestAtomGroupPickle(object):
    def setUp(self):
        """Set up hopefully unique universes."""
        # _n marks named universes/atomgroups/pickled strings
        self.universe = mda.Universe(PDB_small, PDB_small, PDB_small)
        self.universe_n = mda.Universe(PDB_small, PDB_small, PDB_small,
                                       anchor_name="test1")
        self.ag = self.universe.atoms[:20]  # prototypical AtomGroup
        self.ag_n = self.universe_n.atoms[:10]
        self.pickle_str = cPickle.dumps(self.ag,
                                        protocol=cPickle.HIGHEST_PROTOCOL)
        self.pickle_str_n = cPickle.dumps(self.ag_n,
                                          protocol=cPickle.HIGHEST_PROTOCOL)

    def tearDown(self):
        # Some tests might delete attributes.
        try:
            del self.universe
            del self.universe_n
            del self.ag
            del self.ag_n
        except AttributeError:
            pass

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
        # Kill AtomGroup and Universe
        del self.ag_n
        del self.universe_n
        # and make sure they're very dead
        gc.collect()
        # we shouldn't be able to unpickle
        assert_raises(RuntimeError, cPickle.loads, self.pickle_str_n)

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

    def test_pickle_unpickle_empty(self):
        """Test that an empty AtomGroup can be pickled/unpickled (Issue 293)"""
        ag = self.universe.atoms[[]]
        pickle_str = cPickle.dumps(ag, protocol=cPickle.HIGHEST_PROTOCOL)
        newag = cPickle.loads(pickle_str)
        assert_equal(len(newag), 0)


class TestPicklingUpdatingAtomGroups(object):
    def setUp(self):
        self.u = mda.Universe(PDB_small)

    def tearDown(self):
        del self.u

    def test_pickling_uag(self):
        ag = self.u.atoms[:100]
        uag = ag.select_atoms('name C', updating=True)
        pickle_str = cPickle.dumps(uag, protocol=cPickle.HIGHEST_PROTOCOL)
        new_uag = cPickle.loads(pickle_str)

        assert_array_equal(uag.indices, new_uag.indices)

    def test_pickling_uag_of_uag(self):
        uag1 = self.u.select_atoms('name C or name H', updating=True)
        uag2 = uag1.select_atoms('name C', updating=True)
        pickle_str = cPickle.dumps(uag2, protocol=cPickle.HIGHEST_PROTOCOL)
        new_uag2 = cPickle.loads(pickle_str)

        assert_array_equal(uag2.indices, new_uag2.indices)
