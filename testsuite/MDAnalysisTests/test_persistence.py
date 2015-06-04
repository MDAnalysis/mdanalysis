# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

import MDAnalysis
from MDAnalysis.tests.datafiles import PSF, DCD, PDB_small
import MDAnalysis.core.AtomGroup
from MDAnalysis.core.AtomGroup import AtomGroup

import numpy
from numpy.testing import *

import gc
import cPickle

class TestAtomGroupPickle(TestCase):

    def setUp(self):
        """Set up a hopefully unique universe."""
        self.universe = MDAnalysis.Universe(PDB_small, PDB_small, PDB_small)
        self.ag = self.universe.atoms[:20]  # prototypical AtomGroup

    def tearDown(self):
        del self.universe
        del self.ag

    def test_pickle_unpickle(self):
        """Test that an AtomGroup can be pickled/unpickled (Issue 293)"""
        pickle_str = cPickle.dumps(self.ag, protocol=cPickle.HIGHEST_PROTOCOL)
        newag = cPickle.loads(pickle_str)
        # Can unpickle
        assert_array_equal(self.ag.indices(), newag.indices())
        # Shouldn't unpickle if the universe is removed from the anchors
        self.universe.remove_anchor()
        # In the complex (parallel) testing environment there's the risk of other compatible Universes being available
        #  for anchoring even after this one is expressly removed.
        assert_raises(RuntimeError, cPickle.loads, pickle_str)
          # If this fails to raise an exception either:"
          #  1-the anchoring Universe failed to remove_anchor or"
          #  2-another Universe with the same characteristics was created for testing and is being used as anchor."
        # Now it goes back into the anchor list again
        self.universe.make_anchor()
        newag2 = cPickle.loads(pickle_str)
        assert_array_equal(self.ag.indices(), newag2.indices())

def test_pickle_unpickle_empty():
    """Test that an empty AtomGroup can be pickled/unpickled (Issue 293)"""
    import gc
    ag = AtomGroup([])
    pickle_str = cPickle.dumps(ag, protocol=cPickle.HIGHEST_PROTOCOL)
    newag = cPickle.loads(pickle_str)
    assert_equal(len(newag), 0)

#def test_memleak(self):
#    u1 = MDAnalysis.Universe(PSF, DCD)
#    u2 = MDAnalysis.Universe(PSF, DCD, DCD)
#    del u1
#    del u2
#    gc.collect()
#    assert_equal(len(gc.garbage), 0, "Garbage collector can't collect the following: %r" % gc.garbage)
#
