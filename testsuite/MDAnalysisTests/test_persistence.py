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
import cPickle
import warnings
import tempdir


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


class _GromacsReader_offsets(TestCase):
    # This base class assumes same lengths and dt for XTC and TRR test cases!
    filename = None
    ref_unitcell = np.array([80.017, 80.017, 80.017, 60., 60., 90.],
                            dtype=np.float32)
    # computed with Gromacs: 362.26999999999998 nm**3 * 1000 A**3/nm**3
    ref_volume = 362270.0
    ref_offsets = None
    _reader = None

    def setUp(self):
        # since offsets are automatically generated in the same directory
        # as the trajectory, we do everything from a temporary directory
        self.tmpdir = tempdir.TempDir()
        shutil.copy(self.filename, self.tmpdir.name)

        self.traj = os.path.join(self.tmpdir.name,
                                 os.path.basename(self.filename))

        self.trajectory = self._reader(self.traj)
        self.prec = 3
        self.ts = self.trajectory.ts

    def tearDown(self):
        del self.tmpdir
        del self.trajectory

    @dec.slow
    def test_offsets(self):
        self.trajectory._read_offsets(store=True)
        assert_array_almost_equal(self.trajectory._xdr.offsets,
                                  self.ref_offsets,
                                  err_msg="wrong frame offsets")

        outfile_offsets = XDR.offsets_filename(self.traj)
        with open(outfile_offsets) as f:
            saved_offsets = {k: v for k, v in np.load(f).iteritems()}

        assert_array_almost_equal(self.trajectory._xdr.offsets,
                                  saved_offsets['offsets'],
                                  err_msg="error saving frame offsets")
        assert_array_almost_equal(self.ref_offsets, saved_offsets['offsets'],
                                  err_msg="saved frame offsets don't match "
                                  "the known ones")

        self.trajectory._load_offsets()
        assert_array_almost_equal(self.trajectory._xdr.offsets,
                                  self.ref_offsets,
                                  err_msg="error loading frame offsets")
        assert_equal(saved_offsets['ctime'], os.path.getctime(self.traj))
        assert_equal(saved_offsets['size'], os.path.getsize(self.traj))

    @dec.slow
    def test_persistent_offsets_size_mismatch(self):
        self.trajectory._read_offsets(store=True)

        # check that stored offsets are not loaded when trajectory size differs
        # from stored size
        with open(XDR.offsets_filename(self.traj), 'rb') as f:
            saved_offsets = {k: v for k, v in np.load(f).iteritems()}
        saved_offsets['size'] += 1
        with open(XDR.offsets_filename(self.traj), 'wb') as f:
            np.savez(f, **saved_offsets)

        # u = MDAnalysis.Universe(self.top, self.traj)
        # assert_equal((u.trajectory._offsets is None), True)

    @dec.slow
    def test_persistent_offsets_last_frame_wrong(self):
        self.trajectory._read_offsets(store=True)

        # check that stored offsets are not loaded when the offsets themselves
        # appear to be wrong
        with open(XDR.offsets_filename(self.traj), 'rb') as f:
            saved_offsets = {k: v for k, v in np.load(f).iteritems()}
        saved_offsets['offsets'] += 1
        with open(XDR.offsets_filename(self.traj), 'wb') as f:
            np.savez(f, **saved_offsets)

        # with warnings.catch_warnings():
        #     u = MDAnalysis.Universe(self.top, self.traj)
        #     assert_equal((u.trajectory._xdr.offsets is None), True)

    @dec.slow
    def test_persistent_offsets_readonly(self):
        os.remove(XDR.offsets_filename(self.traj))
        assert_equal(os.path.exists(
            XDR.offsets_filename(self.trajectory.filename)), False)

        os.chmod(self.tmpdir.name, 0555)
        self.trajectory._read_offsets(store=True)
        assert_equal(os.path.exists(
            XDR.offsets_filename(self.trajectory.filename)), False)


class TestXTCReader_offsets(_GromacsReader_offsets):
    filename = XTC
    ref_offsets = np.array([0, 165188, 330364, 495520, 660708, 825872, 991044,
                            1156212, 1321384, 1486544])
    _reader = MDAnalysis.coordinates.XTC.XTCReader


class TestTRRReader_offsets(_GromacsReader_offsets):
    filename = TRR
    ref_offsets = np.array([0, 1144464, 2288928, 3433392, 4577856, 5722320,
                            6866784, 8011248, 9155712, 10300176])
    _reader = MDAnalysis.coordinates.TRR.TRRReader
