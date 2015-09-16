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
from MDAnalysis.tests.datafiles import PSF, DCD, PDB_small, GRO, XTC, TRR
import MDAnalysis.core.AtomGroup
from MDAnalysis.core.AtomGroup import AtomGroup

import numpy as np
from numpy.testing import *

import os
import shutil
import cPickle
import tempfile
import warnings

class TestAtomGroupPickle(TestCase):
    def setUp(self):
        """Set up hopefully unique universes."""
        # _n marks named universes/atomgroups/pickled strings
        self.universe = MDAnalysis.Universe(PDB_small, PDB_small, PDB_small)
        self.universe_n = MDAnalysis.Universe(PDB_small, PDB_small, PDB_small, is_anchor=False, anchor_name="test1")
        self.ag = self.universe.atoms[:20]  # prototypical AtomGroup
        self.ag_n = self.universe_n.atoms[:10]
        self.pickle_str = cPickle.dumps(self.ag, protocol=cPickle.HIGHEST_PROTOCOL)
        self.pickle_str_n = cPickle.dumps(self.ag_n, protocol=cPickle.HIGHEST_PROTOCOL)

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
        assert_(newag.universe is self.universe, "Unpickled AtomGroup on wrong Universe.")

    def test_unpickle_named(self):
        """Test that an AtomGroup can be unpickled (Issue 293)"""
        newag = cPickle.loads(self.pickle_str_n)
        # Can unpickle
        assert_array_equal(self.ag_n.indices, newag.indices)
        assert_(newag.universe is self.universe_n, "Unpickled AtomGroup on wrong Universe.")

    def test_unpickle_noanchor(self):
        # Shouldn't unpickle if the universe is removed from the anchors
        self.universe.remove_anchor()
        # In the complex (parallel) testing environment there's the risk of other compatible Universes being available
        #  for anchoring even after this one is expressly removed.
        assert_raises(RuntimeError, cPickle.loads, self.pickle_str)
          # If this fails to raise an exception either:"
          #  1-the anchoring Universe failed to remove_anchor or"
          #  2-another Universe with the same characteristics was created for testing and is being used as anchor."

    def test_unpickle_reanchor(self):
        # universe is removed from the anchors
        self.universe.remove_anchor()
        # now it goes back into the anchor list again
        self.universe.make_anchor()
        newag = cPickle.loads(self.pickle_str)
        assert_array_equal(self.ag.indices, newag.indices)
        assert_(newag.universe is self.universe, "Unpickled AtomGroup on wrong Universe.")

    def test_unpickle_reanchor_other(self):
        # universe is removed from the anchors
        self.universe.remove_anchor()
        # and universe_n goes into the anchor list
        self.universe_n.make_anchor()
        newag = cPickle.loads(self.pickle_str)
        assert_array_equal(self.ag.indices, newag.indices)
        assert_(newag.universe is self.universe_n, "Unpickled AtomGroup on wrong Universe.")

    def test_unpickle_wrongname(self):
        # we change the universe's anchor_name
        self.universe_n.anchor_name = "test2"
        # shouldn't unpickle if no name matches, even if there's a compatible
        #  universe in the unnamed anchor list.
        assert_raises(RuntimeError, cPickle.loads, self.pickle_str_n)

    def test_unpickle_rename(self):
        # we change universe_n's anchor_name
        self.universe_n.anchor_name = "test2"
        # and make universe a named anchor
        self.universe.anchor_name = "test1"
        newag = cPickle.loads(self.pickle_str_n)
        assert_array_equal(self.ag_n.indices, newag.indices)
        assert_(newag.universe is self.universe, "Unpickled AtomGroup on wrong Universe.")


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
    ref_unitcell = np.array([80.017, 80.017, 80.017, 60., 60., 90.], dtype=np.float32)
    ref_volume = 362270.0  # computed with Gromacs: 362.26999999999998 nm**3 * 1000 A**3/nm**3
    ref_offsets = None

    def setUp(self):
        # since offsets are automatically generated in the same directory
        # as the trajectory, we do everything from a temporary directory

        self.tmpdir = tempfile.mkdtemp()
        # loading from GRO is 4x faster than the PDB reader
        shutil.copy(GRO, self.tmpdir)
        shutil.copy(self.filename, self.tmpdir)

        self.top = os.path.join(self.tmpdir, os.path.basename(GRO))
        self.traj = os.path.join(self.tmpdir, os.path.basename(self.filename))

        self.universe = MDAnalysis.Universe(self.top, self.traj, convert_units=True)
        self.trajectory = self.universe.trajectory
        self.prec = 3
        self.ts = self.universe.coord

        # dummy output file
        ext = os.path.splitext(self.filename)[1]
        fd, self.outfile = tempfile.mkstemp(suffix=ext)
        os.close(fd)
        fd, self.outfile_offsets = tempfile.mkstemp(suffix='.pkl')
        os.close(fd)

    def tearDown(self):
        try:
            os.unlink(self.outfile)
            os.unlink(self.outfile_offsets)
            shutil.rmtree(self.tmpdir)
        except:
            pass
        del self.universe

    @dec.slow
    def test_offsets(self):
        if self.trajectory._offsets is None:
            self.trajectory.n_frames
        assert_array_almost_equal(self.trajectory._offsets, self.ref_offsets,
                                  err_msg="wrong frame offsets")

        # Saving
        self.trajectory.save_offsets(self.outfile_offsets)
        with open(self.outfile_offsets, 'rb') as f:
            saved_offsets = cPickle.load(f)
        assert_array_almost_equal(self.trajectory._offsets, saved_offsets['offsets'],
                                  err_msg="error saving frame offsets")
        assert_array_almost_equal(self.ref_offsets, saved_offsets['offsets'],
                                  err_msg="saved frame offsets don't match the known ones")

        # Loading
        self.trajectory.load_offsets(self.outfile_offsets)
        assert_array_almost_equal(self.trajectory._offsets, self.ref_offsets,
                                  err_msg="error loading frame offsets")

    @dec.slow
    def test_persistent_offsets_new(self):
        # check that offsets will be newly generated and not loaded from stored
        # offsets
        assert_equal(self.trajectory._offsets, None)

    @dec.slow
    def test_persistent_offsets_stored(self):
        # build offsets
        self.trajectory.n_frames
        assert_equal((self.trajectory._offsets is None), False)

        # check that stored offsets present
        assert_equal(os.path.exists(self.trajectory._offset_filename()), True)

    @dec.slow
    def test_persistent_offsets_ctime_match(self):
        # build offsets
        self.trajectory.n_frames

        with open(self.trajectory._offset_filename(), 'rb') as f:
            saved_offsets = cPickle.load(f)

        # check that stored offsets ctime matches that of trajectory file
        assert_equal(saved_offsets['ctime'], os.path.getctime(self.traj))

    @dec.slow
    def test_persistent_offsets_size_match(self):
        # build offsets
        self.trajectory.n_frames

        with open(self.trajectory._offset_filename(), 'rb') as f:
            saved_offsets = cPickle.load(f)

        # check that stored offsets size matches that of trajectory file
        assert_equal(saved_offsets['size'], os.path.getsize(self.traj))

    @dec.slow
    def test_persistent_offsets_autoload(self):
        # build offsets
        self.trajectory.n_frames

        # check that stored offsets are loaded for new universe
        u = MDAnalysis.Universe(self.top, self.traj)
        assert_equal((u.trajectory._offsets is not None), True)

    @dec.slow
    def test_persistent_offsets_ctime_mismatch(self):
        # build offsets
        self.trajectory.n_frames

        # check that stored offsets are not loaded when trajectory ctime
        # differs from stored ctime
        with open(self.trajectory._offset_filename(), 'rb') as f:
            saved_offsets = cPickle.load(f)
        saved_offsets['ctime'] = saved_offsets['ctime'] - 1
        with open(self.trajectory._offset_filename(), 'wb') as f:
            cPickle.dump(saved_offsets, f)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # Drop the warnings silently
            u = MDAnalysis.Universe(self.top, self.traj)
            assert_equal((u.trajectory._offsets is None), True)

    @dec.slow
    def test_persistent_offsets_size_mismatch(self):
        # build offsets
        self.trajectory.n_frames

        # check that stored offsets are not loaded when trajectory size differs
        # from stored size
        with open(self.trajectory._offset_filename(), 'rb') as f:
            saved_offsets = cPickle.load(f)
        saved_offsets['size'] += 1
        with open(self.trajectory._offset_filename(), 'wb') as f:
            cPickle.dump(saved_offsets, f)

        u = MDAnalysis.Universe(self.top, self.traj)
        assert_equal((u.trajectory._offsets is None), True)

    @dec.slow
    def test_persistent_offsets_last_frame_wrong(self):
        # build offsets
        self.trajectory.n_frames

        # check that stored offsets are not loaded when the offsets themselves
        # appear to be wrong
        with open(self.trajectory._offset_filename(), 'rb') as f:
            saved_offsets = cPickle.load(f)
        saved_offsets['offsets'] += 1
        with open(self.trajectory._offset_filename(), 'wb') as f:
            cPickle.dump(saved_offsets, f)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # Drop the warnings silently
            u = MDAnalysis.Universe(self.top, self.traj)
            assert_equal((u.trajectory._offsets is None), True)

    @dec.slow
    def test_persistent_offsets_readonly(self):
        # build offsets
        self.trajectory.n_frames

        # check that if directory is read-only offsets aren't stored
        os.unlink(self.trajectory._offset_filename())
        for root, dirs, files in os.walk(self.tmpdir, topdown=False):
            for item in dirs:
                os.chmod(os.path.join(root, item), 0444)
            for item in files:
                os.chmod(os.path.join(root, item), 0444)

        u = MDAnalysis.Universe(self.top, self.traj)
        assert_equal(os.path.exists(self.trajectory._offset_filename()), False)

    @dec.slow
    def test_persistent_offsets_refreshTrue(self):
        # build offsets
        self.trajectory.n_frames

        # check that the *refresh_offsets* keyword ensures stored offsets
        # aren't retrieved
        u = MDAnalysis.Universe(self.top, self.traj, refresh_offsets=True)
        assert_equal((u.trajectory._offsets is None), True)

    @dec.slow
    def test_persistent_offsets_refreshFalse(self):
        # build offsets
        self.trajectory.n_frames

        # check that the *refresh_offsets* keyword as False grabs offsets
        u = MDAnalysis.Universe(self.top, self.traj, refresh_offsets=False)
        assert_equal((u.trajectory._offsets is None), False)


class TestXTCReader_offsets(_GromacsReader_offsets):
    filename = XTC
    ref_offsets = np.array([0,  165188,  330364,  495520,  660708,  825872,  991044, 1156212,
                            1321384, 1486544])


class TestTRRReader_offsets(_GromacsReader_offsets):
    filename = TRR
    ref_offsets = np.array([0,  1144464,  2288928,  3433392,  4577856,  5722320,
                       6866784,  8011248,  9155712, 10300176])
