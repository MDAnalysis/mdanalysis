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
import tempfile
import os
from numpy.testing import *

from MDAnalysisTests.datafiles import mol2_molecules, mol2_molecule, mol2_broken_molecule
from MDAnalysis import Universe


class TestMol2(TestCase):
    def setUp(self):
        fd, self.outfile = tempfile.mkstemp(suffix=".mol2")
        os.close(fd)

    def tearDown(self):
        try:
            os.unlink(self.outfile)
        except OSError:
            pass

    def test_read(self):
        u = Universe(mol2_molecules)
        assert_equal(len(u.atoms), 49)
        assert_equal(u.trajectory.n_frames, 200)

        u.trajectory[199]
        assert_array_almost_equal(u.atoms.positions[0], [1.7240, 11.2730, 14.1200])

    def test_write(self):
        ref = Universe(mol2_molecules)
        ref.atoms.write(self.outfile)
        u = Universe(self.outfile)
        assert_equal(len(u.atoms), len(ref.atoms))
        assert_equal(len(u.trajectory), len(ref.trajectory))
        assert_array_equal(u.atoms.positions, ref.atoms.positions)
        u.trajectory[199]
        ref.trajectory[199]
        assert_array_equal(u.atoms.positions, ref.atoms.positions)

    def test_write_selection(self):
        ref = Universe(mol2_molecule)
        gr0 = ref.select_atoms("name C*")
        gr0.write(self.outfile)
        u = Universe(self.outfile)
        gr1 = u.select_atoms("name C*")
        assert_equal(len(gr0), len(gr1))

    def test_broken_molecule(self):
        assert_raises(ValueError, Universe, mol2_broken_molecule)

        # This doesn't work with 2.6
        # Checks the text of the error message, so it low priority
        #with self.assertRaises(Exception) as context:
        #    u = Universe(mol2_broken_molecule)
        #self.assertEqual("The mol2 block (BrokenMolecule.mol2:0) has no atoms" in context.exception.message,
        # True)


class TestMol2_traj(TestCase):
    def setUp(self):
        self.universe = Universe(mol2_molecules)
        self.traj = self.universe.trajectory
        self.ts = self.universe.coord

    def tearDown(self):
        del self.universe
        del self.traj
        del self.ts

    def test_rewind_traj(self):
        self.traj.rewind()
        assert_equal(self.ts.frame, 0, "rewinding to frame 0")

    def test_next_traj(self):
        self.traj.rewind()
        self.traj.next()
        assert_equal(self.ts.frame, 1, "loading frame 1")

    def test_jump_traj(self):
        self.traj[15]  # index is 0-based and frames are 0-based
        assert_equal(self.ts.frame, 15, "jumping to frame 15")

    def test_jump_lastframe_traj(self):
        self.traj[-1]
        assert_equal(self.ts.frame, 199, "indexing last frame with traj[-1]")

    def test_slice_traj(self):
        frames = [ts.frame for ts in self.traj[5:17:3]]
        assert_equal(frames, [5, 8, 11, 14], "slicing traj [5:17:3]")

    def test_reverse_traj(self):
        frames = [ts.frame for ts in self.traj[20:5:-1]]
        assert_equal(frames, range(20, 5, -1), "reversing traj [20:5:-1]")

    def test_n_frames(self):
        assert_equal(self.universe.trajectory.n_frames, 200, "wrong number of frames in traj")
