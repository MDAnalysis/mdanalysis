import MDAnalysis as mda
import numpy as np

from numpy.testing import assert_equal, assert_raises, assert_almost_equal
from unittest import TestCase

from MDAnalysisTests.coordinates.reference import RefAdKSmall


class _SingleFrameReader(TestCase, RefAdKSmall):
    # see TestPDBReader how to set up!

    def tearDown(self):
        del self.universe

    def test_flag_permissive_pdb_reader(self):
        """test_flag_permissive_pdb_reader: permissive_pdb_reader==True enables
        primitive PDB reader"""
        assert_equal(mda.core.flags['permissive_pdb_reader'], True,
                     "'permissive_pdb_reader' flag should be True as "
                     "MDAnalysis default")

    def test_load_file(self):
        U = self.universe
        assert_equal(len(U.atoms), self.ref_n_atoms,
                     "load Universe from file %s" % U.trajectory.filename)
        assert_equal(U.atoms.select_atoms('resid 150 and name HA2').atoms[0],
                     U.atoms[self.ref_E151HA2_index], "Atom selections")

    def test_n_atoms(self):
        assert_equal(self.universe.trajectory.n_atoms, self.ref_n_atoms,
                     "wrong number of atoms")

    def test_numres(self):
        assert_equal(self.universe.atoms.n_residues, 214,
                     "wrong number of residues")

    def test_n_frames(self):
        assert_equal(self.universe.trajectory.n_frames, 1,
                     "wrong number of frames in pdb")

    def test_time(self):
        assert_equal(self.universe.trajectory.time, 0.0,
                     "wrong time of the frame")

    def test_frame(self):
        assert_equal(
            self.universe.trajectory.frame, 0,
            "wrong frame number (0-based, should be 0 for single frame "
            "readers)")

    def test_frame_index_0(self):
        self.universe.trajectory[0]
        assert_equal(self.universe.trajectory.ts.frame, 0,
                     "frame number for frame index 0 should be 0")

    def test_frame_index_1_raises_IndexError(self):
        def go_to_2(traj=self.universe.trajectory):
            traj[1]

        assert_raises(IndexError, go_to_2)

    def test_dt(self):
        """testing that accessing universe.trajectory.dt gives 1.0
        (the default)"""
        assert_equal(1.0, self.universe.trajectory.dt)

    def test_coordinates(self):
        A10CA = self.universe.atoms.CA[10]
        # restrict accuracy to maximum in PDB files (3 decimals)
        assert_almost_equal(A10CA.pos,
                            self.ref_coordinates['A10CA'],
                            3,
                            err_msg="wrong coordinates for A10:CA")

    def test_distances(self):
        NTERM = self.universe.atoms.N[0]
        CTERM = self.universe.atoms.C[-1]
        d = mda.lib.mdamath.norm(NTERM.position - CTERM.position)
        assert_almost_equal(d,
                            self.ref_distances['endtoend'],
                            self.prec,
                            err_msg="distance between M1:N and G214:C")

    def test_full_slice(self):
        trj_iter = self.universe.trajectory[:]
        frames = [ts.frame for ts in trj_iter]
        assert_equal(frames, np.arange(self.universe.trajectory.n_frames))

    def test_last_slice(self):
        # should be same as above: only 1 frame!
        trj_iter = self.universe.trajectory[-1:]
        frames = [ts.frame for ts in trj_iter]
        assert_equal(frames, np.arange(self.universe.trajectory.n_frames))
