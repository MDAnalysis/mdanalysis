import numpy as np
from numpy.testing import assert_array_almost_equal
import pytest
import os

import MDAnalysis as mda
from MDAnalysisTests.datafiles import GRO, XTC, DihedralsArray, GLYDihedralsArray
from MDAnalysis.analysis import dihedrals


class TestRamachandran(object):

    @pytest.fixture()
    def universe(self):
        return mda.Universe(GRO, XTC)

    def test_ramachandran(self, universe):
        rama = dihedrals.Ramachandran(universe.select_atoms("protein")).run()
        test_rama = np.load(DihedralsArray)

        assert_array_almost_equal(rama.angles, test_rama, 5,
                                  err_msg="error: dihedral angles should "
                                  "match test values")

    def test_ramachandran_single_frame(self, universe):
        rama = dihedrals.Ramachandran(universe.select_atoms("protein"),
                                      start=5, stop=6).run()
        test_rama = [np.load(DihedralsArray)[5]]

        assert_array_almost_equal(rama.angles, test_rama, 5,
                                  err_msg="error: dihedral angles should "
                                  "match test values")

    def test_ramachandran_identical_frames(self, universe, tmpdir):

        outfile = os.path.join(str(tmpdir), "dihedrals.xtc")

        # write a dummy trajectory of all the same frame
        with mda.Writer(outfile, universe.atoms.n_atoms) as W:
            for _ in range(universe.trajectory.n_frames):
                W.write(universe)

        universe = mda.Universe(GRO, outfile)
        rama = dihedrals.Ramachandran(universe.select_atoms("protein")).run()
        test_rama = [np.load(DihedralsArray)[0] for ts in universe.trajectory]

        assert_array_almost_equal(rama.angles, test_rama, 5,
                                  err_msg="error: dihedral angles should "
                                  "match test values")

    def test_ramachandran_residue_selections(self, universe):
        rama = dihedrals.Ramachandran(universe.select_atoms("resname GLY")).run()
        test_rama = np.load(GLYDihedralsArray)

        assert_array_almost_equal(rama.angles, test_rama, 5,
                                  err_msg="error: dihedral angles should "
                                  "match test values")
