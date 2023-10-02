import pytest
import glob
import MDAnalysis as mda

from MDAnalysis.analysis.dssp import DSSP
from MDAnalysisTests.datafiles import DSSP_FOLDER
from MDAnalysisTests.datafiles import XTC_multi_frame


@pytest.mark.parametrize("pdb_filename", glob.glob(f"{DSSP_FOLDER}/*.pdb"))
def test_file_guess_hydrogens(pdb_filename):
    u = mda.Universe(pdb_filename)
    with open(f"{pdb_filename}.dssp", "r") as fin:
        correct_answ = fin.read().strip().split()[0]

    run = DSSP(u, guess_hydrogens=True).run()
    answ = "".join(run.results.dssp[0])
    assert answ == correct_answ


@pytest.mark.parametrize("pdb_filename", glob.glob(f"{DSSP_FOLDER}/*.pdb"))
def test_trajectory_without_guess_hydrogens(pdb_filename):
    u = mda.Universe(pdb_filename)
    DSSP(u, guess_hydrogens=False).run()

    # with_guess_s = "".join(DSSP(u, guess_hydrogens=True).run().results.dssp[0])
    # without_guess_s = "".join(DSSP(u, guess_hydrogens=False).run().results.dssp[0])
    # num_different = sum((s1 != s2 for s1, s2 in zip(with_guess_s, without_guess_s)))
    # assert num_different / len(with_guess_s) <= 0.05
