import pytest
import glob
import MDAnalysis as mda

from MDAnalysis.analysis.dssp import DSSP, translate
from MDAnalysisTests.datafiles import DSSP_FOLDER
from MDAnalysisTests.datafiles import TPR, XTC


@pytest.mark.parametrize("pdb_filename", glob.glob(f"{DSSP_FOLDER}/*.pdb.gz"))
def test_file_guess_hydrogens(pdb_filename):
    u = mda.Universe(pdb_filename)
    with open(f"{pdb_filename.rstrip('.gz')}.dssp", "r") as fin:
        correct_answ = fin.read().strip().split()[0]

    run = DSSP(u, guess_hydrogens=True).run()
    answ = "".join(run.results.dssp[0])
    assert answ == correct_answ


@pytest.mark.parametrize("pdb_filename", glob.glob(f"{DSSP_FOLDER}/*.pdb.gz"))
def test_file_without_guess_hydrogens(pdb_filename):
    u = mda.Universe(pdb_filename)
    DSSP(u, guess_hydrogens=False).run()


def test_trajectory():
    u = mda.Universe(TPR, XTC).select_atoms('protein').universe
    run = DSSP(u).run()
    first_frame = ''.join(run.results.dssp[0])
    last_frame = ''.join(run.results.dssp[-1])
    avg_frame = ''.join(translate(run.results.dssp_ndarray.mean(axis=0)))

    assert first_frame and last_frame and avg_frame
