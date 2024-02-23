import pytest
import glob
import MDAnalysis as mda

from MDAnalysis.analysis.dssp import DSSP, translate
from MDAnalysisTests.datafiles import DSSP as DSSP_FOLDER
from MDAnalysisTests.datafiles import TPR, XTC


@pytest.mark.parametrize("pdb_filename", glob.glob(f"{DSSP_FOLDER}/?????.pdb.gz"))
def test_file_guess_hydrogens(pdb_filename):
    u = mda.Universe(pdb_filename)
    with open(f"{pdb_filename.rstrip('.gz')}.dssp", "r") as fin:
        correct_answ = fin.read().strip().split()[0]

    run = DSSP(u, guess_hydrogens=True).run()
    answ = "".join(run.results.dssp[0])
    assert answ == correct_answ


def test_trajectory():
    u = mda.Universe(TPR, XTC).select_atoms('protein').universe
    run = DSSP(u).run(stop=10)
    first_frame = ''.join(run.results.dssp[0])
    last_frame = ''.join(run.results.dssp[-1])
    avg_frame = ''.join(translate(run.results.dssp_ndarray.mean(axis=0)))

    assert first_frame[:10] != last_frame[:10] == avg_frame[:10] == '-EEEEEE---'


def test_trajectory_with_hydrogens():
    u = mda.Universe(TPR, XTC).select_atoms('protein').universe
    run = DSSP(u, guess_hydrogens=False).run(stop=10)
    first_frame = ''.join(run.results.dssp[0])
    last_frame = ''.join(run.results.dssp[-1])
    avg_frame = ''.join(translate(run.results.dssp_ndarray.mean(axis=0)))

    assert first_frame[:10] == last_frame[:10] == avg_frame[:10] == '-EEEEEE---'


@pytest.mark.parametrize("pdb_filename", glob.glob(f"{DSSP_FOLDER}/2xdgA.pdb.gz"))
def test_trajectory_without_hydrogen_fails(pdb_filename):
    u = mda.Universe(pdb_filename)
    with pytest.raises(ValueError):
        DSSP(u, guess_hydrogens=False)


@pytest.mark.parametrize("pdb_filename", glob.glob(f"{DSSP_FOLDER}/1mr1D_failing.pdb.gz"))
def test_trajectory_with_uneven_number_of_atoms_fails(pdb_filename):
    u = mda.Universe(pdb_filename)
    with pytest.raises(ValueError):
        DSSP(u, guess_hydrogens=False)
