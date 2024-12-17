import pytest
import glob
import MDAnalysis as mda

from MDAnalysis.analysis.dssp import DSSP, translate
from MDAnalysisTests.datafiles import DSSP as DSSP_FOLDER
from MDAnalysisTests.datafiles import TPR, XTC


# Files that match glob pattern '????.pdb.gz' and matching '????.pdb.dssp' files,
# containing the secondary structure assignment string, will be tested automatically.
@pytest.mark.parametrize("pdb_filename", glob.glob(f"{DSSP_FOLDER}/?????.pdb.gz"))
def test_file_guess_hydrogens(pdb_filename, client_DSSP):
    u = mda.Universe(pdb_filename)
    with open(f"{pdb_filename.rstrip('.gz')}.dssp", "r") as fin:
        correct_answ = fin.read().strip().split()[0]

    run = DSSP(u, guess_hydrogens=True).run(**client_DSSP)
    answ = "".join(run.results.dssp[0])
    assert answ == correct_answ


def test_trajectory(client_DSSP):
    u = mda.Universe(TPR, XTC).select_atoms("protein").universe
    run = DSSP(u).run(**client_DSSP, stop=10)
    first_frame = "".join(run.results.dssp[0])
    last_frame = "".join(run.results.dssp[-1])
    avg_frame = "".join(translate(run.results.dssp_ndarray.mean(axis=0)))

    assert first_frame[:10] != last_frame[:10] == avg_frame[:10] == "-EEEEEE---"
    protein = mda.Universe(TPR, XTC).select_atoms("protein")
    run = DSSP(protein).run(**client_DSSP, stop=10)


def test_atomgroup(client_DSSP):
    protein = mda.Universe(TPR, XTC).select_atoms("protein")
    run = DSSP(protein).run(**client_DSSP, stop=10)
    first_frame = "".join(run.results.dssp[0])
    last_frame = "".join(run.results.dssp[-1])
    avg_frame = "".join(translate(run.results.dssp_ndarray.mean(axis=0)))

    assert first_frame[:10] != last_frame[:10] == avg_frame[:10] == "-EEEEEE---"


def test_trajectory_with_hydrogens(client_DSSP):
    u = mda.Universe(TPR, XTC).select_atoms("protein").universe
    run = DSSP(u, guess_hydrogens=False).run(**client_DSSP, stop=10)
    first_frame = "".join(run.results.dssp[0])
    last_frame = "".join(run.results.dssp[-1])
    avg_frame = "".join(translate(run.results.dssp_ndarray.mean(axis=0)))

    assert first_frame[:10] == last_frame[:10] == avg_frame[:10] == "-EEEEEE---"


@pytest.mark.parametrize("pdb_filename", glob.glob(f"{DSSP_FOLDER}/2xdgA.pdb.gz"))
def test_trajectory_without_hydrogen_fails(pdb_filename, client_DSSP):
    u = mda.Universe(pdb_filename)
    with pytest.raises(ValueError):
        DSSP(u, guess_hydrogens=False).run(**client_DSSP)


@pytest.mark.parametrize(
    "pdb_filename", glob.glob(f"{DSSP_FOLDER}/1mr1D_failing.pdb.gz")
)
def test_trajectory_with_uneven_number_of_atoms_fails(pdb_filename,
                                                      client_DSSP):
    u = mda.Universe(pdb_filename)
    with pytest.raises(ValueError):
        DSSP(u, guess_hydrogens=True).run(**client_DSSP)


@pytest.mark.parametrize(
    "pdb_filename", glob.glob(f"{DSSP_FOLDER}/wrong_hydrogens.pdb.gz")
)
def test_exception_raises_with_atom_index(pdb_filename, client_DSSP):
    u = mda.Universe(pdb_filename)
    with pytest.raises(
        ValueError,
        match="Residue <Residue SER, 298> contains*",
    ):
        DSSP(u, guess_hydrogens=False).run(**client_DSSP)
