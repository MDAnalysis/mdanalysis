import pytest
import glob
import MDAnalysis as mda

from MDAnalysis.analysis.dssp import DSSP, translate, HAS_EINOPS
from MDAnalysisTests.datafiles import DSSP as DSSP_FOLDER
from MDAnalysisTests.datafiles import TPR, XTC


@pytest.mark.skipif(HAS_EINOPS, reason="einops present")
def test_without_einops():
    u = mda.Universe(TPR, XTC).select_atoms('protein').universe
    with pytest.raises(ImportError):
        run = DSSP(u).run()

    first_frame = ''.join(run.results.dssp[0])
    last_frame = ''.join(run.results.dssp[-1])
    avg_frame = ''.join(translate(run.results.dssp_ndarray.mean(axis=0)))

    assert first_frame and last_frame and avg_frame


@pytest.mark.skipif(not HAS_EINOPS, reason="einops present")
@pytest.mark.parametrize("pdb_filename", glob.glob(f"{DSSP_FOLDER}/*.pdb.gz"))
def test_file_guess_hydrogens(pdb_filename):
    u = mda.Universe(pdb_filename)
    with open(f"{pdb_filename.rstrip('.gz')}.dssp", "r") as fin:
        correct_answ = fin.read().strip().split()[0]

    run = DSSP(u, guess_hydrogens=True).run()
    answ = "".join(run.results.dssp[0])
    assert answ == correct_answ


@pytest.mark.skipif(not HAS_EINOPS, reason="einops present")
@pytest.mark.parametrize("pdb_filename", glob.glob(f"{DSSP_FOLDER}/*.pdb.gz"))
def test_file_without_guess_hydrogens(pdb_filename):
    u = mda.Universe(pdb_filename)
    DSSP(u, guess_hydrogens=False).run()


@pytest.mark.skipif(not HAS_EINOPS, reason="einops present")
def test_trajectory():
    u = mda.Universe(TPR, XTC).select_atoms('protein').universe
    run = DSSP(u).run()
    first_frame = ''.join(run.results.dssp[0])
    last_frame = ''.join(run.results.dssp[-1])
    avg_frame = ''.join(translate(run.results.dssp_ndarray.mean(axis=0)))

    assert first_frame and last_frame and avg_frame
