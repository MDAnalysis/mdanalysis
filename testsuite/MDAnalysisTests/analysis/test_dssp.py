import MDAnalysis as mda
import pytest
from MDAnalysis.analysis.dssp import DSSP
import glob

from MDAnalysisTests.datafiles import DSSP_FOLDER


@pytest.mark.parametrize("pdb_filename", glob.glob(f"{DSSP_FOLDER}/*.pdb"))
def test_file(pdb_filename):
    u = mda.Universe(pdb_filename)
    with open(f"{pdb_filename}.dssp", "r") as fin:
        correct_answ = fin.read().strip().split()[0]

    run = DSSP(u).run()
    answ = "".join(run.results.dssp[0])
    assert answ == correct_answ
