import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.rmsf_residue import RMSFResidue
from MDAnalysisTests.datafiles import PSF, DCD


def test_rmsf_residue_runs():
    # Load tiny test system included in MDAnalysis
    u = mda.Universe(PSF, DCD)

    # Create the analysis
    rmsf = RMSFResidue(u, select="all")

    # Should run without errors
    rmsf.run()

    # Check that results exist
    assert hasattr(rmsf.results, "residue_rmsf")

    # Number of residues should match universe
    n_res = len(u.residues)
    assert rmsf.results.residue_rmsf.shape[0] == n_res

    # RMSF values must be >= 0
    assert np.all(rmsf.results.residue_rmsf >= 0)
