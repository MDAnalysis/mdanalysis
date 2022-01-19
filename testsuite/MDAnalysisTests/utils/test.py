import MDAnalysis as mda
from openmm import app
from MDAnalysisTests.datafiles import PDBX

u = mda.Universe(app.PDBxFile(PDBX))
