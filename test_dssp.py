from MDAnalysis.analysis.dssp import DSSP
import MDAnalysis as mda
u = mda.Universe('large_data/NhaA_non_water.gro', 'large_data/NhaA_non_water.xtc')
DSSP(u, guess_hydrogens=False).run(stop=20)
