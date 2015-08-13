import MDAnalysis
from MDAnalysis.tests.datafiles import GRO, XTC
universe = MDAnalysis.Universe(GRO, XTC)

#old selection
all_selection = universe.selectAtoms('all')

#additional old selectAtoms selection (this comment shouldn't be modified despite containing the method name)
all_selection.selectAtoms('bynum 1:10')
