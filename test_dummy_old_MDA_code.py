import MDAnalysis
from MDAnalysis.tests.datafiles import GRO, XTC
universe = MDAnalysis.Universe(GRO, XTC)

#old selection
all_selection = universe.selectAtoms('all')

#additional old selectAtoms selection (this comment shouldn't be modified despite containing the method name)
all_selection.selectAtoms('bynum 1:10')

#testing atomgroup methods to properties (and exclusion of comments from conversion):

#all_selection.residues()
all_selection.residues()
#all_selection.charges()
all_selection.charges()
#all_selection.indices()
all_selection.indices()
#all_selection.masses()
all_selection.masses()
#all_selection.names()
all_selection.names()
#all_selection.types()
all_selection.types()
#all_selection.radii()
all_selection.radii()
#all_selection.resids()
all_selection.resids()
#all_selection.resnames()
all_selection.resnames()
#all_selection.resnums()
all_selection.resnums()
#all_selection.segids()
all_selection.segids()

#similarly for atomgroup count method renaming:

#all_selection.numberOfAtoms()
all_selection.numberOfAtoms()

#all_selection.numberOfResidues()
all_selection.numberOfResidues()

#all_selection.numberOfSegments()
all_selection.numberOfSegments()
