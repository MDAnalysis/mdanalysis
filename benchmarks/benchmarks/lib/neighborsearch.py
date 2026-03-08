import numpy as np
from MDAnalysis.lib.NeighborSearch import AtomNeighborSearch
import MDAnalysis as mda
from MDAnalysisTests.datafiles import PSF, DCD


class NeighborSearchBench:

    def setup(self):
        self.u = mda.Universe(PSF, DCD)
        self.atoms = self.u.select_atoms("protein")
        self.ns = AtomNeighborSearch(self.atoms)

    def time_neighbor_search(self):
        self.ns.search(self.atoms, 5.0)