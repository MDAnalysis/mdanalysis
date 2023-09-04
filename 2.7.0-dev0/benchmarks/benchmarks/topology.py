import MDAnalysis
import numpy as np
from MDAnalysis.topology import guessers

try:
    from MDAnalysisTests.datafiles import GRO
    from MDAnalysis.exceptions import NoDataError
except:
    pass

class TopologyGuessBench(object):
    """Benchmarks for individual
    topology functions
    """
    params = (10, 100, 1000, 10000)
    param_names = ['num_atoms']
    
    def setup(self, num_atoms):
        self.u = MDAnalysis.Universe(GRO)
        self.ag = self.u.atoms[:num_atoms]
        self.vdwradii = {'H':1.0,
                         'C':1.0,
                         'N':1.0,
                         'O':1.0,
                         'DUMMY':1.0}

    def time_guessbonds(self, num_atoms):
        """Benchmark for guessing bonds"""
        guessers.guess_bonds(self.ag, self.ag.positions,
                             box=self.ag.dimensions,
                             vdwradii=self.vdwradii)


class BondsBench(object):
    """
    Benchmarks for getting bonds
    """

    params = (1000, 10000, 100000, 1000000)
    param_names = ['num_bonds']

    def setup(self, num_bonds):
        self.u = MDAnalysis.Universe.empty(2*num_bonds)
        bonds = np.arange(2*num_bonds).reshape(num_bonds, 2)
        self.u.add_bonds(bonds)

    def time_bonds(self, num_bonds):
        """Benchmark for calculating bonds"""
        self.u.bonds
