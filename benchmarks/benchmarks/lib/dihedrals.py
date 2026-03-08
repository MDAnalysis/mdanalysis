import numpy as np
from MDAnalysis.lib.distances import calc_dihedrals


class DihedralBench:

    params = [100, 500, 1000, 10000]
    param_names = ["n_dihedrals"]

    def setup(self, n_dihedrals):
        self.coords = np.random.random((n_dihedrals, 4, 3))

    def time_calc_dihedrals(self, n_dihedrals):
        calc_dihedrals(
            self.coords[:,0],
            self.coords[:,1],
            self.coords[:,2],
            self.coords[:,3]
        )