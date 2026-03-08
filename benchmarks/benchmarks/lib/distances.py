import numpy as np
from MDAnalysis.lib.distances import distance_array, self_distance_array


class DistanceArrayBench:

    params = [100, 500, 1000, 10000]
    param_names = ["n_atoms"]

    def setup(self, n_atoms):
        self.coords1 = np.random.random((n_atoms, 3))
        self.coords2 = np.random.random((n_atoms, 3))

    def time_distance_array(self, n_atoms):
        distance_array(self.coords1, self.coords2)

class SelfDistanceArrayBench:

    params = [100, 500, 1000, 10000]
    param_names = ["n_atoms"]

    def setup(self, n_atoms):
        self.coords = np.random.random((n_atoms, 3))

    def time_self_distance_array(self, n_atoms):
        self_distance_array(self.coords)