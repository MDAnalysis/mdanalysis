import numpy as np
from MDAnalysis.lib.nsgrid import FastNS


class FastNSBench:

    params = [100, 500, 1000, 10000]
    param_names = ["n_atoms"]

    def setup(self, n_atoms):
        self.coords = np.random.random((n_atoms, 3)).astype(np.float32)
        self.box = np.array([10, 10, 10, 90, 90, 90], dtype=np.float32)

    def time_fastns(self, n_atoms):
        ns = FastNS(5.0, self.coords, self.box)
        ns.self_search()