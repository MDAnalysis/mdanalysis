from asv_runner.benchmarks.mark import SkipNotImplemented
import numpy as np

try:
    from MDAnalysis.lib import distances, mdamath
except ImportError:
    pass


class DistancesMinimizeVectorsBench(object):
    """Benchmarks for MDAnalysis.lib.distances.minimize_vectors
    function.
    """

    params = (
        [10, 100, 1000, 10000, 100000, 1000000],
        ["ortho", "triclinic"],
        ["serial", "openmp", "distopia"],
    )
    param_names = ["num_atoms", "box", "backend"]

    def setup(self, num_atoms, box, backend):
        np.random.seed(17809)
        self.coords_1 = np.random.random_sample((num_atoms, 3)).astype(np.float32)
        np.random.seed(9008716)
        self.coords_2 = np.random.random_sample((num_atoms, 3)).astype(np.float32)
        if box == "ortho":
            box = np.eye(3) * 10
        elif box == "triclinic":
            box = np.array([[10, 0, 0], [2, 10, 0], [2, 2, 10]])
        else:
            raise ValueError(f"box needs to be ortho or triclinic, passed {box}")
        self.box = mdamath.triclinic_box(*box)

    def time_minimize_vectors(self, num_atoms, box, backend):
        """Benchmark calculation of minimized vectors
        based on the unitcell dimensions box.
        """
        vectors = self.coords_1 - self.coords_2
        try:
            distances.minimize_vectors(vectors, self.box, backend=backend)
        except TypeError:
            raise SkipNotImplemented("minimize_vectors is skipped")
