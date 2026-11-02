import numpy as np
from MDAnalysis.lib.pkdtree import PeriodicKDTree
from MDAnalysis.lib.distances import capped_distance
from scipy.spatial import cKDTree


class NeighborsBench:
    """Benchmarks for neighbor searching functions."""

    params = ([100, 1000, 10000, 100000], [20, 30, 36, 42, 48, 50, 60])
    param_names = ["number_of_atoms", "cutoff"]

    def setup(self, number_of_atoms, cutoff):
        """Setup called before each benchmark with each parameter combination."""
        self.box = np.array(
            [170.0, 70.0, 120.0, 90.0, 90.0, 90.0], dtype=np.float32
        )
        self.positions = (
            np.random.rand(number_of_atoms, 3) * self.box[:3]
        ).astype(np.float32)
        self.centre = (self.box[:3] / 2.0).reshape(1, 3)
        self.cutoff = cutoff

        self.scipy_tree = cKDTree(self.positions, boxsize=self.box[:3])
        self.mda_tree = PeriodicKDTree(box=self.box)
        self.mda_tree.set_coords(self.positions, cutoff=self.cutoff)

    def time_mda_tree_search(self, number_of_atoms, cutoff):
        """Benchmark just the search operation on pre-built tree."""
        self.mda_tree.search(self.centre, self.cutoff)

    def time_scipy_tree_query(self, number_of_atoms, cutoff):
        """Benchmark just the query operation on pre-built tree."""
        self.scipy_tree.query_ball_point(self.centre, self.cutoff)

    def time_mda_PKDtree_with_setup(self, number_of_atoms, cutoff):
        """Benchmark tree construction + search."""
        tree = PeriodicKDTree(box=self.box)
        tree.set_coords(self.positions, cutoff=self.cutoff)
        tree.search(self.centre, self.cutoff)

    def time_scipy_cKDTree_with_setup(self, number_of_atoms, cutoff):
        """Benchmark tree construction + query."""
        tree = cKDTree(self.positions, boxsize=self.box[:3])
        tree.query_ball_point(self.centre, self.cutoff)

    def time_capped_distance_array(self, number_of_atoms, cutoff):
        """Benchmark capped distance calculation."""
        capped_distance(
            self.centre, self.positions, max_cutoff=self.cutoff, box=self.box
        )
