import MDAnalysis as mda
import numpy as np

try:
    from MDAnalysisTests.datafiles import GRO
except:
    pass

try:
    from MDAnalysis.analysis import distances
except:
    pass


class DistancesBench(object):
    """Benchmarks for MDAnalysis.analysis.distances
    functions.
    """

    # TODO: eventually we should include box / pbc
    # unit cell information in the benchmarks
    params = (10, 100, 1000, 10000)
    param_names = ["num_atoms"]

    def setup(self, num_atoms):
        np.random.seed(17809)
        self.coords_1 = np.random.random_sample((num_atoms, 3)).astype(
            np.float32
        )
        np.random.seed(9008716)
        self.coords_2 = np.random.random_sample((num_atoms, 3)).astype(
            np.float32
        )
        np.random.seed(15809)
        self.coords_3 = np.random.random_sample((num_atoms, 3)).astype(
            np.float32
        )
        np.random.seed(871600)
        self.coords_4 = np.random.random_sample((num_atoms, 3)).astype(
            np.float32
        )

        self.allocated_array_2D = np.empty(
            (num_atoms, num_atoms), dtype=np.float64
        )
        self.array_shape_1D = int(num_atoms * (num_atoms - 1) / 2.0)
        self.allocated_array_1D = np.empty(
            self.array_shape_1D, dtype=np.float64
        )
        self.u = mda.Universe(GRO)
        self.ag1 = self.u.atoms[:num_atoms]
        self.ag2 = self.u.atoms[num_atoms: 2 * num_atoms]
        self.ag3 = self.u.atoms[-num_atoms:]

    def time_distance_array(self, num_atoms):
        """Benchmark calculation of all distances
        between two numpy arrays of coordinates,
        using default arguments to distance_array.
        """
        distances.distance_array(
            reference=self.coords_1,
            configuration=self.coords_2,
            box=None,
            result=None,
            backend="serial",
        )

    def time_distance_array_pre_allocated(self, num_atoms):
        """Benchmark calculation of all distances
        between two numpy arrays of coordinates,
        using distance_array with a preallocated
        result array.
        """
        distances.distance_array(
            reference=self.coords_1,
            configuration=self.coords_2,
            box=None,
            result=self.allocated_array_2D,
            backend="serial",
        )

    def time_self_distance_array(self, num_atoms):
        """Benchmark calculation of all distances
        within a single numpy array of coordinates
        using default parameters to self_distance_array.
        """
        distances.self_distance_array(
            reference=self.coords_1, box=None, result=None, backend="serial"
        )

    def time_self_distance_array_pre_allocated(self, num_atoms):
        """Benchmark calculation of all distances
        within a single numpy array of coordinates
        using self_distance_array with preallocated
        result array.
        """
        distances.self_distance_array(
            reference=self.coords_1,
            box=None,
            result=self.allocated_array_1D,
            backend="serial",
        )

    def time_contact_matrix(self, num_atoms):
        """Benchmark calculation of contacts within
        a single numpy array using the default arguments
        to contact_matrix.
        """
        distances.contact_matrix(
            coord=self.coords_1, cutoff=15.0, returntype="numpy", box=None
        )

    def time_contact_matrix_sparse(self, num_atoms):
        """Benchmark calculation of contacts within
        a single numpy array using the slower reduced
        memory implementation of contact_matrix intended
        for larger systems.
        """
        distances.contact_matrix(
            coord=self.coords_1, cutoff=15.0, returntype="sparse", box=None
        )

    def time_dist(self, num_atoms):
        """Benchmark calculation of distances between
        atoms in two atomgroups with no offsets
        to resids.
        """
        distances.dist(A=self.ag1, B=self.ag2, offset=0)

    def time_dist_offsets(self, num_atoms):
        """Benchmark calculation of distances between
        atoms in two atomgroups with offsets
        to resids.
        """
        distances.dist(A=self.ag1, B=self.ag2, offset=20)

    def time_between(self, num_atoms):
        """Benchmark determination of subgroup
        of atomgroup that is within a specific
        distance of two other atomgroups.
        """
        distances.between(
            group=self.ag3, A=self.ag1, B=self.ag2, distance=15.0
        )

    def time_calc_bonds(self, num_atoms):
        """Benchmark calculation of bonds between
        atoms in two atomgroups.
        """
        mda.lib.distances.calc_bonds(self.coords_1, self.coords_2)

    def time_calc_angles(self, num_atoms):
        """Benchmark calculation of angles between
        atoms in three atomgroups.
        """
        mda.lib.distances.calc_angles(
            self.coords_1, self.coords_2, self.coords_3
        )

    def time_calc_dihedrals(self, num_atoms):
        """Benchmark calculation of dihedrals between
        atoms in four atomgroups.
        """
        mda.lib.distances.calc_dihedrals(
            self.coords_1, self.coords_2, self.coords_3, self.coords_4
        )
