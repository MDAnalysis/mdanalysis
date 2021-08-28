import MDAnalysis
import numpy as np


try:
    from MDAnalysis.lib import mdamath
except:
    pass

try:
    from MDAnalysisTests.datafiles import GRO
except:
    pass

try:
    from MDAnalysis.analysis import distances
except:
    pass


class BetweenBench(object):
    """Benchmark for the MDAnalysis.analysis.distances.between
    function.
    """
    u = MDAnalysis.Universe(GRO)
    tri_clinic_dimensions = u.dimensions

    params = ([10, 100, 1000, 10000])
    param_names = (['num_atoms'])

    def setup(self, num_atoms):

        self.u = MDAnalysis.Universe(GRO)

        self.ag1 = self.u.atoms[:num_atoms]
        self.ag2 = self.u.atoms[num_atoms: 2 * num_atoms]
        self.ag3 = self.u.atoms[-num_atoms:]

        np.random.seed(17809)
        self.coords_1 = np.random.random_sample((num_atoms, 3)).astype(np.float32)
        np.random.seed(9008716)
        self.coords_2 = np.random.random_sample((num_atoms, 3)).astype(np.float32)
        self.allocated_array_2D = np.empty((num_atoms, num_atoms),
                                           dtype=np.float64)
        self.array_shape_1D = int(num_atoms * (num_atoms - 1) / 2.)
        self.allocated_array_1D = np.empty(self.array_shape_1D,
                                           dtype=np.float64)

    def time_between(self, num_atoms):
        """Benchmark determination of subgroup
        of atomgroup that is within a specific
        distance of two other atomgroups.
        """
        distances.between(group=self.ag3,
                          A=self.ag1,
                          B=self.ag2,
                          distance=15.0)


class DistancesBench(object):
    """Benchmarks for MDAnalysis.analysis.distances
    functions. Excluding contact matrices.
    """
    u = MDAnalysis.Universe(GRO)
    tri_clinic_dimensions = u.dimensions

    params = ([10, 100, 1000, 10000], [None, 'orthogonal', 'triclinic'])
    param_names = ['num_atoms', 'pbc_type']

    def setup(self, num_atoms, pbc_type):

        self.u = MDAnalysis.Universe(GRO)

        if pbc_type is None:
            self.box_dims = None

        elif pbc_type == 'orthogonal':
            box_shape = self.u.dimensions
            basis_vectors = mdamath.triclinic_vectors(box_shape)

            # this will calculate the furthest point from the origin on
            # the triclinic unit cell.
            furthest_point = np.sum(basis_vectors, axis=0)
            # get full size of cube and also deal with certain distance
            # functions requiring float32
            orthogonal_box_size = np.float32(np.max(furthest_point))
            angle = np.float32(90.0)

            # make dummy cubic box, must be np.array, does not accept list
            self.box_dims = np.array([orthogonal_box_size,
                                      orthogonal_box_size,
                                      orthogonal_box_size,
                                      angle, angle, angle])

        elif pbc_type == 'triclinic':
            self.box_dims = self.u.dimensions

        else:
            raise Exception("Specify only 'orthogonal', 'triclinic', or None for the box type")

        self.ag1 = self.u.atoms[:num_atoms]
        self.ag2 = self.u.atoms[num_atoms: 2 * num_atoms]
        self.ag3 = self.u.atoms[-num_atoms:]

        np.random.seed(17809)
        self.coords_1 = np.random.random_sample((num_atoms, 3)).astype(np.float32)
        np.random.seed(9008716)
        self.coords_2 = np.random.random_sample((num_atoms, 3)).astype(np.float32)
        self.allocated_array_2D = np.empty((num_atoms, num_atoms),
                                           dtype=np.float64)
        self.array_shape_1D = int(num_atoms * (num_atoms - 1) / 2.)
        self.allocated_array_1D = np.empty(self.array_shape_1D,
                                           dtype=np.float64)

    def time_distance_array(self, num_atoms, pbc_type):
        """Benchmark calculation of all distances
        between two numpy arrays of coordinates,
        using default arguments to distance_array.
        """
        distances.distance_array(reference=self.coords_1,
                                 configuration=self.coords_2,
                                 box=self.box_dims,
                                 result=None,
                                 backend='serial')

    def time_distance_array_pre_allocated(self, num_atoms, pbc_type):
        """Benchmark calculation of all distances
        between two numpy arrays of coordinates,
        using distance_array with a preallocated
        result array.
        """
        distances.distance_array(reference=self.coords_1,
                                 configuration=self.coords_2,
                                 box=self.box_dims,
                                 result=self.allocated_array_2D,
                                 backend='serial')

    def time_self_distance_array(self, num_atoms, pbc_type):
        """Benchmark calculation of all distances
        within a single numpy array of coordinates
        using default parameters to self_distance_array.
        """
        distances.self_distance_array(reference=self.coords_1,
                                      box=self.box_dims,
                                      result=None,
                                      backend='serial')

    def time_self_distance_array_pre_allocated(self, num_atoms, pbc_type):
        """Benchmark calculation of all distances
        within a single numpy array of coordinates
        using self_distance_array with preallocated
        result array.
        """
        distances.self_distance_array(reference=self.coords_1,
                                      box=self.box_dims,
                                      result=self.allocated_array_1D,
                                      backend='serial')

    def time_dist(self, num_atoms, box_dims):
        """Benchmark calculation of distances between
        atoms in two atomgroups with no offsets
        to resids.
        """
        distances.dist(A=self.ag1,
                       B=self.ag2,
                       box=self.box_dims,
                       offset=0)

    def time_dist_offsets(self, num_atoms, box_dims):
        """Benchmark calculation of distances between
        atoms in two atomgroups with offsets
        to resids.
        """
        distances.dist(A=self.ag1,
                       B=self.ag2,
                       box=self.box_dims,
                       offset=20)


class ContactsBench(object):
    """Benchmarks for the MDAnalysis.analysis.distances.contact_matrix
    function in both default and sparse settings.
    """
    u = MDAnalysis.Universe(GRO)
    tri_clinic_dimensions = u.dimensions

    params = ([10, 100, 1000], [None, 'orthogonal', 'triclinic'])
    param_names = ['num_atoms', 'pbc_type']

    def setup(self, num_atoms, pbc_type):

        if pbc_type is None:
            self.u = MDAnalysis.Universe(GRO)
            self.box_dims = None

        elif pbc_type == 'orthogonal':
            # Some strange math/type issues. Essentially we create the smallest
            # cubic box which encloses the triclinic unit cell.
            self.u = MDAnalysis.Universe(GRO)
            box_shape = self.u.dimensions
            basis_vectors = mdamath.triclinic_vectors(box_shape)

            # This will calculate the furthest point from the origin on the
            # triclinic unit cell.
            furthest_point = np.sum(basis_vectors, axis=0)

            # Get full size of cube and deal with certain distance functions
            # requiring float32.
            orthogonal_box_size = np.float32(np.max(furthest_point))
            angle = np.float32(90.0)

            # make dummy cubic box, must be np.array, some functions do
            # not accept list
            self.box_dims = np.array([orthogonal_box_size,
                                      orthogonal_box_size,
                                      orthogonal_box_size,
                                      angle, angle, angle])

        elif pbc_type == 'triclinic':
            self.u = MDAnalysis.Universe(GRO)
            print(self.u.dimensions)
            self.box_dims = self.u.dimensions
        else:
            raise Exception("Specify only 'orthogonal', 'triclinic', or None for the box type")

        self.ag1 = self.u.atoms[:num_atoms]
        self.ag2 = self.u.atoms[num_atoms: 2 * num_atoms]
        self.ag3 = self.u.atoms[-num_atoms:]

        np.random.seed(17809)
        self.coords_1 = np.random.random_sample((num_atoms, 3)).astype(np.float32)
        np.random.seed(9008716)
        self.coords_2 = np.random.random_sample((num_atoms, 3)).astype(np.float32)
        self.allocated_array_2D = np.empty((num_atoms, num_atoms),
                                           dtype=np.float64)
        self.array_shape_1D = int(num_atoms * (num_atoms - 1) / 2.)
        self.allocated_array_1D = np.empty(self.array_shape_1D,
                                           dtype=np.float64)

    def time_contact_matrix(self, num_atoms, pbc_type):
        """Benchmark calculation of contacts within
        a single numpy array using the default arguments
        to contact_matrix.
        """
        distances.contact_matrix(coord=self.coords_1,
                                 cutoff=15.0,
                                 returntype='numpy',
                                 box=self.box_dims)

    def time_contact_matrix_sparse(self, num_atoms, box_dims):
        """Benchmark calculation of contacts within
        a single numpy array using the slower reduced
        memory implementation of contact_matrix intended
        for larger systems.
        """
        distances.contact_matrix(coord=self.coords_1,
                                 cutoff=15.0,
                                 returntype='sparse',
                                 box=self.box_dims)
