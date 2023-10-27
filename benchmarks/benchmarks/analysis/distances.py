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
    functions. Excluding contact matrices.
    """
    timeout = 180
    params = ([10, 100, 1000, 10000])
    param_names = ['num_atoms']

    def setup(self, num_atoms):

        # Some functions are performance dependent on a realistic
        # atom density. So we make a new universe with random atom
        # positions, density derived from an example universe.

        # Make fake universe with 3*num_atoms because some of our benchmarks
        # require 3 atom groups.
        universe_num_atoms = num_atoms*3
        self.u = MDAnalysis.Universe.empty(universe_num_atoms, trajectory=True)
        # calculating reasonable density for an atomistic system.
        ex_universe = MDAnalysis.Universe(GRO)
        density = ex_universe.coord.volume/ex_universe.atoms.n_atoms
        desired_box_volume = universe_num_atoms*density
        # making a cube with orthogonal basis vectors
        cube_side_length = np.cbrt(desired_box_volume)
        sides = [float(cube_side_length)]*3
        ortho_angles = [float(90)]*3
        self.u.dimensions = np.array(sides + ortho_angles)

        # filling the cube with random points.
        np.random.seed(17809)
        random_samples = np.random.random_sample((universe_num_atoms,
                                                  3)).astype(np.float32)
        self.u.atoms.positions = cube_side_length * random_samples

        # dealing with missing topology information from empty universe
        # avoids errors
        self.u.add_TopologyAttr('resid', [1])
        self.u.add_TopologyAttr('segid', ["AP1"])
        # splitting atom groups into 2 groups with the same number of atoms
        self.ag1 = self.u.atoms[:num_atoms]
        self.ag2 = self.u.atoms[num_atoms:-num_atoms]
        self.ag3 = self.u.atoms[-num_atoms:]

        np.random.seed(17809)
        self.coords_1 = np.random.random_sample((num_atoms,
                                                 3)).astype(np.float32)
        np.random.seed(9008716)
        self.coords_2 = np.random.random_sample((num_atoms,
                                                 3)).astype(np.float32)
        self.coords_3 = np.random.random_sample((num_atoms,
                                                 3)).astype(np.float32)
        self.coords_4 = np.random.random_sample((num_atoms,
                                                 3)).astype(np.float32)
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

    def time_dist(self, num_atoms):
        """Benchmark calculation of distances between
        atoms in two atomgroups with no offsets
        to resids.
        """
        distances.dist(A=self.ag1,
                       B=self.ag2,
                       box=None,
                       offset=0)

    def time_dist_offsets(self, num_atoms):
        """Benchmark calculation of distances between
        atoms in two atomgroups with offsets
        to resids.
        """
        distances.dist(A=self.ag1,
                       B=self.ag2,
                       box=None,
                       offset=20)

    def time_contact_matrix(self, num_atoms):
        """Benchmark calculation of contacts within
        a single numpy array using the default arguments
        to contact_matrix.
        """
        distances.contact_matrix(
            coord=self.ag1.positions, cutoff=15.0, returntype="numpy", box=None
        )

    def time_contact_matrix_sparse(self, num_atoms):
        """Benchmark calculation of contacts within
        a single numpy array using the slower reduced
        memory implementation of contact_matrix intended
        for larger systems.
        """
        distances.contact_matrix(coord=self.ag1.positions,
                                 cutoff=15.0,
                                 returntype='sparse',
                                 box=None)

    def time_calc_bonds(self, num_atoms):
        """Benchmark calculation of bonds between
        atoms in two atomgroups.
        """
        MDAnalysis.lib.distances.calc_bonds(self.coords_1, self.coords_2)

    def time_calc_angles(self, num_atoms):
        """Benchmark calculation of angles between
        atoms in three atomgroups.
        """
        MDAnalysis.lib.distances.calc_angles(
            self.coords_1, self.coords_2, self.coords_3
        )

    def time_calc_dihedrals(self, num_atoms):
        """Benchmark calculation of dihedrals between
        atoms in four atomgroups.
        """
        MDAnalysis.lib.distances.calc_dihedrals(
            self.coords_1, self.coords_2, self.coords_3, self.coords_4
        )


class DistancesBenchPBC(object):
    """Benchmarks for MDAnalysis.analysis.distances
    functions.
    """
    timeout = 180
    params = ([10, 100, 1000, 10000], ['orthogonal', 'triclinic'])
    param_names = ['num_atoms', 'pbc_type']


    def setup(self, num_atoms, pbc_type):

        # This setup function is virtually identical to the one
        # in the DistancesBench class above, we just have to do some
        # extra work to set up different periodic boundaries.

        # Some functions are performance dependent on a realistic
        # atom density. So we make a new universe with random atom
        # positions, density derived from an example universe.

        # Make fake universe with 3*num_atoms because some of our benchmarks
        # require 3 atom groups.
        universe_num_atoms = num_atoms*3
        self.u = MDAnalysis.Universe.empty(universe_num_atoms, trajectory=True)
        # calculating reasonable density for an atomistic system.
        ex_universe = MDAnalysis.Universe(GRO)
        density = ex_universe.coord.volume/ex_universe.atoms.n_atoms
        desired_box_volume = universe_num_atoms*density
        # making a cube with orthogonal basis vectors
        cube_side_length = np.cbrt(desired_box_volume)
        sides = [float(cube_side_length)]*3
        ortho_angles = [float(90)]*3
        self.u.dimensions = np.array(sides + ortho_angles)

        # filling the cube with random points.
        np.random.seed(17809)
        random_samples = np.random.random_sample((universe_num_atoms,
                                                  3)).astype(np.float32)
        self.u.atoms.positions = cube_side_length * random_samples

        if pbc_type is None:
            self.box_dims = None

        elif pbc_type == 'orthogonal':
            self.box_dims = self.u.dimensions

        elif pbc_type == 'triclinic':
            # making a triclinic box with the same volume as the cube
            alpha = 2*np.pi/3
            beta = np.pi/2
            gamma = np.pi/3

            # change side lengths so that the resulting box has correct
            # volume
            a = cube_side_length
            b = cube_side_length/np.sin(gamma)
            c = np.sqrt(cube_side_length**2/(1 -
                                             ((np.cos(alpha)
                                              - np.cos(beta)
                                              * np.cos(gamma))
                                              / np.sin(gamma))**2
                                             - np.cos(beta)**2))

            self.u.dimensions = [a, b, c,
                                 np.rad2deg(alpha),
                                 np.rad2deg(beta),
                                 np.rad2deg(gamma)]

            self.box_dims = self.u.dimensions
            # wrapping atoms to reflect new triclinic basis
            self.u.atoms.wrap(inplace=True)

        # dealing with missing topology information from empty universe
        # avoids errors
        self.u.add_TopologyAttr('resid', [1])
        self.u.add_TopologyAttr('segid', ["AP1"])
        # splitting atom groups into 2 groups with the same number of atoms
        self.ag1 = self.u.atoms[:num_atoms]
        self.ag2 = self.u.atoms[num_atoms:-num_atoms]
        self.ag3 = self.u.atoms[-num_atoms:]

        np.random.seed(17809)
        self.coords_1 = np.random.random_sample((num_atoms,
                                                 3)).astype(np.float32)
        np.random.seed(9008716)
        self.coords_2 = np.random.random_sample((num_atoms,
                                                 3)).astype(np.float32)
        self.allocated_array_2D = np.empty((num_atoms, num_atoms),
                                           dtype=np.float64)
        self.array_shape_1D = int(num_atoms * (num_atoms - 1) / 2.)
        self.allocated_array_1D = np.empty(self.array_shape_1D,
                                           dtype=np.float64)

    def time_distance_array_pbc(self, num_atoms, pbc_type):
        """Benchmark calculation of all distances
        between two numpy arrays of coordinates,
        using default arguments to distance_array.
        """
        distances.distance_array(reference=self.coords_1,
                                 configuration=self.coords_2,
                                 box=self.box_dims,
                                 result=None,
                                 backend='serial')

    def time_distance_array_pre_allocated_pbc(self, num_atoms, pbc_type):
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

    def time_self_distance_array_pbc(self, num_atoms, pbc_type):
        """Benchmark calculation of all distances
        within a single numpy array of coordinates
        using default parameters to self_distance_array.
        """
        distances.self_distance_array(reference=self.coords_1,
                                      box=self.box_dims,
                                      result=None,
                                      backend='serial')

    def time_self_distance_array_pre_allocated_pbc(self, num_atoms, pbc_type):
        """Benchmark calculation of all distances
        within a single numpy array of coordinates
        using self_distance_array with preallocated
        result array.
        """
        distances.self_distance_array(reference=self.coords_1,
                                      box=self.box_dims,
                                      result=self.allocated_array_1D,
                                      backend='serial')

    def time_dist_pbc(self, num_atoms, pbc_type):
        """Benchmark calculation of distances between
        atoms in two atomgroups with no offsets
        to resids.
        """
        distances.dist(A=self.ag1,
                       B=self.ag2,
                       box=self.box_dims,
                       offset=0)


    def time_dist_offsets_pbc(self, num_atoms, pbc_type):
        """Benchmark calculation of distances between
        atoms in two atomgroups with offsets
        to resids.
        """
        distances.dist(A=self.ag1,
                       B=self.ag2,
                       box=self.box_dims,
                       offset=20)


    def time_contact_matrix_pbc(self, num_atoms, pbc_type):
        """Benchmark calculation of contacts within
        a single numpy array using the default arguments
        to contact_matrix.
        """
        distances.contact_matrix(coord=self.ag1.positions,
                                 cutoff=15.0,
                                 returntype='numpy',
                                 box=self.box_dims)

    def time_contact_matrix_sparse_pbc(self, num_atoms, pbc_type):
        """Benchmark calculation of contacts within
        a single numpy array using the slower reduced
        memory implementation of contact_matrix intended
        for larger systems.
        """
        distances.contact_matrix(coord=self.ag1.positions,
                                 cutoff=15.0,
                                 returntype='sparse',
                                 box=self.box_dims)

