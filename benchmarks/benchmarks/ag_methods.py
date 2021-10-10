import MDAnalysis
import numpy as np

try:
    from MDAnalysisTests.datafiles import (GRO, TPR, XTC,
                                           PSF, DCD,
                                           TRZ_psf, TRZ)
    from MDAnalysis.exceptions import NoDataError
except:
    pass

class AtomGroupMethodsBench(object):
    """Benchmarks for the various MDAnalysis
    atomgroup methods.
    """
    # NOTE: the write() method has been
    # excluded as file writing is considered
    # a separate benchmarking category

    params = (10, 100, 1000, 10000)
    param_names = ['num_atoms']

    def setup(self, num_atoms):
        self.u = MDAnalysis.Universe(GRO)
        self.ag = self.u.atoms[:num_atoms]
        self.weights = np.ones(num_atoms)
        self.vdwradii = {'H':1.0,
                         'C':1.0,
                         'N':1.0,
                         'O':1.0,
                         'DUMMY':1.0}
        self.rot_matrix = np.ones((3,3))
        self.trans = np.ones((4,4))

    def time_bbox_pbc(self, num_atoms):
        """Benchmark bounding box calculation
        with pbc active.
        """
        self.ag.bbox(pbc=True)

    def time_bbox_no_pbc(self, num_atoms):
        """Benchmark bounding box calculation
        with pbc inactive.
        """
        self.ag.bbox(pbc=False)

    def time_bsphere_pbc(self, num_atoms):
        """Benchmark bounding sphere calculation
        with pbc active.
        """
        self.ag.bsphere(pbc=True)

    def time_bsphere_no_pbc(self, num_atoms):
        """Benchmark bounding sphere calculation
        with pbc inactive.
        """
        self.ag.bsphere(pbc=False)

    def time_center_pbc(self, num_atoms):
        """Benchmark center calculation with
        pbc active.
        """
        self.ag.center(weights=self.weights,
                       pbc=True)

    def time_center_no_pbc(self, num_atoms):
        """Benchmark center calculation with
        pbc inactive.
        """
        self.ag.center(weights=self.weights,
                       pbc=False)

    def time_centroid_pbc(self, num_atoms):
        """Benchmark centroid calculation with
        pbc active.
        """
        self.ag.centroid(pbc=True)

    def time_centroid_no_pbc(self, num_atoms):
        """Benchmark centroid calculation with
        pbc inactive.
        """
        self.ag.centroid(pbc=False)

    def time_concatenate(self, num_atoms):
        """Benchmark atomgroup concatenation.
        """
        self.ag.concatenate(self.ag)

    def time_difference(self, num_atoms):
        """Benchmark atomgroup difference
        operation.
        """
        self.ag.difference(self.ag)

    def time_groupby(self, num_atoms):
        """Benchmark atomgroup groupby
        operation.
        """
        self.ag.groupby('resnames')

    def time_guess_bonds(self, num_atoms):
        """Benchmark atomgroup bond guessing
        with artificially-seeded vdw values.
        """
        self.ag.guess_bonds(self.vdwradii)

    def time_intersection(self, num_atoms):
        """Benchmark ag intersection.
        """
        self.ag.intersection(self.ag)

    def time_is_strict_subset(self, num_atoms):
        """Benchmark ag strict subset operation.
        """
        self.ag.is_strict_subset(self.ag)

    def time_is_strict_superset(self, num_atoms):
        """Benchmark ag strict superset operation.
        """
        self.ag.is_strict_superset(self.ag)

    def time_isdisjoint(self, num_atoms):
        """Benchmark disjoint operation between
        atomgroups.
        """
        self.ag.isdisjoint(self.ag)

    def time_issubset(self, num_atoms):
        """Benchmark subset operation between
        atomgroups.
        """
        self.ag.issubset(self.ag)

    def time_issuperset(self, num_atoms):
        """Benchmark superset operation between
        atomgroups.
        """
        self.ag.issuperset(self.ag)

    def time_pack_into_box(self, num_atoms):
        """Benchmark shifting atoms of ag
        into primary unit cell, using
        default parameters.
        """
        self.ag.pack_into_box()

    def time_rotate(self, num_atoms):
        """Benchmark simple rotation operation
        on atomgroup.
        """
        self.ag.rotate(self.rot_matrix)

    def time_rotateby(self, num_atoms):
        """Benchmark rotation by an angle
        of the ag coordinates.
        """
        self.ag.rotateby(angle=45,
                         axis=[1,0,0])

    def time_split(self, num_atoms):
        """Benchmark ag splitting into
        multiple ags based on a simple
        criterion.
        """
        self.ag.split('residue')

    def time_subtract(self, num_atoms):
        """Benchmark ag subtraction.
        """
        self.ag.subtract(self.ag)

    def time_symmetric_difference(self, num_atoms):
        """Benchmark ag symmetric difference
        operation.
        """
        self.ag.symmetric_difference(self.ag)

    def time_transform(self, num_atoms):
        """Benchmark application of transformation
        matrix to atomgroup.
        """
        self.ag.transform(self.trans)

    def time_translate(self, num_atoms):
        """Benchmark the application of a
        translation vector to the ag
        coordinates.
        """
        self.ag.translate([0,0.5,1])

    def time_union(self, num_atoms):
        """Benchmark union operation
        on atomgroups.
        """
        self.ag.union(self.ag)

    def time_wrap(self, num_atoms):
        """Benchmark wrap() operation on
        atomgroup with default params.
        """
        self.ag.wrap()



class AtomGroupAttrsBench(object):
    """Benchmarks for the various MDAnalysis
    atomgroup attributes.
    """

    params = (10, 100, 1000, 10000)
    param_names = ['num_atoms']

    def setup(self, num_atoms):
        self.u = MDAnalysis.Universe(GRO)
        self.ag = self.u.atoms[:num_atoms]

    def time_angle(self, num_atoms):
        """Benchmark simple angle
        calculation. Requires ag
        with three atoms.
        """
        self.ag[:3].angle

    def time_atoms(self, num_atoms):
        """Benchmark returning of identical
        atomgroup.
        """
        self.ag.atoms

    def time_dihedral(self, num_atoms):
        """Benchmark Dihedral object
        creation time. Requires ag of
        size 4.
        """
        self.ag[:4].dihedral

    #TODO: use universe / ag that
    # is suitable for force calc
    def time_forces(self, num_atoms):
        """Benchmark atomgroup force
        calculation.
        """
        try:
            self.ag.forces
        except NoDataError:
            pass

    #TODO: use universe / ag that
    # is suitable for velocity extraction
    def time_velocity(self, num_atoms):
        """Benchmark atomgroup velocity
        values return.
        """
        try:
            self.ag.velocities
        except NoDataError:
            pass

    def time_improper(self, num_atoms):
        """Benchmark improper dihedral
        calculation. Requires ag of size
        4.
        """
        self.ag[:4].improper

    def time_indices(self, num_atoms):
        """Benchmark atom index calculation.
        """
        self.ag.ix

    def time_atomcount(self, num_atoms):
        """Benchmark counting of atoms in
        atomgroup.
        """
        self.ag.n_atoms

    def time_residuecount(self, num_atoms):
        """Benchmark counting of residues in
        atomgroup.
        """
        self.ag.n_residues

    def time_segmentcount(self, num_atoms):
        """Benchmark counting of segments in
        atomgroup.
        """
        self.ag.n_segments

    def time_positions(self, num_atoms):
        """Benchmark returning the positions
        of the atoms in the group.
        """
        self.ag.positions

    def time_residues(self, num_atoms):
        """Benchmark creation of the ResidueGroup
        from the AtomGroup.
        """
        self.ag.residues

    def time_segments(self, num_atoms):
        """Benchmark determination of sorted
        SegmentGroup from AtomGroup.
        """
        self.ag.segments

    def time_ts(self, num_atoms):
        """Benchmark returning of a timestep
        instance from atomgroup.
        """
        self.ag.ts

    def time_unique(self, num_atoms):
        """Benchmark determination of unique
        elements in atomgroup.
        """
        self.ag.unique

    def time_bond(self, num_atoms):
        """Benchmark Bond object creation.
        Requires ag of size 2.
        """
        self.ag[:2].bond


class CompoundSplitting(object):
    """Test how fast can we split compounds into masks and apply them
    
    The benchmark used in Issue #3000. Parameterizes multiple compound number
    and size combinations.
    """
    
    params = [(100, 10000, 1000000),  # n_atoms
              (1, 10, 100),           # n_compounds
              (True, False),          # homogeneous
              (True, False)]          # contiguous

    def setup(self, n_atoms, n_compounds, homogeneous, contiguous):
        rg = np.random.Generator(np.random.MT19937(3000))

        # some parameter screening for nonsensical combinations.
        if n_compounds > n_atoms:
            raise NotImplementedError

        if n_compounds == 1 and not (homogeneous and contiguous):
            raise NotImplementedError
            
        if n_compounds == n_atoms:
            if not (homogeneous and contiguous):
                raise NotImplementedError
            compound_indices = np.arange(n_compounds)
        elif homogeneous:
            ats_per_compound, remainder = divmod(n_atoms, n_compounds)
            if remainder:
                raise NotImplementedError
            compound_indices = np.tile(np.arange(n_compounds),
                                       (ats_per_compound, 1)).T.ravel()
        else:
            compound_indices = np.sort(np.floor(rg.random(n_atoms)
                                                * n_compounds).astype(np.int))
                
        unique_indices = np.unique(compound_indices)
        if len(unique_indices) != n_compounds:
            raise RuntimeError
        
        if not contiguous:
            rg.shuffle(compound_indices)
        
        self.u = MDAnalysis.Universe.empty(n_atoms,
                                           n_residues=n_compounds,
                                           n_segments=1,
                                           atom_resindex=compound_indices,
                                           trajectory=True)
        self.u.atoms.positions = rg.random((n_atoms, 3),
                                           dtype=np.float32) * 100
        self.u.dimensions = [50, 50, 50, 90, 90, 90]

    def time_center_compounds(self, *args):
        self.u.atoms.center(None, compound='residues')


class FragmentFinding(object):
    """Test how quickly we find fragments (distinct molecules from bonds)"""
    # if we try to parametrize over topology &
    # trajectory formats asv will use all
    # possible combinations, so instead handle
    # this in setup()
    params = ('large_fragment_small_solvents',
              'large_fragment',
              'polymer_chains', # 20ish polymer chains
              )
    param_names = ['universe_type']

    def setup(self, universe_type):
        if universe_type == 'large_fragment_small_solvents':
            univ = (TPR, XTC) 
        elif universe_type == 'large_fragment':
            univ = (PSF, DCD)
        else:
            univ = (TRZ_psf, TRZ) 
        self.u = MDAnalysis.Universe(*univ)

    def time_find_fragments(self, universe_type):
        frags = self.u.atoms.fragments


class FragmentCaching(FragmentFinding):
    """Test how quickly we find cached fragments"""
    def setup(self, universe_type):
        super(FragmentCaching, self).setup(universe_type)
        frags = self.u.atoms.fragments  # Priming the cache

    def time_find_cached_fragments(self, universe_type):
        frags = self.u.atoms.fragments
