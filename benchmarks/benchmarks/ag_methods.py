from __future__ import division, absolute_import, print_function

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


class FragmentFinding(object):
    """Test how quickly we find fragments (distinct molecules from bonds)"""
    params = [(TPR, XTC),  # single large fragment, many small solvents
              (PSF, DCD),  # single large fragment
              (TRZ_psf, TRZ)]  # 20ish polymer chains
    param_names = ['universe']

    def setup(self, universe):
        self.u = MDAnalysis.Universe(*universe)

    def time_find_fragments(self, universe):
        frags = self.u.atoms.fragments
