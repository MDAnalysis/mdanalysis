# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""
ParmEd topology parser
======================

Converts a `ParmEd <https://parmed.github.io/ParmEd/html>`_
:class:`parmed.structure.Structure` into a :class:`MDAnalysis.core.Topology`.


Example
-------

If you want to use an MDAnalysis-written ParmEd structure for simulation
in ParmEd, you need to first read your files with ParmEd to include the
necessary topology parameters. ::

    >>> import parmed as pmd
    >>> import MDAnalysis as mda
    >>> from MDAnalysis.tests.datafiles import PRM7_ala2, RST7_ala2
    >>> prm = pmd.load_file(PRM7_ala2, RST7_ala2)
    >>> prm
    <AmberParm 3026 atoms; 1003 residues; 3025 bonds; PBC (orthogonal); parametrized>

We can then convert this to an MDAnalysis structure, select only the
protein atoms, and then convert it back to ParmEd. ::

    >>> u = mda.Universe(prm)
    >>> u
    <Universe with 3026 atoms>
    >>> prot = u.select_atoms('protein')
    >>> prm_prot = prot.convert_to('PARMED')
    >>> prm_prot
    <Structure 23 atoms; 2 residues; 22 bonds; PBC (orthogonal); parametrized>

From here you can create an OpenMM simulation system and minimize the
energy. ::

    >>> import simtk.openmm as mm
    >>> import simtk.openmm.app as app
    >>> from parmed import unit as u
    >>> system = prm_prot.createSystem(nonbondedMethod=app.NoCutoff,
    ...                                constraints=app.HBonds,
    ...                                implicitSolvent=app.GBn2)
    >>> integrator = mm.LangevinIntegrator(
    ...                         300*u.kelvin,       # Temperature of heat bath
    ...                         1.0/u.picoseconds,  # Friction coefficient
    ...                         2.0*u.femtoseconds, # Time step
    ... )
    >>> sim = app.Simulation(prm_prot.topology, system, integrator)
    >>> sim.context.setPositions(prm_prot.positions)
    >>> sim.minimizeEnergy(maxIterations=500)

Now you can continue on and run a simulation, if you wish.

Classes
-------

.. autoclass:: ParmEdParser
   :members:
   :inherited-members:

"""
import logging
import numpy as np

from .base import TopologyReaderBase, change_squash
from .tables import Z2SYMB
from ..core.topologyattrs import (
    Atomids,
    Atomnames,
    AltLocs,
    ChainIDs,
    Atomtypes,
    Occupancies,
    Tempfactors,
    Elements,
    Masses,
    Charges,
    Resids,
    Resnums,
    Resnames,
    Segids,
    GBScreens,
    SolventRadii,
    NonbondedIndices,
    RMins,
    Epsilons,
    RMin14s,
    Epsilon14s,
    Bonds,
    UreyBradleys,
    Angles,
    Dihedrals,
    Impropers,
    CMaps
)
from ..core.topology import Topology

logger = logging.getLogger("MDAnalysis.topology.ParmEdParser")


def squash_identical(values):
    if len(values) == 1:
        return values[0]
    else:
        return tuple(values)


class ParmEdParser(TopologyReaderBase):
    """
    For ParmEd structures
    """
    format = 'PARMED'

    @staticmethod
    def _format_hint(thing):
        """Can this Parser read object *thing*?

        .. versionadded:: 1.0.0
        """
        try:
            import parmed as pmd
        except ImportError:  # if no parmed, probably not parmed
            return False
        else:
            return isinstance(thing, pmd.Structure)

    def parse(self, **kwargs):
        """Parse PARMED into Topology

        Returns
        -------
        MDAnalysis *Topology* object


        .. versionchanged:: 2.0.0
           Elements are no longer guessed, if the elements present in the
           parmed object are not recoginsed (usually given an atomic mass of 0)
           then they will be assigned an empty string.
        """
        structure = self.filename

        #### === ATOMS === ####
        names = []
        masses = []
        charges = []
        types = []
        atomic_numbers = []
        serials = []

        resnames = []
        resids = []
        chainids = []
        segids = []

        altLocs = []
        bfactors = []
        occupancies = []

        screens = []
        solvent_radii = []
        nonbonded_indices = []

        rmins = []
        epsilons = []
        rmin14s = []
        epsilon14s = []

        for atom in structure.atoms:
            names.append(atom.name)
            masses.append(atom.mass)
            charges.append(atom.charge)
            types.append(atom.type)
            atomic_numbers.append(atom.atomic_number)
            serials.append(atom.number)

            resnames.append(atom.residue.name)
            resids.append(atom.residue.number)
            chainids.append(atom.residue.chain)
            segids.append(atom.residue.segid)

            altLocs.append(atom.altloc)
            bfactors.append(atom.bfactor)
            occupancies.append(atom.occupancy)

            screens.append(atom.screen)
            solvent_radii.append(atom.solvent_radius)
            nonbonded_indices.append(atom.nb_idx)

            rmins.append(atom.rmin)
            epsilons.append(atom.epsilon)
            rmin14s.append(atom.rmin_14)
            epsilon14s.append(atom.epsilon_14)

        attrs = []

        n_atoms = len(names)

        elements = []

        for z, name in zip(atomic_numbers, names):
            try:
                elements.append(Z2SYMB[z])
            except KeyError:
                elements.append('')

        # Make Atom TopologyAttrs
        for vals, Attr, dtype in (
                (names, Atomnames, object),
                (masses, Masses, np.float32),
                (charges, Charges, np.float32),
                (types, Atomtypes, object),
                (elements, Elements, object),
                (serials, Atomids, np.int32),
                (chainids, ChainIDs, object),

                (altLocs, AltLocs, object),
                (bfactors, Tempfactors, np.float32),
                (occupancies, Occupancies, np.float32),

                (screens, GBScreens, np.float32),
                (solvent_radii, SolventRadii, np.float32),
                (nonbonded_indices, NonbondedIndices, np.int32),

                (rmins, RMins, np.float32),
                (epsilons, Epsilons, np.float32),
                (rmin14s, RMin14s, np.float32),
                (epsilon14s, Epsilon14s, np.float32),
        ):
            attrs.append(Attr(np.array(vals, dtype=dtype)))

        resids = np.array(resids, dtype=np.int32)
        resnames = np.array(resnames, dtype=object)
        chainids = np.array(chainids, dtype=object)
        segids = np.array(segids, dtype=object)

        residx, (resids, resnames, chainids, segids) = change_squash(
            (resids, resnames, chainids, segids),
            (resids, resnames, chainids, segids))

        n_residues = len(resids)
        attrs.append(Resids(resids))
        attrs.append(Resnums(resids.copy()))
        attrs.append(Resnames(resnames))

        segidx, (segids,) = change_squash((segids,), (segids,))
        n_segments = len(segids)
        attrs.append(Segids(segids))

        #### === OTHERS === ####
        bond_values = {}
        bond_types = []
        bond_orders = []

        ub_values = {}
        ub_types = []

        angle_values = {}
        angle_types = []

        dihedral_values = {}
        dihedral_types = []

        improper_values = {}
        improper_types = []

        cmap_values = {}
        cmap_types = []

        for bond in structure.bonds:
            idx = (bond.atom1.idx, bond.atom2.idx)
            if idx not in bond_values:
                bond_values[idx] = ([bond], [bond.order])
            else:
                bond_values[idx][0].append(bond)
                bond_values[idx][1].append(bond.order)

        try:
            bond_values, values = zip(*list(bond_values.items()))
        except ValueError:
            bond_values, bond_types, bond_orders = [], [], []
        else:
            bond_types, bond_orders = zip(*values)

            bond_types = list(map(squash_identical, bond_types))
            bond_orders = list(map(squash_identical, bond_orders))

        attrs.append(Bonds(bond_values, types=bond_types, guessed=False,
                           order=bond_orders))

        for pmdlist, na, values, types in (
            (structure.urey_bradleys, 2, ub_values, ub_types),
            (structure.angles, 3, angle_values, angle_types),
            (structure.dihedrals, 4, dihedral_values, dihedral_types),
            (structure.impropers, 4, improper_values, improper_types),
            (structure.cmaps, 5, cmap_values, cmap_types),
        ):

            for p in pmdlist:
                atoms = ['atom{}'.format(i) for i in range(1, na+1)]
                idx = tuple(getattr(p, a).idx for a in atoms)
                if idx not in values:
                    values[idx] = [p]
                else:
                    values[idx].append(p)

        for dct, Attr in (
            (ub_values, UreyBradleys),
            (angle_values, Angles),
            (dihedral_values, Dihedrals),
            (improper_values, Impropers),
            (cmap_values, CMaps),
        ):
            try:
                vals, types = zip(*list(dct.items()))
            except ValueError:
                vals, types = [], []

            types = list(map(squash_identical, types))
            attrs.append(Attr(vals, types=types, guessed=False, order=None))

        top = Topology(n_atoms, n_residues, n_segments,
                       attrs=attrs,
                       atom_resindex=residx,
                       residue_segindex=segidx)

        return top
