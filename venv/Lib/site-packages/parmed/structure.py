"""
This module contains the core base class for all of the chemical structures with
various topological and force field features.
"""
import logging
import math
import os
from collections import defaultdict
from copy import copy, deepcopy
import warnings

import numpy as np

from . import unit as u
from . import residue
from .constants import DEG_TO_RAD, SMALL
from .exceptions import ParameterError, ParmedWarning
from .geometry import (STANDARD_BOND_LENGTHS_SQUARED, box_lengths_and_angles_to_vectors,
                       box_vectors_to_lengths_and_angles, distance2)
from .topologyobjects import (AcceptorDonor, Angle, Atom, AtomList, Bond, ChiralFrame, Cmap,
                              Dihedral, DihedralType, DihedralTypeList, DrudeAtom, ExtraPoint, Group, Improper,
                              MultipoleFrame, NonbondedException, NoUreyBradley, OutOfPlaneBend,
                              OutOfPlaneExtraPointFrame, PiTorsion, ResidueList, StretchBend,
                              ThreeParticleExtraPointFrame, TorsionTorsion, TrackedList,
                              TrigonalAngle, TwoParticleExtraPointFrame, UnassignedAtomType,
                              UreyBradley, Link, LocalCoordinatesFrame, QualitativeBondType)
from .utils import PYPY, find_atom_pairs, tag_molecules
from .utils.decorators import needs_openmm
from .vec3 import Vec3
from .rdkit import RDKit

# Try to import the OpenMM modules
try:
    from openmm import app
    import openmm as mm
    from openmm.app.internal.unitcell import reducePeriodicBoxVectors
except ImportError:
    mm = app = None

LOGGER = logging.getLogger(__name__)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Private attributes and methods


def _strip_box_units(args):
    new_args = []
    for arg in args:
        # Handle 3 types of arguments here: units, regular numbers, and
        # iterables. Iterables have units removed recursively. But protect
        # against strings, since you get infinite recursion that way
        if u.is_quantity(arg):
            if arg.unit.is_compatible(u.degree):
                new_args.append(arg.value_in_unit(u.degree))
            else:
                new_args.append(arg.value_in_unit(u.angstroms))
        elif isinstance(arg, str):
            raise TypeError('Unit cell cannot have strings')
        else:
            try:
                iter(arg)
                arg = _strip_box_units(arg)
            except TypeError:
                pass
            new_args.append(arg)
    return new_args


def _bondi(atom):
    if atom.atomic_number == 6:
        return 1.7
    if atom.atomic_number == 1:
        return 1.2
    if atom.atomic_number == 7:
        return 1.55
    if atom.atomic_number == 14:
        return 2.1
    if atom.atomic_number == 15:
        return 1.85
    if atom.atomic_number == 16:
        return 1.8
    return 1.5


def _mbondi(atom):
    if atom.atomic_number == 1:
        bondeds = atom.bond_partners
        if not bondeds or bondeds[0].atomic_number in (6, 7):
            return 1.3
        if bondeds[0].atomic_number in (8, 16):
            return 0.8
        return 1.2
    return _bondi(atom)


def _mbondi2(atom):
    if atom.atomic_number == 1:
        if not atom.bond_partners or atom.bond_partners[0].atomic_number == 7:
            return 1.3
        return 1.2
    return _bondi(atom)


def _mbondi3(atom):
    if atom.residue.name in ('GLU', 'ASP', 'GL4', 'AS4'):
        if atom.name.startswith('OE') or atom.name.startswith('OD'):
            return 1.4
    elif atom.residue.name == 'ARG':
        if atom.name.startswith('HH') or atom.name.startswith('HE'):
            return 1.17
    if atom.name == 'OXT':
        return 1.4
    return _mbondi2(atom)


@needs_openmm
def _gb_rad_screen(atom, model):
    """
    Gets the default GB parameters for a given atom according to a specific
    Generalized Born model

    Parameters
    ----------
    atom : :class:`Atom`
        The atom to get the default GB parameters for
    model : ``app.HCT, app.OBC1, or app.OBC2``
        The GB model to get the default parameters for (app.GBn and app.GBn2 are
        already handled in Structure._get_gb_parameters)

    Returns
    -------
    radius, screen [,alpha, beta, gamma] : ``float, float [,float, float, float]``
        The intrinsic radius of the atom and the screening factor of the atom.
        If the model is GBn2, alpha, beta, and gamma parameters are also
        returned
    """
    if model in (app.OBC1, app.OBC2):
        rad = _mbondi2(atom)
    else:
        rad = _mbondi(atom)
    if atom.atomic_number == 1:
        return rad, 0.85
    if atom.atomic_number == 6:
        return rad, 0.72
    if atom.atomic_number == 7:
        return rad, 0.79
    if atom.atomic_number == 8:
        return rad, 0.85
    if atom.atomic_number == 9:
        return rad, 0.88
    if atom.atomic_number == 15:
        return rad, 0.86
    if atom.atomic_number == 16:
        return rad, 0.96
    return rad, 0.8

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


class Structure:
    """
    A chemical structure composed of atoms, bonds, angles, torsions, and other
    topological features

    Attributes
    ----------
    atoms : :class:`AtomList`
        List of all atoms in the structure
    residues : :class:`ResidueList`
        List of all residues in the structure
    bonds : :class:`TrackedList` (:class:`Bond`)
        List of all bonds in the structure
    angles : :class:`TrackedList` (:class:`Angle`)
        List of all angles in the structure
    dihedrals : :class:`TrackedList` (:class:`Dihedral`)
        List of all dihedrals in the structure
    rb_torsions : :class:`TrackedList` (:class:`Dihedral`)
        List of all Ryckaert-Bellemans torsions in the structure
    urey_bradleys : :class:`TrackedList` (:class:`UreyBradley`)
        List of all Urey-Bradley angle bends in the structure
    impropers : :class:`TrackedList` (:class:`Improper`)
        List of all CHARMM-style improper torsions in the structure
    cmaps : :class:`TrackedList` (:class:`Cmap`)
        List of all CMAP objects in the structure
    trigonal_angles : :class:`TrackedList` (:class:`TrigonalAngle`)
        List of all AMOEBA-style trigonal angles in the structure
    out_of_plane_bends : :class:`TrackedList` (:class:`OutOfPlaneBends`)
        List of all AMOEBA-style out-of-plane bending angles
    pi_torsions : :class:`TrackedList` (:class:`PiTorsion`)
        List of all AMOEBA-style pi-torsion angles
    stretch_bends : :class:`TrackedList` (:class:`StretchBend`)
        List of all AMOEBA-style stretch-bend compound bond/angle terms
    torsion_torsions : :class:`TrackedList` (:class:`TorsionTorsion`)
        List of all AMOEBA-style coupled torsion-torsion terms
    chiral_frames : :class:`TrackedList` (:class:`ChiralFrame`)
        List of all AMOEBA-style chiral frames defined in the structure
    multipole_frames : :class:`TrackedList` (:class:`MultipoleFrame`)
        List of all AMOEBA-style multipole frames defined in the structure
    adjusts : :class:`TrackedList` (:class:`NonbondedException`)
        List of all nonbonded pair-exception rules
    acceptors : :class:`TrackedList` (:class:`AcceptorDonor`)
        List of all H-bond acceptors, if that information is present
    donors : :class:`TrackedList` (:class:`AcceptorDonor`)
        List of all H-bond donors, if that information is present
    groups : :class:`TrackedList` (:class:`Group`)
        List of all CHARMM-style GROUP objects (whatever those are used for)
    box : ``list of 6 floats``
        Box dimensions (a, b, c, alpha, beta, gamma) for the unit cell. If no
        box is defined, `box` is set to `None`
    space_group : ``str``
        The space group of the structure (default is "P 1")
    nrexcl : ``int``
        The number of bonds away that an atom can be in order to be excluded
        from the direct nonbonded summation
    title : ``str``
        Cosmetic only, it is an arbitrary title assigned to the system. Default
        value is an empty string
    positions : u.Quantity(list(Vec3), u.angstroms)
        Unit-bearing atomic coordinates. If not all atoms have coordinates, this
        property is None
    coordinates : np.ndarray of shape (nframes, natom, 3)
        If no coordinates are set, this is set to None. The first frame will
        match the coordinates present on the atoms.
    symmetry : :class:`Symmetry`
        if no symmetry is set, this is set to None.
    links : :class:`TrackedList` (:class:`Link`)
        The list of Link definitions for this Structure
    rdkit_mol : :class:`rdkit.Chem.Mol`
        An RDKit Mol object from the current structure.

    Notes
    -----
    This class also has a handful of type lists for each of the attributes above
    (excluding `atoms`, `residues`, `chiral_frames`, and `multipole_frames`).
    They are all TrackedList instances that are designed to hold the relevant
    parameter type. The list is:
        bond_types, angle_types, dihedral_types, urey_bradley_types,
        improper_types, cmap_types, trigonal_angle_types,
        out_of_plane_bend_types, pi_torsion_types, stretch_bend_types,
        torsion_torsion_types, adjust_types

    dihedral_types _may_ be a list of :class:`DihedralType` instances, since
    torsion profiles are often represented by a Fourier series with multiple
    terms
    """
    # Force groups assigned to each type of force
    BOND_FORCE_GROUP = 0
    ANGLE_FORCE_GROUP = 1
    DIHEDRAL_FORCE_GROUP = 2
    UREY_BRADLEY_FORCE_GROUP = 3
    IMPROPER_FORCE_GROUP = 4
    CMAP_FORCE_GROUP = 5
    TRIGONAL_ANGLE_FORCE_GROUP = 6
    OUT_OF_PLANE_BEND_FORCE_GROUP = 7
    PI_TORSION_FORCE_GROUP = 8
    STRETCH_BEND_FORCE_GROUP = 9
    TORSION_TORSION_FORCE_GROUP = 10
    NONBONDED_FORCE_GROUP = 11
    RB_TORSION_FORCE_GROUP = 12
    DRUDE_FORCE_GROUP = 13

    # ===================================================

    def __init__(self):

        # Topological object lists
        self.atoms = AtomList()
        self.residues = ResidueList()
        self.bonds = TrackedList()
        self.angles = TrackedList()
        self.dihedrals = TrackedList()
        self.rb_torsions = TrackedList()
        self.urey_bradleys = TrackedList()
        self.impropers = TrackedList()
        self.cmaps = TrackedList()
        self.trigonal_angles = TrackedList()
        self.out_of_plane_bends = TrackedList()
        self.pi_torsions = TrackedList()
        self.stretch_bends = TrackedList()
        self.torsion_torsions = TrackedList()
        self.chiral_frames = TrackedList()
        self.multipole_frames = TrackedList()
        self.adjusts = TrackedList()
        # Extraneous information stored in CHARMM PSF files... not used as far
        # as I can tell for much of anything
        self.acceptors = TrackedList()
        self.donors = TrackedList()
        self.groups = TrackedList()

        # Parameter type lists
        self.bond_types = TrackedList()
        self.angle_types = TrackedList()
        self.dihedral_types = TrackedList()
        self.urey_bradley_types = TrackedList()
        self.improper_types = TrackedList()
        self.rb_torsion_types = TrackedList()
        self.cmap_types = TrackedList()
        self.trigonal_angle_types = TrackedList()
        self.out_of_plane_bend_types = TrackedList()
        self.pi_torsion_types = TrackedList()
        self.stretch_bend_types = TrackedList()
        self.torsion_torsion_types = TrackedList()
        self.adjust_types = TrackedList()
        self.links = TrackedList()

        self._box = None
        self._coordinates = None
        self.space_group = "P 1"
        self.unknown_functional = False
        self.nrexcl = 3
        self.title = ''
        self._combining_rule = 'lorentz'
        self.symmetry = None

    # ===================================================

    def __repr__(self):
        natom = len(self.atoms)
        nres = len(self.residues)
        nextra = sum([isinstance(a, ExtraPoint) for a in self.atoms])
        retstr = ['<%s %d atoms' % (type(self).__name__, natom)]
        if nextra > 0:
            retstr.append(' [%d EPs]' % nextra)
        retstr.append('; %d residues' % nres)
        nbond = len(self.bonds)
        retstr.append('; %d bonds' % nbond)
        if self.box is not None:
            retstr.append('; PBC')
            if (abs(self.box[3]-90) > SMALL or abs(self.box[4]-90) > SMALL or
                    abs(self.box[5]-90) > SMALL):
                retstr.append(' (triclinic)')
            else:
                retstr.append(' (orthogonal)')
        # Just assume that if the first atom has a defined type, so does
        # everything else... we don't want __repr__ to be super expensive
        if self.atoms[0].atom_type is UnassignedAtomType and self.atoms[0]._epsilon is None:
            retstr.append('; NOT parameterized>')
        else:
            retstr.append('; parameterized>')
        return ''.join(retstr)

    # ===================================================

    def add_atom(self, atom, resname, resnum, chain='', inscode='', segid=''):
        """
        Adds a new atom to the Structure, adding a new residue to `residues` if
        it has a different name or number as the last residue added and adding
        it to the `atoms` list.

        Parameters
        ----------
        atom : :class:`Atom`
            The atom to add to this residue list
        resname : ``str``
            The name of the residue this atom belongs to
        resnum : ``int``
            The number of the residue this atom belongs to
        chain : ``str``
            The chain ID character for this residue
        inscode : ``str``
            The insertion code ID character for this residue (it is stripped)
        segid : ``str``
            The segment identifier for this residue (it is stripped)

        Notes
        -----
        If the residue name and number differ from the last residue in this
        list, a new residue is added and the atom is added to that residue
        """
        self.residues.add_atom(atom, resname, resnum, chain, inscode, segid)
        self.atoms.append(atom)

    # ===================================================

    def add_atom_to_residue(self, atom, residue):
        """
        Adds a new atom to the Structure at the end if the given residue

        Parameters
        ----------
        atom : :class:`Atom`
            The atom to add to the system
        residue : :class:`Residue`
            The residue to which to add this atom. It MUST be part of this
            Structure instance already or a ValueError is raised

        Notes
        -----
        This atom is added at the end of the residue and is inserted into the
        `atoms` list in such a way that all residues are composed of atoms
        contiguous in the atoms list. For large systems, this may be a
        relatively expensive operation
        """
        # Make sure residue belongs to this list
        if residue.list is not self.residues:
            raise ValueError('Residue is not part of the structure')
        last_atom = residue.atoms[-1]
        residue.add_atom(atom)
        # Special-case if this is going to be the last atom
        if not self.atoms or last_atom is self.atoms[-1]:
            self.atoms.append(atom)
        else:
            self.atoms.insert(last_atom.idx + 1, atom)

    # ===================================================

    def __copy__(self):
        """ A deep copy of the Structure """
        return self.copy(type(self))

    # ===================================================

    def copy(self, cls, split_dihedrals=False):
        """
        Makes a copy of the current structure as an instance of a specified
        subclass

        Parameters
        ----------
        cls : Structure subclass
            The returned object is a copy of this structure as a `cls` instance
        split_dihedrals : ``bool``
            If True, then the Dihedral entries will be split up so that each one
            is paired with a single DihedralType (rather than a
            DihedralTypeList)

        Returns
        -------
        *cls* instance
            The instance of the Structure subclass `cls` with a copy of the
            current Structure's topology information
        """
        c = cls()
        for atom in self.atoms:
            res = atom.residue
            a = deepcopy(atom)
            c.add_atom(a, res.name, res.number, res.chain,
                       res.insertion_code, res.segid)
        # Now copy all of the types
        for bt in self.bond_types:
            c.bond_types.append(deepcopy(bt))
        c.bond_types.claim()
        for at in self.angle_types:
            c.angle_types.append(deepcopy(at))
        c.angle_types.claim()
        if split_dihedrals:
            ndt = 0
            mapdt = {}
            for idt, dt in enumerate(self.dihedral_types):
                if hasattr(dt, '__iter__'):
                    for t in dt:
                        c.dihedral_types.append(deepcopy(t))
                        mapdt.setdefault(idt, []).append(ndt)
                        ndt += 1
                else:
                    mapdt.setdefault(idt, []).append(ndt)
                    ndt += 1
                    c.dihedral_types.append(deepcopy(dt))
        else:
            for dt in self.dihedral_types:
                c.dihedral_types.append(deepcopy(dt))
        c.dihedral_types.claim()
        for ut in self.urey_bradley_types:
            c.urey_bradley_types.append(deepcopy(ut))
        c.urey_bradley_types.claim()
        for it in self.improper_types:
            c.improper_types.append(deepcopy(it))
        c.improper_types.claim()
        for rt in self.rb_torsion_types:
            c.rb_torsion_types.append(deepcopy(rt))
        c.rb_torsion_types.claim()
        for ct in self.cmap_types:
            c.cmap_types.append(deepcopy(ct))
        c.cmap_types.claim()
        for ta in self.trigonal_angle_types:
            c.trigonal_angle_types.append(deepcopy(ta))
        c.trigonal_angle_types.claim()
        for ot in self.out_of_plane_bend_types:
            c.out_of_plane_bend_types.append(deepcopy(ot))
        c.out_of_plane_bend_types.claim()
        for pt in self.pi_torsion_types:
            c.pi_torsion_types.append(deepcopy(pt))
        c.pi_torsion_types.claim()
        for st in self.stretch_bend_types:
            c.stretch_bend_types.append(deepcopy(st))
        c.stretch_bend_types.claim()
        for tt in self.torsion_torsion_types:
            c.torsion_torsion_types.append(deepcopy(tt))
        c.torsion_torsion_types.claim()
        for at in self.adjust_types:
            c.adjust_types.append(deepcopy(at))
        c.adjust_types.claim()
        # Now create the topological objects
        atoms = c.atoms
        for b in self.bonds:
            c.bonds.append(
                Bond(atoms[b.atom1.idx], atoms[b.atom2.idx],
                     b.type and c.bond_types[b.type.idx], order=b.order,
                     qualitative_type=b.qualitative_type)
            )
            c.bonds[-1].funct = b.funct
        for a in self.angles:
            c.angles.append(
                Angle(atoms[a.atom1.idx], atoms[a.atom2.idx],
                      atoms[a.atom3.idx],
                      a.type and c.angle_types[a.type.idx])
            )
            c.angles[-1].funct = a.funct
        if split_dihedrals:
            for d in self.dihedrals:
                if hasattr(d.type, '__iter__'):
                    for i in range(len(d.type)):
                        ie = d.ignore_end or i < len(d.type) - 1
                        ti = mapdt[d.type.idx][i]
                        c.dihedrals.append(
                            Dihedral(atoms[d.atom1.idx],
                                     atoms[d.atom2.idx],
                                     atoms[d.atom3.idx],
                                     atoms[d.atom4.idx],
                                     improper=d.improper, ignore_end=ie,
                                     type=c.dihedral_types[ti])
                        )
                        c.dihedrals[-1]._funct = d._funct
                else:
                    ti = mapdt[d.type.idx][0]
                    c.dihedrals.append(
                        Dihedral(atoms[d.atom1.idx], atoms[d.atom2.idx],
                                 atoms[d.atom3.idx], atoms[d.atom4.idx],
                                 improper=d.improper, ignore_end=d.ignore_end,
                                 type=d.type and c.dihedral_types[ti])
                    )
                    c.dihedrals[-1]._funct = d._funct
        else:
            for d in self.dihedrals:
                if d.type is None:
                    typ = None
                else:
                    typ = c.dihedral_types[d.type.idx]
                c.dihedrals.append(
                    Dihedral(atoms[d.atom1.idx], atoms[d.atom2.idx],
                             atoms[d.atom3.idx], atoms[d.atom4.idx],
                             improper=d.improper, ignore_end=d.ignore_end,
                             type=typ)
                )
                c.dihedrals[-1]._funct = d._funct
        for ub in self.urey_bradleys:
            if ub.type is NoUreyBradley:
                typ = NoUreyBradley
            else:
                typ = ub.type and c.urey_bradley_types[ub.type.idx]
            c.urey_bradleys.append(
                UreyBradley(atoms[ub.atom1.idx], atoms[ub.atom2.idx], typ)
            )
        for i in self.impropers:
            c.impropers.append(
                Improper(atoms[i.atom1.idx], atoms[i.atom2.idx],
                         atoms[i.atom3.idx], atoms[i.atom4.idx],
                         i.type and c.improper_types[i.type.idx])
            )
        for r in self.rb_torsions:
            c.rb_torsions.append(
                Dihedral(atoms[r.atom1.idx], atoms[r.atom2.idx],
                         atoms[r.atom3.idx], atoms[r.atom4.idx],
                         type=r.type and c.rb_torsion_types[r.type.idx])
            )
            c.rb_torsions[-1]._funct = r._funct
        for cm in self.cmaps:
            c.cmaps.append(
                Cmap(atoms[cm.atom1.idx], atoms[cm.atom2.idx],
                     atoms[cm.atom3.idx], atoms[cm.atom4.idx],
                     atoms[cm.atom5.idx],
                     cm.type and c.cmap_types[cm.type.idx])
            )
            c.cmaps[-1].funct = cm.funct
        for t in self.trigonal_angles:
            c.trigonal_angles.append(
                TrigonalAngle(atoms[t.atom1.idx], atoms[t.atom2.idx],
                              atoms[t.atom3.idx], atoms[t.atom4.idx],
                              t.type and c.trigonal_angle_types[t.type.idx])
            )
        for o in self.out_of_plane_bends:
            c.out_of_plane_bends.append(
                OutOfPlaneBend(atoms[o.atom1.idx], atoms[o.atom2.idx],
                               atoms[o.atom3.idx], atoms[o.atom4.idx],
                               o.type and c.out_of_plane_bend_types[o.type.idx])
            )
        for p in self.pi_torsions:
            c.pi_torsions.append(
                PiTorsion(atoms[p.atom1.idx], atoms[p.atom2.idx],
                          atoms[p.atom3.idx], atoms[p.atom4.idx],
                          atoms[p.atom5.idx], atoms[p.atom6.idx],
                          p.type and c.pi_torsion_types[p.type.idx])
            )
        for s in self.stretch_bends:
            c.stretch_bends.append(
                StretchBend(atoms[s.atom1.idx], atoms[s.atom2.idx],
                            atoms[s.atom3.idx],
                            s.type and c.stretch_bend_types[s.type.idx])
            )
        for t in self.torsion_torsions:
            c.torsion_torsions.append(
                TorsionTorsion(atoms[t.atom1.idx], atoms[t.atom2.idx], atoms[t.atom3.idx], atoms[t.atom4.idx],
                               atoms[t.atom5.idx], t.type and c.torsion_torsion_types[t.type.idx])
            )
        for ch in self.chiral_frames:
            c.chiral_frames.append(
                ChiralFrame(atoms[ch.atom1.idx],
                            atoms[ch.atom2.idx], ch.chirality)
            )
        for m in self.multipole_frames:
            c.multipole_frames.append(
                MultipoleFrame(
                    atoms[m.atom.idx], m.frame_pt_num, m.vectail, m.vechead, m.nvec)
            )
        for a in self.adjusts:
            c.adjusts.append(
                NonbondedException(
                    atoms[a.atom1.idx], atoms[a.atom2.idx], a.type and c.adjust_types[a.type.idx])
            )
        for a in self.acceptors:
            c.acceptors.append(
                AcceptorDonor(atoms[a.atom1.idx], atoms[a.atom2.idx])
            )
        for d in self.donors:
            c.donors.append(
                AcceptorDonor(atoms[d.atom1.idx], atoms[d.atom2.idx])
            )
        for g in self.groups:
            c.groups.append(Group(atoms[g.atom.idx], g.type, g.move))
        for l in self.links:
            c.links.append(
                Link(atoms[l.atom1.idx], atoms[l.atom2.idx],
                     l.length, l.symmetry_op1, l.symmetry_op2)
            )
        c._box = deepcopy(self._box)
        c._coordinates = deepcopy(self._coordinates)
        c.combining_rule = self.combining_rule
        # Transfer TER cards
        for r1, r2 in zip(c.residues, self.residues):
            r1.ter = r2.ter
        return c

    # ===================================================

    def to_dataframe(self):
        """ Generates a DataFrame from the current Structure's atomic properties

        Returns
        -------
        df : DataFrame
            DataFrame with all atomic properties

        See Also
        --------
        :func:`parmed.utils.pandautils.create_dataframe`
        """
        from .utils.pandautils import create_dataframe
        return create_dataframe(self)

    # ===================================================

    def load_dataframe(self, df):
        """ Loads atomic properties from an input DataFrame

        Parameters
        ----------
        df : pandas.DataFrame
            A pandas DataFrame with atomic properties that will be used to set
            the properties on the current list of atoms

        See Also
        --------
        :func:`parmed.utils.pandautils.load_dataframe`
        """
        from .utils.pandautils import load_dataframe
        return load_dataframe(self, df)

    # ===================================================

    def is_changed(self):
        """ Determines if any of the topology has changed for this structure """
        return (self.atoms.changed or self.residues.changed or
                self.bonds.changed or self.trigonal_angles.changed or
                self.dihedrals.changed or self.urey_bradleys.changed or
                self.impropers.changed or self.cmaps.changed or
                self.angles.changed or self.out_of_plane_bends.changed or
                self.pi_torsions.changed or self.stretch_bends.changed or
                self.torsion_torsions.changed or self.chiral_frames.changed or
                self.multipole_frames.changed or self.adjusts.changed or
                self.acceptors.changed or self.donors.changed or
                self.groups.changed or self.bond_types.changed or
                self.angle_types.changed or self.dihedral_types.changed or
                self.urey_bradley_types.changed or self.cmap_types.changed or
                self.improper_types.changed or self.adjust_types.changed or
                self.trigonal_angle_types.changed or
                self.out_of_plane_bends.changed or
                self.stretch_bend_types.changed or
                self.torsion_torsion_types.changed or
                self.pi_torsion_types.changed or self.rb_torsions.changed or
                self.rb_torsion_types.changed)

    # ===================================================

    def unchange(self):
        """ Toggles all lists so that they do not indicate any changes """
        self.atoms.changed = False
        self.residues.changed = False
        self.bonds.changed = False
        self.angles.changed = False
        self.dihedrals.changed = False
        self.urey_bradleys.changed = False
        self.impropers.changed = False
        self.cmaps.changed = False
        self.trigonal_angles.changed = False
        self.out_of_plane_bends.changed = False
        self.pi_torsions.changed = False
        self.stretch_bends.changed = False
        self.torsion_torsions.changed = False
        self.chiral_frames.changed = False
        self.multipole_frames.changed = False
        self.adjusts.changed = False
        self.acceptors.changed = False
        self.donors.changed = False
        self.groups.changed = False

        # Parameter type lists
        self.bond_types.changed = False
        self.angle_types.changed = False
        self.dihedral_types.changed = False
        self.urey_bradley_types.changed = False
        self.improper_types.changed = False
        self.cmap_types.changed = False
        self.trigonal_angle_types.changed = False
        self.out_of_plane_bend_types.changed = False
        self.pi_torsion_types.changed = False
        self.stretch_bend_types.changed = False
        self.torsion_torsion_types.changed = False
        self.adjust_types.changed = False

    # ===================================================

    def prune_empty_terms(self):
        """
        Looks through all of the topological lists and gets rid of terms
        in which at least one of the atoms is None or has an `idx` attribute set
        to -1 (indicating that it has been removed from the `atoms` atom list)
        """
        self._prune_empty_bonds()
        self._prune_empty_angles()
        self._prune_empty_dihedrals()
        self._prune_empty_rb_torsions()
        self._prune_empty_ureys()
        self._prune_empty_impropers()
        self._prune_empty_cmaps()
        self._prune_empty_trigonal_angles()
        self._prune_empty_out_of_plane_bends()
        self._prune_empty_pi_torsions()
        self._prune_empty_stretch_bends()
        self._prune_empty_torsion_torsions()
        self._prune_empty_chiral_frames()
        self._prune_empty_multipole_frames()
        self._prune_empty_adjusts()

    # ===================================================

    def update_dihedral_exclusions(self):
        """
        Nonbonded exclusions and exceptions have the following priority:

        bond -> angle -> dihedral

        Since bonds and angles are completely excluded, any ring systems in
        which two atoms are attached by a bond or angle as well as a dihedral
        should be completely excluded as the bond and angle exclusion rules take
        precedence.  If a Bond or Angle was _added_ to the structure between a
        pair of atoms previously connected only by a dihedral term, it's
        possible that those two atoms have both an exclusion *and* an exception
        defined. The result of this scenario is that sander and pmemd will
        happily compute an energy, _including_ the 1-4 nonbonded terms between
        atoms now connected by a bond or an Angle.  OpenMM, on the other hand,
        will complain about an exception specified multiple times. This method
        scans through all of the dihedrals in which `ignore_end` is `False` and
        turns it to `True` if the two end atoms are in the bond or angle
        partners arrays
        """
        set14 = set()
        deferred_dihedrals = []  # to work around pmemd tossing pn=0 dihedrals
        for dihedral in self.dihedrals:
            if dihedral.ignore_end:
                continue
            if (dihedral.atom1 in dihedral.atom4.bond_partners or
                    dihedral.atom1 in dihedral.atom4.angle_partners):
                dihedral.ignore_end = True
            elif (dihedral.atom1.idx, dihedral.atom4.idx) in set14:
                # Avoid double counting of 1-4 in a six-membered ring
                dihedral.ignore_end = True
            elif isinstance(dihedral.type, DihedralType) and dihedral.type.per == 0:
                # This needs to be done to work around a "feature" in pmemd where periodicity=0
                # torsions are thrown away. Since AmberParm never has any DihedralTypeList, we only
                # need to defer periodicity=0 torsions for DihedralType torsions.
                deferred_dihedrals.append(dihedral)
            else:
                set14.add((dihedral.atom1.idx, dihedral.atom4.idx))
                set14.add((dihedral.atom4.idx, dihedral.atom1.idx))
        # Only keep ignore_end = False if we *must*; i.e., if the exclusion it would have
        # added was not added by anybody else
        for dihedral in deferred_dihedrals:
            if (dihedral.atom1.idx, dihedral.atom4.idx) in set14:
                dihedral.ignore_end = True

    # ===================================================

    def strip(self, selection):
        """
        Deletes a subset of the atoms corresponding to an atom-based selection.

        Parameters
        ----------
        selection : :class:`AmberMask`, ``str``, or ``iterable``
            This is the selection of atoms that will be deleted from this
            structure. If it is a string, it will be interpreted as an
            AmberMask. If it is an AmberMask, it will be converted to a
            selection of atoms. If it is an iterable, it must be the same length
            as the `atoms` list.
        """
        from .amber import AmberMask
        if isinstance(selection, AmberMask):
            if selection.parm is not self:
                raise TypeError('passed mask does not belong to Structure')
            sel = selection.Selection()
        elif isinstance(selection, str):
            sel = AmberMask(self, selection).Selection()
        else:
            try:
                sel = list(selection)
            except TypeError:
                raise TypeError('Selection not a supported type [%s]' %
                                type(selection).__name__)
            if len(sel) != len(self.atoms):
                raise ValueError('Selection iterable wrong length')
        atomlist = sorted([i for i, s in enumerate(sel) if s])
        for i in reversed(atomlist):
            del self.atoms[i]
        self.prune_empty_terms()
        self.residues.prune()
        self.unchange()
        # Slice out coordinates if present
        if self._coordinates is not None:
            if PYPY:
                # numpypy does not currently support advanced indexing it seems
                self._coordinates = np.array(  # pragma: no cover
                    [[[x, y, z] for i, (x, y, z) in enumerate(crd)
                      if sel[i] == 0] for crd in self._coordinates]
                )
            else:
                # Pure numpy is faster in CPython, so do that when we can
                self._coordinates = self._coordinates[:, np.array(sel) == 0]

    # ===================================================

    def assign_bonds(self, *reslibs):
        """
        Assigns bonds to all atoms based on the provided residue template
        libraries. Atoms whose names are *not* in the templates, as well as
        those residues for whom no template is found, is assigned to bonds based
        on distances.

        Parameters
        ----------
        reslibs : dict{str: ResidueTemplate}
            Any number of residue template libraries. By default, assign_bonds
            knows about the standard amino acid, RNA, and DNA residues.
        """
        # Import here to avoid circular references
        from .modeller import StandardBiomolecularResidues
        # Build a composite dict of all residue templates
        all_residues = copy(StandardBiomolecularResidues)
        for lib in reslibs:
            all_residues.update(lib)
        # Walk through every residue and assign bonds from the templates
        unassigned_residues = set()
        unassigned_atoms = set()
        cysteine_sg = set()
        for res in self.residues:
            templ = _res_in_templlib(res, all_residues)
            if templ is None:
                # Don't have a template for this residue. Keep a note and go on
                unassigned_atoms.update(res.atoms)
                unassigned_residues.add(res)
                continue
            resatoms = {a.name: a for a in res.atoms}
            for a in res.atoms:
                if a.name not in templ.map:
                    unassigned_atoms.add(a)
                    continue
                # We matched our template atom. Transfer the element
                # information, as it is more reliable
                a.atomic_number = templ.map[a.name].atomic_number
                for bp in templ.map[a.name].bond_partners:
                    if (bp.name in resatoms and resatoms[bp.name] not in a.bond_partners):
                        if a not in resatoms[bp.name].bond_partners:
                            self.bonds.append(Bond(a, resatoms[bp.name]))
        # Now go through each residue and assign heads and tails. This walks
        # through the residues and bonds residue i's tail with residue j's head.
        # So there is nothing to do for the last residue. Omit it from the loop
        # so self.residues[i+1] never raises an IndexError
        for i, res in enumerate(self.residues[:-1]):
            templ = _res_in_templlib(res, all_residues)
            # TER cards and changing chains prevents bonding to the next residue
            if res.ter or (res.chain and (res.chain != self.residues[i+1].chain)):
                continue
            # We should also skip bonding to the next residue if there are detected gaps in the
            # residue sequence numbering
            if res.number != -1 and self.residues[i + 1].number != -1:
                try:
                    if self.residues[i + 1].number != res.number + 1:
                        continue
                except TypeError:
                    # The residue numbers are not numeric... weird, but move on as though no residue numbers were set
                    pass
            ntempl = _res_in_templlib(self.residues[i+1], all_residues)
            if templ is None and ntempl is None:
                # Any cross-link here picked up later
                continue
            if templ is None:
                if ntempl.head is None:
                    continue  # Next residue doesn't bond to the previous one
                # See if any atom in templ is close enough to bond to the head
                # atom of the next residue's template
                for head in self.residues[i+1].atoms:
                    if head.name == ntempl.head.name:
                        break
                else:
                    LOGGER.warning('Could not find the head atom of the next template! Bond '
                                   'pattern may be wrong, which could lead to extra TER cards in a '
                                   'PDB file')
                    continue  # head atom not found!
                for a in res.atoms:
                    maxdist = STANDARD_BOND_LENGTHS_SQUARED[(
                        a.atomic_number, head.atomic_number)]
                    if distance2(a, head) < maxdist:
                        if a not in head.bond_partners:
                            self.bonds.append(Bond(a, head))
                        break
                continue
            if templ.tail is None:
                continue  # This residue does not bond with the next one
            if ntempl is None:
                # See if any atom in the next residue is bonding distance away
                # from my tail
                for tail in res.atoms:
                    if tail.name == templ.tail.name:
                        break
                else:
                    continue  # tail not found
                for a in self.residues[i+1].atoms:
                    maxdist = STANDARD_BOND_LENGTHS_SQUARED[(
                        a.atomic_number, tail.atomic_number)]
                    if distance2(a, tail) < maxdist:
                        if a not in tail.bond_partners:
                            self.bonds.append(Bond(a, tail))
                        break
                continue
            if ntempl.head is None:
                continue  # Next residue does not bond with this one
            # We have templates for both atoms, and both have a head and a tail
            for tail in res.atoms:
                if tail.name == templ.tail.name:
                    break
            else:
                continue  # head could not be found
            for head in self.residues[i+1].atoms:
                if head.name == ntempl.head.name:
                    break
            else:
                continue  # tail could not be found
            if head not in tail.bond_partners:
                self.bonds.append(Bond(head, tail))
        # Now time to find bonds based on distance. First thing to do is to find
        # all CYS residues and pull out those that have an SG atom with only 1
        # bond partner, since those are *very* commonly bonded by
        # post-translational modifications.
        for res in self.residues:
            if len(res.name) != 3 or not residue.AminoAcidResidue.has(res.name):
                continue
            if residue.AminoAcidResidue.get(res.name).abbr != 'CYS':
                continue
            # Cysteine! Find the SG and add it to unassigned atoms if it only
            # has 1 atom bonded to it
            for a in res.atoms:
                if a.name == 'SG' and len(a.bond_partners) < 2:
                    unassigned_atoms.add(a)
                    cysteine_sg.add(a)
                    break
        # Can't do distance-based bond determination without coordinates
        if self.coordinates is None:
            return
        # OK, now go through all atoms that have not been assigned. If an atom
        # has not been assigned, but its residue HAS been assigned, then that
        # means the name is wrong. A likely culprit is bad nomenclature for
        # atoms. So in this case (and this case only), look for bond partners in
        # an 'assigned' residue (but only the same residue we are part of). One
        # thing this misses is when the head or tail atom is misnamed.
        mindist = math.sqrt(max(STANDARD_BOND_LENGTHS_SQUARED.values()))
        pairs = find_atom_pairs(self, mindist, unassigned_atoms)
        for atom in unassigned_atoms:
            for partner in pairs[atom.idx]:
                maxdist = STANDARD_BOND_LENGTHS_SQUARED[(
                    atom.atomic_number, partner.atomic_number)]
                if (distance2(atom, partner) < maxdist and atom not in partner.bond_partners):
                    self.bonds.append(Bond(atom, partner))
            # Now look through all atoms in this template if it's a template
            # that's already been assigned. If it's already in an unassigned
            # residue, then all that residue's atoms were in unassigned_atoms,
            # so there's no need to look through the residue again. Also, we
            # added cysteine SG atoms to the unassigned atoms list to see if it
            # is involved in post-translational bonds. But there *was* a
            # template match for cysteine and SG was paired correctly. So don't
            # try to assign it to bonds with other atoms in the cysteine here.
            if atom.residue in unassigned_residues or atom in cysteine_sg:
                continue
            for partner in atom.residue.atoms:
                if partner is atom:
                    continue
                maxdist = STANDARD_BOND_LENGTHS_SQUARED[(
                    atom.atomic_number, partner.atomic_number)]
                if (distance2(atom, partner) < maxdist and atom not in partner.bond_partners):
                    self.bonds.append(Bond(atom, partner))
        # All reasonable bonds have now been added

    # ===================================================

    def visualize(self, *args, **kwargs):
        """Use nglview for visualization. This only works with Jupyter notebook
        and require to install `nglview`

        Examples
        --------
        >>> import parmed as pmd
        >>> parm = pmd.download_PDB('1tsu')
        >>> parm.visualize()

        Parameters
        ----------
        args and kwargs : positional and keyword arguments given to nglview, optional
        """
        if self.coordinates is None:
            raise ValueError('coordinates must not be None')
        from nglview import show_parmed
        return show_parmed(self, *args, **kwargs)

    # ===================================================

    def __getitem__(self, selection):
        """
        Allows extracting a single atom from the structure or a slice of atoms
        as a new Structure instance. The following syntaxes are allowed:

            - struct[str] : str is interpreted as an Amber selection mask
            - struct[[sel,[sel,]],sel] : sel can be a list of indices (or
                            strings for chain selection) an integer, or a slice

        Parameters
        ----------
        selection : str, or slice|iter|str, slice|iter|int, slice|iter|int
            Which atoms to select for the new structure

        Returns
        -------
        struct : :class:`Structure`
            If more than one atom is selected, the resulting return value is a
            new structure containing all selected atoms and all valence terms
            (and types) in which all involved atoms are present

        or

        atom : :class:`Atom`
            If the selection is a single integer, a tuple of two integers, or a
            tuple of a string and two integers. This is equivalent to selecting
            a particular atom, a particular atom from a particular residue, or a
            particular atom from a particular residue from a particular chain,
            respectively. All other selections return a Structure instance, even
            if that selection happens to only select a single atom.

        Notes
        -----
        When selecting more than 1 atom, this is a costly operation. It is
        currently implemented by making a full copy of the object and then
        stripping the unused atoms.

        Raises
        ------
        ``TypeError`` : if selection (or any part of the selection) is not an
                        allowed type

        ``ValueError`` : if the selection is a boolean-like list and its length
                         is not the same as the number of atoms in the system
        """
        if isinstance(selection, int):
            return self.atoms[selection]

        selection = self._get_selection_array(selection)
        # _get_selection_array might have returned either an Atom or None if
        # there is no selection. Handle those and return, or continue if we got
        # a list of atoms
        if selection is None:
            return type(self)()
        elif isinstance(selection, Atom):
            return selection
        sumsel = sum(selection)
        if sumsel == 0:
            # No atoms selected. Return empty type
            return type(self)()
        # The cumulative sum of selection will give our index + 1 of each
        # selected atom into the new structure
        scan = [selection[0]]
        for i in range(1, len(selection)):
            scan.append(scan[i-1] + selection[i])
        # Zero-out the unselected atoms
        scan = [x * y for x, y in zip(scan, selection)]
        # Copy all parameters
        struct = type(self)()
        res_numbers_are_fake = {res.number for res in self.residues} == {0}
        res_number_attr = "idx" if res_numbers_are_fake else "number"
        for i, atom in enumerate(self.atoms):
            if not selection[i]:
                continue
            res = atom.residue
            num = getattr(res, res_number_attr)
            struct.add_atom(copy(atom), res.name, num, res.chain, res.insertion_code, res.segid)

        def copy_valence_terms(oval, otyp, sval, styp, atom_attrs, extra_attrs=None):
            """ Copies the valence terms from one list to another;
            oval=Other VALence; otyp=Other TYPe; sval=Self VALence;
            styp=Self TYPe; attrlist=ATTRibute LIST (atom1, atom2, ...)
            """
            otypcp = [copy(typ) for typ in styp]
            used_types = [False for typ in otypcp]
            for val in sval:
                ats = [getattr(val, attr) for attr in atom_attrs]
                # Make sure all of our atoms in this valence term is "selected"
                indices = [scan[at.idx] for at in ats if isinstance(at, Atom)]
                if not all(indices):
                    continue
                # Add the type if applicable
                kws = dict()
                if hasattr(val, 'type') and val.type is NoUreyBradley:
                    kws['type'] = NoUreyBradley  # special-case singleton
                elif otypcp and val.type is not None:
                    kws['type'] = otypcp[val.type.idx]
                    used_types[val.type.idx] = True
                for i, at in enumerate(ats):
                    if isinstance(at, Atom):
                        ats[i] = struct.atoms[scan[at.idx]-1]
                if extra_attrs:
                    kws.update({attr: getattr(val, attr)
                               for attr in extra_attrs})
                oval.append(type(val)(*ats, **kws))
                if hasattr(val, 'funct'):
                    oval[-1].funct = val.funct
            # Now tack on the "new" types copied from `other`
            for used, typ in zip(used_types, otypcp):
                if used:
                    otyp.append(typ)
            if hasattr(otyp, 'claim'):
                otyp.claim()
        copy_valence_terms(struct.bonds, struct.bond_types, self.bonds,
                           self.bond_types, ['atom1', 'atom2'], ["order", "qualitative_type"])
        copy_valence_terms(struct.angles, struct.angle_types, self.angles,
                           self.angle_types, ['atom1', 'atom2', 'atom3'])
        copy_valence_terms(struct.dihedrals, struct.dihedral_types,
                           self.dihedrals, self.dihedral_types,
                           ['atom1', 'atom2', 'atom3', 'atom4', 'improper', 'ignore_end'])
        copy_valence_terms(struct.rb_torsions, struct.rb_torsion_types,
                           self.rb_torsions, self.rb_torsion_types,
                           ['atom1', 'atom2', 'atom3', 'atom4', 'improper', 'ignore_end'])
        copy_valence_terms(struct.urey_bradleys, struct.urey_bradley_types, self.urey_bradleys,
                           self.urey_bradley_types, ['atom1', 'atom2'])
        copy_valence_terms(struct.impropers, struct.improper_types, self.impropers,
                           self.improper_types, ['atom1', 'atom2', 'atom3', 'atom4'])
        copy_valence_terms(struct.cmaps, struct.cmap_types, self.cmaps, self.cmap_types,
                           ['atom1', 'atom2', 'atom3', 'atom4', 'atom5'])
        copy_valence_terms(struct.trigonal_angles, struct.trigonal_angle_types,
                           self.trigonal_angles, self.trigonal_angle_types,
                           ['atom1', 'atom2', 'atom3', 'atom4'])
        copy_valence_terms(struct.out_of_plane_bends, struct.out_of_plane_bend_types,
                           self.out_of_plane_bends, self.out_of_plane_bend_types,
                           ['atom1', 'atom2', 'atom3', 'atom4'])
        copy_valence_terms(struct.pi_torsions, struct.pi_torsion_types,
                           self.pi_torsions, self.pi_torsion_types,
                           ['atom1', 'atom2', 'atom3', 'atom4', 'atom5', 'atom6'])
        copy_valence_terms(struct.stretch_bends, struct.stretch_bend_types, self.stretch_bends,
                           self.stretch_bend_types, ['atom1', 'atom2', 'atom3'])
        copy_valence_terms(struct.torsion_torsions, struct.torsion_torsion_types,
                           self.torsion_torsions, self.torsion_torsion_types,
                           ['atom1', 'atom2', 'atom3', 'atom4', 'atom5'])
        copy_valence_terms(struct.chiral_frames, [], self.chiral_frames, [],
                           ['atom1', 'atom2', 'chirality'])
        copy_valence_terms(struct.multipole_frames, [], self.multipole_frames, [],
                           ['atom', 'frame_pt_num', 'vectail', 'vechead', 'nvec'])
        copy_valence_terms(struct.adjusts, struct.adjust_types, self.adjusts,
                           self.adjust_types, ['atom1', 'atom2'])
        copy_valence_terms(struct.donors, [], self.donors,
                           [], ['atom1', 'atom2'])
        copy_valence_terms(struct.acceptors, [],
                           self.acceptors, [], ['atom1', 'atom2'])
        copy_valence_terms(struct.groups, [], self.groups,
                           [], ['atom', 'type', 'move'])
        struct._box = self._box
        struct.symmetry = self.symmetry
        struct.space_group = self.space_group
        return struct

    def _get_selection_array(self, selection):
        """
        Private method to convert a selection into an array -- common use for
        viewing and regular selection

        Parameters
        ----------
        selection : selector (slice, tuple, ... etc.)
            The selection given to the [] operator
        """
        from .amber import AmberMask
        # Now we select a subset of atoms. Convert "selection" into a natom list
        # with 0s and 1s, depending on what the input selection is
        if isinstance(selection, str):
            mask = AmberMask(self, selection)
            selection = mask.Selection()
        elif isinstance(selection, slice):
            sel = [0 for a in self.atoms]
            for idx in list(range(len(self.atoms)))[selection]:
                sel[idx] = 1
            selection = sel
        elif isinstance(selection, tuple) and len(selection) in (2, 3):
            # This is a selection of chains, and/or residues, and/or atoms
            if len(selection) == 2:
                # Select just residues and atoms. Special-case a single-atom
                # selection for speed -- orders of magnitude improvement in
                # efficiency
                ressel, atomsel = selection
                if isinstance(ressel, int) and isinstance(atomsel, int):
                    return self.residues[ressel][atomsel]
                has_chain = False
            elif len(selection) == 3:
                chainsel, ressel, atomsel = selection
                chainmap = defaultdict(TrackedList)
                for r in self.residues:
                    chainmap[r.chain].append(r)
                if (isinstance(chainsel, str) and isinstance(ressel, int) and
                        isinstance(atomsel, int)):
                    # Special-case single-atom selection for efficiency
                    chainmap = dict(chainmap)  # no longer defaultdict
                    try:
                        return chainmap[chainsel][ressel][atomsel]
                    except KeyError:
                        raise IndexError('No chain %s in Structure' % chainsel)
                if isinstance(chainsel, str):
                    if chainsel in chainmap:
                        chainset = set([chainsel])
                    else:
                        # The selected chain is not in the system. Cannot
                        # possibly select *any* atoms. Bail now for speed
                        return None
                elif isinstance(chainsel, slice):
                    # Build an ordered set of chains
                    chains = [self.residues[0].chain]
                    for res in self.residues:
                        if res.chain != chains[-1]:
                            chains.append(res.chain)
                    chainset = set(chains[chainsel])
                else:
                    chainset = set(chainsel)
                for chain in chainset:
                    if chain in chainmap:
                        break
                else:
                    # No requested chain is present. Bail now for speed
                    return None
                has_chain = True
            # Residue selection can either be by name or index
            if isinstance(ressel, slice):
                resset = set(list(range(len(self.residues)))[ressel])
            elif isinstance(ressel, str) or isinstance(ressel, int):
                resset = set([ressel])
            else:
                resset = set(ressel)
            if isinstance(atomsel, slice):
                atomset = set(list(range(len(self.atoms)))[atomsel])
            elif isinstance(atomsel, str) or isinstance(atomsel, int):
                atomset = set([atomsel])
            else:
                atomset = set(atomsel)
            if has_chain:
                try:
                    # To allow us to index the atoms residues from their list of
                    # chains, temporarily have the chainmap lists claim the
                    # residues. This must be reversed or the Structure will be
                    # broken
                    for chain_name, chain in chainmap.items():
                        chain.claim()
                    selection = [
                        (a.residue.chain in chainset) and
                        (a.residue.name in resset or a.residue.idx in resset) and
                        (a.name in atomset or a.idx -
                         a.residue[0].idx in atomset)

                        for a in self.atoms
                    ]
                finally:
                    # Make sure that the residues list *always* reclaims its
                    # contents
                    self.residues.claim()
            else:
                selection = [
                    (a.name in atomset or a.idx-a.residue[0].idx in atomset)
                    and (a.residue.name in resset or a.residue.idx in resset)
                    for a in self.atoms
                ]
        else:
            # Assume it is an iterable. If it is the same length as the atoms,
            # it is a boolean mask array. Otherwise, it is a list of atom
            # indices to select
            sel = [0 for atom in self.atoms]
            selection = list(selection)
            if len(selection) == len(self.atoms):
                for i, val in enumerate(selection):
                    if val:
                        sel[i] = 1
            elif len(selection) > len(self.atoms):
                raise ValueError('Selection iterable is too long')
            else:
                try:
                    for val in selection:
                        sel[val] = 1
                except IndexError:
                    raise ValueError('Selected atom out of range')
            selection = sel
        # Make sure we have an integer array of 0s and 1s
        return [int(bool(x)) for x in selection]

    # ===================================================

    @property
    def view(self):
        """
        Returns an indexable object that can be indexed like a standard
        Structure, but returns a *view* rather than a copy

        See Also
        --------
        Structure.__getitem__
        """
        return _StructureViewerCreator(self)

    # ===================================================

    def split(self):
        """
        Split the current Structure into separate Structure instances for each
        unique molecule. A molecule is defined as all atoms connected by a graph
        of covalent bonds.

        Returns
        -------
        [structs, counts] : list of (:class:`Structure`, list) tuples
            List of all molecules in the order that they appear first in the
            parent structure accompanied by the list of the molecule numbers
            in which that molecule appears in the Structure
        """
        tag_molecules(self)
        mollist = [atom.marked for atom in self.atoms]
        nmol = max(mollist)
        structs = []
        counts = []
        res_molecules = dict()
        # Divvy up the atoms into their separate molecules
        molatoms = [[] for i in range(nmol)]
        for atom in self.atoms:
            molatoms[atom.marked-1].append(atom)
        for i in range(nmol):
            sel = molatoms[i]
            involved_residues = set(atom.residue.idx for atom in sel)
            # Shortcut -- keep names of single-residue molecules and the number
            # of atoms present in those residues in a dict and use that to see
            # if two molecules are unique -- this gives drastic speedup for
            # systems with many duplicated single-residue molecules
            # (i.e., just about every solvent out there)
            if len(involved_residues) == 1:
                res = sel[0].residue
                names = tuple(a.name for a in res)
                charges = tuple('%.6f' % a.charge for a in res)
                rmins = tuple('%.6f' % a.rmin for a in res)
                epsilons = tuple('%.6f' % a.epsilon for a in res)
                oneres_key = (res.name, len(res), names,
                              charges, rmins, epsilons)
                if oneres_key in res_molecules:
                    counts[res_molecules[oneres_key]].add(i)
                    continue
            is_duplicate = False
            for j, struct in enumerate(structs):
                if len(struct.atoms) == len(sel):
                    for a1, a2 in zip(struct.atoms, sel):
                        assert None not in (
                            a1.residue, a2.residue), 'Residues must all be set'
                        if a1.residue.name != a2.residue.name:
                            break
                        if a1.name != a2.name:
                            break
                        if '%.6f' % a1.charge != '%.6f' % a2.charge:
                            break
                        if '%.6f' % a1.rmin != '%.6f' % a2.rmin:
                            break
                        if '%.6f' % a1.epsilon != '%.6f' % a2.epsilon:
                            break
                    else:
                        counts[j].add(i)
                        is_duplicate = True
                        break
            if not is_duplicate:
                if len(involved_residues) == 1:
                    # Cache single-residue molecules
                    res_molecules[oneres_key] = len(structs)
                mol = self[[atom.marked == i+1 for atom in self.atoms]]
                structs.append(mol)
                counts.append(set([i]))
        return list(zip(structs, counts))

    # ===================================================

    def save(self, fname, format=None, overwrite=False, **kwargs):
        """
        Saves the current Structure in the requested file format. Supported
        formats can be specified explicitly or determined by file-name
        extension. The following formats are supported, with the recognized
        suffix and ``format`` keyword shown in parentheses:

            - PDB (.pdb, pdb)
            - PDBx/mmCIF (.cif, cif)
            - PQR (.pqr, pqr)
            - Amber topology file (.prmtop/.parm7, amber)
            - CHARMM PSF file (.psf, psf)
            - CHARMM coordinate file (.crd, charmmcrd)
            - Gromacs topology file (.top, gromacs)
            - Gromacs GRO file (.gro, gro)
            - DLPOLY topology file (.field, field)
            - DLPOLY coordinate file (.config, config)
            - Mol2 file (.mol2, mol2)
            - Mol3 file (.mol3, mol3)
            - Amber ASCII restart (.rst7/.inpcrd/.restrt, rst7)
            - Amber NetCDF restart (.ncrst, ncrst)

        Parameters
        ----------
        fname : str or file-like object
            Name of the file or file-like object to save. If ``format`` is
            ``None`` (see below), the file type will be determined based on
            the filename extension. If ``fname`` is file-like object,  ``format``
            must be  provided. If the type cannot be determined, a ValueError is raised.
        format : str, optional
            The case-insensitive keyword specifying what type of file ``fname``
            should be saved as. If ``None`` (default), the file type will be
            determined from filename extension of ``fname``
        overwrite : bool, optional
            If True, allow the target file to be overwritten. Otherwise, an
            IOError is raised if the file exists. Default is False
        kwargs : keyword-arguments
            Remaining arguments are passed on to the file writing routines that
            are called by this function

        Raises
        ------
        ValueError if either filename extension or ``format`` are not recognized
        TypeError if the structure cannot be converted to the desired format for
        whatever reason
        IOError if the file cannot be written either because it exists and
        ``overwrite`` is ``False``, the filesystem is read-only, or write
        permissions are not granted for the user
        """
        from . import amber, charmm, formats, gromacs, dlpoly
        extmap = {
            '.pdb': 'PDB',
            '.pqr': 'PQR',
            '.cif': 'CIF',
            '.pdbx': 'CIF',
            '.parm7': 'AMBER',
            '.prmtop': 'AMBER',
            '.psf': 'PSF',
            '.top': 'GROMACS',
            '.gro': 'GRO',
            '.field': 'FIELD',
            '.config': 'CONFIG',
            '.mol2': 'MOL2',
            '.mol3': 'MOL3',
            '.crd': 'CHARMMCRD',
            '.rst7': 'RST7',
            '.inpcrd': 'RST7',
            '.restrt': 'RST7',
            '.ncrst': 'NCRST',
        }
        # Basically everybody uses atom type names instead of type indexes. So
        # convert to atom type names and switch back if need be
        if not hasattr(fname, 'write'):
            if os.path.exists(fname) and not overwrite:
                raise IOError('%s exists; not overwriting' % fname)
        else:
            if format is None:
                raise RuntimeError(
                    'Must provide supported format if using file-like object')
        all_ints = True
        for atom in self.atoms:
            if (isinstance(atom.type, int) and atom.atom_type is not UnassignedAtomType):
                atom.type = str(atom.atom_type)
            else:
                all_ints = False
        try:
            if format is not None:
                format = format.upper()
            else:
                base, ext = os.path.splitext(fname)
                if ext in ('.bz2', '.gz'):
                    ext = os.path.splitext(base)[1]
                try:
                    format = extmap[ext]
                except KeyError:
                    raise ValueError(
                        'Could not determine file type of %s' % fname)
            # Dispatch
            if format == 'PDB':
                self.write_pdb(fname, **kwargs)
            elif format == 'CIF':
                self.write_cif(fname, **kwargs)
            elif format == 'PQR':
                formats.PQRFile.write(self, fname, **kwargs)
            elif format == 'PSF':
                s = charmm.CharmmPsfFile.from_structure(self)
                s.write_psf(fname, **kwargs)
            elif format == 'GRO':
                gromacs.GromacsGroFile.write(self, fname, **kwargs)
            elif format == 'CONFIG':
                dlpoly.DlpolyConfigFile.write(self, fname, **kwargs)
            elif format == 'MOL2':
                formats.Mol2File.write(self, fname, **kwargs)
            elif format == 'MOL3':
                formats.Mol2File.write(self, fname, mol3=True, **kwargs)
            elif format == 'GROMACS':
                s = gromacs.GromacsTopologyFile.from_structure(self)
                s.write(fname, **kwargs)
            elif format == 'FIELD':
                s = dlpoly.DlpolyFieldFile.from_structure(self)
                s.write(fname, **kwargs)
            elif format == 'CHARMMCRD':
                charmm.CharmmCrdFile.write(self, fname, **kwargs)
            elif format == 'AMBER':
                if (self.trigonal_angles or self.out_of_plane_bends or self.torsion_torsions or
                        self.pi_torsions or self.stretch_bends or self.chiral_frames or
                        self.multipole_frames):
                    s = amber.AmoebaParm.from_structure(self)
                    s.write_parm(fname, **kwargs)
                elif self.urey_bradleys or self.impropers:
                    s = amber.ChamberParm.from_structure(self)
                    s.write_parm(fname, **kwargs)
                else:
                    try:
                        s = amber.AmberParm.from_structure(self)
                    except TypeError as e:
                        if 'Cannot translate exceptions' in str(e):
                            s = amber.ChamberParm.from_structure(self)
                        else:
                            raise  # pragma: no cover
                    s.write_parm(fname, **kwargs)
            elif format in ('RST7', 'NCRST'):
                rst7 = amber.Rst7(natom=len(self.atoms), **kwargs)
                rst7.coordinates = self.coordinates
                rst7.vels = self.velocities
                rst7.box = self.box
                rst7.write(fname, netcdf=(format == 'NCRST'))
            else:
                raise ValueError('No file type matching %s' % format)
        finally:
            if all_ints:
                for atom in self.atoms:
                    atom.type = int(atom.atom_type)

    # ===================================================

    def join_dihedrals(self):
        """
        Joins multi-term torsions into a single term and makes all of the
        parameters DihedralTypeList instances. If any dihedrals are *already*
        DihedralTypeList instances, or any are not parametrized, or there are no
        dihedral_types, this method returns without doing anything
        """
        if not self.dihedral_types:
            return  # nothing to do
        if any(isinstance(t, DihedralTypeList) for t in self.dihedral_types):
            return  # already done
        if any(d.type is None for d in self.dihedrals):
            return  # Not fully parametrized
        dihedrals_to_delete = list()
        dihedrals_processed = dict()
        new_dihedral_types = TrackedList()
        for i, d in enumerate(self.dihedrals):
            if d.atom1 < d.atom4:
                key = (d.atom1, d.atom2, d.atom3, d.atom4)
            else:
                key = (d.atom4, d.atom3, d.atom2, d.atom1)
            if key in dihedrals_processed:
                dihedrals_processed[key].append(d.type)
                dihedrals_to_delete.append(i)
            else:
                dihedrals_processed[key] = dtl = DihedralTypeList()
                dtl.append(d.type)
                new_dihedral_types.append(dtl)
                d.type = dtl
        # Now drop the new dihedral types into place
        self.dihedral_types = new_dihedral_types
        # Remove the "duplicate" dihedrals
        for i in reversed(dihedrals_to_delete):
            self.dihedrals[i].delete()
            del self.dihedrals[i]

    # ===================================================

    @property
    def combining_rule(self):
        return self._combining_rule

    @combining_rule.setter
    def combining_rule(self, thing):
        if thing not in ('lorentz', 'geometric'):
            raise ValueError("combining_rule must be 'lorentz' or 'geometric'")
        self._combining_rule = thing

    # ===================================================

    @property
    def rdkit_mol(self):
        if not self.is_changed:
            try:
                return self._rdkit_mol
            except AttributeError:
                pass
        else:
            self.prune_empty_terms()
            self.unchange()
        self._rdkit_mol = RDKit.to_mol(self)
        return self._rdkit_mol

    # ===================================================

    @property
    @needs_openmm
    def topology(self):
        """
        The OpenMM Topology object. Cached when possible, but any changes to the
        Structure instance results in the topology being deleted and rebuilt

        Notes
        -----
        This function calls ``prune_empty_terms`` if any topology lists have
        changed.
        """
        if not self.is_changed():
            try:
                return self._topology
            except AttributeError:
                pass
        else:
            self.prune_empty_terms()
            self.unchange()

        self._topology = top = app.Topology()
        try:
            last_chain = self.residues[0].chain
            last_residue = None
            last_omm_residue = None
        except IndexError:
            # Empty topology
            return self._topology
        chain = top.addChain(id=last_chain)
        # Add the atoms
        for i, atom in enumerate(self.atoms):
            # See if we need to add a new residue
            if atom.residue is not last_residue:
                # See if we need a new chain
                if last_chain != atom.residue.chain:
                    last_chain = atom.residue.chain
                    chain = top.addChain(id=last_chain)
                last_residue = atom.residue
                last_omm_residue = top.addResidue(
                    atom.residue.name,
                    chain,
                    id=str(atom.residue.number),
                    insertionCode=atom.residue.insertion_code,
                )
            try:
                elem = app.element.Element.getByAtomicNumber(
                    atom.atomic_number)
            except KeyError:
                elem = None
            top.addAtom(atom.name, elem, last_omm_residue)
        # Add the bonds
        bond_orders_by_openmm_type = {
            app.Single: 1.0,
            app.Double: 2.0,
            app.Triple: 3.0,
            app.Amide: 1.25,
            app.Aromatic: 1.5,
        }
        atoms = list(top.atoms())
        for bond in self.bonds:
            # Determine the bond type
            type = None
            for t, o in bond_orders_by_openmm_type.items():
                if math.isclose(bond.order, o):
                    type = t
            # Determine the bond order (an integer)
            order = None
            n = round(bond.order)
            if math.isclose(bond.order, n):
                order = n
            top.addBond(
                atoms[bond.atom1.idx], atoms[bond.atom2.idx], type=type, order=order
            )
        # Set the unit cell dimensions
        if self.box is not None:
            top.setPeriodicBoxVectors(
                reducePeriodicBoxVectors(
                    box_lengths_and_angles_to_vectors(*self.box))
            )
        return top

    # ===================================================

    @property
    def positions(self):
        """
        A list of 3-element Quantity tuples of dimension length representing the
        atomic positions for every atom in the system. If set with unitless
        numbers, those numbers are assumed to be in angstroms. If any atoms do
        not have coordinates, this is simply ``None``.
        """
        try:
            return [Vec3(a.xx, a.xy, a.xz) for a in self.atoms] * u.angstroms
        except AttributeError:
            return None

    @positions.setter
    def positions(self, value):
        """
        A list of 3-element Quantity tuples of dimension length representing the
        atomic positions for every atom in the system. If set with unitless
        numbers, those numbers are assumed to be in angstroms.

        Raises
        ------
        ValueError if the positions are not either a (natom, 3)-shape iterable
        of iterables or a (natom*3)-length iterable.
        """
        if u.is_quantity(value):
            value = value.value_in_unit(u.angstroms)
        # See if the array is flattened
        if len(value) == len(self.atoms):
            # It had better all be 3-length iterables
            for i, atom in enumerate(self.atoms):
                atom.xx, atom.xy, atom.xz = value[i]
        elif len(value) == 3 * len(self.atoms):
            for i, atom in enumerate(self.atoms):
                i3 = i * 3
                atom.xx, atom.xy, atom.xz = value[i3:i3+3]
        else:
            raise ValueError('Wrong shape for position array')

    @property
    def coordinates(self):
        try:
            coords = [[a.xx, a.xy, a.xz] for a in self.atoms]
        except AttributeError:
            return None
        return np.array(coords)

    @coordinates.setter
    def coordinates(self, value):
        """ Setting coordinates will also set xx, xy, and xz on the atoms """
        if value is None:
            # Wipe out coordinates
            self._coordinates = None
            for atom in self.atoms:
                try:
                    del atom.xx, atom.xy, atom.xz
                except AttributeError:
                    pass
        else:
            if u.is_quantity(value):
                value = value.value_in_unit(u.angstroms)
            value = list(value)
            coords = np.asanyarray(value, dtype=np.float64)
            coords = coords.reshape((-1, len(self.atoms), 3))
            if len(coords) > 0:
                for a, xyz in zip(self.atoms, coords[0]):
                    a.xx, a.xy, a.xz = xyz
                self._coordinates = coords
            else:
                # We set coordinates to an empty iterable. Wipe them out here
                for a in self.atoms:
                    try:
                        del a.xx, a.xy, a.xz
                    except AttributeError:
                        pass
                self._coordinates = None

    def get_coordinates(self, frame='all'):
        """
        In some cases, multiple conformations may be stored in the Structure.
        This function retrieves a particular frame's coordinates

        Parameters
        ----------
        frame : int or 'all', optional
            The frame number whose coordinates should be retrieved. Default is
            'all'

        Returns
        -------
        coords : np.ndarray, shape([#,] natom, 3) or None
            If frame is 'all', all coordinates are returned with shape
            (#, natom, 3). Otherwise the requested frame is returned with shape
            (natom, 3). If no coordinates exist and 'all' is requested, None is
            returned

        Raises
        ------
        IndexError if there are fewer than ``frame`` coordinates
        """
        try:
            coords = [[a.xx, a.xy, a.xz] for a in self.atoms]
        except AttributeError:
            coords = None
        else:
            coords = np.array(coords)
        # Wipe out our existing coordinates if we think anything might have changed
        if self._coordinates is not None:
            if coords is None:
                self._coordinates = None
            elif coords.shape != self._coordinates.shape[1:]:
                self._coordinates = None
            elif np.abs(coords - self._coordinates[0]).max() > SMALL:
                self._coordinates = None

        if frame == 'all':
            if self._coordinates is not None:
                return self._coordinates
            elif coords is not None:
                return coords.reshape((1, len(self.atoms), 3))
            return None
        elif self._coordinates is None:
            if frame == 0 and coords is not None:
                return coords
            # Requested NOT the first frame
            raise IndexError('No coordinate frames present')
        return self._coordinates[frame]

    @property
    def box(self):
        if self._box is None:
            return None
        return self._box[0]

    @box.setter
    def box(self, value):
        if value is None:
            self._box = None
        else:
            if isinstance(value, np.ndarray):
                box = value
            else:
                box = _strip_box_units(list(value))
            box = np.asanyarray(box, dtype=np.float64)
            if box.shape != (6,):
                if len(box.shape) != 2 or box.shape[-1] != 6:
                    raise ValueError('Box information must be 6 floats')
            self._box = box.reshape((-1, 6))

    def get_box(self, frame='all'):
        """
        In some cases, multiple conformations may be stored in the Structure.
        This function retrieves a particular frame's unit cell (box) dimensions

        Parameters
        ----------
        frame : int or 'all', optional
            The frame number whose unit cell should be retrieved. Default is
            'all'

        Returns
        -------
        box : np.ndarray, shape([#,] 6) or None
            If frame is 'all', all unit cells are returned with shape
            (#, 6). Otherwise the requested frame is returned with shape
            (6,). If no unit cell exist and 'all' is requested, None is
            returned

        Raises
        ------
        IndexError if there are fewer than ``frame`` unit cell dimensions
        """
        if self._box is None:
            if frame != 'all':
                raise IndexError('No unit cell frames present')
            return None
        elif frame == 'all':
            return self._box
        else:
            return self._box[frame]

    # ===================================================

    @property
    def velocities(self):
        """
        A (natom, 3)-shape numpy array with atomic velocities for every atom in
        the system (in units of angstrom/picosecond), or None if there are no
        velocities
        """
        try:
            return np.array([[a.vx, a.vy, a.vz] for a in self.atoms])
        except AttributeError:
            return None

    @velocities.setter
    def velocities(self, value):
        """
        A list of 3-element Quantity tuples of dimension length representing the
        atomic velocities for every atom in the system
        """
        if u.is_quantity(value):
            value = value.value_in_unit(u.angstroms/u.picoseconds)
        if value is None:
            for atom in self.atoms:
                try:
                    del atom.vx, atom.vy, atom.vz
                except AttributeError:
                    pass
        else:
            value = np.asarray(value).reshape(
                (-1, len(self.atoms), 3))
            for atom, xyz in zip(self.atoms, value[0]):
                atom.vx, atom.vy, atom.vz = xyz

    # ===================================================

    def has_NBFIX(self):
        """
        Returns whether or not any pairs of atom types have their LJ
        interactions modified by an NBFIX definition

        Returns
        -------
        has_nbfix : bool
            If True, at least two atom types have NBFIX mod definitions
        """
        typemap = dict()
        for a in self.atoms:
            if a.atom_type is UnassignedAtomType:
                continue
            typemap[str(a.atom_type)] = a.atom_type
        # Now we have a map of all atom types that we have defined in our
        # system. Look through all of the atom types and see if any of their
        # NBFIX definitions are also keys in typemap
        for key, type in typemap.items():
            for key in type.nbfix:
                if key in typemap:
                    return True
        return False

    # ===================================================

    @property
    def box_vectors(self):
        """
        3, 3-element tuple of unit cell vectors that are Quantity objects of
        dimension length
        """
        if self._box is None:
            return None
        return box_lengths_and_angles_to_vectors(*self.box)

    @box_vectors.setter
    def box_vectors(self, value):
        """
        3, 3-element tuple of unit cell vectors that are Quantity objects of
        dimension length
        """
        (a, b, c), (A, B, G) = box_vectors_to_lengths_and_angles(*value)
        a = a.value_in_unit(u.angstroms)
        b = b.value_in_unit(u.angstroms)
        c = c.value_in_unit(u.angstroms)
        A = A.value_in_unit(u.degrees)
        B = B.value_in_unit(u.degrees)
        G = G.value_in_unit(u.degrees)
        self._box = np.array([[a, b, c, A, B, G]], dtype=np.float64)

    # ===================================================

    @needs_openmm
    def createSystem(self, nonbondedMethod=None,
                     nonbondedCutoff=8.0*u.angstroms,
                     switchDistance=0.0*u.angstroms,
                     constraints=None,
                     rigidWater=True,
                     implicitSolvent=None,
                     implicitSolventKappa=None,
                     implicitSolventSaltConc=0.0*u.moles/u.liters,
                     temperature=298.15*u.kelvin,
                     soluteDielectric=1.0,
                     solventDielectric=78.5,
                     useSASA=False,
                     removeCMMotion=True,
                     hydrogenMass=None,
                     ewaldErrorTolerance=0.0005,
                     flexibleConstraints=True,
                     splitDihedrals=False,
                     drudeMass=0.4*u.dalton,
                     **kwargs):
        """
        Construct an OpenMM System representing the topology described by the
        prmtop file.

        Parameters
        ----------
        nonbondedMethod : cutoff method
            This is the cutoff method. It can be either the NoCutoff,
            CutoffNonPeriodic, CutoffPeriodic, PME, LJPME, or Ewald objects
            from the openmm.app namespace
        nonbondedCutoff : float or distance Quantity
            The nonbonded cutoff must be either a floating point number
            (interpreted as nanometers) or a Quantity with attached units. This
            is ignored if nonbondedMethod is NoCutoff.
        switchDistance : float or distance Quantity
            The distance at which the switching function is turned on for van
            der Waals interactions. This is ignored when no cutoff is used, and
            no switch is used if switchDistance is 0, negative, or greater than
            the cutoff
        constraints : None, app.HBonds, app.HAngles, or app.AllBonds
            Which type of constraints to add to the system (e.g., SHAKE). None
            means no bonds are constrained. HBonds means bonds with hydrogen are
            constrained
        rigidWater : bool=True
            If True, water is kept rigid regardless of the value of constraints.
            A value of False is ignored if constraints is not None.
        implicitSolvent : None, app.HCT, app.OBC1, app.OBC2, app.GBn, app.GBn2
            The Generalized Born implicit solvent model to use.
        implicitSolventKappa : float or 1/distance Quantity = None
            This is the Debye kappa property related to modeling saltwater
            conditions in GB. It should have units of 1/distance (1/nanometers
            is assumed if no units present). A value of None means that kappa
            will be calculated from implicitSolventSaltConc (below)
        implicitSolventSaltConc : float or amount/volume Quantity=0 moles/liter
            If implicitSolventKappa is None, the kappa will be computed from the
            salt concentration. It should have units compatible with mol/L
        temperature : float or temperature Quantity = 298.15 kelvin
            This is only used to compute kappa from implicitSolventSaltConc
        soluteDielectric : float=1.0
            The dielectric constant of the protein interior used in GB
        solventDielectric : float=78.5
            The dielectric constant of the water used in GB
        useSASA : bool=False
            If True, use the ACE non-polar solvation model. Otherwise, use no
            SASA-based nonpolar solvation model.
        removeCMMotion : bool=True
            If True, the center-of-mass motion will be removed periodically
            during the simulation. If False, it will not.
        hydrogenMass : float or mass quantity = None
            If not None, hydrogen masses will be changed to this mass and the
            difference subtracted from the attached heavy atom (hydrogen mass
            repartitioning)
        ewaldErrorTolerance : float=0.0005
            When using PME or Ewald, the Ewald parameters will be calculated
            from this value
        flexibleConstraints : bool=True
            If False, the energies and forces from the constrained degrees of
            freedom will NOT be computed. If True, they will (but those degrees
            of freedom will *still* be constrained).
        splitDihedrals : bool=False
            If True, the dihedrals will be split into two forces -- proper and
            impropers. This is primarily useful for debugging torsion parameter
            assignments.
        drudeMass : float or mass quantity = 0.4 dalton
            The mass to assign to the Drude particles. Default is 0.4 daltons

        Notes
        -----
        This function calls prune_empty_terms if any Topology lists have changed
        """
        if self.unknown_functional:
            raise ParameterError('Cannot createSystem: unknown functional')
        # Establish defaults
        if nonbondedMethod is None:
            nonbondedMethod = app.NoCutoff
        system = mm.System()
        # Make sure periodic simulations have a box
        if nonbondedMethod in (app.CutoffPeriodic, app.PME, app.Ewald, app.LJPME):
            if self.box is None:
                raise ValueError('No periodic boundary conditions detected')
        # Do hydrogen mass repartitioning if necessary
        masses = [atom.mass for atom in self.atoms]
        if hydrogenMass is not None:
            if u.is_quantity(hydrogenMass):
                hydrogenMass = hydrogenMass.value_in_unit(u.dalton)
            if hydrogenMass <= 0:
                raise ValueError('Hydrogen mass must be positive')
            for atom in self.atoms:
                if atom.element == 1:
                    heavy_atom = None
                    for a2 in atom.bond_partners:
                        if a2.element != 1:
                            heavy_atom = a2
                            break
                    if heavy_atom is not None:
                        masses[atom.idx] = hydrogenMass
                        masses[heavy_atom.idx] -= hydrogenMass - atom.mass
        for mass in masses:
            system.addParticle(mass)
        self.omm_add_constraints(system, constraints, rigidWater)
        # Prune empty terms if we have changed
        if self.is_changed():
            self.prune_empty_terms()
            self.unchange()
        # Add the various types of forces
        LOGGER.info('Adding bonds...')
        self._add_force_to_system(system, self.omm_bond_force(
            constraints, rigidWater, flexibleConstraints))
        LOGGER.info('Adding angles...')
        self._add_force_to_system(system, self.omm_angle_force(
            constraints, flexibleConstraints))
        LOGGER.info('Adding dihedrals...')
        self._add_force_to_system(
            system, self.omm_dihedral_force(splitDihedrals))
        LOGGER.info('Adding Ryckaert-Bellemans torsions...')
        self._add_force_to_system(system, self.omm_rb_torsion_force())
        LOGGER.info('Adding Urey-Bradleys...')
        self._add_force_to_system(system, self.omm_urey_bradley_force())
        LOGGER.info('Adding improper torsions...')
        self._add_force_to_system(system, self.omm_improper_force())
        LOGGER.info('Adding CMAP torsions...')
        self._add_force_to_system(system, self.omm_cmap_force())
        LOGGER.info('Adding trigonal angle terms...')
        self._add_force_to_system(system, self.omm_trigonal_angle_force())
        LOGGER.info('Adding out-of-plane bends...')
        self._add_force_to_system(system, self.omm_out_of_plane_bend_force())
        LOGGER.info('Adding pi-torsions...')
        self._add_force_to_system(system, self.omm_pi_torsion_force())
        LOGGER.info('Adding stretch-bends...')
        self._add_force_to_system(system, self.omm_stretch_bend_force())
        LOGGER.info('Adding torsion-torsions...')
        self._add_force_to_system(system, self.omm_torsion_torsion_force())
        LOGGER.info('Adding Nonbonded force...')
        if implicitSolvent is not None:
            rf_dielc = 1.0
        else:
            rf_dielc = 78.5
        self._add_force_to_system(system,
                                  self.omm_nonbonded_force(
                                      nonbondedMethod, nonbondedCutoff, switchDistance, ewaldErrorTolerance, rf_dielc)
                                  )
        if implicitSolvent is not None:
            LOGGER.info('Adding GB force...')
            self._add_force_to_system(system,
                                      self.omm_gbsa_force(implicitSolvent, nonbondedMethod, nonbondedCutoff,
                                                          soluteDielectric, solventDielectric, implicitSolventKappa,
                                                          implicitSolventSaltConc, temperature, useSASA)
                                      )
        if removeCMMotion:
            system.addForce(mm.CMMotionRemover())
        if self.box is not None:
            system.setDefaultPeriodicBoxVectors(
                *reducePeriodicBoxVectors(self.box_vectors))
        self.omm_set_virtual_sites(system)
        if any(isinstance(atom, DrudeAtom) for atom in self.atoms):
            self._add_force_to_system(system, self.omm_drude_force(drudeMass))
        return system

    # ===================================================

    @needs_openmm
    def omm_add_constraints(self, system, constraints, rigidWater):
        """ Adds constraints to a given system

        Parameters
        ----------
        system : mm.System
            The OpenMM system for which constraints should be added
        constraints : None, app.HBonds, app.AllBonds, or app.HAngles
            Which kind of constraints should be used
        rigidWater : bool
            If True, water bonds are constrained regardless of whether
            constrains is None
        """

        if constraints is None and not rigidWater:
            return
        if constraints not in (None, app.HBonds, app.AllBonds, app.HAngles):
            raise ValueError(
                f"Unrecognized constraints option ({constraints})")

        length_conv = u.angstrom.conversion_factor_to(u.nanometer)

        constraint_bond_set = set()
        constraint_angle_set = set()
        is_water = _settler(self)

        if constraints is app.AllBonds or constraints is app.HAngles:
            for bond in self.bonds:
                # Skip all extra points... don't constrain those
                if isinstance(bond.atom1, ExtraPoint):
                    continue
                if isinstance(bond.atom2, ExtraPoint):
                    continue
                atom1, atom2 = min(bond.atom1, bond.atom2), max(
                    bond.atom1, bond.atom2)
                constraint_bond_set.add(frozenset((atom1.idx, atom2.idx)))
        elif constraints is app.HBonds:
            for bond in self.bonds:
                if isinstance(bond.atom1, ExtraPoint):
                    continue
                if isinstance(bond.atom2, ExtraPoint):
                    continue
                if bond.atom1.element == 1 or bond.atom2.element == 1:
                    atom1, atom2 = min(bond.atom1, bond.atom2), max(
                        bond.atom1, bond.atom2)
                    constraint_bond_set.add(frozenset((atom1.idx, atom2.idx)))
        if rigidWater:
            for bond in self.bonds:
                if isinstance(bond.atom1, ExtraPoint):
                    continue
                if isinstance(bond.atom2, ExtraPoint):
                    continue
                if is_water[bond.atom1.residue.idx]:
                    atom1, atom2 = min(bond.atom1, bond.atom2), max(
                        bond.atom1, bond.atom2)
                    constraint_bond_set.add(frozenset((atom1.idx, atom2.idx)))

        # Add bond constraints
        for bond in self.bonds:
            if frozenset((bond.atom1.idx, bond.atom2.idx)) in constraint_bond_set:
                system.addConstraint(
                    bond.atom1.idx, bond.atom2.idx, bond.type.req * length_conv)

        if constraints is app.HAngles:
            for angle in self.angles:
                numH = 0
                if angle.atom1.element == 1:
                    numH += 1
                if angle.atom3.element == 1:
                    numH += 1
                if numH == 2 or (numH == 1 and angle.atom2.element == 8):
                    atom1, atom2, atom3 = min(angle.atom1, angle.atom3), angle.atom2, max(
                        angle.atom1, angle.atom3)
                    constraint_angle_set.add((atom1.idx, atom2.idx, atom3.idx))
        if rigidWater:
            for angle in self.angles:
                if is_water[angle.atom1.residue.idx]:
                    atom1, atom2, atom3 = min(angle.atom1, angle.atom3), angle.atom2, max(
                        angle.atom1, angle.atom3)
                    constraint_angle_set.add(
                        (angle.atom1.idx, angle.atom2.idx, angle.atom3.idx))
        # Add angle constraints
        for angle in self.angles:
            if (angle.atom1.idx, angle.atom2.idx, angle.atom3.idx) in constraint_angle_set:
                if frozenset((angle.atom1.idx, angle.atom3.idx)) in constraint_bond_set:
                    continue
                if frozenset((angle.atom3.idx, angle.atom1.idx)) in constraint_bond_set:
                    continue
                # Constrain this angle
                l1 = l2 = None
                for bond in angle.atom2.bonds:
                    if bond in angle and angle.atom1 in bond:
                        l1 = bond.type.req * length_conv
                    elif bond in angle and angle.atom3 in bond:
                        l2 = bond.type.req * length_conv
                # Law of cosines to find the constraint distance
                if l1 is None or l2 is None:
                    continue  # no bonds found...
                cost = math.cos(angle.type.theteq * DEG_TO_RAD)
                length = math.sqrt(l1 * l1 + l2 * l2 - 2 * l1 * l2 * cost)
                system.addConstraint(angle.atom1.idx, angle.atom3.idx, length)
    # ===================================================

    @needs_openmm
    def omm_set_virtual_sites(self, system):
        """
        Sets the virtual sites in a given OpenMM `System` object from the extra
        points defined in this system

        Parameters
        ----------
        system : mm.System
            The system for which the virtual sites will be set. All particles
            must have already been added to this System before calling this
            method
        """
        if system.getNumParticles() != len(self.atoms):
            raise ValueError('OpenMM System does not correspond to Structure')
        for atom in self.atoms:
            if not isinstance(atom, ExtraPoint):
                continue
            # This is a virtual site... get its frame type
            typ = atom.frame_type
            weights = typ.get_weights()
            refatoms = typ.get_atoms()
            if isinstance(typ, TwoParticleExtraPointFrame):
                a1, a2 = refatoms
                w1, w2 = weights
                system.setVirtualSite(
                    atom.idx, mm.TwoParticleAverageSite(a1.idx, a2.idx, w1, w2))
            elif isinstance(typ, ThreeParticleExtraPointFrame):
                a1, a2, a3 = refatoms
                w1, w2, w3 = weights
                system.setVirtualSite(
                    atom.idx, mm.ThreeParticleAverageSite(
                        a1.idx, a2.idx, a3.idx, w1, w2, w3)
                )
            elif isinstance(typ, OutOfPlaneExtraPointFrame):
                a1, a2, a3 = refatoms
                w1, w2, w3 = weights
                system.setVirtualSite(
                    atom.idx, mm.OutOfPlaneSite(
                        a1.idx, a2.idx, a3.idx, w1, w2, w3)
                )
            elif isinstance(typ, LocalCoordinatesFrame):
                a1, a2, a3 = refatoms
                origin_weights, x_weights, y_weights, local_position = weights
                # Convert positions to nanometers
                position_nm = [0.1 * position for position in local_position]
                system.setVirtualSite(
                    atom.idx, mm.LocalCoordinatesSite(
                        a1.idx, a2.idx, a3.idx, origin_weights, x_weights, y_weights, position_nm)
                )
            else:
                raise RuntimeError(
                    f"Did not recognize frame type {type(typ).__name__}")

    # ===================================================

    @needs_openmm
    def omm_bond_force(self, constraints=None, rigidWater=True, flexibleConstraints=True):
        """
        Creates an OpenMM Bond Force object (or AmoebaBondForce if the bonds are
        for an Amoeba-parametrized system)

        Parameters
        ----------
        constraints : None, app.HBonds, app.AllBonds, or app.HAngles
            The types of constraints that are on the system. If
            flexibleConstraints is False, then the constrained bonds will not be
            added to the resulting Force
        rigidWater : bool=True
            Should water-H bonds be constrained regardless of `constraints`?
        flexibleConstraints : bool=True
            If True, all bonds are added to the force regardless of
            `constraints`

        Returns
        -------
        force
            HarmonicBondForce (or AmoebaBondForce if this is an Amoeba system),
            or None if there are no bonds to add
        """
        if not flexibleConstraints and constraints in (app.HAngles,
                                                       app.AllBonds) or not self.bonds:
            return None  # No bonds to add
        length_conv = u.angstroms.conversion_factor_to(u.nanometers)
        _ambfrc = u.kilocalorie_per_mole/u.angstrom**2
        _ommfrc = u.kilojoule_per_mole/u.nanometer**2
        frc_conv = _ambfrc.conversion_factor_to(_ommfrc)
        # See if we need to add Amoeba bonds or regular bonds
        if hasattr(self.bond_types, 'degree') and hasattr(self.bond_types, 'coeffs'):
            energy = "k*(d^2 + %s*d^3 + %s*d^4); d=r-r0" % (
                self.bond_types.coeffs[3]/length_conv, self.bond_types.coeffs[4]/length_conv**2)
            force = mm.CustomBondForce(energy)
            force.addPerBondParameter('r0')
            force.addPerBondParameter('k')
            force.setName('AmoebaBond')
        else:
            force = mm.HarmonicBondForce()
        force.setForceGroup(self.BOND_FORCE_GROUP)
        # Add the bonds
        if rigidWater:
            is_water = _settler(self)
        else:
            is_water = [False for r in self.residues]
        for bond in self.bonds:
            if 1 in (bond.atom1.element, bond.atom2.element) and (
                    not flexibleConstraints and constraints is app.HBonds):
                continue
            if not flexibleConstraints and is_water[bond.atom1.residue.idx]:
                # is_water is False even for waters if rigidWater is False, so
                # we can rely on is_water to be correct for our needs
                # regardless
                continue
            if bond.type is None:
                raise ParameterError('Cannot find necessary parameters')
            if isinstance(force, mm.HarmonicBondForce):
                force.addBond(bond.atom1.idx, bond.atom2.idx, bond.type.req*length_conv,
                              2*bond.type.k*frc_conv)
            else:
                force.addBond(bond.atom1.idx, bond.atom2.idx, (bond.type.req*length_conv,
                              2*bond.type.k*frc_conv))
        # Done adding the force
        if force.getNumBonds() == 0:
            return None
        return force

    # ===================================================

    @needs_openmm
    def omm_angle_force(self, constraints=None, flexibleConstraints=True):
        """
        Creates an OpenMM HarmonicAngleForce object (or AmoebaAngleForce if the
        angles are for an Amoeba-parametrized system)

        Parameters
        ----------
        constraints : None, app.HBonds, app.AllBonds, or app.HAngles
            The types of constraints that are on the system. If
            flexibleConstraints is False, then the constrained bonds will not be
            added to the resulting Force
        flexibleConstraints : bool=True
            If True, all bonds are added to the force regardless of
            `constraints`

        Returns
        -------
        force
            HarmonicAngleForce (or AmoebaAngleForce if this is an Amoeba
            system), or None if there are no angles to add
        """
        if not self.angles:
            return None
        frc_conv = u.kilocalories.conversion_factor_to(u.kilojoules)
        if (hasattr(self.angle_types, 'degree') and
                hasattr(self.angle_types, 'coeffs')):
            c = self.angle_types.coeffs
            energy = "k*(d^2 + %s*d^3 + %s*d^4 + %s*d^5 + %s*d^6); d=%.15g*theta-theta0" % (
                c[3], c[4], c[5], c[6], 180/math.pi)
            force = mm.CustomAngleForce(energy)
            force.addPerAngleParameter('theta0')
            force.addPerAngleParameter('k')
            force.setName('AmoebaAngle')
        else:
            force = mm.HarmonicAngleForce()
        force.setForceGroup(self.ANGLE_FORCE_GROUP)
        for angle in self.angles:
            num_h = (angle.atom1.element == 1) + (angle.atom3.element == 1)
            if constraints is app.HAngles and (num_h == 2 or (num_h == 1 and
                                                              angle.atom2.element == 8) and not flexibleConstraints):
                continue  # pragma: no cover
            if angle.type is None:
                raise ParameterError('Cannot find angle parameters')
            if isinstance(force, mm.HarmonicAngleForce):
                force.addAngle(angle.atom1.idx, angle.atom2.idx, angle.atom3.idx,
                               angle.type.theteq*DEG_TO_RAD, 2*angle.type.k*frc_conv)
            else:
                force.addAngle(angle.atom1.idx, angle.atom2.idx, angle.atom3.idx,
                               (angle.type.theteq*DEG_TO_RAD, 2*angle.type.k*frc_conv))
        if force.getNumAngles() == 0:
            return None
        return force

    # ===================================================

    @needs_openmm
    def omm_dihedral_force(self, split=False):
        """ Creates the OpenMM PeriodicTorsionForce modeling dihedrals

        Parameters
        ----------
        split : bool, optional, default=False
            If True, separate PeriodicTorsionForce instances with the propers in
            the first and impropers in the second return item. If no impropers
            or propers are present, the instances with zero terms are not
            returned.

        Returns
        -------
        PeriodicTorsionForce[, PeriodicTorsionForce]
            Or None if no torsions are present in this system
        """
        if not self.dihedrals:
            return None
        frc_conv = u.kilocalories.conversion_factor_to(u.kilojoules)
        proper = mm.PeriodicTorsionForce()
        improper = mm.PeriodicTorsionForce()
        proper.setForceGroup(self.DIHEDRAL_FORCE_GROUP)
        improper.setForceGroup(self.IMPROPER_FORCE_GROUP)
        for tor in self.dihedrals:
            if tor.type is None:
                raise ParameterError('Cannot find torsion parameters')
            if split and tor.improper:
                force = improper
            else:
                force = proper
            if isinstance(tor.type, DihedralTypeList):
                for typ in tor.type:
                    if typ.per == 0:
                        warnings.warn(
                            "Torsion with 0-periodicity is being omitted. This will shift the energy by a constant, "
                            "but will not affect forces", ParmedWarning
                        )
                        continue
                    force.addTorsion(tor.atom1.idx, tor.atom2.idx, tor.atom3.idx, tor.atom4.idx,
                                     abs(int(typ.per)), typ.phase*DEG_TO_RAD, typ.phi_k*frc_conv)
            else:
                if tor.type.per == 0:
                    warnings.warn(
                        "Torsion with 0-periodicity is being omitted. This will shift the energy by a constant, "
                        "but will not affect forces", ParmedWarning
                    )
                    continue
                force.addTorsion(tor.atom1.idx, tor.atom2.idx, tor.atom3.idx, tor.atom4.idx,
                                 abs(int(tor.type.per)), tor.type.phase*DEG_TO_RAD,
                                 tor.type.phi_k*frc_conv)
        if proper.getNumTorsions() == 0:
            return improper
        if improper.getNumTorsions() == 0:
            return proper
        return proper, improper

    # ===================================================

    @needs_openmm
    def omm_rb_torsion_force(self):
        """ Creates the OpenMM RBTorsionForce for Ryckaert-Bellemans torsions

        Returns
        -------
        RBTorsionForce
            Or None if no torsions are present in this system
        """
        if not self.rb_torsions:
            return None
        conv = u.kilocalories.conversion_factor_to(u.kilojoules)
        force = mm.RBTorsionForce()
        force.setForceGroup(self.RB_TORSION_FORCE_GROUP)
        for tor in self.rb_torsions:
            if tor.type is None:
                raise ParameterError('Cannot find R-B torsion parameters')
            force.addTorsion(tor.atom1.idx, tor.atom2.idx, tor.atom3.idx, tor.atom4.idx,
                             tor.type.c0*conv, tor.type.c1*conv, tor.type.c2*conv,
                             tor.type.c3*conv, tor.type.c4*conv, tor.type.c5*conv)
        return force

    # ===================================================

    @needs_openmm
    def omm_urey_bradley_force(self):
        """ Creates the OpenMM Urey-Bradley force

        Returns
        -------
        HarmonicBondForce
            Or None, if no urey-bradleys are present
        """
        if not self.urey_bradleys:
            return None
        length_conv = u.angstroms.conversion_factor_to(u.nanometers)
        _ambfrc = u.kilocalorie_per_mole/u.angstrom**2
        _ommfrc = u.kilojoule_per_mole/u.nanometer**2
        frc_conv = _ambfrc.conversion_factor_to(_ommfrc)
        force = mm.HarmonicBondForce()
        force.setForceGroup(self.UREY_BRADLEY_FORCE_GROUP)
        for urey in self.urey_bradleys:
            if urey.type is None:
                raise ParameterError('Cannot find urey-bradley parameters')
            force.addBond(urey.atom1.idx, urey.atom2.idx, urey.type.req*length_conv,
                          2*urey.type.k*frc_conv)
        return force

    # ===================================================

    @needs_openmm
    def omm_improper_force(self):
        """ Creates the OpenMM improper torsion force (quadratic bias)

        Returns
        -------
        CustomTorsionForce
            With the formula k*(phi-phi0)^2, or None if there are no impropers
        """
        if not self.impropers:
            return None
        frc_conv = u.kilocalories.conversion_factor_to(u.kilojoules)
        energy_function = 'k*dtheta_torus^2;'
        energy_function += 'dtheta_torus = dtheta - floor(dtheta/(2*pi)+0.5)*(2*pi);'
        energy_function += 'dtheta = theta - theta0;'
        energy_function += 'pi = %f;' % math.pi
        force = mm.CustomTorsionForce(energy_function)
        force.addPerTorsionParameter('k')
        force.addPerTorsionParameter('theta0')
        force.setForceGroup(self.IMPROPER_FORCE_GROUP)
        for imp in self.impropers:
            if imp.type is None:
                raise ParameterError('Cannot find improper torsion parameters')
            force.addTorsion(imp.atom1.idx, imp.atom2.idx, imp.atom3.idx,
                             imp.atom4.idx, (imp.type.psi_k*frc_conv,
                                             imp.type.psi_eq*DEG_TO_RAD))
        return force

    # ===================================================

    @needs_openmm
    def omm_cmap_force(self):
        """ Creates the OpenMM CMAP torsion force

        Returns
        -------
        CMAPTorsionForce
            Or None, if no CMAP terms are present
        """
        if not self.cmaps:
            return None
        frc_conv = u.kilocalories.conversion_factor_to(u.kilojoules)
        force = mm.CMAPTorsionForce()
        force.setForceGroup(self.CMAP_FORCE_GROUP)
        # First get the list of cmap maps we're going to use. Just store the IDs
        # so we have simple integer comparisons to do later
        cmap_type_list = []
        cmap_map = dict()
        for cmap in self.cmaps:
            if cmap.type is None:
                raise ParameterError('Cannot find CMAP torsion parameters')
            if not id(cmap.type) in cmap_type_list:
                ct = cmap.type
                cmap_type_list.append(id(ct))
                # OpenMM takes the correction maps in the range 0 to 360
                # degrees, but we store it in -180 -- 180 (and in the transpose
                # of what OpenMM expects).
                grid = ct.grid.switch_range().T
                m = force.addMap(ct.resolution, [x*frc_conv for x in grid])
                cmap_map[id(ct)] = m
        # Now add all of the cmaps
        for cmap in self.cmaps:
            force.addTorsion(cmap_map[id(cmap.type)], cmap.atom1.idx, cmap.atom2.idx,
                             cmap.atom3.idx, cmap.atom4.idx, cmap.atom2.idx, cmap.atom3.idx,
                             cmap.atom4.idx, cmap.atom5.idx)
        return force

    # ===================================================

    @needs_openmm
    def omm_nonbonded_force(self, nonbondedMethod=None, nonbondedCutoff=8*u.angstroms,
                            switchDistance=0*u.angstroms, ewaldErrorTolerance=0.0005,
                            reactionFieldDielectric=78.5):
        """ Creates the OpenMM NonbondedForce instance

        Parameters
        ----------
        nonbondedMethod : cutoff method
            This is the cutoff method. It can be either the NoCutoff, CutoffNonPeriodic,
            CutoffPeriodic, PME, or Ewald objects from the openmm.app namespace
        nonbondedCutoff : float or distance Quantity
            The nonbonded cutoff must be either a floating point number
            (interpreted as nanometers) or a Quantity with attached units. This
            is ignored if nonbondedMethod is NoCutoff.
        switchDistance : float or distance Quantity
            The distance at which the switching function is turned on for van
            der Waals interactions. This is ignored when no cutoff is used, and
            no switch is used if switchDistance is 0, negative, or greater than the cutoff
        ewaldErrorTolerance : float=0.0005
            When using PME or Ewald, the Ewald parameters will be calculated from this value
        reactionFieldDielectric : float=78.5
            If the nonbondedMethod is CutoffPeriodic or CutoffNonPeriodic, the
            region beyond the cutoff is treated using a reaction field method
            with this dielectric constant. It should be set to 1 if another
            implicit solvent model is being used (e.g., GB)

        Returns
        -------
        NonbondedForce
            This just implements the very basic NonbondedForce with the typical
            charge-charge and 12-6 Lennard-Jones interactions with the
            Lorentz-Berthelot combining rules.

        Notes
        -----
        Subclasses of Structure for which this nonbonded treatment is inadequate
        should override this method to implement what is needed.

        If nrexcl is set to 3 and no exception parameters are stored in the
        adjusts list, the 1-4 interactions are determined from the list of
        dihedrals
        """
        if not self.atoms:
            return None
        length_conv = u.angstrom.conversion_factor_to(u.nanometer)
        ene_conv = u.kilocalories.conversion_factor_to(u.kilojoules)
        force = mm.NonbondedForce()
        force.setForceGroup(self.NONBONDED_FORCE_GROUP)
        if u.is_quantity(nonbondedCutoff):
            nonbondedCutoff = nonbondedCutoff.value_in_unit(u.nanometers)
        if nonbondedMethod is None or nonbondedMethod is app.NoCutoff:
            force.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
        elif nonbondedMethod is app.CutoffNonPeriodic:
            force.setNonbondedMethod(mm.NonbondedForce.CutoffNonPeriodic)
            force.setCutoffDistance(nonbondedCutoff)
        elif nonbondedMethod is app.CutoffPeriodic:
            force.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
            force.setCutoffDistance(nonbondedCutoff)
        elif nonbondedMethod is app.PME:
            force.setNonbondedMethod(mm.NonbondedForce.PME)
            force.setCutoffDistance(nonbondedCutoff)
            force.setEwaldErrorTolerance(ewaldErrorTolerance)
        elif nonbondedMethod is app.LJPME:
            force.setNonbondedMethod(mm.NonbondedForce.LJPME)
            force.setCutoffDistance(nonbondedCutoff)
            force.setEwaldErrorTolerance(ewaldErrorTolerance)
        elif nonbondedMethod is app.Ewald:
            force.setNonbondedMethod(mm.NonbondedForce.Ewald)
            force.setCutoffDistance(nonbondedCutoff)
            force.setEwaldErrorTolerance(ewaldErrorTolerance)
        else:
            raise ValueError(
                f'Unrecognized nonbondedMethod ({nonbondedMethod})')
        force.setReactionFieldDielectric(reactionFieldDielectric)
        # Now add the particles
        sigma_scale = length_conv * 2 * 2**(-1/6)
        for atom in self.atoms:
            force.addParticle(atom.charge, atom.sigma *
                              length_conv, abs(atom.epsilon*ene_conv))
        # Add exclusions from the bond graph out to nrexcl-1 bonds away (atoms
        # nrexcl bonds away will be exceptions defined later)

        def exclude_to(origin, atom, level, end):
            if level >= end:
                return
            for partner in atom.bond_partners:
                if partner is origin:
                    continue
                if atom is not origin:
                    force.addException(origin.idx, atom.idx,
                                       0.0, 0.5, 0.0, True)
                # Exclude EP children, too
                for child in origin.children:
                    force.addException(
                        partner.idx, child.idx, 0.0, 0.5, 0.0, True)
                for child in partner.children:
                    force.addException(child.idx, origin.idx,
                                       0.0, 0.5, 0.0, True)
                for child in origin.children:
                    for child2 in partner.children:
                        force.addException(
                            child.idx, child2.idx, 0.0, 0.5, 0.0, True)
                exclude_to(origin, partner, level+1, end)
        for atom in self.atoms:
            exclude_to(atom, atom, 0, self.nrexcl)
        # Add the exceptions from the dihedral list IFF no explicit exceptions
        # (or *adjusts*) have been specified. If dihedral.ignore_end is False, a
        # 1-4 with the appropriate scaling factor is used as the exception.
        sigma_scale = 2**(-1/6) * length_conv
        if self.combining_rule == 'lorentz':
            def comb_sig(sig1, sig2): return 0.5 * (sig1 + sig2)
        elif self.combining_rule == 'geometric':
            def comb_sig(sig1, sig2): return math.sqrt(sig1 * sig2)
        if not self.adjusts:
            for dih in self.dihedrals + self.rb_torsions:
                if dih.ignore_end:
                    continue
                if isinstance(dih.type, DihedralTypeList):
                    scee = scnb = 0
                    i = 0
                    while (scee == 0 or scnb == 0) and i < len(dih.type):
                        scee = dih.type[i].scee
                        scnb = dih.type[i].scnb
                        i += 1
                    if scee == 0 or scnb == 0:
                        raise ValueError('Detected scaling constants of 0 for '
                                         'dihedral containing 1-4 info!')
                else:
                    scee = dih.type.scee
                    scnb = dih.type.scnb
                try:
                    rij, wdij, rij14, wdij14 = dih.atom1.atom_type.nbfix[str(
                        dih.atom4.atom_type)]
                except (KeyError, AttributeError):
                    epsprod = abs(dih.atom1.epsilon_14 * dih.atom4.epsilon_14)
                    epsprod = math.sqrt(epsprod) * ene_conv / scnb
                    sigprod = comb_sig(dih.atom1.sigma_14, dih.atom4.sigma_14)
                    sigprod *= length_conv
                else:
                    epsprod = wdij14 * ene_conv / scnb
                    sigprod = rij14 * length_conv * sigma_scale
                chgprod = dih.atom1.charge * dih.atom4.charge / scee
                force.addException(dih.atom1.idx, dih.atom4.idx,
                                   chgprod, sigprod, epsprod, True)
                for child in dih.atom1.children:
                    epsprod = abs(child.epsilon_14 * dih.atom4.epsilon_14)
                    epsprod = math.sqrt(epsprod) * ene_conv / scnb
                    sigprod = comb_sig(child.sigma_14, dih.atom4.sigma_14)
                    sigprod *= length_conv
                    chgprod = (child.charge * dih.atom4.charge) / scee
                    force.addException(
                        child.idx, dih.atom4.idx, chgprod, sigprod, epsprod, True)
                for child in dih.atom4.children:
                    epsprod = abs(child.epsilon_14 * dih.atom1.epsilon_14)
                    epsprod = math.sqrt(epsprod) * ene_conv / scnb
                    sigprod = comb_sig(child.sigma_14, dih.atom1.sigma_14)
                    sigprod *= length_conv
                    chgprod = child.charge * dih.atom1.charge / scee
                    force.addException(
                        child.idx, dih.atom1.idx, chgprod, sigprod, epsprod, True)
                for c1 in dih.atom1.children:
                    for c2 in dih.atom2.children:
                        epsprod = abs(c1.epsilon_14 * c2.epsilon_14)
                        epsprod = math.sqrt(epsprod) * ene_conv / scnb
                        sigprod = comb_sig(c1.sigma_14, c2.sigma_14)
                        sigprod *= length_conv
                        chgprod = c1.charge * c2.charge / scee
                        force.addException(
                            c1.idx, c2.idx, chgprod, sigprod, epsprod, True)
        # Allow our specific exceptions (in adjusts) to override anything that
        # came before
        for pair in self.adjusts:
            chgprod = pair.atom1.charge * pair.atom2.charge * pair.type.chgscale
            force.addException(pair.atom1.idx, pair.atom2.idx, chgprod, pair.type.sigma*length_conv,
                               pair.type.epsilon*ene_conv, True)

        # Any exclusion partners we already added will zero-out existing
        # exclusions/exceptions
        for atom in self.atoms:
            for a2 in atom.exclusion_partners:
                force.addException(atom.idx, a2.idx, 0.0, 0.5, 0.0, True)

        if switchDistance and nonbondedMethod is not app.NoCutoff:
            if u.is_quantity(switchDistance):
                switchDistance = switchDistance.value_in_unit(u.nanometers)
            if 0 < switchDistance < nonbondedCutoff:
                force.setUseSwitchingFunction(True)
                force.setSwitchingDistance(switchDistance)

        # Now see if we have any NBFIXes that we need to implement via a
        # CustomNonbondedForce
        if self.combining_rule == 'geometric':
            if self.has_NBFIX():
                return (force, self._omm_nbfixed_force(force, nonbondedMethod),
                        self._omm_geometric_force(force, nonbondedMethod))
            else:
                return force, self._omm_geometric_force(force, nonbondedMethod)

        if self.has_NBFIX():
            return force, self._omm_nbfixed_force(force, nonbondedMethod)

        return force

    # ===================================================

    def omm_drude_force(self, drude_mass):
        force = mm.DrudeForce()
        force.setForceGroup(self.DRUDE_FORCE_GROUP)
        return force

    # ===================================================

    def _omm_nbfixed_force(self, nonbfrc, nonbondedMethod):
        """ Private method for creating a CustomNonbondedForce with a lookup
        table. This should not be called by users -- you have been warned.

        Parameters
        ----------
        nonbfrc : NonbondedForce
            NonbondedForce for the "standard" nonbonded interactions. This will
            be modified (specifically, L-J ixns will be zeroed)
        nonbondedMethod : Nonbonded Method (e.g., NoCutoff, PME, etc.)
            The nonbonded method to apply here. Ewald and PME will be
            interpreted as CutoffPeriodic for the CustomNonbondedForce.

        Returns
        -------
        force : CustomNonbondedForce
            The L-J force with NBFIX implemented as a lookup table
        """
        length_conv = u.angstroms.conversion_factor_to(u.nanometers)
        ene_conv = u.kilocalories.conversion_factor_to(u.kilojoules)
        # We need a CustomNonbondedForce to implement the NBFIX functionality.
        # First derive the type lookup tables
        lj_idx_list = [0 for atom in self.atoms]
        lj_radii, lj_depths = [], []
        num_lj_types = 0
        lj_type_list = []
        for i, atom in enumerate(self.atoms):
            atype = atom.atom_type
            if lj_idx_list[i]:
                continue  # already assigned
            num_lj_types += 1
            lj_idx_list[i] = num_lj_types
            ljtype = (atype.rmin, abs(atype.epsilon))
            lj_type_list.append(atype)
            lj_radii.append(atype.rmin)
            lj_depths.append(abs(atype.epsilon))
            for j in range(i+1, len(self.atoms)):
                atype2 = self.atoms[j].atom_type
                if lj_idx_list[j] > 0:
                    continue  # already assigned
                if atype2 is atype:
                    lj_idx_list[j] = num_lj_types
                elif not atype.nbfix:
                    # Only non-NBFIXed atom types can be compressed
                    ljtype2 = (atype2.rmin, abs(atype2.epsilon))
                    if ljtype == ljtype2:
                        lj_idx_list[j] = num_lj_types
        # Now everything is assigned. Create the A- and B-coefficient arrays
        acoef = [0 for i in range(num_lj_types*num_lj_types)]
        bcoef = acoef[:]
        for i in range(num_lj_types):
            for j in range(num_lj_types):
                namej = lj_type_list[j].name
                try:
                    rij, wdij, rij14, wdij14 = lj_type_list[i].nbfix[namej]
                except KeyError:
                    rij = (lj_radii[i] + lj_radii[j]) * length_conv
                    wdij = math.sqrt(lj_depths[i] * lj_depths[j]) * ene_conv
                else:
                    rij *= length_conv
                    wdij *= ene_conv
                rij6 = rij**6
                acoef[i+num_lj_types*j] = math.sqrt(wdij) * rij6
                bcoef[i+num_lj_types*j] = 2 * wdij * rij6
        force = mm.CustomNonbondedForce('(a/r6)^2-b/r6; r6=r2*r2*r2; r2=r^2; '
                                        'a=acoef(type1, type2); b=bcoef(type1, type2)')
        force.addTabulatedFunction('acoef', mm.Discrete2DFunction(
            num_lj_types, num_lj_types, acoef))
        force.addTabulatedFunction('bcoef', mm.Discrete2DFunction(
            num_lj_types, num_lj_types, bcoef))
        force.addPerParticleParameter('type')
        force.setForceGroup(self.NONBONDED_FORCE_GROUP)
        force.setUseLongRangeCorrection(True)
        if nonbondedMethod in (app.PME, app.Ewald, app.CutoffPeriodic, app.LJPME):
            force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        elif nonbondedMethod is app.NoCutoff:
            force.setNonbondedMethod(mm.CustomNonbondedForce.NoCutoff)
        elif nonbondedMethod is app.CutoffNonPeriodic:
            force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic)
        else:
            raise ValueError(
                f'Unrecognized nonbonded method [{nonbondedMethod}]')
        # Add the particles
        for i in lj_idx_list:
            force.addParticle((i-1,))
        # Now wipe out the L-J parameters in the nonbonded force
        for i in range(nonbfrc.getNumParticles()):
            chg, sig, eps = nonbfrc.getParticleParameters(i)
            nonbfrc.setParticleParameters(i, chg, 0.5, 0.0)
        # Now transfer the exclusions
        for ii in range(nonbfrc.getNumExceptions()):
            i, j, qq, ss, ee = nonbfrc.getExceptionParameters(ii)
            force.addExclusion(i, j)
        # Now transfer the other properties (cutoff, switching function, etc.)
        force.setCutoffDistance(nonbfrc.getCutoffDistance())
        if nonbfrc.getUseSwitchingFunction():
            force.setUseSwitchingFunction(True)
            force.setSwitchingDistance(nonbfrc.getSwitchingDistance())

        return force

    # ===================================================

    def _omm_geometric_force(self, nonbfrc, nonbondedMethod):
        """ Private method for creating a CustomNonbondedForce with a lookup
        table. This should not be called by users -- you have been warned.

        Parameters
        ----------
        nonbfrc : NonbondedForce
            NonbondedForce for the "standard" nonbonded interactions. This will
            be modified (specifically, L-J ixns will be zeroed)
        nonbondedMethod : Nonbonded Method (e.g., NoCutoff, PME, etc.)
            The nonbonded method to apply here. Ewald and PME will be
            interpreted as CutoffPeriodic for the CustomNonbondedForce.

        Returns
        -------
        force : CustomNonbondedForce
            The L-J force with NBFIX implemented as a lookup table
        """
        length_conv = u.angstroms.conversion_factor_to(u.nanometers)
        ene_conv = u.kilocalories.conversion_factor_to(u.kilojoules)
        # We need a CustomNonbondedForce to implement the geometric combining
        # rules
        force = mm.CustomNonbondedForce('epsilon1*epsilon2*(sigr6^2-sigr6); sigr6=sigr2*sigr2*sigr2; '
                                        'sigr2=(sigc/r)^2; sigc=sigma1*sigma2')
        force.addPerParticleParameter('epsilon')
        force.addPerParticleParameter('sigma')
        force.setForceGroup(self.NONBONDED_FORCE_GROUP)
        if (nonbondedMethod is app.PME or nonbondedMethod is app.Ewald or
                nonbondedMethod is app.CutoffPeriodic):
            force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        elif nonbondedMethod is app.NoCutoff:
            force.setNonbondedMethod(mm.CustomNonbondedForce.NoCutoff)
        elif nonbondedMethod is app.CutoffNonPeriodic:
            force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic)
        else:
            raise ValueError(
                f'Unrecognized nonbonded method [{nonbondedMethod}]')
        # Add the particles
        for atom in self.atoms:
            eps = math.sqrt(atom.epsilon*ene_conv) * 2
            sig = math.sqrt(atom.sigma*length_conv)
            force.addParticle((eps, sig))
        # Now wipe out the L-J parameters in the nonbonded force
        for i in range(nonbfrc.getNumParticles()):
            chg, sig, eps = nonbfrc.getParticleParameters(i)
            nonbfrc.setParticleParameters(i, chg, 0.5, 0.0)
        # Now transfer the exclusions
        for ii in range(nonbfrc.getNumExceptions()):
            i, j, qq, ss, ee = nonbfrc.getExceptionParameters(ii)
            force.addExclusion(i, j)
        # Now transfer the other properties (cutoff, switching function, etc.)
        force.setUseLongRangeCorrection(True)
        if nonbondedMethod is app.NoCutoff:
            force.setNonbondedMethod(mm.CustomNonbondedForce.NoCutoff)
        elif nonbondedMethod is app.CutoffNonPeriodic:
            force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic)
        elif nonbondedMethod in (app.PME, app.Ewald, app.CutoffPeriodic):
            force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        else:
            raise AssertionError(
                'Unsupported nonbonded method %s' % nonbondedMethod)
        force.setCutoffDistance(nonbfrc.getCutoffDistance())
        if nonbfrc.getUseSwitchingFunction():
            force.setUseSwitchingFunction(True)
            force.setSwitchingDistance(nonbfrc.getSwitchingDistance())

        return force

    # ===================================================

    @needs_openmm
    def omm_gbsa_force(self, implicitSolvent,
                       nonbondedMethod=None,
                       nonbondedCutoff=30.0*u.angstroms,
                       soluteDielectric=1.0,
                       solventDielectric=78.5,
                       implicitSolventKappa=None,
                       implicitSolventSaltConc=0.0*u.moles/u.liter,
                       temperature=298.15*u.kelvin,
                       useSASA=True):
        """
        Creates a Generalized Born force for running implicit solvent
        calculations

        Parameters
        ----------
        implicitSolvent : app.HCT, app.OBC1, app.OBC2, app.GBn, app.GBn2
            The Generalized Born implicit solvent model to use.
        nonbondedMethod : cutoff method
            This is the cutoff method. It can be either the NoCutoff,
            CutoffNonPeriodic, CutoffPeriodic, PME, or Ewald objects from the
            openmm.app namespace. Default is NoCutoff
        nonbondedCutoff : float or distance Quantity
            The nonbonded cutoff must be either a floating opint number
            (interpreted as nanometers) or a Quantity with attached units. This
            is ignored if nonbondedMethod is NoCutoff
        implicitSolventKappa : float or 1/distance Quantity = None
            This is the Debye kappa property related to modeling saltwater
            conditions in GB. It should have units of 1/distance (1/nanometers
            is assumed if no units present). A value of None means that kappa
            will be calculated from implicitSolventSaltConc (below)
        implicitSolventSaltConc : float or amount/volume Quantity=0 moles/liter
            If implicitSolventKappa is None, the kappa will be computed from the
            salt concentration. It should have units compatible with mol/L
        temperature : float or temperature Quantity = 298.15 kelvin
            This is only used to compute kappa from implicitSolventSaltConc
        soluteDielectric : float=1.0
            The dielectric constant of the protein interior used in GB
        solventDielectric : float=78.5
            The dielectric constant of the water used in GB
        """
        from openmm.app.internal.customgbforces import (
            GBSAHCTForce, GBSAOBC1Force, GBSAOBC2Force, GBSAGBnForce, GBSAGBn2Force
        )
        if implicitSolvent is None:
            return None
        if useSASA:
            sasa = 'ACE'
        else:
            sasa = None
        if nonbondedMethod is None:
            nonbondedMethod = app.NoCutoff
        if implicitSolvent not in (app.HCT, app.OBC1, app.OBC2, app.GBn, app.GBn2):
            raise ValueError('Unrecognized implicit solvent model')
        gb_parms = self._get_gb_parameters(implicitSolvent)
        if implicitSolventKappa is None:
            if u.is_quantity(implicitSolventSaltConc):
                sc = implicitSolventSaltConc.value_in_unit(u.moles/u.liter)
                implicitSolventSaltConc = sc
            if u.is_quantity(temperature):
                temperature = temperature.value_in_unit(u.kelvin)
            # The constant is 1 / sqrt(eps_0 * kB / (2*NA*q^2*1000)) where NA is
            # Avogadro's number, eps_0 is the permittivity of free space, q is
            # the charge (this # matches Amber's conversion factor)
            scale = math.sqrt(implicitSolventSaltConc /
                              solventDielectric / temperature)
            implicitSolventKappa = 50.33355 * scale
            # Multiply by 0.73 to account for ion exclusions, and multiply by 10
            # to convert to 1/nm from 1/angstroms
            implicitSolventKappa *= 7.3
        elif u.is_quantity(implicitSolventKappa):
            implicitSolventKappa = implicitSolventKappa.value_in_unit(
                u.nanometer**-1)

        if nonbondedMethod is app.NoCutoff:
            cutoff = None
        elif u.is_quantity(nonbondedCutoff):
            cutoff = nonbondedCutoff.value_in_unit(u.nanometers)
        else:
            cutoff = nonbondedCutoff
        if implicitSolvent is app.HCT:
            force = GBSAHCTForce(
                solventDielectric, soluteDielectric, sasa, cutoff, kappa=implicitSolventKappa)
        elif implicitSolvent is app.OBC1:
            force = GBSAOBC1Force(
                solventDielectric, soluteDielectric, sasa, cutoff, kappa=implicitSolventKappa)
        elif implicitSolvent is app.OBC2:
            force = GBSAOBC2Force(
                solventDielectric, soluteDielectric, sasa, cutoff, kappa=implicitSolventKappa)
        elif implicitSolvent is app.GBn:
            force = GBSAGBnForce(
                solventDielectric, soluteDielectric, sasa, cutoff, kappa=implicitSolventKappa)
        elif implicitSolvent is app.GBn2:
            force = GBSAGBn2Force(
                solventDielectric, soluteDielectric, sasa, cutoff, kappa=implicitSolventKappa)
        else:
            raise AssertionError(
                'Unexpected implicit solvent model... should not be here')
        for atom, parms in zip(self.atoms, gb_parms):
            force.addParticle([atom.charge] + list(parms))
        try:
            force.finalize()
        except AttributeError:
            # only new versions of omm require calling finalize
            pass
        # Set cutoff method
        if nonbondedMethod is app.NoCutoff:
            force.setNonbondedMethod(mm.CustomGBForce.NoCutoff)
        elif nonbondedMethod is app.CutoffNonPeriodic:
            force.setNonbondedMethod(mm.CustomGBForce.CutoffNonPeriodic)
            force.setCutoffDistance(cutoff)
        else:  # cutoff periodic (PME, CutoffPeriodic, Ewald)
            force.setNonbondedMethod(mm.CustomGBForce.CutoffPeriodic)
            force.setCutoffDistance(cutoff)
        force.setForceGroup(self.NONBONDED_FORCE_GROUP)

        return force

    # ===================================================

    # Amoeba-specific forces below

    @needs_openmm
    def omm_trigonal_angle_force(self):
        """ Creates the Amoeba trigonal-angle force

        Returns
        -------
        AmoebaInPlaneAngleForce
            The trigonal in-plane Angle force
        """
        if not self.trigonal_angles:
            return None
        frc_conv = u.kilocalories.conversion_factor_to(u.kilojoules)
        if (not hasattr(self.trigonal_angle_types, 'degree') or
                not hasattr(self.trigonal_angle_types, 'coeffs')):
            raise ParameterError(
                'Do not have the trigonal angle force table parameters')
        c = self.trigonal_angle_types.coeffs
        energy = """k*(d^2 + %s*d^3 + %s*d^4 + %s*d^5 + %s*d^6); d=theta-theta0;
                    theta = %.15g*pointangle(x1, y1, z1, projx, projy, projz, x3, y3, z3);
                    projx = x2-nx*dot; projy = y2-ny*dot; projz = z2-nz*dot;
                    dot = nx*(x2-x3) + ny*(y2-y3) + nz*(z2-z3);
                    nx = px/norm; ny = py/norm; nz = pz/norm;
                    norm = sqrt(px*px + py*py + pz*pz);
                    px = (d1y*d2z-d1z*d2y); py = (d1z*d2x-d1x*d2z); pz = (d1x*d2y-d1y*d2x);
                    d1x = x1-x4; d1y = y1-y4; d1z = z1-z4;
                    d2x = x3-x4; d2y = y3-y4; d2z = z3-z4""" % (c[3], c[4], c[5], c[6], 180/math.pi)
        force = mm.CustomCompoundBondForce(4, energy)
        force.addPerBondParameter("theta0")
        force.addPerBondParameter("k")
        force.setName('AmoebaInPlaneAngle')
        force.setForceGroup(self.TRIGONAL_ANGLE_FORCE_GROUP)
        for ang in self.trigonal_angles:
            if ang.type is None:
                raise ParameterError('Missing trigonal angle parameters')
            force.addBond((ang.atom1.idx, ang.atom2.idx, ang.atom3.idx, ang.atom4.idx),
                          (ang.type.theteq, ang.type.k*frc_conv))
        return force

    # ===================================================

    @needs_openmm
    def omm_out_of_plane_bend_force(self):
        """ Creates the Amoeba out-of-plane bend force

        Returns
        -------
        AmoebaOutOfPlaneBendForce
            The out-of-plane bend Angle force
        """
        if not self.out_of_plane_bends:
            return None
        frc_conv = u.kilocalories.conversion_factor_to(u.kilojoules)
        if (not hasattr(self.out_of_plane_bend_types, 'degree') or
                not hasattr(self.out_of_plane_bend_types, 'coeffs')):
            raise ParameterError(
                'Do not have the trigonal angle force table parameters')
        c = self.out_of_plane_bend_types.coeffs
        energy = """k*(theta^2 + %s*theta^3 + %s*theta^4 + %s*theta^5 + %s*theta^6);
                    theta = %.15g*pointangle(x2, y2, z2, x4, y4, z4, projx, projy, projz);
                    projx = x2-nx*dot; projy = y2-ny*dot; projz = z2-nz*dot;
                    dot = nx*(x2-x3) + ny*(y2-y3) + nz*(z2-z3);
                    nx = px/norm; ny = py/norm; nz = pz/norm;
                    norm = sqrt(px*px + py*py + pz*pz);
                    px = (d1y*d2z-d1z*d2y); py = (d1z*d2x-d1x*d2z); pz = (d1x*d2y-d1y*d2x);
                    d1x = x1-x4; d1y = y1-y4; d1z = z1-z4;
                    d2x = x3-x4; d2y = y3-y4; d2z = z3-z4""" % (c[3], c[4], c[5], c[6], 180/math.pi)
        force = mm.CustomCompoundBondForce(4, energy)
        force.addPerBondParameter("k")
        force.setName('AmoebaOutOfPlaneBend')
        force.setForceGroup(self.OUT_OF_PLANE_BEND_FORCE_GROUP)
        for ang in self.out_of_plane_bends:
            if ang.type is None:
                raise ParameterError('Missing out-of-plane bend parameters')
            force.addBond((ang.atom1.idx, ang.atom2.idx, ang.atom3.idx, ang.atom4.idx),
                          (2*ang.type.k*frc_conv,))
        return force

    # ===================================================

    @needs_openmm
    def omm_pi_torsion_force(self):
        """ Creates the Amoeba pi-torsion force

        Returns
        -------
        AmoebaPiTorsionForce
            The pi-torsion force
        """
        if not self.pi_torsions:
            return None
        frc_conv = u.kilocalories.conversion_factor_to(u.kilojoules)
        energy = """2*k*sin(phi)^2;
                    phi = pointdihedral(x3+c1x, y3+c1y, z3+c1z, x3, y3, z3, x4, y4, z4, x4+c2x, y4+c2y, z4+c2z);
                    c1x = (d14y*d24z-d14z*d24y); c1y = (d14z*d24x-d14x*d24z); c1z = (d14x*d24y-d14y*d24x);
                    c2x = (d53y*d63z-d53z*d63y); c2y = (d53z*d63x-d53x*d63z); c2z = (d53x*d63y-d53y*d63x);
                    d14x = x1-x4; d14y = y1-y4; d14z = z1-z4;
                    d24x = x2-x4; d24y = y2-y4; d24z = z2-z4;
                    d53x = x5-x3; d53y = y5-y3; d53z = z5-z3;
                    d63x = x6-x3; d63y = y6-y3; d63z = z6-z3"""
        force = mm.CustomCompoundBondForce(6, energy)
        force.addPerBondParameter('k')
        force.setName('AmoebaPiTorsion')
        force.setForceGroup(self.PI_TORSION_FORCE_GROUP)
        for ang in self.pi_torsions:
            if ang.type is None:
                raise ParameterError('Missing pi-torsion parameters')
            force.addBond((ang.atom1.idx, ang.atom2.idx, ang.atom3.idx, ang.atom4.idx,
                           ang.atom5.idx, ang.atom6.idx), (ang.type.phi_k*frc_conv,))
        return force

    # ===================================================

    @needs_openmm
    def omm_stretch_bend_force(self):
        """ Create the OpenMM Amoeba stretch-bend force for this system

        Returns
        -------
        AmoebaStretchBendForce
            The stretch-bend force containing all terms in this system
        """
        if not self.stretch_bends:
            return None
        # Conversion factor taken from pyopenmm/processTinkerForceField.py
        frc_conv = math.pi / 180 * 41.84  # 4.184 * 10
        length_conv = u.angstroms.conversion_factor_to(u.nanometers)
        energy = "(k1*(distance(p1,p2)-r12) + k2*(distance(p2,p3)-r23))*(%.15g*(angle(p1,p2,p3)-theta0))" % (180/math.pi)
        force = mm.CustomCompoundBondForce(3, energy)
        force.addPerBondParameter("r12")
        force.addPerBondParameter("r23")
        force.addPerBondParameter("theta0")
        force.addPerBondParameter("k1")
        force.addPerBondParameter("k2")
        force.setName('AmoebaStretchBend')
        force.setForceGroup(self.STRETCH_BEND_FORCE_GROUP)
        for strbnd in self.stretch_bends:
            if strbnd.type is None:
                raise ParameterError("Missing stretch-bend parameters")
            force.addBond((strbnd.atom1.idx, strbnd.atom2.idx, strbnd.atom3.idx),
                          (strbnd.type.req1*length_conv, strbnd.type.req2*length_conv,
                          strbnd.type.theteq*math.pi/180, strbnd.type.k1*frc_conv,
                          strbnd.type.k2*frc_conv))
        return force

    # ===================================================

    @needs_openmm
    def omm_torsion_torsion_force(self):
        """ Create the OpenMM Amoeba coupled-torsion (CMAP) force

        Returns
        -------
        AmoebaTorsionTorsionForce
            The torsion-torsion (CMAP) force with all coupled-torsion parameters
            for this system
        """
        if not self.torsion_torsions:
            return None
        # Not implemented yet...
        raise NotImplementedError("Torsion-torsions not yet implemented")

    # ===================================================

    def _prune_empty_bonds(self):
        """ Gets rid of any empty bonds """
        for i in reversed(range(len(self.bonds))):
            bond = self.bonds[i]
            if bond.atom1 is None and bond.atom2 is None:
                del self.bonds[i]
            elif bond.atom1.idx == -1 or bond.atom2.idx == -1:
                bond.delete()
                del self.bonds[i]

    # ===================================================

    def _prune_empty_angles(self):
        """ Gets rid of any empty angles """
        for i in reversed(range(len(self.angles))):
            angle = self.angles[i]
            if (angle.atom1 is None and angle.atom2 is None and
                    angle.atom3 is None):
                del self.angles[i]
            elif -1 in (angle.atom1.idx, angle.atom2.idx, angle.atom3.idx):
                angle.delete()
                del self.angles[i]

    # ===================================================

    def _prune_empty_dihedrals(self, dlist='dihedrals'):
        """ Gets rid of any empty dihedrals """
        for i in reversed(range(len(getattr(self, dlist)))):
            dihed = getattr(self, dlist)[i]
            if (dihed.atom1 is None and dihed.atom2 is None and
                    dihed.atom3 is None and dihed.atom4 is None):
                del getattr(self, dlist)[i]
            elif -1 in (dihed.atom1.idx, dihed.atom2.idx, dihed.atom3.idx, dihed.atom4.idx):
                dihed.delete()
                del getattr(self, dlist)[i]

    # ===================================================

    def _prune_empty_rb_torsions(self):
        """ Gets rid of any empty R-B torsions """
        self._prune_empty_dihedrals('rb_torsions')

    # ===================================================

    def _prune_empty_ureys(self):
        """ Gets rid of any empty Urey-Bradley terms """
        for i in reversed(range(len(self.urey_bradleys))):
            ub = self.urey_bradleys[i]
            if ub.atom1 is None and ub.atom2 is None:
                del self.urey_bradleys[i]
            elif ub.atom1.idx == -1 or ub.atom2.idx == -1:
                ub.delete()
                del self.urey_bradleys[i]

    # ===================================================

    def _prune_empty_impropers(self):
        """ Gets rid of any empty improper torsions """
        for i in reversed(range(len(self.impropers))):
            imp = self.impropers[i]
            if imp.atom1 is None and imp.atom2 is None and imp.atom3 is None and imp.atom4 is None:
                del self.impropers[i]
            elif -1 in (imp.atom1.idx, imp.atom2.idx, imp.atom3.idx, imp.atom4.idx):
                imp.delete()
                del self.impropers[i]

    # ===================================================

    def _prune_empty_cmaps(self):
        """ Gets rid of any empty CMAP terms """
        for i in reversed(range(len(self.cmaps))):
            cmap = self.cmaps[i]
            if (cmap.atom1 is None and cmap.atom2 is None and cmap.atom3 is None
                    and cmap.atom4 is None and cmap.atom5 is None):
                del self.cmaps[i]
            elif -1 in (cmap.atom1.idx, cmap.atom2.idx, cmap.atom3.idx,
                        cmap.atom4.idx, cmap.atom5.idx):
                cmap.delete()
                del self.cmaps[i]

    # ===================================================

    def _prune_empty_trigonal_angles(self):
        """ Gets rid of any empty trigonal angles """
        self._prune_empty_dihedrals('trigonal_angles')

    # ===================================================

    def _prune_empty_out_of_plane_bends(self):
        """ Gets rid of any empty out-of-plane bends """
        self._prune_empty_dihedrals('out_of_plane_bends')

    # ===================================================

    def _prune_empty_pi_torsions(self):
        """ Gets rid of any empty pi-torsions """
        for i in reversed(range(len(self.pi_torsions))):
            pit = self.pi_torsions[i]
            if (pit.atom1 is None and pit.atom2 is None and pit.atom3 is None and
                    pit.atom4 is None and pit.atom5 is None and pit.atom6 is None):
                del self.pi_torsions[i]
            elif -1 in (pit.atom1.idx, pit.atom2.idx, pit.atom3.idx, pit.atom4.idx,
                        pit.atom5.idx, pit.atom6.idx):
                # Not stored anywhere, no need to call delete()
                del self.pi_torsions[i]

    # ===================================================

    def _prune_empty_stretch_bends(self):
        """ Gets rid of any empty stretch-bend terms """
        for i in reversed(range(len(self.stretch_bends))):
            sb = self.stretch_bends[i]
            if sb.atom1 is None and sb.atom2 is None and sb.atom3 is None:
                del self.stretch_bends[i]
            elif sb.atom1.idx == -1 or sb.atom2.idx == -1 or sb.atom3.idx == -1:
                # Not stored anywhere, no need to call delete()
                del self.stretch_bends[i]

    # ===================================================

    def _prune_empty_torsion_torsions(self):
        """ Gets rid of any empty torsion-torsion terms """
        for i in reversed(range(len(self.torsion_torsions))):
            tt = self.torsion_torsions[i]
            if (tt.atom1 is None and tt.atom2 is None and tt.atom3 is None and
                    tt.atom4 is None and tt.atom5 is None):
                del self.torsion_torsions[i]
            elif -1 in (tt.atom1.idx, tt.atom2.idx, tt.atom3.idx, tt.atom4.idx, tt.atom5.idx):
                tt.delete()
                del self.torsion_torsions[i]

    # ===================================================

    def _prune_empty_chiral_frames(self):
        """ Gets rid of any empty chiral frame terms """
        for i in reversed(range(len(self.chiral_frames))):
            cf = self.chiral_frames[i]
            if cf.atom1 is None or cf.atom2 is None:
                del self.chiral_frames[i]
            elif cf.atom1.idx == -1 or cf.atom2.idx == -1:
                del self.chiral_frames[i]

    # ===================================================

    def _prune_empty_multipole_frames(self):
        """ Gets rid of any empty multipole frame terms """
        for i in reversed(range(len(self.multipole_frames))):
            mf = self.multipole_frames[i]
            if mf.atom is None or mf.atom.idx == -1:
                del self.multipole_frames[i]

    # ===================================================

    def _prune_empty_adjusts(self):
        """ Gets rid of any empty nonbonded exception adjustments """
        for i in reversed(range(len(self.adjusts))):
            adj = self.adjusts[i]
            if adj.atom1 is None or adj.atom2 is None:
                del self.adjusts[i]
            elif adj.atom1.idx == -1 or adj.atom2.idx == -1:
                del self.adjusts[i]

    # ===================================================

    @needs_openmm
    def _get_gb_parameters(self, implicitSolvent):
        """ Gets the GB parameters for the requested GB model used by OpenMM

        Parameters
        ----------
        implicitSolvent : app.HCT, app.OBC1, app.OBC2, app.GBn, or app.GBn2
            The object specifying a particular GB model in OpenMM

        Returns
        -------
        parameters : list of float
            List of parameters for the requested GB model
        """
        if implicitSolvent is app.GBn:
            screen = [0.5 for atom in self.atoms]
            radii = [atom.solvent_radius for atom in self.atoms]
            for i, atom in enumerate(self.atoms):
                if atom.element == 6:
                    screen[i] = 0.48435382330
                elif atom.element == 1:
                    screen[i] = 1.09085413633
                elif atom.element == 7:
                    screen[i] = 0.700147318409
                elif atom.element == 8:
                    screen[i] = 1.06557401132
                elif atom.element == 16:
                    screen[i] = 0.602256336067
                if radii[i] == 0:
                    radii[i] = _bondi(atom)
        elif implicitSolvent is app.GBn2:
            # Add non-optimized values as defaults
            alpha = [1.0 for i in self.atoms]
            beta = [0.8 for i in self.atoms]
            gamma = [4.85 for i in self.atoms]
            screen = [0.5 for i in self.atoms]
            radii = [atom.solvent_radius for atom in self.atoms]
            for i, atom in enumerate(self.atoms):
                if atom.element == 6:
                    screen[i] = 1.058554
                    alpha[i] = 0.733756
                    beta[i] = 0.506378
                    gamma[i] = 0.205844
                elif atom.element == 1:
                    screen[i] = 1.425952
                    alpha[i] = 0.788440
                    beta[i] = 0.798699
                    gamma[i] = 0.437334
                elif atom.element == 7:
                    screen[i] = 0.733599
                    alpha[i] = 0.503364
                    beta[i] = 0.316828
                    gamma[i] = 0.192915
                elif atom.element == 8:
                    screen[i] = 1.061039
                    alpha[i] = 0.867814
                    beta[i] = 0.876635
                    gamma[i] = 0.387882
                elif atom.element == 16:
                    screen[i] = -0.703469
                    alpha[i] = 0.867814
                    beta[i] = 0.876635
                    gamma[i] = 0.387882
                if not radii[i]:
                    radii[i] = _mbondi3(atom)
        else:
            radii = [atom.solvent_radius for atom in self.atoms]
            screen = [atom.screen for atom in self.atoms]
            for i, atom in enumerate(self.atoms):
                if not radii[i] or not screen[i]:
                    # Replace with defaults
                    radii[i], screen[i] = _gb_rad_screen(atom, implicitSolvent)

        length_conv = u.angstrom.conversion_factor_to(u.nanometer)
        radii = [x * length_conv for x in radii]

        if implicitSolvent is app.GBn2:
            return list(zip(radii, screen, alpha, beta, gamma))
        return list(zip(radii, screen))

    # ===================================================

    def __str__(self):
        if hasattr(self, 'name') and self.name:
            return self.name
        return repr(self)

    # ===================================================

    # Support concatenating and replicating the Structure

    def __add__(self, other):
        cp = copy(self)
        cp += other
        return cp

    def __iadd__(self, other):
        if not isinstance(other, Structure):
            return NotImplemented
        # Cache my coordinates, since changing the structure will destroy the
        # coordinate list
        mycrd = self.get_coordinates('all')
        # Make a copy of the "other" structure so we can change the atom types if
        # there are any clashes with those in `self`
        other = copy(other)
        other._deduplicate_atom_type_names(self)
        # The basic approach taken here is to extend the atom list, then scan
        # through all of the valence terms in `other`, adding them to the
        # corresponding arrays of `self`, using an offset to look into the atom
        # array. The types are done a little differently. The types are all
        # copied into a "temporary" array so they have the same index as the
        # original Structure `other`, and as the valence terms are created the
        # reference to that type is added. Afterwards, those types are tacked on
        # to the end of the types list on `self`
        aoffset = len(self.atoms)
        try:
            roffset = self.residues[-1].number + 1
        except IndexError:
            # No residues... must be an empty structure
            roffset = 0
        for atom in other.atoms:
            res = atom.residue
            self.add_atom(copy(atom), res.name, res.idx+roffset,
                          res.chain, res.insertion_code, res.segid)

        def copy_valence_terms(oval, otyp, sval, styp, attrlist):
            """ Copies the valence terms from one list to another;
            oval=Other VALence; otyp=Other TYPe; sval=Self VALence;
            styp=Self TYPe; attrlist=ATTRibute LIST (atom1, atom2, ...)
            """
            otypcp = [copy(typ) for typ in otyp]
            for val in oval:
                ats = [getattr(val, attr) for attr in attrlist]
                # Replace any atoms with the ones from self.atoms
                for i, at in enumerate(ats):
                    if isinstance(at, Atom):
                        ats[i] = self.atoms[at.idx+aoffset]
                kws = dict()
                if hasattr(val, 'type') and val.type is NoUreyBradley:
                    kws['type'] = NoUreyBradley  # special-case singleton
                elif otypcp and val.type is not None:
                    kws['type'] = otypcp[val.type.idx]
                sval.append(type(val)(*ats, **kws))
                if hasattr(val, 'funct'):
                    sval[-1].funct = val.funct
            # Now tack on the "new" types copied from `other`
            styp.extend(otypcp)
            if hasattr(styp, 'claim'):
                styp.claim()
        copy_valence_terms(other.bonds, other.bond_types, self.bonds,
                           self.bond_types, ['atom1', 'atom2'])
        copy_valence_terms(other.angles, other.angle_types, self.angles,
                           self.angle_types, ['atom1', 'atom2', 'atom3'])
        copy_valence_terms(other.dihedrals, other.dihedral_types,
                           self.dihedrals, self.dihedral_types,
                           ['atom1', 'atom2', 'atom3', 'atom4', 'improper', 'ignore_end'])
        copy_valence_terms(other.rb_torsions, other.rb_torsion_types,
                           self.rb_torsions, self.rb_torsion_types,
                           ['atom1', 'atom2', 'atom3', 'atom4', 'improper', 'ignore_end'])
        copy_valence_terms(other.urey_bradleys, other.urey_bradley_types, self.urey_bradleys,
                           self.urey_bradley_types, ['atom1', 'atom2'])
        copy_valence_terms(other.impropers, other.improper_types, self.impropers,
                           self.improper_types, ['atom1', 'atom2', 'atom3', 'atom4'])
        copy_valence_terms(other.cmaps, other.cmap_types, self.cmaps, self.cmap_types,
                           ['atom1', 'atom2', 'atom3', 'atom4', 'atom5'])
        copy_valence_terms(other.trigonal_angles, other.trigonal_angle_types,
                           self.trigonal_angles, self.trigonal_angle_types,
                           ['atom1', 'atom2', 'atom3', 'atom4'])
        copy_valence_terms(other.out_of_plane_bends, other.out_of_plane_bend_types,
                           self.out_of_plane_bends, self.out_of_plane_bend_types,
                           ['atom1', 'atom2', 'atom3', 'atom4'])
        copy_valence_terms(other.pi_torsions, other.pi_torsion_types,
                           self.pi_torsions, self.pi_torsion_types,
                           ['atom1', 'atom2', 'atom3', 'atom4', 'atom5', 'atom6'])
        copy_valence_terms(other.stretch_bends, other.stretch_bend_types, self.stretch_bends,
                           self.stretch_bend_types, ['atom1', 'atom2', 'atom3'])
        copy_valence_terms(other.torsion_torsions, other.torsion_torsion_types,
                           self.torsion_torsions, self.torsion_torsion_types,
                           ['atom1', 'atom2', 'atom3', 'atom4', 'atom5'])
        copy_valence_terms(other.chiral_frames, [], self.chiral_frames, [],
                           ['atom1', 'atom2', 'chirality'])
        copy_valence_terms(other.multipole_frames, [], self.multipole_frames, [],
                           ['atom', 'frame_pt_num', 'vectail', 'vechead', 'nvec'])
        copy_valence_terms(other.adjusts, other.adjust_types, self.adjusts,
                           self.adjust_types, ['atom1', 'atom2'])
        copy_valence_terms(other.donors, [], self.donors,
                           [], ['atom1', 'atom2'])
        copy_valence_terms(other.acceptors, [],
                           self.acceptors, [], ['atom1', 'atom2'])
        copy_valence_terms(other.groups, [], self.groups,
                           [], ['atom', 'type', 'move'])
        if mycrd is None or other._coordinates is None:
            self._coordinates = None
        else:
            ocrd = other.get_coordinates('all')
            nframes = min(ocrd.shape[0], mycrd.shape[0])
            self._coordinates = np.concatenate(
                (mycrd[:nframes, :, :], ocrd[:nframes, :, :]), axis=1)
        return self

    def __mul__(self, ncopies):
        """ Replicates the current Structure `ncopies` times """
        cp = copy(self)
        return cp.__imul__(ncopies, self)

    __rmul__ = __mul__

    def __imul__(self, ncopies, other=None):
        if not isinstance(ncopies, int):
            return NotImplemented
        # The basic approach here is similar to what we used in __iadd__, except
        # we don't have to extend the type arrays at all -- we just point to the
        # same one that the first copy pointed to.

        def copy_valence_terms(oval, aoffset, sval, styp, attrlist):
            """ Copies the valence terms from one list to another;
            oval=Other VALence; otyp=Other TYPe; sval=Self VALence;
            styp=Self TYPe; attrlist=ATTRibute LIST (atom1, atom2, ...)
            """
            for val in oval:
                ats = [getattr(val, attr) for attr in attrlist]
                # Replace any atoms with the ones from self.atoms
                for i, at in enumerate(ats):
                    if isinstance(at, Atom):
                        ats[i] = self.atoms[at.idx+aoffset]
                kws = dict()
                if hasattr(val, 'type') and val.type is NoUreyBradley:
                    kws['type'] = NoUreyBradley  # special-case singleton
                elif styp and val.type is not None:
                    kws['type'] = styp[val.type.idx]
                sval.append(type(val)(*ats, **kws))
                if hasattr(val, 'funct'):
                    sval[-1].funct = val.funct
        if other is None:
            other = copy(self)
        for i in range(ncopies-1):
            aoffset = len(self.atoms)
            roffset = self.residues[-1].number + 1
            for atom in other.atoms:
                res = atom.residue
                self.add_atom(copy(atom), res.name, res.idx+roffset, res.chain,
                              res.insertion_code, res.segid)
            copy_valence_terms(other.bonds, aoffset, self.bonds,
                               self.bond_types, ['atom1', 'atom2'])
            copy_valence_terms(other.angles, aoffset, self.angles,
                               self.angle_types, ['atom1', 'atom2', 'atom3'])
            copy_valence_terms(other.dihedrals, aoffset, self.dihedrals,
                               self.dihedral_types, ['atom1', 'atom2', 'atom3',
                                                     'atom4', 'improper', 'ignore_end'])
            copy_valence_terms(other.rb_torsions, aoffset,
                               self.rb_torsions, self.rb_torsion_types,
                               ['atom1', 'atom2', 'atom3', 'atom4', 'improper', 'ignore_end'])
            copy_valence_terms(other.urey_bradleys, aoffset,
                               self.urey_bradleys, self.urey_bradley_types, ['atom1', 'atom2'])
            copy_valence_terms(other.impropers, aoffset, self.impropers, self.improper_types,
                               ['atom1', 'atom2', 'atom3', 'atom4'])
            copy_valence_terms(other.cmaps, aoffset, self.cmaps, self.cmap_types,
                               ['atom1', 'atom2', 'atom3', 'atom4', 'atom5'])
            copy_valence_terms(other.trigonal_angles, aoffset, self.trigonal_angles,
                               self.trigonal_angle_types, ['atom1', 'atom2', 'atom3', 'atom4'])
            copy_valence_terms(other.out_of_plane_bends, aoffset, self.out_of_plane_bends,
                               self.out_of_plane_bend_types, ['atom1', 'atom2', 'atom3', 'atom4'])
            copy_valence_terms(other.pi_torsions, aoffset, self.pi_torsions, self.pi_torsion_types,
                               ['atom1', 'atom2', 'atom3', 'atom4', 'atom5', 'atom6'])
            copy_valence_terms(other.stretch_bends, aoffset, self.stretch_bends,
                               self.stretch_bend_types, ['atom1', 'atom2', 'atom3'])
            copy_valence_terms(other.torsion_torsions, aoffset,
                               self.torsion_torsions, self.torsion_torsion_types,
                               ['atom1', 'atom2', 'atom3', 'atom4', 'atom5'])
            copy_valence_terms(other.chiral_frames, aoffset, self.chiral_frames,
                               [], ['atom1', 'atom2', 'chirality'])
            copy_valence_terms(other.multipole_frames, aoffset, self.multipole_frames, [],
                               ['atom', 'frame_pt_num', 'vectail', 'vechead', 'nvec'])
            copy_valence_terms(other.adjusts, aoffset, self.adjusts,
                               self.adjust_types, ['atom1', 'atom2'])
            copy_valence_terms(other.donors, aoffset,
                               self.donors, [], ['atom1', 'atom2'])
            copy_valence_terms(other.acceptors, aoffset,
                               self.acceptors, [], ['atom1', 'atom2'])
            copy_valence_terms(other.groups, aoffset, self.groups, [], [
                               'atom', 'type', 'move'])
        if self._coordinates is not None:
            self._coordinates = np.hstack(
                [self._coordinates for i in range(ncopies)])
        return self

    # ===================================================

    # For the bool-ness of Structure. An empty structure evaluates to
    # boolean-false, but this means that the structure must have no atoms,
    # residues, topological features, or parameter types at all.

    def __bool__(self):
        return bool(self.atoms or self.residues or self.bonds or self.angles or self.dihedrals or
                    self.impropers or self.rb_torsions or self.cmaps or self.torsion_torsions or
                    self.stretch_bends or self.out_of_plane_bends or self.trigonal_angles or
                    self.torsion_torsions or self.pi_torsions or self.urey_bradleys or
                    self.chiral_frames or self.multipole_frames or self.adjusts or self.acceptors or
                    self.donors or self.groups or self.bond_types or self.angle_types or
                    self.dihedral_types or self.urey_bradley_types or self.improper_types or
                    self.rb_torsion_types or self.cmap_types or self.trigonal_angle_types or
                    self.out_of_plane_bend_types or self.pi_torsion_types or
                    self.torsion_torsion_types or self.adjust_types)

    __nonzero__ = __bool__  # for Python 2

    # ===================================================

    @staticmethod
    def _add_force_to_system(system, force):
        """ Adds an OpenMM force to a system IFF the force is not None """
        if force is None:
            return
        if isinstance(force, tuple) or isinstance(force, list):
            # It's possible we got multiple forces to add
            for f in force:
                system.addForce(f)
            return
        system.addForce(force)

    # ===================================================

    def __getstate__(self):
        """ Serializes a structure """
        retdict = dict(residues=self.residues,
                       bond_types=self.bond_types,
                       angle_types=self.angle_types,
                       dihedral_types=self.dihedral_types,
                       urey_bradley_types=self.urey_bradley_types,
                       improper_types=self.improper_types,
                       rb_torsion_types=self.rb_torsion_types,
                       cmap_types=self.cmap_types,
                       trigonal_angle_types=self.trigonal_angle_types,
                       out_of_plane_bend_types=self.out_of_plane_bend_types,
                       pi_torsion_types=self.pi_torsion_types,
                       stretch_bend_types=self.stretch_bend_types,
                       torsion_torsion_types=self.torsion_torsion_types,
                       adjust_types=self.adjust_types,
                       groups=self.groups,
                       _coordinates=self._coordinates,
                       _box=self._box,
                       nrexcl=self.nrexcl,
                       _combining_rule=self._combining_rule,
                       unknown_functional=self.unknown_functional,
                       space_group=self.space_group,
                       )

        def idx(thing):
            return thing.idx if thing is not None else None
        retdict['bonds'] = [(b.atom1.idx, b.atom2.idx, idx(b.type), b.order, b.qualitative_type)
                            for b in self.bonds]
        retdict['angles'] = [(a.atom1.idx, a.atom2.idx, a.atom3.idx,
                              idx(a.type)) for a in self.angles]
        retdict['dihedrals'] = [(d.atom1.idx, d.atom2.idx, d.atom3.idx, d.atom4.idx, d.improper,
                                 d.ignore_end, idx(d.type)) for d in self.dihedrals]
        retdict['impropers'] = [(d.atom1.idx, d.atom2.idx, d.atom3.idx, d.atom4.idx, idx(d.type))
                                for d in self.impropers]
        retdict['rb_torsions'] = [(d.atom1.idx, d.atom2.idx, d.atom3.idx, d.atom4.idx, idx(d.type))
                                  for d in self.rb_torsions]
        retdict['urey_bradleys'] = [(u.atom1.idx, u.atom2.idx, idx(u.type))
                                    for u in self.urey_bradleys]
        retdict['cmaps'] = [(c.atom1.idx, c.atom2.idx, c.atom3.idx, c.atom4.idx, c.atom5.idx,
                             idx(c.type)) for c in self.cmaps]
        retdict['trigonal_angles'] = [(t.atom1.idx, t.atom2.idx, t.atom3.idx, t.atom4.idx,
                                       idx(t.type)) for t in self.trigonal_angles]
        retdict['out_of_plane_bends'] = [(o.atom1.idx, o.atom2.idx, o.atom3.idx, o.atom4.idx,
                                          idx(o.type)) for o in self.out_of_plane_bends]
        retdict['pi_torsions'] = [(p.atom1.idx, p.atom2.idx, p.atom3.idx, p.atom4.idx, p.atom5.idx,
                                   p.atom6.idx, idx(p.type)) for p in self.pi_torsions]
        retdict['stretch_bends'] = [(s.atom1.idx, s.atom2.idx, s.atom3.idx,
                                     idx(s.type)) for s in self.stretch_bends]
        retdict['torsion_torsions'] = [(t.atom1.idx, t.atom2.idx, t.atom3.idx, t.atom4.idx,
                                        t.atom5.idx, idx(t.type)) for t in self.torsion_torsions]
        retdict['chiral_frames'] = [(c.atom1.idx, c.atom2.idx, c.chirality)
                                    for c in self.chiral_frames]
        retdict['multipole_frames'] = [(f.atom.idx, f.frame_pt_num, f.vectail,
                                        f.vechead, f.nvec) for f in self.multipole_frames]
        retdict['adjusts'] = [(e.atom1.idx, e.atom2.idx, idx(e.type))
                              for e in self.adjusts]
        retdict['acceptors'] = [(a.atom1.idx, a.atom2.idx)
                                for a in self.acceptors]
        retdict['donors'] = [(d.atom1.idx, d.atom2.idx) for d in self.donors]
        retdict['exclusions'] = [
            tuple(e.idx for e in a._exclusion_partners) for a in self.atoms]
        retdict['links'] = [(link.atom1.idx, link.atom2.idx, link.length,
                             link.symmetry_op1, link.symmetry_op2) for link in self.links]

        # Now the metadata stuff, if applicable
        for key in ('experimental', 'journal', 'authors', 'keywords', 'doi',
                    'pmid', 'journal_authors', 'volume', 'title', 'year',
                    'resolution', 'related_entries'):
            try:
                retdict[key] = getattr(self, key)
            except AttributeError:
                continue

        return retdict

    def __setstate__(self, d):
        # Assign the type lists we need to just copy over, and make sure the
        # attributes are claimed.
        for attr in ('residues', 'bond_types', 'angle_types', 'dihedral_types',
                     'urey_bradley_types', 'improper_types', 'rb_torsion_types',
                     'cmap_types', 'trigonal_angle_types', 'adjust_types',
                     'out_of_plane_bend_types', 'pi_torsion_types', 'groups',
                     'stretch_bend_types', 'torsion_torsion_types'):
            setattr(self, attr, d[attr])
            getattr(self, attr).claim()
        # Assign the possible metadata
        for key in ('experimental', 'journal', 'authors', 'keywords', 'doi',
                    'pmid', 'journal_authors', 'volume', 'title', 'year',
                    'resolution', 'related_entries', '_coordinates', '_box',
                    'nrexcl', '_combining_rule', 'unknown_functional',
                    'space_group'):
            if key in d:
                setattr(self, key, d[key])

        self.atoms = AtomList()
        for r in self.residues:
            self.atoms.extend(r.atoms)

        def assign_type(typelist, idx):
            if idx is None:
                return None
            return typelist[idx]

        def assign_order(bond_tuple):
            if len(bond_tuple) < 4:
                return 1.0
            return bond_tuple[3]

        def assign_qualitative_type(bond_tuple):
            if len(bond_tuple) < 5:
                return None
            qualitative_type = bond_tuple[4]
            return QualitativeBondType(qualitative_type) if qualitative_type is not None else None

        # Set the topology arrays
        self.bonds = TrackedList(
            Bond(
                self.atoms[it[0]], self.atoms[it[1]], type=assign_type(self.bond_types, it[2]),
                order=assign_order(it), qualitative_type=assign_qualitative_type(it),
            )
            for it in d['bonds']
        )
        self.angles = TrackedList(
            Angle(self.atoms[it[0]], self.atoms[it[1]], self.atoms[it[2]], type=assign_type(self.angle_types, it[3]))
            for it in d['angles']
        )
        self.dihedrals = TrackedList(
            Dihedral(self.atoms[it[0]], self.atoms[it[1]], self.atoms[it[2]], self.atoms[it[3]],
                     improper=it[4], ignore_end=it[5], type=assign_type(self.dihedral_types, it[6]))
            for it in d['dihedrals']
        )
        self.impropers = TrackedList(
            Improper(self.atoms[it[0]], self.atoms[it[1]], self.atoms[it[2]], self.atoms[it[3]],
                     type=assign_type(self.improper_types, it[4])) for it in d['impropers']
        )
        self.urey_bradleys = TrackedList(
            UreyBradley(self.atoms[it[0]], self.atoms[it[1]],
                        type=assign_type(self.urey_bradley_types, it[2]))
            for it in d['urey_bradleys']
        )
        self.rb_torsions = TrackedList(
            Dihedral(self.atoms[it[0]], self.atoms[it[1]], self.atoms[it[2]], self.atoms[it[3]],
                     type=assign_type(self.rb_torsion_types, it[4]))
            for it in d['rb_torsions']
        )
        self.cmaps = TrackedList(
            Cmap(self.atoms[it[0]], self.atoms[it[1]], self.atoms[it[2]], self.atoms[it[3]],
                 self.atoms[it[4]], type=assign_type(self.cmap_types, it[5]))
            for it in d['cmaps']
        )
        self.trigonal_angles = TrackedList(
            TrigonalAngle(self.atoms[it[0]], self.atoms[it[1]], self.atoms[it[2]],
                          self.atoms[it[3]], type=assign_type(self.trigonal_angle_types, it[4]))
            for it in d['trigonal_angles']
        )
        self.out_of_plane_bends = TrackedList(
            OutOfPlaneBend(self.atoms[it[0]], self.atoms[it[1]], self.atoms[it[2]],
                           self.atoms[it[3]], type=assign_type(self.out_of_plane_bend_types, it[4]))
            for it in d['out_of_plane_bends']
        )
        self.pi_torsions = TrackedList(
            PiTorsion(self.atoms[it[0]], self.atoms[it[1]], self.atoms[it[2]], self.atoms[it[3]],
                      self.atoms[it[4]], self.atoms[it[5]], type=assign_type(self.pi_torsion_types, it[6]))
            for it in d['pi_torsions']
        )
        self.stretch_bends = TrackedList(
            StretchBend(self.atoms[it[0]], self.atoms[it[1]], self.atoms[it[2]],
                        type=assign_type(self.stretch_bend_types, it[3]))
            for it in d['stretch_bends']
        )
        self.torsion_torsions = TrackedList(
            TorsionTorsion(self.atoms[it[0]], self.atoms[it[1]], self.atoms[it[2]],
                           self.atoms[it[3]], self.atoms[it[4]], type=assign_type(self.torsion_torsion_types, it[5]))
            for it in d['torsion_torsions']
        )
        self.chiral_frames = TrackedList(
            ChiralFrame(self.atoms[it[0]], self.atoms[it[1]], it[2]) for it in d['chiral_frames']
        )
        self.multipole_frames = TrackedList(
            MultipoleFrame(self.atoms[it[0]], *it[1:]) for it in d['multipole_frames']
        )
        self.adjusts = TrackedList(
            NonbondedException(self.atoms[it[0]], self.atoms[it[1]], assign_type(self.adjust_types, it[2]))
            for it in d['adjusts']
        )
        self.acceptors = TrackedList(
            AcceptorDonor(self.atoms[it[0]], self.atoms[it[1]]) for it in d['acceptors']
        )
        self.donors = TrackedList(
            AcceptorDonor(self.atoms[it[0]], self.atoms[it[1]]) for it in d['donors']
        )
        self.links = TrackedList(
            Link(self.atoms[it[0]], self.atoms[it[1]], it[2], it[3], it[4]) for it in d["links"]
        )
        # Transfer the exclusions
        for atom, excl in zip(self.atoms, d['exclusions']):
            for idx in excl:
                atom.exclude(self.atoms[idx])

    def _deduplicate_atom_type_names(self, other):
        """ De-duplicates atom type names from another structure. Useful if adding two Structures """
        other_atom_types = {
            atom.atom_type.name: atom.atom_type for atom in other.atoms if atom.atom_type is not UnassignedAtomType
        }
        copied_atom_type_map = dict()
        for atom in self.atoms:
            if atom.atom_type is UnassignedAtomType:
                continue
            if other_atom_types.get(atom.atom_type.name, atom.atom_type) != atom.atom_type:
                if id(atom.atom_type) in copied_atom_type_map:
                    new_atom_type = copied_atom_type_map[id(atom.atom_type)]
                else:
                    new_atom_type = copy(atom.atom_type)
                    new_atom_type.name += "D"
                    copied_atom_type_map[id(atom.atom_type)] = new_atom_type
                atom.type += "D"
                atom.atom_type = new_atom_type

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


class _StructureViewerCreator(object):
    """ Class responsible for creating a StructureView when indexed """

    def __init__(self, struct):
        self.struct = struct

    def __getitem__(self, selection):
        struct = self.struct
        if isinstance(selection, int):
            return struct.atoms[selection]

        view = StructureView()
        selection = struct._get_selection_array(selection)
        # _get_selection_array might have returned either an Atom or None if
        # there is no selection. Handle those and return, or continue if we got
        # a list of atoms
        if selection is None:
            return view
        elif isinstance(selection, Atom):
            return selection
        sumsel = sum(selection)
        if sumsel == 0:
            # No atoms selected. Return None
            return view
        # The cumulative sum of selection will give our index + 1 of each
        # selected atom into the new structure
        scan = [selection[0]]
        for i in range(1, len(selection)):
            scan.append(scan[i-1] + selection[i])
        # Zero-out the unselected atoms
        scan = [x * y for x, y in zip(scan, selection)]
        # Copy all parameters
        sel_res = set()
        for i, atom in enumerate(struct.atoms):
            if not selection[i]:
                continue
            view.atoms.append(atom)
            if atom.residue in sel_res:
                continue
            view.residues.append(atom.residue)
            sel_res.add(atom.residue)

        def add_valence_terms(oval, sval, attrlist):
            """ Adds the valence terms from one list to another;
            oval=Other VALence; otyp=Other TYPe; sval=Self VALence;
            styp=Self TYPe; attrlist=ATTRibute LIST (atom1, atom2, ...)
            """
            for val in sval:
                ats = [getattr(val, attr) for attr in attrlist]
                # Make sure all of our atoms in this valence term is "selected"
                indices = [scan[at.idx] for at in ats if isinstance(at, Atom)]
                if not all(indices):
                    continue
                # Add the type if applicable
                oval.append(val)
        add_valence_terms(view.bonds, struct.bonds, ['atom1', 'atom2'])
        add_valence_terms(view.angles, struct.angles,
                          ['atom1', 'atom2', 'atom3'])
        add_valence_terms(view.dihedrals, struct.dihedrals, [
                          'atom1', 'atom2', 'atom3', 'atom4'])
        add_valence_terms(view.rb_torsions, struct.rb_torsions,
                          ['atom1', 'atom2', 'atom3', 'atom4'])
        add_valence_terms(view.urey_bradleys,
                          struct.urey_bradleys, ['atom1', 'atom2'])
        add_valence_terms(view.impropers, struct.impropers, [
                          'atom1', 'atom2', 'atom3', 'atom4'])
        add_valence_terms(view.cmaps, struct.cmaps, [
                          'atom1', 'atom2', 'atom3', 'atom4', 'atom5'])
        add_valence_terms(view.trigonal_angles, struct.trigonal_angles,
                          ['atom1', 'atom2', 'atom3', 'atom4'])
        add_valence_terms(view.out_of_plane_bends, struct.out_of_plane_bends,
                          ['atom1', 'atom2', 'atom3', 'atom4'])
        add_valence_terms(view.pi_torsions, struct.pi_torsions,
                          ['atom1', 'atom2', 'atom3', 'atom4', 'atom5', 'atom6'])
        add_valence_terms(view.stretch_bends, struct.stretch_bends, [
                          'atom1', 'atom2', 'atom3'])
        add_valence_terms(view.torsion_torsions, struct.torsion_torsions,
                          ['atom1', 'atom2', 'atom3', 'atom4', 'atom5'])
        add_valence_terms(view.chiral_frames,
                          struct.chiral_frames, ['atom1', 'atom2'])
        add_valence_terms(view.multipole_frames,
                          struct.multipole_frames, ['atom'])
        add_valence_terms(view.adjusts, struct.adjusts, ['atom1', 'atom2'])
        add_valence_terms(view.donors, struct.donors, ['atom1', 'atom2'])
        add_valence_terms(view.acceptors, struct.acceptors, ['atom1', 'atom2'])
        return view

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


class StructureView(object):
    """
    A view of a Structure. In many cases, this can serve as a duck-typed
    Structure and it has many of the same attributes. However, none of its lists
    *own* their objects, and the lists of atoms, residues, and
    parameters/topologies are regular lists, rather than TrackedList instances.
    Therefore, the indexes correspond to the indexes from the *original*
    Structure from which this view was taken. Furthermore, there are no "type"
    lists, as they would be exactly equivalent to the type lists of the parent
    Structure instance.

    Attributes
    ----------
    atoms : :class:`AtomList`
        List of all atoms in the structure
    residues : :class:`ResidueList`
        List of all residues in the structure
    bonds : :class:`TrackedList` (:class:`Bond`)
        List of all bonds in the structure
    angles : :class:`TrackedList` (:class:`Angle`)
        List of all angles in the structure
    dihedrals : :class:`TrackedList` (:class:`Dihedral`)
        List of all dihedrals in the structure
    rb_torsions : :class:`TrackedList` (:class:`RBTorsion`)
        List of all Ryckaert-Bellemans torsions in the structure
    urey_bradleys : :class:`TrackedList` (:class:`UreyBradley`)
        List of all Urey-Bradley angle bends in the structure
    impropers : :class:`TrackedList` (:class:`Improper`)
        List of all CHARMM-style improper torsions in the structure
    cmaps : :class:`TrackedList` (:class:`Cmap`)
        List of all CMAP objects in the structure
    trigonal_angles : :class:`TrackedList` (:class:`TrigonalAngle`)
        List of all AMOEBA-style trigonal angles in the structure
    out_of_plane_bends : :class:`TrackedList` (:class:`OutOfPlaneBends`)
        List of all AMOEBA-style out-of-plane bending angles
    pi_torsions : :class:`TrackedList` (:class:`PiTorsion`)
        List of all AMOEBA-style pi-torsion angles
    stretch_bends : :class:`TrackedList` (:class:`StretchBend`)
        List of all AMOEBA-style stretch-bend compound bond/angle terms
    torsion_torsions : :class:`TrackedList` (:class:`TorsionTorsion`)
        List of all AMOEBA-style coupled torsion-torsion terms
    chiral_frames : :class:`TrackedList` (:class:`ChiralFrame`)
        List of all AMOEBA-style chiral frames defined in the structure
    multipole_frames : :class:`TrackedList` (:class:`MultipoleFrame`)
        List of all AMOEBA-style multipole frames defined in the structure
    adjusts : :class:`TrackedList` (:class:`NonbondedException`)
        List of all nonbonded pair-exception rules
    acceptors : :class:`TrackedList` (:class:`AcceptorDonor`)
        List of all H-bond acceptors, if that information is present
    donors : :class:`TrackedList` (:class:`AcceptorDonor`)
        List of all H-bond donors, if that information is present
    positions : u.Quantity(list(Vec3), u.angstroms)
        Unit-bearing atomic coordinates. If not all atoms have coordinates, this
        property is None
    coordinates : np.ndarray of shape (nframes, natom, 3)
        If no coordinates are set, this is set to None. The first frame will
        match the coordinates present on the atoms.
    """

    def __init__(self):
        self.atoms = []
        self.residues = []
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.rb_torsions = []
        self.urey_bradleys = []
        self.impropers = []
        self.cmaps = []
        self.trigonal_angles = []
        self.out_of_plane_bends = []
        self.pi_torsions = []
        self.stretch_bends = []
        self.torsion_torsions = []
        self.chiral_frames = []
        self.multipole_frames = []
        self.adjusts = []
        self.acceptors = []
        self.donors = []

    def to_dataframe(self):
        """ Generates a DataFrame from the current Structure's atomic properties

        Returns
        -------
        df : DataFrame
            DataFrame with all atomic properties

        See Also
        --------
        :func:`parmed.utils.pandautils.create_dataframe`
        """
        from .utils.pandautils import create_dataframe
        return create_dataframe(self)

    def load_dataframe(self, df):
        """ Loads atomic properties from an input DataFrame

        Parameters
        ----------
        df : pandas.DataFrame
            A pandas DataFrame with atomic properties that will be used to set
            the properties on the current list of atoms

        See Also
        --------
        :func:`parmed.utils.pandautils.load_dataframe`
        """
        from .utils.pandautils import load_dataframe
        return load_dataframe(self, df)

    @property
    def coordinates(self):
        try:
            return np.array([[a.xx, a.xy, a.xz] for a in self.atoms])
        except AttributeError:
            return None

    @property
    def positions(self):
        try:
            return [Vec3(a.xx, a.xy, a.xz) for a in self.atoms] * u.angstroms
        except AttributeError:
            return None

    def __bool__(self):
        return bool(self.atoms or self.residues or self.bonds or self.angles or
                    self.dihedrals or self.impropers or self.rb_torsions or
                    self.cmaps or self.torsion_torsions or self.stretch_bends or
                    self.out_of_plane_bends or self.trigonal_angles or
                    self.torsion_torsions or self.pi_torsions or
                    self.urey_bradleys or self.chiral_frames or
                    self.multipole_frames or self.adjusts or self.acceptors or
                    self.donors)
        # Ignore types here because it's just a view

    def __repr__(self):
        natom = len(self.atoms)
        nres = len(self.residues)
        nextra = sum([isinstance(a, ExtraPoint) for a in self.atoms])
        retstr = ['<%s %d atoms' % (type(self).__name__, natom)]
        if nextra > 0:
            retstr.append(' [%d EPs]' % nextra)
        retstr.append('; %d residues' % nres)
        nbond = len(self.bonds)
        retstr.append('; %d bonds>' % nbond)
        return ''.join(retstr)

    __nonzero__ = __bool__  # For Python 2

    def __iter__(self):
        return iter(self.atoms)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


def _settler(parm):
    """ Identifies residues that can have SETTLE applied to it """
    is_water = []
    for r in parm.residues:
        na = sum(1 for a in r if not isinstance(a, ExtraPoint))
        if na != 3:
            is_water.append(False)
            continue
        # Make sure we are not bonded to any other residues
        for a in r:
            for a2 in a.bond_partners:
                if a2.residue is not r:
                    is_water.append(False)
                    break
            else:
                # this atom is OK
                continue
            # We only hit here if we hit the break above. So break here
            break
        else:
            # Don't check elements, since they may not *always* be accurate.
            # It's more likely that a 3-atom residue not bonded to any other
            # residue is water.
            is_water.append(True)
    assert len(is_water) == len(parm.residues), 'Incorrect length'
    return is_water


def _res_in_templlib(res, lib):
    """ Returns the residue template inside lib that matches res """
    if res.name in lib:
        return lib[res.name]
    if len(res.name) == 3 and residue.AminoAcidResidue.has(res.name):
        return lib[residue.AminoAcidResidue.get(res.name).abbr]
    if residue.DNAResidue.has(res.name):
        return lib[residue.DNAResidue.get(res.name).abbr]
    if residue.RNAResidue.has(res.name) and residue.RNAResidue.get(res.name).abbr != 'T':
        return lib[residue.RNAResidue.get(res.name).abbr]
    # Not present
    return None
