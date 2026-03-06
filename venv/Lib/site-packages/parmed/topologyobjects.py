"""
This module contains objects that deal with system topology and parameters, such
as atoms, residues, bonds, angles, etc.

by Jason Swails
"""
import enum
import math
from functools import cached_property
from typing import Optional
import warnings
from abc import ABC, abstractmethod
from copy import copy
from functools import wraps

from . import unit as u
from .constants import TINY_DIGITS as _TINY_DIGITS, DEG_TO_RAD, RAD_TO_DEG, TINY
from .exceptions import MoleculeError, ParameterError, ParameterWarning
from .geometry import angle, dihedral, distance2
from .periodic_table import Element

__all__ = ['Angle', 'AngleType', 'Atom', 'AtomList', 'Bond', 'BondType', 'ChiralFrame', 'Cmap',
           'CmapType', 'Dihedral', 'DihedralType', 'DihedralTypeList', 'Improper', 'ImproperType',
           'MultipoleFrame', 'OutOfPlaneBend', 'PiTorsion', 'Residue', 'ResidueList', 'StretchBend',
           'StretchBendType', 'TorsionTorsion', 'TorsionTorsionType', 'TrigonalAngle',
           'TrackedList', 'UreyBradley', 'OutOfPlaneBendType', 'NonbondedException',
           'NonbondedExceptionType', 'AmoebaNonbondedExceptionType', 'AcceptorDonor', 'Group',
           'AtomType', 'NoUreyBradley', 'ExtraPoint', 'TwoParticleExtraPointFrame',
           'ThreeParticleExtraPointFrame', 'OutOfPlaneExtraPointFrame', 'RBTorsionType',
           'UnassignedAtomType', 'Link', 'DrudeAtom', 'DrudeAnisotropy', 'Hybridization', 'QualitativeBondType']

def _strip_units(value, unit=None):
    """
    Strips any units from the given value by casting them into the AKMA unit
    system (or the requested unit)
    """
    if u.is_quantity(value):
        # special-case angles, since pure angles are always in degrees
        if unit is None:
            if value.unit.is_compatible(u.degrees):
                return value.value_in_unit(u.degrees)
            return value.value_in_unit_system(u.akma_unit_system)
        else:
            return value.value_in_unit(unit)
    return value

def _exception_to_notimplemented(func):
    """
    Wraps comparison operators to return NotImplemented instead of raising an AttributeError
    """
    @wraps(func)
    def wrapper(self, other):
        try:
            return func(self, other)
        except AttributeError:
            return NotImplemented
    return wrapper

def _getstate_with_exclusions(exclusions=None):
    """ Serializes based on all attributes except requested exclusions

    Parameters
    ----------
    exclusions : list of str, optional
        List of all attributes to exclude from serialization (should be
        descriptors and 'list'). Default is None

    Notes
    -----
    If exclusions is None, it defaults to excluding 'list'. If this is not
    desired, set exclusions to the empty list.
    """
    if exclusions is None:
        exclusions = ['list']
    def __getstate__(self):
        return {key : val for (key, val) in self.__dict__.items() if key not in exclusions}
    return __getstate__

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Private classes and methods

class _ListItem:
    """
    Helpful methods for items that appear in some kind of tracked list.
    Subclasses of this method interact with their `list` attribute, if they have
    one, in order to determine _where_ in the list they occur.

    Attributes
    ----------
    idx : ``int``
        This is intended to be a read-only variable that determines where in the
        list this particular object is. If there is no `list` attribute for this
        object, or the item is not in the list at all, `idx` is -1

    Notes
    -----
    For lists that support indexing its members and tracking when that list is
    changed (so we know when to update indexes), this is a fast lookup, as
    indexing only needs to be done once, at most, for the entire list (until
    changes are made that could change the indexing).

    The ``idx`` lookup in a TrackedList is therefore an O(1) operation, while
    it is O(N) for standard containers (O(N^2) for looking up *all* indexes).
    """

    @property
    def idx(self):
        try:
            mylist = self.list
        except AttributeError:
            return -1

        if mylist is None:
            return -1

        try:
            needs_indexing = mylist.needs_indexing
        except AttributeError:
            # This isn't a tracked list, so just look through the list
            for i, item in enumerate(mylist):
                if item is self: return i
            return -1
        else:
            if needs_indexing or self._idx == -1:
                try:
                    self._idx = -1
                    mylist.index_members()
                    return self._idx
                except AttributeError:
                    for i, item in enumerate(mylist):
                        if item is self: return i
                    return -1
            else:
                return self._idx

class _FourAtomTerm:
    """
    A base class for a parameter that spans 4 atoms

    Parameters
    ----------
    atom1 : :class:`Atom`
        The first atom in the term
    atom2 : :class:`Atom`
        The second atom in the term
    atom3 : :class:`Atom`
        The third atom in the term
    atom4 : :class:`Atom`
        The fourth atom in the term

    Notes
    -----
    This is a base class that should be overridden by terms consisting of 4
    atoms (like torsions). There are a lot of such terms in the Amoeba force
    field, so this functionality is broken out into a new form
    """
    def __init__(self, atom1, atom2, atom3, atom4):
        if len({atom1, atom2, atom3, atom4}) != 4:
            raise MoleculeError(
                f'4-atom term cannot have duplicate atoms! Atoms are: {atom1} {atom2} {atom3} {atom4}'
            )
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4

    def __contains__(self, obj):
        return self.atom1 is obj or self.atom2 is obj or self.atom3 is obj or self.atom4 is obj

    def delete(self):
        """ Sets all atoms in this term to None, and its type if it exists """
        self.atom1 = self.atom2 = self.atom3 = self.atom4 = None
        if hasattr(self, 'type'):
            self.type = None

class _ParameterType:
    """
    A parameter type that defines the nature of a particular molecular
    interaction. This is a base class that simply indicates whether a particular
    parameter type is "used", which in turn is used to determine whether or not
    this parameter type needs to be printed

    Attributes
    ----------
    used : ``bool``
        If ``True``, then this parameter type should be considered *used*. If
        ``False``, it is not being used and does not need to be printed.
    penalty : float or None
        If this is assigned from a database, there might be a penalty assigned
        to the determination of this parameter
    """

    def __init__(self):
        self.used = False
        self.penalty = None

    def __ne__(self, other):
        return not self == other

def _delete_from_list(list, item):
    """
    Deletes a requested item from a list. If the item does not exist in the
    list, a ValueError is raised

    Parameters
    ----------
    list : ``list``
        The list from which an item will be deleted
    item : ``object``
        The object to delete from the list
    """
    list.pop(list.index(item))

def _safe_assigns(dest, source, attrs):
    """
    Shallow-copies all requested attributes from `source` to `dest` if they
    are present. If not present, nothing is done
    """
    for attr in attrs:
        if not hasattr(source, attr): continue
        myattr = getattr(source, attr)
        setattr(dest, attr, myattr)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Hybridization(enum.IntEnum):
    """This represents atom hybridization, and mirrors the values in RDKit
    See https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html#rdkit.Chem.rdchem.HybridizationType
    """
    UNSPECIFIED = 0
    S = 1
    SP = 2
    SP2 = 3
    SP3 = 4
    SP2D = 5
    SP3D = 6
    SP3D2 = 7
    OTHER = 8


class QualitativeBondType(enum.IntEnum):
    """These are bond types as understood in the cheminformatics sense (e.g., single, double, etc.).
    The definitions here adopt the RDKit values seen at https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html#rdkit.Chem.rdchem.BondType
    """
    UNSPECIFIED = 0
    SINGLE = 1
    DOUBLE = 2
    TRIPLE = 3
    QUADRUPLE = 4
    QUINTUPLE = 5
    HEXTUPLE = 6
    ONEANDAHALF = 7
    TWOANDAHALF = 8
    THREEANDAHALF = 9
    FOURANDAHALF = 10
    FIVEANDAHALF = 11
    AROMATIC = 12
    IONIC = 13
    HYDROGEN = 14
    THREECENTER = 15
    DATIVEONE = 16
    DATIVE = 17
    DATIVEL = 18
    DATIVER = 19
    OTHER = 20
    ZERO = 21

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Atom(_ListItem):
    """
    An atom. Only use these as elements in AtomList instances, since AtomList
    will keep track of when indexes and other stuff needs to be updated. All
    parameters are optional.

    Parameters
    ----------
    atomic_number : ``int``
        The atomic number of this atom
    name : ``str``
        The name of this atom
    type : ``str``
        The type name of this atom
    charge : ``float``
        The partial atomic charge of this atom in fractions of an electron
    mass : ``float``
        The atomic mass of this atom in daltons
    nb_idx : ``int``
        The nonbonded index. This is a pointer that is relevant in the context
        of an Amber topology file and identifies its Lennard-Jones atom type
    solvent_radius : ``float``
        The intrinsic solvation radius of this atom.
    screen : ``float``
        The Generalized Born screening factor for this atom.
    occupancy : ``float``
        The occupancy of the atom (see PDB file)
    bfactor : ``float``
        The B-factor of the atom (see PDB file)
    altloc : ``str``
        Alternate location indicator (see PDB file)

    Other Parameters
    ----------------
    list : :class:`AtomList`
        The AtomList that this atom belongs to. If None, this atom does not
        belong to any list. This can be any iterable, but should typically be an
        AtomList or None
    tree : ``str``
        The tree chain identifier assigned to this atom. Relevant in the context
        of an Amber topology file, and not used for very much.
    join : ``int``
        The 'join` property of atoms stored in the Amber topology file. At the
        time of writing this class, `join` is unused, but still oddly required.
        Add support for future-proofing
    irotat : ``int``
        The `irotat` property of atoms stored in the Amber topology file.
        Unused, but included for future-proofing.
    number : ``int``
        The serial number given to the atom (see PDB file)
    rmin : ``float``
        The Rmin/2 Lennard-Jones parameter for this atom. Default evaluates to 0
    epsilon : ``float``
        The epsilon (well depth) Lennard-Jones parameter for this atom. Default
        evaluates to 0
    rmin14 : ``float``
        The Rmin/2 Lennard-Jones parameter for this atom in 1-4 interactions.
        Default evaluates to 0
    epsilon14 : ``float``
        The epsilon (well depth) Lennard-Jones parameter for this atom in 1-4
        interactions. Default evaluates to 0
    formal_charge : ``int``
        The formal charge of this atom. The default value is `None`
    hybridization : ``HybridizationType``
        The hybridization type of this atom. Default is ``None``.
    aromatic : ``bool``
        If set, indicates whether or not the atom is aromatic

    Other Attributes
    ----------------
    element : ``int``
        This is an alias for atomic_number
    atom_type : :class:`AtomType`
        In some cases, "type" is an ambiguous choice of an integer serial number
        or a string descriptor. In this case, atom_type is an AtomType instance
        that disambiguates this discrepancy.
    anisou : numpy.ndarray(float64) (or list of floats)
        Anisotropic temperature scaling factors. This is a 6-element numpy array
        They are the 3x3 symmetric matrix elements U(1,1), U(2,2), U(3,3),
        U(1,2), U(1,3), U(2,3). If no factors available, it is None.
    idx : ``int``
        The index of this atom in the list. Set to -1 if this atom is not part
        of a list or the index cannot otherwise be determined (i.e., if the
        containing list does not support indexing its members)
    residue : :class:`Residue`
        The Residue that this atom belongs to. This is assigned when this atom
        is passed to `Residue.add_atom` -- see below for more information. Until
        it is set there, it is None
    other_locations : ``dict`` of :class:`Atom`
        A dict of Atom instances that represent alternate conformers of this
        atom. The keys are the `altloc` characters for those Atoms.
    bonds : ``list`` of :class:`Bond`
        list of Bond objects in which this atom is a member. This attribute
        should not be modified.
    angles : ``list`` of :class:`Angle`
        list of Angle objects in which this atom is a member. This attribute
        should not be modified.
    dihedrals : ``list`` of :class:`Dihedral`
        list of Dihedral objects in which this atom is a member. This attribute
        should not be modified.
    urey_bradleys : ``list`` of :class:`UreyBradley`
        list of UreyBradley objects in which this atom is a member (CHARMM,
        AMOEBA). This attribute should not be modified.
    impropers : ``list`` of :class:`Improper`
        list of Improper objects in which the atom is a member (CHARMM). This
        attribute should not be modified.
    cmaps : ``list`` of :class:`Cmap`
        list of Cmap objects in which the atom is a member (CHARMM, AMOEBA).
        This attribute should not be modified.
    tortors : ``list`` of :class:`TorsionTorsion`
        list of TorsionTorsion objects in which the atom is a member (AMOEBA).
        This attribute should not be modified.
    bond_partners : ``list`` of :class:`Atom`
        list of Atoms to which this atom is bonded. Do not modify this
        attribute -- it will almost certainly not do what you think it will
    angle_partners : ``list`` of :class:`Atom`
        list of Atoms to which this atom forms an angle, but not a bond. Do not
        modify this attribute -- it will almost certainly not do what you think
        it will
    dihedral_partners : ``list`` of :class:`Atom`
        list of Atoms to which this atom forms an dihedral, but not a bond or
        angle. Do not modify this attribute -- it will almost certainly not do
        what you think it will
    tortor_partners : ``list`` of :class:`Atom`
        list of Atoms to which this atom forms a coupled Torsion-Torsion, but
        not a bond or angle (AMOEBA). Do not modify this attribute -- it will
        almost certainly not do what you think it will
    exclusion_partners : ``list`` of :class:`Atom`
        list of Atoms with which this atom is excluded, but not bonded, angled,
        or dihedraled to. Do not modify this attribute -- it will almost
        certainly not do what you think it will
    marked : ``int``
        Mainly for internal use, it is used to indicate when certain atoms have
        been "marked" when traversing the bond network identifying topological
        features (like molecules and rings)
    children : ``list`` of :class:`ExtraPoint`
        This is the list of "child" ExtraPoint objects bonded to this atom
    number : ``int``
        The serial number of the atom in the input structure (e.g., PDB file).
        If not present in the original structure, a default value of -1 is used
    rmin : ``float``
        The Rmin/2 Lennard-Jones parameter for this atom. Default value is 0.
        If not set, it is taken from the `atom_type` attribute (if available)
    epsilon : ``float``
        The epsilon (well depth) Lennard-Jones parameter for this atom. Default
        value is 0. If not set, it is taken from the `atom_type` attribute (if
        available)
    rmin_14 : ``float``
        The Rmin/2 L-J parameter for 1-4 pairs. Default value is `rmin` (see
        above). If not set, it is taken from the `atom_type` attribute (if
        available).
    epsilon_14 : ``float``
        The epsilon L-J parameter for 1-4 pairs. Default value is `epsilon` (see
        above). If not set, it is taken from the `atom_type` attribute (if
        available).

    Possible Attributes
    -------------------
    xx : ``float``
        The X-component of the position of this atom. Only present if
        coordinates have been loaded. Otherwise raises AttributeError
    xy : ``float``
        The Y-component of the position of this atom. Only present if
        coordinates have been loaded. Otherwise raises AttributeError
    xz : ``float``
        The Z-component of the position of this atom. Only present if
        coordinates have been loaded. Otherwise raises AttributeError
    vx : ``float``
        The X-component of the velocity of this atom. Only present if the
        velocities have been loaded. Otherwise raises AttributeError
    vy : ``float``
        The Y-component of the velocity of this atom. Only present if the
        velocities have been loaded. Otherwise raises AttributeError
    vz : ``float``
        The Z-component of the velocity of this atom. Only present if the
        velocities have been loaded. Otherwise raises AttributeError
    type_idx : ``int``
        The AMOEBA atom type index. Only present if initialized with the AMOEBA
        force field. Otherwise raises AttributeError.
    class_idx : ``int``
        The AMOEBA class type index Only present if initialized with the AMOEBA
        force field. Otherwise raises AttributeError.
    multipoles : ``list(float)``
        The list of the 10 multipole moments up through quadrupoles Only present
        if initialized with the AMOEBA force field. Otherwise raises
        AttributeError.
    polarizability : ``list(float)``
        The polarizability of the atom. Only present if initialized with the
        AMOEBA force field. Otherwise raises AttributeError.
    vdw_parent : :class:`Atom`
        In the AMOEBA force field, this is the parent atom for the van der Waals
        term
    vdw_weight : ``float``
        In the AMOEBA force field, this is the weight of the van der Waals
        interaction on the parent atom

    Notes
    -----
    The bond_partners, angle_partners, dihedral_partners, and exclusion_partners
    arrays are actually generated as properties by taking differences of sets
    and sorting them. As a result, accessing this attribute constantly can be
    less efficient than you would expect. Iterating over them in a loop requires
    minimal overhead. But if frequent access is needed and these sets are
    guaranteed not to change, you should save a reference to the object and use
    that instead.

    Binary comparisons are done by atom index and are used primarily for sorting
    sublists of atoms. The == operator is not defined, so Atom equality should
    be done using the `is` operator. This allows Atom instances to be hashable
    (and so used as `dict` keys and put in `set`s)

    Examples
    --------
    >>> a1 = Atom(name='CO', type='C', charge=0.5, mass=12.01)
    >>> a2 = Atom(name='OC', type='O', charge=-0.5, mass=12.01)
    >>> a1.bond_to(a2)
    >>> a1 in a2.bond_partners and a2 in a1.bond_partners
    True
    >>> a1.idx # Not part of a container
    -1

    This object also supports automatic indexing when it is part of a container

    >>> atom_list = []
    >>> atom_list.append(Atom(list=atom_list, name='CO', charge=0.5))
    >>> atom_list.append(Atom(list=atom_list, name='OC', charge=-0.5))
    >>> atom_list[0].idx
    0
    >>> atom_list[1].idx
    1
    """
    #===================================================

    def __init__(self, list=None, atomic_number=0, name='', type='',
                 charge=None, mass=0.0, nb_idx=0, solvent_radius=0.0,
                 screen=0.0, tree='BLA', join=0.0, irotat=0.0, occupancy=0.0,
                 bfactor=0.0, altloc='', number=-1, rmin=None, epsilon=None,
                 rmin14=None, epsilon14=None, formal_charge=None, hybridization=None, aromatic=None):
        self.list = list
        self._idx = -1
        self.atomic_number = atomic_number
        self.name = name.strip()
        try:
            self.type = type.strip()
        except AttributeError:
            self.type = type
        self._charge = _strip_units(charge, u.elementary_charge)
        self.mass = _strip_units(mass, u.dalton)
        self.nb_idx = nb_idx
        self.solvent_radius = _strip_units(solvent_radius, u.angstrom)
        self.screen = screen
        self.tree = tree
        self.join = join
        self.irotat = irotat
        self.bfactor = bfactor
        self.altloc = altloc
        self.occupancy = occupancy
        self._bond_partners = []
        self._angle_partners = []
        self._dihedral_partners = []
        self._tortor_partners = []
        self._exclusion_partners = [] # For arbitrary/other exclusions
        self.residue = None
        self.marked = 0 # For setting molecules
        self.bonds, self.angles, self.dihedrals = [], [], []
        self.urey_bradleys, self.impropers, self.cmaps = [], [], []
        self.tortors = []
        self.other_locations = {} # A dict of Atom instances
        self.atom_type = UnassignedAtomType
        self.number = number
        self.anisou = None
        self._rmin = _strip_units(rmin, u.angstroms)
        self._epsilon = _strip_units(epsilon, u.kilocalories_per_mole)
        self._rmin14 = _strip_units(rmin14, u.angstroms)
        self._epsilon14 = _strip_units(epsilon14, u.kilocalories_per_mole)
        self.children = []
        self.aromatic = aromatic
        self.formal_charge = formal_charge
        self.hybridization = Hybridization(hybridization) if hybridization is not None else None

    #===================================================

    @classmethod
    def _copy(cls, item):
        new = cls(atomic_number=item.atomic_number, name=item.name, type=item.type,
                  charge=item.charge, mass=item.mass, nb_idx=item.nb_idx,
                  solvent_radius=item.solvent_radius, screen=item.screen, tree=item.tree,
                  join=item.join, irotat=item.irotat, occupancy=item.occupancy,
                  bfactor=item.bfactor, altloc=item.altloc, formal_charge=item.formal_charge,
                  hybridization=item.hybridization, aromatic=item.aromatic)
        new.atom_type = item.atom_type
        new.anisou = copy(item.anisou)
        for key in item.other_locations:
            new.other_locations[key] = copy(item.other_locations[key])
        _safe_assigns(new, item, ('xx', 'xy', 'xz', 'vx', 'vy', 'vz', 'type_idx', 'class_idx',
                                  'multipoles', 'polarizability', 'vdw_parent', 'vdw_weight'))
        return new

    def __copy__(self):
        """ Returns a deep copy of this atom, but not attached to any list """
        return type(self)._copy(self)

    #===================================================

    @property
    def bond_partners(self):
        """ Go through all bonded partners """
        bp = set(self._bond_partners)
        for p in self._bond_partners:
            for c in p.children:
                bp.add(c)
        return sorted(list(bp))

    @property
    def angle_partners(self):
        """ List of all angle partners that are NOT bond partners """
        bp = set(self._bond_partners)
        ap = set(self._angle_partners) - bp
        toadd = set()
        for p in ap:
            for c in p.children:
                toadd.add(c)
        ap |= toadd
        return sorted(list(ap))

    @property
    def dihedral_partners(self):
        " List of all dihedral partners that are NOT angle or bond partners "
        bp = set(self._bond_partners)
        ap = set(self._angle_partners)
        dp = set(self._dihedral_partners) - ap - bp
        toadd = set()
        for p in dp:
            for c in p.children:
                toadd.add(c)
        dp |= toadd
        return sorted(list(dp))

    @property
    def tortor_partners(self):
        """
        List of all 1-5 partners that are NOT in angle or bond partners. This is
        *only* used in the Amoeba force field
        """
        bp = set(self._bond_partners)
        ap = set(self._angle_partners)
        dp = set(self._dihedral_partners)
        tp = set(self._tortor_partners) - dp - ap - bp
        toadd = set()
        for p in tp:
            for c in p.children:
                toadd.add(c)
        tp |= toadd
        return sorted(list(tp))

    @property
    def exclusion_partners(self):
        """
        List of all exclusions not otherwise excluded by bonds/angles/torsions
        """
        # A little expensive, but the only way to ensure this is completely
        # correct easily
        bp = set(self.bond_partners)
        ap = set(self.angle_partners)
        dp = set(self.dihedral_partners)
        tp = set(self.tortor_partners)
        ep = set(self._exclusion_partners) - tp - dp - ap - bp
        toadd = set()
        for p in ep:
            for c in p.children:
                toadd.add(c)
        ep |= toadd
        return sorted(list(ep))

    #===================================================

    # Various parameters that can be taken from the AtomType if not set on the
    # atom directly.

    @property
    def charge(self):
        if self._charge is None:
            if self.atom_type is UnassignedAtomType or self.atom_type.charge is None:
                return 0.0
            return self.atom_type.charge
        return self._charge

    @charge.setter
    def charge(self, value):
        self._charge = _strip_units(value, unit=u.elementary_charge)

    @property
    def ucharge(self):
        """ Charge with units """
        return self.charge * u.elementary_charge

    @property
    def rmin(self):
        """ Lennard-Jones Rmin/2 parameter (the Lennard-Jones radius) """
        if self._rmin is None:
            if self.atom_type is UnassignedAtomType or self.atom_type.rmin is None:
                return 0.0
            return self.atom_type.rmin
        return self._rmin

    @rmin.setter
    def rmin(self, value):
        """ Lennard-Jones Rmin/2 parameter (the Lennard-Jones radius) """
        self._rmin = _strip_units(value, unit=u.angstroms)

    @property
    def urmin(self):
        """ Lennard-Jones Rmin/2 parameter with units """
        return self.rmin * u.angstrom

    @property
    def sigma(self):
        """ Lennard-Jones sigma parameter -- directly related to Rmin """
        return self.rmin * 2**(-1/6) * 2

    @sigma.setter
    def sigma(self, value):
        self._rmin = _strip_units(value, unit=u.angstroms) * 2**(1/6) / 2

    @property
    def usigma(self):
        """ Lennard-Jones sigma parameter with units """
        return self.sigma * u.angstroms

    @property
    def epsilon(self):
        """ Lennard-Jones epsilon parameter (the Lennard-Jones well depth) """
        if self._epsilon is None:
            if self.atom_type is UnassignedAtomType or self.atom_type.epsilon is None:
                return 0.0
            return self.atom_type.epsilon
        return self._epsilon

    @epsilon.setter
    def epsilon(self, value):
        """ Lennard-Jones epsilon parameter (the Lennard-Jones well depth) """
        self._epsilon = _strip_units(value, unit=u.kilocalorie_per_mole)

    @property
    def uepsilon(self):
        """ Lennard-Jones epsilon parameter with units """
        return self.epsilon * u.kilocalories_per_mole

    @property
    def rmin_14(self):
        """ The 1-4 Lennard-Jones Rmin/2 parameter """
        if self._rmin14 is None:
            if self.atom_type is UnassignedAtomType or self.atom_type.rmin_14 is None:
                return self.rmin
            return self.atom_type.rmin_14
        return self._rmin14

    @rmin_14.setter
    def rmin_14(self, value):
        """ The 1-4 Lennard-Jones Rmin/2 parameter """
        self._rmin14 = _strip_units(value, u.angstroms)

    @property
    def urmin_14(self):
        """ The 1-4 Lennard-Jones Rmin/2 parameter with units """
        return self.rmin_14 * u.angstroms

    @property
    def sigma_14(self):
        """ Lennard-Jones sigma parameter -- directly related to Rmin """
        if self._rmin14 is None:
            if self.atom_type is UnassignedAtomType or self.atom_type.rmin_14 is None:
                return self.sigma
            return self.atom_type.rmin_14 * 2**(-1/6) * 2
        return self._rmin14 * 2**(-1/6) * 2

    @sigma_14.setter
    def sigma_14(self, value):
        self._rmin14 = _strip_units(value, u.angstroms) * 2**(1/6) / 2

    @property
    def usigma_14(self):
        """ The 1-4 Lennard-Jones sigma parameter with units """
        return self.sigma_14 * u.angstroms

    @property
    def epsilon_14(self):
        """ The 1-4 Lennard-Jones epsilon parameter """
        if self._epsilon14 is None:
            if self.atom_type is UnassignedAtomType or self.atom_type.epsilon_14 is None:
                return self.epsilon
            return self.atom_type.epsilon_14
        return self._epsilon14

    @epsilon_14.setter
    def epsilon_14(self, value):
        """ The 1-4 Lennard-Jones epsilon parameter """
        self._epsilon14 = _strip_units(value, u.kilocalories_per_mole)

    @property
    def uepsilon_14(self):
        """ The 1-4 Lennard-Jones epsilon parameter """
        return self.epsilon_14 * u.kilocalories_per_mole

    @property
    def umass(self):
        return self.mass * u.daltons

    @property
    def usolvent_radius(self):
        """ Solvation radius with units attached """
        return self.solvent_radius * u.angstroms

    #===================================================

    def nonbonded_exclusions(self, only_greater=True, index_from=0):
        """
        Returns the total number of nonbonded atom exclusions for this atom. The
        general rules for building the exclusion list is to include both
        exceptions AND exclusions (i.e., the Amber scaling of 1-4 interactions
        means that the 1-4 terms are excluded and a special pairlist is built to
        handle those exceptions).

        All atoms in the `_partners` arrays are nonbonded exclusions.

        Parameters
        ----------
        only_greater : ``bool``, optional
            If True (default), only atoms whose `idx` value is greater than this
            `Atom`s `idx` will be counted as an exclusion (to avoid double-
            counting exclusions). If False, all exclusions will be counted.
        index_from : ``int``, optional
            This is the index of the first atom, and is intended to be 0 (for C-
            and Python-style numbering, default) or 1 (for Fortran-style
            numbering, such as that used in the Amber and CHARMM topology files)

        Returns
        -------
        ``list of int``
            The returned list will be the atom indexes of the exclusion partners
            for this atom (indexing starts from ``index_from``)

        Notes
        -----
        If this instance's `idx` attribute evaluates to -1 -- meaning it is not
        in an AtomList -- a IndexError will be raised. If you have two extra
        points (i.e., those with atomic numbers of 0) bonded to each other, this
        routine may raise a ``RuntimeError`` if the recursion depth is exceeded.
        """
        if self.idx < 0:
            raise IndexError('Cannot find exclusions of an unindexed Atom')
        if only_greater:
            baseline = self.idx
        else:
            baseline = -1
        excl = []
        for atm in self.bond_partners:
            i = atm.idx + index_from
            if i > baseline:
                excl.append(i)
        for atm in self.angle_partners:
            i = atm.idx + index_from
            if i > baseline:
                excl.append(i)
        for atm in self.dihedral_partners:
            i = atm.idx + index_from
            if i > baseline:
                excl.append(i)
        for atm in self.tortor_partners:
            i = atm.idx + index_from
            if i > baseline:
                excl.append(i)
        for atm in self.exclusion_partners:
            i = atm.idx + index_from
            if i > baseline:
                excl.append(i)
        return sorted(excl)

    #===================================================

    # Make 'element' an alias for 'atomic_number'

    @property
    def element(self):
        return self.atomic_number
    @element.setter
    def element(self, value):
        self.atomic_number = value
    @property
    def element_name(self):
        return Element[self.atomic_number]

    #===================================================

    def bond_to(self, other):
        """
        Log this atom as bonded to another atom.

        Parameters
        ----------
        other : :class:`Atom`
            An atom that will be added to `bond_partners`

        Notes
        -----
        This action adds `self` to `other.bond_partners`. Raises
        :class:`MoleculeError` if `other is self`
        """
        if isinstance(other, ExtraPoint):
            self.children.append(other)
        elif isinstance(self, ExtraPoint):
            other.children.append(self)
        if self is other:
            raise MoleculeError(f"Cannot bond atom to itself! Atoms are: {self} {other}")
        self._bond_partners.append(other)
        other._bond_partners.append(self)

    #===================================================

    def angle_to(self, other):
        """
        Log this atom as angled to another atom.

        Parameters
        ----------
        other : :class:`Atom`
            An atom that will be added to `angle_partners`

        Notes
        -----
        This action adds `self` to `other.angle_partners`. Raises
        :class:`MoleculeError` if `other is self`
        """
        if self is other:
            raise MoleculeError(f"Cannot angle an atom with itself! Atoms are: {self} {other}")
        self._angle_partners.append(other)
        other._angle_partners.append(self)

    #===================================================

    def dihedral_to(self, other):
        """
        Log this atom as dihedral-ed to another atom.

        Parameters
        ----------
        other : :class:`Atom`
            An atom that will be added to `dihedral_partners`

        Notes
        -----
        This action adds `self` to `other.dihedral_partners`. Raises
        :class:`MoleculeError` if `other is self`
        """
        if self is other:
            raise MoleculeError(f"Cannot dihedral an atom with itself! Atoms are: {self} {other}")
        self._dihedral_partners.append(other)
        other._dihedral_partners.append(self)

    #===================================================

    def tortor_to(self, other):
        """
        Log this atom as 1-5 partners to another atom

        Parameters
        ----------
        other : :class:`Atom`
            An atom that will be added to `tortor_partners`

        Notes
        -----
        This action adds `self` to `other.tortor_partners`. Raises
        :class:`MoleculeError` if `other is self`
        """
        if self is other:
            raise MoleculeError(f'Cannot coupled-dihedral atom to itself Atoms are: {self} {other}')

        self._tortor_partners.append(other)
        other._tortor_partners.append(self)

    #===================================================

    def exclude(self, other):
        """
        Add one atom to my arbitrary exclusion list

        Parameters
        ----------
        other : :class:`Atom`
            An atom that will be added to `exclusion_partners`.

        Notes
        -----
        This action adds `self` to `other.exclusion_partners`. Raises
        :class:`MoleculeError` if `other is self`
        """
        if self is other:
            raise MoleculeError(f"Cannot exclude an atom from itself! Atoms are: {self} {other}")
        self._exclusion_partners.append(other)
        other._exclusion_partners.append(self)

    #===================================================

    def prune_exclusions(self):
        """
        For extra points, the exclusion partners may be filled before the bond,
        angle, dihedral, and tortor partners. Since we don't want memory of
        these exclusions if any of those topological features were to break, we
        want to *remove* those from the exclusion list. This function makes sure
        that nothing in the bond, angle, dihedral, and tortor lists appears in
        the exclusion list.
        """
        excludes = (set(self._exclusion_partners) - set(self._tortor_partners) -
                    set(self._dihedral_partners) - set(self._angle_partners) -
                    set(self._bond_partners))
        self._exclusion_partners = sorted(list(excludes))

    #===================================================

    # Comparisons are done by comparing indexes

    def __gt__(self, other):
        return self.idx > other.idx

    def __lt__(self, other):
        return self.idx < other.idx

    def __ge__(self, other):
        return not self < other

    def __le__(self, other):
        return not self > other

    def __repr__(self):
        start = f'<{self.__class__.__name__} {self.name} [{self.idx}]'
        if self.residue is not None and hasattr(self.residue, 'idx'):
            return start + f'; In {self.residue.name} {self.residue.idx}>'
        elif self.residue is not None:
            return start + f'; In object {repr(self.residue)}>'
        return start + '>'

    #===================================================

    # For pickleability

    def __getstate__(self):
        retval = dict(name=self.name, type=self.type, atom_type=self.atom_type,
                      _charge=self._charge, mass=self.mass, nb_idx=self.nb_idx,
                      solvent_radius=self.solvent_radius, screen=self.screen,
                      tree=self.tree, join=self.join, irotat=self.irotat, bfactor=self.bfactor,
                      altloc=self.altloc, occupancy=self.occupancy, number=self.number,
                      anisou=self.anisou, _rmin=self._rmin, _epsilon=self._epsilon,
                      _rmin14=self._rmin14, _epsilon14=self._epsilon14, children=self.children,
                      atomic_number=self.atomic_number, formal_charge=self.formal_charge,
                      hybridization=self.hybridization, aromatic=self.aromatic)
        for key in ('xx', 'xy', 'xz', 'vx', 'vy', 'vz', 'multipoles', 'type_idx', 'class_idx',
                    'polarizability', 'vdw_weight', 'weights', '_frame_type'):
            try:
                retval[key] = getattr(self, key)
            except AttributeError:
                continue

        return retval

    def __setstate__(self, d):
        self._bond_partners = []
        self._angle_partners = []
        self._dihedral_partners = []
        self._tortor_partners = []
        self._exclusion_partners = []
        self.residue = None
        self.marked = 0
        self.bonds, self.angles, self.dihedrals = [], [], []
        self.urey_bradleys, self.impropers, self.cmaps = [], [], []
        self.tortors = []
        self.other_locations = {}
        self.__dict__.update(d)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ExtraPoint(Atom):
    """
    An extra point is a massless, virtual site that is used to extend the
    flexibility of partial atomic charge fitting. They can be used to, for
    instance, give atoms permanent multipoles in a fixed-charge format.

    However, virtual sites are massless particles whose position is fixed with
    respect to a particular frame of reference. As a result, they must be
    treated specially when running dynamics. This class extends the `Atom` class
    with extra functionality specific to these "Extra points" or "virtual
    sites". See the documentation for the :class:`Atom` class for more
    information.

    Parameters
    ----------
    weights : list of float, optional, keyword-only
        This is the list of weights defining its frame.

    See Also
    --------
    :class:`Atom`: The ExtraPoint constructor also takes all `Atom` arguments
                   as well, and shares all of the same properties
    """
    def __init__(self, *args, **kwargs):
        if 'weights' in kwargs:
            self.weights = kwargs.pop('weights')
        else:
            self.weights = None
        super().__init__(*args, **kwargs)
        self._frame_type = None

    #===================================================

    @property
    def parent(self):
        """
        The parent atom of this extra point is the atom it is bonded to (there
        should only be 1 bond partner, but the parent is defined as its first)
        """
        try:
            return self._bond_partners[0]
        except IndexError:
            return None

    #===================================================

    @property
    def bond_partners(self):
        """ List of all atoms bonded to this atom and its parent """
        try:
            return sorted([self.parent] + [x for x in self.parent.bond_partners if x is not self])
        except AttributeError:
            return []

    @property
    def angle_partners(self):
        """ List of all atoms angled to this atom and its parent """
        try:
            return self.parent.angle_partners
        except AttributeError:
            return []

    @property
    def dihedral_partners(self):
        try:
            return self.parent.dihedral_partners
        except AttributeError:
            return []

    @property
    def tortor_partners(self):
        try:
            return self.parent.tortor_partners
        except AttributeError:
            return []

    @property
    def exclusion_partners(self):
        try:
            return self.parent.exclusion_partners
        except AttributeError:
            return []

    #===================================================

    @property
    def frame_type(self):
        """
        The type of frame used for this extra point. Whether the EP lies
        beyond (outside) or between (inside) the is determined from the
        positions of each atom in the frame and the extra point. If those are
        not set, the EP is assumed to be between the two points in the frame
        """
        if self._frame_type is not None:
            return self._frame_type
        if len(self.parent.bonds) == 2:
            mybond = None
            otherbond = None
            for bond in self.parent.bonds:
                if self in bond:
                    mybond = bond
                else:
                    otherbond = bond
            assert mybond is not None and otherbond is not None, 'Strange bond pattern detected'
            bonddist = otherbond.type.req
            if otherbond.atom1 is self.parent:
                otheratom = otherbond.atom2
            else:
                otheratom = otherbond.atom1
            try:
                x1, y1, z1 = self.xx, self.xy, self.xz
                x2, y2, z2 = otheratom.xx, otheratom.xy, otheratom.xz
            except AttributeError:
                self._frame_type = TwoParticleExtraPointFrame(self, True)
            else:
                dx, dy, dz = x1-x2, y1-y2, z1-z2
                if dx*dx + dy*dy + dz*dz > bonddist:
                    self._frame_type = TwoParticleExtraPointFrame(self, False)
                else:
                    self._frame_type = TwoParticleExtraPointFrame(self, True)
            return self._frame_type
        elif len(self.parent.bonds) == 3:
            # Just take 1 other atom -- an acute angle is "inside", an obtuse
            # angle is "outside"
            other_atom = None
            for bond in self.parent.bonds:
                if not self in bond:
                    if bond.atom1 is self.parent:
                        other_atom = bond.atom2
                    else:
                        other_atom = bond.atom1
                    break
            assert other_atom is not None, 'Strange bond pattern detected'
            try:
                x1, y1, z1 = other_atom.xx, other_atom.xy, other_atom.xz
                x2, y2, z2 = self.parent.xx, self.parent.xy, self.parent.xz
                x3, y3, z3 = self.xx, self.xy, self.xz
            except AttributeError:
                # See if both other atoms are hydrogens and the current atom is
                # an oxygen. If so, this is TIP4P and we go inside. Otherwise,
                # we go outside
                other_atom2 = None
                for bond in self.parent.bonds:
                    if not self in bond and not other_atom in bond:
                        if bond.atom1 is self.parent:
                            other_atom2 = bond.atom2
                        else:
                            other_atom2 = bond.atom1
                        break
                if other_atom2 is None:
                    raise RuntimeError('Strange bond pattern detected')
                if (self.parent.element == 8 and other_atom.element == 1 and
                        other_atom2.element == 1): # TIP4P!
                    self._frame_type = ThreeParticleExtraPointFrame(self, True)
                else:
                    self._frame_type = ThreeParticleExtraPointFrame(self, False)
            else:
                vec1 = [x1-x2, y1-y2, z1-z2]
                vec2 = [x3-x2, y3-y2, z3-z2]
                dot11 = sum([x*y for x, y in zip(vec1, vec1)])
                dot22 = sum([x*y for x, y in zip(vec2, vec2)])
                dot12 = sum([x*y for x, y in zip(vec1, vec2)])
                angle = math.acos(dot12 / (math.sqrt(dot11)*math.sqrt(dot22)))
                if angle < math.pi / 2:
                    self._frame_type = ThreeParticleExtraPointFrame(self, True)
                else:
                    self._frame_type = ThreeParticleExtraPointFrame(self, False)
        elif len(self.parent.bonds) == 4:
            # Only supported option here is the OOP one (e.g., TIP5P) -- always
            # assume tetrahedral angle
            self._frame_type = OutOfPlaneExtraPointFrame(self)

        return self._frame_type

    @frame_type.setter
    def frame_type(self, value):
        self._frame_type = value

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ExtraPointFrame(ABC):

    @abstractmethod
    def get_weights(self):
        pass

    @abstractmethod
    def get_atoms(self):
        pass

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class TwoParticleExtraPointFrame(ExtraPointFrame):
    r"""
    This class defines a frame of reference for a given extra point with a frame
    of reference defined by 2 particles

    Parameters
    ----------
    ep : :class:`ExtraPoint`
        The extra point defined by this frame
    inside : ``bool``
        If True, the extra point is contained inside the curve connecting all
        points in the frame of reference. If False, the point is outside that
        curve. See figures below for example

    Figures
    -------
    In the figures, x marks the real particles, e marks the extra point
    (properly numbered). The numbers refer to the ordering this class expects
    those atoms to be in. The real particle that is the parent of the extra
    point is shown in upper-case

    Inside  : x1-----e--X2

    Outside : x1--------X2--e
    """
    def __init__(self, ep, inside=False):
        # Don't do error-checking now, since we want to give it time for the
        # bonds to be generated -- only error check when trying to get weights
        self.ep = ep
        self.inside = inside

    def get_atoms(self):
        """
        Returns the particles involved in the frame

        Returns
        -------
        a1, a2 : :class:`Atom`, :class:`Atom`
            a1 is the parent atom of the extra point. a2 is the other atom
            bonded to the parent atom
        """
        try:
            mybond, = [bond for bond in self.ep.parent.bonds if self.ep not in bond]
        except (ValueError, TypeError):
            raise RuntimeError("Bad bond pattern in EP frame")

        if mybond.atom1 is self.ep.parent:
            return self.ep.parent, mybond.atom2
        else:
            return self.ep.parent, mybond.atom1

    def get_weights(self):
        """
        Returns the weights for the two particles

        Returns
        -------
        w1, w2 : ``float, float``
            w1 is the weight of particle 1 (the parent atom), whereas w2 is the
            weight of particle 2 (the atom bonded to the parent atom)
        """
        ep = self.ep
        if len(ep.parent.bonds) != 2:
            raise ValueError('EP parent bond pattern inconsistent with 2-point virtual site frame')
        b1, b2 = ep.parent.bonds
        if ep in b1:
            r1 = b2.type.req
            r2 = b1.type.req
        else:
            r1 = b1.type.req
            r2 = b2.type.req

        if self.inside:
            # It is closer to atom 1, but both weights are positive and add to 1
            return ((r1 - r2) / r1), (r2 / r1)
        else:
            return ((r1 + r2) / r1), -(r2 / r1)

class ThreeParticleExtraPointFrame(ExtraPointFrame):
    r"""
    This class defines a frame of reference for a given extra point with a frame
    of reference defined by 3 particles

    Parameters
    ----------
    ep : :class:`ExtraPoint`
        The extra point defined by this frame
    inside : ``bool``
        If True, the extra point is contained inside the curve connecting all
        points in the frame of reference. If False, the point is outside that
        curve. See figures below for example

    Figures
    -------
    In the figures, x marks the real particles, e marks the extra point
    (properly numbered). The numbers refer to the ordering this class expects
    those atoms to be in. The real particle that is the parent of the extra
    point is shown in upper-case

    Inside :         X
                    /|\
                   / e \
                  x     x
    Outside :
                     e
                     |
                     X
                    / \
                   /   \
                  x     x
    """
    def __init__(self, ep, inside=True):
        # Don't do error-checking now, since we want to give it time for the
        # bonds to be generated -- only error check when trying to get weights
        self.ep = ep
        self.inside = inside

    def get_atoms(self):
        """
        Returns the particles involved in the frame

        Returns
        -------
        a1, a2, a3 : :class:`Atom`, :class:`Atom`, :class:`Atom`
            a1 is the parent atom of the extra point. a2 and a3 are the other
            atoms that are both bonded to the parent atom
        """
        try:
            b1, b2 = [bond for bond in self.ep.parent.bonds if self.ep not in bond]
        except (ValueError, TypeError):
            raise RuntimeError('Unsupported bonding pattern in EP frame')
        if b1.atom1 is self.ep.parent:
            oatom1 = b1.atom2
        else:
            oatom1 = b1.atom1
        if b2.atom1 is self.ep.parent:
            oatom2 = b2.atom2
        else:
            oatom2 = b2.atom1
        return self.ep.parent, oatom1, oatom2

    @staticmethod
    def from_weights(parent, a1, a2, w1, w2, dp1=None, dp2=None, theteq=None, d12=None):
        """
        This function determines the necessary bond length between an ExtraPoint
        and its parent atom from the weights that are calculated in
        ``get_weights``.

        Parameters
        ----------
        parent : :class:`Atom`
            The parent atom to the ExtraPoint
        a1 : :class:`Atom`
            The first atom in the frame bonded to the parent
        a2 : :class:`Atom`
            The second atom in the frame bonded to the parent
        w1 : float
            The first weight defining the ExtraPoint position wrt ``a1``
        w2 : float
            The second weight defining the ExtraPoint position wrt ``a2``
        dp1 : float, optional
            Equilibrium distance between parent and a1. If None, a bond with a
            bond type must be defined between parent and a1, and the distance
            will be taken from there. Units must be Angstroms. Default is None
        dp2 : float, optional
            Same as dp1 above, but between atoms parent and a2. Either both dp1
            and dp2 should be None or neither should be. If one is None, the
            other will be ignored if not None.
        theteq : float, optional
            Angle between bonds parent-a1 and parent-a2. If None, d12 (below)
            will be used to calculate the angle. If both are None, the angle
            or d12 distance will be calculated from parametrized bonds or
            angles between a1, parent, and a2. Units of degrees. Default None.
        d12 : float, optional
            Distance between a1 and a2 as an alternative way to specify theteq.
            Default is None.

        Returns
        -------
        dist : float
            The distance between the ExtraPoint and its parent atom that
            satisfies the weights

        Notes
        -----
        parent must form a Bond with both a1 and a2.  Then, if a1-parent-a2
        forms an angle and it has an assigned type, that equilibrium value is
        used. Otherwise, the a1-a2 bond distance is used.  If neither of those
        are defined, a ValueError is raised.

        Raises
        ------
        ValueError if the necessary geometry requirements are not set or if the
        two weights are different
        """
        if a1 not in parent.bond_partners or a2 not in parent.bond_partners:
            raise ValueError('Parent atom not bound to other 2 atoms in frame')
        if a2 not in a1.angle_partners and a2 not in a1.bond_partners:
            raise ValueError('No geometry defined between 2 non-parent atoms')
        if w1 != w2:
            raise ValueError('Currently only equal weights are supported')
        # distance(OV) = w1 * cos(EP-Parent-a1) * parent-a1 dist
        if dp1 is None or dp2 is None:
            for bond in parent.bonds:
                if a1 in bond:
                    if bond.type is None:
                        raise ParameterError('Could not determine virtual site geometry')
                    dp1 = bond.type.req
                if a2 in bond:
                    if bond.type is None:
                        raise ParameterError('Could not determine virtual site geometry')
                    dp2 = bond.type.req
        if theteq is None and d12 is None:
            if a2 not in a1.angle_partners:
                for bond in a1.bonds:
                    if a2 not in bond: continue
                    if bond.type is None:
                        raise ParameterError('Could not determine virtual site geometry')
                    d12 = bond.type.req
                # Get angle from law of cosines
                theteq = math.acos((dp1*dp1+dp2*dp2-d12*d12)/(2*dp1*dp2))
            else:
                for ang in a1.angles:
                    if a2 in ang and a2 is not ang.atom2: #TODO test angle.type is None
                        theteq = ang.type.theteq * DEG_TO_RAD
                        break
                else:
                    assert False, "Could not find matching angle"
        elif theteq is None:
            theteq = math.acos((dp1*dp1+dp2*dp2-d12*d12)/(2*dp1*dp2))
        else:
            theteq *= DEG_TO_RAD

        if abs(dp1 - dp2) > TINY:
            raise ValueError('Cannot deal with asymmetry in EP frame')

        return abs(w1 * 2 * math.cos(theteq * 0.5) * dp1)

    def get_weights(self):
        """
        Returns the weights for the three particles

        Returns
        -------
        w1, w2, w3 : ``float, float, float``
            w1 is the weight of particle 1 (the parent atom), whereas w2 and w3
            are the weights of particles 2 and 3 (the atoms bonded to the
            parent atom)
        """
        ep = self.ep
        if len(ep.parent.bonds) != 3:
            raise ValueError('EP parent bond pattern inconsistent with 3-point virtual site frame')
        b1, b2, b3 = ep.parent.bonds
        # There are 2 possibilities here -- there is an angle between the 3
        # atoms in the frame OR there is a triangle of bonds (e.g., TIP4P).
        # Compute the 'ideal' distance between the 3 frame-of-ref. atoms
        if ep in b1:
            b1, b3 = b3, b1
        elif ep in b2:
            b2, b3 = b3, b2
        # See if there is an angle with both b1 and b2 in it
        found = False
        for ang in ep.parent.angles:
            if b1 in ang and b2 in ang:
                found = True
                break
        if found:
            # Compute the 2-3 distance from the two bond lengths and the angles
            # using law of cosines
            r1 = b1.type.req
            r2 = b2.type.req
            theta = ang.type.theteq * DEG_TO_RAD
            req23 = math.sqrt(r1*r1 + r2*r2 - 2*r1*r2*math.cos(theta))
        else:
            # See if there is a bond between particles 2 and 3
            if b1.atom1 is ep.parent:
                a1 = b1.atom2
            else:
                a1 = b1.atom1
            if b2.atom1 is ep.parent:
                a2 = b2.atom2
            else:
                a2 = b2.atom1
            if a1 not in a2.bond_partners:
                raise RuntimeError('EP frame definition incomplete for 3-point virtual site... '
                                   'cannot determine distance between particles 2 and 3')
            req23 = None
            for bond in a1.bonds:
                if a2 in bond:
                    req23 = bond.type.req
            if req23 is None:
                raise RuntimeError('Cannot determine 2-3 particle distance in '
                                   'three-site virtual site frame')
        req12 = b1.type.req
        req13 = b2.type.req
        weight = b3.type.req / math.sqrt(req12*req13 - 0.25*req23*req23)

        if self.inside:
            return 1 - weight, weight / 2, weight / 2
        else:
            return 1 + weight, -weight / 2, -weight / 2

class OutOfPlaneExtraPointFrame(ExtraPointFrame):
    r"""
    This class defines a frame of reference for a given extra point with a frame
    of reference defined by 3 particles, but with the virtual site out of the
    plane of those 3 particles. For example, TIP5P

    Parameters
    ----------
    ep : :class:`ExtraPoint`
        The extra point defined by this frame
    angle : ``float``
        The angle out-of-plane that the extra point is. By default, it is half
        of the angle of a tetrahedron. Given in degrees

    Figures
    -------
    In the figures, x marks the real particles, e marks the extra point
    (properly numbered). The numbers refer to the ordering this class expects
    those atoms to be in. The real particle that is the parent of the extra
    point is shown in upper-case

                   e   e
                    \ /
                     X
                    / \
                   /   \
                  x     x
    """
    def __init__(self, ep, angle=54.735):
        # Don't do error-checking now, since we want to give it time for the
        # bonds to be generated -- only error check when trying to get weights
        self.ep = ep
        self.angle = angle

    def get_atoms(self):
        """
        Returns the particles involved in the frame

        Returns
        -------
        a1, a2, a3 : :class:`Atom`, :class:`Atom`, :class:`Atom`
            a1 is the parent atom of the extra point. a2 and a3 are the other
            atoms that are both bonded to the parent atom but are not EPs
        """
        try:
            b1, b2 = [bond for bond in self.ep.parent.bonds
                      if (not isinstance(bond.atom1, ExtraPoint) and
                          not isinstance(bond.atom2, ExtraPoint))
                     ]
        except (ValueError, TypeError):
            raise RuntimeError('Unsupported bonding pattern in EP frame')
        if b1.atom1 is self.ep.parent:
            oatom1 = b1.atom2
        else:
            oatom1 = b1.atom1
        if b2.atom1 is self.ep.parent:
            oatom2 = b2.atom2
        else:
            oatom2 = b2.atom1
        return self.ep.parent, oatom2, oatom1

    def get_weights(self):
        """
        Returns the weights for the three particles

        Returns
        -------
        w1, w2, w3 : ``float, float, float``
            w1 and w2 are the weights with respect to the second two particles
            in the frame (i.e., NOT the parent atom). w3 is the weight of the
            cross-product
        """
        ep = self.ep
        # Find the two bonds that do not involve any extra points
        regbonds = []
        mybond = None
        for bond in ep.parent.bonds:
            if (not isinstance(bond.atom1, ExtraPoint) and not isinstance(bond.atom2, ExtraPoint)):
                regbonds.append(bond)
            elif ep in bond:
                if mybond is not None:
                    raise ValueError('multiple bonds detected to extra point')
                mybond = bond
        #FIXME -- is this necessary?
        if len(regbonds) != 2:
            raise ValueError('EP parent bond pattern inconsistent with an '
                             'out-of-plane, 3-point virtual site frame')
        if mybond is None:
            raise RuntimeError('No EP bond found... should not be here')
        b1, b2 = regbonds[:2]
        req12 = b1.type.req
        req13 = b2.type.req
        # See if there is an angle with both b1 and b2 in it
        found = False
        for ang in ep.parent.angles:
            if b1 in ang and b2 in ang:
                found = True
                break
        if found:
            # Compute the 2-3 distance from the two bond lengths and the angles
            # using law of cosines
            t213 = ang.theteq
            r1 = b1.type.req
            r2 = b2.type.req
            theta = ang.type.theteq * DEG_TO_RAD
            req23 = math.sqrt(r1*r1 + r2*r2 - 2*r1*r2*math.cos(theta))
        else:
            # See if there is a bond between particles 2 and 3
            if b1.atom1 is ep.parent:
                a1 = b1.atom2
            else:
                a1 = b1.atom1
            if b2.atom1 is ep.parent:
                a2 = b2.atom2
            else:
                a2 = b2.atom1
            if a1 not in a2.bond_partners:
                raise RuntimeError('EP frame definition incomplete for 3-point virtual site... '
                                   'cannot determine distance between particles 2 and 3')
            req23 = None
            for bond in a1.bonds:
                if a2 in bond:
                    req23 = bond.type.req
            if req23 is None:
                raise RuntimeError('Cannot determine 2-3 particle distance in '
                                   'three-site virtual site frame')
            # Now calculate the angle
            t213 = 2 * math.asin(req23 / (req12 + req13)) * RAD_TO_DEG
        # Some necessary constants
        length_conv = u.angstroms.conversion_factor_to(u.nanometers)
        req12 *= length_conv
        req13 *= length_conv
        req23 *= length_conv
        sinOOP = math.sin(self.angle * DEG_TO_RAD)
        cosOOP = math.cos(self.angle * DEG_TO_RAD)
        sin213 = math.sin(t213 * DEG_TO_RAD)
        # Find how big the cross product is
        lenCross = req12 * req13 * sin213
        # Find our weights (assume symmetric)
        weightCross = sinOOP * mybond.type.req * length_conv / lenCross
        weight = cosOOP * mybond.type.req * length_conv / math.sqrt(req12*req13 - 0.25*req23*req23)
        # Find out of the cross product weight should be positive or negative.
        a1, a2, a3 = self.get_atoms()
        try:
            v12 = (a2.xx-a1.xx, a2.xy-a1.xy, a2.xz-a1.xz)
            v13 = (a3.xx-a1.xx, a3.xy-a1.xy, a3.xz-a1.xz)
            v1e = (self.ep.xx-a1.xx, self.ep.xy-a1.xy, self.ep.xz-a1.xz)
        except AttributeError:
            # No coordinates... we have to guess. The first EP will have a
            # positive weight, the second will have a negative (this matches
            # what happens with Amber prmtop files...)
            for a in self.ep.parent.bond_partners:
                if a is self.ep:
                    break
                if isinstance(a, ExtraPoint):
                    weightCross = -weightCross
                    break
        else:
            # Take the cross product of v12 and v13, then dot that with the EP
            # vector. An acute angle is a positive weightCross. An obtuse one is
            # negative
            cross = (v12[1]*v13[2] - v12[2]*v13[1],
                     v12[2]*v13[0] - v12[0]*v13[2],
                     v12[0]*v13[1] - v12[1]*v13[0])
            lencross = math.sqrt(sum([cross[i]*cross[i] for i in range(3)]))
            lenv1e = math.sqrt(sum([v1e[i]*v1e[i] for i in range(3)]))
            v1edotcross = sum([v1e[i]*cross[i] for i in range(3)])
            costheta = v1edotcross / (lenv1e*lencross)
            weightCross = -weightCross if costheta < 0 else weightCross
        return weight / 2, weight / 2, weightCross

class LocalCoordinatesFrame(ExtraPointFrame):
    """
    Frame based on the local coordinates

    Parameters
    ----------
    atom1: :class:`Atom`
        The first atom defining the frame
    atom2: :class:`Atom`
        The second atom defining the frame
    atom3: :class:`Atom`
        The third atom defining the frame
    origin_weights: list
        weights (must sum to 1) for the three particles when computing the origin location
    x_weights: list
        the weight factors for the three particles when computing x-direction
    y_weights: list
        the weight factors for the three particles when computing y-direction
    distance: float
        the distance in the local coordinate frame
    angle: float
        the angle in the local coordinate frame
    dihedral: float
        the dihedral angle in the local coordinate frame
    frame_size: int
        The number of immediately bonded atoms defining the frame

    Notes
    -----
    All of the lists must have 3 elements.

    Raises
    ------
    ValueError if any of the weights or position are invalid
    """
    def __init__(self, atom1, atom2, atom3, origin_weights, x_weights, y_weights, distance, angle, dihedral, frame_size):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.origin_weights = self._validate_weights(origin_weights, lambda x: abs(sum(x) - 1) > 1e-5)
        self.x_weights = self._validate_weights(x_weights)
        self.y_weights = self._validate_weights(y_weights)
        self.distance = _strip_units(distance, u.angstrom)
        self.angle = _strip_units(angle, u.degrees)
        self.dihedral = _strip_units(dihedral, u.degrees)
        self.frame_size = frame_size

    @staticmethod
    def _validate_weights(weights, does_violate=None) -> list:
        my_weights = list(weights)
        if len(my_weights) != 3:
            raise ValueError(f"Weights must have 3 elements, not {len(my_weights)}")
        if does_violate is not None and does_violate(weights):
            raise ValueError("Weights failed condition requirement")
        return [float(wt) for wt in my_weights]

    def get_weights(self):
        return self.origin_weights, self.x_weights, self.y_weights, self.local_position

    def get_atoms(self):
        return self.atom1, self.atom2, self.atom3

    @cached_property
    def local_position(self):
        def small_to_zero(num: float) -> float:
            if abs(num) < 1e-10:
                return 0.0
            return num
        angle = self.angle * DEG_TO_RAD
        dihedral = (180 - self.dihedral) * DEG_TO_RAD
        distance = abs(self.distance)
        part = distance * math.sin(angle)
        return [
            small_to_zero(distance * math.cos(angle)),
            small_to_zero(part * math.cos(dihedral)),
            small_to_zero(part * math.sin(dihedral)),
        ]

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Bond:
    """
    A covalent bond connecting two atoms.

    Parameters
    ----------
    atom1 : :class:`Atom`
        The first atom involved in the bond
    atom2 : :class:`Atom`
        The other atom involved in the bond
    type : :class:`BondType` or None, optional
        The bond type that defines the parameters for this bond. Default is None
    order : float, optional
        The bond order of this bond. Bonds are classified as follows:
            1.0 -- single bond
            2.0 -- double bond
            3.0 -- triple bond
            1.5 -- aromatic bond
            1.25 -- amide bond
        Default is 1.0
    qualitative_type : :class:`QualitativeBondType`
        An optional qualitative bond type

    Notes
    -----
    You can test whether an :class:`Atom` is contained within the bond using the
    `in` operator. A `MoleculeError` is raised if `atom1` and `atom2` are identical.
    This bond instance is `append`ed to the `bonds` list for both `atom1` and
    `atom2` and is automatically removed from those lists upon garbage
    collection

    Examples
    --------
    >>> a1, a2 = Atom(), Atom()
    >>> bond = Bond(a1, a2)
    >>> a1 in bond and a2 in bond
    True
    """

    def __init__(self, atom1, atom2, type=None, order=1.0, qualitative_type=None):
        """ Bond constructor """
        # Make sure we're not bonding me to myself
        if atom1 is atom2:
            raise MoleculeError(f'Cannot bond atom to itself! Atoms are: {atom1} {atom2}')
        if isinstance(atom1, ExtraPoint) and isinstance(atom2, ExtraPoint):
            raise MoleculeError('Cannot bond two virtual sites/extra points together')
        # Order the atoms so the lowest atom # is first
        self.atom1 = atom1
        self.atom2 = atom2
        # Load this struct into the atoms
        self.atom1.bonds.append(self)
        self.atom2.bonds.append(self)
        atom1.bond_to(atom2)
        self.type = type
        self.funct = 1
        self._order = None
        self.order = order
        self.qualitative_type = QualitativeBondType(qualitative_type) if qualitative_type is not None else None

    def __contains__(self, thing):
        """ Quick and easy way to see if an Atom is in this Bond """
        return thing is self.atom1 or thing is self.atom2

    def delete(self):
        """
        Deletes this bond from the atoms that make it up. This method removes
        this bond from the `bonds` list for both atom1 and atom2 as well as
        removing atom1 from atom2.bond_partners and vice versa.
        """
        _delete_from_list(self.atom1.bonds, self)
        _delete_from_list(self.atom2.bonds, self)
        _delete_from_list(self.atom1._bond_partners, self.atom2)
        _delete_from_list(self.atom2._bond_partners, self.atom1)

        self.atom1 = self.atom2 = self.type = None

    @property
    def order(self):
        """ Bond order. See description in :class:`Bond` argument list """
        return self._order

    @order.setter
    def order(self, value):
        """
        Order of the bond. Must be a float, or ValueError or TypeError will be
        raised
        """
        self._order = float(value)

    def measure(self):
        """ Measures the current bond

        Returns
        -------
        measurement : float or None
            If the atoms have coordinates, returns the distance between the two
            atoms. If any coordinates are missing, returns None
        """
        if None in (self.atom1, self.atom2):
            return None
        try:
            return math.sqrt(distance2(self.atom1, self.atom2))
        except AttributeError:
            return None

    def umeasure(self):
        """ Same as "measure", but with units """
        m = self.measure()
        return m * u.angstroms if m is not None else None

    def energy(self):
        """ Measures the current bond energy

        Returns
        -------
        energy : float or None
            Bond strain energy in kcal/mol. Return value is None if either the
            coordinates of either atom is not set or the bond type is not set
        """
        d = self.measure()
        if self.type is None or d is None:
            return None
        dx = d - self.type.req
        return self.type.k * dx * dx

    def uenergy(self):
        """ Same as energy(), but with units """
        ene = self.energy()
        return ene * u.kilocalories_per_mole if ene is not None else None

    def __repr__(self):
        return f'<{self.__class__.__name__} {repr(self.atom1)}--{repr(self.atom2)}; type={repr(self.type)}>'

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class BondType(_ParameterType, _ListItem):
    """
    A bond type with a set of bond parameters

    Parameters
    ----------
    k : ``float``
        Force constant in kcal/mol/Angstrom^2
    req : ``float``
        Equilibrium bond distance in Angstroms
    list : :class:`TrackedList`
        A list of :class:`BondType`s in which this is a member

    Inherited Attributes
    --------------------
    idx : int
        The index of this BondType inside its containing list

    Notes
    -----
    Two `BondType`s are equal if their `k` and `req` attributes are equal

    Examples
    --------
    >>> bt1 = BondType(10.0, 1.0)
    >>> bt2 = BondType(10.0, 1.0)
    >>> bt1 is bt2
    False
    >>> bt1 == bt2
    True
    >>> bt1.idx # Not in a list or container
    -1

    As part of a list, they can be indexed

    >>> bond_list = []
    >>> bond_list.append(BondType(10.0, 1.0, list=bond_list))
    >>> bond_list.append(BondType(10.0, 1.0, list=bond_list))
    >>> bond_list[0].idx
    0
    >>> bond_list[1].idx
    1
    """

    def __init__(self, k, req, list=None):
        _ParameterType.__init__(self)
        self.k = _strip_units(k, u.kilocalories_per_mole/u.angstrom**2)
        self.req = _strip_units(req, u.angstrom)
        self.list = list
        self._idx = -1

    @_exception_to_notimplemented
    def __eq__(self, other):
        return abs(self.k - other.k) < TINY and abs(self.req - other.req) < TINY

    def __repr__(self):
        return '<%s; k=%.3f, req=%.3f>' % (type(self).__name__, self.k, self.req)

    def __copy__(self):
        """ Not bound to any list """
        # Hack to keep NoUreyBradley as a singleton
        if self is NoUreyBradley:
            return self
        return BondType(self.k, self.req)

    __getstate__ = _getstate_with_exclusions()

    def __reduce__(self):
        """
        Special-case NoUreyBradley, which should be a singleton. So if it's
        NoUreyBradley, return the same object. Otherwise, create a new one
        """
        if self is NoUreyBradley:
            return 'NoUreyBradley'
        return super().__reduce__()

    def __hash__(self):
        return hash((round(self.k, _TINY_DIGITS), round(self.req, _TINY_DIGITS)))

    @property
    def uk(self):
        return self.k * u.kilocalories_per_mole / u.angstroms**2

    @property
    def ureq(self):
        return self.req * u.angstroms

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Angle:
    """
    A valence angle between 3 atoms separated by two covalent bonds.

    Parameters
    ----------
    atom1 : :class:`Atom`
        An atom one end of the valence angle
    atom2 : :class:`Atom`
        The atom in the middle of the valence angle bonded to both other atoms
    atom3 : :class:`Atom`
        An atom on the other end of the valence angle to atom1
    type : :class:`AngleType`
        The AngleType object containing the parameters for this angle

    Notes
    -----
    An Angle can contain bonds or atoms. A bond is contained if it exists
    between atoms 1 and 2 or atoms 2 and 3.

    Examples
    --------
    >>> a1, a2, a3 = Atom(), Atom(), Atom()
    >>> angle = Angle(a1, a2, a3)
    >>> Bond(a1, a2) in angle and Bond(a3, a2) in angle
    True
    >>> a1 in angle and a2 in angle and a3 in angle
    True
    >>> Bond(a1, a3) in angle # this is not part of the angle definition
    False
    """

    def __init__(self, atom1, atom2, atom3, type=None):
        # Make sure we're not angling me to myself
        if atom1 is atom2 or atom1 is atom3 or atom2 is atom3:
            raise MoleculeError(f'Cannot angle atom to itself! Atoms are: {atom1} {atom2} {atom3}')
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        # Log these angles in each atom
        atom1.angles.append(self)
        atom2.angles.append(self)
        atom3.angles.append(self)
        # Load the force constant and equilibrium angle
        self.type = type
        atom1.angle_to(atom2)
        atom1.angle_to(atom3)
        atom2.angle_to(atom3)
        self.funct = 1

    def __contains__(self, obj):
        """ Quick and easy way to see if an Atom or a Bond is in this Angle """
        if isinstance(obj, Atom):
            return obj is self.atom1 or obj is self.atom2 or obj is self.atom3
        return ((self.atom1 in obj and self.atom2 in obj) or
                (self.atom2 in obj and self.atom3 in obj))

    def delete(self):
        """
        Deletes this angle from the atoms that make it up. This method removes
        this angle from the `angles` list for atom1, atom2, and atom3 as well as
        removing each atom form each others' angle partner lists
        """
        _delete_from_list(self.atom1.angles, self)
        _delete_from_list(self.atom2.angles, self)
        _delete_from_list(self.atom3.angles, self)

        _delete_from_list(self.atom1._angle_partners, self.atom2)
        _delete_from_list(self.atom1._angle_partners, self.atom3)
        _delete_from_list(self.atom2._angle_partners, self.atom1)
        _delete_from_list(self.atom2._angle_partners, self.atom3)
        _delete_from_list(self.atom3._angle_partners, self.atom1)
        _delete_from_list(self.atom3._angle_partners, self.atom2)

        self.atom1 = self.atom2 = self.atom3 = self.type = None

    def measure(self):
        """ Measures the current angle

        Returns
        -------
        measurement : float or None
            If the atoms have coordinates, returns the angle between the three
            atoms. If any coordinates are missing, returns None
        """
        if None in (self.atom1, self.atom2, self.atom3):
            return None
        try:
            return angle(self.atom1, self.atom2, self.atom3)
        except AttributeError:
            return None

    def umeasure(self):
        """ Same as measure(), but with units """
        m = self.measure()
        return m * u.degrees if m is not None else None

    def energy(self):
        """ Measures the current angle energy

        Returns
        -------
        energy : float or None
            Angle strain energy in kcal/mol. Return value is None if either the
            coordinates of either atom is not set or the angle type is not set
        """
        a = self.measure()
        if self.type is None or a is None:
            return None
        da = (a - self.type.theteq) * DEG_TO_RAD
        return self.type.k * da * da

    def uenergy(self):
        """ Same as energy(), but with units """
        ene = self.energy()
        return ene * u.kilocalories_per_mole if ene is not None else None

    def __repr__(self):
        return f'<{self.__class__.__name__} {repr(self.atom1)}--{repr(self.atom2)}--{repr(self.atom3)}; type={repr(self.type)}>'

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AngleType(_ParameterType, _ListItem):
    """
    An angle type with a set of angle parameters

    Parameters
    ----------
    k : ``float``
        Force constant in kcal/mol/radians^2
    theteq : ``float``
        Equilibrium angle in Degrees
    list : :class:`TrackedList`
        A list of `AngleType`s in which this is a member

    Inherited Attributes
    --------------------
    idx : ``int``
        The index of this AngleType inside its containing list

    Notes
    -----
    Two `AngleType`s are equal if their `k` and `theteq` attributes are equal

    Examples
    --------
    >>> at1 = AngleType(10.0, 180.0)
    >>> at2 = AngleType(10.0, 180.0)
    >>> at1 is at2
    False
    >>> at1 == at2
    True
    >>> at1.idx # not part of any list or iterable
    -1

    As part of a list, they can be indexed

    >>> angle_list = []
    >>> angle_list.append(AngleType(10.0, 180.0, list=angle_list))
    >>> angle_list.append(AngleType(10.0, 180.0, list=angle_list))
    >>> angle_list[0].idx
    0
    >>> angle_list[1].idx
    1
    """
    def __init__(self, k, theteq, list=None):
        _ParameterType.__init__(self)
        self.k = _strip_units(k, u.kilocalories_per_mole/u.radians**2)
        self.theteq = _strip_units(theteq, u.degrees)
        self._idx = -1
        self.list = list

    @_exception_to_notimplemented
    def __eq__(self, other):
        return abs(self.k - other.k) < TINY and abs(self.theteq - other.theteq) < TINY

    def __repr__(self):
        return f'<{self.__class__.__name__}; k={self.k:.3f}, theteq={self.theteq:.3f}>'

    def __copy__(self):
        return AngleType(self.k, self.theteq)

    __getstate__ = _getstate_with_exclusions()

    def __hash__(self):
        return hash((round(self.k, _TINY_DIGITS), round(self.theteq, _TINY_DIGITS)))

    @property
    def uk(self):
        return self.k * u.kilocalories_per_mole / u.radians**2

    @property
    def utheteq(self):
        return  self.theteq * u.degrees

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Dihedral(_FourAtomTerm):
    """ A valence dihedral between 4 atoms separated by three covalent bonds.

    Parameters
    ----------
    atom1 : :class:`Atom`
        An atom on one end of the valence dihedral bonded to atom 2
    atom2 : :class:`Atom`
        An atom in the middle of the valence dihedral bonded to atom1 and atom3
    atom3 : :class:`Atom`
        An atom in the middle of the valence dihedral bonded to atom2 and atom4
    atom4 : :class:`Atom`
        An atom on the other end of the valence dihedral bonded to atom 3
    improper : ``bool``
        If True, this is an Amber-style improper torsion, where atom3 is the
        "central" atom bonded to atoms 1, 2, and 4 (atoms 1, 2, and 4 are *only*
        bonded to atom 3 in this instance)
    ignore_end : ``bool``
        If True, the end-group interactions for this torsion are ignored, either
        because it is involved in a ring where the end-group atoms are excluded
        or because it is one term in a multi-term dihedral expansion
    type : :class:`DihedralType`
        The DihedralType object containing the parameters for this dihedral

    Notes
    -----
    A Dihedral can contain bonds or atoms. A bond is contained if it exists
    between atoms 1 and 2, between atoms 2 and 3, or between atoms 3 and 4.

    Examples
    --------
    >>> a1, a2, a3, a4 = Atom(), Atom(), Atom(), Atom()
    >>> dihed = Dihedral(a1, a2, a3, a4)
    >>> Bond(a1,a2) in dihed and Bond(a3,a2) in dihed and Bond(a3,a4) in dihed
    True
    >>> a1 in dihed and a2 in dihed and a3 in dihed and a4 in dihed
    True
    >>> Bond(a1, a4) in dihed # this is not part of the angle definition
    False

    For improper torsions, the bond pattern is different

    >>> a1, a2, a3, a4 = Atom(), Atom(), Atom(), Atom()
    >>> dihed = Dihedral(a1, a2, a3, a4, improper=True)
    >>> Bond(a1,a3) in dihed and Bond(a2,a3) in dihed and Bond(a3,a4) in dihed
    True
    >>> Bond(a1,a2) in dihed # Not like a normal dihedral!
    False
    >>> a1 in dihed and a2 in dihed and a3 in dihed and a4 in dihed
    True
    """

    def __init__(self, atom1, atom2, atom3, atom4, improper=False, ignore_end=False, type=None):
        _FourAtomTerm.__init__(self, atom1, atom2, atom3, atom4)
        # improper _implies_ ignore_end
        ignore_end = improper or ignore_end
        # Log these dihedrals in each atom
        atom1.dihedrals.append(self)
        atom2.dihedrals.append(self)
        atom3.dihedrals.append(self)
        atom4.dihedrals.append(self)
        self.type = type
        self.improper = improper
        self.ignore_end = ignore_end
        if improper:
            self._funct = 4
        else:
            atom1.dihedral_to(atom2)
            atom1.dihedral_to(atom3)
            atom1.dihedral_to(atom4)
            atom2.dihedral_to(atom3)
            atom2.dihedral_to(atom4)
            atom3.dihedral_to(atom4)
            self._funct = None

    @property
    def funct(self):
        if self._funct is None:
            # see https://manual.gromacs.org/documentation/current/reference-manual/topologies/topology-file-formats.html
            # for more info on funct types and their meaning
            if self.type is not None and isinstance(self.type, DihedralTypeList):
                return 9 # multiple proper dihedral
            elif self.type is not None and isinstance(self.type, DihedralType):
                return 1 # proper dihedral
            elif self.type is not None and isinstance(self.type, RBTorsionType):
                return 3 # Ryckaert-Bellemans dihedral
            elif self.type is not None:
                raise NotImplementedError(
                    f"Dihedral.type (self.type) must be: DihedralTypeList, DihedralType, or RBTorsionType."
                )
        return self._funct

    @funct.setter
    def funct(self, value):
        self._funct = value

    def __contains__(self, thing):
        """ Quick and easy way to find out if an Atom or Bond is in this Dihedral """
        if isinstance(thing, Atom):
            return _FourAtomTerm.__contains__(self, thing)
        # A dihedral is made up of 3 bonds
        if self.improper:
            # An improper is different... Atom 3 is the central atom
            return ((self.atom1 in thing and self.atom3 in thing) or
                    (self.atom2 in thing and self.atom3 in thing) or
                    (self.atom4 in thing and self.atom3 in thing))
        return ((self.atom1 in thing and self.atom2 in thing) or
                (self.atom2 in thing and self.atom3 in thing) or
                (self.atom3 in thing and self.atom4 in thing))

    def same_atoms(self, thing):
        """
        Determines if this dihedral has the same atoms (or atom indexes) as
        another object

        Parameters
        ----------
        thing : :class:`Dihedral` or other ``iterable``
            A Dihedral or an iterable with 4 indexes that will be used to
            identify whether this torsion has the same atoms

        Returns
        -------
        is_same : ``bool``
            True if ``thing`` and ``self`` have the same atoms. ``False``
            otherwise

        Notes
        -----
        This raises a ``TypeError`` if thing is not a Dihedral and is not
        iterable

        Raises
        ------
        TypeError
        """
        if isinstance(thing, Dihedral):
            # I'm comparing with another Dihedral here
            return ((self.atom1 is thing.atom1 and self.atom2 is thing.atom2 and
                     self.atom3 is thing.atom3 and self.atom4 is thing.atom4) or
                    (self.atom1 is thing.atom4 and self.atom2 is thing.atom3 and
                     self.atom4 is thing.atom1)
                   )
        thing = list(thing)
        # Here, atoms are expected to index from 0 (Python standard) if we
        # are comparing with a list or tuple
        if len(thing) != 4:
            raise TypeError(f'comparative {thing.__class__.__name__} has {len(thing)} elements! Expect 4.')
        # Compare starting_index, since we may not have an index right now
        return ( (self.atom1.idx == thing[0] and
                self.atom2.idx == thing[1] and
                self.atom3.idx == thing[2] and
                self.atom4.idx == thing[3]) or
                (self.atom1.idx == thing[3] and
                self.atom2.idx == thing[2] and
                self.atom3.idx == thing[1] and
                self.atom4.idx == thing[0]) )

    def delete(self):
        """
        Deletes this dihedral from the atoms that make it up. This method
        removes this dihedral from the `dihedrals` list for atom1, atom2, atom3,
        and atom4 as well as removing each atom form each others' dihedral
        partner lists
        """
        _delete_from_list(self.atom1.dihedrals, self)
        _delete_from_list(self.atom2.dihedrals, self)
        _delete_from_list(self.atom3.dihedrals, self)
        _delete_from_list(self.atom4.dihedrals, self)

        if not self.improper:
            _delete_from_list(self.atom1._dihedral_partners, self.atom2)
            _delete_from_list(self.atom1._dihedral_partners, self.atom3)
            _delete_from_list(self.atom1._dihedral_partners, self.atom4)
            _delete_from_list(self.atom2._dihedral_partners, self.atom1)
            _delete_from_list(self.atom2._dihedral_partners, self.atom3)
            _delete_from_list(self.atom2._dihedral_partners, self.atom4)
            _delete_from_list(self.atom3._dihedral_partners, self.atom1)
            _delete_from_list(self.atom3._dihedral_partners, self.atom2)
            _delete_from_list(self.atom3._dihedral_partners, self.atom4)
            _delete_from_list(self.atom4._dihedral_partners, self.atom1)
            _delete_from_list(self.atom4._dihedral_partners, self.atom2)
            _delete_from_list(self.atom4._dihedral_partners, self.atom3)

        self.atom1 = self.atom2 = self.atom3 = self.atom4 = self.type = None

    def measure(self):
        """ Measures the current dihedral angle

        Returns
        -------
        measurement : float or None
            If the atoms have coordinates, returns the dihedral angle between the four atoms in
            degrees. If any coordinates are missing, returns None
        """
        if None in (self.atom1, self.atom2, self.atom3, self.atom4):
            return None
        try:
            return dihedral(self.atom1, self.atom2, self.atom3, self.atom4)
        except AttributeError:
            return None

    def umeasure(self):
        """ Same as measure(), but with units """
        m = self.measure()
        return m * u.degrees if m is not None else None

    def energy(self):
        """ Measures the current dihedral angle energy

        Returns
        -------
        energy : float or None
            Dihedral angle energy in kcal/mol. Return value is None if either
            the coordinates of either atom is not set or the angle type is not
            set
        """
        phi = self.measure()
        if phi is None or self.type is None:
            return None
        phi *=  DEG_TO_RAD
        if isinstance(self.type, DihedralType):
            return self.type.phi_k * (1 + math.cos(self.type.per*phi - self.type.phase*DEG_TO_RAD))
        elif isinstance(self.type, DihedralTypeList):
            e = 0
            for term in self.type:
                e += term.phi_k * (1 + math.cos(term.per*phi - term.phase*DEG_TO_RAD))
            e *= 0.5
            return e
        else:
            raise NotImplementedError('Only DihedralType and DihedralTypeList '
                                      'are supported for energy calculations')

    def uenergy(self):
        """ Same as energy(), but with units """
        ene = self.energy()
        return ene * u.kilocalories_per_mole if ene is not None else None

    def __repr__(self):
        name = self.__class__.__name__
        if self.improper:
            name += ' [imp]'
        elif self.ignore_end:
            name += ' [ign]'
        return f'<{name}; {repr(self.atom1)}--{repr(self.atom2)}--{repr(self.atom3)}--{repr(self.atom4)}; type={repr(self.type)}>'

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class DihedralType(_ParameterType, _ListItem):
    """
    A dihedral type with a set of dihedral parameters

    Parameters
    ----------
    phi_k : ``float``
        The force constant in kcal/mol
    per : ``int``
        The dihedral periodicity
    phase : ``float``
        The dihedral phase in degrees
    scee : ``float``, optional
        1-4 electrostatic scaling factor. Default is 1.0
    scnb : ``float``, optional
        1-4 Lennard-Jones scaling factor. Default is 1.0
    list : :class:`TrackedList`
        A list of `DihedralType`s in which this is a member

    Inherited Attributes
    --------------------
    idx : ``int``
        The index of this DihedralType inside its containing list

    Notes
    -----
    Two :class:`DihedralType`s are equal if their `phi_k`, `per`, `phase`,
    `scee`, and `scnb` attributes are equal

    Examples
    --------
    >>> dt1 = DihedralType(10.0, 2, 180.0, 1.2, 2.0)
    >>> dt2 = DihedralType(10.0, 2, 180.0, 1.2, 2.0)
    >>> dt1 is dt2
    False
    >>> dt1 == dt2
    True
    >>> dt1.idx # not part of any list or iterable
    -1

    As part of a list, they can be indexed

    >>> dihedral_list = []
    >>> dihedral_list.append(DihedralType(10.0, 2, 180.0, 1.2, 2.0,
    ...                                   list=dihedral_list))
    >>> dihedral_list.append(DihedralType(10.0, 2, 180.0, 1.2, 2.0,
    ...                                   list=dihedral_list))
    >>> dihedral_list[0].idx
    0
    >>> dihedral_list[1].idx
    1
    """

    #===================================================

    def __init__(self, phi_k, per, phase, scee=1.0, scnb=1.0, list=None):
        """ DihedralType constructor """
        _ParameterType.__init__(self)
        self.locked = True
        self.phi_k = _strip_units(phi_k, u.kilocalories_per_mole)
        self.per = int(per)
        self.phase = _strip_units(phase, u.degrees)
        self.scee = scee
        self.scnb = scnb
        self.list = list
        self._idx = -1
        self.locked = True

    #===================================================

    @_exception_to_notimplemented
    def __eq__(self, other):
        return (abs(self.phi_k - other.phi_k) < TINY and self.per == other.per and
                abs(self.phase - other.phase) < TINY and abs(self.scee - other.scee) < TINY and
                abs(self.scnb - other.scnb) < TINY)

    def __repr__(self):
        retstr = [f'<{self.__class__.__name__}; phi_k={self.phi_k:.3f}, per={self.per}, phase={self.phase:.3f}, ']
        retstr.append(f' scee={self.scee:.3f}, scnb={self.scnb:.3f}>')
        return ''.join(retstr)

    def __copy__(self):
        return DihedralType(self.phi_k, self.per, self.phase, self.scee, self.scnb)

    __getstate__ = _getstate_with_exclusions()

    def __hash__(self):
        return hash((round(self.phi_k, _TINY_DIGITS), round(self.per, _TINY_DIGITS),
                     round(self.phase, _TINY_DIGITS), round(self.scee, _TINY_DIGITS),
                     round(self.scnb, _TINY_DIGITS)))

    @property
    def uphi_k(self):
        return self.phi_k * u.kilocalories_per_mole

    @property
    def uphase(self):
        return self.phase * u.degrees

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class RBTorsionType(_ParameterType, _ListItem):
    """
    A Ryckaert-Bellemans type with a set of dihedral parameters

    Parameters (and Attributes)
    ---------------------------
    c0 : float
        The coefficient of the constant term in kcal/mol
    c1 : float
        The coefficient of the linear term in kcal/mol
    c2 : float
        The coefficient of the quadratic term in kcal/mol
    c3 : float
        The coefficient of the cubic term in kcal/mol
    c4 : float
        The coefficient of the quartic term in kcal/mol
    c5 : float
        The coefficient of the quintic term in kcal/mol
    scee : float, optional
        1-4 electrostatic scaling factor. Default is 1.0
    scnb : float, optional
        1-4 van der Waals scaling factor. Default is 1.0
    list : TrackedList=None
        A list of `RBTorsionType`s in which this is a member

    Inherited Attributes
    --------------------
    idx : int
        The index of this RBTorsionType inside its containing list

    Notes
    -----
    Two `RBTorsionType`s are equal if their coefficients are all equal

    Examples
    --------
    >>> dt1 = RBTorsionType(10.0, 20.0, 30.0, 40.0, 50.0, 60.0)
    >>> dt2 = RBTorsionType(10.0, 20.0, 30.0, 40.0, 50.0, 60.0)
    >>> dt1 is dt2
    False
    >>> dt1 == dt2
    True
    >>> dt1.idx # not part of any list or iterable
    -1

    As part of a list, they can be indexed

    >>> rb_torsion_list = []
    >>> rb_torsion_list.append(RBTorsionType(10.0, 20.0, 30.0, 40.0, 50.0, 60.0,
    ...                                      list=rb_torsion_list))
    >>> rb_torsion_list.append(RBTorsionType(10.0, 20.0, 30.0, 40.0, 50.0, 60.0,
    ...                                      list=rb_torsion_list))
    >>> rb_torsion_list[0].idx
    0
    >>> rb_torsion_list[1].idx
    1
    """

    #===================================================

    def __init__(self, c0, c1, c2, c3, c4, c5, scee=1.0, scnb=1.0, list=None):
        _ParameterType.__init__(self)
        self.c0 = _strip_units(c0, u.kilocalories_per_mole)
        self.c1 = _strip_units(c1, u.kilocalories_per_mole)
        self.c2 = _strip_units(c2, u.kilocalories_per_mole)
        self.c3 = _strip_units(c3, u.kilocalories_per_mole)
        self.c4 = _strip_units(c4, u.kilocalories_per_mole)
        self.c5 = _strip_units(c5, u.kilocalories_per_mole)
        self.scee = scee
        self.scnb = scnb
        self.list = list
        self._idx = -1

    #===================================================

    @_exception_to_notimplemented
    def __eq__(self, other):
        return (abs(self.c0 - other.c0) < TINY and
                abs(self.c1 - other.c1) < TINY and
                abs(self.c2 - other.c2) < TINY and
                abs(self.c3 - other.c3) < TINY and
                abs(self.c4 - other.c4) < TINY and
                abs(self.c5 - other.c5) < TINY)

    def __copy__(self):
        return RBTorsionType(self.c0, self.c1, self.c2, self.c3, self.c4,
                             self.c5, self.scee, self.scnb)

    def __repr__(self):
        return (
            f'<{self.__class__.__name__}; c0={self.c0:.3f}; c1={self.c1:.3f}; c2={self.c2:.3f}; '
            f'c3={self.c3:.3f}; c4={self.c4:.3f}; c5={self.c5:.3f}; scee={self.scee}; scnb={self.scnb}>'
        )

    __getstate__ = _getstate_with_exclusions()

    def __hash__(self):
        return hash((round(self.c0, _TINY_DIGITS), round(self.c1, _TINY_DIGITS),
                     round(self.c2, _TINY_DIGITS), round(self.c3, _TINY_DIGITS),
                     round(self.c4, _TINY_DIGITS), round(self.c5, _TINY_DIGITS)))

    @property
    def uc0(self):
        return self.c0 * u.kilocalories_per_mole

    @property
    def uc1(self):
        return self.c2 * u.kilocalories_per_mole

    @property
    def uc2(self):
        return self.c2 * u.kilocalories_per_mole

    @property
    def uc3(self):
        return self.c3 * u.kilocalories_per_mole

    @property
    def uc4(self):
        return self.c4 * u.kilocalories_per_mole

    @property
    def uc5(self):
        return self.c5 * u.kilocalories_per_mole

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class DihedralTypeList(list, _ListItem):
    """
    Dihedral types are a Fourier expansion of terms. In some cases, they are
    stored in a list like this one inside another TrackedList. In other cases,
    each term is a separate entry in the :class:`TrackedList`.

    In cases where `DihedralType`s are stored with every term in the same
    container, this object supports list assignment and indexing like
    :class:`DihedralType`.

    Parameters
    ----------
    *args : objects
        Any arguments that ``list`` would take.
    list : TrackedList, optional
        A list that "contains" this DihedralTypeList instance. This is a
        keyword-only argument. Default is ``None`` (i.e., belonging to no list)
    **kwargs : keyword argument list
        All other keyword arguments passed directly to the ``list`` constructor
    """
    def __init__(self, *args, **kwargs):
        if 'list' in kwargs:
            self.list = kwargs.pop('list')
        else:
            self.list = None
        list.__init__(self, *args, **kwargs)
        self._idx = -1
        self.used = False

    @classmethod
    def from_rbtorsion(cls, rbtorsion):
        """
        Creates a Fourier series of proper torsions from a Ryckaerts-Bellemans
        torsion.

        Parameters
        ----------
        rbtorsion : RBTorsionType
            The R-B torsion type to convert to a series of proper torsions

        Raises
        ------
        ValueError if all terms in rbtorsion are zero except for c0

        """
        c0 = rbtorsion.c0
        c1 = rbtorsion.c1
        c2 = rbtorsion.c2
        c3 = rbtorsion.c3
        c4 = rbtorsion.c4
        c5 = rbtorsion.c5

        phi = 0 * u.degrees
        fc0 = (4 * c0 + 4 * c1 + 4 * c3 - c4 + 4 * c5) / 8.
        fc1 = (-8 * c1 - 6 * c3 - 5 * c5) / 8
        fc2 = (c2 + c4) / 2
        fc3 = (-4 * c3 - 5 * c5) / 16
        fc4 = c4 / 8
        fc5 = -c5 / 16

        inst = cls()
        for i, f in enumerate((fc0, fc1, fc2, fc3, fc4, fc5)):
            if abs(f) > TINY:
                inst.append(DihedralType(f, i, phi, scee=rbtorsion.scee, scnb=rbtorsion.scnb))
        if len(inst) == 0:
            # All force constants were zeros:
            if abs(fc0) < TINY:
                inst.append(DihedralType(0.0, 0, phi, scee=rbtorsion.scee, scnb=rbtorsion.scnb))
            else:
                raise ValueError('Unable to convert RB torsion to propers.')
        return inst

    @_exception_to_notimplemented
    def __eq__(self, other):
        if len(self) != len(other): return False
        for t1, t2 in zip(self, other):
            if not t1 == t2:
                return False
        return True

    def append(self, other, override=False):
        """ Adds a DihedralType to the DihedralTypeList

        Parameters
        ----------
        other : :class:`DihedralType`
            The DihedralType instance to add to this list. It cannot have the
            same periodicity as any other DihedralType in the list
        override : bool, optional, default=False
            If True, this will override an existing torsion, but will raise a
            warning if it does so

        Raises
        ------
        TypeError if other is not an instance of :class:`DihedralTypeList`

        ParameterError if other has the same periodicity as another member in
        this list and override is False

        ParameterWarning if other has same periodicit as another member in this
        list and override is True
        """
        if not isinstance(other, DihedralType):
            raise TypeError('Can only add DihedralType to DihedralTypeList')
        for i, existing in enumerate(self):
            if existing.per == other.per:
                # Do not add duplicate periodicities
                if existing == other: return
                if override:
                    warnings.warn('Overriding DihedralType in DihedralTypeList '
                                  'with same periodicity', ParameterWarning)
                    self[i] = other
                else:
                    raise ParameterError('Cannot add two DihedralType instances with the same '
                                         'periodicity to the same DihedralTypeList')
        list.append(self, other)

    @property
    def penalty(self):
        penalty = None
        for dt in self:
            if dt.penalty is not None:
                if penalty is None:
                    penalty = dt.penalty
                else:
                    penalty = max(dt.penalty, penalty)
        return penalty

    def __repr__(self):
        return '<DihedralTypes %s>' % (list.__repr__(self))

    def __copy__(self):
        return DihedralTypeList([copy(x) for x in self])

    __getstate__ = _getstate_with_exclusions(['list', 'penalty'])

    def __hash__(self):
        return hash(tuple(self))

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class UreyBradley(Bond):
    """
    A Urey-Bradley angle type with a set of parameters. It has the same
    functional form as a bond, but it is defined between two atoms forming a
    valence angle separated by two bonds.

    Parameters
    ----------
    atom1 : :class:`Atom`
        The first atom involved in the Urey-Bradley bond
    atom2 : :class:`Atom`
        The other atom involved in the Urey-Bradley bond
    type : :class:`BondType`
        The Urey-Bradley bond type that defines the parameters for this bond

    Notes
    -----
    You can test whether an :class:`Atom` is contained within the bond using the
    `in` operator. A :class:`MoleculeError` is raised if `atom1` and `atom2` are
    identical.  You can also test that a :class:`Bond` is contained in this
    Urey-Bradley valence angle

    Examples
    --------
    >>> a1, a2, a3 = Atom(), Atom(), Atom()
    >>> b1 = Bond(a1, a2)
    >>> b2 = Bond(a2, a3)
    >>> angle = Angle(a1, a2, a3)
    >>> urey = UreyBradley(a1, a3)
    >>> a1 in urey and a3 in urey
    True
    >>> b1 in urey and b2 in urey
    True
    """
    def __init__(self, atom1, atom2, type=None):
        """ Bond constructor """
        # Make sure we're not bonding me to myself
        if atom1 is atom2:
            raise MoleculeError(f'Cannot bond atom to itself! Atoms are {atom1} {atom2}')
        # Order the atoms so the lowest atom # is first
        self.atom1 = atom1
        self.atom2 = atom2
        # Log this urey-bradley in the atoms
        atom1.urey_bradleys.append(self)
        atom2.urey_bradleys.append(self)
        # Load the force constant and equilibrium distance
        self.type = type

    def __contains__(self, thing):
        " Quick and easy way to see if an Atom or Bond is in this Urey-Bradley "
        if isinstance(thing, Atom):
            return thing is self.atom1 or thing is self.atom2
        # If this is a bond things are a bit more complicated since we don't
        # know the central atom of this Urey-Bradley. We need to make sure that
        # one of the atoms of the bond is either atom1 or atom2 and that the
        # OTHER atom in Bond "thing" has the OTHER atom in this Urey-Bradley in
        # its list of bonded partners.
        if not thing.atom1 in self:
            if not thing.atom2 in self:
                # Neither atom is in this Urey-Bradley...
                return False
            # If we are here, thing.atom2 is in self
            end1 = thing.atom2
            cent = thing.atom1
        else:
            # If we are here, thing.atom1 is in self
            end1 = thing.atom1
            cent = thing.atom2
        # If we are here, end1 and cent are set. Look through the bond
        # partners of cent(er) and see if any of them is in this
        # Urey-Bradley (but ONLY if that atom is not the original end1)
        for atm in cent.bond_partners:
            if atm is end1: continue
            if atm in self:
            # If we got here, we found both atoms in this Urey-Bradley
            # separated by 2 bonds, so "thing" IS in this Urey-Bradley
                return True
        # If we got here, we could not find the other atom connected by 2
        # bonds
        return False

    def delete(self):
        """
        Deletes this Urey-Bradley from the atoms that make it up. This method
        removes this urey-bradley from the `urey_bradleys` list for atom1 and
        atom2.
        """
        _delete_from_list(self.atom1.urey_bradleys, self)
        _delete_from_list(self.atom2.urey_bradleys, self)

        self.atom1 = self.atom2 = self.type = None

    def __repr__(self):
        return f'<{self.__class__.__name__} {repr(self.atom1)}--{repr(self.atom2)}; type={repr(self.type)}>'

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Improper(_FourAtomTerm):
    """
    A CHARMM-style improper torsion between 4 atoms. The first atom must be the
    central atom, as shown in the schematic below

                      A3
                      |
                      |
             A4 ----- A1 ----- A2

    Parameters
    ----------
    atom1 : :class:`Atom`
        The central atom A1 in the schematic above
    atom2 : :class:`Atom`
        An atom in the improper, A2 in the schematic above
    atom3 : :class:`Atom`
        An atom in the improper, A3 in the schematic above
    atom4 : :class:`Atom`
        An atom in the improper, A4 in the schematic above
    type : :class:`ImproperType`
        The ImproperType object containing the parameters for this improper
        torsion

    Notes
    -----
    An Improper torsion can contain bonds or atoms. A bond is contained if it
    exists between atom 1 and any other atom. Raises :class:`MoleculeError` if any
    of the atoms are duplicates.

    Examples
    --------
    >>> a1, a2, a3, a4 = Atom(), Atom(), Atom(), Atom()
    >>> imp = Improper(a1, a2, a3, a4)
    >>> Bond(a1, a2) in imp and Bond(a1, a3) in imp and Bond(a1, a4) in imp
    True
    >>> Bond(a2, a3) in imp
    False
    """
    def __init__(self, atom1, atom2, atom3, atom4, type=None):
        _FourAtomTerm.__init__(self, atom1, atom2, atom3, atom4)
        # Log these impropers in each atom
        atom1.impropers.append(self)
        atom2.impropers.append(self)
        atom3.impropers.append(self)
        atom4.impropers.append(self)
        # Load the force constant and equilibrium angle
        self.type = type
        self.funct = 2

    def __contains__(self, thing):
        """
        Quick and easy way to find out if an Atom or Bond is in this Improper
        """
        if isinstance(thing, Atom):
            return _FourAtomTerm.__contains__(self, thing)
        # Treat it like a Bond
        return ((self.atom1 in thing and self.atom2 in thing) or
                (self.atom1 in thing and self.atom3 in thing) or
                (self.atom1 in thing and self.atom4 in thing))

    def same_atoms(self, thing):
        """
        An improper has the same 4 atoms if atom1 is the same for both and the
        other 3 can be in any order
        """
        if isinstance(thing, Improper):
            # I'm comparing with another Improper here. Central atom must be
            # the same. Others can be in any order
            if self.atom1 is not thing.atom1:
                return False
            # Make a set with the remaining atoms. If they are equal, the
            # impropers are equivalent
            selfset = set([self.atom2, self.atom3, self.atom4])
            otherset = set([thing.atom2, thing.atom3, thing.atom4])
            return selfset == otherset
        thing = list(thing)
        # Here, atoms are expected to index from 0 (Python standard) if we
        # are comparing with a list or tuple
        if len(thing) != 4:
            raise MoleculeError(f'Impropers have 4 atoms, not {len(thing)}')
        if self.atom1.idx != thing[0]:
            return False
        selfset = set([self.atom2.idx, self.atom3.idx, self.atom4.idx])
        otherset = set([thing[1], thing[2], thing[3]])
        return selfset == otherset

    def measure(self):
        """ Measures the current torsional angle

        Returns
        -------
        measurement : float or None
            If the atoms have coordinates, returns the torsional angle between
            the four atoms. If any coordinates are missing, returns None
        """
        if None in (self.atom1, self.atom2, self.atom3, self.atom4):
            return None
        try:
            return dihedral(self.atom1, self.atom2, self.atom3, self.atom4)
        except AttributeError:
            return None

    def umeasure(self):
        """ Same as measure(), but with units """
        m = self.measure()
        return m * u.degrees if m is not None else None

    def energy(self):
        """ Measures the current dihedral angle energy

        Returns
        -------
        energy : float or None
            Dihedral angle energy in kcal/mol. Return value is None if either
            the coordinates of either atom is not set or the angle type is not
            set
        """
        phi = self.measure()
        if phi is None or self.type is None:
            return None
        dphi = (self.type.psi_eq - phi) * DEG_TO_RAD
        return self.type.psi_k * dphi * dphi

    def uenergy(self):
        """ Same as energy(), but with units """
        ene = self.energy()
        return ene * u.kilocalories_per_mole if ene is not None else None

    def delete(self):
        """
        Deletes this Improper from the atoms that make it up. This method
        removes this Improper from the `impropers` list for atom1, atom2, atom3,
        and atom4
        """
        _delete_from_list(self.atom1.impropers, self)
        _delete_from_list(self.atom2.impropers, self)
        _delete_from_list(self.atom3.impropers, self)
        _delete_from_list(self.atom4.impropers, self)

        self.atom1 = self.atom2 = self.atom3 = self.atom4 = self.type = None

    def __repr__(self):
        return (
            f'<{self.__class__.__name__}; {repr(self.atom1)}--({repr(self.atom2)},'
            f'{repr(self.atom3)},{repr(self.atom4)}); type={repr(self.type)}>'
        )

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ImproperType(_ParameterType, _ListItem):
    """
    An improper type with a set of improper torsion parameters

    Parameters
    ----------
    psi_k : ``float``
        Force constant in kcal/mol/radians^2
    psi_eq : ``float``
        Equilibrium torsion angle in Degrees
    list : :class:`TrackedList`
        A list of `ImproperType`s in which this is a member

    Inherited Attributes
    --------------------
    idx : ``int``
        The index of this ImproperType inside its containing list

    Notes
    -----
    Two `ImproperType`s are equal if their `psi_k` and `psi_eq` attributes are
    equal

    Examples
    --------
    >>> it1 = ImproperType(10.0, 180.0)
    >>> it2 = ImproperType(10.0, 180.0)
    >>> it1 is it2
    False
    >>> it1 == it2
    True
    >>> it1.idx # Not part of any list or iterable
    -1

    As part of a list, they can be indexed

    >>> improper_list = []
    >>> improper_list.append(ImproperType(10.0, 180.0, list=improper_list))
    >>> improper_list.append(ImproperType(10.0, 180.0, list=improper_list))
    >>> improper_list[0].idx
    0
    >>> improper_list[1].idx
    1
    """
    def __init__(self, psi_k, psi_eq, list=None):
        _ParameterType.__init__(self)
        self.psi_k = _strip_units(psi_k, u.kilocalories_per_mole/u.radians**2)
        self.psi_eq = _strip_units(psi_eq, u.degrees)
        self.list = list
        self._idx = -1

    @_exception_to_notimplemented
    def __eq__(self, other):
        return (abs(self.psi_k - other.psi_k) < TINY and
                abs(self.psi_eq - other.psi_eq) < TINY)

    def __repr__(self):
        return f'<{self.__class__.__name__}; psi_k={self.psi_k:.3f}, psi_eq={self.psi_eq:.3f}>'

    def __copy__(self):
        return ImproperType(self.psi_k, self.psi_eq)

    __getstate__ = _getstate_with_exclusions()

    def __hash__(self):
        return hash((round(self.psi_k, _TINY_DIGITS), round(self.psi_eq, _TINY_DIGITS)))

    @property
    def upsi_k(self):
        return self.psi_k * u.kilocalories_per_mole / u.radians ** 2

    @property
    def upsi_eq(self):
        return self.psi_eq * u.degrees

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Cmap:
    """
    A coupled-torsion correction map term defined between 5 atoms connected by
    four covalent bonds. This is a coupled-torsion potential in which the
    torsions are consecutive.

    Parameters
    ----------
    atom1 : :class:`Atom`
        An atom on one end of the valence coupled-torsion bonded to atom2
    atom2 : :class:`Atom`
        An atom in the middle of the CMAP bonded to atoms 1 and 3
    atom3 : :class:`Atom`
        An atom in the middle of the CMAP bonded to atoms 2 and 4
    atom4 : :class:`Atom`
        An atom in the middle of the CMAP bonded to atoms 3 and 5
    atom5 : :class:`Atom`
        An atom in the middle of the CMAP bonded to atom 4
    type : :class:`CmapType`
        The CmapType object containing the parameter map for this term

    Notes
    -----
    A CMAP can contain bonds or atoms. A bond is contained if it exists between
    atoms 1 and 2, between atoms 2 and 3, between atoms 3 and 4, or between
    atoms 4 and 5.

    Examples
    --------
    >>> a1, a2, a3, a4, a5 = Atom(), Atom(), Atom(), Atom(), Atom()
    >>> cmap = Cmap(a1, a2, a3, a4, a5)
    >>> Bond(a1, a2) in cmap and Bond(a2, a3) in cmap
    True
    >>> Bond(a1, a3) in cmap
    False
    """
    def __init__(self, atom1, atom2, atom3, atom4, atom5, type=None):
        # Make sure we're not CMAPping me to myself
        atmlist = [atom1, atom2, atom3, atom4, atom5]
        for i in range(len(atmlist)):
            for j in range(i+1, len(atmlist)):
                if atmlist[i] is atmlist[j]:
                    raise MoleculeError(f'Cannot cmap atom to itself! Atoms are {atmlist[i]} {atmlist[j]}')
        # Set up instances
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.atom5 = atom5
        # Log these cmaps in each atom
        atom1.cmaps.append(self)
        atom2.cmaps.append(self)
        atom3.cmaps.append(self)
        atom4.cmaps.append(self)
        atom5.cmaps.append(self)
        # Load the CMAP interpolation table
        self.type = type
        self.funct = 1

    @classmethod
    def extended(cls, atom1, atom2, atom3, atom4, atom5, atom6, atom7, atom8, type=None):
        """
        Alternative constructor for correction maps defined with 8 atoms (each
        torsion being separately specified). Correction maps are, to the best of
        my knowledge, used to parametrized "consecutive" torsions (i.e., atoms
        2, 3, and 4 are common to both torsions). However, the CHARMM definition
        specifies the torsions separately.

        If the torsions are _not_ consecutive (i.e., if atoms 2, 3, and 4 are
        not the same as atoms 5, 6, and 7), NotImplementedError is raised.

        Parameters
        ----------
        atom1 : :class:`Atom`
            The first atom of the first torsion
        atom2 : :class:`Atom`
            The second atom of the first torsion
        atom3 : :class:`Atom`
            The third atom of the first torsion
        atom4 : :class:`Atom`
            The fourth atom of the first torsion
        atom5 : :class:`Atom`
            The first atom of the second torsion
        atom6 : :class:`Atom`
            The second atom of the second torsion
        atom7 : :class:`Atom`
            The third atom of the second torsion
        atom8 : :class:`Atom`
            The fourth atom of the second torsion
        type : :class:`CmapType`
            The CmapType object containing the parameter map for this term
        """
        if atom2 is not atom5 or atom3 is not atom6 or atom7 is not atom4:
            raise NotImplementedError('Only consecutive coupled-torsions are '
                                      'supported by CMAPs currently')
        return cls(atom1, atom2, atom3, atom4, atom8, type=type)

    def __contains__(self, thing):
        """
        Quick and easy way to find out an atom or bond is in this
        coupled-torsion
        """
        if isinstance(thing, Atom):
            return (thing is self.atom1 or thing is self.atom2 or
                    thing is self.atom3 or thing is self.atom4 or
                    thing is self.atom5)
        return ((self.atom1 in thing and self.atom2 in thing) or
                (self.atom2 in thing and self.atom3 in thing) or
                (self.atom3 in thing and self.atom4 in thing) or
                (self.atom4 in thing and self.atom5 in thing))

    def same_atoms(self, thing):
        """
        A coupled-torsion is equivalent if the 5 atoms are in the same or
        reverse order Allow comparison with another type of cmap or with a
        sequence of 5 indexes
        """
        if isinstance(thing, Cmap):
            # I'm comparing with another Improper here. Central atom must be the
            # same. Others can be in any order
            return ((self.atom1 is thing.atom1 and self.atom2 is thing.atom2 and
                    self.atom3 is thing.atom3 and self.atom4 is thing.atom4 and
                    self.atom5 is thing.atom5) or
                    (self.atom1 is thing.atom5 and self.atom2 is thing.atom4 and
                    self.atom3 is thing.atom3 and self.atom4 is thing.atom2 and
                    self.atom5 is thing.atom1))
        thing = list(thing)
        # Here, atoms are expected to index from 0 (Python standard) if we
        # are comparing with a list or tuple
        if len(thing) != 5:
            raise MoleculeError(f'CMAP can compare to 5 elements, not {len(thing)}')
        return ((self.atom1.idx == thing[0] and self.atom2.idx == thing[1] and
                 self.atom3.idx == thing[2] and self.atom4.idx == thing[3] and
                 self.atom5.idx == thing[4]) or
                (self.atom1.idx == thing[4] and self.atom2.idx == thing[3] and
                 self.atom3.idx == thing[2] and self.atom4.idx == thing[1] and
                 self.atom5.idx == thing[0]))

    def delete(self):
        """
        Deletes this Cmap from the atoms that make it up. This method removes
        the Cmap from the `cmaps` list for atom1, atom2, atom3, atom4, and atom5
        """
        _delete_from_list(self.atom1.cmaps, self)
        _delete_from_list(self.atom2.cmaps, self)
        _delete_from_list(self.atom3.cmaps, self)
        _delete_from_list(self.atom4.cmaps, self)
        _delete_from_list(self.atom5.cmaps, self)

        self.atom1 = self.atom2 = self.atom3 = self.atom4 = self.atom5 = self.type = None

    def __repr__(self):
        return (
            f'<{self.__class__.__name__}; {repr(self.atom1)}--{repr(self.atom2)}--'
            f'{repr(self.atom3)}--{repr(self.atom4)}--{repr(self.atom5)}; type={repr(self.type)}>'
        )

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class CmapType(_ParameterType, _ListItem):
    """
    A CMAP type with a potential energy interpoloation grid mapping out the 2-D
    potential of coupled torsions.

    Parameters
    ----------
    resolution : ``int``
        The number of grid points in the correction map potential in both
        torsion dimensions
    grid : ``iterable of floats``
        This must be a 1-dimensional list of grid values. The dimension must be
        `resolution*resolution`, and must be row-major (i.e., the second
        dimension changes the fastest) when indexed with 2 indices.
    list : :class:`TrackedList`
        A list of `CmapType`s in which this is a member

    Attributes
    ----------
    resolution : ``int``
        Potential grid resolution (see description in Parameters)
    grid : :class:`_CmapGrid`
        A _CmapGrid object defining the interpolating potential energy grid,
        with each point having the units kcal/mol
    comments : ``list(str)``
        List of strings that represent comments about this parameter type
    list : :class:`TrackedList`
        If not None, this is a list in which this instance _may_ be a member
    idx : ``int``
        The index of this CmapType inside its containing list

    Notes
    -----
    Two `CmapType`s are equal if their resolution is the same and each grid
    point is the same to within 1E-8

    See the docs for `_CmapGrid` for information on how to access the
    interpolating grid data if necessary.

    Raises
    ------
    TypeError is raised if the grid does not have the correct number of elements
    for the given resolution

    Examples
    --------
    >>> ct = CmapType(2, [0, 1, 2, 3])
    >>> ct.grid[0,0], ct.grid[0]
    (0, 0)
    >>> ct.grid[0,1], ct.grid[1]
    (1, 1)
    >>> ct.grid[1,0], ct.grid[2]
    (2, 2)
    >>> ct.grid[1,1], ct.grid[3]
    (3, 3)
    >>> ct.idx # not part of a list or iterable
    -1
    """
    def __init__(self, resolution, grid, comments=None, list=None):
        _ParameterType.__init__(self)
        self.resolution = resolution
        self.grid = _CmapGrid(resolution, grid)
        if len(grid) != self.resolution * self.resolution:
            raise TypeError('CMAP grid does not match expected resolution')
        if comments is None:
            self.comments = []
        else:
            self.comments = comments
        self._idx = -1
        self.list = list

    @_exception_to_notimplemented
    def __eq__(self, other):
        return (self.resolution == other.resolution and
                all(abs(i - j) < TINY for i, j in zip(self.grid, other.grid)))

    def __repr__(self):
        return f'<{self.__class__.__name__}; resolution={self.resolution}>'

    def __copy__(self):
        return CmapType(self.resolution, copy(self.grid._data), self.comments[:])

    __getstate__ = _getstate_with_exclusions()

    def __hash__(self):
        return hash((self.resolution, self.grid))

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class _CmapGrid:
    """
    A grid object for storing Correction map data. Data can be accessed in one
    of two ways; either with 1 or 2 indexes. If 2 indexes [i,j] are given, the
    index into the flattened array is i*resolution+j. Indexing starts from 0.

    The _CmapGrid usually has ranges for the two angles from -180 to 180. Some
    places will expect the range to be 0-360 degrees (e.g., OpenMM). The
    switch_range method returns a _CmapGrid with the "other" range than the
    current object.  i.e., if the current range is from 0 to 360 degrees, the
    _CmapGrid returned by `switch_range` will be from -180 to 180, and vice
    versa.

    Example:
    >>> g = _CmapGrid(2, [0, 1, 2, 3])
    >>> print('%s %s' % (g[0], g[0,0]))
    0 0
    >>> print('%s %s' % (g[1], g[0,1]))
    1 1
    >>> print(g[1,0])
    2
    >>> g[1,1] = 10
    >>> print(g[3])
    10
    >>> print(g.switch_range())
    [10.0000, 2.0000
     1.0000, 0.0000]
    """

    def __init__(self, resolution, data=None):
        self.resolution = resolution
        if data is None:
            self._data = [0 for i in range(self.resolution*self.resolution)]
        else:
            self._data = _strip_units(data, u.kilocalories_per_mole)

    @property
    def transpose(self):
        """ The transpose of the potential grid """
        if hasattr(self, '_transpose'):
            return self._transpose
        _transpose = []
        size = len(self._data)
        for i in range(self.resolution):
            piece = [self[j] for j in range(i, size, self.resolution)]
            _transpose += piece
        self._transpose = _CmapGrid(self.resolution, _transpose)
        return self._transpose

    T = transpose

    def __getitem__(self, idx):
        if isinstance(idx, tuple):
            if idx[0] >= self.resolution or idx[1] >= self.resolution:
                raise IndexError('_CmapGrid: Index out of range')
            return self._data[self.resolution*idx[0]+idx[1]]
        return self._data[idx]

    def __setitem__(self, idx, val):
        if isinstance(idx, tuple):
            if idx[0] >= self.resolution or idx[1] >= self.resolution:
                raise IndexError('_CmapGrid: Index out of range')
            self._data[self.resolution*idx[0]+idx[1]] = _strip_units(val, u.kilocalories_per_mole)
        else:
            try:
                indices = range(*idx.indices(len(self._data)))
            except AttributeError:
                self._data[idx] = _strip_units(val, u.kilocalories_per_mole)
            else:
                try:
                    lenval = len(val)
                except TypeError:
                    lenval = 1
                if lenval == 1:
                    for x in indices:
                        self._data[x] = _strip_units(val, u.kilocalories_per_mole)
                elif lenval != len(indices):
                    raise ValueError('Wrong number of values setting a slice')
                else:
                    for x, y in zip(indices, val):
                        self._data[x] = _strip_units(y, u.kilocalories_per_mole)

    def __len__(self):
        return len(self._data)

    def __iter__(self):
        return iter(self._data)

    def __repr__(self):
        return f'<_CmapGrid: {self.resolution}x{self.resolution}>'

    def __str__(self):
        retstr = f'[{self._data[0]:.4f},'
        fmt = ' {0:.4f}'
        for i, val in enumerate(self):
            if i == 0: continue
            retstr += fmt.format(val)
            if (i+1) % self.resolution == 0 and i != len(self._data) - 1:
                retstr += '\n'
            elif i != len(self) - 1:
                retstr += ','
        return retstr + ']'

    @_exception_to_notimplemented
    def __eq__(self, other):
        if self.resolution != other.resolution:
            return False
        for x, y in zip(self, other):
            if abs(x - y) > TINY:
                return False
        return True

    def switch_range(self):
        """
        Returns a grid object whose range is 0 to 360 degrees in both dimensions
        instead of -180 to 180 degrees (or -180 to 180 degrees if the range is
        already 0 to 360 degrees)
        """
        res = self.resolution
        mid = res // 2
        newgrid = _CmapGrid(res)
        for i in range(res):
            ii = (i + mid) % res
            for j in range(res):
                jj = (j + mid) % res
                # Start from the middle
                newgrid[i, j] = self[ii, jj]
        return newgrid

    def __copy__(self):
        return _CmapGrid(self.resolution, copy(self._data))

    def __hash__(self):
        return hash(tuple([round(x, _TINY_DIGITS) for x in self]))

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class TrigonalAngle(_FourAtomTerm):
    """
    A trigonal-angle term in the AMOEBA force field. It exists in a pattern like
    the one shown below

                                A1
                                |
                                |
                          A4----A2----A3

    Parameters
    ----------
    atom1 : :class:`Atom`
        The first atom involved in the trigonal angle
    atom2 : :class:`Atom`
        The central atom involved in the trigonal angle
    atom3 : :class:`Atom`
        The third atom involved in the trigonal angle
    atom4 : :class:`Atom`
        The fourth atom involved in the trigonal angle
    type : :class:`AngleType`
        The angle type containing the parameters

    Notes
    -----
    Either `Atom`s or `Bond`s can be contained within this trigonal angle
    """
    def __init__(self, atom1, atom2, atom3, atom4, type=None):
        _FourAtomTerm.__init__(self, atom1, atom2, atom3, atom4)
        self.type = type

    def __contains__(self, thing):
        if isinstance(thing, Atom):
            return _FourAtomTerm.__contains__(self, thing)
        return ((self.atom1 in thing and self.atom2 in thing) or
                (self.atom2 in thing and self.atom3 in thing) or
                (self.atom2 in thing and self.atom4 in thing))

    def __repr__(self):
        return (
            f'<{self.__class__.__name__}; {repr(self.atom2)}--({repr(self.atom1)},'
            f'{repr(self.atom3)},{repr(self.atom4)}); type={repr(self.type)}>'
        )

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class OutOfPlaneBend(_FourAtomTerm):
    """
    Out-of-plane bending term in the AMOEBA force field. The bond pattern is the
    same as :class:`TrigonalAngle`

    Parameters
    ----------
    atom1 : :class:`Atom`
        The first atom involved in the trigonal angle
    atom2 : :class:`Atom`
        The central atom involved in the trigonal angle
    atom3 : :class:`Atom`
        The third atom involved in the trigonal angle
    atom4 : :class:`Atom`
        The fourth atom involved in the trigonal angle
    type : :class:`OutOfPlaneBendType`
        The angle type containing the parameters

    Notes
    -----
    Either `Atom`s or `Bond`s can be contained within this trigonal angle
    """
    def __init__(self, atom1, atom2, atom3, atom4, type=None):
        super().__init__(atom1, atom2, atom3, atom4)
        self.type = type

    def __contains__(self, thing):
        if isinstance(thing, Atom):
            return _FourAtomTerm.__contains__(self, thing)
        return ((self.atom1 in thing and self.atom2 in thing) or
                (self.atom2 in thing and self.atom3 in thing) or
                (self.atom2 in thing and self.atom4 in thing))

    def __repr__(self):
        return (
            f'<{self.__class__.__name__}; {repr(self.atom2)}--({repr(self.atom1)},'
            f'{repr(self.atom3)},{repr(self.atom4)}); type={repr(self.type)}>'
        )

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class OutOfPlaneBendType(_ParameterType, _ListItem):
    """
    An angle type with a set of angle parameters

    Parameters
    ----------
    k : ``float``
        Force constant in kcal/mol/radians^2
    list : :class:`TrackedList`
        A list of `OutOfPlaneBendType`s in which this is a member

    Inherited Attributes
    --------------------
    idx : ``int``
        The index of this OutOfPlaneBendType inside its containing list

    Notes
    -----
    Two `OutOfPlaneBendType`s are equal if their `k` attribute is equal

    Examples
    --------
    >>> ot1 = OutOfPlaneBendType(10.0)
    >>> ot2 = OutOfPlaneBendType(10.0)
    >>> ot1 is ot2
    False
    >>> ot1 == ot2
    True
    >>> ot1.idx # not part of any list or iterable
    -1

    As part of a list, they can be indexed

    >>> oopbend_list = []
    >>> oopbend_list.append(OutOfPlaneBendType(10.0, list=oopbend_list))
    >>> oopbend_list.append(OutOfPlaneBendType(10.0, list=oopbend_list))
    >>> oopbend_list[0].idx
    0
    >>> oopbend_list[1].idx
    1
    """
    def __init__(self, k, list=None):
        _ParameterType.__init__(self)
        self.k = _strip_units(k, u.kilocalories_per_mole/u.radians**2)
        self._idx = -1
        self.list = list

    @_exception_to_notimplemented
    def __eq__(self, other):
        return abs(self.k - other.k) < TINY

    def __repr__(self):
        return '<%s; k=%.3f>' % (type(self).__name__, self.k)

    def __copy__(self):
        return OutOfPlaneBendType(self.k)

    __getstate__ = _getstate_with_exclusions()

    def __hash__(self):
        return hash(round(self.k, _TINY_DIGITS))

    @property
    def uk(self):
        return self.k * u.kilocalories_per_mole / u.radians**2

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class PiTorsion:
    r"""
    Defines a pi-torsion term in the AMOEBA force field. The Pi-torsion is
    defined around a sp2-hybridized pi-delocalized orbital (like an amide) by 6
    atoms, as shown in the schematic below.


         A2           A5-AAA
           \         /
            \       /
             A3---A4
            /       \
           /         \
     AAA-A1           A6

    In the above schematic, A3 and A4 are sp2-hybridized, and atoms A2 and A6
    are bonded *only* to A3 and A4, respectively. Atoms A1 and A5 are each
    bonded to 3 other atoms.

    Parameters
    ----------
    atom1 : :class:`Atom`
        atom A1 in the schematic above
    atom2 : :class:`Atom`
        atom A2 in the schematic above
    atom3 : :class:`Atom`
        atom A3 in the schematic above
    atom4 : :class:`Atom`
        atom A4 in the schematic above
    atom5 : :class:`Atom`
        atom A5 in the schematic above
    atom6 : :class:`Atom`
        atom A6 in the schematic above
    type : :class:`DihedralType`
        The parameters for this Pi-torsion

    Notes
    -----
    Both :class:`Bond`s and :class:`Atom`s can be contained in a pi-torsion
    """
    def __init__(self, atom1, atom2, atom3, atom4, atom5, atom6, type=None):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.atom5 = atom5
        self.atom6 = atom6
        self.type = type

    def __contains__(self, thing):
        if isinstance(thing, Atom):
            return (thing is self.atom1 or thing is self.atom2 or
                    thing is self.atom3 or thing is self.atom4 or
                    thing is self.atom5 or thing is self.atom6)
        # Assume Bond
        return ((self.atom2 in thing and self.atom3 in thing) or
                (self.atom1 in thing and self.atom3 in thing) or
                (self.atom3 in thing and self.atom4 in thing) or
                (self.atom4 in thing and self.atom5 in thing) or
                (self.atom4 in thing and self.atom6 in thing))

    def __repr__(self):
        return (
            f'<{self.__class__.__name__}; ({repr(self.atom1)},{repr(self.atom2)})--'
            f'{repr(self.atom3)}--{repr(self.atom4)}--({repr(self.atom5)},{repr(self.atom6)}); '
            f'type={repr(self.type)}>'
        )

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class StretchBend:
    """
    This term models the stretching and bending of a standard valence angle, and
    is used in the AMOEBA force field

    Parameters
    ----------
    atom1 : :class:`Atom`
        The first atom on one end of the angle
    atom2 : :class:`Atom`
        The central atom in the angle
    atom3 : :class:`Atom`
        The atom on the other end of the angle
    type : :class:`StretchBendType`
        The type containing the stretch-bend parameters

    Notes
    -----
    Both :class:`Bond`s and :class:`Atom`s can be contained in a stretch-bend term
    """
    def __init__(self, atom1, atom2, atom3, type=None):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.type = type

    def __contains__(self, thing):
        if isinstance(thing, Atom):
            return (self.atom1 is thing or self.atom2 is thing or
                    self.atom3 is thing)
        return ((self.atom1 in thing and self.atom2 in thing) or
                (self.atom2 in thing and self.atom3 in thing))

    def __repr__(self):
        return f'<{self.__class__.__name__}; {repr(self.atom1)}--{repr(self.atom2)}--{repr(self.atom3)}; type={repr(self.type)}>'

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class StretchBendType(_ParameterType, _ListItem):
    """
    A stretch-bend type with two distances and an angle in AMOEBA

    Parameters
    ----------
    k1 : ``float``
        First force constant in kcal/mol/(radians*angstroms)
    k2 : ``float``
        Second force constant in kcal/mol/(radians*angstroms)
    req1 : ``float``
        Equilibrium bond distance for bond between the first and second atoms in
        Angstroms
    req2 : ``float``
        Equilibrium bond distance for bond between the second and third atoms in
        Angstroms
    theteq : ``float``
        Equilibrium angle in degrees
    list : :class:`TrackedList`
        A list of `StretchBendType`s in which this is a member

    Inherited Attributes
    --------------------
    idx : ``int``
        The index of this StretchBendType inside its containing list

    Notes
    -----
    Two `StretchBendType`s are equal if their `req1`, `req2`, `theteq`, and `k`
    attributes are equal

    Examples
    --------
    >>> sbt1 = StretchBendType(10.0, 10.0, 1.0, 1.0, 180.0)
    >>> sbt2 = StretchBendType(10.0, 10.0, 1.0, 1.0, 180.0)
    >>> sbt1 is sbt2
    False
    >>> sbt1 == sbt2
    True
    >>> sbt1.idx # Not part of any list or iterable
    -1

    As part of a list, they can be indexed

    >>> strbnd_list = []
    >>> strbnd_list.append(StretchBendType(10.0, 10.0, 1.0, 1.0, 180.0, strbnd_list))
    >>> strbnd_list.append(StretchBendType(10.0, 10.0, 1.0, 1.0, 180.0, strbnd_list))
    >>> strbnd_list[0].idx
    0
    >>> strbnd_list[1].idx
    1
    """
    def __init__(self, k1, k2, req1, req2, theteq, list=None):
        _ParameterType.__init__(self)
        self.k1 = _strip_units(k1, u.kilocalories_per_mole/(u.radian * u.angstroms))
        self.k2 = _strip_units(k2, u.kilocalories_per_mole/(u.radian * u.angstroms))
        self.req1 = _strip_units(req1, u.angstrom)
        self.req2 = _strip_units(req2, u.angstrom)
        self.theteq = _strip_units(theteq, u.degrees)
        self._idx = -1
        self.list = list

    @_exception_to_notimplemented
    def __eq__(self, other):
        return (abs(self.k1 - other.k1) < TINY and
                abs(self.k2 - other.k2) < TINY and
                abs(self.req1 - other.req1) < TINY and
                abs(self.req2 - other.req2) < TINY and
                abs(self.theteq - other.theteq) < TINY)

    def __repr__(self):
        return (
            f'<{self.__class__.__name__}; req1={self.req1:.3f}, req2={self.req2:.3f}, '
            f'theteq={self.theteq:.3f}, k1={self.k1:.3f}, k2={self.k2:.3f}>'
        )

    def __copy__(self):
        return StretchBendType(self.k1, self.k2, self.req1, self.req2, self.theteq)

    __getstate__ = _getstate_with_exclusions()

    def __hash__(self):
        return hash((round(self.k1, _TINY_DIGITS), round(self.k2, _TINY_DIGITS),
                     round(self.req1, _TINY_DIGITS), round(self.req2, _TINY_DIGITS),
                     round(self.theteq, _TINY_DIGITS)))

    @property
    def uk1(self):
        return self.k1 * u.kilocalories_per_mole / (u.radian * u.angstroms)

    @property
    def uk2(self):
        return self.k2 * u.kilocalories_per_mole / (u.radian * u.angstroms)

    @property
    def ureq1(self):
        return self.req1 * u.angstroms

    @property
    def ureq2(self):
        return self.req2 * u.angstroms

    @property
    def utheteq(self):
        return self.theteq * u.degrees

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class TorsionTorsion(Cmap):
    """
    This is a coupled-torsion map used in the AMOEBA force field similar to the
    correction-map (CMAP) potential used by the CHARMM force field

    Parameters
    ----------
    atom1 : :class:`Atom`
        An atom on one end of the valence torsion-torsion bonded to atom2
    atom2 : :class:`Atom`
        An atom in the middle of the torsion-torsion bonded to atoms 1 and 3
    atom3 : :class:`Atom`
        An atom in the middle of the torsion-torsion bonded to atoms 2 and 4
    atom4 : :class:`Atom`
        An atom in the middle of the torsion-torsion bonded to atoms 3 and 5
    atom5 : :class:`Atom`
        An atom in the middle of the torsion-torsion bonded to atom 4
    type : :class:`TorsionTorsionType`
        The TorsionTorsionType object containing the parameter map for this term

    Notes
    -----
    A TorsionTorsion can contain bonds or atoms. A bond is contained if it
    exists between atoms 1 and 2, between atoms 2 and 3, between atoms 3 and 4,
    or between atoms 4 and 5.

    Examples
    --------
    >>> a1, a2, a3, a4, a5 = Atom(), Atom(), Atom(), Atom(), Atom()
    >>> tortor = TorsionTorsion(a1, a2, a3, a4, a5)
    >>> Bond(a1, a2) in tortor and Bond(a2, a3) in tortor
    True
    >>> Bond(a1, a3) in tortor
    False
    """
    def __init__(self, atom1, atom2, atom3, atom4, atom5, type=None):
        Cmap.__init__(self, atom1, atom2, atom3, atom4, atom5, type)
        atom1.tortors.append(self)
        atom2.tortors.append(self)
        atom3.tortors.append(self)
        atom4.tortors.append(self)
        atom5.tortors.append(self)
        atom1.tortor_to(atom2)
        atom1.tortor_to(atom3)
        atom1.tortor_to(atom4)
        atom1.tortor_to(atom5)
        atom2.tortor_to(atom3)
        atom2.tortor_to(atom4)
        atom2.tortor_to(atom5)
        atom3.tortor_to(atom4)
        atom3.tortor_to(atom5)
        atom4.tortor_to(atom5)

    def delete(self):
        """
        Deletes this TorsionTorsion from the atoms that make it up. This method
        removes the TorsionTorsion from the `tortors` list for atom1, atom2,
        atom3, atom4, and atom5, and removes each atom from the others'
        tortor_partners list.
        """
        _delete_from_list(self.atom1.tortors, self)
        _delete_from_list(self.atom2.tortors, self)
        _delete_from_list(self.atom3.tortors, self)
        _delete_from_list(self.atom4.tortors, self)
        _delete_from_list(self.atom5.tortors, self)

        _delete_from_list(self.atom1._tortor_partners, self.atom2)
        _delete_from_list(self.atom1._tortor_partners, self.atom3)
        _delete_from_list(self.atom1._tortor_partners, self.atom4)
        _delete_from_list(self.atom1._tortor_partners, self.atom5)
        _delete_from_list(self.atom2._tortor_partners, self.atom1)
        _delete_from_list(self.atom2._tortor_partners, self.atom3)
        _delete_from_list(self.atom2._tortor_partners, self.atom4)
        _delete_from_list(self.atom2._tortor_partners, self.atom5)
        _delete_from_list(self.atom3._tortor_partners, self.atom1)
        _delete_from_list(self.atom3._tortor_partners, self.atom2)
        _delete_from_list(self.atom3._tortor_partners, self.atom4)
        _delete_from_list(self.atom3._tortor_partners, self.atom5)
        _delete_from_list(self.atom4._tortor_partners, self.atom1)
        _delete_from_list(self.atom4._tortor_partners, self.atom2)
        _delete_from_list(self.atom4._tortor_partners, self.atom3)
        _delete_from_list(self.atom4._tortor_partners, self.atom5)
        _delete_from_list(self.atom5._tortor_partners, self.atom1)
        _delete_from_list(self.atom5._tortor_partners, self.atom2)
        _delete_from_list(self.atom5._tortor_partners, self.atom3)
        _delete_from_list(self.atom5._tortor_partners, self.atom4)

        self.atom1 = self.atom2 = self.atom3 = self.atom4 = self.atom5 = None
        self.type = None

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class _TorTorTable:
    """
    Contains an interpolating potential grid for a coupled-torsion in the AMOEBA
    force field.

    Parameters
    ----------
    ang1 : ``list of floats``
        Angles in the first dimension of the interpolation table
    ang2 : ``list of floats``
        Angles in the second dimension of the interpolation table
    data : ``list of floats``
        Value of the potential grid at each point (ang2 varies fastest)

    Notes
    -----
    Raises `TypeError` if the dimension of the data array does not match the
    number of data required by ang1 and ang2

    Elements from the table are obtained and/or set by using angles as indexes.
    If the pair of angles was not one of the original angles passed to this
    table, a KeyError is raised.

    Examples
    --------
    >>> table = _TorTorTable([1.0, 2.0], [1.0, 2.0, 3.0], [1, 2, 3, 4, 5, 6])
    >>> table[1.0,1.0]
    1
    >>> table[1.0,2.0]
    2
    >>> table[2.0,2.0]
    5
    >>> table[2.0,2.0] = 10
    >>> table.data
    [1, 2, 3, 4, 10, 6]
    """
    def __init__(self, ang1, ang2, data):
        if len(data) != len(ang1) * len(ang2):
            raise TypeError(
                f'Coupled torsion parameter size mismatch. {len(ang1)}x{len(ang2)} '
                f'grid expects {len(ang1) * len(ang2)} elements (got {len(data)})'
            )
        self.data = _strip_units(data, u.kilocalories_per_mole)
        self._indexes = dict()
        i = 0
        for a1 in ang1:
            for a2 in ang2:
                self._indexes[(a1, a2)] = i
                i += 1

    def __getitem__(self, idx):
        return self.data[self._indexes[idx]]

    def __setitem__(self, idx, value):
        idx = self._indexes[_strip_units(idx, u.degrees)]
        self.data[idx] = _strip_units(value, u.kilocalories_per_mole)

    @_exception_to_notimplemented
    def __eq__(self, other):
        try:
            for idx in self._indexes.keys():
                if abs(self[idx] - other[idx]) > TINY:
                    return False
        except KeyError:
            return False
        else:
            return True

    def __ne__(self, other):
        return not self.__eq__(other)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class TorsionTorsionType(_ParameterType, _ListItem):
    """
    The type containing the parameter maps for the Amoeba torsion-torsion
    potentials. It contains the original potential as well as interpolated first
    and second derivatives for the AMOEBA force field.

    Parameters
    ----------
    dims : ``tuple of 2 ints``
        The table dimensions
    ang1 : ``list of floats``
        The list of angles in the first dimension
    ang2 : ``list of floats``
        The list of angles in the second dimension
    f : ``list of floats``
        The interpolation table for the energy
    dfda1 : ``list of floats``
        The interpolation table of the gradient w.r.t. angle 1
    dfda2 : ``list of floats``
        The interpolation table of the gradient w.r.t. angle 2
    d2fda1da2 : ``list of floats``
        The interpolation table of the 2nd derivative w.r.t. both angles
    list : :class:`TrackedList`
        The list containing this coupled torsion-torsion map

    Attributes
    ----------
    dims : ``tuple of 2 ints``
        The table dimensions
    ang1 : ``list of floats``
        The list of angles in the first dimension
    ang2 : ``list of floats``
        The list of angles in the second dimension
    f : :class:`_TorTorTable`
        The interpolation table for the energy as a _TorTorTable
    dfda1 : :class:`_TorTorTable`
        The interpolation table for the first gradient as a _TorTorTable
    dfda2 : :class:`_TorTorTable`
        The interpolation table for the second gradient as a _TorTorTable
    d2fda1da2 : :class:`_TorTorTable`
        The interpolation table for the second derivative as a _TorTorTable
    list : :class:`TrackedList`
        The list that may, or may not, contain this TorsionTorsionType
    idx : ``int``
        The index of this item in the list or iterable defined by `list`
    """
    def __init__(self, dims, ang1, ang2, f,
                 dfda1=None, dfda2=None, d2fda1da2=None, list=None):
        _ParameterType.__init__(self)
        if len(dims) != 2:
            raise ValueError('dims must be a 2-dimensional iterable')
        if len(ang1) != dims[0] or len(ang2) != dims[1]:
            raise ValueError('dims does match the angle definitions')
        self.dims = tuple(dims)
        self.ang1 = _strip_units(ang1, u.degrees)
        self.ang2 = _strip_units(ang2, u.degrees)
        self.f = _TorTorTable(ang1, ang2, f)
        if dfda1 is None:
            self.dfda1 = None
        else:
            self.dfda1 = _TorTorTable(ang1, ang2, dfda1)
        if dfda2 is None:
            self.dfda2 = None
        else:
            self.dfda2 = _TorTorTable(ang1, ang2, dfda2)
        if d2fda1da2 is None:
            self.d2fda1da2 = None
        else:
            self.d2fda1da2 = _TorTorTable(ang1, ang2, d2fda1da2)
        self._idx = -1
        self.list = list

    @_exception_to_notimplemented
    def __eq__(self, other):
        if self.dims != other.dims: return False
        if self.ang1 != other.ang1: return False
        if self.ang2 != other.ang2: return False
        return (self.f == other.f and self.dfda1 == other.dfda1 and
                self.dfda2 == other.dfda2 and self.d2fda1da2 == other.d2fda1da2)

    def __repr__(self):
        return f'<{self.__class__.__name__}; {self.dims[0]}x{self.dims[1]}>'

    def __copy__(self):
        f = copy(self.f.data)
        # dfda1
        if self.dfda1 is None:
            dfda1 = None
        else:
            dfda1 = copy(self.dfda1.data)
        # dfda2
        if self.dfda2 is None:
            dfda2 = None
        else:
            dfda2 = copy(self.dfda2.data)
        # d2fda1da2
        if self.d2fda1da2 is None:
            d2fda1da2 = None
        else:
            d2fda1da2 = copy(self.d2fda1da2.data)
        # Copy
        return TorsionTorsionType(self.dims, self.ang1, self.ang2, f, dfda1, dfda2, d2fda1da2)

    __getstate__ = _getstate_with_exclusions()

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ChiralFrame:
    """
    A chiral frame as defined in the AMOEBA force field. It defines the frame of
    reference for a chiral center

    Parameters
    ----------
    atom1 : :class:`Atom`
        The first atom defined in the chiral frame
    atom2 : :class:`Atom`
        The second atom defined in the chiral frame
    chirality : ``int``
        Either 1 or -1 to identify directionality. A ValueError is raised if a
        different value is provided

    Notes
    -----
    A chiral frame can only contain atoms.
    """
    def __init__(self, atom1, atom2, chirality):
        self.atom1 = atom1
        self.atom2 = atom2
        if chirality != 1 and chirality != -1:
            raise ValueError('chirality must be 1 or -1')
        self.chirality = chirality

    def __contains__(self, thing):
        return thing is self.atom1 or thing is self.atom2

    def __repr__(self):
        return f'<{self.__class__.__name__}; {repr(self.atom1)}--{repr(self.atom2)}, direction={self.chirality}>'

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class MultipoleFrame:
    """
    This defines the frame of reference for computing multipole interactions in
    the AMOEBA force field.

    Parameters
    ----------
    atom : :class:`Atom`
        The atom for which the frame of reference is defined
    frame_pt_num : ``int``
        The frame point number
    vectail : ``int``
        The vector tail index
    vechead : ``int``
        The vector head index
    nvec : ``int``
        The number of vectors

    Examples
    --------
    >>> atom = Atom()
    >>> mf = MultipoleFrame(atom, 0, 1, 2, 3)
    >>> atom in mf
    True
    >>> mf.frame_pt_num
    0
    >>> mf.vectail
    1
    >>> mf.vechead
    2
    >>> mf.nvec
    3
    """
    def __init__(self, atom, frame_pt_num, vectail, vechead, nvec):
        self.atom = atom
        self.frame_pt_num = frame_pt_num
        self.vectail = vectail
        self.vechead = vechead
        self.nvec = nvec

    def __contains__(self, thing):
        return self.atom is thing

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class DrudeAtom(Atom):
    """
    An Atom that has a Drude particle attached to it.  This is a subclass of
    Atom, so it also has all the properties defined for regular Atoms.

    Parameters
    ----------
    alpha : ``float``
        the atomic polarizability
    thole : ``float``
        the Thole damping facior
    drude_type : ``str``
        the atom type to use for the Drude particle.

    Other Attributes
    ----------------
    anisotropy : :class:`DrudeAnisotropy`
        describes how this atom is anisotropically polarizable.  For isotropic
        atoms, this is None.
    """

    def __init__(self, alpha=0.0, thole=1.3, drude_type='DRUD', parent=None, **kwargs):
        super().__init__(**kwargs)
        self.alpha = alpha
        self.thole = thole
        self.drude_type = drude_type
        self.anisotropy: Optional["DrudeAnisotropy"] = None
        self.parent = parent

    @property
    def drude_charge(self):
        sign = (-1 if self.alpha < 0 else 1)
        alpha = abs(self.alpha)*u.angstrom**3/(138.935456*u.kilojoules_per_mole*u.nanometer)
        return sign*math.sqrt(alpha*2*(500*u.kilocalories_per_mole/u.angstrom**2))

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class DrudeAnisotropy:
    """
    A description of an anisotropically polarizable atom.

    Atom 1 is a :class:`DrudeAtom` whose polarizability is anisotropic.  The
    other three atoms define the coordinate frame.

    Parameters
    ----------
    atom1 : :class:`Atom`
        the parent atom to the Drude atom
    atom2 : :class:`Atom`
        the second atom defining the coordinate frame
    atom3 : :class:`Atom`
        the third atom defining the coordinate frame
    atom4 : :class:`Atom`
        the fourth atom defining the coordinate frame
    a11 : ``float``
        the scale factor for the polarizability along the direction defined by atom1 and atom2
    a22 : ``float``
        the scale factor for the polarizability along the direction defined by atom3 and atom4
    **params : Any
        arbitrary parameter list that can be used to store additional information about this
        Drude anisotropy (useful for supporting PSF writes)
    """

    def __init__(self, atom1, atom2, atom3, atom4, a11, a22, **params):
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.a11 = a11
        self.a22 = a22
        self.params = params

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Residue(_ListItem):
    """
    A single residue that is composed of a small number of atoms

    Parameters
    ----------
    name : ``str``
        Name of the residue. Typical convention is to choose a name that is 4
        characters or shorter
    number : ``int``, optional
        Residue number assigned in the input structure. Default is -1
    chain : ``str``, optional
        The 1-letter chain identifier for this residue. Default is empty string
    insertion_code : ``str``, optional
        The insertion code (used in PDB files) for this residue. Default is
        empty string
    segid : ``str``, optional
        The segment identifier, used by CHARMM in a way similar to chain. Dfault
        is empty string
    list : :class:`TrackedList`
        List of residues in which this residue is a member

    Attributes
    ----------
    name : ``str``
        The name of this residue
    number : ``int``
        The number of this residue in the input structure
    idx : ``int``
        The index of this residue inside the container. If this residue has no
        container, or it is not present in the container, idx is -1
    chain : ``str``
        The 1-letter chain identifier for this residue
    insertion_code : ``str``
        The insertion code (used in PDB files) for this residue
    ter : ``bool``
        If True, there is a TER card directly after this residue (i.e., a
        molecule or chain ends). By default, it is False
    list : :class:`TrackedList`
        The container that _may_ have this residue contained inside
    atoms : ``list of`` :`class`Atom` ``instances``
        This is the list of `Atom`s that make up this residue

    Notes
    -----
    - Iterating over a residue will iterate over the atoms. It is exactly
      equivalent to iterating over the `atoms` attribute
    - Supports testing if an Atom instance is contained `in` this residue
    - `len()` returns the number of atoms in this residue
    """

    def __init__(self, name, number=-1, chain='', insertion_code='',
                 segid='', list=None):
        self.name = name.strip()
        self.number = number
        self.chain = chain.strip()
        self.insertion_code = insertion_code.strip()
        self.list = list
        self._idx = -1
        self.atoms = []
        self.ter = False
        self.segid = segid

    def add_atom(self, atom):
        """ Adds an atom to this residue

        Parameters
        ----------
        atom : :class:`Atom`
            The atom to add to this residue

        Notes
        -----
        This action assigns the `residue` attribute to `atom`
        """
        atom.residue = self
        self.atoms.append(atom)

    def delete_atom(self, atom):
        """
        If an atom is present in this residue, delete it from the list of
        atoms. No change if an atom is not present in this residue.
        """
        for a in self.atoms:
            if atom.residue is self:
                atom.residue = None
                self.atoms = [a for a in self.atoms if a is not atom]

    # Implement some container methods over the list of atoms
    def __contains__(self, thing):
        """ True if an atom is present in this residue """
        return thing in self.atoms

    def __len__(self):
        return len(self.atoms)

    def __iter__(self):
        return iter(self.atoms)

    def is_empty(self):
        """
        Determines if there are any atoms in this residue

        Returns
        -------
        empty: ``bool``
            ``True`` if there are no atoms left. ``False`` otherwise.
        """
        return len(self) == 0

    def sort(self):
        """ Sorts the atoms in this list by atom index """
        self.atoms.sort()

    def __getitem__(self, idx):
        return self.atoms.__getitem__(idx)

    # Sort by atom indices
    def __lt__(self, other):
        return self.atoms[0].idx < other.atoms[0].idx
    def __gt__(self, other):
        return self.atoms[0].idx > other.atoms[0].idx
    def __le__(self, other):
        return not self.atoms[0].idx > other.atoms[0].idx
    def __ge__(self, other):
        return not self.atoms[0].idx < other.atoms[0].idx

    def __repr__(self):
        if self.number == -1:
            num = self.idx
        else:
            num = self.number
        rep = f'<{self.__class__.__name__} {self.name}[{num}]'
        if self.chain:
            rep += f'; chain={self.chain}'
        if self.insertion_code:
            rep += f'; insertion_code={self.insertion_code}'
        if self.segid:
            rep += f'; segid={self.segid}'
        return rep + '>'

    __getstate__ = _getstate_with_exclusions()

    def __setstate__(self, d):
        self.__dict__.update(d)
        for a in self:
            a.residue = self

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def _changes(func):
    """ Decorator to indicate the list has changed """
    def new_func(self, *args, **kwargs):
        self.changed = True
        self.needs_indexing = True
        return func(self, *args, **kwargs)
    return new_func

class TrackedList(list):
    """
    This creates a list type that allows you to see if anything has changed

    Attributes
    ----------
    changed : ``bool``
        Determines if something has been done to fundamentally change the
        underlying topology defined by this list such that the topology needs to
        be rebuilt
    needs_indexing : ``bool``
        A flag to determine whether or not the items in a tracked list need to
        be indexed or not.

    Examples
    --------
    >>> tl = TrackedList()
    >>> tl.append(Atom())
    >>> tl.append(Atom())
    >>> tl.append(Atom())
    >>> tl.needs_indexing, tl.changed
    (True, True)
    >>> tl.index_members()
    >>> tl.needs_indexing, tl.changed
    (False, True)
    >>> tl.changed = False # Must do when changes have been incorporated
    >>> tl.needs_indexing, tl.changed
    (False, False)
    """
    def __init__(self, *args):
        self.changed = False
        self.needs_indexing = False
        return list.__init__(self, *args)

    def __repr__(self):
        retstr = ['%s([\n' % (type(self).__name__)]
        if len(self) > 30:
            retstr.extend(f'\t{repr(self[i])}\n' for i in range(24))
            retstr.append('\t...\n')
            retstr.extend(f'\t{repr(self[i])}\n' for i in range(-5, 0))
        else:
            retstr.extend(f'\t{repr(i)}\n' for i in self)
        retstr.append('])')
        return ''.join(retstr)

    @_changes
    def __delitem__(self, item):
        """ Deletes items and slices. Make sure all items """
        try:
            indices = range(*item.indices(len(self)))
        except AttributeError:
            indices = [item]

        for index in indices:
            try:
                self[index]._idx = -1
            except AttributeError:
                # If we can't set _idx attribute on this object, don't fret
                pass
            try:
                self[index].list = None
            except AttributeError:
                # If we can't set list attribute on this object, don't fret
                pass

        return list.__delitem__(self, item)

    @_changes
    def __delslice__(self, start, stop):
        """ Python 2 still uses __delslice__... """
        self.__delitem__(slice(start, stop)) # pragma: no cover

    @_changes
    def pop(self, idx=-1):
        item = list.pop(self, idx)
        if hasattr(item, '_idx'):
            item._idx = -1
        return item

    @_changes
    def remove(self, thing):
        list.remove(self, thing)
        # If this did not raise an exception, it was part of the list, so
        # de-index it
        if hasattr(thing, '_idx'):
            thing._idx = -1

    append = _changes(list.append)
    extend = _changes(list.extend)
    insert = _changes(list.insert)
    __setitem__ = _changes(list.__setitem__)
    __iadd__ = _changes(list.__iadd__)
    __imul__ = _changes(list.__imul__)

    # Type-safe methods that return another instance
    def __add__(self, other):
        return TrackedList(list.__add__(self, other))

    def __mul__(self, fac):
        return TrackedList(list.__mul__(self, fac))

    def index_members(self):
        """
        Assigns the idx variable for every member of this list to its place in
        the list, if the members of this list permit
        """
        for i, item in enumerate(self):
            try:
                item._idx = i
            except AttributeError:
                # Must be some kind of immutable type... don't worry
                pass
        self.needs_indexing = False

    def claim(self):
        """
        This method causes this list to "claim" all of the items it contains and
        subsequently indexes all of its items.
        """
        for i, item in enumerate(self):
            try:
                item.list = self
            except AttributeError:
                # Must be some kind of immutable type... don't worry
                pass
        self.index_members()

    def prune_unused(self):
        """
        This method inspects the `used` attribute of all of its members, if it
        has one, and deletes any item in which it is set to `False`
        """
        for i in reversed(range(len(self))):
            try:
                if not self[i].used:
                    del self[i]
            except AttributeError:
                # Don't worry if we don't have a `used` attribute
                pass

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ResidueList(TrackedList):
    """ Array of `Residue` instances """

    def add_atom(self, atom, resname, resnum, chain='', inscode='', segid=''):
        """
        Adds a new atom to the ResidueList, adding a new residue to this list if
        it has a different name or number as the last residue

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
        inscode = inscode.strip()
        segid = segid.strip()
        try:
            last = self[-1]
        except IndexError:
            # Empty list -- add our first residue
            new_res = Residue(resname, resnum, chain, inscode, segid, list=self)
            new_res.add_atom(atom)
            self.append(new_res)
        else:
            if (last.number != resnum or last.name != resname.strip()
                or last.chain != chain.strip()
                or last.segid != segid.strip()
                or last.insertion_code != inscode.strip()):
                new_res = Residue(resname, resnum, chain, inscode, segid, list=self)
                new_res.add_atom(atom)
                self.append(new_res)
            else:
                last.add_atom(atom)

    def prune(self):
        """
        This function goes through the residue list and removes all empty
        residues from the list. This isn't done automatically when atoms are
        deleted, since it will become very slow. You must remember to do this to
        avoid including empty residues
        """
        # Delete from the back to avoid indexes changing as we iterate
        for i in reversed(range(len(self))):
            res = self[i]
            if res.is_empty():
                del self[i]

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AtomList(TrackedList):
    """
    Array of Atoms

    Notes
    -----
    Deleting an atom from the AtomList also deletes that atom from the residue
    it belongs to.
    """

    @_changes
    def __delitem__(self, idx):
        """ Deleting an atom also needs to delete it from the residue """
        try:
            indices = range(*idx.indices(len(self)))
        except AttributeError:
            indices = [idx]

        for index in indices:
            atom = self[index]
            atom._idx = -1
            atom.list = None
            # Make sure we delete this atom from its respective residue
            if atom.residue is not None: atom.residue.delete_atom(atom)

        list.__delitem__(self, idx)

    @_changes
    def pop(self, idx=-1):
        atom = list.pop(self, idx)
        atom._idx = -1
        atom.list = None
        if atom.residue is not None: atom.residue.delete_atom(atom)
        return atom

    @_changes
    def remove(self, atom):
        if atom.list is not self:
            raise ValueError(f'{repr(atom)} is not in list')
        if atom.residue is not None: atom.residue.delete_atom(atom)
        atom._idx = -1
        atom.list = None
        list.remove(self, atom)

    def unmark(self):
        """ Unmark all atoms in this list """
        for atm in self: atm.marked = 0

    @_changes
    def append(self, item):
        """
        Add an Atom to the end of the list and have this list claim ownership of
        the item.

        Parameters
        ----------
        item : :class:`Atom`
            The atom to add to this list

        Notes
        -----
        Only Atom objects should be added here, so if `item` does not have a
        `list` attribute, an AttributeError will be raised. This action assigns
        this list as the `list` attribute to the passed Atom.
        """
        item.list = self
        return list.append(self, item)

    @_changes
    def extend(self, items):
        """
        Add an iterable of `Atom`s to the end of the list and have this list
        claim ownership of the items.

        Parameters
        ----------
        items : iterable of :class:`Atom`s
            The iterable containing Atom instances to add to the end of this
            list

        Notes
        -----
        Any generator passed here will be exhausted. The `list` attribute for
        every object in `items` will have their `list` attribute set to this
        list, and an AttributeError will be raised if this is not possible.
        """
        for item in items:
            item.list = self
        return list.extend(self, items)

    @_changes
    def insert(self, idx, item):
        """
        Insert an Atom into the atom list

        Parameters
        ----------
        idx : ``int``
            The index in front of (i.e., before) which to insert the item
        item : :class:`Atom`
            The atom to insert in the desired index. This atom will be claimed
            by the AtomList
        """
        item.list = self
        return list.insert(self, idx, item)

    def assign_nbidx_from_types(self):
        """
        Assigns the nb_idx attribute of every atom inside here from the
        atom_type definition. If the atom_type is not assigned, RuntimeError is
        raised.

        Returns
        -------
        ``list of dict``
            Each element is a `set` of the `nb_idx` indices for which NBFIX
            alterations are defined for the type with that given that index
            (minus 1, to adjust for indexing from 0 and nb_idx starting from 1)
        """
        idx = 1
        natoms = len(self)
        atom_type_lookups = dict() # For fast lookups
        atom_type_list = []
        for i, atom in enumerate(self):
            if atom.atom_type is UnassignedAtomType:
                raise RuntimeError('atom types are not assigned')
            atom.atom_type._idx = -1

        for i, atom in enumerate(self):
            type1 = atom.atom_type
            # Skip atom types that have already been assigned
            if type1._idx != -1: continue
            type1._idx = idx
            atom_type_lookups[str(type1)] = type1
            atom_type_list.append(type1)
            for j in range(i+1, natoms):
                atom2 = self[j]
                type2 = atom2.atom_type
                # Skip atom types that have already been assigned
                if type2._idx != -1: continue
                if type1 == type2:
                    type2._idx = idx
            idx += 1
        # Now go back through and assign nb_idx from type._idx
        for atom in self:
            atom.nb_idx = atom.atom_type._idx
        # Now collect the nbfixes
        nbfix_list = [set() for i in range(idx-1)]
        # Look through all of the atom types and add the nbfixes
        for i, type in enumerate(atom_type_list):
            for key in type.nbfix:
                rmin, eps, rmin14, eps14 = type.nbfix[key]
                try:
                    otype = atom_type_lookups[key]
                except KeyError:
                    continue
                else:
                    obj = (otype._idx, rmin, eps, rmin14, eps14)
                    nbfix_list[i].add(obj)

        return nbfix_list

    def find_original_index(self, idx):
        """
        Finds an atom with the given original index. Cannot assume that the
        original indexes are in order, since reordering may have been necessary.
        As a result, the complexity of this algorithm is O(N)

        Parameters
        ----------
        idx : int
            The integer corresponding to the original index

        Returns
        -------
        atom : :class:`Atom`
            The atom with the original index ``idx``

        Raises
        ------
        IndexError
            If no atom has the original index ``idx``
        """
        for atom in self:
            if atom.number == idx: return atom
        raise IndexError(f'No atom found with index {idx}')

    def __iadd__(self, other):
        return NotImplemented

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class NonbondedException:
    """
    The AMOEBA force field has complex exclusion and exception rules (referred
    to as "adjustments" in the Amber-converted files). This class stores
    per-particle exceptions.

    Parameters
    ----------
    atom1 : :class:`Atom`
        One of the atoms in the exclusion pair
    atom2 : :class:`Atom`
        The other atom in the exclusion pair
    type : :class:`NonbondedExceptionType`
        The nonbonded exception type that describes how the various nonbonded
        interactions between these two atoms should work

    Notes
    -----
    NonbondedException objects "contain" two atoms and will return True when
    used with the binary `in` operator
    """
    def __init__(self, atom1, atom2, type=None):
        self.atom1 = atom1
        self.atom2 = atom2
        self.type = type
        self.funct = 1

    def __contains__(self, thing):
        return thing is self.atom1 or thing is self.atom2

    def __repr__(self):
        retstr = [f'<{self.__class__.__name__}; {repr(self.atom1)} and {repr(self.atom2)}']
        if self.type is not None:
            retstr.append(f', type={repr(self.type)}>')
        else:
            retstr.append('>')
        return ''.join(retstr)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class NonbondedExceptionType(_ParameterType, _ListItem):
    """
    A parameter describing how the various nonbonded interactions between a
    particular pair of atoms behaves in a specified nonbonded exception (e.g.,
    in 1-4 interacting terms)

    Parameters
    ----------
    rmin : float
        The combined Rmin value for this particular pair of atom types
        (dimension length, default units are Angstroms)
    epsilon : float
        The combined well-depth value for this particular pair of atom types
        (dimension energy, default units are kcal/mol)
    chgscale : float, optional
        The scaling factor by which to multiply the product of the charges for
        this pair. Default is 1.0.
    list : :class:`TrackedList`
        The list containing this nonbonded exception
    """

    def __init__(self, rmin, epsilon, chgscale=1.0, list=None):
        _ParameterType.__init__(self)
        self.rmin = _strip_units(rmin, u.angstroms)
        self.epsilon = _strip_units(epsilon, u.kilocalories_per_mole)
        self.chgscale = chgscale
        self._idx = None
        self.list = list

    @property
    def sigma(self):
        return self.rmin * 2**(-1/6)

    @sigma.setter
    def sigma(self, value):
        self.rmin = value * 2**(1/6)

    @property
    def usigma(self):
        return self.sigma * u.angstroms

    @property
    def urmin(self):
        return self.rmin * u.angstroms

    @property
    def uepsilon(self):
        return self.epsilon * u.kilocalories_per_mole

    def __repr__(self):
        return f'<{self.__class__.__name__}; rmin={self.rmin:.4f}, epsilon={self.epsilon:.4f}, chgscale={self.chgscale:.4f}>'

    @_exception_to_notimplemented
    def __eq__(self, other):
        return (abs(self.rmin - other.rmin) < TINY and abs(self.epsilon - other.epsilon) < TINY and
                abs(self.chgscale - other.chgscale) < TINY)

    __getstate__ = _getstate_with_exclusions()

    def __hash__(self):
        return hash((round(self.rmin, _TINY_DIGITS), round(self.epsilon, _TINY_DIGITS),
                     round(self.chgscale, _TINY_DIGITS)))

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AmoebaNonbondedExceptionType(NonbondedExceptionType):
    """
    A parameter describing how the various nonbonded interactions between a
    particular pair of atoms is scaled in the AMOEBA force field

    Parameters
    ----------
    vdw_weight : ``float``
        The scaling factor by which van der Waals interactions are multiplied
    multipole_weight : ``float``
        The scaling factor by which multipole interactions are multiplied
    direct_weight : ``float``
        The scaling factor by which direct-space interactions are multiplied
    polar_weight : ``float``
        The scaling factor by which polarization interactions are multiplied
    mutual_weight : ``float``
        The scaling factor by which mutual interactions are multiplied
    list : :class:`TrackedList`
        The list containing this nonbonded exception

    Other Attributes
    ----------------
    idx : ``int``
        The index of this term in the list that contains it
    """
    def __init__(self, vdw_weight, multipole_weight, direct_weight,
                 polar_weight, mutual_weight, list=None):
        self.vdw_weight = vdw_weight
        self.multipole_weight = multipole_weight
        self.direct_weight = direct_weight
        self.polar_weight = polar_weight
        self.mutual_weight = mutual_weight
        self._idx = -1
        self.list = list

    @_exception_to_notimplemented
    def __eq__(self, other):
        return (abs(self.vdw_weight - other.vdw_weight) < TINY and
                abs(self.multipole_weight - other.multipole_weight) < TINY and
                abs(self.direct_weight - other.direct_weight) < TINY and
                abs(self.polar_weight - other.polar_weight) < TINY and
                abs(self.mutual_weight - other.mutual_weight) < TINY)

    def __copy__(self):
        return AmoebaNonbondedExceptionType(
            self.vdw_weight, self.multipole_weight, self.direct_weight, self.polar_weight,
            self.mutual_weight
        )

    __getstate__ = _getstate_with_exclusions()

    def __hash__(self):
        return hash((round(self.vdw_weight, _TINY_DIGITS),
                     round(self.multipole_weight, _TINY_DIGITS),
                     round(self.direct_weight, _TINY_DIGITS),
                     round(self.polar_weight, _TINY_DIGITS),
                     round(self.mutual_weight, _TINY_DIGITS)))

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AtomType:
    """
    Atom types can either be compared by indexes or names. Can be assigned with
    a string, integer, (string is not automatically cast) or with both.

    Parameters
    ----------
    name : ``str``
        The name of the atom type
    number : ``int``
        The serial index of the atom type
    mass : ``float``
        The mass of the atom type
    atomic_number : ``int``, optional
        The atomic number of the element of the atom type. Default -1
    bond_type : ``str``, optional
        If defined, this is the type name used to look up bonded parameters.
        Default is None (which falls back to ``name``)
    charge : ``float``, optional
        If defined, this is the partial atomic charge in elementary charge
        units. Default is None

    Other Attributes
    ----------------
    epsilon : ``float``
        If set, it is the Lennard-Jones well depth of this atom type
    rmin : ``float``
        If set, it is the Lennard-Jones Rmin/2 parameter of this atom type
    epsilon_14 : ``float``
        If set, it is the Lennard-Jones well depth of this atom type in 1-4
        nonbonded interactions
    rmin_14 : ``float``
        If set, it is the Lennard-Jones Rmin/2 parameter of this atom type in
        1-4 nonbonded interactions
    sigma : ``float``
        This is the sigma parameter, which is just equal to Rmin*2^(1/6)
    sigma_14 : ``float``
        This is the sigma parameter corresponding to rmin_14, which is just equal to Rmin_14*2^(1/6)
    nbfix : ``dict(str:tuple)``
        A hash that maps atom type names of other atom types with which _this_
        atom type has a defined NBFIX with a tuple containing the terms
        (Rmin, epsilon, Rmin14, Epsilon14)
    _bond_type : str or None
        If an explicit value was given to bond_type, _bond_type will be set to a
        value. Otherwise, _bond_type will be None (and bond_type will be the
        same as ``name``). This can be used to determine if a bond_type was
        specified by comparing to None.

    Notes
    -----
    This object is primarily used to build parameter databases from parameter
    files. Also, sigma is related to Rmin, but rmin is Rmin/2, so there is an
    extra factor of 2 in the sigma for this reason.

    Examples
    --------
    >>> at = AtomType('HA', 1, 1.008, 1)
    >>> at.name, at.number
    ('HA', 1)
    >>> at2 = AtomType('CA', 2, 12.01, 6)
    >>> at2.name, at2.number
    ('CA', 2)
    >>> print("%s: %d" % (str(at), int(at)))
    HA: 1
    """
    _digits_tol = 6

    def __init__(self, name, number, mass, atomic_number=-1, bond_type=None, charge=0.0):
        if number is None and name is not None:
            # If we were given an integer, assign it to number. Otherwise,
            # assign it to the name
            if isinstance(name, int):
                self.number = name
                self.name = None
            else:
                self.name = name
                self.number = None
        else:
            self.name = name
            self.number = int(number)
        self.mass = _strip_units(mass, u.daltons)
        self.atomic_number = atomic_number
        # We have no LJ parameters as of yet
        self.epsilon = self.rmin = self.epsilon_14 = self.rmin_14 = None
        # Store each NBFIX term as a dict with the atom type string matching to
        # a 2-element tuple that is rmin, epsilon
        self.nbfix = dict()
        # Also store the NBTHOLE term similarly, except it is just a single number
        self.nbthole = dict()
        self._idx = -1 # needed for some internal bookkeeping
        self._bond_type = bond_type
        self.charge = charge

    def __repr__(self):
        return (f"<{self.__class__.__name__} {self.name} [# {self.number}]; q={self.charge:.4f} rmin={self.rmin:.4f} "
                f"epsilon={self.epsilon:.4f} rmin_14={self.rmin_14:.4f} epsilon_14={self.epsilon_14:.4f}>")

    @_exception_to_notimplemented
    def __eq__(self, other):
        """
        Compares based on available properties (name and number, just name,
        or just number)
        """
        if isinstance(other, AtomType):
            if self.name != other.name or self.number != other.number:
                return False
            # Now check if LJ parameters are defined, and make sure those are
            # also equal
            has_all = True
            has_none = True
            for attr in ('epsilon', 'rmin', 'epsilon_14', 'rmin_14'):
                if getattr(self, attr) is None and getattr(other, attr) is None:
                    has_all = False
                    continue
                elif getattr(self, attr) is None or getattr(other, attr) is None:
                    return False # can't be equal
                else:
                    has_none = False
            assert has_all ^ has_none, 'Please report this bug!'
            if not has_all:
                return True
            # Check charges
            if self.charge is None or other.charge is None:
                if self.charge is not other.charge:
                    return False
            elif abs(self.charge - other.charge) > TINY:
                return False
            # At this point, we have all the attributes we need to compare
            return (abs(self.epsilon - other.epsilon) < TINY and
                    abs(self.rmin - other.rmin) < TINY and
                    abs(self.epsilon_14 - other.epsilon_14) < TINY and
                    abs(self.rmin_14 - other.rmin_14) < TINY and
                    self.nbfix == other.nbfix)
        if isinstance(other, str):
            return self.name == other
        if isinstance(other, int):
            return self.number == other
        return other == (self.number, self.name)

    def set_lj_params(self, eps, rmin, eps14=None, rmin14=None):
        """ Sets Lennard-Jones parameters on this atom type """
        if u.is_quantity(eps):
            eps = eps.value_in_unit(u.kilocalories_per_mole)
        if u.is_quantity(rmin):
            rmin = rmin.value_in_unit(u.angstroms)
        if eps14 is None:
            eps14 = eps
        if rmin14 is None:
            rmin14 = rmin
        self.epsilon = eps
        self.rmin = rmin
        self.epsilon_14 = eps14
        self.rmin_14 = rmin14

    @property
    def bond_type(self):
        if self._bond_type is None:
            return self.name
        return self._bond_type

    @bond_type.setter
    def bond_type(self, value):
        self._bond_type = value

    def __int__(self):
        """ The integer representation of an AtomType is its index """
        return self.number

    def add_nbfix(self, typename, rmin, epsilon, rmin14=None, epsilon14=None):
        """ Adds a new NBFIX exclusion for this atom type

        Parameters
        ----------
        typename : str
            The name of the *other* type with which this NBFIX is defined
        rmin : float
            The combined Rmin value for this NBFIXed pair. If no units, assumed
            to be in Angstroms
        epsilon : float
            The combined epsilon value for this NBFIXed pair. If no units,
            assumed to be in kcal/mol
        rmin14 : float, optional
            Same as rmin, but for 1-4 interactions. If None (default), it is
            given the same value as rmin
        epsilon14 : float, optional
            Same as epsilon, but for 1-4 interactions. If None (default), it is
            given the same value as rmin
        """
        if rmin14 is None: rmin14 = rmin
        if epsilon14 is None: epsilon14 = epsilon
        self.nbfix[typename] = (rmin, epsilon, rmin14, epsilon14)

    def add_nbthole(self, typename, nbt):
        """ Adds a new NBTHOLE screening factor for this atom type """
        self.nbthole[typename] = nbt

    @property
    def sigma(self):
        """ Sigma is Rmin / 2^(1/6) """
        return self.rmin * 2**(-1/6) * 2

    @sigma.setter
    def sigma(self, value):
        self.rmin = value * 2**(1/6) / 2

    @property
    def usigma(self):
        return self.sigma * u.angstroms

    @property
    def sigma_14(self):
        """ Sigma is Rmin / 2^(1/6) """
        return self.rmin_14 * 2**(-1/6) * 2

    @sigma_14.setter
    def sigma_14(self, value):
        self.rmin_14 = value * 2**(1/6) / 2

    @property
    def usigma_14(self):
        return self.sigma_14 * u.angstroms

    @property
    def urmin(self):
        return self.rmin * u.angstroms

    @property
    def uepsilon(self):
        return self.epsilon * u.kilocalories_per_mole

    @property
    def urmin_14(self):
        return self.rmin_14 * u.angstroms

    @property
    def uepsilon_14(self):
        return self.epsilon_14 * u.kilocalories_per_mole

    def __str__(self):
        return self.name

    def __copy__(self):
        cp = AtomType(self.name, self.number, self.mass, self.atomic_number,
                      bond_type=self._bond_type, charge=self.charge)
        cp.epsilon = self.epsilon
        cp.rmin = self.rmin
        cp.epsilon_14 = self.epsilon_14
        cp.rmin_14 = self.rmin_14
        cp.nbfix = self.nbfix.copy()
        return cp

    def __hash__(self):
        return hash((self.name, self._round_trunc(self.mass), self.atomic_number, self.bond_type,
                     self._round_trunc(self.charge), self._round_trunc(self.epsilon),
                     self._round_trunc(self.rmin), self._round_trunc(self.epsilon_14),
                     self._round_trunc(self.rmin_14), tuple(self.nbfix.items())))

    @classmethod
    def _round_trunc(cls, value):
        return round(value, _TINY_DIGITS) if value is not None else None

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class _UnassignedAtomType:
    """
    This raises the appropriate exceptions (ParameterError) when you try to
    access its properties
    """
    _OBJ = None

    def __new__(cls):
        if cls._OBJ is None:
            cls._OBJ = super(_UnassignedAtomType, cls).__new__(cls)
        return cls._OBJ

    def __int__(self):
        raise ParameterError('Atom type is not defined')

    def __str__(self):
        raise ParameterError('Atom type is not defined')

    def __eq__(self, other):
        return isinstance(other, _UnassignedAtomType) # Behave like a singleton

    def __reduce__(self):
        return 'UnassignedAtomType'

UnassignedAtomType = _UnassignedAtomType()
# Make sure it's a singleton
assert UnassignedAtomType is _UnassignedAtomType(), "Not a singleton"

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class AcceptorDonor:
    """ Just a holder for donors and acceptors in CHARMM speak

    Parameters
    ----------
    atom1 : :class:`Atom`
        First atom in the donor/acceptor group
    atom2 : :class:`Atom`
        Second atom in the donor/acceptor group
    """
    def __init__(self, atom1, atom2):
        self.atom1 = atom1
        self.atom2 = atom2

    def __repr__(self):
        return f'<AcceptorDonor; {repr(self.atom1)} {repr(self.atom2)}>'

    def __contains__(self, thing):
        """ See if the atom is in this donor/acceptor """
        return thing is self.atom1 or thing is self.atom2

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Group:
    """ An 'interacting' group defined by CHARMM PSF files

    Parameters
    ----------
    atom : :class:`Atom`
        The first atom within a group
    type : ``int``
        Flag for group information; 0 when all atoms have zero charge,
        1 when group has a net zero charge but at least one atom has a non-zero
        partial charge, 2 when the net charge of the group is not zero
    move : ``int``
        0 if the atoms are not fixed, 1 when they are

    Notes
    -----
    See the discussion on Github for the source of the meanings of these
    variables: https://github.com/ParmEd/ParmEd/pull/307#issuecomment-128244134
    """
    def __init__(self, atom, type, move):
        self.atom = atom
        self.type = type
        self.move = move

    @_exception_to_notimplemented
    def __eq__(self, other):
        return (self.atom is other.atom and self.type == other.type and self.move == other.move)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Link:
    """ An intra-residue "Link" as defined by the PDB standard:

    See http://www.wwpdb.org/documentation/file-format-content/format33/sect6.html#LINK for more
    information

    Parameters
    ----------
    atom1 : :class:`Atom`
        The first Atom involved in the Link
    atom2 : :class:`Atom`
        The other atom to which ``atom1`` is bonded in this link
    length : float
        The length of the link
    symmetry_op1 : str, optional
        The first symmetry operator for the link
    symmetry_op2 : str, optional
        The second symmetry operator for the link
    """

    def __init__(self, atom1, atom2, length, symmetry_op1='1555', symmetry_op2='1555'):
        self.atom1 = atom1
        self.atom2 = atom2
        self.length = length
        self.symmetry_op1 = symmetry_op1
        self.symmetry_op2 = symmetry_op2

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

NoUreyBradley = BondType(0.0, 0.0) # singleton representing lack of a U-B term
