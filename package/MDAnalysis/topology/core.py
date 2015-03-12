# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
"""
Common functions for topology building --- :mod:`MDAnalysis.topology.core`
==========================================================================

The various topology parsers make use of functions and classes in this
module. They are mostly of use to developers.

.. SeeAlso:: :mod:`MDAnalysis.topology.tables` for some hard-coded atom
   information that is used by functions such as :func:`guess_atom_type` and
   :func:`guess_atom_mass`.

"""
from __future__ import print_function
# Global imports
import os.path
from math import sqrt
import numpy
from collections import defaultdict
from itertools import izip

# Local imports
from . import tables
from ..core import distances
from ..core.util import norm, dihedral, cached
from ..core.util import angle as slowang
from ..core import AtomGroup


def build_segments(atoms):
    """Create all :class:`~MDAnalysis.core.AtomGroup.Segment` instancess from
    a list of :class:`~MDAnalysis.core.AtomGroup.Atom` instances.

    The function also builds the :class:`~MDAnalysis.core.AtomGroup.Residue`
    instances by tracking residue numbers.

    Updating segments also changes the underlying
    :class:`~MDAnalysis.core.AtomGroup.Atom` instances, which record
    to which residue and segment an atom belongs.

    :Returns: structure dict, which associates a segname with a
              :class:`~MDAnalysis.core.AtomGroup.Segment`

    .. versionchanged:: 0.9.0
       Now allows resids in a given Segment to be given in non sequential order.
    """
    struc = {}
    resatomlist = defaultdict(list)
    curr_segname = atoms[0].segid

    for a in atoms:
        if a.segid == curr_segname:  # if still in same Segment
            resatomlist[a.resid].append(a)
        else:
            # We've come to a new segment
            # Build the Segment we just left
            residues = [AtomGroup.Residue(ats[0].resname, k, ats)
                        for k, ats in resatomlist.iteritems()]
            struc[curr_segname] = AtomGroup.Segment(curr_segname, residues)

            # Reset things and start again
            resatomlist = defaultdict(list)
            resatomlist[a.resid].append(a)
            curr_segname = a.segid

    # Add the last segment
    residues = [AtomGroup.Residue(ats[0].resname, k, ats)
                for k, ats in resatomlist.iteritems()]
    struc[curr_segname] = AtomGroup.Segment(curr_segname, residues)
    return struc


def build_residues(atoms):
    """Create a list :class:`~MDAnalysis.core.AtomGroup.Residue` instances from
    a list of :class:`~MDAnalysis.core.AtomGroup.Atom` instances.

    Updating residues also changes the underlying
    :class:`~MDAnalysis.core.AtomGroup.Atom` instances, which record
    to which residue an atom belongs.

    :Returns: List of :class:`~MDAnalysis.core.AtomGroup.Residue` instances

    .. versionadded:: 0.8
    .. versionchanged:: 0.9.0
       Now allows resids to be given in non sequential order
    """
    resatomlist = defaultdict(list)
    for a in atoms:
        resatomlist[a.resid].append(a)

    residues = [AtomGroup.Residue(ats[0].resname, k, ats)
                for k, ats in resatomlist.iteritems()]

    return residues


class TopologyObject(object):
    """Base class for all Topology items

    Defines the behaviour by which Bonds/Angles/etc in MDAnalysis should
    behave.

    .. versionadded 0.9.0
    """
    __slots__ = ("atoms")

    def __init__(self, atoms):
        self.atoms = tuple(atoms)

    @property
    def type(self):
        """Type of the bond as a tuple

        .. Note::
             When comparing types, it is important to consider the
             reverse of the type too, ie::
               a.type == b.type or a.type == b.type[::-1]
        """
        return tuple([a.type for a in self.atoms])

    def __repr__(self):
        return "<{cname} between: {conts}>".format(
            cname = self.__class__.__name__,
            conts = ", ".join([
                "Atom {} ({} of {}-{})".format(
                    a.number+1, a.name, a.resname, a.resid)
                for a in self.atoms]))

    def __contains__(self, other):
        """Check whether an atom is in this TopologyObject"""
        return other in self.atoms

    def __eq__(self, other):
        """Check whether two bonds have identical contents"""
        my_tup = tuple([a.number for a in self.atoms])
        ot_tup = tuple([a.number for a in other.atoms])

        return (my_tup == ot_tup) or (my_tup == ot_tup[::-1])

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):  # so bondlists can be sorted
        return self.atoms < other.atoms
        
    def __gt__(self, other):
        return self.atoms > other.atoms

    def __getitem__(self, item):
        """Can retrieve a given Atom from within"""
        return self.atoms[item]

    def __iter__(self):
        return iter(self.atoms)

    def __len__(self):
        return len(self.atoms)


class Bond(TopologyObject):
    """A bond between two :class:`~MDAnalysis.core.AtomGroup.Atom` instances.

    Two :class:`Bond` instances can be compared with the ``==`` and
    ``!=`` operators. A bond is equal to another if the same atom
    numbers are connected and they have the same bond order. The
    ordering of the two atom numbers is ignored as is the fact that a
    bond was guessed.

    The presence of a particular atom can also be queried::
       >>> Atom in Bond
       True / False

    .. versionchanged:: 0.9.0
       * Now a subclass of :class:`TopologyObject`
       * Changed class to use __slots__
       * Changed class to store atoms in .atoms attribute.
    """
    __slots__ = ("atoms", "order", "__is_guessed")

    def __init__(self, atoms, order=None):
        self.atoms = tuple(atoms)
        self.order = order
        self.__is_guessed = False

    def partner(self, atom):
        """Bond.partner(Atom)

        :Returns: the other :class:`~MDAnalysis.core.AtomGroup.Atom` in this
        bond
        """
        if atom is self.atoms[0]:
            return self.atoms[1]
        elif atom is self.atoms[1]:
            return self.atoms[0]
        else:
            raise ValueError("Unrecognised Atom")

    @property
    def is_guessed(self):
        """``True`` if the bond was guessed.

        .. SeeAlso:: :func:`guess_bonds`
        """
        return self.__is_guessed

    @is_guessed.setter
    def is_guessed(self, b):
        self.__is_guessed = b

    def length(self):
        """Length of the bond."""
        bond = self.atoms[0].pos - self.atoms[1].pos
        bond2 = bond * bond
        return sqrt(bond2[0] + bond2[1] + bond2[2])

    def __repr__(self):
        a1, a2 = self.atoms
        s_id = "<Bond between: Atom {0:d} ({1.name} of {1.resname} {1.resid}"\
               " {1.altLoc}) and Atom {2:d} ({3.name} of {3.resname}"\
               "{3.resid} {3.altLoc})".format(
                   a1.number + 1, a1, a2.number + 1, a2)
        try:
            s_length = ", length {0:.2f} A".format(self.length())
        except AttributeError:
            s_length = ""  # no trajectory/coordinates available
        return s_id + s_length + ">"


class Angle(TopologyObject):
    """An angle between three :class:`~MDAnalysis.core.AtomGroup.Atom` instances.
    Atom 2 is the apex of the angle

    .. versionadded:: 0.8
    .. versionchanged:: 0.9.0
       * Now a subclass of :class:`TopologyObject`
       * Changed class to use __slots__
       * Changed class to store atoms in .atoms attribute.
    """
    def angle(self):
        """Returns the angle in degrees of this Angle.

        Angle between atoms 0 and 2 with apex at 1::

              2
             /
            /
           1------0

        .. versionadded:: 0.9.0
        """
        a = self[0].pos - self[1].pos
        b = self[2].pos - self[1].pos
        return numpy.rad2deg(numpy.arccos(numpy.dot(a, b) / (numpy.linalg.norm(a)*numpy.linalg.norm(b))))


class Torsion(TopologyObject):
    """Torsion (dihedral angle) between four
    :class:`~MDAnalysis.core.AtomGroup.Atom` instances.

    The torsion is defined as the angle between the planes formed by
    Atoms (1, 2, 3) and (2, 3, 4).

    .. versionadded:: 0.8
    .. versionchanged:: 0.9.0
       * Now a subclass of :class:`TopologyObject`
       * Changed class to use __slots__
       * Changed class to store atoms in .atoms attribute.
    """
    # http://cbio.bmt.tue.nl/pumma/uploads/Theory/dihedral.png
    def torsion(self):
        """Calculate the dihedral angle in degrees.

        Dihedral angle around axis connecting atoms 1 and 2 (i.e. the angle
        between the planes spanned by atoms (0,1,2) and (1,2,3))::

                  3
                  |
            1-----2
           /
          0

        .. versionadded:: 0.9.0
        """   
        A, B, C, D = self.atoms
        ab = A.position - B.position
        bc = B.position - C.position
        cd = C.position - D.position
        return numpy.rad2deg(dihedral(ab, bc, cd))


class Improper_Torsion(Torsion):  # subclass Torsion to inherit torsion method
    """
    Improper Torsion (improper dihedral angle) between four
    :class:`~MDAnalysis.core.AtomGroup.Atom` instances.

    MDAnalysis treats the improper torsion angle as the angle between 
    the planes formed by Atoms (1, 2, 3) and (2, 3, 4).

    .. warning:: Definitions of Atom ordering in improper torsions
                 can change. Check the definitions here against 
                 your software.

    .. versionadded 0.9.0
    """
    # http://cbio.bmt.tue.nl/pumma/uploads/Theory/improper.png
    def improper(self):
        """Improper dihedral angle in degrees"""
        return self.torsion()


def get_parser_for(filename, permissive=False, tformat=None):
    """Return the appropriate topology parser for *filename*.

    Automatic detection is disabled when an explicit *format* is
    provided.
    """
    from . import _topology_parsers

    tformat = guess_format(filename, format=tformat)
    if tformat == 'PDB' and permissive:
        return _topology_parsers['Permissive_PDB']
    else:
        return _topology_parsers[tformat]


def guess_format(filename, format=None):
    """Returns the type of topology file *filename*.

    The current heuristic simply looks at the filename extension but
    more complicated probes could be implemented here or in the
    individual packages (e.g. as static methods).

    If *format* is supplied then it overrides the auto detection.
    """
    from . import _topology_parsers

    if format is None:
        # simple extension checking... something more complicated is left
        # for the ambitious
        try:
            root, ext = os.path.splitext(filename)
            if ext.startswith('.'):
                ext = ext[1:]
            format = ext.upper()
        except:
            raise TypeError("Cannot determine topology type for %r" % filename)
    else:
        # internally, formats are all uppercase
        format = str(format).upper()

    # sanity check
    if format not in _topology_parsers:
        raise TypeError("Unknown topology format %r for %r; "
                        "only %r are implemented in MDAnalysis." %
                        (format, filename, _topology_parsers.keys()))
    return format


# following guess_* used by PDB parser
def guess_atom_type(atomname):
    """Guess atom type from the name.

    At the moment, this function simply returns the element, as
    guessed by :func:`guess_atom_element`.

    .. SeeAlso:: :func:`guess_atom_element` and
                 :mod:`MDAnalysis.topology.tables`
    """
    return guess_atom_element(atomname)


def guess_atom_element(atomname):
    """Guess the element of the atom from the name.

    Looks in dict to see if element is found, otherwise it uses the first
    character in the atomname. The table comes from CHARMM and AMBER atom
    types, where the first character is not sufficient to determine the atom
    type. Some GROMOS ions have also been added.

    .. Warning: The translation table is incomplete. This will probably result
                in some mistakes, but it still better than nothing!

    .. SeeAlso:: :func:`guess_atom_type` and
                 :mod:`MDAnalysis.topology.tables` (where the data are stored)
    """
    try:
        return tables.atomelements[atomname]
    except KeyError:
        if atomname[0].isdigit():
            # catch 1HH etc
            return atomname[1]
        return atomname[0]


def guess_bonds(atoms, coords, **kwargs):
    """Guess if bonds exist between two atoms based on their distance.

    Bond between two atoms is created, if the two atoms are within

    .. math::

          d < f * (R_1 + R_2)

    of each other, where :math:`R_1` and :math:`R_2` are the VdW radii
    of the atoms and :math:`f` is an ad-hoc *fudge_factor*. This is
    the `same algorithm that VMD uses`_.

    :Keywords:

      *fudge_factor*
        The factor by which atoms must overlap eachother to be considered a 
        bond.  Larger values will increase the number of bonds found. [0.72]

      *vdwradii*
        To supply custom vdwradii for atoms in the algorithm. Must be a dict
        of format {type:radii}. The default table of van der Waals radii is 
        hard-coded as :data:`MDAnalysis.topology.tables.vdwradii`.  Any user
        defined vdwradii passed as an argument will supercede the table
        values. [``None``]

      *lower_bound*
        The minimum bond length. All bonds found shorter than this length will
        be ignored. This is useful for parsing PDB with altloc records where
        atoms with altloc A and B maybe very close together and there should be
        no chemical bond between them. [0.1]

      *box*
        Bonds are found using a distance search, if unit cell information is
        given, periodic boundary conditions will be considered in the distance
        search. [``None``]

    :Returns:
       List of tuples suitable for use in Universe topology building.

    .. warning::

       No check is done after the bonds are guessed to see if Lewis
       structure is correct. This is wrong and will burn somebody.

    .. _`same algorithm that VMD uses`:
       http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.1/ug/node26.html

    .. versionadded:: 0.7.7
    .. versionchanged:: 0.9.0
       * Updated method internally to use more numpy, should work faster.
       * Should also use less memory, previously scaled O(n^2)
       * vdwradii argument now augments table list rather than replacing entirely
    """
    # why not just use atom.positions?
    if len(atoms) != len(coords):
        raise ValueError("'atoms' and 'coord' must be the same length")

    fudge_factor = kwargs.get('fudge_factor', 0.72)

    vdwradii = tables.vdwradii.copy()  # so I don't permanently change it
    user_vdwradii = kwargs.get('vdwradii', None)
    if user_vdwradii:  # this should make algo use their values over defaults
        vdwradii.update(user_vdwradii)

    try:
        atomtypes = set(atoms.types())
    except AttributeError:  # sometimes atoms is just list of atoms not AG
        atomtypes = set([a.type for a in atoms])
    # check that all types have a defined vdw
    if not all([val in vdwradii for val in atomtypes]):
        raise ValueError(("vdw radii for types: " + 
                          ", ".join([t for t in atomtypes if
                                     not t in vdwradii]) + 
                          ". These can be defined manually using the" +
                          " keyword 'vdwradii'"))

    lower_bound = kwargs.get('lower_bound', 0.1)
    
    box = kwargs.get('box', None)

    # to speed up checking, calculate what the largest possible bond
    # atom that would warrant attention.
    # then use this to quickly mask distance results later
    max_vdw = max([vdwradii[t] for t in atomtypes])

    bonds = []

    for i, atom in enumerate(atoms[:-1]):
        vdw_i = vdwradii[atom.type]
        max_d = (vdw_i + max_vdw) * fudge_factor

        # using self_distance_array scales O(n^2)
        # 20,000 atoms = 1.6 Gb memory
        dist = distances.distance_array(coords[i][None, :], coords[i+1:],
                                        box=box)[0]
        idx = numpy.where((dist > lower_bound) & (dist <= max_d))[0]

        for a in idx:
            atom_j = atoms[i + 1 + a]

            if dist[a] < (vdw_i + vdwradii[atom_j.type]) * fudge_factor:
                # because of method used, same bond won't be seen twice, 
                # so don't need to worry about duplicates
                bonds.append((atom.number, atom_j.number))

    return tuple(bonds)


def guess_angles(bonds):
    """Given a list of Bonds, find all angles that exist between atoms.

    Works by assuming that if atoms 1 & 2 are bonded, and 2 & 3 are bonded,
    then (1,2,3) must be an angle.

    :Returns:
      List of tuples defining the angles.  Suitable for use in u._psf

    .. seeAlso:: :meth:`guess_bonds`

    .. versionadded 0.9.0
    """
    angles_found = set()

    for b in bonds:
        for atom in b:
            other_a = b.partner(atom)  # who's my friend currently in Bond
            for other_b in atom.bonds:
                if other_b != b:  # if not the same bond I start as
                    third_a = other_b.partner(atom)
                    desc = tuple([other_a.number, atom.number, third_a.number])
                    if desc[0] > desc[-1]:  # first index always less than last
                        desc = desc[::-1]
                    angles_found.add(desc)

    return tuple(angles_found)


def guess_torsions(angles):
    """Given a list of Angles, find all torsions that exist between atoms.

    Works by assuming that if (1,2,3) is an angle, and 3 & 4 are bonded,
    then (1,2,3,4) must be a torsion.

    :Returns:
      List of tuples defining the torsions.  Suitable for use in u._psf

    .. versionadded 0.9.0
    """
    torsions_found = set()

    for b in angles:
        a_tup = tuple([a.number for a in b])  # angle as tuple of numbers
        # if searching with b[0], want tuple of (b[2], b[1], b[0], +new)
        # search the first and last atom of each angle
        for atom, prefix in zip([b.atoms[0], b.atoms[-1]], 
                                [a_tup[::-1], a_tup]):
            for other_b in atom.bonds:
                if not other_b.partner(atom) in b:
                    third_a = other_b.partner(atom)
                    desc = prefix + (third_a.number,)
                    if desc[0] > desc[-1]:
                        desc = desc[::-1]
                    torsions_found.add(desc)

    return tuple(torsions_found)


def guess_improper_torsions(angles):
    """Given a list of Angles, find all improper torsions that exist between 
    atoms.

    Works by assuming that if (1,2,3) is an angle, and 2 & 4 are bonded,
    then (2, 1, 3, 4) must be an improper torsion.
    ie the improper torsion is the angle between the planes formed by
    (1, 2, 3) and (1, 3, 4)

    :Returns:
      List of tuples defining the improper torsions.  Suitable for use in 
      u._psf

    .. versionadded 0.9.0
    """
    torsions_found = set()

    for b in angles:
        atom = b[1]  # select middle atom in angle
        a_tup = tuple([b[a].number for a in [1, 2, 0]])  # start of improper tuple
        # if searching with b[1], want tuple of (b[1], b[2], b[0], +new)
        # search the first and last atom of each angle
        for other_b in atom.bonds:
            other_atom = other_b.partner(atom)
            if not other_atom in b:  # if this atom isn't in the angle I started with
                desc = a_tup + (other_atom.number,)
                if desc[0] > desc[-1]:
                    desc = desc[::-1]
                torsions_found.add(desc)

    return tuple(torsions_found)


def get_atom_mass(element):
    """Return the atomic mass in u for *element*.

    Masses are looked up in :data:`MDAnalysis.topology.tables.masses`.

    .. Warning:: Unknown masses are set to 0.00
    """
    try:
        return tables.masses[element]
    except KeyError:
        return 0.000


def guess_atom_mass(atomname):
    """Guess a mass based on the atom name.

    :func:`guess_atom_element` is used to determine the kind of atom.

    .. warning:: Anything not recognized is simply set to 0; if you rely on the
                 masses you might want to double check.
    """
    return get_atom_mass(guess_atom_element(atomname))


def guess_atom_charge(atomname):
    """Guess atom charge from the name.

    .. Warning:: Not implemented; simply returns 0.
    """
    # TODO: do something slightly smarter, at least use name/element
    return 0.0


class TopologyDict(object):
    """A customised dictionary designed for sorting the bonds, angles and
    torsions present in a group of atoms.

    Usage::

      topologydict = TopologyDict(members)

    :Arguments:
        *members*
            A list of :class:`TopologyObject` instances

    :Returns:
        *topologydict*
            A specialised dictionary of the topology instances passed to it

    TopologyDicts are also built lazily from a :class:`TopologyGroup.topDict`
    attribute.

    The TopologyDict collects all the selected topology type from the
    atoms and categorises them according to the types of the atoms within.
    A :class:`TopologyGroup` containing all of a given bond type can
    be made by querying with the appropriate key.  The keys to the
    topologyDict are a tuple of the atom types that the bond represents
    and can be viewed using the keys() method.

    For example, from a system containing pure ethanol ::

      >>> td = u.bonds.topDict
      >>> td.keys()
      [('C', 'C'),
       ('C', 'H'),
       ('O', 'H'),
       ('C', 'O')]
      >>> td['C', 'O']
      < TopologyGroup containing 912 bonds >

    .. Note::

       The key for a bond is taken from the type attribute of the atoms.

       Getting and setting types of bonds is done smartly, so a C-C-H
       angle is considered identical to a H-C-C angle.

    Duplicate entries are automatically removed upon creation and
    combination of different Dicts.  This means a bond between atoms
    1 and 2 will only ever appear once in a dict despite both atoms 1
    and 2 having the bond in their .bond attribute.

    Two TopologyDicts can be combined using addition, this will not
    create any duplicate bonds in the process.

    .. versionadded:: 0.8
    .. versionchanged:: 0.9.0
       * Changed initialisation to use a list of TopologyObjects
       instead of list of atoms
       * Now used from within TopologyGroup instead of accessed from AtomGroup
    """
    def __init__(self, members):
        self.dict = dict()
        # Detect what I've been given
        if isinstance(members[0], TopologyObject):
            self.toptype = members[0].__class__.__name__
        else:  # Throw error if not given right thing
            raise TypeError('Unrecognised input')

        for b in members:
            btype = b.type
            try:
                self.dict[btype] += [b]
            except KeyError:
                self.dict[btype] = [b]

        self._removeDupes()

    def _removeDupes(self):
        """Sorts through contents and makes sure that there are no duplicate keys
        (through type reversal)
        """
        newdict = dict()

        # Go through all keys, if the reverse of the key exists add this to
        # that entry else make a new entry
        for k in self.dict:
            if not k[::-1] in newdict:
                newdict[k] = self.dict[k]
            else:
                newdict[k[::-1]] += self.dict[k]

        self.dict = newdict

    def __len__(self):
        """Returns the number of types of bond in the topology dictionary"""
        return len(self.dict.keys())

    def keys(self):
        """Returns a list of the different types of available bonds"""
        return self.dict.keys()

    def __iter__(self):
        """Iterator over keys in this dictionary"""
        return iter(self.dict)

    def __repr__(self):
        return "<TopologyDict with {num} unique {type}s>".format(
            num=len(self), type=self.toptype)

    def __getitem__(self, key):
        """Returns a TopologyGroup matching the criteria if possible,
        otherwise returns ``None``
        """
        if key in self:
            if key in self.dict:
                selection = self.dict[key]
            else:
                selection = self.dict[key[::-1]]

            return TopologyGroup(selection)
        else:
            raise KeyError(key)

    def __contains__(self, other):
        """
        Returns boolean on whether a given type exists within this dictionary

        For topology groups the key (1,2,3) is considered the same as (3,2,1)
        """
        return other in self.dict or other[::-1] in self.dict


class TopologyGroup(object):
    """A container for a groups of bonds.

    All bonds of a certain types can be retrieved from within the 
    :class:`TopologyGroup` by querying with a tuple of types::

      tg2 = tg.selectBonds([key])

    Where *key* describes the desired bond as a tuple of the involved 
    :class:`~MDAnalysis.AtomGroup.Atom` types, as defined by the .type Atom 
    attribute). A list of available keys can be displayed using the 
    :meth:`types` method.

    Alternatively, all the bonds which are in a given 
    :class:`~MDAnalysis.AtomGroup.AtomGroup` can be extracted using 
    :meth:`atomgroup_intersection`::

      tg2 = tg.atomgroup_intersection(ag)

    This allows the keyword *strict* to be given, which forces all members of
    all bonds to be inside the AtomGroup passed to it.

    Finally, a TopologyGroup can be sliced similarly to AtomGroups::

      tg2 = tg[5:10]

    The :meth:`bonds`, :meth:`angles` and :meth:`torsions` methods offer
    a "shortcut" to the Cython distance calculation functions in
    :class:`MDAnalysis.core.distances`.

    TopologyGroups can be combined with TopologyGroups of the same bond
    type (ie can combine two angle containing TopologyGroups).

    .. versionadded:: 0.8
    .. versionchanged:: 0.9.0
       Overhauled completely
       * Added internal TopologyDict accessible by the topDict attribute
       * :meth:`selectBonds` allows the topDict to be queried with tuple of 
         types
       * Added :meth:`atomgroup_intersection` to allow bonds which are in 
         a given AtomGroup to be retrieved.
    """
    def __init__(self, bondlist):
        if len(bondlist) == 0:
            raise ValueError("Can't make empty TopologyGroup")
        if isinstance(bondlist[0], TopologyObject):
            self.toptype = bondlist[0].__class__.__name__
        else:
            raise TypeError("Input not recognised")
        # Would be nice to make everything work internally using sets, BUT
        # sets can't be indexed, so couldn't work backward from .angles()
        # results to find which angle is a certain value.
        # Sorted so that slicing returns sensible results
        self.bondlist = tuple(sorted(set(bondlist)))

        self._cache = dict()  # used for topdict saving

    def selectBonds(self, selection):
        """Retrieves a selection from this topology group based on types.

        .. seeAlso :meth:`types`

        .. versionadded 0.9.0
        """
        return self.topDict[selection]

    def types(self):
        """Return a list of the bond types in this TopologyGroup

        .. versionadded 0.9.0
        """
        return self.topDict.keys()

    @property
    @cached('dict')
    def topDict(self):
        """
        Returns the TopologyDict for this topology group.

        This is used for the selectBonds method when fetching a certain type
        of bond.

        This is a cached property so will be generated the first time it is
        accessed.

        .. versionadded 0.9.0
        """
        return TopologyDict(self.bondlist)

    def atomgroup_intersection(self, ag, **kwargs):
        """Retrieve all bonds from within this TopologyGroup that are within
        the AtomGroup which is passed.

        :Keywords:
          *strict*
            Only retrieve bonds which are completely contained within the
            AtomGroup. [``False``]
        
        .. versionadded:: 0.9.0
        """
        strict = kwargs.get('strict', False)
        if strict:
            return self._strict_intersection(ag)
        else:
            return self._loose_intersection(ag)

    def _loose_intersection(self, other):
        """Copies bonds if it appears even once in an AtomGroup

        This means that some bonds might extend out of the defined AtomGroup

        .. SeeAlso:: :meth:`_strict_intersection` for including bonds
                     more strictly

        .. versionadded 0.9.0
        """
        if self.toptype == 'Bond':
            other_set = set(other.bonds)
        elif self.toptype == 'Angle':
            other_set = set(other.angles)
        elif self.toptype == 'Torsion':
            other_set = set(other.torsions)
        elif self.toptype == 'Improper_Torsion':
            other_set = set(other.impropers)
        else:
            raise ValueError("Unsupported intersection")

        newlist = list(set(self.bondlist).intersection(other_set))

        if len(newlist) > 0:
            return TopologyGroup(newlist)
        else:
            return None

    def _strict_intersection(self, other):
        """Copies bonds only if all members of the bond appear in the AtomGroup

        This means that all bonds will be contained within the AtomGroup

        .. SeeAlso:: :meth:`_loose_intersection` for including bonds
                     less strictly

        .. versionadded 0.9.0
        """
        # Create a dictionary of all bonds within this TG, initialised to 0
        # for each
        # 
        # Then go through all TopObjs in AtomGroup and count their appearances
        # keeping track using the dict
        # 
        # Then see how many times each TopObj was spotted in the AtomGroup's bonds
        #
        # If this count is equal to crit, (bond=2, angle=3, torsion=4) then
        # the TopObj was seen enough for it to have to be completely be 
        # present in the AtomGroup

        # each bond starts with 0 appearances
        # I'm only interested in intersection, so if its not in tg then
        # i'll get keyerrors which i'll pass
        count_dict = {b: 0 for b in self.bondlist}

        # then go through ag and count appearances of bonds
# This seems to benchmark slow, because __getattribute__ is slower than a.bonds
#        for atom in other:
#            for b in atom.__getattribute__(req_attr):
#                try:
#                    count_dict[b] += 1
#                except KeyError:  # if he's not in dict then meh
#                    pass
# So I'll bruteforce here, despite it being fugly
        if self.toptype == 'Bond':
            crit = 2
            for atom in other:
                for b in atom.bonds:
                    try:
                        count_dict[b] += 1
                    except KeyError:
                        pass
        elif self.toptype == 'Angle':
            crit = 3
            for atom in other:
                for b in atom.angles:
                    try:
                        count_dict[b] += 1
                    except KeyError:
                        pass
        elif self.toptype == 'Torsion':
            crit = 4
            for atom in other:
                for b in atom.torsions:
                    try:
                        count_dict[b] += 1
                    except KeyError:
                        pass
        elif self.toptype == 'Improper_Torsion':
            crit = 4
            for atom in other:
                for b in atom.impropers:
                    try:
                        count_dict[b] += 1
                    except KeyError:
                        pass

        # now make new list, which only includes bonds with enough appearances
        newlist = [b for b in self.bondlist if count_dict[b] == crit]

        if len(newlist) > 0:
            return TopologyGroup(newlist)
        else:
            return None

    def dump_contents(self):
        """Return a tuple of tuples which define the contents of this
        TopologyGroup in terms of the atom numbers,
        (0 based index within u.atoms)

        This format should be identical to the original contents of the
        entries in universe._psf.
        Note that because bonds are sorted as they are initialised, the order
        that atoms are defined in each entry might be reversed.

        .. versionadded 0.9.0
        """
        # should allow topology information to be pickled even if it is
        # substantially changed from original input, 
        # eg through merging universes or defining new bonds manually.
        bondlist = tuple([tuple([a.number for a in b.atoms]) 
                          for b in self.bondlist])

        return bondlist

    @property
    @cached('atom1')
    def atom1(self):
        return AtomGroup.AtomGroup([b[0] for b in self.bondlist])

    @property
    @cached('atom2')
    def atom2(self):
        return AtomGroup.AtomGroup([b[1] for b in self.bondlist])

    @property
    @cached('atom3')
    def atom3(self):
        return AtomGroup.AtomGroup([b[2] for b in self.bondlist])

    @property
    @cached('atom4')
    def atom4(self):
        return AtomGroup.AtomGroup([b[3] for b in self.bondlist])

    def __len__(self):
        """Number of bonds in the topology group"""
        return len(self.bondlist)

    def __add__(self, other):
        """Combine two TopologyGroups together.

        Can combined two TopologyGroup of the same type, or add a single
        TopologyObject to a TopologyGroup.
        """
        if isinstance(other, TopologyObject):
            if type(other) != type(self.bondlist[0]):
                raise TypeError("Cannot add different types of "
                                "TopologyObjects together")
            else:
                return TopologyGroup(self.bondlist + (other,))
        if not isinstance(other, TopologyGroup):
            raise TypeError("Can only combine two TopologyGroups")
        elif self.toptype != other.toptype:
            raise TypeError("Can only combine TopologyGroups of the same type")

        return TopologyGroup(self.bondlist + other.bondlist)

    def __getitem__(self, item):
        """Returns a particular bond as single object or a subset of
        this TopologyGroup as another TopologyGroup
        """
        if numpy.dtype(type(item)) == numpy.dtype(int):
            return self.bondlist[item]  # single TopObj
        elif type(item) == slice:
            return TopologyGroup(self.bondlist[item])  # new TG
        elif isinstance(item, (numpy.ndarray, list)):
            return TopologyGroup([self.bondlist[i] for i in item])

    def __iter__(self):
        """Iterator over all bonds"""
        return iter(self.bondlist)

    def __contains__(self, item):
        """Tests if this TopologyGroup contains a bond"""
        return item in set(self.bondlist)

    def __repr__(self):
        return "<TopologyGroup containing {num} {type}s>".format(
            num=len(self), type=self.toptype)

    def __eq__(self, other):
        """Test if contents of TopologyGroups are equal"""
        return set(self) == set(other)

    def __ne__(self, other):
        return not self.__eq__(other)

    # Distance calculation methods below
    # "Slow" versions exist as a way of testing the Cython implementations
    def _bondsSlow(self, pbc=False):  # pragma: no cover
        """Slow version of bond (numpy implementation)"""
        if not self.toptype == 'Bond':
            return TypeError("TopologyGroup is not of type 'Bond'")
        else:
            bond_dist = self.atom1.coordinates() - self.atom2.coordinates()
            if pbc:
                box = self.atom1.dimensions
                # orthogonal and divide by zero check
                if (box[6:9] == 90.).all() and not (box[0:3] == 0).any():
                    bond_dist -= numpy.rint(bond_dist / box[0:3]) * box[0:3]
                else:
                    raise ValueError("Only orthogonal boxes supported")

            return numpy.array([norm(a) for a in bond_dist])

    def bonds(self, pbc=False, result=None):
        """Calculates the distance between all bonds in this TopologyGroup

        :Keywords:
           *pbc*
              apply periodic boundary conditions when calculating distance
              [False]
           *result*
              allows a predefined results array to be used,
              note that this will be overwritten

        Uses cython implementation
        """
        if not self.toptype == 'Bond':
            raise TypeError("TopologyGroup is not of type 'Bond'")
        if not result:
            result = numpy.zeros(len(self), numpy.float64)
        if pbc:
            return distances.calc_bonds(self.atom1.coordinates(),
                                        self.atom2.coordinates(),
                                        box=self.atom1.dimensions,
                                        result=result)
        else:
            return distances.calc_bonds(self.atom1.coordinates(),
                                        self.atom2.coordinates(),
                                        result=result)

    def _anglesSlow(self):  # pragma: no cover
        """Slow version of angle (numpy implementation)"""
        if not self.toptype == 'Angle':
            raise TypeError("TopologyGroup is not of type 'Angle'")

        vec1 = self.atom1.coordinates() - self.atom2.coordinates()
        vec2 = self.atom3.coordinates() - self.atom2.coordinates()

        angles = numpy.array([slowang(a, b) for a, b in izip(vec1, vec2)])
        return angles

    def angles(self, result=None, pbc=False):
        """Calculates the angle in radians formed between a bond
        between atoms 1 and 2 and a bond between atoms 2 & 3

        :Keywords:
           *result*
              allows a predefined results array to be used, note that this
              will be overwritten
           *pbc*
              apply periodic boundary conditions when calculating angles
              [False] this is important when connecting vectors between atoms
              might require minimum image convention

        Uses cython implementation

        .. versionchanged :: 0.9.0
           Added pbc option (default False)
        """
        if not self.toptype == 'Angle':
            raise TypeError("TopologyGroup is not of type 'Angle'")
        if not result:
            result = numpy.zeros(len(self), numpy.float64)
        if pbc:
            return distances.calc_angles(self.atom1.coordinates(),
                                         self.atom2.coordinates(),
                                         self.atom3.coordinates(),
                                         box=self.atom1.dimensions,
                                         result=result)
        else:
            return distances.calc_angles(self.atom1.coordinates(),
                                         self.atom2.coordinates(),
                                         self.atom3.coordinates(),
                                         result=result)

    def _torsionsSlow(self):  # pragma: no cover
        """Slow version of torsion (numpy implementation)"""
        if self.toptype not in ['Torsion', 'Improper_Torsion']:
            raise TypeError("TopologyGroup is not of type 'Torsion' or "
                            "'Improper_Torsion'")

        vec1 = self.atom2.coordinates() - self.atom1.coordinates()
        vec2 = self.atom3.coordinates() - self.atom2.coordinates()
        vec3 = self.atom4.coordinates() - self.atom3.coordinates()

        return numpy.array([dihedral(a, b, c)
                            for a, b, c in izip(vec1, vec2, vec3)])

    def torsions(self, result=None, pbc=False):
        """Calculate the torsional angle in radians for this topology
        group.

        Defined as the angle between a plane formed by atoms 1, 2 and
        3 and a plane formed by atoms 2, 3 and 4.

        :Keywords:
           *result*
              allows a predefined results array to be used, note that this
              will be overwritten
           *pbc*
              apply periodic boundary conditions when calculating angles
              [False] this is important when connecting vectors between
              atoms might require minimum image convention

        Uses cython implementation.

        .. versionchanged:: 0.9.0
           Added pbc option (default False)
        """
        if self.toptype not in ['Torsion', 'Improper_Torsion']:
            raise TypeError("TopologyGroup is not of type 'Torsion' or "
                            "'Improper_Torsion'")
        if not result:
            result = numpy.zeros(len(self), numpy.float64)
        if pbc:
            return distances.calc_torsions(self.atom1.coordinates(),
                                           self.atom2.coordinates(),
                                           self.atom3.coordinates(),
                                           self.atom4.coordinates(),
                                           box=self.atom1.dimensions,
                                           result=result)
        else:
            return distances.calc_torsions(self.atom1.coordinates(),
                                           self.atom2.coordinates(),
                                           self.atom3.coordinates(),
                                           self.atom4.coordinates(),
                                           result=result)
