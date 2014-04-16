# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2011 Naveen Michaud-Agrawal,
#               Elizabeth J. Denning, Oliver Beckstein,
#               and contributors (see website for details)
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
#     N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
#     O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
#     Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
#     doi:10.1002/jcc.21787
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
import itertools
import os.path
from math import sqrt
import numpy
import sys

# Local imports
from . import tables
from ..core import distances
from ..core.util import norm
from ..core import AtomGroup


def build_segments(atoms):
    """Create all :class:`~MDAnalysis.core.AtomGroup.Segment` instancess from a list 
    of :class:`~MDAnalysis.core.AtomGroup.Atom` instances.

    The function also builds the :class:`~MDAnalysis.core.AtomGroup.Residue`
    instances by tracking residue numbers.

    Updating segments also changes the underlying
    :class:`~MDAnalysis.core.AtomGroup.Atom` instances, which record
    to which residue and segment an atom belongs.

    :Returns: structure dict, which associates a segname with a
              :class:`~MDAnalysis.core.AtomGroup.Segment`

    """
    struc = {}
    residues = []
    resatomlist = []
    curr_segname = atoms[0].segid
    curr_resnum = atoms[0].resid
    curr_resname = atoms[0].resname
    for a in atoms:
        if a.segid == curr_segname:
            if a.resid == curr_resnum:
                resatomlist.append(a)
            else:
                # New residue
                residues.append(AtomGroup.Residue(curr_resname, curr_resnum, resatomlist))
                resatomlist = [a]
                curr_resnum = a.resid
                curr_resname = a.resname
        else:
            # We've come to a new segment
            residues.append(AtomGroup.Residue(curr_resname, curr_resnum, resatomlist))
            struc[curr_segname] = AtomGroup.Segment(curr_segname, residues)
            residues = []
            resatomlist = [a]
            curr_resnum = a.resid
            curr_resname = a.resname
            curr_segname = a.segid
    # Add the last segment
    residues.append(AtomGroup.Residue(curr_resname, curr_resnum, resatomlist))
    struc[curr_segname] = AtomGroup.Segment(curr_segname, residues)
    return struc


def build_residues(atoms):
    """Create a list :class:`~MDAnalysis.core.AtomGroup.Residue` instances from a list 
    of :class:`~MDAnalysis.core.AtomGroup.Atom` instances.

    Updating residues also changes the underlying
    :class:`~MDAnalysis.core.AtomGroup.Atom` instances, which record
    to which residue an atom belongs.

    :Returns: List of :class:`~MDAnalysis.core.AtomGroup.Residue` instances

    .. versionadded:: 0.8
    """
    struc = {}
    residues = []
    resatomlist = []
    curr_resnum = atoms[0].resid
    curr_resname = atoms[0].resname
    for a in atoms:
        if a.resid == curr_resnum:
            resatomlist.append(a)
        else:
            # New residue
            residues.append(AtomGroup.Residue(curr_resname, curr_resnum, resatomlist))
            resatomlist = [a]
            curr_resnum = a.resid
            curr_resname = a.resname
    # Add the last residue
    residues.append(AtomGroup.Residue(curr_resname, curr_resnum, resatomlist))
    return residues


class Bond(object):
    """A bond between two :class:`~MDAnalysis.core.AtomGroup.Atom` instances.

    Two :class:`Bond` instances can be compared with the ``==`` and
    ``!=`` operators. A bond is equal to another if the same atom
    numbers are connected and they have the same bond order. The
    ordering of the two atom numbers is ignored as is the fact that a
    bond was guessed.

    The presence of a particular atom can also be queried::
       >>> Atom in Bond
       True / False
    """

    def __init__(self, a1, a2, order=None):
        self.atom1 = a1
        self.atom2 = a2
        a1.bonds.append(self)
        a2.bonds.append(self)
        self.order = order
        self.__is_guessed = False

    def partner(self, atom):
        """Bond.partner(Atom)

        :Returns: the other :class:`~MDAnalysis.core.AtomGroup.Atom` in this bond
        """
        if atom is self.atom1:
            return self.atom2
        else:
            return self.atom1

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
        bond = self.atom1.pos - self.atom2.pos
        bond2 = bond * bond
        return sqrt(bond2[0] + bond2[1] + bond2[2])

    def __repr__(self):
        a1, a2 = self.atom1, self.atom2
        s_id = "< Bond between: Atom {0:d} ({1.name} of {1.resname} {1.resid} {1.altLoc}) and " \
               "Atom {2:d} ({3.name} of {3.resname} {3.resid} {3.altLoc})".format(a1.number + 1, a1, a2.number + 1, a2)
        try:
            s_length = ", length {0:.2f} A".format(self.length())
        except AttributeError:
            s_length = ""  # no trajectory/coordinates available
        return s_id + s_length + ">"

    def __contains__(self, other):
        """Can check if a particular atom is present in this bond"""
        return other == self.atom1 or other == self.atom2

    def __eq__(self, other):
        """This bond is equal to *other* if the same atom numbers are connected.

        The bond order must also be the same. Only two :class:`Bond`
        instances can be compared with each other.

        The ordering of the two atom numbers is ignored as is the fact
        that a bond was guessed.
        """
        if type(other) is type(self):
            return self.order == other.order and \
                   set((self.atom1.number, self.atom2.number)) == set((other.atom1.number, other.atom2.number))
        else:
            return False

    def __ne__(self, other):
        """This bond is not equal to *other* if different atom  numbers are connected.

        .. SeeAlso:: This is the logical opposite of :meth:`Bond.__eq__`.
        """
        return not self.__eq__(other)


class Angle(object):
    """An angle between three :class:`~MDAnalysis.core.AtomGroup.Atom` instances.
    Atom 2 is the apex of the angle

    .. versionadded:: 0.8
    """
    #currently a stub of a class, maybe add to this to make it as useful as bonds
    def __init__(self, a1, a2, a3):
        self.atom1 = a1
        self.atom2 = a2 # middle atom in angle
        self.atom3 = a3

    def __repr__(self):
        a1 = self.atom1
        a2 = self.atom2
        a3 = self.atom3
        return "< Angle between: Atom %d (%s of %s-%d), Atom %d (%s of %s-%d) and Atom %d (%s of %s-%d) >" % \
               (a1.number + 1, a1.name, a1.resname, a1.resid,
                a2.number + 1, a2.name, a2.resname, a2.resid,
                a3.number + 1, a3.name, a3.resname, a3.resid)

    def __contains__(self, other):
        return other == self.atom1 or other == self.atom2 or other == self.atom3

class Torsion(object):
    """Torsion (dihedral angle) between four :class:`~MDAnalysis.core.AtomGroup.Atom` instances.

    The torsion is defined as the angle between the planes formed by (1, 2, 3) and (2, 3, 4).

    .. versionadded:: 0.8
    """

    def __init__(self, a1, a2, a3, a4):
        self.atom1 = a1
        self.atom2 = a2
        self.atom3 = a3
        self.atom4 = a4

    def __repr__(self):
        a1 = self.atom1
        a2 = self.atom2
        a3 = self.atom3
        a4 = self.atom4
        return "< Torsion between Atom %d (%s of %s-%d), Atom %d (%s of %s-%d), Atom %d (%s of %s-%d) and Atom %d (%s of %s-%d) >" % \
               (a1.number + 1, a1.name, a1.resname, a1.resid,
                a2.number + 1, a2.name, a2.resname, a2.resid,
                a3.number + 1, a3.name, a3.resname, a3.resid,
                a4.number + 1, a4.name, a4.resname, a4.resid)

    def __contains__(self, other):
        return other == self.atom1 or other == self.atom2 or \
            other == self.atom3 or other == self.atom4


def build_bondlists(atoms, bonds=None, angles=None, torsions=None):
    """Construct the topology lists of each :class:`~MDAnalysis.core.AtomGroup.Atom`.

    The lists are stored in the attributes
    :attr:`MDAnalysis.core.AtomGroup.Atom.bonds`
    :attr:`MDAnalysis.core.AtomGroup.Atom.angles`
    :attr:`MDAnalysis.core.AtomGroup.Atom.torsions`
    and consist of a list of
    :class:`Bond` :class:`Angle` and :class:`Torsion` instances respectively

    .. versionchanged:: 0.8
    """
    if bonds:
        for a in atoms:
            a.bonds = []
        for a1, a2 in bonds:
            atom1 = atoms[a1]
            atom2 = atoms[a2]
            b = Bond(atom1, atom2)
            atom1.bonds.append(b)
            atom2.bonds.append(b)
    if angles:
        for a in atoms:
            a.angles = []
        for a1, a2, a3 in angles:
            atom1 = atoms[a1]
            atom2 = atoms[a2]
            atom3 = atoms[a3]
            b = Angle(atom1, atom2, atom3)
            atom1.angles.append(b)
            atom2.angles.append(b)
            atom3.angles.append(b)
    if torsions:
        for a in atoms:
            a.torsions = []
        for a1, a2, a3, a4 in torsions:
            atom1 = atoms[a1]
            atom2 = atoms[a2]
            atom3 = atoms[a3]
            atom4 = atoms[a4]
            b = Torsion(atom1, atom2, atom3, atom4)
            atom1.torsions.append(b)
            atom2.torsions.append(b)
            atom3.torsions.append(b)
            atom4.torsions.append(b)


def get_parser_for(filename, permissive=False, bonds=False, format=None):
    """Return the appropriate topology parser for *filename*.

    Automatic detection is disabled when an explicit *format* is
    provided.
    """
    from . import _topology_parsers_permissive, _topology_parsers, _topology_parsers_bonds
    format = guess_format(filename, format=format)
    if permissive and not bonds:
        return _topology_parsers_permissive[format]
    if permissive and bonds:
        return _topology_parsers_bonds[format]
    return _topology_parsers[format]


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
        root, ext = os.path.splitext(filename)
        try:
            if ext.startswith('.'):
                ext = ext[1:]
            format = ext.upper()
        except:
            raise TypeError("Cannot determine topology type for %r" % filename)
    else:
        # internally, formats are all uppercase
        format = str(format).upper()

    # sanity check
    if not format in _topology_parsers:
        raise TypeError("Unknown topology format %r for %r; only %r are implemented in MDAnalysis." %
                        (format, filename, _topology_parsers.keys()))
    return format

# following guess_* used by PDB parser

def guess_atom_type(atomname):
    """Guess atom type from the name.

    At the moment, this function simply returns the element, as
    guessed by :func:`guess_atom_element`.

    .. SeeAlso:: :func:`guess_atom_element` and :mod:`MDAnalysis.topology.tables`
    """
    return guess_atom_element(atomname)


def guess_atom_element(atomname):
    """Guess the element of the atom from the name.

    Looks in dict to see if element is found, otherwise it uses the first character in the atomname.
    The table comes from CHARMM and AMBER atom types, where the first character is not sufficient to
    determine the atom type. Some GROMOS ions have also been added.

    .. Warning: The translation table is incomplete. This will probably result
                in some mistakes, but it still better than nothing!

    .. SeeAlso:: :func:`guess_atom_type` and
                 :mod:`MDAnalysis.topology.tables` (where the data are stored)
    """
    try:
        return tables.atomelements[atomname]
    except KeyError:
        if atomname.startswith(('1', '2', '3', '4', '5', '6', '7', '8', '9')):
            # catch 1HH etc
            return atomname[1]
        return atomname[0]


def guess_bonds(atoms, coords, fudge_factor=0.72, vdwradii=None, lower_bound=0.1):
    """Guess if bonds exist between two atoms based on their distance.

    Bond between two atoms is created, if the two atoms are within

    .. math::

          d < f * (R_1 + R_2)

    of each other, where :math:`R_1` and :math:`R_2` are the VdW radii
    of the atoms and :math:`f` is an ad-hoc *fudge_factor*. This is
    the `same algorithm that VMD uses`_.

    The table of van der Waals radii is hard-coded as
    :data:`MDAnalysis.topology.tables.vdwradii` or a user-specified
    table can be provided as a dictionary in the keyword argument
    *vdwradii*. Atoms are found by their :attr:`Atom.type`.
    
    *lower_bound* defines a heuristic cutoff below which a bond is too short to 
    exist. This is useful for parsing PDB with altloc records where atoms with 
    altloc A and B maybe very close together and there should be no chemical 
    bond between them. 

    .. warning::

       No check is done after the bonds are guessed to see if Lewis
       structure is correct. This is wrong and will burn somebody.

    The code is also in pure python now, so it's slow.

    .. _`same algorithm that VMD uses`:
       http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.1/ug/node26.html

    .. versionadded:: 0.7.7
    """
    # Taken from GROMACS gromacs/top/vdwradii.dat
    # FIXME by JD: these should be stored in an asset file, rather than in
    # source code.
    # FIXME this is not the whole periodic table... (eg halogens are missing)
    if vdwradii is None:
        vdwradii = tables.vdwradii

    assert len([a for a in atoms if a]) == coords.shape[0]

    bonds = set()
    for a in atoms:
        if not a.type in vdwradii.keys():
            # SB: logger is not defined. printting to stderr instead
            print("ERROR: guess_bonds(): %s has no defined vdw radius, cannot guess bonds" % a.type, file=sys.stderr)
            return bonds
    # 1-D vector of the upper-triangle of all-to-all distance matrix
    dist = distances.self_distance_array(coords)
    N = len(coords)

    pairs = ((i, j) for i in xrange(N) for j in xrange(i + 1, N))
    for x, (d, (i, j)) in enumerate(itertools.izip(dist, pairs)):
        a1, a2 = atoms[i], atoms[j]
        r1, r2 = vdwradii[a1.type], vdwradii[a2.type]
        
        # ideal bond distance smaller or equal distance
        if (r1 + r2) * fudge_factor <= dist[x]:
            continue
          
        # filter out unusually short bonds - like the ones detected between identical atoms with different altloc records
        if dist[x] <= lower_bound: continue
        #print "BOND", ((r1 + r2) * 10 * fudge_factor), dist[i,j]
        bonds.add(frozenset([a1.number, a2.number]))  # orbeckst: Atom.number are 0-based, do not add 1.

    return bonds


def get_atom_mass(element):
    """Return the atomic mass in u for *element*.

    Masses are looked up in :data:`MDAnalysis.topology.tables.masses`.

    .. Warning:: Unknown masses are set to 0.
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
    """A customised dictionary designed for sorting the bonds, angles and torsions
    present in a group of atoms.

    Usage::

      topologydict = TopologyDict(topologytype, atoms)

    :Arguments:
        *topologytype* 
            one of 'bond' 'angle' or 'torsion'; a single
            TopologyDict can only handle one type of topology.

        *atoms* 
            a list of :class:`MDAnalysis.core.AtomGroup.Atom` objects.

    :Returns:
        *topologydict*
            A dictionary of the selected topology type
    
    TopologyDicts are also built lazily from a :class:`MDAnalysis.core.AtomGroup.AtomGroup`
    using the :meth:`MDAnalysis.core.AtomGroup.bondDict` :meth:`MDAnalysis.core.AtomGroup.angleDict`
    or :meth:`MDAnalysis.core.AtomGroup.torsionDict` methods.

    The TopologyDict collects all the selected topology type from the
    atoms and categorises them according to the types of the atoms within.
    A :class:`TopologyGroup` containing all of a given bond type can
    be made by querying with the appropriate key.  The keys to the
    topologyDict are a tuple of the atom types that the bond represents
    and can be viewed using the keys() method.

    For example, from a system containing pure ethanol ::

      >>> td = u.atoms.bondDict
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
    """

    def __init__(self, toptype, members):
        self.dict = dict()
        self.toptype = toptype
        if isinstance(members[0], dict): #combining two topologyDicts
            for d in members:
                for k in d.keys():
                    if not k in self.dict:
                        self.dict[k] = d[k]
                    else:
                        self.dict[k] += d[k]
        elif toptype == 'bond':
            for atom in members: #loop through all atoms
                for b in atom.bonds: #loop through all bonds this atom has
                    btype = (b.atom1.type, b.atom2.type)
                    try:
                        self.dict[btype] += [b]
                    except KeyError:
                        self.dict[btype] = [b]
        elif toptype == 'angle':
            for atom in members:
                for b in atom.angles:
                    btype = (b.atom1.type, b.atom2.type, b.atom3.type)
                    try:
                        self.dict[btype] += [b]
                    except KeyError:
                        self.dict[btype] = [b]
        elif toptype == 'torsion':
            for atom in members:
                for b in atom.torsions:
                    btype = (b.atom1.type, b.atom2.type, b.atom3.type, b.atom4.type)
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
        #First remove duplicate keys
        for k in self.dict:
            if not tuple(reversed(k)) in newdict:#newdict starts blank, so having k already is impossible
                newdict[k] = []
            else:
                self.dict[k] += self.dict[tuple(reversed(k))] #if the reverse exists already, pile those values onto the reverse in old dict
        #newdict now has unique set of keys but no values (yet!)
        for k in newdict: #loop over only newdict's keys as all values are in these key values in the old dict
            for v in self.dict[k]:
                if not v in newdict[k]:
                    newdict[k] += [v]
        self.dict = newdict

    def __add__(self, other):
        if not isinstance(other, TopologyDict):
            raise TypeError("Can only combine topologyDicts")
        if not self.toptype == other.toptype:
            raise TypeError("topologyDicts are of different type, cannot combine")
        return TopologyDict(self.toptype, [self.dict, other.dict])

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
        return '<' + self.__class__.__name__ + ' with ' + repr(len(self)) + ' unique ' + self.toptype + 's >'

    def __getitem__(self, key):
        """Returns a TopologyGroup matching the criteria if possible, otherwise returns None"""
        if key in self:
            if key in self.dict:
                selection = self.dict[key]
            else:
                selection = self.dict[tuple(reversed(key))]

            return TopologyGroup(selection)
        else:
            raise KeyError(key)

    def __contains__(self, other):
        """Returns boolean on whether a topology group exists within this dictionary"""
        # For topology groups, 1-2-3 is considered the same as 3-2-1
        return other in self.dict or tuple(reversed(other)) in self.dict


class TopologyGroup(object):
    """ A container for a group of bonds (either bonds, angles or torsions)::

      tg = atomGroup.selectBonds(key)
      tg = TopologyDict[key]

    *key* describes the desired bond as a tuple of the involved atom
    types (as defined by the .type atom attribute). A list of available
    topology keys can be displayed using the .keys() method.

    The TopologyGroup contains :class:`MDAnalysis.core.AtomGroup.AtomGroup`
    instances which correspond to the components of the bonds the
    TopologyGroup contains. Ie a bond has 2 AtomGroups whereas an angle
    has 3.

    The :meth:`bonds`, :meth:`angles` and :meth:`torsions` methods offer
    a "shortcut" to the Cython distance calculation functions in
    :class:`MDAnalysis.core.distances`.

    TopologyGroups can be combined with TopologyGroups of the same bond
    type (ie can combine two angle containing TopologyGroups).

    TopologyGroups can be indexed to return a single :class:`Bond`
    :class:`Angle` or :class:`Torsion` ::

      tg[0], tg[-2]

    Or sliced to return a TopologyGroup containing a subset of the original::

      tg[4:-4]

    .. versionadded:: 0.8
    """

    def __init__(self, bondlist):
        self.bondlist = bondlist
        if isinstance(self.bondlist[0], Bond):
            self.toptype = "bond"
        elif isinstance(self.bondlist[0], Angle):
            self.toptype = "angle"
        elif isinstance(self.bondlist[0], Torsion):
            self.toptype = "torsion"
        else:
            raise ValueError("Input not recognised")

        self._removeDupes()

        self._buildAtomGroups()

        self.len = len(self.atom1)

    def _removeDupes(self):
        """Removes duplicate bonds from a TopologyGroup

        Will rearrange the order of the bondlist within the TopologyGroup after use.
        """
        self.bondlist = list(set(self.bondlist))

    def _buildAtomGroups(self):
        """Builds the "vertical" AtomGroups which are used by the coordinate methods"""
        self.atom1 = AtomGroup.AtomGroup([bond.atom1 for bond in self.bondlist])
        self.atom2 = AtomGroup.AtomGroup([bond.atom2 for bond in self.bondlist])
        if not self.toptype == 'bond':
            self.atom3 = AtomGroup.AtomGroup([bond.atom3 for bond in self.bondlist])
        if self.toptype == 'torsion':
            self.atom4 = AtomGroup.AtomGroup([bond.atom4 for bond in self.bondlist])

    def __len__(self):
        """Number of bonds in the topology group"""
        return self.len

    def __add__(self, other):
        """
        Combine two TopologyGroups together.

        Currently only TopologyGroups of the same type can be combined in such a way
        """
        if not isinstance(other, TopologyGroup):
            raise TypeError("Can only combine two TopologyGroups")
        elif self.toptype != other.toptype:
            raise TypeError("Can only combine TopologyGroups of the same type")

        return TopologyGroup(self.bondlist + other.bondlist)

    def __getitem__(self, item):
        """Returns a particular bond as single object or a subset of
        this TopologyGroup as another TopologyGroup
        """
        if numpy.dtype(type(item)) == numpy.dtype(int): #return a single bond
            return self.bondlist[item]
        elif type(item) == slice: #return a subset of this TopologyGroup
            return TopologyGroup(self.bondlist[item])

    def __iter__(self):
        """Iterator over all bonds"""
        return iter(self.bondlist)

    def __contains__(self, item):
        """Tests if this TopologyGroup contains a bond"""
        return item in self.bondlist

    def __repr__(self):
        return "< " + self.__class__.__name__ + " containing " + str(len(self)) + " " + self.toptype + "s >"

    def _bondsSlow(self, pbc=False):
        """Slow version of bond (numpy implementation)"""
        if not self.toptype == 'bond':
            return "This TopologyGroup is not a bond group!"
        else:
            bond_dist = self.atom1.coordinates() - self.atom2.coordinates()
            if pbc:
                box = self.atom1.dimensions
                if (box[6:9] == 90.).all() and not (box[0:3] == 0).any():  # orthogonal and divide by zero check
                    bond_dist -= numpy.rint(bond_dist / box[0:3]) * box[0:3]
                else:
                    raise ValueError("Only orthogonal boxes supported")

            return numpy.array([norm(a) for a in bond_dist])

    def bonds(self, pbc=False, result=None):
        """Calculates the distance between all bonds in this TopologyGroup

        :Keywords:
           *pbc*
              apply periodic boundary conditions when calculating distance [False]
           *result*
              allows a predefined results array to be used, note that this will be overwritten

        Uses cython implementation
        """
        if not self.toptype == 'bond':
            raise TypeError("TopologyGroup is not of 'bond' type")
        if not result:
            result = numpy.zeros((self.len,), numpy.float64)
        if pbc:
            return distances.calc_bonds(self.atom1.coordinates(), self.atom2.coordinates(),
                                        box=self.atom1.dimensions, result=result)
        else:
            return distances.calc_bonds(self.atom1.coordinates(), self.atom2.coordinates(),
                                        result=result)

    def _anglesSlow(self):
        """Slow version of angle (numpy implementation)"""
        if not self.toptype == 'angle':
            raise TypeError("TopologyGroup is not of type 'angle'")
        from MDAnalysis.core.util import angle as slowang
        from itertools import izip

        vec1 = self.atom1.coordinates() - self.atom2.coordinates()
        vec2 = self.atom3.coordinates() - self.atom2.coordinates()

        angles = numpy.array([slowang(a, b) for a, b in izip(vec1, vec2)])
        return angles

    def angles(self, result=None, pbc=False):
        """Calculates the angle in radians formed between a bond
        between atoms 1 and 2 and a bond between atoms 2 & 3

        :Keywords:
           *result*
              allows a predefined results array to be used, note that this will be overwritten
           *pbc*
              apply periodic boundary conditions when calculating angles [False]
              this is important when connecting vectors between atoms might require 
              minimum image convention

        Uses cython implementation

        .. versionchanged :: 0.8.2
           Added pbc option (default False)
        """
        if not self.toptype == 'angle':
            raise TypeError("topology group is not of type 'angle'")
        if not result:
            result = numpy.zeros((self.len,), numpy.float64)
        if pbc:
            return distances.calc_angles(self.atom1.coordinates(), self.atom2.coordinates(),
                                         self.atom3.coordinates(), box=self.atom1.dimensions, result=result)
        else:
            return distances.calc_angles(self.atom1.coordinates(), self.atom2.coordinates(),
                                         self.atom3.coordinates(), result=result)

    def _torsionsSlow(self):
        """Slow version of torsion (numpy implementation)"""
        if not self.toptype == 'torsion':
            raise TypeError("topology group is not of type 'torsion'")
        from MDAnalysis.core.util import dihedral
        from itertools import izip

        vec1 = self.atom2.coordinates() - self.atom1.coordinates()
        vec2 = self.atom3.coordinates() - self.atom2.coordinates()
        vec3 = self.atom4.coordinates() - self.atom3.coordinates()

        return numpy.array([dihedral(a, b, c) for a, b, c in izip(vec1, vec2, vec3)])

    def torsions(self, result=None, pbc=False):
        """Calculate the torsional angle in radians for this topology
        group.

        Defined as the angle between a plane formed by atoms 1, 2 and
        3 and a plane formed by atoms 2, 3 and 4.

        :Keywords:
           *result*
              allows a predefined results array to be used, note that this will be overwritten
           *pbc*
              apply periodic boundary conditions when calculating angles [False]
              this is important when connecting vectors between atoms might require 
              minimum image convention

        Uses cython implementation.
        
        .. versionchanged:: 0.8.2
           Added pbc option (default False)           
        """
        if not self.toptype == 'torsion':
            raise TypeError("topology group is not of type 'torsion'")
        if not result:
            result = numpy.zeros((self.len,), numpy.float64)
        if pbc:
            return distances.calc_torsions(self.atom1.coordinates(), self.atom2.coordinates(),
                                           self.atom3.coordinates(), self.atom4.coordinates(),
                                           box=self.atom1.dimensions, result=result)
        else:
            return distances.calc_torsions(self.atom1.coordinates(), self.atom2.coordinates(),
                                           self.atom3.coordinates(), self.atom4.coordinates(),
                                           result=result)
