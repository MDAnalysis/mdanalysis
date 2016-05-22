# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.MDAnalysis.org
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
Core Topology Objects --- :mod:`MDAnalysis.core.topologyobjects`
================================================================

The building blocks for MDAnalysis' description of topology

"""
from __future__ import print_function, absolute_import

from six.moves import zip
import numpy as np

from ..lib.mdamath import norm, dihedral
from ..lib.mdamath import angle as slowang
from ..lib.util import cached
from ..lib import distances

class TopologyObject(object):

    """Base class for all Topology items.

    Defines the behaviour by which Bonds/Angles/etc in MDAnalysis should
    behave.

    .. versionadded:: 0.9.0
    .. versionchanged:: 0.10.0
       All TopologyObject now keep track of if they were guessed or not
       via the ``is_guessed`` managed property.
    .. versionadded:: 0.11.0
       Added the `value` method to return the size of the object
    """
    __slots__ = ("atoms", "_is_guessed")

    def __init__(self, atoms, is_guessed=False):
        self.atoms = tuple(atoms)
        self._is_guessed = is_guessed

    @property
    def indices(self):
        """Tuple of indices describing this object

        .. versionadded:: 0.10.0
        """
        return tuple([a.index for a in self.atoms])

    @property
    def type(self):
        """Type of the bond as a tuple

        When comparing types, it is important to consider the reverse
        of the type too, i.e.::

            a.type == b.type or a.type == b.type[::-1]

        """
        return tuple([a.type for a in self.atoms])

    @property
    def is_guessed(self):  # is property so it gets nice docs?
        """``True`` if the bond was guessed.

        .. SeeAlso:: :func:`guess_bonds` :func:`guess_angles` and
                     :func:`guess_dihedrals`
        """
        return self._is_guessed

    @is_guessed.setter
    def is_guessed(self, b):
        self._is_guessed = b

    def __repr__(self):
        return "<{cname} between: {conts}>".format(
            cname=self.__class__.__name__,
            conts=", ".join([
                "Atom {num} ({name} of {resname}-{resid})".format(
                    num=a.index + 1,
                    name=a.name,
                    resname=a.resname,
                    resid=a.resid)
                for a in self.atoms]))

    def __contains__(self, other):
        """Check whether an atom is in this :class:`TopologyObject`"""
        return other in self.atoms

    def __eq__(self, other):
        """Check whether two bonds have identical contents"""
        my_tup = tuple([a.index for a in self.atoms])
        ot_tup = tuple([a.index for a in other.atoms])

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

    def _cmp_key(self):
        """Unique key for the object to be used to generate the object hash"""
        # This key must be equal for two object considered as equal by __eq__
        return self.__class__, tuple(sorted(self.indices))

    def __hash__(self):
        """Makes the object hashable"""
        return hash(self._cmp_key())


class Bond(TopologyObject):

    """A bond between two :class:`~MDAnalysis.core.AtomGroup.Atom` instances.

    Two :class:`Bond` instances can be compared with the ``==`` and
    ``!=`` operators. A bond is equal to another if the same atom
    numbers are connected and they have the same bond order. The
    ordering of the two atom numbers is ignored as is the fact that a
    bond was guessed.

    The presence of a particular atom can also be queried::

      >>> Atom in Bond

    will return either ``True`` or ``False``.

    .. versionchanged:: 0.9.0
       Now a subclass of :class:`TopologyObject`. Changed class to use
       :attr:`__slots__` and stores atoms in :attr:`atoms` attribute.
    """
    __slots__ = ("atoms", "order", "_is_guessed")

    def __init__(self, atoms, order=None, is_guessed=False):
        self.atoms = tuple(atoms)
        self.order = order
        self._is_guessed = is_guessed

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

    def length(self, pbc=False):
        """Length of the bond.

        .. versionchanged:: 0.11.0
           Added pbc keyword
        """
        if pbc:
            box = self.atoms[0].universe.dimensions
            return distances.self_distance_array(
                np.array([self.atoms[0].position, self.atoms[1].position]),
                box=box)
        else:
            return norm(self.atoms[0].position - self.atoms[1].position)

    value = length

    def __repr__(self):
        a1, a2 = self.atoms
        s_id = "<Bond between: Atom {0:d} ({1.name} of {1.resname} {1.resid}"\
               " {1.altLoc}) and Atom {2:d} ({3.name} of {3.resname}"\
               "{3.resid} {3.altLoc})".format(
                   a1.index + 1, a1, a2.index + 1, a2)
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
       Now a subclass of :class:`TopologyObject`; now uses
       :attr:`__slots__` and stores atoms in :attr:`atoms` attribute
    """

    def angle(self):
        """Returns the angle in degrees of this Angle.

        Angle between atoms 0 and 2 with apex at 1::

              2
             /
            /
           1------0

        .. Note:: The numerical precision is typically not better than
                  4 decimals (and is only tested to 3 decimals).

        .. versionadded:: 0.9.0
        """
        a = self[0].position - self[1].position
        b = self[2].position - self[1].position
        return np.rad2deg(
            np.arccos(np.dot(a, b) / (norm(a) * norm(b))))

    value = angle


class Dihedral(TopologyObject):

    """Dihedral (dihedral angle) between four
    :class:`~MDAnalysis.core.AtomGroup.Atom` instances.

    The dihedral is defined as the angle between the planes formed by
    Atoms (1, 2, 3) and (2, 3, 4).

    .. versionadded:: 0.8
    .. versionchanged:: 0.9.0
       Now a subclass of :class:`TopologyObject`; now uses :attr:`__slots__` and
       stores atoms in :attr:`atoms` attribute.
    .. versionchanged:: 0.11.0
       Renamed to Dihedral (was Torsion)
    """
    # http://cbio.bmt.tue.nl/pumma/uploads/Theory/dihedral.png

    def dihedral(self):
        """Calculate the dihedral angle in degrees.

        Dihedral angle around axis connecting atoms 1 and 2 (i.e. the angle
        between the planes spanned by atoms (0,1,2) and (1,2,3))::

                  3
                  |
            1-----2
           /
          0


        .. Note:: The numerical precision is typically not better than
                  4 decimals (and is only tested to 3 decimals).

        .. versionadded:: 0.9.0
        """
        A, B, C, D = self.atoms
        ab = A.position - B.position
        bc = B.position - C.position
        cd = C.position - D.position
        return np.rad2deg(dihedral(ab, bc, cd))

    value = dihedral


class ImproperDihedral(Dihedral):  # subclass Dihedral to inherit dihedral method

    """
    Improper Dihedral (improper dihedral angle) between four
    :class:`~MDAnalysis.core.AtomGroup.Atom` instances.

    MDAnalysis treats the improper dihedral angle as the angle between
    the planes formed by Atoms (1, 2, 3) and (2, 3, 4).

    .. warning:: Definitions of Atom ordering in improper dihedrals
                 can change. Check the definitions here against
                 your software.

    .. versionadded:: 0.9.0
    .. versionchanged:: 0.11.0
       Renamed to ImproperDihedral (was Improper_Torsion)
    """
    # http://cbio.bmt.tue.nl/pumma/uploads/Theory/improper.png

    def improper(self):
        """Improper dihedral angle in degrees.

        .. Note:: The numerical precision is typically not better than
                  4 decimals (and is only tested to 3 decimals).
        """
        return self.dihedral()


class TopologyDict(object):

    """A customised dictionary designed for sorting the bonds, angles and
    dihedrals present in a group of atoms.

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

    The :class:`TopologyDict` collects all the selected topology type from the
    atoms and categorises them according to the types of the atoms within.
    A :class:`TopologyGroup` containing all of a given bond type can
    be made by querying with the appropriate key.  The keys to the
    :class:`TopologyDict` are a tuple of the atom types that the bond represents
    and can be viewed using the :meth:`keys` method.

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
    and 2 having the bond in their :attr:`bond` attribute.

    Two :class:`TopologyDict` instances can be combined using
    addition and it will not create any duplicate bonds in the process.

    .. versionadded:: 0.8
    .. versionchanged:: 0.9.0
       Changed initialisation to use a list of :class:`TopologyObject`
       instances instead of list of atoms; now used from within
       :class:`TopologyGroup` instead of accessed from :class:`AtomGroup`.
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

      tg2 = tg.select_bonds([key])

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

    The :meth:`bonds`, :meth:`angles` and :meth:`dihedrals` methods offer
    a "shortcut" to the Cython distance calculation functions in
    :class:`MDAnalysis.lib.distances`.

    TopologyGroups can be combined with TopologyGroups of the same bond
    type (ie can combine two angle containing TopologyGroups).

    .. versionadded:: 0.8
    .. versionchanged:: 0.9.0
       Overhauled completely: (1) Added internal :class:`TopologyDict`
       accessible by the :attr:`topDict` attribute. (2)
       :meth:`selectBonds` allows the :attr:`topDict` to be queried
       with tuple of types. (3) Added :meth:`atomgroup_intersection`
       to allow bonds which are in a given :class:`AtomGroup` to be retrieved.
    .. versionchanged:: 0.10.0
       Added :func:`from_indices` constructor, allowing class to be created
       from indices.
       Can now create empty Group.
       Renamed :meth:`dump_contents` to :meth:`to_indices`
    .. versionchanged:: 0.11.0
       Added `values` method to return the size of each object in this group
       Deprecated selectBonds method in favour of select_bonds
    """

    def __init__(self, bondlist):
        try:
            self.toptype = bondlist[0].__class__.__name__
        except IndexError:
            self.toptype = None
        else:
            # only check type if list has length
            if not isinstance(bondlist[0], TopologyObject):
                raise TypeError("Input must be TopologyObject")
        # Would be nice to make everything work internally using sets, BUT
        # sets can't be indexed, so couldn't work backward from .angles()
        # results to find which angle is a certain value.
        # Sorted so that slicing returns sensible results
        self.bondlist = tuple(sorted(set(bondlist)))

        self._cache = dict()  # used for topdict saving

    @classmethod
    def from_indices(cls, bondlist, atomgroup,
                     bondclass=None, guessed=True,
                     remove_duplicates=False):
        """Initiate from a list of indices.

        :Arguments:
          *bondlist*
            A list of lists of indices.  For example `[(0, 1), (1, 2)]`
            Note that these indices refer to the index of the Atoms
            within the supplied AtomGroup, not their global index.
          *atomgroup*
            An AtomGroup which the indices from bondlist will be used on.

        :Keywords:
          *bondclass*
            The Class of the topology object to be made.
            If missing this will try and be guessed according to the number
            of indices in each record.
          *guessed*
            Whether or not the bonds were guessed. [``True``]
          *remove_duplicates*
            Sort through items and make sure that no duplicates are created [``False``]

        .. versionadded:: 0.10.0
        """
        if remove_duplicates:
            # always have first index less than last
            bondlist = {b if b[0] < b[-1] else b[::-1] for b in bondlist}

        if bondclass is None:  # try and guess
            try:
                bondclass = {
                    2:Bond,
                    3:Angle,
                    4:Dihedral
                }[len(bondlist[0])]
            except KeyError:
                raise ValueError("Can't detect bondclass for provided indices")

        bonds = [bondclass([atomgroup[a] for a in entry], is_guessed=guessed)
                 for entry in bondlist]

        return cls(bonds)

    def select_bonds(self, selection):
        """Retrieves a selection from this topology group based on types.

        .. seeAlso :meth:`types`

        .. versionadded 0.9.0
        """
        return self.topDict[selection]

    selectBonds = select_bonds

    def types(self):
        """Return a list of the bond types in this TopologyGroup

        .. versionadded 0.9.0
        """
        return list(self.topDict.keys())

    @property
    @cached('dict')
    def topDict(self):
        """
        Returns the TopologyDict for this topology group.

        This is used for the select_bonds method when fetching a certain type
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
        elif self.toptype == 'Dihedral':
            other_set = set(other.dihedrals)
        elif self.toptype == 'ImproperDihedral':
            other_set = set(other.impropers)
        else:
            raise ValueError("Unsupported intersection")

        newlist = list(set(self.bondlist).intersection(other_set))

        return TopologyGroup(newlist)

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
        # If this count is equal to crit, (bond=2, angle=3, dihedral=4) then
        # the TopObj was seen enough for it to have to be completely be
        # present in the AtomGroup

        # each bond starts with 0 appearances
        # I'm only interested in intersection, so if its not in tg then
        # i'll get keyerrors which i'll pass
        count_dict = dict.fromkeys(self.bondlist, 0)

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
        elif self.toptype == 'Dihedral':
            crit = 4
            for atom in other:
                for b in atom.dihedrals:
                    try:
                        count_dict[b] += 1
                    except KeyError:
                        pass
        elif self.toptype == 'ImproperDihedral':
            crit = 4
            for atom in other:
                for b in atom.impropers:
                    try:
                        count_dict[b] += 1
                    except KeyError:
                        pass

        # now make new list, which only includes bonds with enough appearances
        newlist = [b for b in self.bondlist if count_dict[b] == crit]

        return TopologyGroup(newlist)

    def to_indices(self):
        """Return a tuple of tuples which define the contents of this
        TopologyGroup in terms of the atom numbers,
        (0 based index within u.atoms)

        This format should be identical to the original contents of the
        entries in universe._topology.
        Note that because bonds are sorted as they are initialised, the order
        that atoms are defined in each entry might be reversed.

        .. versionadded:: 0.9.0
        .. versionchanged:: 0.10.0
           Renamed from "dump_contents" to "to_indices"
        """
        # should allow topology information to be pickled even if it is
        # substantially changed from original input,
        # eg through merging universes or defining new bonds manually.
        bondlist = tuple([b.indices for b in self.bondlist])

        return bondlist

    dump_contents = to_indices

    @property
    @cached('atom1')
    def atom1(self):
        from .AtomGroup import AtomGroup
        return AtomGroup([b[0] for b in self.bondlist])

    @property
    @cached('atom2')
    def atom2(self):
        from .AtomGroup import AtomGroup
        return AtomGroup([b[1] for b in self.bondlist])

    @property
    @cached('atom3')
    def atom3(self):
        from .AtomGroup import AtomGroup
        return AtomGroup([b[2] for b in self.bondlist])

    @property
    @cached('atom4')
    def atom4(self):
        from .AtomGroup import AtomGroup
        return AtomGroup([b[3] for b in self.bondlist])

    def __len__(self):
        """Number of bonds in the topology group"""
        return len(self.bondlist)

    def __add__(self, other):
        """Combine two TopologyGroups together.

        Can combined two TopologyGroup of the same type, or add a single
        TopologyObject to a TopologyGroup.
        """
        # check addition is sane
        if not (isinstance(other, TopologyObject)
                or isinstance(other, TopologyGroup)):
            raise TypeError("Can only combine TopologyObject or TopologyGroup to"
                            " TopologyGroup, not {0}".format(type(other)))

        # cases where either other or self is empty TG
        if not other:  # adding empty TG to me
            return self
        if not self:
            if isinstance(other, TopologyObject):
                return TopologyGroup([other])
            else:
                return TopologyGroup(other.bondlist)

        # add TO to me
        if isinstance(other, TopologyObject):
            if not isinstance(other, type(self.bondlist[0])):
                raise TypeError("Cannot add different types of "
                                "TopologyObjects together")
            else:
                return TopologyGroup(self.bondlist + (other,))

        # add TG to me
        if self.toptype != other.toptype:
            raise TypeError("Can only combine TopologyGroups of the same type")
        else:
            return TopologyGroup(self.bondlist + other.bondlist)

    def __getitem__(self, item):
        """Returns a particular bond as single object or a subset of
        this TopologyGroup as another TopologyGroup

        .. versionchanged:: 0.10.0
           Allows indexing via boolean numpy array
        """
        if np.dtype(type(item)) == np.dtype(int):
            return self.bondlist[item]  # single TopObj
        elif isinstance(item, slice):
            return TopologyGroup(self.bondlist[item])  # new TG
        elif isinstance(item, (np.ndarray, list)):
            try:
                if isinstance(item[0], np.bool_):
                    item = np.arange(len(item))[item]
            except IndexError:  # zero length item
                pass
            return TopologyGroup([self.bondlist[i] for i in item])

    def __iter__(self):
        """Iterator over all bonds"""
        return iter(self.bondlist)

    def __contains__(self, item):
        """Tests if this TopologyGroup contains a bond"""
        return item in self.bondlist

    def __repr__(self):
        return "<TopologyGroup containing {num} {type}s>".format(
            num=len(self), type=self.toptype)

    def __eq__(self, other):
        """Test if contents of TopologyGroups are equal"""
        return set(self) == set(other)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __nonzero__(self):
        return not len(self) == 0

    # Distance calculation methods below
    # "Slow" versions exist as a way of testing the Cython implementations
    def values(self, **kwargs):
        """Return the size of each object in this Group

        :Keywords:
           *pbc*
              apply periodic boundary conditions when calculating distance
              [``False``]
           *result*
              allows a predefined results array to be used,
              note that this will be overwritten

        .. versionadded:: 0.11.0
        """
        if self.toptype == 'Bond':
            return self.bonds(**kwargs)
        elif self.toptype == 'Angle':
            return self.angles(**kwargs)
        elif self.toptype == 'Dihedral':
            return self.dihedrals(**kwargs)
        elif self.toptype == 'ImproperDihedral':
            return self.dihedrals(**kwargs)

    def _bondsSlow(self, pbc=False):  # pragma: no cover
        """Slow version of bond (numpy implementation)"""
        if not self.toptype == 'Bond':
            return TypeError("TopologyGroup is not of type 'Bond'")
        else:
            bond_dist = self.atom1.positions - self.atom2.positions
            if pbc:
                box = self.atom1.dimensions
                # orthogonal and divide by zero check
                if (box[6:9] == 90.).all() and not (box[0:3] == 0).any():
                    bond_dist -= np.rint(bond_dist / box[0:3]) * box[0:3]
                else:
                    raise ValueError("Only orthogonal boxes supported")

            return np.array([norm(a) for a in bond_dist])

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
            result = np.zeros(len(self), np.float64)
        if pbc:
            return distances.calc_bonds(self.atom1.positions,
                                        self.atom2.positions,
                                        box=self.atom1.dimensions,
                                        result=result)
        else:
            return distances.calc_bonds(self.atom1.positions,
                                        self.atom2.positions,
                                        result=result)

    def _anglesSlow(self):  # pragma: no cover
        """Slow version of angle (numpy implementation)"""
        if not self.toptype == 'Angle':
            raise TypeError("TopologyGroup is not of type 'Angle'")

        vec1 = self.atom1.positions - self.atom2.positions
        vec2 = self.atom3.positions - self.atom2.positions

        angles = np.array([slowang(a, b) for a, b in zip(vec1, vec2)])
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
              [``False``] this is important when connecting vectors between atoms
              might require minimum image convention

        Uses cython implementation

        .. versionchanged :: 0.9.0
           Added *pbc* option (default ``False``)
        """
        if not self.toptype == 'Angle':
            raise TypeError("TopologyGroup is not of type 'Angle'")
        if not result:
            result = np.zeros(len(self), np.float64)
        if pbc:
            return distances.calc_angles(self.atom1.positions,
                                         self.atom2.positions,
                                         self.atom3.positions,
                                         box=self.atom1.dimensions,
                                         result=result)
        else:
            return distances.calc_angles(self.atom1.positions,
                                         self.atom2.positions,
                                         self.atom3.positions,
                                         result=result)

    def _dihedralsSlow(self):  # pragma: no cover
        """Slow version of dihedral (numpy implementation)"""
        if self.toptype not in ['Dihedral', 'ImproperDihedral']:
            raise TypeError("TopologyGroup is not of type 'Dihedral' or "
                            "'ImproperDihedral'")

        vec1 = self.atom2.positions - self.atom1.positions
        vec2 = self.atom3.positions - self.atom2.positions
        vec3 = self.atom4.positions - self.atom3.positions

        return np.array([dihedral(a, b, c)
                         for a, b, c in zip(vec1, vec2, vec3)])

    def dihedrals(self, result=None, pbc=False):
        """Calculate the dihedralal angle in radians for this topology
        group.

        Defined as the angle between a plane formed by atoms 1, 2 and
        3 and a plane formed by atoms 2, 3 and 4.

        :Keywords:
           *result*
              allows a predefined results array to be used, note that this
              will be overwritten
           *pbc*
v              apply periodic boundary conditions when calculating angles
              [``False``] this is important when connecting vectors between
              atoms might require minimum image convention

        Uses cython implementation.

        .. versionchanged:: 0.9.0
           Added *pbc* option (default ``False``)
        """
        if self.toptype not in ['Dihedral', 'ImproperDihedral']:
            raise TypeError("TopologyGroup is not of type 'Dihedral' or "
                            "'ImproperDihedral'")
        if not result:
            result = np.zeros(len(self), np.float64)
        if pbc:
            return distances.calc_dihedrals(self.atom1.positions,
                                            self.atom2.positions,
                                            self.atom3.positions,
                                            self.atom4.positions,
                                            box=self.atom1.dimensions,
                                            result=result)
        else:
            return distances.calc_dihedrals(self.atom1.positions,
                                            self.atom2.positions,
                                            self.atom3.positions,
                                            self.atom4.positions,
                                            result=result)
