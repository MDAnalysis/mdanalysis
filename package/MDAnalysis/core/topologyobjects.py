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
Core Topology Objects --- :mod:`MDAnalysis.core.topologyobjects`
================================================================

The building blocks for MDAnalysis' description of topology

"""
import numbers
import numpy as np
import functools

from ..lib import mdamath
from ..lib.util import cached
from ..lib import util
from ..lib import distances


@functools.total_ordering
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
    .. versionchanged:: 2.6.0
       Updated Atom ID representation order to match that of AtomGroup indices
    """
    __slots__ = ("_ix", "_u", "btype", "_bondtype", "_guessed", "order")

    def __init__(self, ix, universe, type=None, guessed=False, order=None):
        """Create a topology object

        Parameters
        ----------
        ix : numpy array
          indices of the Atoms
        universe : MDAnalysis.Universe
        type : optional
          Type of the bond
        guessed : optional
          If the Bond is guessed
        """
        self._ix = ix
        self._u = universe
        self._bondtype = type
        self._guessed = guessed
        self.order = order

    @property
    def atoms(self):
        """Atoms within this Bond"""
        return self._u.atoms[self._ix]

    @property
    def indices(self):
        """Tuple of indices describing this object

        .. versionadded:: 0.10.0
        """
        return self._ix

    @property
    def universe(self):
        return self._u

    @property
    def type(self):
        """Type of the bond as a tuple

        Note
        ----
        When comparing types, it is important to consider the reverse
        of the type too, i.e.::

            a.type == b.type or a.type == b.type[::-1]

        """
        if self._bondtype is not None:
            return self._bondtype
        else:
            return tuple(self.atoms.types)

    @property
    def is_guessed(self):
        return bool(self._guessed)

    def __hash__(self):
        return hash((self._u, tuple(self.indices)))

    def __repr__(self):
        """Return representation in same order of AtomGroup indices"""
        return "<{cname} between: {conts}>".format(
            cname=self.__class__.__name__,
            conts=", ".join([
                "Atom {0}".format(i)
                for i in self.indices]))

    def __contains__(self, other):
        """Check whether an atom is in this :class:`TopologyObject`"""
        return other in self.atoms

    def __eq__(self, other):
        """Check whether two bonds have identical contents"""
        if not self.universe == other.universe:
            return False
        return (np.array_equal(self.indices, other.indices) or
                np.array_equal(self.indices[::-1], other.indices))

    def __ne__(self, other):
        return not self == other

    def __lt__(self, other):
        return tuple(self.indices) < tuple(other.indices)

    def __getitem__(self, item):
        """Can retrieve a given Atom from within"""
        return self.atoms[item]

    def __iter__(self):
        return iter(self.atoms)

    def __len__(self):
        return len(self._ix)


class Bond(TopologyObject):

    """A bond between two :class:`~MDAnalysis.core.groups.Atom` instances.

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
    btype = 'bond'

    def partner(self, atom):
        """Bond.partner(Atom)

        Returns
        -------
        the other :class:`~MDAnalysis.core.groups.Atom` in this
        bond
        """
        if atom == self.atoms[0]:
            return self.atoms[1]
        elif atom == self.atoms[1]:
            return self.atoms[0]
        else:
            raise ValueError("Unrecognised Atom")

    def length(self, pbc=True):
        """Length of the bond.

        .. versionchanged:: 0.11.0
           Added pbc keyword
        .. versionchanged:: 0.19.0
           Changed default of pbc to True
        """
        box = self.universe.dimensions if pbc else None

        return distances.calc_bonds(self[0].position, self[1].position, box)

    value = length


class Angle(TopologyObject):

    """An angle between three :class:`~MDAnalysis.core.groups.Atom` instances.
    Atom 2 is the apex of the angle

    .. versionadded:: 0.8
    .. versionchanged:: 0.9.0
       Now a subclass of :class:`TopologyObject`; now uses
       :attr:`__slots__` and stores atoms in :attr:`atoms` attribute
    """
    btype = 'angle'

    def angle(self, pbc=True):
        """Returns the angle in degrees of this Angle.

        Angle between atoms 0 and 2 with apex at 1::

              2
             /
            /
           1------0

        Note
        ----
        The numerical precision is typically not better than
        4 decimals (and is only tested to 3 decimals).

        .. versionadded:: 0.9.0
        .. versionchanged:: 0.17.0
           Fixed angles close to 180 giving NaN
        .. versionchanged:: 0.19.0
           Added pbc keyword, default True
        """
        box = self.universe.dimensions if pbc else None

        return np.rad2deg(distances.calc_angles(
            self[0].position, self[1].position, self[2].position, box))

    value = angle


class Dihedral(TopologyObject):

    """Dihedral (dihedral angle) between four
    :class:`~MDAnalysis.core.groups.Atom` instances.

    The dihedral is defined as the angle between the planes formed by
    Atoms (1, 2, 3) and (2, 3, 4).

    .. versionadded:: 0.8
    .. versionchanged:: 0.9.0
       Now a subclass of :class:`TopologyObject`; now uses :attr:`__slots__`
       and stores atoms in :attr:`atoms` attribute.
    .. versionchanged:: 0.11.0
       Renamed to Dihedral (was Torsion)

    """
    # http://cbio.bmt.tue.nl/pumma/uploads/Theory/dihedral.png
    btype = 'dihedral'

    def dihedral(self, pbc=True):
        """Calculate the dihedral angle in degrees.

        Dihedral angle around axis connecting atoms 1 and 2 (i.e. the angle
        between the planes spanned by atoms (0,1,2) and (1,2,3))::

                  3
                  |
            1-----2
           /
          0


        Note
        ----
        The numerical precision is typically not better than
        4 decimals (and is only tested to 3 decimals).

        .. versionadded:: 0.9.0
        .. versionchanged:: 0.19.0
           Added pbc keyword, default True
        """
        box = self.universe.dimensions if pbc else None
        A, B, C, D = self.atoms

        return np.rad2deg(distances.calc_dihedrals(
            A.position, B.position, C.position, D.position, box))

    value = dihedral


# subclass Dihedral to inherit dihedral method
class ImproperDihedral(Dihedral):
    """
    Improper Dihedral (improper dihedral angle) between four
    :class:`~MDAnalysis.core.groups.Atom` instances.

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
    btype = 'improper'

    def improper(self):
        """Improper dihedral angle in degrees.

        Note
        ----
        The numerical precision is typically not better than
        4 decimals (and is only tested to 3 decimals).
        """
        return self.dihedral()


class UreyBradley(TopologyObject):

    """A Urey-Bradley angle between two :class:`~MDAnalysis.core.groups.Atom` instances.
    Two :class:`UreyBradley` instances can be compared with the ``==`` and
    ``!=`` operators. A UreyBradley angle is equal to another if the same atom
    numbers are involved.

    .. versionadded:: 1.0.0
    """
    btype = 'ureybradley'

    def partner(self, atom):
        """UreyBradley.partner(Atom)
        Returns
        -------
        the other :class:`~MDAnalysis.core.groups.Atom` in this
        interaction
        """
        if atom == self.atoms[0]:
            return self.atoms[1]
        elif atom == self.atoms[1]:
            return self.atoms[0]
        else:
            raise ValueError("Unrecognised Atom")

    def distance(self, pbc=True):
        """Distance between the atoms.
        """
        box = self.universe.dimensions if pbc else None
        return distances.calc_bonds(self[0].position, self[1].position, box)

    value = distance


class CMap(TopologyObject):
    """
    Coupled-torsion correction map term between five 
    :class:`~MDAnalysis.core.groups.Atom` instances.

    .. versionadded:: 1.0.0
    """
    btype = 'cmap'


class TopologyDict(object):

    """A customised dictionary designed for sorting the bonds, angles and
    dihedrals present in a group of atoms.

    Usage::

      topologydict = TopologyDict(members)


    TopologyDicts are also built lazily from a :class:`TopologyGroup.topDict`
    attribute.

    The :class:`TopologyDict` collects all the selected topology type from the
    atoms and categorises them according to the types of the atoms within. A
    :class:`TopologyGroup` containing all of a given bond type can be made by
    querying with the appropriate key. The keys to the :class:`TopologyDict`
    are a tuple of the atom types that the bond represents and can be viewed
    using the :meth:`keys` method.

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

    Arguments
    ---------
    members :
      A list of :class:`TopologyObject` instances


    .. versionadded:: 0.8
    .. versionchanged:: 0.9.0
       Changed initialisation to use a list of :class:`TopologyObject`
       instances instead of list of atoms; now used from within
       :class:`TopologyGroup` instead of accessed from :class:`AtomGroup`.

    """

    def __init__(self, topologygroup):
        if not isinstance(topologygroup, TopologyGroup):
            raise TypeError("Can only construct from TopologyGroup")

        self.dict = dict()
        self._u = topologygroup.universe
        self.toptype = topologygroup.btype

        for b in topologygroup:
            btype = b.type
            try:
                self.dict[btype] += [b]
            except KeyError:
                self.dict[btype] = [b]

        # Some force field types define bonds with a type
        # (Ex. 1: 12 or 21), while others define with a tuple of atom types
        # (Ex. 2: ("H", "O")  or ("O", "H")). If the bond type is a tuple
        # then the bond types in our second example are equivalent and one
        # should be removed. If the bonds are defined as an integer then
        # our first example would also be combined if `_removeDupes()`
        # is run.
        if self.dict and isinstance(list(self.dict.keys())[0], tuple):
            self._removeDupes()

    def _removeDupes(self):
        """Sorts through contents and makes sure that there are
        no duplicate keys (through type reversal)
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

    @property
    def universe(self):
        return self._u

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

            bix = np.vstack([s.indices for s in selection])

            return TopologyGroup(bix, self._u, btype=self.toptype)
        else:
            raise KeyError(key)

    def __contains__(self, other):
        """
        Returns boolean on whether a given type exists within this dictionary

        For topology groups the key (1,2,3) is considered the same as (3,2,1)
        """
        return other in self.dict or other[::-1] in self.dict


_BTYPE_TO_SHAPE = {'bond': 2, 'ureybradley': 2, 'angle': 3, 
                   'dihedral': 4, 'improper': 4, 'cmap': 5}


class TopologyGroup(object):

    """A container for a groups of bonds.

    All bonds of a certain types can be retrieved from within the
    :class:`TopologyGroup` by querying with a tuple of types::

      tg2 = tg.select_bonds([key])

    Where *key* describes the desired bond as a tuple of the involved
    :class:`~MDAnalysis.core.groups.Atom` types, as defined by the .type Atom
    attribute). A list of available keys can be displayed using the
    :meth:`types` method.

    Alternatively, all the bonds which are in a given
    :class:`~MDAnalysis.core.groups.AtomGroup` can be extracted using
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
    .. versionchanged:: 0.19.0
       Empty TopologyGroup now returns correctly shaped empty array via
       indices property and to_indices()
    .. versionchanged::1.0.0
       ``type``, ``guessed``, and ``order`` are no longer reshaped to arrays
       with an extra dimension
    """
    def __init__(self, bondidx, universe, btype=None, type=None, guessed=None,
                 order=None):
        if btype is None:
            # guess what I am
            # difference between dihedral and improper
            # not really important
            self.btype = {2: 'bond',
                          3: 'angle',
                          4: 'dihedral'}[len(bondidx[0])]
        elif btype in _BTYPE_TO_SHAPE:
            self.btype = btype
        else:
            raise ValueError("Unsupported btype, use one of '{}'"
                             "".format(', '.join(_BTYPE_TO_SHAPE)))

        bondidx = np.asarray(bondidx)
        nbonds = len(bondidx)
        # remove duplicate bonds
        if type is None:
            type = np.repeat(None, nbonds)
        if guessed is None:
            guessed = np.repeat(True, nbonds)
        elif guessed is True or guessed is False:
            guessed = np.repeat(guessed, nbonds)
        else:
            guessed = np.asarray(guessed, dtype=bool)
        if order is None:
            order = np.repeat(None, nbonds)

        if nbonds > 0:
            uniq, uniq_idx = util.unique_rows(bondidx, return_index=True)

            self._bix = uniq
            self._bondtypes = type[uniq_idx]
            self._guessed = guessed[uniq_idx]
            self._order = order[uniq_idx]

            # Create vertical AtomGroups
            self._ags = [universe.atoms[self._bix[:, i]]
                         for i in range(self._bix.shape[1])]
        else:
            # Empty TopologyGroup
            self._bix = np.array([])
            self._bondtypes = np.array([])
            self._guessed = np.array([])
            self._order = np.array([])
            self._ags = []
        self._u = universe

        self._cache = dict()  # used for topdict saving

    @property
    def universe(self):
        """The Universe that we belong to"""
        return self._u

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
        return TopologyDict(self)

    def atomgroup_intersection(self, ag, **kwargs):
        """Retrieve all bonds from within this TopologyGroup that are within
        the AtomGroup which is passed.

        Parameters
        ----------
        ag : AtomGroup
            The `:class:~MDAnalysis.core.groups.AtomGroup` to intersect
            with.
        strict : bool
            Only retrieve bonds which are completely contained within the
            AtomGroup. [``False``]


        .. versionadded:: 0.9.0
        """
        # Issue #780 - if self is empty, return self to avoid invalid mask
        if not self:
            return self

        # Strict requires all items in a row to be seen,
        # otherwise any item in a row
        func = np.all if kwargs.get('strict', False) else np.any

        atom_idx = ag.indices
        # Create a list of boolean arrays,
        # each representing a column of bond indices.
        seen = [np.isin(col, atom_idx) for col in self._bix.T]

        # Create final boolean mask by summing across rows
        mask = func(seen, axis=0)

        return self[mask]

    @property
    def indices(self):
        """all bond indices

        See Also
        --------
        to_indices : function that just returns `indices`
        """
        if not self:
            # empty TG
            shape = _BTYPE_TO_SHAPE[self.btype]
            return np.zeros((0, shape), dtype=np.int32)
        else:
            return self._bix

    def to_indices(self):
        """Return a data structure with atom indices describing the bonds.

        This format should be identical to the original contents of the
        entries in universe._topology.
        Note that because bonds are sorted as they are initialised, the order
        that atoms are defined in each entry might be reversed.

        Returns
        -------
        indices : tuple
            A tuple of tuples which define the contents of this
            TopologyGroup in terms of the atom numbers.  (0 based
            index within u.atoms)

        .. versionadded:: 0.9.0
        .. versionchanged:: 0.10.0
           Renamed from "dump_contents" to "to_indices"

        """
        return self.indices

    dump_contents = to_indices

    def __len__(self):
        """Number of bonds in the topology group"""
        return self._bix.shape[0]

    def __add__(self, other):
        """Combine two TopologyGroups together.

        Can combined two TopologyGroup of the same type, or add a single
        TopologyObject to a TopologyGroup.
        """
        # check addition is sane
        if not isinstance(other, (TopologyObject, TopologyGroup)):
            raise TypeError("Can only combine TopologyObject or "
                            "TopologyGroup to TopologyGroup, not {0}"
                            "".format(type(other)))

        # cases where either other or self is empty TG
        if not other:  # adding empty TG to me
            return self
        if not self:
            if isinstance(other, TopologyObject):
                # Reshape indices to be 2d array
                return TopologyGroup(other.indices[None, :],
                                     other.universe,
                                     btype=other.btype,
                                     type=np.array([other._bondtype]),
                                     guessed=np.array([other.is_guessed]),
                                     order=np.array([other.order]),
                                     )
            else:
                return TopologyGroup(other.indices,
                                     other.universe,
                                     btype=other.btype,
                                     type=other._bondtypes,
                                     guessed=other._guessed,
                                     order=other._order,
                                     )
        else:
            if not other.btype == self.btype:
                raise TypeError("Cannot add different types of "
                                "TopologyObjects together")
            if isinstance(other, TopologyObject):
                # add TO to me
                return TopologyGroup(
                    np.concatenate([self.indices, other.indices[None, :]]),
                    self.universe,
                    btype=self.btype,
                    type=np.concatenate([self._bondtypes,
                                         np.array([other._bondtype])]),
                    guessed=np.concatenate([self._guessed,
                                            np.array([other.is_guessed])]),
                    order=np.concatenate([self._order,
                                          np.array([other.order])]),
                )
            else:
                # add TG to me
                return TopologyGroup(
                    np.concatenate([self.indices, other.indices]),
                    self.universe,
                    btype=self.btype,
                    type=np.concatenate([self._bondtypes, other._bondtypes]),
                    guessed=np.concatenate([self._guessed, other._guessed]),
                    order=np.concatenate([self._order, other._order]),
                )

    def __getitem__(self, item):
        """Returns a particular bond as single object or a subset of
        this TopologyGroup as another TopologyGroup

        .. versionchanged:: 0.10.0
           Allows indexing via boolean numpy array
        """
        # Grab a single Item, similar to Atom/AtomGroup relationship
        if isinstance(item, numbers.Integral):
            outclass = {'bond': Bond,
                        'angle': Angle,
                        'dihedral': Dihedral,
                        'improper': ImproperDihedral,
                        'ureybradley': UreyBradley,
                        'cmap': CMap}[self.btype]
            return outclass(self._bix[item],
                            self._u,
                            type=self._bondtypes[item],
                            guessed=self._guessed[item],
                            order=self._order[item])
        else:
            # Slice my index array with the item
            return self.__class__(self._bix[item],
                                  self._u,
                                  btype=self.btype,
                                  type=self._bondtypes[item],
                                  guessed=self._guessed[item],
                                  order=self._order[item],)

    def __contains__(self, item):
        """Tests if this TopologyGroup contains a bond"""
        return item.indices in self._bix

    def __repr__(self):
        return "<TopologyGroup containing {num} {type}s>".format(
            num=len(self), type=self.btype)

    def __eq__(self, other):
        """Test if contents of TopologyGroups are equal"""
        return np.array_equal(self.indices, other.indices)

    def __ne__(self, other):
        return not self == other

    def __nonzero__(self):
        return not len(self) == 0

    @property
    def atom1(self):
        """The first atom in each TopologyObject in this Group"""
        return self._ags[0]

    @property
    def atom2(self):
        """The second atom in each TopologyObject in this Group"""
        return self._ags[1]

    @property
    def atom3(self):
        """The third atom in each TopologyObject in this Group"""
        try:
            return self._ags[2]
        except IndexError:
            nvert = _BTYPE_TO_SHAPE[self.btype]
            errmsg = (f"TopologyGroup of {self.btype}s only has {nvert} "
                      f"vertical AtomGroups")
            raise IndexError(errmsg) from None

    @property
    def atom4(self):
        """The fourth atom in each TopologyObject in this Group"""
        try:
            return self._ags[3]
        except IndexError:
            nvert = _BTYPE_TO_SHAPE[self.btype]
            errmsg = (f"TopologyGroup of {self.btype}s only has {nvert} "
                      f"vertical AtomGroups")
            raise IndexError(errmsg) from None

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
        if self.btype == 'bond':
            return self.bonds(**kwargs)
        elif self.btype == 'angle':
            return self.angles(**kwargs)
        elif self.btype == 'dihedral':
            return self.dihedrals(**kwargs)
        elif self.btype == 'improper':
            return self.dihedrals(**kwargs)

    def _calc_connection_values(self, func, *btypes, result=None, pbc=False):
        if not any(self.btype == btype for btype in btypes):
            strbtype = "' or '".join(btypes)
            raise TypeError(f"TopologyGroup is not of type '{strbtype}'")
        if not result:
            result = np.zeros(len(self), np.float64)
        box = None if not pbc else self._ags[0].dimensions
        positions = [ag.positions for ag in self._ags]
        return func(*positions, box=box, result=result)

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
        return self._calc_connection_values(distances.calc_bonds, "bond",
                                            pbc=pbc, result=result)

    def angles(self, result=None, pbc=False):
        """Calculates the angle in radians formed between a bond
        between atoms 1 and 2 and a bond between atoms 2 & 3

        Parameters
        ----------
        result : array_like
            allows a predefined results array to be used, note that this
            will be overwritten
        pbc : bool
            apply periodic boundary conditions when calculating angles
            [``False``] this is important when connecting vectors between
            atoms might require minimum image convention

        Returns
        -------
        angles : ndarray

        .. versionchanged :: 0.9.0
           Added *pbc* option (default ``False``)

        """
        return self._calc_connection_values(distances.calc_angles, "angle",
                                            pbc=pbc, result=result)

    def dihedrals(self, result=None, pbc=False):
        """Calculate the dihedral angle in radians for this topology
        group.

        Defined as the angle between a plane formed by atoms 1, 2 and
        3 and a plane formed by atoms 2, 3 and 4.

        Parameters
        ----------
        result : array_like
            allows a predefined results array to be used, note that this
            will be overwritten
        pbc : bool
            apply periodic boundary conditions when calculating angles
            [``False``] this is important when connecting vectors between
            atoms might require minimum image convention

        Returns
        -------
        angles : ndarray

        .. versionchanged:: 0.9.0
           Added *pbc* option (default ``False``)
        """
        return self._calc_connection_values(distances.calc_dihedrals,
                                            "dihedral", "improper",
                                            pbc=pbc, result=result)
