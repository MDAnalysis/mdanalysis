# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""\
==========================================================
Core objects: Containers --- :mod:`MDAnalysis.core.groups`
==========================================================

The :class:`~MDAnalysis.core.universe.Universe` instance contains all the
particles in the system (which MDAnalysis calls :class:`Atom`). Groups of
:class:`atoms<Atom>` are handled as :class:`AtomGroup` instances. The
:class:`AtomGroup` is probably the most important object in MDAnalysis because
virtually everything can be accessed through it. :class:`AtomGroup` instances
can be easily created (e.g., from an :meth:`AtomGroup.select_atoms` selection or
simply by slicing).

For convenience, chemically meaningful groups of :class:`Atoms<Atom>` such as a
:class:`Residue` or a :class:`Segment` (typically a whole molecule or all of the
solvent) also exist as containers, as well as groups of these units
(:class:`ResidueGroup`, :class:`SegmentGroup`).


Classes
=======

Collections
-----------

.. autoclass:: AtomGroup
   :members:
   :inherited-members:
.. autoclass:: ResidueGroup
   :members:
   :inherited-members:
.. autoclass:: SegmentGroup
   :members:
   :inherited-members:
.. autoclass:: UpdatingAtomGroup
   :members:

Chemical units
--------------

.. autoclass:: Atom
   :members:
   :inherited-members:
.. autoclass:: Residue
   :members:
   :inherited-members:
.. autoclass:: Segment
   :members:
   :inherited-members:

Levels
------

Each of the above classes has a *level* attribute.  This can be used to verify
that two objects are of the same level, or to access a particular class::

   u = mda.Universe()

   ag = u.atoms[:10]
   at = u.atoms[11]

   ag.level == at.level  # Returns True

   ag.level.singular  # Returns Atom class
   at.level.plural  # Returns AtomGroup class

"""
from __future__ import absolute_import
from six.moves import zip
from six import string_types

from collections import namedtuple
import numpy as np
import functools
import itertools
import numbers
import os
import warnings

from numpy.lib.utils import deprecate

from .. import _ANCHOR_UNIVERSES
from ..lib import util
from ..lib.util import cached, warn_if_not_unique, unique_int_1d
from ..lib import distances
from ..lib import transformations
from ..selections import get_writer as get_selection_writer_for
from . import selection
from . import flags
from ..exceptions import NoDataError
from . import topologyobjects
from ._get_readers import get_writer_for


def _unpickle(uhash, ix):
    try:
        u = _ANCHOR_UNIVERSES[uhash]
    except KeyError:
        # doesn't provide as nice an error message as before as only hash of universe is stored
        # maybe if we pickled the filename too we could do better...
        raise RuntimeError(
            "Couldn't find a suitable Universe to unpickle AtomGroup onto "
            "with Universe hash '{}'.  Available hashes: {}"
            "".format(uhash, ', '.join([str(k)
                                        for k in _ANCHOR_UNIVERSES.keys()])))
    return u.atoms[ix]

def _unpickle_uag(basepickle, selections, selstrs):
    bfunc, bargs = basepickle[0], basepickle[1:][0]
    basegroup = bfunc(*bargs)
    return UpdatingAtomGroup(basegroup, selections, selstrs)


def make_classes():
    """Make a fresh copy of all classes

    Returns
    -------
    Two dictionaries. One with a set of :class:`_TopologyAttrContainer` classes
    to serve as bases for :class:`~MDAnalysis.core.universe.Universe`\ -specific
    MDA container classes. Another with the final merged versions of those
    classes. The classes themselves are used as hashing keys.

    """
    bases = {}
    classes = {}
    groups = (AtomGroup, ResidueGroup, SegmentGroup)
    components = (Atom, Residue, Segment)

    # The 'GBase' middle man is needed so that a single topologyattr
    #  patching applies automatically to all groups.
    GBase = bases[GroupBase] = _TopologyAttrContainer._subclass(is_group=True)
    for cls in groups:
        bases[cls] = GBase._subclass(is_group=True)
    # CBase for patching all components
    CBase = bases[ComponentBase] = _TopologyAttrContainer._subclass(
        is_group=False)
    for cls in components:
        bases[cls] = CBase._subclass(is_group=False)

    # Initializes the class cache.
    for cls in groups + components:
        classes[cls] = bases[cls]._mix(cls)

    return bases, classes


class _TopologyAttrContainer(object):
    """Class factory for receiving sets of :class:`TopologyAttr` objects.

    :class:`_TopologyAttrContainer` is a convenience class to encapsulate the
    functions that deal with:
    * the import and namespace transplant of
      :class:`~MDAnalysis.core.topologyattrs.TopologyAttr` objects;
    * the copying (subclassing) of itself to create distinct bases for the
      different container classes (:class:`AtomGroup`, :class:`ResidueGroup`,
      :class:`SegmentGroup`, :class:`Atom`, :class:`Residue`, :class:`Segment`,
      and subclasses thereof);
    * the mixing (subclassing and co-inheritance) with the container classes.
      The mixed subclasses become the final container classes specific to each
      :class:`~MDAnalysis.core.universe.Universe`.
    """
    @classmethod
    def _subclass(cls, is_group):
        """Factory method returning :class:`_TopologyAttrContainer` subclasses.

        Parameters
        ----------
        is_group : bool
            The :attr:`_is_group` of the returned class will be set to
            `is_group`. This is used to distinguish between Groups
            (:class:`AtomGroup` etc.) and Components (:class:`Atom` etc.) in
            internal methods when considering actions such as addition of
            objects or adding
            :class:`TopologyAttributes<MDAnalysis.core.topologyattrs.TopologyAttr>`
            to them.

        Returns
        -------
        type
            A subclass of :class:`_TopologyAttrContainer`, with the same name.
        """
        newcls = type(cls.__name__, (cls,), {'_is_group': bool(is_group)})
        if is_group:
            newcls._SETATTR_WHITELIST = {
                'positions', 'velocities', 'forces', 'dimensions',
                'atoms', 'residue', 'residues', 'segment', 'segments',
            }
        else:
            newcls._SETATTR_WHITELIST = {
                'position', 'velocity', 'force', 'dimensions',
                'atoms', 'residue', 'residues', 'segment',
            }

        return newcls

    @classmethod
    def _mix(cls, other):
        """Creates a subclass with ourselves and another class as parents.

        Classes mixed at this point override :meth:`__new__`, causing further
        instantiations to shortcut to :meth:`~object.__new__` (skipping the
        cache-fetch process for :class:`_MutableBase` subclasses).

        The new class will have an attribute `_derived_class` added, pointing
        to itself. This pointer instructs which class to use when
        slicing/adding instances of the new class. At initialization time, the
        new class may choose to point `_derived_class` to another class (as is
        done in the initialization of :class:`UpdatingAtomGroup`).

        Parameters
        ----------
        other : type
            The class to mix with ourselves.

        Returns
        -------
        type
            A class of parents :class:`_ImmutableBase`, *other* and this class.
            Its name is the same as *other*'s.
        """
        newcls = type(other.__name__, (_ImmutableBase, other, cls), {})
        newcls._derived_class = newcls
        return newcls

    @classmethod
    def _add_prop(cls, attr):
        """Add `attr` into the namespace for this class

        Parameters
        ----------
        attr : A :class:`TopologyAttr` object
        """
        def getter(self):
            return attr.__getitem__(self)

        def setter(self, values):
            return attr.__setitem__(self, values)

        if cls._is_group:
            setattr(cls, attr.attrname,
                    property(getter, setter, None, attr.groupdoc))
            cls._SETATTR_WHITELIST.add(attr.attrname)
        else:
            setattr(cls, attr.singular,
                    property(getter, setter, None, attr.singledoc))
            cls._SETATTR_WHITELIST.add(attr.singular)

    def __setattr__(self, attr, value):
        # `ag.this = 42` calls setattr(ag, 'this', 42)
        if not (attr.startswith('_') or  # 'private' allowed
                attr in self._SETATTR_WHITELIST or  # known attributes allowed
                hasattr(self, attr)):  # preexisting (eg properties) allowed
            raise AttributeError(
                "Cannot set arbitrary attributes to a {}".format(
                    'Group' if self._is_group else 'Component'))
        # if it is, we allow the setattr to proceed by deferring to the super
        # behaviour (ie do it)
        super(_TopologyAttrContainer, self).__setattr__(attr, value)


class _MutableBase(object):
    """
    Base class that merges appropriate :class:`_TopologyAttrContainer` classes.

    Implements :meth:`__new__`. In it the instantiating class is fetched from
    :attr:`~MDAnalysis.core.universe.Universe._classes`. If there is a cache
    miss, a merged class is made
    with a base from :attr:`~MDAnalysis.core.universe.Universe._class_bases`
    and cached.

    The classes themselves are used as the cache dictionary keys for simplcity
    in cache retrieval.

    """
    def __new__(cls, *args, **kwargs):
        # This pre-initialization wrapper must be pretty generic to
        # allow for different initialization schemes of the possible classes.
        # All we really need here is to fish a universe out of the arg list.
        # The AtomGroup cases get priority and are fished out first.
        try:
            u = args[-1].universe
        except (IndexError, AttributeError):
            try:
                # older AtomGroup init method..
                u = args[0][0].universe
            except (TypeError, IndexError, AttributeError):
                from .universe import Universe
                # Let's be generic and get the first argument that's either a
                # Universe, a Group, or a Component, and go from there.
                # This is where the UpdatingAtomGroup args get matched.
                for arg in args+tuple(kwargs.values()):
                    if isinstance(arg, (Universe, GroupBase,
                                        ComponentBase)):
                        u = arg.universe
                        break
                else:
                    raise TypeError("No universe, or universe-containing "
                                   "object passed to the initialization of "
                                    "{}".format(cls.__name__))
        try:
            return object.__new__(u._classes[cls])
        except KeyError:
            # Cache miss. Let's find which kind of class this is and merge.
            try:
                parent_cls = next(u._class_bases[parent]
                                  for parent in cls.mro()
                                  if parent in u._class_bases)
            except StopIteration:
                raise TypeError("Attempted to instantiate class '{}' but "
                                "none of its parents are known to the "
                                "universe. Currently possible parent "
                                "classes are: {}".format(cls.__name__,
                                    str(sorted(u._class_bases.keys()))))
            newcls = u._classes[cls] = parent_cls._mix(cls)
            return object.__new__(newcls)


class _ImmutableBase(object):
    """Class used to shortcut :meth:`__new__` to :meth:`object.__new__`.

    """
    # When mixed via _TopologyAttrContainer._mix this class has MRO priority.
    #  Setting __new__ like this will avoid having to go through the
    #  cache lookup if the class is reused (as in ag._derived_class(...)).
    __new__ = object.__new__



def _only_same_level(function):
    @functools.wraps(function)
    def wrapped(self, other):
        if not isinstance(other, (ComponentBase, GroupBase)):  # sanity check
            raise TypeError("Can't perform '{}' between objects:"
                            " '{}' and '{}'".format(
                                function.__name__,
                                type(self).__name__,
                                type(other).__name__))
        if self.level != other.level:
            raise TypeError("Can't perform '{}' on different level objects"
                            "".format(function.__name__))
        if self.universe is not other.universe:
            raise ValueError(
                "Can't operate on objects from different Universes")
        return function(self, other)
    return wrapped


class GroupBase(_MutableBase):
    """Base class from which a :class:`<~MDAnalysis.core.universe.Universe`\ 's
    Group class is built.

    Instances of :class:`GroupBase` provide the following operations that
    conserve element repetitions and order:

    +-------------------------------+------------+----------------------------+
    | Operation                     | Equivalent | Result                     |
    +===============================+============+============================+
    | ``len(s)``                    |            | number of elements (atoms, |
    |                               |            | residues or segment) in    |
    |                               |            | the group                  |
    +-------------------------------+------------+----------------------------+
    | ``s == t``                    |            | test if ``s`` and ``t``    |
    |                               |            | contain the same elements  |
    |                               |            | in the same order          |
    +-------------------------------+------------+----------------------------+
    | ``x in s``                    |            | test if component ``x`` is |
    |                               |            | part of group ``s``        |
    +-------------------------------+------------+----------------------------+
    | ``s.concatenate(t)``          | ``s + t``  | new Group with elements    |
    |                               |            | from ``s`` and from ``t``  |
    +-------------------------------+------------+----------------------------+
    | ``s.subtract(t)``             |            | new Group with elements    |
    |                               |            | from ``s`` that are not    |
    |                               |            | in ``t``                   |
    +-------------------------------+------------+----------------------------+

    The following operations treat the Group as set. Any result will have any
    duplicate entries removed and the Group will be sorted.

    +-------------------------------+------------+----------------------------+
    | Operation                     | Equivalent | Result                     |
    +===============================+============+============================+
    | ``s.isdisjoint(t)``           |            | ``True`` if ``s`` and      |
    |                               |            | ``t`` do not share         |
    |                               |            | elements                   |
    +-------------------------------+------------+----------------------------+
    | ``s.issubset(t)``             |            | test if all elements of    |
    |                               |            | ``s`` are part of ``t``    |
    +-------------------------------+------------+----------------------------+
    | ``s.is_strict_subset(t)``     |            | test if all elements of    |
    |                               |            | ``s`` are part of ``t``,   |
    |                               |            | and ``s != t``             |
    +-------------------------------+------------+----------------------------+
    | ``s.issuperset(t)``           |            | test if all elements of    |
    |                               |            | ``t`` are part of ``s``    |
    +-------------------------------+------------+----------------------------+
    | ``s.is_strict_superset(t)``   |            | test if all elements of    |
    |                               |            | ``t`` are part of ``s``,   |
    |                               |            | and ``s != t``             |
    +-------------------------------+------------+----------------------------+
    | ``s.union(t)``                | ``s | t``  | new Group with elements    |
    |                               |            | from both ``s`` and ``t``  |
    +-------------------------------+------------+----------------------------+
    | ``s.intersection(t)``         | ``s & t``  | new Group with elements    |
    |                               |            | common to ``s`` and ``t``  |
    +-------------------------------+------------+----------------------------+
    | ``s.difference(t)``           | ``s - t``  | new Group with elements of |
    |                               |            | ``s`` that are not in ``t``|
    +-------------------------------+------------+----------------------------+
    | ``s.symmetric_difference(t)`` | ``s ^ t``  | new Group with elements    |
    |                               |            | that are part of ``s`` or  |
    |                               |            | ``t`` but not both         |
    +-------------------------------+------------+----------------------------+
    """
    def __init__(self, *args):
        try:
            if len(args) == 1:
                # list of atoms/res/segs, old init method
                ix = [at.ix for at in args[0]]
                u = args[0][0].universe
            else:
                # current/new init method
                ix, u = args
        except (AttributeError,  # couldn't find ix/universe
                TypeError):  # couldn't iterate the object we got
            raise TypeError(
                "Can only initialise a Group from an iterable of Atom/Residue/"
                "Segment objects eg: AtomGroup([Atom1, Atom2, Atom3]) "
                "or an iterable of indices and a Universe reference "
                "eg: AtomGroup([0, 5, 7, 8], u).")

        # indices for the objects I hold
        self._ix = np.asarray(ix, dtype=np.intp)
        self._u = u
        self._cache = dict()

    def __hash__(self):
        return hash((self._u, self.__class__, tuple(self.ix.tolist())))

    def __len__(self):
        return len(self.ix)

    def __getitem__(self, item):
        # supports
        # - integer access
        # - boolean slicing
        # - fancy indexing
        # because our _ix attribute is a numpy array
        # it can be sliced by all of these already,
        # so just return ourselves sliced by the item
        if isinstance(item, numbers.Integral):
            return self.level.singular(self.ix[item], self.universe)
        else:
            if isinstance(item, list) and item:  # check for empty list
                # hack to make lists into numpy arrays
                # important for boolean slicing
                item = np.array(item)
            # We specify _derived_class instead of self.__class__ to allow
            # subclasses, such as UpdatingAtomGroup, to control the class
            # resulting from slicing.
            return self._derived_class(self.ix[item], self.universe)

    def __repr__(self):
        name = self.level.name
        return ("<{}Group with {} {}{}>"
                "".format(name.capitalize(), len(self), name,
                "s"[len(self)==1:])) # Shorthand for a conditional plural 's'.

    def __str__(self):
        name = self.level.name
        if len(self) <= 10:
            return '<{}Group {}>'.format(name.capitalize(), repr(list(self)))
        else:
            return '<{}Group {}, ..., {}>'.format(name.capitalize(),
                                                  repr(list(self)[:3])[:-1],
                                                  repr(list(self)[-3:])[1:])

    def __add__(self, other):
        """Concatenate the Group with another Group or Component of the same
        level.

        Parameters
        ----------
        other : Group or Component
            Group or Component with `other.level` same as `self.level`

        Returns
        -------
        Group
            Group with elements of `self` and `other` concatenated

        """
        return self.concatenate(other)

    def __radd__(self, other):
        """Using built-in sum requires supporting 0 + self. If other is
        anything other 0, an exception will be raised.

        Parameters
        ----------
        other : int
            Other should be 0, or else an exception will be raised.

        Returns
        -------
        self
            Group with elements of `self` reproduced

        """
        if other == 0:
            return self._derived_class(self.ix, self.universe)
        else:
            raise TypeError("unsupported operand type(s) for +:"
                            " '{}' and '{}'".format(type(self).__name__,
                                                    type(other).__name__))
    def __sub__(self, other):
        return self.difference(other)

    @_only_same_level
    def __eq__(self, other):
        """Test group equality.

        Two groups are equal if they contain the same indices in
        the same order. Groups that are not at the same level or that belong
        to different :class:`Universes<MDAnalysis.core.universe.Universe>`
        cannot be compared.
        """
        o_ix = other.ix
        return np.array_equal(self.ix, o_ix)

    def __contains__(self, other):
        if not other.level == self.level:
            # maybe raise TypeError instead?
            # eq method raises Error for wrong comparisons
            return False
        return other.ix in self.ix

    def __or__(self, other):
        return self.union(other)

    def __and__(self, other):
        return self.intersection(other)

    def __xor__(self, other):
        return self.symmetric_difference(other)

    @property
    def universe(self):
        """The underlying :class:`~MDAnalysis.core.universe.Universe` the group
        belongs to.
        """
        return self._u

    @property
    def ix(self):
        """Unique indices of the components in the Group.

        - If this Group is an :class:`AtomGroup`, these are the
          indices of the :class:`Atom` instances.
        - If it is a :class:`ResidueGroup`, these are the indices of
          the :class:`Residue` instances.
        - If it is a :class:`SegmentGroup`, these are the indices of
          the :class:`Segment` instances.

        """
        return self._ix

    @property
    def ix_array(self):
        """Unique indices of the components in the Group.

        For a Group, :attr:`ix_array` is the same as :attr:`ix`. This method
        gives a consistent API between components and groups.

        See Also
        --------
        :attr:`ix`
        """
        return self._ix

    @property
    def dimensions(self):
        """Obtain a copy of the dimensions of the currently loaded Timestep"""
        return self.universe.trajectory.ts.dimensions.copy()

    @dimensions.setter
    def dimensions(self, dimensions):
        self.universe.trajectory.ts.dimensions = dimensions

    @property
    @cached('isunique')
    def isunique(self):
        """Boolean indicating whether all components of the group are unique,
        i.e., the group contains no duplicates.

        Examples
        --------

           >>> ag = u.atoms[[2, 1, 2, 2, 1, 0]]
           >>> ag
           <AtomGroup with 6 atoms>
           >>> ag.isunique
           False
           >>> ag2 = ag.unique
           >>> ag2
           <AtomGroup with 3 atoms>
           >>> ag2.isunique
           True


        .. versionadded:: 0.19.0
        """
        if len(self) <= 1:
            return True
        # Fast check for uniqueness
        # 1. get sorted array of component indices:
        s_ix = np.sort(self._ix)
        # 2. If the group's components are unique, no pair of adjacent values in
        #    the sorted indices array are equal. We therefore compute a boolean
        #    mask indicating equality of adjacent sorted indices:
        mask = s_ix[1:] == s_ix[:-1]
        # 3. The group is unique if all elements in the mask are False. We could
        #    return ``not np.any(mask)`` here but using the following is faster:
        return not np.count_nonzero(mask)

    @warn_if_not_unique
    def center(self, weights, pbc=None, compound='group'):
        """Weighted center of (compounds of) the group

        Computes the weighted center of :class:`Atoms<Atom>` in the group.
        Weighted centers per :class:`Residue` or per :class:`Segment` can be
        obtained by setting the `compound` parameter accordingly.

        Parameters
        ----------
        weights : array_like or None
            Weights to be used. Setting `weights=None` is equivalent to passing
            identical weights for all atoms of the group.
        pbc : bool or None, optional
            If ``True`` and `compound` is ``'group'``, move all atoms to the
            primary unit cell before calculation. If ``True`` and `compound` is
            ``'segments'`` or ``'residues'``, the center of each compound will
            be calculated without moving any :class:`Atoms<Atom>` to keep the
            compounds intact. Instead, the resulting position vectors will be
            moved to the primary unit cell after calculation.
        compound : {'group', 'segments', 'residues'}, optional
            If ``'group'``, the weighted center of all atoms in the group will
            be returned as a single position vector. Else, the weighted centers
            of each :class:`Segment` or :class:`Residue` will be returned as an
            array of position vectors, i.e. a 2d array. Note that, in any case,
            *only* the positions of :class:`Atoms<Atom>` *belonging to the
            group* will be taken into account.

        Returns
        -------
        center : numpy.ndarray
            Position vector(s) of the weighted center(s) of the group.
            If `compound` was set to ``'group'``, the output will be a single
            position vector.
            If `compound` was set to ``'segments'`` or ``'residues'``, the
            output will be a 2d array of shape ``(n, 3)`` where ``n`` is the
            number of compounds.

        Examples
        --------

        To find the center of charge of a given :class:`AtomGroup`::

            >>> sel = u.select_atoms('prop mass > 4.0')
            >>> sel.center(sel.charges)

        To find the centers of mass per residue of all CA :class:`Atoms<Atom>`::

            >>> sel = u.select_atoms('name CA')
            >>> sel.center(sel.masses, compound='residues')

        Notes
        -----
        If the :class:`MDAnalysis.core.flags` flag *use_pbc* is set to
        ``True`` then the `pbc` keyword is used by default.


        .. versionchanged:: 0.19.0 Added `compound` parameter
        """

        if pbc is None:
            pbc = flags['use_pbc']

        atoms = self.atoms

        # enforce calculations in double precision:
        dtype = np.float64

        if compound.lower() == 'group':
            if pbc:
                coords = atoms.pack_into_box(inplace=False)
            else:
                coords = atoms.positions
            # promote coords or weights to dtype if required:
            if weights is None:
                coords = coords.astype(dtype, copy=False)
            else:
                weights = weights.astype(dtype, copy=False)
            return np.average(coords, weights=weights, axis=0)
        elif compound.lower() == 'residues':
            compound_indices = atoms.resindices
            n_compounds = atoms.n_residues
        elif compound.lower() == 'segments':
            compound_indices = atoms.segindices
            n_compounds = atoms.n_segments
        else:
            raise ValueError("Unrecognized compound definition: {}\nPlease use"
                             " one of 'group', 'residues', or 'segments'."
                             "".format(compound))

        # Sort positions and weights by compound index and promote to dtype if
        # required:
        sort_indices = np.argsort(compound_indices)
        compound_indices = compound_indices[sort_indices]
        coords = atoms.positions[sort_indices]
        if weights is None:
            coords = coords.astype(dtype, copy=False)
        else:
            weights = weights.astype(dtype, copy=False)
            weights = weights[sort_indices]
        # Allocate output array:
        centers = np.zeros((n_compounds, 3), dtype=dtype)
        # Get sizes of compounds:
        unique_compound_indices, compound_sizes = np.unique(compound_indices,
                                                            return_counts=True)
        unique_compound_sizes = unique_int_1d(compound_sizes)
        # Compute centers per compound for each compound size:
        for compound_size in unique_compound_sizes:
            compound_mask = compound_sizes == compound_size
            _compound_indices = unique_compound_indices[compound_mask]
            atoms_mask = np.in1d(compound_indices, _compound_indices)
            _coords = coords[atoms_mask].reshape((-1, compound_size, 3))
            if weights is None:
                _centers = _coords.mean(axis=1)
            else:
                _weights = weights[atoms_mask].reshape((-1, compound_size))
                _centers = (_coords * _weights[:, :, None]).sum(axis=1)
                _centers /= _weights.sum(axis=1)[:, None]
            centers[compound_mask] = _centers
        if pbc:
            centers = distances.apply_PBC(centers, atoms.dimensions)
        return centers

    @warn_if_not_unique
    def center_of_geometry(self, pbc=None, compound='group'):
        """Center of geometry (also known as centroid) of the group.

        Computes the center of geometry (a.k.a. centroid) of
        :class:`Atoms<Atom>` in the group. Centers of geometry per
        :class:`Residue` or per :class:`Segment` can be obtained by setting the
        `compound` parameter accordingly.

        Parameters
        ----------
        pbc : bool or None, optional
            If ``True`` and `compound` is ``'group'``, move all atoms to the
            primary unit cell before calculation. If ``True`` and `compound` is
            ``'segments'`` or ``'residues'``, the center of each compound will
            be calculated without moving any :class:`Atoms<Atom>` to keep the
            compounds intact. Instead, the resulting position vectors will be
            moved to the primary unit cell after calculation.
        compound : {'group', 'segments', 'residues'}, optional
            If ``'group'``, the center of geometry of all :class:`Atoms<Atom>`
            in the group will be returned as a single position vector. Else, the
            centers of geometry of each :class:`Segment` or :class:`Residue`
            will be returned as an array of position vectors, i.e. a 2d array.
            Note that, in any case, *only* the positions of :class:`Atoms<Atom>`
            *belonging to the group* will be taken into account.

        Returns
        -------
        center : numpy.ndarray
            Position vector(s) of the geometric center(s) of the group.
            If `compound` was set to ``'group'``, the output will be a single
            position vector.
            If `compound` was set to ``'segments'`` or ``'residues'``, the
            output will be a 2d array of shape ``(n, 3)`` where ``n`` is the
            number of compounds.

        Notes
        -----
        If the :class:`MDAnalysis.core.flags` flag *use_pbc* is set to
        ``True`` then the `pbc` keyword is used by default.


        .. versionchanged:: 0.8 Added `pbc` keyword
        .. versionchanged:: 0.19.0 Added `compound` parameter
        """
        return self.center(None, pbc=pbc, compound=compound)

    centroid = center_of_geometry

    def bbox(self, **kwargs):
        """Return the bounding box of the selection.

        The lengths A,B,C of the orthorhombic enclosing box are ::

          L = AtomGroup.bbox()
          A,B,C = L[1] - L[0]

        Parameters
        ----------
        pbc : bool, optional
            If ``True``, move all :class:`Atoms<Atom>` to the primary unit cell
            before calculation. [``False``]

        Returns
        -------
         corners : numpy.ndarray
            2x3 array giving corners of bounding box as
            ``[[xmin, ymin, zmin], [xmax, ymax, zmax]]``.

        Note
        ----
        The :class:`MDAnalysis.core.flags` flag *use_pbc* when set to
        ``True`` allows the *pbc* flag to be used by default.


        .. versionadded:: 0.7.2
        .. versionchanged:: 0.8 Added *pbc* keyword
        """
        atomgroup = self.atoms
        pbc = kwargs.pop('pbc', flags['use_pbc'])

        if pbc:
            x = atomgroup.pack_into_box(inplace=False)
        else:
            x = atomgroup.positions

        return np.array([x.min(axis=0), x.max(axis=0)])

    def bsphere(self, **kwargs):
        """Return the bounding sphere of the selection.

        The sphere is calculated relative to the
        :meth:`center of geometry<center_of_geometry>`.

        Parameters
        ----------
        pbc : bool, optional
            If ``True``, move all atoms to the primary unit cell before
            calculation. [``False``]

        Returns
        -------
        R : float
            Radius of the bounding sphere.
        center : numpy.ndarray
            Coordinates of the sphere center as ``[xcen, ycen, zcen]``.

        Note
        ----
        The :class:`MDAnalysis.core.flags` flag *use_pbc* when set to
        ``True`` allows the *pbc* flag to be used by default.


        .. versionadded:: 0.7.3
        .. versionchanged:: 0.8 Added *pbc* keyword
        """
        atomgroup = self.atoms.unique
        pbc = kwargs.pop('pbc', flags['use_pbc'])

        if pbc:
            x = atomgroup.pack_into_box(inplace=False)
            centroid = atomgroup.center_of_geometry(pbc=True)
        else:
            x = atomgroup.positions
            centroid = atomgroup.center_of_geometry(pbc=False)

        R = np.sqrt(np.max(np.sum(np.square(x - centroid), axis=1)))

        return R, centroid

    def transform(self, M):
        r"""Apply homogenous transformation matrix `M` to the coordinates.

        :class:`Atom` coordinates are rotated and translated in-place.

        Parameters
        ----------
        M : array_like
            4x4 matrix with the rotation in ``R = M[:3, :3]`` and the
            translation in ``t = M[:3, 3]``.

        Returns
        -------
        self

        See Also
        --------
        MDAnalysis.lib.transformations : module of all coordinate transforms

        Notes
        -----
        The rotation :math:`\mathsf{R}` is about the origin and is applied
        before the translation :math:`\mathbf{t}`:

        .. math::

           \mathbf{x}' = \mathsf{R}\mathbf{x} + \mathbf{t}

        """
        M = np.asarray(M)
        R = M[:3, :3]
        t = M[:3, 3]
        return self.rotate(R, [0, 0, 0]).translate(t)

    def translate(self, t):
        r"""Apply translation vector `t` to the selection's coordinates.

        :class:`Atom` coordinates are translated in-place.

        Parameters
        ----------
        t : array_like
            vector to translate coordinates with

        Returns
        -------
        self

        See Also
        --------
        MDAnalysis.lib.transformations : module of all coordinate transforms

        Notes
        -----
        The method applies a translation to the :class:`AtomGroup`
        from current coordinates :math:`\mathbf{x}` to new coordinates
        :math:`\mathbf{x}'`:

        .. math::

            \mathbf{x}' = \mathbf{x} + \mathbf{t}

        """
        atomgroup = self.atoms.unique
        vector = np.asarray(t)
        # changes the coordinates in place
        atomgroup.universe.trajectory.ts.positions[atomgroup.indices] += vector
        return self

    def rotate(self, R, point=(0, 0, 0)):
        r"""Apply a rotation matrix `R` to the selection's coordinates.
        :math:`\mathsf{R}` is a 3x3 orthogonal matrix that transforms a vector
        :math:`\mathbf{x} \rightarrow \mathbf{x}'`:

        .. math::

            \mathbf{x}' = \mathsf{R}\mathbf{x}

        :class:`Atom` coordinates are rotated in-place.

        Parameters
        ----------
        R : array_like
            3x3 rotation matrix
        point : array_like, optional
            Center of rotation

        Returns
        -------
        self

        Notes
        -----
        By default, rotates about the origin ``point=(0, 0, 0)``. To rotate
        a group ``g`` around its center of geometry, use
        ``g.rotate(R, point=g.center_of_geometry())``.

        See Also
        --------
        rotateby : rotate around given axis and angle
        MDAnalysis.lib.transformations : module of all coordinate transforms

        """
        R = np.asarray(R)
        point = np.asarray(point)

        # changes the coordinates (in place)
        atomgroup = self.atoms.unique
        require_translation = bool(np.count_nonzero(point))
        if require_translation:
            atomgroup.translate(-point)
        x = atomgroup.universe.trajectory.ts.positions
        idx = atomgroup.indices
        x[idx] = np.dot(x[idx], R.T)
        if require_translation:
            atomgroup.translate(point)

        return self

    def rotateby(self, angle, axis, point=None):
        r"""Apply a rotation to the selection's coordinates.

        Parameters
        ----------
        angle : float
            Rotation angle in degrees.
        axis : array_like
            Rotation axis vector.
        point : array_like, optional
            Center of rotation. If ``None`` then the center of geometry of this
            group is used.

        Returns
        -------
        self

        Notes
        -----
        The transformation from current coordinates :math:`\mathbf{x}`
        to new coordinates :math:`\mathbf{x}'` is

        .. math::

          \mathbf{x}' = \mathsf{R}\,(\mathbf{x}-\mathbf{p}) + \mathbf{p}

        where :math:`\mathsf{R}` is the rotation by `angle` around the
        `axis` going through `point` :math:`\mathbf{p}`.

        See Also
        --------
        MDAnalysis.lib.transformations.rotation_matrix :
            calculate :math:`\mathsf{R}`

        """
        alpha = np.radians(angle)
        axis = np.asarray(axis)
        if point is None:
            point = self.center_of_geometry()
        point = np.asarray(point)
        M = transformations.rotation_matrix(alpha, axis, point=point)
        return self.transform(M)

    def pack_into_box(self, box=None, inplace=True):
        r"""Shift all :class:`Atoms<Atom>` in this group to the primary unit
        cell.

        Parameters
        ----------
        box : array_like
            Box dimensions, can be either orthogonal or triclinic information.
            Cell dimensions must be in an identical to format to those returned
            by :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`,
            ``[lx, ly, lz, alpha, beta, gamma]``. If ``None``, uses these
            timestep dimensions.
        inplace : bool
            ``True`` to change coordinates in place.

        Returns
        -------
        coords : numpy.ndarray
            Shifted atom coordinates.

        Notes
        -----
        All atoms will be moved so that they lie between 0 and boxlength
        :math:`L_i` in all dimensions, i.e. the lower left corner of the
        simulation box is taken to be at (0,0,0):

        .. math::

           x_i' = x_i - \left\lfloor\frac{x_i}{L_i}\right\rfloor

        The default is to take unit cell information from the underlying
        :class:`~MDAnalysis.coordinates.base.Timestep` instance. The optional
        argument `box` can be used to provide alternative unit cell information
        (in the MDAnalysis standard format
        ``[Lx, Ly, Lz, alpha, beta, gamma]``).

        Works with either orthogonal or triclinic box types.

        .. note::
            :meth:`AtomGroup.pack_into_box` is currently faster than
            :meth:`ResidueGroup.pack_into_box` or
            :meth:`SegmentGroup.pack_into_box`.


        .. versionadded:: 0.8
        """
        # Try and auto detect box dimensions:
        if box is None:
            box = self.dimensions

        # For a vector representation, the diagonal must not be zero, for a
        # [x, y, z, alpha, beta, gamma] representation, no element must be zero:
        if (box.shape == (3, 3) and (box.diagonal() == 0.0).any()) \
            or (box == 0).any():
                raise ValueError("One or more box dimensions is zero."
                                 "  You can specify a box with 'box ='")
        # no matter what kind of group we have, we need to work on its atoms:
        ag = self.atoms
        # If self is unique, so are its atoms. Thus, if the group is unique, we
        # can use its atoms as is.
        if self.isunique:
            packed_coords = distances.apply_PBC(ag.positions, box)
            if inplace:
                ag.universe.trajectory.ts.positions[ag.ix] = packed_coords
                return ag.universe.trajectory.ts.positions[ag.ix]
            return packed_coords
        else:
            unique_ag = ag.unique
            unique_ix = unique_ag.ix
            restore_mask = ag._unique_restore_mask
            coords = ag.universe.trajectory.ts.positions[unique_ix]
            packed_coords = distances.apply_PBC(coords, box)
            if inplace:
                ag.universe.trajectory.ts.positions[unique_ix] = packed_coords
                return ag.universe.trajectory.ts.positions[unique_ix[restore_mask]]
            return packed_coords[restore_mask]

    def wrap(self, compound="atoms", center="com", box=None):
        """Shift the contents of this group back into the unit cell.

        This is a more powerful version of :meth:`pack_into_box`, allowing
        groups of :class:`Atoms<Atom>` to be kept together through the process.

        Parameters
        ----------
        compound : {'atoms', 'group', 'residues', 'segments', 'fragments'}
            The group which will be kept together through the shifting process.
        center : {'com', 'cog'}
            How to define the center of a given group of atoms.
        box : array_like
            Box dimensions, can be either orthogonal or triclinic information.
            Cell dimensions must be in an identical to format to those returned
            by :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`,
            ``[lx, ly, lz, alpha, beta, gamma]``. If ``None``, uses these
            timestep dimensions.

        Notes
        -----
        When specifying a `compound`, the translation is calculated based on
        each compound. The same translation is applied to all atoms
        within this compound, meaning it will not be broken by the shift.
        This might however mean that not all atoms from the compound are
        inside the unit cell, but rather the center of the compound is.

        `center` allows the definition of the center of each group to be
        specified. This can be either ``'com'`` for center of mass, or ``'cog'``
        for center of geometry.

        `box` allows a unit cell to be given for the transformation. If not
        specified, the :attr:`~MDAnalysis.coordinates.base.Timestep.dimensions`
        information from the current
        :class:`~MDAnalysis.coordinates.base.Timestep` will be used.

        .. note::
           :meth:`wrap` with all default keywords is identical to
           :meth:`pack_into_box`


        .. versionadded:: 0.9.2
        """
        atomgroup = self.atoms.unique
        if compound.lower() == "atoms":
            return atomgroup.pack_into_box(box=box)

        if compound.lower() == 'group':
            objects = [atomgroup.atoms]
        elif compound.lower() == 'residues':
            objects = atomgroup.residues
        elif compound.lower() == 'segments':
            objects = atomgroup.segments
        elif compound.lower() == 'fragments':
            objects = atomgroup.fragments
        else:
            raise ValueError("Unrecognized compound definition: {0}"
                             "Please use one of 'group' 'residues' 'segments'"
                             "or 'fragments'".format(compound))

        # TODO: ADD TRY-EXCEPT FOR MASSES PRESENCE
        if center.lower() in ('com', 'centerofmass'):
            centers = np.vstack([o.atoms.center_of_mass() for o in objects])
        elif center.lower() in ('cog', 'centroid', 'centerofgeometry'):
            centers = np.vstack([o.atoms.center_of_geometry() for o in objects])
        else:
            raise ValueError("Unrecognized center definition: {0}"
                             "Please use one of 'com' or 'cog'".format(center))
        centers = centers.astype(np.float32)

        if box is None:
            box = atomgroup.dimensions

        # calculate shift per object center
        dests = distances.apply_PBC(centers, box=box)
        shifts = dests - centers

        for o, s in zip(objects, shifts):
            # Save some needless shifts
            if not all(s == 0.0):
                o.atoms.translate(s)

    def copy(self):
        """Get another group identical to this one.


        .. versionadded:: 0.19.0
        """
        group = self[:]
        # Try to fill the copied group's uniqueness caches:
        try:
            group._cache['isunique'] = self._cache['isunique']
            if group._cache['isunique']:
                group._cache['unique'] = group
        except KeyError:
            pass
        return group

    def groupby(self, topattrs):
        """Group together items in this group according to values of *topattr*

        Parameters
        ----------
        topattrs: str or list
           One or more topology attributes to group components by.
           Single arguments are passed as a string. Multiple arguments
           are passed as a list.

        Returns
        -------
        dict
            Unique values of the multiple combinations of topology attributes
            as keys, Groups as values.

        Example
        -------
        To group atoms with the same mass together:

        >>> ag.groupby('masses')
        {12.010999999999999: <AtomGroup with 462 atoms>,
         14.007: <AtomGroup with 116 atoms>,
         15.999000000000001: <AtomGroup with 134 atoms>}

        To group atoms with the same residue name and mass together:

          >>> ag.groupby(['resnames', 'masses'])
          {('ALA', 1.008): <AtomGroup with 95 atoms>,
           ('ALA', 12.011): <AtomGroup with 57 atoms>,
           ('ALA', 14.007): <AtomGroup with 19 atoms>,
           ('ALA', 15.999): <AtomGroup with 19 atoms>},
           ('ARG', 1.008): <AtomGroup with 169 atoms>,
           ...}

          >>> ag.groupby(['resnames', 'masses'])('ALA', 15.999)
           <AtomGroup with 19 atoms>


        .. versionadded:: 0.16.0
        .. versionchanged:: 0.18.0 The function accepts multiple attributes
        """

        res = dict()

        if isinstance(topattrs, (string_types, bytes)):
            attr=topattrs
            if isinstance(topattrs, bytes):
                attr = topattrs.decode('utf-8')
            ta = getattr(self, attr)

            return {i: self[ta == i] for i in set(ta)}

        else:
            attr = topattrs[0]
            ta = getattr(self, attr)
            for i in set(ta):
                if len(topattrs) == 1:
                    res[i] = self[ta == i]
                else:
                    res[i] = self[ta == i].groupby(topattrs[1:])

            return util.flatten_dict(res)

    @_only_same_level
    def concatenate(self, other):
        """Concatenate with another Group or Component of the same level.

        Duplicate entries and original order is preserved. It is synomymous to
        the `+` operator.

        Parameters
        ----------
        other : Group or Component
            Group or Component with `other.level` same as `self.level`

        Returns
        -------
        Group
            Group with elements of `self` and `other` concatenated

        Example
        -------
        The order of the original contents (including duplicates)
        are preserved when performing a concatenation.

        >>> ag1 = u.select_atoms('name O')
        >>> ag2 = u.select_atoms('name N')
        >>> ag3 = ag1 + ag2  # or ag1.concatenate(ag2)
        >>> ag3[:3].names
        array(['O', 'O', 'O'], dtype=object)
        >>> ag3[-3:].names
        array(['N', 'N', 'N'], dtype=object)


        .. versionadded:: 0.16.0
        """
        o_ix = other.ix_array
        return self._derived_class(np.concatenate([self.ix, o_ix]),
                                   self.universe)

    @_only_same_level
    def union(self, other):
        """Group of elements either in this Group or another

        On the contrary to concatenation, this method sort the elements and
        removes duplicate ones. It is synomymous to the `|` operator.

        Parameters
        ----------
        other : Group or Component
            Group or Component with `other.level` same as `self.level`

        Returns
        -------
        Group
            Group with the combined elements of `self` and `other`, without
            duplicate elements

        Example
        -------
        In contrast to :meth:`concatenate`, any duplicates are dropped
        and the result is sorted.

        >>> ag1 = u.select_atoms('name O')
        >>> ag2 = u.select_atoms('name N')
        >>> ag3 = ag1 | ag2  # or ag1.union(ag2)
        >>> ag3[:3].names
        array(['N', 'O', 'N'], dtype=object)

        See Also
        --------
        concatenate, intersection


        .. versionadded:: 0.16
        """
        o_ix = other.ix_array
        return self._derived_class(np.union1d(self.ix, o_ix), self.universe)

    @_only_same_level
    def intersection(self, other):
        """Group of elements which are in both this Group and another

        This method removes duplicate elements and sorts the result. It is
        synomymous to the `&` operator.

        Parameters
        ----------
        other : Group or Component
            Group or Component with `other.level` same as `self.level`

        Returns
        -------
        Group
            Group with the common elements of `self` and `other`, without
            duplicate elements

        Example
        -------
        Intersections can be used when the select atoms string would
        become too complicated.  For example to find the water atoms
        which are within 4.0A of two segments:

        >>> water = u.select_atoms('resname SOL')
        >>> shell1 = water.select_atoms('around 4.0 segid 1')
        >>> shell2 = water.select_atoms('around 4.0 segid 2')
        >>> common = shell1 & shell2  # or shell1.intersection(shell2)

        See Also
        --------
        union


        .. versionadded:: 0.16
        """
        o_ix = other.ix_array
        return self._derived_class(np.intersect1d(self.ix, o_ix), self.universe)

    @_only_same_level
    def subtract(self, other):
        """Group with elements from this Group that don't appear in other

        The original order of this group is kept, as well as any duplicate
        elements. If an element of this Group is duplicated and appears in
        the other Group or Component, then all the occurences of that element
        are removed from the returned Group.

        Parameters
        ----------
        other : Group or Component
            Group or Component with `other.level` same as `self.level`

        Returns
        -------
        Group
            Group with the elements of `self` that are not in  `other`,
            conserves order and duplicates.

        Example
        -------
        Unlike :meth:`difference` this method will not sort or remove
        duplicates.

        >>> ag1 = u.atoms[[3, 3, 2, 2, 1, 1]]
        >>> ag2 = u.atoms[2]
        >>> ag3 = ag1 - ag2  # or ag1.subtract(ag2)
        >>> ag1.indices
        array([3, 3, 1, 1])

        See Also
        --------
        concatenate, difference


        .. versionadded:: 0.16
        """
        o_ix = other.ix_array
        in_other = np.in1d(self.ix, o_ix)  # mask of in self.ix AND other
        return self[~in_other]  # ie inverse of previous mask

    @_only_same_level
    def difference(self, other):
        """Elements from this Group that do not appear in another

        This method removes duplicate elements and sorts the result. As such,
        it is different from :meth:`subtract`. :meth:`difference` is synomymous
        to the `-` operator.

        Parameters
        ----------
        other : Group or Component
            Group or Component with `other.level` same as `self.level`

        Returns
        -------
        Group
            Group with the elements of `self` that are not in  `other`, without
            duplicate elements

        See Also
        --------
        subtract, symmetric_difference


        .. versionadded:: 0.16
        """
        o_ix = other.ix_array
        return self._derived_class(np.setdiff1d(self._ix, o_ix), self._u)

    @_only_same_level
    def symmetric_difference(self, other):
        """Group of elements which are only in one of this Group or another

        This method removes duplicate elements and the result is sorted. It is
        synomym to the `^` operator.

        Parameters
        ----------
        other : Group or Component
            Group or Component with `other.level` same as `self.level`

        Returns
        -------
        Group
            Group with the elements that are in `self` or in `other` but not in
            both, without duplicate elements

        Example
        -------

        >>> ag1 = u.atoms[[0, 1, 5, 3, 3, 2]]
        >>> ag2 = u.atoms[[4, 4, 6, 2, 3, 5]]
        >>> ag3 = ag1 ^ ag2  # or ag1.symmetric_difference(ag2)
        >>> ag3.indices  # 0 and 1 are only in ag1, 4 and 6 are only in ag2
        [0, 1, 4, 6]

        See Also
        --------
        difference


        .. versionadded:: 0.16
        """
        o_ix = other.ix_array
        return self._derived_class(np.setxor1d(self._ix, o_ix), self._u)

    def isdisjoint(self, other):
        """If the Group has no elements in common with the other Group

        Parameters
        ----------
        other : Group or Component
            Group or Component with `other.level` same as `self.level`

        Returns
        -------
        bool
            ``True`` if the two Groups do not have common elements


        .. versionadded:: 0.16
        """
        return len(self.intersection(other)) == 0

    @_only_same_level
    def issubset(self, other):
        """If all elements of this Group are part of another Group

        Note that an empty group is a subset of any group of the same level.

        Parameters
        ----------
        other : Group or Component
            Group or Component with `other.level` same as `self.level`

        Returns
        -------
        bool
            ``True`` if this Group is a subset of the other one


        .. versionadded:: 0.16
        """
        o_ix = set(other.ix_array)
        s_ix = set(self.ix)
        return s_ix.issubset(o_ix)

    def is_strict_subset(self, other):
        """If this Group is a subset of another Group but not identical

        Parameters
        ----------
        other : Group or Component
            Group or Component with `other.level` same as `self.level`

        Returns
        -------
        bool
            ``True`` if this Group is a strict subset of the other one


        .. versionadded:: 0.16
        """
        return self.issubset(other) and not self == other

    @_only_same_level
    def issuperset(self, other):
        """If all elements of another Group are part of this Group

        Parameters
        ----------
        other : Group or Component
            Group or Component with `other.level` same as `self.level`

        Returns
        -------
        bool
            ``True`` if this Group is a subset of the other one


        .. versionadded:: 0.16
        """
        o_ix = set(other.ix_array)
        s_ix = set(self.ix)
        return s_ix.issuperset(o_ix)

    def is_strict_superset(self, other):
        """If this Group is a superset of another Group but not identical

        Parameters
        ----------
        other : Group or Component
            Group or Component with `other.level` same as `self.level`

        Returns
        -------
        bool
            ``True`` if this Group is a strict superset of the other one


        .. versionadded:: 0.16
        """
        return self.issuperset(other) and not self == other


class AtomGroup(GroupBase):
    """An ordered array of atoms.

    Can be initiated from an iterable of :class:`Atoms<Atom>`::

        ag = AtomGroup([Atom1, Atom2, Atom3])

    Or from providing a list of indices and the
    :class:`~MDAnalysis.core.universe.Universe` it should belong to::

        ag = AtomGroup([72, 14, 25], u)

    Alternatively, an :class:`AtomGroup` is generated by indexing/slicing
    another :class:`AtomGroup`, such as the group of all :class:`Atoms<Atom>` in
    the :class:`~MDAnalysis.core.universe.Universe` at
    :attr:`MDAnalysis.core.universe.Universe.atoms`.

    An :class:`AtomGroup` can be indexed and sliced like a list::

        ag[0], ag[-1]

    will return the first and the last :class:`Atom` in the group whereas the
    slice::

        ag[0:6:2]

    returns an :class:`AtomGroup` of every second element, corresponding to
    indices 0, 2, and 4.

    It also supports "advanced slicing" when the argument is a
    :class:`numpy.ndarray` or a :class:`list`::

        aslice = [0, 3, -1, 10, 3]
        ag[aslice]

    will return a new :class:`AtomGroup` of :class:`Atoms<Atom>` with those
    indices in the old :class:`AtomGroup`.

    Finally, :class:`AtomGroups<AtomGroup>` can be created from a selection.
    See :meth:`select_atoms`.

    .. note::

        :class:`AtomGroups<AtomGroup>` originating from a selection are sorted
        and duplicate elements are removed. This is not true for
        :class:`AtomGroups<AtomGroup>` produced by slicing. Thus, slicing can be
        used when the order of atoms is crucial (for instance, in order to
        define angles or dihedrals).

    :class:`AtomGroups<AtomGroup>` can be compared and combined using group
    operators. For instance, :class:`AtomGroups<AtomGroup>` can be concatenated
    using `+` or :meth:`concatenate`::

        ag_concat = ag1 + ag2  # or ag_concat = ag1.concatenate(ag2)

    When groups are concatenated, the order of the :class:`Atoms<Atom>` is
    conserved. If :class:`Atoms<Atom>` appear several times in one of the
    groups, the duplicates are kept in the resulting group. On the contrary to
    :meth:`concatenate`, :meth:`union` treats the :class:`AtomGroups<AtomGroup>`
    as sets so that duplicates are removed from the resulting group, and
    :class:`Atoms<Atom>` are ordered. The `|` operator is synomymous to
    :meth:`union`::

        ag_union = ag1 | ag2  # or ag_union = ag1.union(ag2)

    The opposite operation to :meth:`concatenate` is :meth:`subtract`. This
    method creates a new group with all the :class:`Atoms<Atom>` of the group
    that are not in a given other group; the order of the :class:`Atoms<Atom>`
    is kept, and so are duplicates. :meth:`difference` is the set version of
    :meth:`subtract`. The resulting group is sorted and deduplicated.

    All set methods are listed in the table below. These methods treat the
    groups as sorted and deduplicated sets of :class:`Atoms<Atom>`.

    +-------------------------------+------------+----------------------------+
    | Operation                     | Equivalent | Result                     |
    +===============================+============+============================+
    | ``s.isdisjoint(t)``           |            | ``True`` if ``s`` and      |
    |                               |            | ``t`` do not share         |
    |                               |            | elements                   |
    +-------------------------------+------------+----------------------------+
    | ``s.issubset(t)``             |            | test if all elements of    |
    |                               |            | ``s`` are part of ``t``    |
    +-------------------------------+------------+----------------------------+
    | ``s.is_strict_subset(t)``     |            | test if all elements of    |
    |                               |            | ``s`` are part of ``t``,   |
    |                               |            | and ``s != t``             |
    +-------------------------------+------------+----------------------------+
    | ``s.issuperset(t)``           |            | test if all elements of    |
    |                               |            | ``t`` are part of ``s``    |
    +-------------------------------+------------+----------------------------+
    | ``s.is_strict_superset(t)``   |            | test if all elements of    |
    |                               |            | ``t`` are part of ``s``,   |
    |                               |            | and ``s != t``             |
    +-------------------------------+------------+----------------------------+
    | ``s.union(t)``                | ``s | t``  | new Group with elements    |
    |                               |            | from both ``s`` and ``t``  |
    +-------------------------------+------------+----------------------------+
    | ``s.intersection(t)``         | ``s & t``  | new Group with elements    |
    |                               |            | common to ``s`` and ``t``  |
    +-------------------------------+------------+----------------------------+
    | ``s.difference(t)``           | ``s - t``  | new Group with elements of |
    |                               |            | ``s`` that are not in ``t``|
    +-------------------------------+------------+----------------------------+
    | ``s.symmetric_difference(t)`` | ``s ^ t``  | new Group with elements    |
    |                               |            | that are part of ``s`` or  |
    |                               |            | ``t`` but not both         |
    +-------------------------------+------------+----------------------------+

    The following methods keep the order of the atoms as well as duplicates.

    +-------------------------------+------------+----------------------------+
    | Operation                     | Equivalent | Result                     |
    +===============================+============+============================+
    | ``len(s)``                    |            | number of elements (atoms, |
    |                               |            | residues or segment) in    |
    |                               |            | the group                  |
    +-------------------------------+------------+----------------------------+
    | ``s == t``                    |            | test if ``s`` and ``t``    |
    |                               |            | contain the same elements  |
    |                               |            | in the same order          |
    +-------------------------------+------------+----------------------------+
    | ``s.concatenate(t)``          | ``s + t``  | new Group with elements    |
    |                               |            | from ``s`` and from ``t``  |
    +-------------------------------+------------+----------------------------+
    | ``s.subtract(t)``             |            | new Group with elements    |
    |                               |            | from ``s`` that are not    |
    |                               |            | in ``t``                   |
    +-------------------------------+------------+----------------------------+

    The `in` operator allows to test if an :class:`Atom` is in the
    :class:`AtomGroup`.

    :class:`AtomGroup` instances are always bound to a
    :class:`MDAnalysis.core.universe.Universe`. They cannot exist in isolation.

    .. rubric:: Deprecated functionality

    *Instant selectors* will be removed in the 1.0 release.  See issue `#1377
    <https://github.com/MDAnalysis/mdanalysis/issues/1377>`_ for more details.

    :class:`Atoms<Atom>` can also be accessed in a Pythonic fashion by using the
    :class:`Atom` name as an attribute. For instance, ::

        ag.CA

    will provide an :class:`AtomGroup` of all CA :class:`Atoms<Atom>` in the
    group. These *instant selector* attributes are auto-generated for
    each atom name encountered in the group.

    Notes
    -----
    The name-attribute instant selector access to :class:`Atoms<Atom>` is mainly
    meant for quick interactive work. Thus it either returns a single
    :class:`Atom` if there is only one matching :class:`Atom`, *or* a
    new :class:`AtomGroup` for multiple matches.  This makes it difficult to use
    the feature consistently in scripts.


    See Also
    --------
    :class:`MDAnalysis.core.universe.Universe`


    .. deprecated:: 0.16.2
       *Instant selectors* of :class:`AtomGroup` will be removed in the 1.0
       release. See :ref:`Instant selectors <instance-selectors>` for details
       and alternatives.
    """
    def __getitem__(self, item):
        # DEPRECATED in 0.16.2
        # REMOVE in 1.0
        #
        # u.atoms['HT1'] access, otherwise default
        if isinstance(item, string_types):
            try:
                return self._get_named_atom(item)
            except (AttributeError, selection.SelectionError):
                pass
        return super(AtomGroup, self).__getitem__(item)

    def __getattr__(self, attr):
        # DEPRECATED in 0.16.2
        # REMOVE in 1.0
        #
        # is this a known attribute failure?
        if attr in ('fragments',):  # TODO: Generalise this to cover many attributes
            # eg:
            # if attr in _ATTR_ERRORS:
            # raise NDE(_ATTR_ERRORS[attr])
            raise NoDataError("AtomGroup has no fragments; this requires Bonds")
        elif hasattr(self.universe._topology, 'names'):
            # Ugly hack to make multiple __getattr__s work
            try:
                return self._get_named_atom(attr)
            except selection.SelectionError:
                pass
        raise AttributeError("{cls} has no attribute {attr}".format(
            cls=self.__class__.__name__, attr=attr))

    def __reduce__(self):
        return (_unpickle, (self.universe.anchor_name, self.ix))

    @property
    def atoms(self):
        """The :class:`AtomGroup` itself.

        See Also
        --------
        copy : return a true copy of the :class:`AtomGroup`


        .. versionchanged:: 0.19.0
           In previous versions, this returned a copy, but now
           the :class:`AtomGroup` itself is returned. This should
           not affect any code but only speed up calculations.

        """
        return self

    @property
    def n_atoms(self):
        """Number of atoms in the :class:`AtomGroup`.

        Equivalent to ``len(self)``.
        """
        return len(self)

    @property
    def residues(self):
        """A sorted :class:`ResidueGroup` of the unique
        :class:`Residues<Residue>` present in the :class:`AtomGroup`.
        """
        rg = self.universe.residues[unique_int_1d(self.resindices)]
        rg._cache['isunique'] = True
        rg._cache['unique'] = rg
        return rg

    @residues.setter
    def residues(self, new):
        # Can set with Res, ResGroup or list/tuple of Res
        if isinstance(new, Residue):
            r_ix = itertools.cycle((new.resindex,))
        elif isinstance(new, ResidueGroup):
            r_ix = new.resindices
        else:
            try:
                r_ix = [r.resindex for r in new]
            except AttributeError:
                raise TypeError("Can only set AtomGroup residues to Residue "
                                "or ResidueGroup not {}".format(
                                    ', '.join(type(r) for r in new
                                              if not isinstance(r, Residue))
                                ))
        if not isinstance(r_ix, itertools.cycle) and len(r_ix) != len(self):
            raise ValueError("Incorrect size: {} for AtomGroup of size: {}"
                             "".format(len(new), len(self)))
        # Optimisation TODO:
        # This currently rebuilds the tt len(self) times
        # Ideally all changes would happen and *afterwards* tables are built
        # Alternatively, if the changes didn't rebuild table, this list
        # comprehension isn't terrible.
        for at, r in zip(self, r_ix):
            self.universe._topology.tt.move_atom(at.ix, r)

    @property
    def n_residues(self):
        """Number of unique :class:`Residues<Residue>` present in the
        :class:`AtomGroup`.

        Equivalent to ``len(self.residues)``.

        """
        return len(self.residues)

    @property
    def segments(self):
        """A sorted :class:`SegmentGroup` of the unique segments present in the
        :class:`AtomGroup`.
        """
        sg = self.universe.segments[unique_int_1d(self.segindices)]
        sg._cache['isunique'] = True
        sg._cache['unique'] = sg
        return sg

    @segments.setter
    def segments(self, new):
        raise NotImplementedError("Cannot assign Segments to AtomGroup. "
                                  "Segments are assigned to Residues")

    @property
    def n_segments(self):
        """Number of unique segments present in the :class:`AtomGroup`.

        Equivalent to ``len(self.segments)``.
        """
        return len(self.segments)

    @property
    @cached('unique_restore_mask')
    def _unique_restore_mask(self):
        # The _unique_restore_mask property's cache is populated whenever the
        # AtomGroup.unique property of a *non-unique* AtomGroup is accessed.
        # If _unique_restore_mask is not cached, it is *definitely* used in the
        # wrong place, so we raise an exception here. In principle, the
        # exception should be an AttributeError, but the error message would
        # then be replaced by the __getattr__() error message. To prevent the
        # message from being overridden, we raise a RuntimeError instead.
        if self.isunique:
            msg = ("{0}._unique_restore_mask is not available if the {0} is "
                   "unique. ".format(self.__class__.__name__))
        else:
            msg = ("{0}._unique_restore_mask is only available after "
                   "accessing {0}.unique. ".format(self.__class__.__name__))
        msg += ("If you see this error message in an unmodified release "
                "version of MDAnalysis, this is almost certainly a bug!")
        raise RuntimeError(msg)

    @_unique_restore_mask.setter
    def _unique_restore_mask(self, mask):
        self._cache['unique_restore_mask'] = mask

    @property
    @cached('unique')
    def unique(self):
        """An :class:`AtomGroup` containing sorted and unique
        :class:`Atoms<Atom>` only.

        If the :class:`AtomGroup` is unique, this is the group itself.

        Examples
        --------

           >>> ag = u.atoms[[2, 1, 2, 2, 1, 0]]
           >>> ag
           <AtomGroup with 6 atoms>
           >>> ag.ix
           array([2, 1, 2, 2, 1, 0])
           >>> ag2 = ag.unique
           >>> ag2
           <AtomGroup with 3 atoms>
           >>> ag2.ix
           array([0, 1, 2])
           >>> ag2.unique is ag2
           True


        .. versionadded:: 0.16.0
        .. versionchanged:: 0.19.0 If the :class:`AtomGroup` is already unique,
            :attr:`AtomGroup.unique` now returns the group itself instead of a
            copy.
        """
        if self.isunique:
            return self
        unique_ix, restore_mask = np.unique(self.ix, return_inverse=True)
        _unique = self.universe.atoms[unique_ix]
        self._unique_restore_mask = restore_mask
        # Since we know that _unique is a unique AtomGroup, we set its
        # uniqueness caches from here:
        _unique._cache['isunique'] = True
        _unique._cache['unique'] = _unique
        return _unique

    @property
    def positions(self):
        """Coordinates of the :class:`Atoms<Atom>` in the :class:`AtomGroup`.

        A :class:`numpy.ndarray` with
        :attr:`~numpy.ndarray.shape`\ ``=(``\ :attr:`~AtomGroup.n_atoms`\ ``, 3)``
        and :attr:`~numpy.ndarray.dtype`\ ``=numpy.float32``.

        The positions can be changed by assigning an array of the appropriate
        shape, i.e., either ``(``\ :attr:`~AtomGroup.n_atoms`\ ``, 3)`` to
        assign individual coordinates, or ``(3,)`` to assign the *same*
        coordinate to all :class:`Atoms<Atom>` (e.g.,
        ``ag.positions = array([0,0,0])`` will move all :class:`Atoms<Atom>`
        to the origin).

        .. note:: Changing positions is not reflected in any files; reading any
                  frame from the
                  :attr:`~MDAnalysis.core.universe.Universe.trajectory` will
                  replace the change with that from the file *except* if the
                  :attr:`~MDAnalysis.core.universe.Universe.trajectory` is held
                  in memory, e.g., when the
                  :meth:`~MDAnalysis.core.universe.Universe.transfer_to_memory`
                  method was used.

        Raises
        ------
        ~MDAnalysis.exceptions.NoDataError
            If the underlying :class:`~MDAnalysis.coordinates.base.Timestep`
            does not contain
            :attr:`~MDAnalysis.coordinates.base.Timestep.positions`.
        """
        return self.universe.trajectory.ts.positions[self.ix]

    @positions.setter
    def positions(self, values):
        ts = self.universe.trajectory.ts
        ts.positions[self.ix, :] = values

    @property
    def velocities(self):
        """Velocities of the :class:`Atoms<Atom>` in the :class:`AtomGroup`.

        A :class:`numpy.ndarray` with
        :attr:`~numpy.ndarray.shape`\ ``=(``\ :attr:`~AtomGroup.n_atoms`\ ``, 3)``
        and :attr:`~numpy.ndarray.dtype`\ ``=numpy.float32``.

        The velocities can be changed by assigning an array of the appropriate
        shape, i.e. either ``(``\ :attr:`~AtomGroup.n_atoms`\ ``, 3)`` to assign
        individual velocities or ``(3,)`` to assign the *same* velocity to all
        :class:`Atoms<Atom>` (e.g. ``ag.velocities = array([0,0,0])`` will give
        all :class:`Atoms<Atom>` zero :attr:`~Atom.velocity`).

        Raises
        ------
        ~MDAnalysis.exceptions.NoDataError
            If the underlying :class:`~MDAnalysis.coordinates.base.Timestep`
            does not contain
            :attr:`~MDAnalysis.coordinates.base.Timestep.velocities`.
        """
        ts = self.universe.trajectory.ts
        try:
            return np.array(ts.velocities[self.ix])
        except (AttributeError, NoDataError):
            raise NoDataError("Timestep does not contain velocities")

    @velocities.setter
    def velocities(self, values):
        ts = self.universe.trajectory.ts
        try:
            ts.velocities[self.ix, :] = values
        except (AttributeError, NoDataError):
            raise NoDataError("Timestep does not contain velocities")

    @property
    def forces(self):
        """Forces on each :class:`Atom` in the :class:`AtomGroup`.

        A :class:`numpy.ndarray` with
        :attr:`~numpy.ndarray.shape`\ ``=(``\ :attr:`~AtomGroup.n_atoms`\ ``, 3)``
        and :attr:`~numpy.ndarray.dtype`\ ``=numpy.float32``.

        The forces can be changed by assigning an array of the appropriate
        shape, i.e. either ``(``\ :attr:`~AtomGroup.n_atoms`\ ``, 3)`` to assign
        individual forces or ``(3,)`` to assign the *same* force to all
        :class:`Atoms<Atom>` (e.g. ``ag.forces = array([0,0,0])`` will give all
        :class:`Atoms<Atom>` a zero :attr:`~Atom.force`).

        Raises
        ------
        ~MDAnalysis.exceptions.NoDataError
            If the :class:`~MDAnalysis.coordinates.base.Timestep` does not
            contain :attr:`~MDAnalysis.coordinates.base.Timestep.forces`.
        """
        ts = self.universe.trajectory.ts
        try:
            return ts.forces[self.ix]
        except (AttributeError, NoDataError):
            raise NoDataError("Timestep does not contain forces")

    @forces.setter
    def forces(self, values):
        ts = self.universe.trajectory.ts
        try:
            ts.forces[self.ix, :] = values
        except (AttributeError, NoDataError):
            raise NoDataError("Timestep does not contain forces")

    @property
    def ts(self):
        """Temporary Timestep that contains the selection coordinates.

        A :class:`~MDAnalysis.coordinates.base.Timestep` instance,
        which can be passed to a trajectory writer.

        If :attr:`~AtomGroup.ts` is modified then these modifications
        will be present until the frame number changes (which
        typically happens when the underlying
        :attr:`~MDAnalysis.core.universe.Universe.trajectory` frame changes).

        It is not possible to assign a new
        :class:`~MDAnalysis.coordinates.base.Timestep` to the
        :attr:`AtomGroup.ts` attribute; change attributes of the object.
        """
        trj_ts = self.universe.trajectory.ts  # original time step

        return trj_ts.copy_slice(self.indices)

    # As with universe.select_atoms, needing to fish out specific kwargs
    # (namely, 'updating') doesn't allow a very clean signature.
    def select_atoms(self, sel, *othersel, **selgroups):
        """Select :class:`Atoms<Atom>` using a selection string.

        Returns an :class:`AtomGroup` with :class:`Atoms<Atom>` sorted according
        to their index in the topology (this is to ensure that there are no
        duplicates, which can happen with complicated selections).

        Raises
        ------
        TypeError
            If the arbitrary groups passed are not of type
            :class:`MDAnalysis.core.groups.AtomGroup`

        Examples
        --------
        All simple selection listed below support multiple arguments which are
        implicitly combined with an or operator. For example

           >>> sel = universe.select_atoms('resname MET GLY')

        is equivalent to

           >>> sel = universe.select_atoms('resname MET or resname GLY')

        Will select all atoms with a residue name of either MET or GLY.

        Subselections can be grouped with parentheses.

           >>> sel = universe.select_atoms("segid DMPC and not ( name H* O* )")
           >>> sel
           <AtomGroup with 3420 atoms>


        Existing :class:`AtomGroup` objects can be passed as named arguments,
        which will then be available to the selection parser.

           >>> universe.select_atoms("around 10 group notHO", notHO=sel)
           <AtomGroup with 1250 atoms>

        Selections can be set to update automatically on frame change, by
        setting the `updating` keyword argument to `True`.  This will return
        a :class:`UpdatingAtomGroup` which can represent the solvation shell
        around another object.

           >>> universe.select_atoms("resname SOL and around 2.0 protein", updating=True)
           <Updating AtomGroup with 100 atoms>

        Notes
        -----

        If exact ordering of atoms is required (for instance, for
        :meth:`~AtomGroup.angle` or :meth:`~AtomGroup.dihedral` calculations)
        then one supplies selections *separately* in the required order. Also,
        when multiple :class:`AtomGroup` instances are concatenated with the
        ``+`` operator, then the order of :class:`Atom` instances is preserved
        and duplicates are *not* removed.


        See Also
        --------
        :ref:`selection-commands-label` for further details and examples.


        .. rubric:: Selection syntax


        The selection parser understands the following CASE SENSITIVE
        *keywords*:

        **Simple selections**

            protein, backbone, nucleic, nucleicbackbone
                selects all atoms that belong to a standard set of residues;
                a protein is identfied by a hard-coded set of residue names so
                it  may not work for esoteric residues.
            segid *seg-name*
                select by segid (as given in the topology), e.g. ``segid 4AKE``
                or ``segid DMPC``
            resid *residue-number-range*
                resid can take a single residue number or a range of numbers. A
                range consists of two numbers separated by a colon (inclusive)
                such as ``resid 1:5``. A residue number ("resid") is taken
                directly from the topology.
                If icodes are present in the topology, then these will be
                taken into account.  Ie 'resid 163B' will only select resid
                163 with icode B while 'resid 163' will select only residue 163.
                Range selections will also respect icodes, so 'resid 162-163B'
                will select all residues in 162 and those in 163 up to icode B.
            resnum *resnum-number-range*
                resnum is the canonical residue number; typically it is set to
                the residue id in the original PDB structure.
            resname *residue-name*
                select by residue name, e.g. ``resname LYS``
            name *atom-name*
                select by atom name (as given in the topology). Often, this is
                force field dependent. Example: ``name CA`` (for C&alpha; atoms)
                or ``name OW`` (for SPC water oxygen)
            type *atom-type*
                select by atom type; this is either a string or a number and
                depends on the force field; it is read from the topology file
                (e.g. the CHARMM PSF file contains numeric atom types). It has
                non-sensical values when a PDB or GRO file is used as a topology
            atom *seg-name*  *residue-number*  *atom-name*
                a selector for a single atom consisting of segid resid atomname,
                e.g. ``DMPC 1 C2`` selects the C2 carbon of the first residue of
                the DMPC segment
            altloc *alternative-location*
                a selection for atoms where alternative locations are available,
                which is often the case with high-resolution crystal structures
                e.g. `resid 4 and resname ALA and altloc B` selects only the
                atoms of ALA-4 that have an altloc B record.
            moltype *molecule-type*
                select by molecule type, e.g. ``moltype Protein_A``. At the
                moment, only the TPR format defines the molecule type.

        **Boolean**

            not
                all atoms not in the selection, e.g. ``not protein`` selects
                all atoms that aren't part of a protein

            and, or
                combine two selections according to the rules of boolean
                algebra, e.g. ``protein and not resname ALA LYS``
                selects all atoms that belong to a protein, but are not in a
                lysine or alanine residue

        **Geometric**

            around *distance*  *selection*
                selects all atoms a certain cutoff away from another selection,
                e.g. ``around 3.5 protein`` selects all atoms not belonging to
                protein that are within 3.5 Angstroms from the protein
            point *x* *y* *z*  *distance*
                selects all atoms within a cutoff of a point in space, make sure
                coordinate is separated by spaces,
                e.g. ``point 5.0 5.0 5.0  3.5`` selects all atoms within 3.5
                Angstroms of the coordinate (5.0, 5.0, 5.0)
            prop [abs] *property*  *operator*  *value*
                selects atoms based on position, using *property*  **x**, **y**,
                or **z** coordinate. Supports the **abs** keyword (for absolute
                value) and the following *operators*: **<, >, <=, >=, ==, !=**.
                For example, ``prop z >= 5.0`` selects all atoms with z
                coordinate greater than 5.0; ``prop abs z <= 5.0`` selects all
                atoms within -5.0 <= z <= 5.0.
            sphzone *radius* *selection*
                Selects all atoms that are within *radius* of the center of
                geometry of *selection*
            sphlayer *inner radius* *outer radius* *selection*
                Similar to sphzone, but also excludes atoms that are within
                *inner radius* of the selection COG
            cyzone *externalRadius* *zMax* *zMin* *selection*
                selects all atoms within a cylindric zone centered in the
                center of geometry (COG) of a given selection,
                e.g. ``cyzone 15 4 -8 protein and resid 42`` selects the
                center of geometry of protein and resid 42, and creates a
                cylinder of external radius 15 centered on the COG. In z, the
                cylinder extends from 4 above the COG to 8 below. Positive
                values for *zMin*, or negative ones for *zMax*, are allowed.
            cylayer *innerRadius* *externalRadius* *zMax* *zMin* *selection*
                selects all atoms within a cylindric layer centered in the
                center of geometry (COG) of a given selection,
                e.g. ``cylayer 5 10 10 -8 protein`` selects the center of
                geometry of protein, and creates a cylindrical layer of inner
                radius 5, external radius 10 centered on the COG. In z, the
                cylinder extends from 10 above the COG to 8 below. Positive
                values for *zMin*, or negative ones for *zMax*, are allowed.

        **Connectivity**

            byres *selection*
                selects all atoms that are in the same segment and residue as
                selection, e.g. specify the subselection after the byres keyword
            bonded *selection*
                selects all atoms that are bonded to selection
                eg: ``select name H and bonded name O`` selects only hydrogens
                bonded to oxygens

        **Index**

            bynum *index-range*
                selects all atoms within a range of (1-based) inclusive indices,
                e.g. ``bynum 1`` selects the first atom in the universe;
                ``bynum 5:10`` selects atoms 5 through 10 inclusive. All atoms
                in the :class:`~MDAnalysis.core.universe.Universe` are
                consecutively numbered, and the index runs from 1 up to the
                total number of atoms.

        **Preexisting selections**

            group `group-name`
                selects the atoms in the :class:`AtomGroup` passed to the
                function as an argument named `group-name`. Only the atoms
                common to `group-name` and the instance
                :meth:`~MDAnalysis.core.groups.AtomGroup.select_atoms`
                was called from will be considered, unless ``group`` is
                preceded by the ``global`` keyword. `group-name` will be
                included in the parsing just by comparison of atom indices.
                This means that it is up to the user to make sure the
                `group-name` group was defined in an appropriate
                :class:`~MDAnalysis.core.universe.Universe`.

            global *selection*
                by default, when issuing
                :meth:`~MDAnalysis.core.groups.AtomGroup.select_atoms` from an
                :class:`~MDAnalysis.core.groups.AtomGroup`, selections and
                subselections are returned intersected with the atoms of that
                instance. Prefixing a selection term with ``global`` causes its
                selection to be returned in its entirety.  As an example, the
                ``global`` keyword allows for
                ``lipids.select_atoms("around 10 global protein")`` --- where
                ``lipids`` is a group that does not contain any proteins. Were
                ``global`` absent, the result would be an empty selection since
                the ``protein`` subselection would itself be empty. When issuing
                :meth:`~MDAnalysis.core.groups.AtomGroup.select_atoms` from a
                :class:`~MDAnalysis.core.universe.Universe`, ``global`` is
                ignored.

        **Dynamic selections**
            If :meth:`~MDAnalysis.core.groups.AtomGroup.select_atoms` is
            invoked with named argument `updating` set to `True`, an
            :class:`~MDAnalysis.core.groups.UpdatingAtomGroup` instance will be
            returned, instead of a regular
            :class:`~MDAnalysis.core.groups.AtomGroup`. It behaves just like
            the latter, with the difference that the selection expressions are
            re-evaluated every time the trajectory frame changes (this happens
            lazily, only when the
            :class:`~MDAnalysis.core.groups.UpdatingAtomGroup` is accessed so
            that there is no redundant updating going on).
            Issuing an updating selection from an already updating group will
            cause later updates to also reflect the updating of the base group.
            A non-updating selection or a slicing operation made on an
            :class:`~MDAnalysis.core.groups.UpdatingAtomGroup` will return a
            static :class:`~MDAnalysis.core.groups.AtomGroup`, which will no
            longer update across frames.


        .. versionchanged:: 0.7.4 Added *resnum* selection.
        .. versionchanged:: 0.8.1 Added *group* and *fullgroup* selections.
        .. deprecated:: 0.11 The use of *fullgroup* has been deprecated in favor
            of the equivalent *global group* selections.
        .. versionchanged:: 0.13.0 Added *bonded* selection.
        .. versionchanged:: 0.16.0 Resid selection now takes icodes into account
            where present.
        .. versionchanged:: 0.16.0 Updating selections now possible by setting
            the `updating` argument.
        .. versionchanged:: 0.17.0 Added *moltype* and *molnum* selections.
        .. versionchanged:: 0.18.1 Added strict type checking for passed groups.
        """
        updating = selgroups.pop('updating', False)
        sel_strs = (sel,) + othersel

        for group, thing in selgroups.items():
            if not isinstance(thing, AtomGroup):
                raise TypeError("Passed groups must be AtomGroups. "
                                "You provided {} for group '{}'".format(
                                    thing.__class__.__name__, group))

        selections = tuple((selection.Parser.parse(s, selgroups)
                            for s in sel_strs))
        if updating:
            atomgrp = UpdatingAtomGroup(self, selections, sel_strs)
        else:
            # Apply the first selection and sum to it
            atomgrp = sum([sel.apply(self) for sel in selections[1:]],
                          selections[0].apply(self))
        return atomgrp

    def split(self, level):
        """Split :class:`AtomGroup` into a :class:`list` of
        :class:`AtomGroups<AtomGroup>` by `level`.

        Parameters
        ----------
        level : {'atom', 'residue', 'molecule', 'segment'}


        .. versionadded:: 0.9.0
        .. versionchanged:: 0.17.0 Added the 'molecule' level.
        """
        accessors = {'segment': 'segindices',
                     'residue': 'resindices',
                     'molecule': 'molnums'}

        if level == "atom":
            return [self.universe.atoms[[a.ix]] for a in self]

        # higher level groupings
        try:
            levelindices = getattr(self, accessors[level])
        except AttributeError:
            raise AttributeError('This universe does not have {} '
                             'information. Maybe it is not provided in the '
                             'topology format in use.'.format(level))
        except KeyError:
            raise ValueError("level = '{0}' not supported, "
                             "must be one of {1}".format(level,
                                                         accessors.keys()))

        return [self[levelindices == index] for index in
                unique_int_1d(levelindices)]

    def guess_bonds(self, vdwradii=None):
        """Guess bonds that exist within this :class:`AtomGroup` and add them to
        the underlying :attr:`~AtomGroup.universe`.

        Parameters
        ----------
        vdwradii : dict, optional
            Dict relating atom types: vdw radii


        See Also
        --------
        :func:`MDAnalysis.topology.guessers.guess_bonds`


        .. versionadded:: 0.10.0
        """
        from ..topology.core import guess_bonds, guess_angles, guess_dihedrals
        from .topologyattrs import Bonds, Angles, Dihedrals

        def get_TopAttr(u, name, cls):
            """either get *name* or create one from *cls*"""
            try:
                return getattr(u._topology, name)
            except AttributeError:
                attr = cls([])
                u.add_TopologyAttr(attr)
                return attr

        # indices of bonds
        b = guess_bonds(self.atoms, self.atoms.positions, vdwradii=vdwradii)
        bondattr = get_TopAttr(self.universe, 'bonds', Bonds)
        bondattr.add_bonds(b, guessed=True)

        a = guess_angles(self.bonds)
        angleattr = get_TopAttr(self.universe, 'angles', Angles)
        angleattr.add_bonds(a, guessed=True)

        d = guess_dihedrals(self.angles)
        diheattr = get_TopAttr(self.universe, 'dihedrals', Dihedrals)
        diheattr.add_bonds(d)

    @property
    def bond(self):
        """This :class:`AtomGroup` represented as a
        :class:`MDAnalysis.core.topologyobjects.Bond` object

        Raises
        ------
        ValueError
            If the :class:`AtomGroup` is not length 2


        .. versionadded:: 0.11.0
        """
        if len(self) != 2:
            raise ValueError(
                "bond only makes sense for a group with exactly 2 atoms")
        return topologyobjects.Bond(self.ix, self.universe)

    @property
    def angle(self):
        """This :class:`AtomGroup` represented as an
        :class:`MDAnalysis.core.topologyobjects.Angle` object

        Raises
        ------
        ValueError
            If the :class:`AtomGroup` is not length 3


        .. versionadded:: 0.11.0
        """
        if len(self) != 3:
            raise ValueError(
                "angle only makes sense for a group with exactly 3 atoms")
        return topologyobjects.Angle(self.ix, self.universe)

    @property
    def dihedral(self):
        """This :class:`AtomGroup` represented as a
        :class:`~MDAnalysis.core.topologyobjects.Dihedral` object

        Raises
        ------
        ValueError
            If the :class:`AtomGroup` is not length 4


        .. versionadded:: 0.11.0
        """
        if len(self) != 4:
            raise ValueError(
                "dihedral only makes sense for a group with exactly 4 atoms")
        return topologyobjects.Dihedral(self.ix, self.universe)

    @property
    def improper(self):
        """This :class:`AtomGroup` represented as an
        :class:`MDAnalysis.core.topologyobjects.ImproperDihedral` object

        Raises
        ------
        ValueError
            If the :class:`AtomGroup` is not length 4


        .. versionadded:: 0.11.0
        """
        if len(self) != 4:
            raise ValueError(
                "improper only makes sense for a group with exactly 4 atoms")
        return topologyobjects.ImproperDihedral(self.ix, self.universe)

    def write(self, filename=None, file_format="PDB",
              filenamefmt="{trjname}_{frame}", frames=None, **kwargs):
        """Write `AtomGroup` to a file.

        The output can either be a coordinate file or a selection, depending on
        the format. 

        Examples
        --------

        >>> ag = u.atoms
        >>> ag.write('selection.ndx')  # Write a gromacs index file
        >>> ag.write('coordinates.pdb')  # Write the current frame as PDB
        >>> # Write the trajectory in XTC format
        >>> ag.write('trajectory.xtc', frames='all')
        >>> # Write every other frame of the trajectory in PBD format
        >>> ag.write('trajectory.pdb', frames=u.trajectory[::2])

        Parameters
        ----------
        filename : str, optional
            ``None``: create TRJNAME_FRAME.FORMAT from filenamefmt [``None``]
        file_format : str, optional
            The name or extension of a coordinate, trajectory, or selection
            file format such as PDB, CRD, GRO, VMD (tcl), PyMol (pml), Gromacs
            (ndx) CHARMM (str) or Jmol (spt); case-insensitive [PDB]
        filenamefmt : str, optional
            format string for default filename; use substitution tokens
            'trjname' and 'frame' ["%(trjname)s_%(frame)d"]
        bonds : str, optional
            how to handle bond information, especially relevant for PDBs.
            ``"conect"``: write only the CONECT records defined in the original
            file. ``"all"``: write out all bonds, both the original defined and
            those guessed by MDAnalysis. ``None``: do not write out bonds.
            Default is ``"conect"``.
        frames: array-like or slice or FrameIteratorBase or str, optional
            An ensemble of frames to write. The ensemble can be an list or
            array of frame indices, a mask of booleans, an instance of
            :class:`slice`, or the value returned when a trajectory is indexed.
            By default, `frames` is set to ``None`` and only the current frame
            is written. If `frames` is set to "all", then all the frame from
            trajectory are written.


        .. versionchanged:: 0.9.0 Merged with write_selection. This method can
            now write both selections out.
        .. versionchanged:: 0.19.0
            Can write multiframe trajectories with the 'frames' argument.
        """
        # TODO: Add a 'verbose' option alongside 'frames'.

        # check that AtomGroup actually has any atoms (Issue #434)
        if len(self.atoms) == 0:
            raise IndexError("Cannot write an AtomGroup with 0 atoms")

        trj = self.universe.trajectory  # unified trajectory API
        if frames is None or frames == 'all':
            trj_frames = trj[::]
        elif isinstance(frames, numbers.Integral):
            # We accept everything that indexes a trajectory and returns a
            # subset of it. Though, numbers return a Timestep instead.
            raise TypeError('The "frames" argument cannot be a number.')
        else:
            try:
                test_trajectory = frames.trajectory
            except AttributeError:
                trj_frames = trj[frames]
            else:
                if test_trajectory is not trj:
                    raise ValueError(
                        'The trajectory of {} provided to the frames keyword '
                        'attribute is different from the trajectory of the '
                        'AtomGroup.'.format(frames)
                    )
                trj_frames = frames

        if filename is None:
            trjname, ext = os.path.splitext(os.path.basename(trj.filename))
            filename = filenamefmt.format(trjname=trjname, frame=trj.frame)
        filename = util.filename(filename, ext=file_format.lower(), keep=True)

        # Some writer behave differently when they are given a "multiframe"
        # argument. It is the case of the PDB writer tht writes models when
        # "multiframe" is True.
        # We want to honor what the user provided with the argument if
        # provided explicitly. If not, then we need to figure out if we write
        # multiple frames or not.
        multiframe = kwargs.pop('multiframe', None)
        if len(trj_frames) > 1 and multiframe == False:
            raise ValueError(
                'Cannot explicitely set "multiframe" to False and request '
                'more than 1 frame with the "frames" keyword argument.'
            ) 
        elif multiframe is None:
            if frames is None:
                # By default we only write the current frame.
                multiframe = False
            else:
                multiframe = len(trj_frames) > 1

        # From the following blocks, one must pass.
        # Both can't pass as the extensions don't overlap.
        # Try and select a Class using get_ methods (becomes `writer`)
        # Once (and if!) class is selected, use it in with block
        try:
            # format keyword works differently in get_writer and get_selection_writer
            # here it overrides everything, in get_sel it is just a default
            # apply sparingly here!
            format = os.path.splitext(filename)[1][1:]  # strip initial dot!
            format = format or file_format
            format = format.strip().upper()

            writer = get_writer_for(filename, format=format, multiframe=multiframe)
        except (ValueError, TypeError):
            pass
        else:
            with writer(filename, n_atoms=self.n_atoms, **kwargs) as w:
                if frames is None:
                    w.write(self.atoms)
                else:
                    current_frame = trj.ts.frame
                    try:
                        for _ in trj_frames:
                            w.write(self.atoms)
                    finally:
                        trj[current_frame]
            return

        try:
            # here `file_format` is only used as default,
            # anything pulled off `filename` will be used preferentially
            writer = get_selection_writer_for(filename, file_format)
        except (TypeError, NotImplementedError):
            pass
        else:
            with writer(filename, n_atoms=self.n_atoms, **kwargs) as w:
                w.write(self.atoms)
            return

        raise ValueError("No writer found for format: {}".format(filename))


class ResidueGroup(GroupBase):
    """ResidueGroup base class.

    This class is used by a :class:`~MDAnalysis.core.universe.Universe` for
    generating its Topology-specific :class:`ResidueGroup` class. All the
    :class:`~MDAnalysis.core.topologyattrs.TopologyAttr` components are obtained
    from :class:`GroupBase`, so this class only includes ad-hoc methods
    specific to :class:`ResidueGroups<ResidueGroup>`.

    ResidueGroups can be compared and combined using group operators. See the
    list of these operators on :class:`GroupBase`.

    .. deprecated:: 0.16.2
       *Instant selectors* of Segments will be removed in the 1.0 release.
       See :ref:`Instant selectors <instance-selectors>` for details and
       alternatives.
    """

    @property
    def atoms(self):
        """An :class:`AtomGroup` of :class:`Atoms<Atom>` present in this
        :class:`ResidueGroup`.

        The :class:`Atoms<Atom>` are ordered locally by :class:`Residue` in the
        :class:`ResidueGroup`.  Duplicates are *not* removed.
        """
        ag = self.universe.atoms[np.concatenate(self.indices)]
        # If the ResidueGroup is known to be unique, this also holds for the
        # atoms therein, since atoms can only belong to one residue at a time.
        # On the contrary, if the ResidueGroup is not unique, this does not
        # imply non-unique atoms, since residues might be empty.
        try:
            if self._cache['isunique']:
                ag._cache['isunique'] = True
                ag._cache['unique'] = ag
        except KeyError:
            pass
        return ag

    @property
    def n_atoms(self):
        """Number of :class:`Atoms<Atom>` present in this :class:`ResidueGroup`,
        including duplicate residues (and thus, duplicate atoms).

        Equivalent to ``len(self.atoms)``.
        """
        return len(self.atoms)

    @property
    def residues(self):
        """The :class:`ResidueGroup` itself.

        See Also
        --------
        copy : return a true copy of the :class:`ResidueGroup`


        .. versionchanged:: 0.19.0
           In previous versions, this returned a copy, but now
           the :class:`ResidueGroup` itself is returned. This should
           not affect any code but only speed up calculations.

        """
        return self

    @property
    def n_residues(self):
        """Number of residues in the :class:`ResidueGroup`.

        Equivalent to ``len(self)``.
        """
        return len(self)

    @property
    def segments(self):
        """Get sorted :class:`SegmentGroup` of the unique segments present in
        the :class:`ResidueGroup`.
        """
        sg = self.universe.segments[unique_int_1d(self.segindices)]
        sg._cache['isunique'] = True
        sg._cache['unique'] = sg
        return sg

    @segments.setter
    def segments(self, new):
        # Can set with Seg, SegGroup or list/tuple of Seg
        if isinstance(new, Segment):
            s_ix = itertools.cycle((new.segindex,))
        elif isinstance(new, SegmentGroup):
            s_ix = new.segindices
        else:
            try:
                s_ix = [s.segindex for s in new]
            except AttributeError:
                raise TypeError("Can only set ResidueGroup segments to Segment "
                                "or SegmentGroup, not {}".format(
                                    ', '.join(type(r) for r in new
                                              if not isinstance(r, Segment))
                                ))
        if not isinstance(s_ix, itertools.cycle) and len(s_ix) != len(self):
            raise ValueError("Incorrect size: {} for ResidueGroup of size: {}"
                             "".format(len(new), len(self)))
        # Optimisation TODO:
        # This currently rebuilds the tt len(self) times
        # Ideally all changes would happen and *afterwards* tables are built
        # Alternatively, if the changes didn't rebuild table, this list
        # comprehension isn't terrible.
        for r, s in zip(self, s_ix):
            self.universe._topology.tt.move_residue(r.ix, s)

    @property
    def n_segments(self):
        """Number of unique segments present in the ResidueGroup.

        Equivalent to ``len(self.segments)``.
        """
        return len(self.segments)

    @property
    @cached('unique')
    def unique(self):
        """Return a :class:`ResidueGroup` containing sorted and unique
        :class:`Residues<Residue>` only.

        If the :class:`ResidueGroup` is unique, this is the group itself.

        Examples
        --------

           >>> rg = u.residues[[2, 1, 2, 2, 1, 0]]
           >>> rg
           <ResidueGroup with 6 residues>
           >>> rg.ix
           array([2, 1, 2, 2, 1, 0])
           >>> rg2 = rg.unique
           >>> rg2
           <ResidueGroup with 3 residues>
           >>> rg2.ix
           array([0, 1, 2])
           >>> rg2.unique is rg2
           True


        .. versionadded:: 0.16.0
        .. versionchanged:: 0.19.0 If the :class:`ResidueGroup` is already
            unique, :attr:`ResidueGroup.unique` now returns the group itself
            instead of a copy.
        """
        if self.isunique:
            return self
        _unique = self.universe.residues[unique_int_1d(self.ix)]
        # Since we know that _unique is a unique ResidueGroup, we set its
        # uniqueness caches from here:
        _unique._cache['isunique'] = True
        _unique._cache['unique'] = _unique
        return _unique


class SegmentGroup(GroupBase):
    """:class:`SegmentGroup` base class.

    This class is used by a :class:`~MDAnalysis.core.universe.Universe` for
    generating its Topology-specific :class:`SegmentGroup` class. All the
    :class:`~MDAnalysis.core.topologyattrs.TopologyAttr` components are obtained
    from :class:`GroupBase`, so this class only includes ad-hoc methods specific
    to :class:`SegmentGroups<SegmentGroup>`.

    :class:`SegmentGroups<SegmentGroup>` can be compared and combined using
    group operators. See the list of these operators on :class:`GroupBase`.

    .. deprecated:: 0.16.2
       *Instant selectors* of Segments will be removed in the 1.0 release.
       See :ref:`Instant selectors <instance-selectors>` for details and
       alternatives.
    """

    @property
    def atoms(self):
        """An :class:`AtomGroup` of :class:`Atoms<Atom>` present in this
        :class:`SegmentGroup`.

        The :class:`Atoms<Atom>` are ordered locally by :class:`Residue`, which
        are further ordered by :class:`Segment` in the :class:`SegmentGroup`.
        Duplicates are *not* removed.
        """
        ag = self.universe.atoms[np.concatenate(self.indices)]
        # If the SegmentGroup is known to be unique, this also holds for the
        # residues therein, and thus, also for the atoms in those residues.
        # On the contrary, if the SegmentGroup is not unique, this does not
        # imply non-unique atoms, since segments or residues might be empty.
        try:
            if self._cache['isunique']:
                ag._cache['isunique'] = True
                ag._cache['unique'] = ag
        except KeyError:
            pass
        return ag

    @property
    def n_atoms(self):
        """Number of atoms present in the :class:`SegmentGroup`, including
        duplicate segments (and thus, duplicate atoms).

        Equivalent to ``len(self.atoms)``.
        """
        return len(self.atoms)

    @property
    def residues(self):
        """A :class:`ResidueGroup` of :class:`Residues<Residue>` present in this
        :class:`SegmentGroup`.

        The :class:`Residues<Residue>` are ordered locally by
        :class:`Segment` in the :class:`SegmentGroup`. Duplicates are *not*
        removed.
        """
        rg = self.universe.residues[np.concatenate(self.resindices)]
        # If the SegmentGroup is known to be unique, this also holds for the
        # residues therein. On the contrary, if the SegmentGroup is not unique,
        # this does not imply non-unique residues, since segments might be
        # empty.
        try:
            if self._cache['isunique']:
                rg._cache['isunique'] = True
                rg._cache['unique'] = rg
        except KeyError:
            pass
        return rg

    @property
    def n_residues(self):
        """Number of residues present in this :class:`SegmentGroup`, including
        duplicate segments (and thus, residues).

        Equivalent to ``len(self.residues)``.
        """
        return len(self.residues)

    @property
    def segments(self):
        """The :class:`SegmentGroup` itself.

        See Also
        --------
        copy : return a true copy of the :class:`SegmentGroup`


        .. versionchanged:: 0.19.0
           In previous versions, this returned a copy, but now
           the :class:`SegmentGroup` itself is returned. This should
           not affect any code but only speed up calculations.

        """
        return self

    @property
    def n_segments(self):
        """Number of segments in the :class:`SegmentGroup`.

        Equivalent to ``len(self)``.
        """
        return len(self)

    @property
    @cached('unique')
    def unique(self):
        """Return a :class:`SegmentGroup` containing sorted and unique
        :class:`Segments<Segment>` only.

        If the :class:`SegmentGroup` is unique, this is the group itself.

        Examples
        --------

           >>> sg = u.segments[[2, 1, 2, 2, 1, 0]]
           >>> sg
           <SegmentGroup with 6 segments>
           >>> sg.ix
           array([2, 1, 2, 2, 1, 0])
           >>> sg2 = sg.unique
           >>> sg2
           <SegmentGroup with 3 segments>
           >>> sg2.ix
           array([0, 1, 2])
           >>> sg2.unique is sg2
           True


        .. versionadded:: 0.16.0
        .. versionchanged:: 0.19.0 If the :class:`SegmentGroup` is already
            unique, :attr:`SegmentGroup.unique` now returns the group itself
            instead of a copy.
        """
        if self.isunique:
            return self
        _unique = self.universe.segments[unique_int_1d(self.ix)]
        # Since we know that _unique is a unique SegmentGroup, we set its
        # uniqueness caches from here:
        _unique._cache['isunique'] = True
        _unique._cache['unique'] = _unique
        return _unique


@functools.total_ordering
class ComponentBase(_MutableBase):
    """Base class from which a :class:`~MDAnalysis.core.universe.Universe`\ 's
    Component class is built.

    Components are the individual objects that are found in Groups.
    """
    def __init__(self, ix, u):
        # index of component
        self._ix = ix
        self._u = u

    def __lt__(self, other):
        if self.level != other.level:
            raise TypeError("Can't compare different level objects")
        return self.ix < other.ix

    def __eq__(self, other):
        if self.level != other.level:
            raise TypeError("Can't compare different level objects")
        return self.ix == other.ix

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        return hash(self.ix)

    @_only_same_level
    def __add__(self, other):
        """Concatenate the Component with another Component or Group of the
        same level.

        Parameters
        ----------
        other : Component or Group
            Component or Group with `other.level` same as `self.level`

        Returns
        -------
        Group
            Group with elements of `self` and `other` concatenated
        """
        o_ix = other.ix_array

        return self.level.plural(
                np.concatenate([self.ix_array, o_ix]), self.universe)

    def __radd__(self, other):
        """Using built-in sum requires supporting 0 + self. If other is
        anything other 0, an exception will be raised.

        Parameters
        ----------
        other : int
            Other should be 0, or else an exception will be raised.

        Returns
        -------
        self
            Group with elements of `self` reproduced
        """
        if other == 0:
            return self.level.plural(self.ix_array, self.universe)
        else:
            raise TypeError("unsupported operand type(s) for +:"
                            " '{}' and '{}'".format(type(self).__name__,
                                                    type(other).__name__))

    @property
    def universe(self):
        return self._u

    @property
    def ix(self):
        """Unique index of this component.

        If this component is an :class:`Atom`, this is the index of the
        :class:`Atom`.
        If it is a :class:`Residue`, this is the index of the :class:`Residue`.
        If it is a :class:`Segment`, this is the index of the :class:`Segment`.
        """
        return self._ix

    @property
    def ix_array(self):
        """Unique index of this component as an array.

        This method gives a consistent API between components and groups.

        See Also
        --------
        ix
        """
        return np.array([self.ix], dtype=np.intp)


class Atom(ComponentBase):
    """:class:`Atom` base class.

    This class is used by a :class:`~MDAnalysis.core.universe.Universe` for
    generating its Topology-specific :class:`Atom` class. All the
    :class:`~MDAnalysis.core.topologyattrs.TopologyAttr` components are obtained
    from :class:`ComponentBase`, so this class only includes ad-hoc methods
    specific to :class:`Atoms<Atom>`.
    """
    def __getattr__(self, attr):
        """Try and catch known attributes and give better error message"""
        if attr in ('fragment',):
            raise NoDataError("Atom has no fragment data, this requires Bonds")
        else:
            raise AttributeError("{cls} has no attribute {attr}".format(
                cls=self.__class__.__name__, attr=attr))

    def __repr__(self):
        me = '<Atom {}:'.format(self.ix + 1)
        if hasattr(self, 'name'):
            me += ' {}'.format(self.name)
        if hasattr(self, 'type'):
            me += ' of type {}'.format(self.type)
        if hasattr(self, 'resname'):
            me += ' of resname {},'.format(self.resname)
        if hasattr(self, 'resid'):
            me += ' resid {}'.format(self.resid)
        if hasattr(self, 'segid'):
            me += ' and segid {}'.format(self.segid)
        if hasattr(self, 'altLoc'):
            me += ' and altLoc {}'.format(self.altLoc)
        return me + '>'

    @property
    def residue(self):
        return self.universe.residues[self.universe._topology.resindices[self]]

    @residue.setter
    def residue(self, new):
        if not isinstance(new, Residue):
            raise TypeError("Can only set Atom residue to Residue, not {}"
                            "".format(type(new)))
        self.universe._topology.tt.move_atom(self.ix, new.resindex)

    @property
    def segment(self):
        return self.universe.segments[self.universe._topology.segindices[self]]

    @segment.setter
    def segment(self, new):
        raise NotImplementedError("Cannot set atom segment.  "
                                  "Segments are assigned to Residues")

    @property
    def position(self):
        """Coordinates of the atom.

        The position can be changed by assigning an array of length (3,).

        .. note:: changing the position is not reflected in any files; reading
                  any frame from the trajectory will replace the change with
                  that from the file

        Raises
        ------
        ~MDAnalysis.exceptions.NoDataError
            If the underlying :class:`~MDAnalysis.coordinates.base.Timestep`
            does not contain
            :attr:`~MDAnalysis.coordinates.base.Timestep.positions`.
        """
        return self.universe.trajectory.ts.positions[self.ix].copy()

    @position.setter
    def position(self, values):
        self.universe.trajectory.ts.positions[self.ix, :] = values

    @property
    def velocity(self):
        """Velocity of the atom.

        The velocity can be changed by assigning an array of shape ``(3,)``.

        .. note:: changing the velocity is not reflected in any files; reading
                  any frame from the trajectory will replace the change with
                  that from the file

        Raises
        ------
        ~MDAnalysis.exceptions.NoDataError
            If the underlying :class:`~MDAnalysis.coordinates.base.Timestep`
            does not contain
            :attr:`~MDAnalysis.coordinates.base.Timestep.velocities`.
        """
        ts = self.universe.trajectory.ts
        try:
            return ts.velocities[self.ix].copy()
        except (AttributeError, NoDataError):
            raise NoDataError("Timestep does not contain velocities")

    @velocity.setter
    def velocity(self, values):
        ts = self.universe.trajectory.ts
        try:
            ts.velocities[self.ix, :] = values
        except (AttributeError, NoDataError):
            raise NoDataError("Timestep does not contain velocities")

    @property
    def force(self):
        """Force on the atom.

        The force can be changed by assigning an array of shape ``(3,)``.

        .. note:: changing the force is not reflected in any files; reading any
                  frame from the trajectory will replace the change with that
                  from the file

        Raises
        ------
        ~MDAnalysis.exceptions.NoDataError
            If the underlying :class:`~MDAnalysis.coordinates.base.Timestep`
            does not contain
            :attr:`~MDAnalysis.coordinates.base.Timestep.forces`.
        """
        ts = self.universe.trajectory.ts
        try:
            return ts.forces[self.ix].copy()
        except (AttributeError, NoDataError):
            raise NoDataError("Timestep does not contain forces")

    @force.setter
    def force(self, values):
        ts = self.universe.trajectory.ts
        try:
            ts.forces[self.ix, :] = values
        except (AttributeError, NoDataError):
            raise NoDataError("Timestep does not contain forces")


class Residue(ComponentBase):
    """:class:`Residue` base class.

    This class is used by a :class:`~MDAnalysis.core.universe.Universe` for
    generating its Topology-specific :class:`Residue` class. All the
    :class:`~MDAnalysis.core.topologyattrs.TopologyAttr` components are obtained
    from :class:`ComponentBase`, so this class only includes ad-hoc methods
    specific to :class:`Residues<Residue>`.
    """
    def __repr__(self):
        me = '<Residue'
        if hasattr(self, 'resname'):
            me += ' {},'.format(self.resname)
        if hasattr(self, 'resid'):
            me += ' {}'.format(self.resid)

        return me + '>'

    @property
    def atoms(self):
        """An :class:`AtomGroup` of :class:`Atoms<Atom>` present in this
        :class:`Residue`.
        """
        ag = self.universe.atoms[self.universe._topology.indices[self][0]]
        ag._cache['isunique'] = True
        ag._cache['unique'] = ag
        return ag

    @property
    def segment(self):
        """The :class:`Segment` this :class:`Residue` belongs to.
        """
        return self.universe.segments[self.universe._topology.segindices[self]]

    @segment.setter
    def segment(self, new):
        if not isinstance(new, Segment):
            raise TypeError("Can only set Residue segment to Segment, not {}"
                            "".format(type(new)))
        self.universe._topology.tt.move_residue(self.ix, new.segindex)


class Segment(ComponentBase):
    """:class:`Segment` base class.

    This class is used by a :class:`~MDAnalysis.core.universe.Universe` for
    generating its Topology-specific :class:`Segment` class. All the
    :class:`~MDAnalysis.core.topologyattrs.TopologyAttr` components are obtained
    from :class:`ComponentBase`, so this class only includes ad-hoc methods
    specific to :class:`Segments<Segment>`.

    .. deprecated:: 0.16.2
       *Instant selectors* of :class:`Segments<Segment>` will be removed in the
       1.0 release. See :ref:`Instant selectors <instance-selectors>` for
       details and alternatives.
    """
    def __repr__(self):
        me = '<Segment'
        if hasattr(self, 'segid'):
            me += ' {}'.format(self.segid)
        return me + '>'

    @property
    def atoms(self):
        """An :class:`AtomGroup` of :class:`Atoms<Atom>` present in this
        :class:`Segment`.
        """
        ag = self.universe.atoms[self.universe._topology.indices[self][0]]
        ag._cache['isunique'] = True
        ag._cache['unique'] = ag
        return ag

    @property
    def residues(self):
        """A :class:`ResidueGroup` of :class:`Residues<Residue>` present in this
        :class:`Segment`.
        """
        rg = self.universe.residues[self.universe._topology.resindices[self][0]]
        rg._cache['isunique'] = True
        rg._cache['unique'] = rg
        return rg

    def __getattr__(self, attr):
        # DEPRECATED in 0.16.2
        # REMOVE in 1.0
        #
        # Segment.r1 access
        if attr.startswith('r') and attr[1:].isdigit():
            resnum = int(attr[1:])
            rg = self.residues[resnum - 1]  # convert to 0 based
            warnings.warn("Instant selectors Segment.r<N> will be removed in "
                          "1.0. Use Segment.residues[N-1] instead.",
                          DeprecationWarning)
            return rg
        # Resname accesss
        if hasattr(self.residues, 'resnames'):
            try:
                return self.residues._get_named_residue(attr)
            except selection.SelectionError:
                pass
        raise AttributeError("{cls} has no attribute {attr}"
                             "".format(cls=self.__class__.__name__, attr=attr))

# Accessing these attrs doesn't trigger an update. The class and instance
# methods of UpdatingAtomGroup that are used during __init__ must all be
# here, otherwise we get __getattribute__ infinite loops.
_UAG_SHORTCUT_ATTRS = {
    # Class information of the UAG
    "__class__", "_derived_class",
    # Metadata of the UAG
    "_base_group", "_selections", "_lastupdate",
    "level", "_u", "universe",
    # Methods of the UAG
    "_ensure_updated",
    "is_uptodate",
    "update_selection",
}

class UpdatingAtomGroup(AtomGroup):
    """:class:`AtomGroup` subclass that dynamically updates its selected atoms.

    Accessing any attribute/method of an :class:`UpdatingAtomGroup` instance
    triggers a check for the last frame the group was updated. If the last
    frame matches the current trajectory frame, the attribute is returned
    normally; otherwise the group is updated (the stored selections are
    re-applied), and only then is the attribute returned.


    .. versionadded:: 0.16.0
    """
    # WARNING: This class has __getattribute__ and __getattr__ methods (the
    # latter inherited from AtomGroup). Because of this bugs introduced in the
    # class that cause an AttributeError may be very hard to diagnose and
    # debug: the most obvious symptom is an infinite loop going through both
    # __getattribute__ and __getattr__, and a solution might be to add said
    # attribute to _UAG_SHORTCUT_ATTRS.

    def __init__(self, base_group, selections, strings):
        """

        Parameters
        ----------
        base_group : :class:`AtomGroup`
            group on which *selections* are to be applied.
        selections : a tuple of :class:`~MDAnalysis.core.selection.Selection`
            instances selections ready to be applied to *base_group*.
        """
        # Because we're implementing __getattribute__, which needs _u for
        # its check, no self.attribute access can be made before this line
        self._u = base_group.universe
        self._selections = selections
        self._selection_strings = strings
        self._base_group = base_group
        self._lastupdate = None
        self._derived_class = base_group._derived_class
        if self._selections:
            # Allows the creation of a cheap placeholder UpdatingAtomGroup
            # by passing an empty selection tuple.
            self._ensure_updated()

    def update_selection(self):
        """
        Forces the reevaluation and application of the group's selection(s).

        This method is triggered automatically when accessing attributes, if
        the last update occurred under a different trajectory frame.
        """
        bg = self._base_group
        sels = self._selections
        if sels:
            # As with select_atoms, we select the first sel and then sum to it.
            ix = sum([sel.apply(bg) for sel in sels[1:]],
                     sels[0].apply(bg)).ix
        else:
            ix = np.array([], dtype=np.intp)
        # Run back through AtomGroup init with this information to remake
        # ourselves
        super(UpdatingAtomGroup, self).__init__(ix, self.universe)
        self.is_uptodate = True

    @property
    def is_uptodate(self):
        """
        Checks whether the selection needs updating based on frame number only.

        Modifications to the coordinate data that render selections stale are
        not caught, and in those cases :attr:`is_uptodate` may return an
        erroneous value.

        Returns
        -------
        bool
            ``True`` if the group's selection is up-to-date, ``False``
            otherwise.
        """
        try:
            return self.universe.trajectory.frame == self._lastupdate
        except AttributeError: # self.universe has no trajectory
            return self._lastupdate == -1

    @is_uptodate.setter
    def is_uptodate(self, value):
        if value:
            try:
                self._lastupdate = self.universe.trajectory.frame
            except AttributeError: # self.universe has no trajectory
                self._lastupdate = -1
        else:
            # This always marks the selection as outdated
            self._lastupdate = None

    def _ensure_updated(self):
        """
        Checks whether the selection needs updating and updates it if needed.

        Returns
        -------
        bool
            ``True`` if the group was already up-to-date, ``False`` otherwise.
        """
        status = self.is_uptodate
        if not status:
            self.update_selection()
        return status

    def __getattribute__(self, name):
        # ALL attribute access goes through here
        # If the requested attribute is public (not starting with '_') and
        # isn't in the shortcut list, update ourselves
        if not (name.startswith('_') or name in _UAG_SHORTCUT_ATTRS):
            self._ensure_updated()
        # Going via object.__getattribute__ then bypasses this check stage
        return object.__getattribute__(self, name)

    def __reduce__(self):
        # strategy for unpickling is:
        # - unpickle base group
        # - recreate UAG as created through select_atoms (basegroup and selstrs)
        # even if base_group is a UAG this will work through recursion
        return (_unpickle_uag,
                (self._base_group.__reduce__(), self._selections,
                 self._selection_strings))

    def __repr__(self):
        basestr = super(UpdatingAtomGroup, self).__repr__()
        if not self._selection_strings:
            return basestr
        sels = "'{}'".format("' + '".join(self._selection_strings))
        # Cheap comparison. Might fail for corner cases but this is
        # mostly cosmetic.
        if self._base_group is self.universe.atoms:
            basegrp = "the entire Universe."
        else:
            basegrp = "another AtomGroup."
        # With a shorthand to conditionally append the 's' in 'selections'.
        return "{}, with selection{} {} on {}>".format(basestr[:-1],
                    "s"[len(self._selection_strings)==1:], sels, basegrp)

    @property
    def atoms(self):
        """Get a *static* :class:`AtomGroup` identical to the group of currently
        selected :class:`Atoms<Atom>` in the :class:`UpdatingAtomGroup`.


        By returning a *static* :class:`AtomGroup` it becomes possible to
        compare the contents of the group *between* trajectory frames. See the
        Example below.


        Note
        ----
        The :attr:`atoms` attribute of an :class:`UpdatingAtomGroup` behaves
        differently from :attr:`AtomGroup.atoms`: the latter returns the
        :class:`AtomGroup` itself whereas the former returns a
        :class:`AtomGroup` and not an :class:`UpdatingAtomGroup` (for this, use
        :meth:`UpdatingAtomGroup.copy`).


        Example
        -------
        The static :attr:`atoms` allows comparison of groups of atoms between
        frames. For example, track water molecules that move in and out of a
        solvation shell of a protein::

          u = mda.Universe(TPR, XTC)
          water_shell = u.select_atoms("name OW and around 3.5 protein", updating=True)
          water_shell_prev = water_shell.atoms

          for ts in u.trajectory:
              exchanged = water_shell - water_shell_prev

              print(ts.time, "waters in shell =", water_shell.n_residues)
              print(ts.time, "waters that exchanged = ", exchanged.n_residues)
              print(ts.time, "waters that remained bound = ",
                    water_shell.n_residues - exchanged.n_residues)

              water_shell_prev = water_shell.atoms

        By remembering the atoms of the current time step in
        `water_shell_prev`, it becomes possible to use the :meth:`subtraction
        of AtomGroups<AtomGroup.subtract>` to find the water molecules that
        changed.


        See Also
        --------
        copy : return a true copy of the :class:`UpdatingAtomGroup`

        """
        return self[:]

    def copy(self):
        """Get another :class:`UpdatingAtomGroup` identical to this one.


        .. versionadded:: 0.19.0
        """
        return UpdatingAtomGroup(self._base_group, self._selections,
                                 self._selection_strings)


# Define relationships between these classes
# with Level objects
_Level = namedtuple('Level', ['name', 'singular', 'plural'])
ATOMLEVEL = _Level('atom', Atom, AtomGroup)
RESIDUELEVEL = _Level('residue', Residue, ResidueGroup)
SEGMENTLEVEL = _Level('segment', Segment, SegmentGroup)

Atom.level = ATOMLEVEL
AtomGroup.level = ATOMLEVEL
Residue.level = RESIDUELEVEL
ResidueGroup.level = RESIDUELEVEL
Segment.level = SEGMENTLEVEL
SegmentGroup.level = SEGMENTLEVEL

def requires(*attrs):
    """Decorator to check if all :class:`AtomGroup` arguments have certain
    attributes

    Example
    -------
    When used to wrap a function, will check all :class:`AtomGroup` arguments
    for the listed requirements

    @requires('masses', 'charges')
    def mass_times_charge(atomgroup):
        return atomgroup.masses * atomgroup.charges
    """
    def require_dec(func):
        @functools.wraps(func)
        def check_args(*args, **kwargs):
            for a in args:  # for each argument
                if isinstance(a, AtomGroup):
                    # Make list of missing attributes
                    missing = [attr for attr in attrs
                               if not hasattr(a, attr)]
                    if missing:
                        raise NoDataError(
                            "{funcname} failed. "
                            "AtomGroup is missing the following required "
                            "attributes: {attrs}".format(
                                funcname=func.__name__,
                                attrs=', '.join(missing)))
            return func(*args, **kwargs)
        return check_args
    return require_dec
