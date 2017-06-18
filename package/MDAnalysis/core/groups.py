# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.mdanalysis.org
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

The :class:`~MDAnalysis.core.universe.Universe` instance contains all
the particles in the system (which MDAnalysis calls
:class:`Atom`). Groups of atoms are handled as :class:`AtomGroup`
instances. The :class:`AtomGroup` is probably the most important
object in MDAnalysis because virtually everything can be accessed
through it. `AtomGroup` instances can be easily created (e.g., from a
:meth:`AtomGroup.select_atoms` selection or just by slicing).

For convenience, chemically meaningful groups of atoms such as a
:class:`Residue` or a :class:`Segment` (typically a whole molecule or
all of the solvent) also exist as containers, as well as groups of
these units ((:class:`ResidueGroup`, :class:`SegmentGroup`).


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

Each of the above classes has a level attribute.  This can be used to
verify that two objects are of the same level, or to access a particular
class::

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

import MDAnalysis
from .. import _ANCHOR_UNIVERSES
from ..lib import util
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
    """Make a fresh copy of all Classes

    Returns
    -------
    Two dictionaries. One with a set of :class:`_TopologyAttrContainer` classes
    to serve as bases for universe-specific MDA container classes. Another with
    the final merged versions of those classes. The classes themselves are used
    as hashing keys.

    """
    bases = {}
    classes = {}
    groups = (AtomGroup, ResidueGroup, SegmentGroup)
    components = (Atom, Residue, Segment)

    # The 'GBase' middle man is needed so that a single topologyattr
    #  patching applies automatically to all groups.
    GBase = bases[GroupBase] = _TopologyAttrContainer._subclass()
    for cls in groups:
        bases[cls] = GBase._subclass()

    # In the current Group-centered topology scheme no attributes apply only
    #  to ComponentBase, so no need to have a 'CB' middle man.
    #CBase = _TopologyAttrContainer(singular=True)
    for cls in components:
        bases[cls] = _TopologyAttrContainer._subclass(singular=True)

    # Initializes the class cache.
    for cls in groups + components:
        classes[cls] = bases[cls]._mix(cls)

    return bases, classes


class _TopologyAttrContainer(object):
    """Class factory for receiving sets of :class:`TopologyAttr` objects.

    :class:`_TopologyAttrContainer` is a convenience class to encapsulate the
    functions that deal with:
    * the import and namespace transplant of :class:`TopologyAttr` objects;
    * the copying (subclassing) of itself to create distinct bases for the
      different container classes (:class:`AtomGroup`, :class:`ResidueGroup`,
      :class:`SegmentGroup`, :class:`Atom`, :class:`Residue`, :class:`Segment`,
      and subclasses thereof);
    * the mixing (subclassing and co-inheritance) with the container classes.
      The mixed subclasses become the final container classes specific to each
      :class:`Universe`.
    """
    _singular = False

    @classmethod
    def _subclass(cls, singular=None):
        """Factory method returning :class:`_TopologyAttrContainer` subclasses.

        Parameters
        ----------
        singular : bool
            The :attr:`_singular` of the returned class will be set to
            *singular*. It controls the type of :class:`TopologyAttr` addition.

        Returns
        -------
        type
            A subclass of :class:`_TopologyAttrContainer`, with the same name.
        """
        if singular is not None:
            return type(cls.__name__, (cls,), {'_singular': bool(singular)})
        else:
            return type(cls.__name__, (cls,), {})

    @classmethod
    def _mix(cls, other):
        """Creates a subclass with ourselves and another class as parents.

        Classes mixed at this point override :meth:`__new__`, causing further
        instantiations to shortcut to :meth:`~object.__new__` (skipping the
        cache-fetch process for :class:`_MutableBase` subclasses).

        The new class will have an attribute `_derived_class` added, pointing
        to itself. This pointer instructs which class to use when
        slicing/adding instances of the new class. At initialization time the
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
        """Add attr into the namespace for this class

        Parameters
        ----------
        attr : A TopologyAttr object
        """
        def getter(self):
            return attr.__getitem__(self)

        def setter(self, values):
            return attr.__setitem__(self, values)

        if cls._singular:
            setattr(cls, attr.singular,
                    property(getter, setter, None, attr.singledoc))
        else:
            setattr(cls, attr.attrname,
                    property(getter, setter, None, attr.groupdoc))


class _MutableBase(object):
    """
    Base class that merges appropriate :class:`_TopologyAttrContainer` classes.

    Implements :meth:`__new__`. In it the instantiating class is fetched from
    :attr:`Universe._classes`. If there is a cache miss, a merged class is made
    with a base from :attr:`Universe._class_bases` and cached.

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
                # deprecated AtomGroup init method..
                u = args[0][0].universe
            except (IndexError, AttributeError):
                # Let's be generic and get the first argument that's either a
                # Universe, a Group, or a Component, and go from there.
                # This is where the UpdatingAtomGroup args get matched.
                for arg in args+tuple(kwargs.values()):
                    if isinstance(arg, (MDAnalysis.Universe, GroupBase,
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
    """Base class from which a Universe's Group class is built.

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
        if len(args) == 1:
            warnings.warn("Using deprecated init method for Group. "
                          "In the future use `Group(indices, universe)`. "
                          "This init method will be removed in version 1.0.",
                          DeprecationWarning)
            # list of atoms/res/segs, old init method
            ix = [at.ix for at in args[0]]
            u = args[0][0].universe
        else:
            # current/new init method
            ix, u = args

        # indices for the objects I hold
        self._ix = np.asarray(ix, dtype=np.intp)
        self._u = u
        self._cache = dict()

    def __hash__(self):
        return hash((self._u, self.__class__, tuple(self.ix.tolist())))

    def __len__(self):
        return len(self._ix)

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
        to different universe cannot be compared.
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

        For a Group, ix_array is the same as ix. This method gives a
        consistent API between components and groups.

        See Also
        --------
        ix
        """
        return self._ix

    @property
    def dimensions(self):
        return self.universe.trajectory.ts.dimensions

    @dimensions.setter
    def dimensions(self, dimensions):
        self.universe.trajectory.ts.dimensions = dimensions

    def center(self, weights, pbc=None):
        """Calculate center of group given some weights

        Parameters
        ----------
        weights : array_like
            weights to be used
        pbc : boolean, optional
            ``True``: Move all atoms within the primary unit cell
            before calculation [``False``]

        Returns
        -------
        center : ndarray
            weighted center of group

        Examples
        --------

        To find the charge weighted center of a given Atomgroup::

            >>> sel = u.select_atoms('prop mass > 4.0')
            >>> sel.center(sel.charges)


        Notes
        -----
        If the :class:`MDAnalysis.core.flags` flag *use_pbc* is set to
        ``True`` then the `pbc` keyword is used by default.

        """
        atoms = self.atoms
        if pbc is None:
            pbc = flags['use_pbc']
        if pbc:
            xyz = atoms.pack_into_box(inplace=False)
        else:
            xyz = atoms.positions

        return np.average(xyz, weights=weights, axis=0)

    def center_of_geometry(self, pbc=None):
        """Center of geometry (also known as centroid) of the selection.

        Parameters
        ----------
        pbc : boolean, optional
            ``True``: Move all atoms within the primary unit cell
            before calculation [``False``]

        Returns
        -------
        center : ndarray
            geometric center of group

        Notes
        -----
        If the :class:`MDAnalysis.core.flags` flag *use_pbc* is set to
        ``True`` then the `pbc` keyword is used by default.


        .. versionchanged:: 0.8 Added `pbc` keyword
        """
        return self.center(None, pbc=pbc)

    centroid = center_of_geometry

    def bbox(self, **kwargs):
        """Return the bounding box of the selection.

        The lengths A,B,C of the orthorhombic enclosing box are ::

          L = AtomGroup.bbox()
          A,B,C = L[1] - L[0]

        Parameters
        ----------
        pbc : bool, optional
            If ``True``, move all atoms within the primary unit cell before
            calculation. [``False``]

        Returns
        -------
         corners : array
            2x3 array giving corners of bounding box as
            [[xmin, ymin, zmin], [xmax, ymax, zmax]].

        Note
        ----
        The :class:`MDAnalysis.core.flags` flag *use_pbc* when set to
        ``True`` allows the *pbc* flag to be used by default.


        .. versionadded:: 0.7.2
        .. versionchanged:: 0.8 Added *pbc* keyword
        """
        atomgroup = self.atoms
        pbc = kwargs.pop('pbc', MDAnalysis.core.flags['use_pbc'])

        if pbc:
            x = atomgroup.pack_into_box(inplace=False)
        else:
            x = atomgroup.positions

        return np.array([x.min(axis=0), x.max(axis=0)])

    def bsphere(self, **kwargs):
        """Return the bounding sphere of the selection.

        The sphere is calculated relative to the centre of geometry.

        Parameters
        ----------
        pbc : bool, optional
            If ``True``, move all atoms within the primary unit cell before
            calculation. [``False``]

        Returns
        -------
        R : float
            Radius of bounding sphere.
        center : array
            Coordinates of sphere center as ``[xcen,ycen,zcen]``.

        Note
        ----
        The :class:`MDAnalysis.core.flags` flag *use_pbc* when set to
        ``True`` allows the *pbc* flag to be used by default.


        .. versionadded:: 0.7.3
        .. versionchanged:: 0.8 Added *pbc* keyword
        """
        atomgroup = self.atoms
        pbc = kwargs.pop('pbc', MDAnalysis.core.flags['use_pbc'])

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

        Parameters
        ----------
        M : array
            4x4 matrix, with the rotation in ``R = M[:3,:3]`` and the
            translation in ``t = M[:3,3]``.

        Returns
        -------
        self

        See Also
        --------
        MDAnalysis.lib.transformations : module of all coordinate transforms

        Notes
        -----
        The rotation :math:`\mathsf{R}` is applied before the translation
        :math:`\mathbf{t}`:

        .. math::

           \mathbf{x}' = \mathsf{R}\mathbf{x} + \mathbf{t}

        """
        R = M[:3, :3]
        t = M[:3, 3]
        return self.rotate(R, [0, 0, 0]).translate(t)

    def translate(self, t):
        r"""Apply translation vector `t` to the selection's coordinates.

        Atom coordinates are translated in-place.

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

        Parameters
        ----------
        R : array_like
            3x3 rotation matrix to use for applying rotation.
        point : array_like, optional
            Center of rotation

        Returns
        -------
        self : AtomGroup

        Notes
        -----
        By default rotates around center of origin ``point=(0, 0, 0)``. To
        rotate around center of geometry of the atomgroup use ``ag.rotate(R,
        point=ag.centroid)``.

        See Also
        --------
        rotateby : rotate around given axis and angle
        MDAnalysis.lib.transformations : module of all coordinate transforms

        """
        R = np.asarray(R)
        point = np.asarray(point)

        # changes the coordinates (in place)
        self.translate(-point)
        x = self.atoms.unique.universe.trajectory.ts.positions
        idx = self.atoms.unique.indices
        x[idx] = np.dot(x[idx], R.T)
        self.translate(point)

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
        self : AtomGroup

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
        MDAnalysis.lib.transformations.rotation_matrix : calculate :math:`\mathsf{R}`

        """
        alpha = np.radians(angle)
        axis = np.asarray(axis)
        point = np.asarray(point) if point is not None else self.centroid()
        M = transformations.rotation_matrix(alpha, axis, point=point)
        return self.transform(M)

    def pack_into_box(self, box=None, inplace=True):
        r"""Shift all atoms in this group to be within the primary unit cell.

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
        coords : array
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
        (in the MDAnalysis standard format ``[Lx, Ly, Lz, alpha, beta,
        gamma]``).

        Works with either orthogonal or triclinic box types.


        .. versionadded:: 0.8

        """
        atomgroup = self.atoms.unique
        if box is None:  # Try and auto detect box dimensions
            box = atomgroup.dimensions  # Can accept any box

        if box.shape == (3, 3):
            # for a vector representation, diagonal cannot be zero
            if (box.diagonal() == 0.0).any():
                raise ValueError("One or more box dimensions is zero."
                                 "  You can specify a boxsize with 'box ='")
        else:
            if (box == 0).any():  # Check that a box dimension isn't zero
                raise ValueError("One or more box dimensions is zero."
                                 "  You can specify a boxsize with 'box='")

        coords = atomgroup.universe.coord.positions[atomgroup.indices]
        if not inplace:
            return distances.apply_PBC(coords, box)

        atomgroup.universe.coord.positions[atomgroup.indices] = distances.apply_PBC(coords, box)

        return atomgroup.universe.coord.positions[atomgroup.indices]

    def wrap(self, compound="atoms", center="com", box=None):
        """Shift the contents of this Group back into the unit cell.

        This is a more powerful version of :meth:`pack_into_box`, allowing
        groups of atoms to be kept together through the process.

        Parameters
        ----------
        compound : {'atoms', 'group', 'residues', 'segments', 'fragments'}
            The group which will be kept together through the shifting process.
        center : {'com', 'cog'}
            How to define the center of a given group of atoms.
        box : array
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
        This might however mean that all atoms from the compound are not
        inside the unit cell, but rather the center of the compound is.

        `center` allows the definition of the center of each group to be
        specified. This can be either 'com' for center of mass, or 'cog' for
        center of geometry.

        `box` allows a unit cell to be given for the transformation. If not
        specified, an the dimensions information from the current Timestep will
        be used.

        .. note::
           wrap with all default keywords is identical to :meth:`pack_into_box`


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
            raise ValueError("Unrecognised compound definition: {0}"
                             "Please use one of 'group' 'residues' 'segments'"
                             "or 'fragments'".format(compound))

        # TODO: ADD TRY-EXCEPT FOR MASSES PRESENCE
        if center.lower() in ('com', 'centerofmass'):
            centers = np.vstack([o.atoms.center_of_mass() for o in objects])
        elif center.lower() in ('cog', 'centroid', 'centerofgeometry'):
            centers = np.vstack([o.atoms.center_of_geometry() for o in objects])
        else:
            raise ValueError("Unrecognised center definition: {0}"
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

    def groupby(self, topattr):
        """Group together items in this group according to values of *topattr*

        Parameters
        ----------
        topattr: str
           Topology attribute to group components by.

        Returns
        -------
        dict
            Unique values of the topology attribute as keys, Groups as values.

        Example
        -------
        To group atoms with the same mass together::

          >>> ag.groupby('masses')
          {12.010999999999999: <AtomGroup with 462 atoms>,
          14.007: <AtomGroup with 116 atoms>,
          15.999000000000001: <AtomGroup with 134 atoms>}

        .. versionadded:: 0.16.0
        """
        ta = getattr(self, topattr)
        return {i: self[ta == i] for i in set(ta)}

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
        Bool
            True if the two Groups do not have common elements

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
        Bool
            True if this Group is a subset of the other one

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
        Bool
            True if this Group is a strict subset of the other one

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
        Bool
            True if this Group is a subset of the other one

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
        Bool
            True if this Group is a strict superset of the other one

        .. versionadded:: 0.16
        """
        return self.issuperset(other) and not self == other


class AtomGroup(GroupBase):
    """A group of atoms.

    An :class:`AtomGroup` is an ordered collection of atoms. Typically, an
    :class:`AtomGroup` is generated from a selection, or by indexing/slicing
    the :class:`AtomGroup` of all atoms in the :class:`Universe` at
    :attr:`MDAnalysis.core.universe.Universe.atoms`.

    An AtomGroup can be indexed and sliced like a list::

        ag[0], ag[-1]

    will return the first and the last :class:`Atom` in the group whereas the
    slice ::

        ag[0:6:2]

    returns an AtomGroup of every second element, corresponding to indices 0,
    2, and 4.

    It also supports "advanced slicing" when the argument is a
    :class:`numpy.ndarray` or a :class:`list`::

        aslice = [0, 3, -1, 10, 3]
        ag[aslice]

    will return a new AtomGroup of atoms with those indices in the old
    AtomGroup.

    .. note::

        AtomGroups originating from a selection are sorted and
        duplicate elements are removed. This is not true for AtomGroups
        produced by slicing. Thus slicing can be used when the order of
        atoms is crucial (for instance, in order to define angles or
        dihedrals).

    AtomGroups can be compared and combined using group operators. For
    instance, AtomGroups can be concatenated using `+` or :meth:`concatenate`::

        ag_concat = ag1 + ag2  # or ag_concat = ag1.concatenate(ag2)

    When groups are concatenated, the order of the atoms is conserved. If atoms
    appear several times in one of the groups, the duplicates are kept in the
    resulting group. On the contrary to :meth:`concatenate`, :meth:`union`
    treats the AtomGroups as sets, duplicates are removed from the resulting
    group, and atoms are ordered. The `|` operator is synomymous to
    :meth:`union`::

        ag_union = ag1 | ag2  # or ag_union = ag1.union(ag2)

    The opposite operation to :meth:`concatenate` is :meth:`subtract`. This
    method creates a new group with all the atoms of the group that are not in
    a given other group; the order of the atoms is kept, so as duplicates.
    :meth:`difference` is the set version of :meth:`subtract`. The resulting
    group is sorted and deduplicated.

    All set methods are listed in the table below. These methods treat the
    groups as sorted and deduplicated sets of atoms.

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

    The following methods keep the order of the atoms, and keep duplicated
    atoms.

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

    The `in` operator allows to test if an :class:`Atom` is in the AtomGroup.

    AtomGroup instances are always bound to a
    :class:`MDAnalysis.core.universe.Universe`. They cannot exist in isolation.


    .. rubric:: Deprecated functionality

    *Instant selectors* will be removed in the 1.0 release.  See issue `#1377
    <https://github.com/MDAnalysis/mdanalysis/issues/1377>`_ for more details.

    Atoms can also be accessed in a Pythonic fashion by using the atom name as
    an attribute. For instance, ::

        ag.CA

    will provide a :class:`AtomGroup` of all CA atoms in the
    group. These *instant selector* attributes are auto-generated for
    each atom name encountered in the group.

    Notes
    -----
    The name-attribute instant selector access to atoms is mainly
    meant for quick interactive work. Thus it either returns a
    single :class:`Atom` if there is only one matching atom, *or* a
    new :class:`AtomGroup` for multiple matches.  This makes it
    difficult to use the feature consistently in scripts.


    See Also
    --------
    :class:`MDAnalysis.core.universe.Universe`


    .. deprecated:: 0.16.2
       *Instant selectors* of AtomGroup will be removed in the 1.0 release.
       See :ref:`Instant selectors <instance-selectors>` for details and alternatives.

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
        """Get another AtomGroup identical to this one."""
        return self.universe.atoms[self.ix]

    @property
    def n_atoms(self):
        """Number of atoms in AtomGroup.

        Equivalent to ``len(self)``."""
        return len(self)

    @property
    def residues(self):
        """Get sorted :class:`ResidueGroup` of the (unique) residues
        represented in the AtomGroup."""
        return self.universe.residues[np.unique(self.resindices)]

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
        """Number of unique residues represented in the AtomGroup.

        Equivalent to ``len(self.residues)``.

        """
        return len(self.residues)

    @property
    def segments(self):
        """Get sorted :class:`SegmentGroup` of the (unique) segments
        represented in the AtomGroup."""
        return self.universe.segments[np.unique(self.segindices)]

    @segments.setter
    def segments(self, new):
        raise NotImplementedError("Cannot assign Segments to AtomGroup. "
                                  "Segments are assigned to Residues")

    @property
    def n_segments(self):
        """Number of unique segments represented in the AtomGroup.

        Equivalent to ``len(self.segments)``.

        """
        return len(self.segments)

    @property
    def unique(self):
        """Return an AtomGroup containing sorted and unique atoms only.

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

        .. versionadded:: 0.16.0
        """
        return self.universe.atoms[np.unique(self.ix)]

    @property
    def positions(self):
        """Coordinates of the atoms in the AtomGroup.

        The positions can be changed by assigning an array of the appropriate
        shape, i.e. either Nx3 to assign individual coordinates or 3, to assign
        the *same* coordinate to all atoms (e.g. ``ag.positions =
        array([0,0,0])`` will move all particles to the origin).

        .. note:: Changing the position is not reflected in any files;
                  reading any frame from the trajectory will replace
                  the change with that from the file *except* if the
                  trajectory is held in memory, e.g., when the
                  :class:`~MDAnalysis.core.universe.Universe.transfer_to_memory`
                  method was used.

        """
        return self.universe.trajectory.ts.positions[self.ix]

    @positions.setter
    def positions(self, values):
        ts = self.universe.trajectory.ts
        ts.positions[self.ix, :] = values

    @property
    def velocities(self):
        """Velocities of the atoms in the AtomGroup.

        The velocities can be changed by assigning an array of the appropriate
        shape, i.e. either Nx3 to assign individual velocities or 3 to assign
        the *same* velocity to all atoms (e.g. ``ag.velocity = array([0,0,0])``
        will give all particles zero velocity).

        Raises a :exc:`NoDataError` if the underlying
        :class:`~MDAnalysis.coordinates.base.Timestep` does not contain
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
        """Forces on each atom in the AtomGroup.

        The velocities can be changed by assigning an array of the appropriate
        shape, i.e. either Nx3 to assign individual velocities or 3 to assign
        the *same* velocity to all atoms (e.g. ``ag.velocity = array([0,0,0])``
        will give all particles zero velocity).

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
        typically happens when the underlying trajectory frame
        changes).

        It is not possible to assign a new
        :class:`~MDAnalysis.coordinates.base.Timestep` to the
        :attr:`AtomGroup.ts` attribute; change attributes of the object.
        """
        trj_ts = self.universe.trajectory.ts  # original time step

        return trj_ts.copy_slice(self.indices)

    # As with universe.select_atoms, needing to fish out specific kwargs
    # (namely, 'updating') doesn't allow a very clean signature.
    def select_atoms(self, sel, *othersel, **selgroups):
        """Select atoms using a selection string.

        Returns an :class:`AtomGroup` with atoms sorted according to
        their index in the topology (this is to ensure that there
        are not any duplicates, which can happen with complicated
        selections).

        Existing :class:`AtomGroup` objects can be passed as named arguments,
        which will then be available to the selection parser.

        Subselections can be grouped with parentheses.

        Selections can be set to update automatically on frame change, by
        setting the `updating` named argument to `True`.

        Examples
        --------

           >>> sel = universe.select_atoms("segid DMPC and not ( name H* or name O* )")
           >>> sel
           <AtomGroup with 3420 atoms>

           >>> universe.select_atoms("around 10 group notHO", notHO=sel)
           <AtomGroup with 1250 atoms>

        Notes
        -----

        If exact ordering of atoms is required (for instance, for
        :meth:`~AtomGroup.angle` or :meth:`~AtomGroup.dihedral`
        calculations) then one supplies selections *separately* in the
        required order. Also, when multiple :class:`AtomGroup`
        instances are concatenated with the ``+`` operator then the
        order of :class:`Atom` instances is preserved and duplicates
        are not removed.


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

        **Boolean**

            not
                all atoms not in the selection, e.g. ``not protein`` selects
                all atoms that aren't part of a protein

            and, or
                combine two selections according to the rules of boolean
                algebra, e.g. ``protein and not (resname ALA or resname LYS)``
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

        **Connectivity**

            byres *selection*
                selects all atoms that are in the same segment and residue as
                selection, e.g. specify the subselection after the byres keyword
            bonded *selection*
                selects all atoms that are bonded to selection
                eg: ``select name H bonded name O`` selects only hydrogens
                bonded to oxygens

        **Index**

            bynum *index-range*
                selects all atoms within a range of (1-based) inclusive indices,
                e.g. ``bynum 1`` selects the first atom in the universe;
                ``bynum 5:10`` selects atoms 5 through 10 inclusive. All atoms
                in the :class:`MDAnalysis.Universe` are consecutively numbered,
                and the index runs from 1 up to the total number of atoms.

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
                :class:`Universe`.

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
                :class:`~MDAnalysis.core.universe.Universe`, ``global`` is ignored.

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


        .. versionchanged:: 0.7.4
           Added *resnum* selection.
        .. versionchanged:: 0.8.1
           Added *group* and *fullgroup* selections.
        .. deprecated:: 0.11
           The use of ``fullgroup`` has been deprecated in favor of the equivalent
           ``global group``.
        .. versionchanged:: 0.13.0
           Added *bonded* selection
        .. versionchanged:: 0.16.0
           Resid selection now takes icodes into account where present.
        .. versionadded:: 0.16.0
           Updating selections now possible by setting the ``updating`` argument.

        """
        updating = selgroups.pop('updating', False)
        sel_strs = (sel,) + othersel
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
        """Split AtomGroup into a list of atomgroups by `level`.

        Parameters
        ----------
        level : {'atom', 'residue', 'segment'}


        .. versionadded:: 0.9.0
        """
        accessors = {'segment': 'segindices',
                     'residue': 'resindices'}

        if level == "atom":
            return [self.universe.atoms[[a.ix]] for a in self]

        # higher level groupings
        try:
            levelindices = getattr(self, accessors[level])
        except KeyError:
            raise ValueError("level = '{0}' not supported, "
                             "must be one of {1}".format(level,
                                                         accessors.keys()))

        return [self[levelindices == index] for index in
                np.unique(levelindices)]

    def guess_bonds(self, vdwradii=None):
        """Guess bonds that exist within this AtomGroup and add to Universe

        Parameters
        ----------
        vdwradii : dict, optional
          Dict relating atom type: vdw radii


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
        """This AtomGroup represented as a Bond object

        Returns
        -------
          A :class:`MDAnalysis.core.topologyobjects.Bond` object

        Raises
        ------
          `ValueError` if the AtomGroup is not length 2


        .. versionadded:: 0.11.0
        """
        if len(self) != 2:
            raise ValueError(
                "bond only makes sense for a group with exactly 2 atoms")
        return topologyobjects.Bond(self.ix, self.universe)

    @property
    def angle(self):
        """This AtomGroup represented as an Angle object

        Returns
        -------
          A :class:`MDAnalysis.core.topologyobjects.Angle` object

        Raises
        ------
          `ValueError` if the AtomGroup is not length 3


        .. versionadded:: 0.11.0
        """
        if len(self) != 3:
            raise ValueError(
                "angle only makes sense for a group with exactly 3 atoms")
        return topologyobjects.Angle(self.ix, self.universe)

    @property
    def dihedral(self):
        """This AtomGroup represented as a Dihedral object

        Returns
        -------
          A :class:`MDAnalysis.core.topologyobjects.Dihedral` object

        Raises
        ------
          `ValueError` if the AtomGroup is not length 4


        .. versionadded:: 0.11.0
        """
        if len(self) != 4:
            raise ValueError(
                "dihedral only makes sense for a group with exactly 4 atoms")
        return topologyobjects.Dihedral(self.ix, self.universe)

    @property
    def improper(self):
        """This AtomGroup represented as an ImproperDihedral object

        Returns
        -------
          A :class:`MDAnalysis.core.topologyobjects.ImproperDihedral` object

        Raises
        ------
          `ValueError` if the AtomGroup is not length 4


        .. versionadded:: 0.11.0
        """
        if len(self) != 4:
            raise ValueError(
                "improper only makes sense for a group with exactly 4 atoms")
        return topologyobjects.ImproperDihedral(self.ix, self.universe)

    def write(self, filename=None, file_format="PDB",
              filenamefmt="{trjname}_{frame}", **kwargs):
        """Write `AtomGroup` to a file.

        The output can either be a coordinate file or a selection, depending on
        the `format`. Only single-frame coordinate files are supported. If you
        need to write out a trajectory, see :mod:`MDAnalysis.coordinates`.

        Parameters
        ----------
        filename : str, optional
            ``None``: create TRJNAME_FRAME.FORMAT from filenamefmt [``None``]
        file_format : str, optional
            PDB, CRD, GRO, VMD (tcl), PyMol (pml), Gromacs (ndx) CHARMM (str)
            Jmol (spt); case-insensitive and can also be supplied as the
            filename extension [PDB]
        filenamefmt : str, optional
            format string for default filename; use substitution tokens
            'trjname' and 'frame' ["%(trjname)s_%(frame)d"]
        bonds : str, optional
            how to handle bond information, especially relevant for PDBs.
            ``"conect"``: write only the CONECT records defined in the original
            file. ``"all"``: write out all bonds, both the original defined and
            those guessed by MDAnalysis. ``None``: do not write out bonds.
            Default os ``"conect"``.


        .. versionchanged:: 0.9.0
           Merged with write_selection.  This method can now write both
           selections out.

        """
        # check that AtomGroup actually has any atoms (Issue #434)
        if len(self.atoms) == 0:
            raise IndexError("Cannot write an AtomGroup with 0 atoms")

        trj = self.universe.trajectory  # unified trajectory API

        if trj.n_frames == 1:
            kwargs.setdefault("multiframe", False)

        if filename is None:
            trjname, ext = os.path.splitext(os.path.basename(trj.filename))
            filename = filenamefmt.format(trjname=trjname, frame=trj.frame)
        filename = util.filename(filename, ext=file_format.lower(), keep=True)

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

            multiframe = kwargs.pop('multiframe', None)

            writer = get_writer_for(filename, format=format, multiframe=multiframe)
            #MDAnalysis.coordinates.writer(filename, **kwargs)
            coords = True
        except (ValueError, TypeError):
            coords = False

        try:
            # here `file_format` is only used as default,
            # anything pulled off `filename` will be used preferentially
            writer = get_selection_writer_for(filename, file_format)
            selection = True
        except (TypeError, NotImplementedError):
            selection = False

        if not (coords or selection):
            raise ValueError("No writer found for format: {}".format(filename))
        else:
            with writer(filename, n_atoms=self.n_atoms, **kwargs) as w:
                w.write(self.atoms)


class ResidueGroup(GroupBase):
    """ResidueGroup base class.

    This class is used by a :class:`Universe` for generating its
    Topology-specific :class:`ResidueGroup` class. All the
    :class:`TopologyAttr` components are obtained from
    :class:`GroupBase`, so this class only includes ad-hoc methods
    specific to ResidueGroups.

    ResidueGroups can be compared and combined using group operators. See the
    list of these operators on :class:`GroupBase`.

    .. deprecated:: 0.16.2
       *Instant selectors* of Segments will be removed in the 1.0 release.
       See :ref:`Instant selectors <instance-selectors>` for details and alternatives.

    """

    @property
    def atoms(self):
        """Get an :class:`AtomGroup` of atoms represented in this
        :class:`ResidueGroup`.

        The atoms are ordered locally by residue in the
        :class:`ResidueGroup`.  No duplicates are removed.

        """
        return self.universe.atoms[np.concatenate(self.indices)]

    @property
    def n_atoms(self):
        """Number of atoms represented in :class:`ResidueGroup`, including
        duplicate residues.

        Equivalent to ``len(self.atoms)``.

        """
        return len(self.atoms)

    @property
    def residues(self):
        """Get another :class:`ResidueGroup` identical to this one.

        """
        return self.universe.residues[self.ix]

    @property
    def n_residues(self):
        """Number of residues in ResidueGroup. Equivalent to ``len(self)``.

        """
        return len(self)

    @property
    def segments(self):
        """Get sorted SegmentGroup of the (unique) segments represented in the
        ResidueGroup.

        """
        return self.universe.segments[np.unique(self.segindices)]

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
                raise TypeError("Can only set ResidueGroup residues to Segment "
                                "or ResidueGroup not {}".format(
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
        """Number of unique segments represented in the ResidueGroup.

        Equivalent to ``len(self.segments)``.

        """
        return len(self.segments)

    @property
    def unique(self):
        """Return a ResidueGroup containing sorted and unique residues only.

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

        .. versionadded:: 0.16.0
        """
        return self.universe.residues[np.unique(self.ix)]


class SegmentGroup(GroupBase):
    """SegmentGroup base class.

    This class is used by a Universe for generating its Topology-specific
    SegmentGroup class. All the TopologyAttr components are obtained from
    GroupBase, so this class only includes ad-hoc methods specific to
    SegmentGroups.

    SegmentGroups can be compared and combined using group operators. See the
    list of these operators on :class:`GroupBase`.

    .. deprecated:: 0.16.2
       *Instant selectors* of Segments will be removed in the 1.0 release.
       See :ref:`Instant selectors <instance-selectors>` for details and alternatives.

    """

    @property
    def atoms(self):
        """Get an AtomGroup of atoms represented in this SegmentGroup.

        The atoms are ordered locally by residue, which are further ordered by
        segment in the SegmentGroup. No duplicates are removed.

        """
        return self.universe.atoms[np.concatenate(self.indices)]

    @property
    def n_atoms(self):
        """Number of atoms represented in SegmentGroup, including duplicate
        segments.

        Equivalent to ``len(self.atoms)``.

        """
        return len(self.atoms)

    @property
    def residues(self):
        """Get a ResidueGroup of residues represented in this SegmentGroup.

        The residues are ordered locally by segment in the SegmentGroup.
        No duplicates are removed.

        """
        return self.universe.residues[np.concatenate(self.resindices)]

    @property
    def n_residues(self):
        """Number of residues represented in SegmentGroup, including duplicate
        segments.

        Equivalent to ``len(self.residues)``.

        """
        return len(self.residues)

    @property
    def segments(self):
        """Get another SegmentGroup identical to this one.

        """
        return self.universe.segments[self.ix]

    @property
    def n_segments(self):
        """Number of segments in SegmentGroup. Equivalent to ``len(self)``.

        """
        return len(self)

    @property
    def unique(self):
        """Return a SegmentGroup containing sorted and unique segments only.

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

        .. versionadded:: 0.16.0

        """
        return self.universe.segments[np.unique(self.ix)]


@functools.total_ordering
class ComponentBase(_MutableBase):
    """Base class from which a Universe's Component class is built.

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

        If this component is an Atom, this is the index of the atom.
        If it is a Residue, this is the index of the residue.
        If it is a Segment, this is the index of the segment.

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
    """Atom base class.

    This class is used by a Universe for generating its Topology-specific Atom
    class. All the TopologyAttr components are obtained from ComponentBase, so
    this class only includes ad-hoc methods specific to Atoms.

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
            raise TypeError(
                "Can only set Atom residue to Residue, not {}".format(type(new)))
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

        .. note:: changing the position is not reflected in any files; reading any
                  frame from the trajectory will replace the change with that
                  from the file
        """
        return self.universe.trajectory.ts.positions[self.ix].copy()

    @position.setter
    def position(self, values):
        self.universe.trajectory.ts.positions[self.ix, :] = values

    @property
    def velocity(self):
        """Velocity of the atom.

        The velocity can be changed by assigning an array of shape (3,).

        .. note:: changing the velocity is not reflected in any files; reading any
                  frame from the trajectory will replace the change with that
                  from the file

        A :exc:`~MDAnalysis.NoDataError` is raised if the trajectory
        does not contain velocities.

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

        The force can be changed by assigning an array of shape (3,).

        .. note:: changing the force is not reflected in any files; reading any
                  frame from the trajectory will replace the change with that
                  from the file

        A :exc:`~MDAnalysis.NoDataError` is raised if the trajectory
        does not contain forces.

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
    """Residue base class.

    This class is used by a Universe for generating its Topology-specific
    Residue class. All the TopologyAttr components are obtained from
    ComponentBase, so this class only includes ad-hoc methods specific to
    Residues.

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
        return self.universe.atoms[self.universe._topology.indices[self][0]]

    @property
    def segment(self):
        return self.universe.segments[self.universe._topology.segindices[self]]

    @segment.setter
    def segment(self, new):
        if not isinstance(new, Segment):
            raise TypeError(
                "Can only set Residue segment to Segment, not {}".format(type(new)))
        self.universe._topology.tt.move_residue(self.ix, new.segindex)


class Segment(ComponentBase):
    """Segment base class.

    This class is used by a Universe for generating its Topology-specific
    Segment class. All the TopologyAttr components are obtained from
    ComponentBase, so this class only includes ad-hoc methods specific to
    Segments.

    .. deprecated:: 0.16.2
       *Instant selectors* of Segments will be removed in the 1.0 release.
       See :ref:`Instant selectors <instance-selectors>` for details and alternatives.

    """
    def __repr__(self):
        me = '<Segment'
        if hasattr(self, 'segid'):
            me += ' {}'.format(self.segid)
        return me + '>'

    @property
    def atoms(self):
        return self.universe.atoms[self.universe._topology.indices[self][0]]

    @property
    def residues(self):
        return self.universe.residues[self.universe._topology.resindices[self][0]]

    def __getattr__(self, attr):
        # DEPRECATED in 0.16.2
        # REMOVE in 1.0
        #
        # Segment.r1 access
        if attr.startswith('r') and attr[1:].isdigit():
            resnum = int(attr[1:])
            rg = self.residues[resnum - 1]  # convert to 0 based
            warnings.warn("Instant selectors Segment.r<N> will be removed in 1.0. "
                          "Use Segment.residues[N-1] instead.",
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
        self.selection_strings = strings
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
        # Run back through AtomGroup init with this information to remake ourselves
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
            `True` if the group's selection is up-to-date, `False` otherwise.

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
        Checks whether the selection needs updating and updates it, if needed.

        Returns
        -------
        bool
            `True` if the group was already up-to-date, `False` otherwise.

        """
        status = self.is_uptodate
        if not status:
            self.update_selection()
        return status

    def __getattribute__(self, name):
        # ALL attribute access goes through here
        # If the requested attribute isn't in the shortcut list, update ourselves
        if not name in _UAG_SHORTCUT_ATTRS:
            self._ensure_updated()
        # Going via object.__getattribute__ then bypasses this check stage
        return object.__getattribute__(self, name)

    def __reduce__(self):
        # strategy for unpickling is:
        # - unpickle base group
        # - recreate UAG as created through select_atoms (basegroup and selstrs)
        # even if base_group is a UAG this will work through recursion
        return (_unpickle_uag,
                (self._base_group.__reduce__(), self._selections, self.selection_strings))

    def __repr__(self):
        basestr = super(UpdatingAtomGroup, self).__repr__()
        if not self.selection_strings:
            return basestr
        sels = "'{}'".format("' + '".join(self.selection_strings))
        # Cheap comparison. Might fail for corner cases but this is
        # mostly cosmetic.
        if self._base_group is self.universe.atoms:
            basegrp = "the entire Universe."
        else:
            basegrp = "another AtomGroup."
        # With a shorthand to conditionally append the 's' in 'selections'.
        return "{}, with selection{} {} on {}>".format(basestr[:-1],
                    "s"[len(self.selection_strings)==1:], sels, basegrp)

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
    """Decorator to check if all AtomGroup arguments have certain attributes

    Example
    -------
    When used to wrap a function, will check all AtomGroup arguments for the
    listed requirements

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
