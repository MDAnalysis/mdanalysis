# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.mdanalysis.org
# Copyright (c) 2006-2016 The MDAnalysis Development Team and contributors
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

"""
from six.moves import zip

import numpy as np
import functools
import itertools
import os
import warnings

import MDAnalysis
from ..lib import util
from ..lib import distances
from ..lib import transformations
from . import selection
from . import flags
from . import levels
from ..exceptions import NoDataError
from . import topologyobjects


def make_classes():
    """Make a fresh copy of all Classes

    Returns
    -------
    Two dictionaries. One with a set of :class:`_TopologyAttrContainer` classes
    to serve as bases for universe-specific MDA container classes. Another,
    with the final merged versions of those classes. The classes themselves are
    used as hashing keys.
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

        Classes mixed at this point always override :meth:`__new__`, causing
        further instantiations to shortcut to :meth:`object.__new__` (skipping
        the cache-fetch process for :class:`_MutableBase` subclasses).

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
        return type(other.__name__, (_ImmutableBase, other, cls), {})

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
    Base class that merges appropriate :class:`TopologyAttrContainer` classes.

    Implements :meth:`__new__`. In it the instantiating class is fetched from
    :attr:`Universe._classes`. If there is a cache miss, a merged class is made
    with a base from :attr:`Universe._class_bases` and cached.

    The classes themselves are used as the cache dictionary keys for simplcity
    in cache retrieval.
    """
    # This signature must be kept in sync with the __init__ signature of
    #  GroupBase and ComponentBase.
    def __new__(cls, *args):
        try:
            ix, u = args
        except ValueError:
            # deprecated AtomGroup init method..
            u = args[0][0].universe
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
                                "classes are: {}".format(
                                    cls.__name__, str(sorted(u._class_bases.keys()))))
            newcls = u._classes[cls] = parent_cls._mix(cls)
            return object.__new__(newcls)


class _ImmutableBase(object):
    """Class used to shortcut :meth:`__new__` to :meth:`object.__new__`.

    """
    # When mixed via _TopologyAttrContainer._mix this class has MRO priority.
    #  Setting __new__ like this will avoid having to go through the
    #  cache lookup if the class is reused (as in ag.__class__(...)).
    __new__ = object.__new__


class GroupBase(_MutableBase):
    """Base class from which a Universe's Group class is built.

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
        self._ix = np.asarray(ix, dtype=np.int64)
        self._u = u
        self._cache = dict()

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
        if isinstance(item, (int, np.int_)):
            return self.level.singular(self._ix[item], self._u)
        else:
            if isinstance(item, list) and item:  # check for empty list
                # hack to make lists into numpy arrays
                # important for boolean slicing
                item = np.array(item)
            return self.__class__(self._ix[item], self._u)

    def __repr__(self):
        name = self.level.name
        return ("<{}Group with {} {}s>"
                "".format(name.capitalize(), len(self), name))

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
        if self.level != other.level:
            raise TypeError("Can't add different level objects")
        if self._u is not other._u:
            raise ValueError("Can't add objects from different Universes")

        # for the case where other is a Component, and so other._ix is an
        # integer
        if isinstance(other._ix, int):
            o_ix = np.array([other._ix])
        else:
            o_ix = other._ix

        return self.__class__(np.concatenate([self._ix, o_ix]), self._u)

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
            return self.__class__(self._ix, self._u)
        else:
            raise TypeError("unsupported operand type(s) for +:"
                            " '{}' and '{}'".format(type(self).__name__,
                                                    type(other).__name__))

    def __contains__(self, other):
        if not other.level == self.level:
            # maybe raise TypeError instead?
            # eq method raises Error for wrong comparisons
            return False
        return other.ix in self._ix

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
    def dimensions(self):
        return self._u.trajectory.ts.dimensions

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

        .. note::
            The :class:`MDAnalysis.core.flags` flag *use_pbc* when set to
            ``True`` allows the *pbc* flag to be used by default.

        Returns
        -------
         corners : array
            2x3 array giving corners of bounding box as
            [[xmin, ymin, zmin], [xmax, ymax, zmax]].


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

        .. note::
            The :class:`MDAnalysis.core.flags` flag *use_pbc* when set to
            ``True`` allows the *pbc* flag to be used by default.

        Returns
        -------
        R : float
            Radius of bounding sphere.
        center : array
            Coordinates of sphere center as ``[xcen,ycen,zcen]``.


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

    def rotate(self, R, point=None):
        r"""Apply a rotation matrix `R` to the selection's coordinates.

        Parameters
        ----------
        R : array_like
            3x3 rotation matrix to use for applying rotation.
        point : array_like, optional
            Center of rotation. If ``None`` then the center of geometry of this
            group is used.

        Returns
        -------
        self : AtomGroup

        Notes
        -----
        By default (``point=None``) the rotation is performed around
        the centroid of the group (:meth:`center_of_geometry`). In
        order to perform a rotation around, say, the origin, use
        ``point=[0, 0, 0]``.

        :math:`\mathsf{R}` is a 3x3 orthogonal matrix that transforms a vector
        :math:`\mathbf{x} \rightarrow \mathbf{x}'`:

        .. math::

            \mathbf{x}' = \mathsf{R}\mathbf{x}

        See Also
        --------
        rotateby : rotate around given axis and angle
        MDAnalysis.lib.transformations : module of all coordinate transforms

        """
        R = np.asarray(R)
        point = np.asarray(point) if point is not None else self.centroid()

        self.translate(-point)
        # changes the coordinates (in place)
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

        for o, s in itertools.izip(objects, shifts):
            # Save some needless shifts
            if not all(s == 0.0):
                o.atoms.translate(s)


class AtomGroup(GroupBase):
    """A group of atoms.

    An AtomGroup is an ordered collection of atoms. Typically, an AtomGroup is
    generated from a selection, or by indexing/slcing the AtomGroup of all
    atoms in the Universe at :attr:`MDAnalysis.core.universe.Universe.atoms`.

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

    Atoms can also be accessed in a Pythonic fashion by using the atom name as
    an attribute. For instance, ::

        ag.CA

    will provide a :class:`AtomGroup` of all CA atoms in the
    group. These *instant selector* attributes are auto-generated for
    each atom name encountered in the group.

    .. note::

        The name-attribute instant selector access to atoms is mainly
        meant for quick interactive work. Thus it either returns a
        single :class:`Atom` if there is only one matching atom, *or* a
        new :class:`AtomGroup` for multiple matches.  This makes it
        difficult to use the feature consistently in scripts.

    AtomGroup instances are always bound to a
    :class:`MDAnalysis.core.universe.Universe`. They cannot exist in isolation.

    .. SeeAlso:: :class:`MDAnalysis.core.universe.Universe`

    """
    def __getattr__(self, attr):
        # is this a known attribute failure?
        if attr in ('fragments',):  # TODO: Generalise this to cover many attributes
            # eg:
            # if attr in _ATTR_ERRORS:
            # raise NDE(_ATTR_ERRORS[attr])
            raise NoDataError("AtomGroup has no fragments this requires Bonds")
        elif hasattr(self.universe._topology, 'names'):
            # Ugly hack to make multiple __getattr__s work
            try:
                return self._get_named_atom(attr)
            except selection.SelectionError:
                pass
        raise AttributeError("{cls} has no attribute {attr}".format(
            cls=self.__class__.__name__, attr=attr))

    @property
    def atoms(self):
        """Get another AtomGroup identical to this one."""
        return self._u.atoms[self.ix]

    @property
    def n_atoms(self):
        """Number of atoms in AtomGroup.

        Equivalent to ``len(self)``."""
        return len(self)

    @property
    def residues(self):
        """Get sorted :class:`ResidueGroup` of the (unique) residues
        represented in the AtomGroup."""
        return self._u.residues[np.unique(self.resindices)]

    @residues.setter
    def residues(self, new):
        # Can set with Res, ResGroup or list/tuple of Res
        if isinstance(new, Residue):
            rix = itertools.cycle((new.resindex,))
        elif isinstance(new, ResidueGroup):
            rix = new.resindices
        else:
            try:
                rix = [r.resindex for r in new]
            except AttributeError:
                raise TypeError("Can only set AtomGroup residues to Residue "
                                "or ResidueGroup not {}".format(
                                    ', '.join(type(r) for r in new
                                              if not isinstance(r, Residue))
                                ))
        if not isinstance(rix, itertools.cycle) and len(rix) != len(self):
            raise ValueError("Incorrect size: {} for AtomGroup of size: {}"
                             "".format(len(new), len(self)))
        # Optimisation TODO:
        # This currently rebuilds the tt len(self) times
        # Ideally all changes would happen and *afterwards* tables are built
        # Alternatively, if the changes didn't rebuild table, this list
        # comprehension isn't terrible.
        for at, r in zip(self, rix):
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
        return self._u.segments[np.unique(self.segindices)]

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
        """Return an AtomGroup containing sorted and unique atoms only."""
        return self._u.atoms[np.unique(self.ix)]

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
        return self._u.trajectory.ts.positions[self._ix]

    @positions.setter
    def positions(self, values):
        ts = self._u.trajectory.ts
        ts.positions[self._ix, :] = values

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
        ts = self._u.trajectory.ts
        try:
            return np.array(ts.velocities[self._ix])
        except (AttributeError, NoDataError):
            raise NoDataError("Timestep does not contain velocities")

    @velocities.setter
    def velocities(self, values):
        ts = self._u.trajectory.ts
        try:
            ts.velocities[self._ix, :] = values
        except AttributeError:
            raise NoDataError("Timestep does not contain velocities")

    @property
    def forces(self):
        """Forces on each atom in the AtomGroup.

        The velocities can be changed by assigning an array of the appropriate
        shape, i.e. either Nx3 to assign individual velocities or 3 to assign
        the *same* velocity to all atoms (e.g. ``ag.velocity = array([0,0,0])``
        will give all particles zero velocity).

        """
        ts = self._u.trajectory.ts
        try:
            return ts.forces[self._ix]
        except (AttributeError, NoDataError):
            raise NoDataError("Timestep does not contain forces")

    @forces.setter
    def forces(self, values):
        ts = self._u.trajectory.ts
        try:
            ts.forces[self._ix, :] = values
        except AttributeError:
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

    def select_atoms(self, sel, *othersel, **selgroups):
        """Select atoms using a selection string.

        Returns an :class:`AtomGroup` with atoms sorted according to
        their index in the topology (this is to ensure that there
        are not any duplicates, which can happen with complicated
        selections).

        Existing :class:`AtomGroup` objects can be passed as named arguments,
        which will then be available to the selection parser.

        Subselections can be grouped with parentheses.

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

            group *group-name*
                selects the atoms in the :class:`AtomGroup` passed to the
                function as an argument named *group-name*. Only the atoms
                common to *group-name* and the instance :meth:`~select_atoms`
                was called from will be considered. *group-name* will be
                 included in the parsing just by comparison of atom indices.
                This means that it is up to the user to make sure they were
                defined in an appropriate :class:`Universe`.

            fullgroup *group-name*
                just like the ``group`` keyword with the difference that all the
                atoms of *group-name* are included. The resulting selection may
                therefore have atoms that were initially absent from the
                instance :meth:`~select_atoms` was called from.


        .. versionchanged:: 0.7.4
           Added *resnum* selection.
        .. versionchanged:: 0.8.1
           Added *group* and *fullgroup* selections.
        .. versionchanged:: 0.13.0
           Added *bonded* selection
        .. versionchanged:: 0.16.0
           Resid selection now takes icodes into account where present.

        """
        atomgrp = selection.Parser.parse(sel, selgroups).apply(self)
        if othersel:
            # Generate a selection for each selection string
            for sel in othersel:
                atomgrp += selection.Parser.parse(sel, selgroups).apply(self)
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
            return [self._u.atoms[[a.ix]] for a in self]

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
        return topologyobjects.Bond(self._ix, self.universe)

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
        return topologyobjects.Angle(self._ix, self.universe)

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
        return topologyobjects.Dihedral(self._ix, self.universe)

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
        return topologyobjects.ImproperDihedral(self._ix, self.universe)

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

        format : str, optional
            PDB, CRD, GRO, VMD (tcl), PyMol (pml), Gromacs (ndx) CHARMM (str)
            Jmol (spt); case-insensitive and can also be supplied as the
            filename extension [PDB]

        filenamefmt : str, optional
            format string for default filename; use substitution tokens
            'trjname' and 'frame' ["%(trjname)s_%(frame)d"]

        bonds : str, optional
           how to handle bond information, especially relevant for PDBs;
           default is ``"conect"``.

           * ``"conect"``: write only the CONECT records defined in the original
             file
           * ``"all"``: write out all bonds, both the original defined and those
             guessed by MDAnalysis
           * ``None``: do not write out bonds


        .. versionchanged:: 0.9.0
           Merged with write_selection.  This method can now write both
           selections out.

        """
        import MDAnalysis.coordinates
        import MDAnalysis.selections

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
        try:
            writer = MDAnalysis.coordinates.writer(filename, **kwargs)
            coords = True
        except TypeError as e:
            # might be selections format
            if 'got an unexpected keyword argument' in str(e):
                raise
            coords = False

        try:
            SelectionWriter = MDAnalysis.selections.get_writer(filename,
                                                               file_format)
            writer = SelectionWriter(filename, **kwargs)
            selection = True
        except (TypeError, NotImplementedError):
            selection = False

        if not (coords or selection):
            raise ValueError("No writer found for format: {}".format(filename))
        else:
            writer.write(self.atoms)
            if coords:  # only these writers have a close method
                writer.close()


class ResidueGroup(GroupBase):
    """ResidueGroup base class.

    This class is used by a :class:`Universe` for generating its
    Topology-specific :class:`ResidueGroup` class. All the
    :class:`TopologyAttr` components are obtained from
    :class:`GroupBase`, so this class only includes ad-hoc methods
    specific to ResidueGroups.

    """
    @property
    def atoms(self):
        """Get an :class:`AtomGroup` of atoms represented in this
        :class:`ResidueGroup`.

        The atoms are ordered locally by residue in the
        :class:`ResidueGroup`.  No duplicates are removed.

        """
        return self._u.atoms[np.concatenate(self.indices)]

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
        return self._u.residues[self.ix]

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
        return self._u.segments[np.unique(self.segindices)]

    @segments.setter
    def segments(self, new):
        # Can set with Seg, SegGroup or list/tuple of Seg
        if isinstance(new, Segment):
            six = itertools.cycle((new.segindex,))
        elif isinstance(new, SegmentGroup):
            six = new.segindices
        else:
            try:
                six = [s.segindex for s in new]
            except AttributeError:
                raise TypeError("Can only set ResidueGroup residues to Segment "
                                "or ResidueGroup not {}".format(
                                    ', '.join(type(r) for r in new
                                              if not isinstance(r, Segment))
                                ))
        if not isinstance(six, itertools.cycle) and len(six) != len(self):
            raise ValueError("Incorrect size: {} for ResidueGroup of size: {}"
                             "".format(len(new), len(self)))
        # Optimisation TODO:
        # This currently rebuilds the tt len(self) times
        # Ideally all changes would happen and *afterwards* tables are built
        # Alternatively, if the changes didn't rebuild table, this list
        # comprehension isn't terrible.
        for r, s in zip(self, six):
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

        """
        return self._u.residues[np.unique(self.ix)]


class SegmentGroup(GroupBase):
    """SegmentGroup base class.

    This class is used by a Universe for generating its Topology-specific
    SegmentGroup class. All the TopologyAttr components are obtained from
    GroupBase, so this class only includes ad-hoc methods specific to
    SegmentGroups.

    """
    @property
    def atoms(self):
        """Get an AtomGroup of atoms represented in this SegmentGroup.

        The atoms are ordered locally by residue, which are further ordered by
        segment in the SegmentGroup. No duplicates are removed.

        """
        return self._u.atoms[np.concatenate(self.indices)]

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
        return self._u.residues[np.concatenate(self.resindices)]

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
        return self._u.segments[self.ix]

    @property
    def n_segments(self):
        """Number of segments in SegmentGroup. Equivalent to ``len(self)``.

        """
        return len(self)

    @property
    def unique(self):
        """Return a SegmentGroup containing sorted and unique segments only.

        """
        return self._u.segments[np.unique(self.ix)]


@functools.total_ordering
class ComponentBase(_MutableBase):
    """Base class from which a Universe's Component class is built.

    Components are the individual objects that are found in Groups.
    """
    def __init__(self, ix, u):
        # index of component
        self._ix = ix
        self._u = u

    def __repr__(self):
        return ("<{} {}>"
                "".format(self.level.name.capitalize(), self._ix))

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
        if self.level != other.level:
            raise TypeError('Can only add {0}s or {1}s (not {2}s/{3}s)'
                            ' to {0}'.format(self.level.singular.__name__,
                                             self.level.plural.__name__,
                                             other.level.singular.__name__,
                                             other.level.plural.__name__))

        if self.universe is not other.universe:
            raise ValueError("Can only add objects from the same Universe")

        if isinstance(other.ix, int):
            o_ix = np.array([other.ix])
        else:
            o_ix = other.ix

        return self.level.plural(
                np.concatenate((np.array([self.ix]), o_ix)), self.universe)

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
            return self.level.plural(np.array([self._ix]), self._u)
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

    @property
    def residue(self):
        residueclass = self.level.parent.singular
        return residueclass(self._u._topology.resindices[self],
                            self._u)

    @residue.setter
    def residue(self, new):
        if not isinstance(new, Residue):
            raise TypeError(
                "Can only set Atom residue to Residue, not {}".format(type(new)))
        self.universe._topology.tt.move_atom(self.ix, new.resindex)

    @property
    def segment(self):
        segmentclass = self.level.parent.parent.singular
        return segmentclass(self._u._topology.segindices[self],
                            self._u)

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
        return self._u.trajectory.ts.positions[self._ix].copy()

    @position.setter
    def position(self, values):
        self._u.trajectory.ts.positions[self._ix, :] = values

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
        ts = self._u.trajectory.ts
        try:
            return ts.velocities[self._ix].copy()
        except (AttributeError, NoDataError):
            raise NoDataError("Timestep does not contain velocities")

    @velocity.setter
    def velocity(self, values):
        ts = self._u.trajectory.ts
        try:
            ts.velocities[self.index, :] = values
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
        ts = self._u.trajectory.ts
        try:
            return ts.forces[self._ix].copy()
        except (AttributeError, NoDataError):
            raise NoDataError("Timestep does not contain forces")

    @force.setter
    def force(self, values):
        ts = self._u.trajectory.ts
        try:
            ts.forces[self._ix, :] = values
        except (AttributeError, NoDataError):
            raise NoDataError("Timestep does not contain forces")


class Residue(ComponentBase):
    """Residue base class.

    This class is used by a Universe for generating its Topology-specific
    Residue class. All the TopologyAttr components are obtained from
    ComponentBase, so this class only includes ad-hoc methods specific to
    Residues.

    """
    @property
    def atoms(self):
        atomsclass = self.level.child.plural
        return atomsclass(self._u._topology.indices[self][0],
                          self._u)

    @property
    def segment(self):
        segmentclass = self.level.parent.singular
        return segmentclass(self._u._topology.segindices[self],
                            self._u)

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

    """
    @property
    def atoms(self):
        atomsclass = self.level.child.child.plural
        return atomsclass(self._u._topology.indices[self][0],
                          self._u)

    @property
    def residues(self):
        residuesclass = self.level.child.plural
        return residuesclass(self._u._topology.resindices[self][0],
                             self._u)

    def __getattr__(self, attr):
        # Segment.r1 access
        if attr.startswith('r') and attr[1:].isdigit():
            resnum = int(attr[1:])
            return self.residues[resnum - 1]  # convert to 0 based
        # Resname accesss
        if hasattr(self.residues, 'resnames'):
            try:
                return self.residues._get_named_residue(attr)
            except selection.SelectionError:
                pass
        raise AttributeError("{cls} has no attribute {attr}"
                             "".format(cls=self.__class__.__name__, attr=attr))

# Define relationships between these classes
# with Level objects
ATOMLEVEL = levels.Level('atom', Atom, AtomGroup)
RESIDUELEVEL = levels.Level('residue', Residue, ResidueGroup)
SEGMENTLEVEL = levels.Level('segment', Segment, SegmentGroup)

ATOMLEVEL.parent = RESIDUELEVEL
ATOMLEVEL.child = None
RESIDUELEVEL.parent = SEGMENTLEVEL
RESIDUELEVEL.child = ATOMLEVEL
SEGMENTLEVEL.parent = None
SEGMENTLEVEL.child = RESIDUELEVEL

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
