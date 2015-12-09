# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#


"""
Fundamental building blocks --- :mod:`MDAnalysis.core.AtomGroup`
================================================================

The most important data structure in MDAnalysis is the
:class:`AtomGroup`, which contains :class:`Atom` instances.

A :class:`Universe` is the user-visible entry point and collects all
information needed to analyze a structure or a whole trajectory.

Segments and residues are a way to refer to a collection of atoms. By
convention, a :class:`Residue` is a single amino acid, or a water
molecule, ion, or ligand. A :class:`Segment` is a collection of
residues such as a whole protein or a chain in a protein or all the
water in the system.

Class Hierarchy
---------------

A :class:`Universe` contains Segments, which contain Residues, which
contain Atoms; all containers are derived from :class:`AtomGroup`, and
thus one can always analyze them as a collection of atoms, independent
of the hierarchical level.

Each :class:`Atom` can only belong to a single :class:`Residue`, and a
:class:`Residue` belongs to one specific :class:`Segment`. This
hierarchy can be described as ::

    Segment > Residue > Atom

Depending on the use case, it can be more convenient to access data
on, for instance, the basis of residues than atoms, or to write out
individual chains (segments) of a protein. MDAnalysis simply provides
access to these groupings and keeps track of where an atom
belongs. Each object provides three attributes (:attr:`~AtomGroup.atoms`,
:attr:`~AtomGroup.residues` or :attr:`~Atom.residue`, :attr:`~AtomGroup.segments` or
:attr:`~Atom.segment`) that give access to the tiers in the hierarchy
that the object belongs to.


Manipulating atoms, residues, and segments
------------------------------------------

When working with MDAnalysis it is useful to remember that the fundamental
object is the :class:`Atom`. Each particle in the topology is represented by
exactly one :class:`Atom` instance. One :class:`Atom`, however, can be a member
of multiple :class:`AtomGroup` collections, for instance from different
selections even though they all refer to the *same* :class:`Atom` object. Thus,
changing a property of a specific and :class:`Atom` in one :class:`AtomGroup`
changes it "everywhere".

The same is mostly true for :class:`Residue` instances although they are
derived from :class:`Atom` instances: all :class:`Atom` objects with the same
:attr:`Atom.resid` are bundled into a single :class:`Residue` with
:class:`Residue.id` = *resid*. This means that just changing, say, the residue
name with a command such as ::

  >>> r = u.select_atoms("resid 99").residues[0]
  >>> print(r)
  <Residue 'ALA', 99>
  >>> r.name = "UNK"
  >>> print(r)
  <Residue 'UNK', 99>
  >>> rnew = u.select_atoms("resid 99").residues[0]
  >>> print(rnew)
  <Residue 'UNK', 99>

will typically work as expected. When working with collections such as
:class:`AtomGroup` or :class:`ResidueGroup` it is generally better to use
properties such as :attr:`AtomGroup.resnames` or :attr:`ResidueGroup.resnames`
to modify their items' attributes.

There are two cases when it is very important to use these collective
properties:

* changing *resid*: :attr:`AtomGroup.resids` and :attr:`ResidueGroup.resids`
* changing *segid*: :attr:`AtomGroup.segids` and :attr:`ResidueGroup.segids`

Because residues are determined by the :attr:`Atom.resid` and segments by
:attr:`Atom.segid`, the above properties take extra care to rebuild the list of
segments and residues. Alternatively, the same effect can be obtained using
the corresponding setter method, e.g. :meth:`AtomGroup.set_resids`.

.. Note::

   Setting any of
   :attr:`AtomGroup.resids`, :attr:`ResidueGroup.resids`,
   :attr:`AtomGroup.segids`, :attr:`ResidueGroup.segids`
   can change the topology: they can split or merge residues or segments.

Splitting/merging of residues is probably not very useful because no chemical
rearrangements are carried out. Manipulating segments might be more useful in
order to add additional structure to a :class:`Universe` and provide instant
segment selectors for interactive work::

  u.select_atoms("protein").set_segids("protein")
  u.select_atoms("resname POPE or resname POPC").set_segids("lipids")
  u.select_atoms("resname SOL").set_segids("water")
  u.select_atoms("resname NA or resname CL").set_segids("ions")

  u.protein.n_residues
  water_oxygens = u.water.OW

The setter methods have the additional advantage that they can assign
lists. For instance, many MD codes number residues consecutively starting from
1. However, the original structure might be missing a few atoms at the
N-terminus. Let's say that the first residue is really residue 10. In order to
store the canonical residue IDs ("resnum") one could the use ::

  protein = u.select_atoms("protein").residues
  protein.set_resnums(protein.resnums + 9)

One can then use ``protein.select("resnum 42")`` to select the residue that has
the canonical residue id 42 (instead of ``resid 33``).

One can also read the resids directly from an original PDB file::

  orig = MDAnalysis.Universe("2jln.pdb")
  protein.set_resnums(orig.select_atoms("protein").resids)


Working with Topologies
-----------------------

If the topology file given to the Universe had bonding information then this
will have been loaded into the Universe as :attr:`Universe.bonds`
:attr:`Universe.angles` :attr:`Universe.dihedrals` and :attr:`Universe.impropers`.


If your topology file does not have this information, it is still possible
to construct it based on the positions of the atoms and assumed vdw radii
for these atoms.  See :meth:`MDAnalysis.AtomGroup.guess_bonds` and
:func:`MDAnalysis.topology.core.guess_bonds` for details.

This Topology information is stored in :class:`MDAnalysis.core.topologyobjects.TopologyGroup`
objects.  These are designed to be analogous to the AtomGroup container for Atoms.

For examples working with a box of ethanol::

    >>> import MDAnalysis as mda
    >>> u = mda.Universe('ethanol.gro', guess_bonds=True)
    >>> u.bonds
    <TopologyGroup containing 11784 Bonds>
    >>> u.bonds.types()  # view available types
    [('O', 'H'), ('C', 'O'), ('C', 'H'), ('C', 'C')]
    >>> u.bonds.select_bonds(('C', 'O'))  # return all C-O bonds from the group
    <TopologyGroup containing 1473 Bonds>

Bonds are categorised based on the types of atoms.  This is done in a way
so that type (a, b, c) is equivalent to (c, b, a) ie. bonds are reversible.
For example::

    >>> u.angles.types()
    [('C', 'C', 'H'),
     ('H', 'C', 'H'),
     ('C', 'O', 'H'),
     ('C', 'C', 'O'),
     ('H', 'C', 'O')]

There are only C-C-H bonds and no H-C-C bonds.  Selection however is
aware that sometimes types are reversed::

    >>> u.angles.select_bonds(('H', 'C', 'C'))  # note reversal of type
    <TopologyGroup containing 7365 Angles>

TopologyGroups can be combined and indexed::

    >>> tg = u.angles.select_bonds(('C', 'C', 'O')) + u.angles.select_bonds(('C', 'O', 'H'))
    >>> tg.types()
    [('C', 'O', 'H'), ('C', 'C', 'O')]
    >>> tg[:100]
    <TopologyGroup containing 100 Angles>

Finally, TopologyGroups are linked to some fast Cython calculation methods to
determine bond lengths and angle sizes::

    >>> tg.values()
    array([ 1.88042373,  1.95928987,  1.74770012, ...,  1.79306789,
            1.95522678,  1.88881045])


Combining objects: system building
----------------------------------

It is often convenient to combined multiple groups of atoms into a single
object. If they are contained in a single :class:`Universe` then the methods
described above (especially manipulating the segments) might be
useful. However, if the atoms reside in different universes, the :func:`Merge`
function can be used.

Merging
~~~~~~~

In the following example for :func:`Merge`, protein, ligand, and solvent were
externally prepared in three different PDB files. They are loaded into separate
:class:`Universe` objects (where they could be further manipulated,
e.g. renumbered, relabeled, rotated, ...) The :func:`Merge` command is used to
combine all of them together::

    import MDAnalysis
    u1 = MDAnalysis.Universe("protein.pdb")
    u2 = MDAnalysis.Universe("ligand.pdb")
    u3 = MDAnalysis.Universe("solvent.pdb")
    u = MDAnalysis.Merge(u1.select_atoms("protein"), u2.atoms, u3.atoms)
    u.atoms.write("system.pdb")

The complete system is then written out to a new PDB file.

Replicating
~~~~~~~~~~~

It is also possible to replicate a molecule to build a system with
multiple copies of the same molecule. In the example, we replicate an
AdK molecule and then translate and rotate the second copy::

    import MDAnalysis; from MDAnalysis.tests.datafiles import *
    u = MDAnalysis.Universe(PSF, DCD)
    p = u.select_atoms("protein")
    m = MDAnalysis.Merge(p,p)

    # now renumber resids and segids for each copy

    # first copy of the protein (need to use atom indices because currently that's the only reliable property in the
    merged universe)
    p1 = m.select_atoms("bynum 1:3341")
    # second copy
    p2 = m.select_atoms("bynum 3342:6682")

    p1.set_segid("A")
    p2.set_segid("B")
    p2.residues.set_resid(p2.residues.resids + p1.residues.resids[-1])  # increment resids for p2 with the last
    resid from p1

    # you must regenerate the selections after modifying them (see notes in the docs!)
    # because the changed resids are not reflected in the selection (due to how residues are referenced internally)
    p1 = m.select_atoms("segid A")       # or as before:  m.select_atoms("bynum 1:3341")
    p2 = m.select_atoms("segid B")

    # rotate and translate
    p2.rotateby(180, [0,0,1])
    p2.translate([50,0,0])

Note that we have to manually set the residue numbers (resids) and
segment identifies because :func:`Merge` simply concatenates the
existing atoms and only ensures that all data structures are contained
in the new merged universe.


Classes and functions
---------------------

.. autoclass:: Universe
   :members:
.. autofunction:: Merge
.. autoclass:: AtomGroup
   :members:

   .. attribute:: _atoms

      immutable list of references to the atoms :class:`Atom` in the
      group

   .. automethod:: _rebuild_caches

   .. automethod:: _clear_caches

   .. automethod:: _fill_cache

.. autoclass:: Atom
   :members:

   .. attribute::     index

      atom index

   .. attribute::     segid

      name of the segment

   .. attribute::     resid

      residue number

   .. attribute::     resnum

      canonical residue number as, for instance, used in the original
      PDB file

      .. versionadded:: 0.7.4

   .. attribute::        resname

      residue name

   .. attribute::        residue

      :class:`Residue` object containing the atoms

   .. attribute::     id

      atom number inside the residue

   .. attribute::       name

      string, short name

   .. attribute::        type

      string or number (from force field), describing the atom type;
      stored as a string.

      .. versionchanged:: 0.7.6
         The :attr:`Atom.type` attribute is always stored as a string.

   .. attribute::        mass

      float, in `atomic mass units`_ (u)

   .. attribute::        charge

      float, in `electron charges`_ (*e*)

   .. attribute::        radius

      Born-radius for electrostatic calculations. (Only if read from a
      PQR file with the
      :class:`~MDAnalysis.coordinates.PQR.PQRReader`.)

   .. attribute::        altLoc

      Alternate location indicator (as used in `ATOM_` records in PDB
      files)

      .. versionadded:: 0.9.0

   .. attribute:: bonds

      A :class:`MDAnalysis.core.topologyobjects.TopologyGroup` of all
      :class:`~MDAnalysis.topogology.objects.Bond` instances that
      contains all the bonds that this atom participates in.

      .. versionadded:: 0.8.1

   .. attribute:: angles

      A :class:`MDAnalysis.core.topologyobjects.TopologyGroup` of all
      :class:`~MDAnalysis.topogology.objects.Angle` instances
      that contains all the angles that this atom participates in.

      .. versionadded:: 0.8.1

   .. attribute:: dihedrals

      A :class:`MDAnalysis.core.topologyobjects.TopologyGroup` of all
      :class:`~MDAnalysis.topogology.objects.Dihedral` instances
      that contains all the dihedrals that this atom
      participates in.

      .. versionadded:: 0.8.1

   .. attribute:: impropers

      A :class:`MDAnalysis.core.topologyobjects.TopologyGroup` of all
      :class:`MDAnalysis.core.topologyobjects.ImproperDihedral` instances
      that this atom is present in.

      .. versionadded:: 0.8.1

.. autoclass:: Residue
   :members:
.. autoclass:: ResidueGroup
   :members:
.. autoclass:: Segment
   :members:
.. autoclass:: SegmentGroup
   :members:

.. autofunction:: as_Universe
.. autoexception:: SelectionError
   :no-members:
.. autoexception:: SelectionWarning
   :no-members:
.. autoexception:: NoDataError
   :no-members:

.. _atomic mass units: http://physics.nist.gov/cgi-bin/cuu/Value?u
.. _electron charges: http://physics.nist.gov/cgi-bin/cuu/Value?e
.. _ATOM: http://www.wwpdb.org/documentation/format23/sect9.html#ATOM

"""
from __future__ import print_function, absolute_import

# Global imports
import warnings
import numpy as np
from numpy.linalg import eig
from numpy.lib.utils import deprecate
import itertools
from collections import defaultdict
import copy
import logging
import os.path
import weakref
import gc
import functools

# Local imports
import MDAnalysis
from .. import SelectionError, NoDataError, SelectionWarning
from ..lib import util
from ..lib import distances
from ..lib import mdamath
from ..lib import transformations
from ..lib.util import cached
from . import topologyobjects as top
from .universe import Universe

logger = logging.getLogger("MDAnalysis.core.AtomGroup")


# Constant to translate the name of an Atom's property
# to the plural version, as found in AtomGroup
_PLURAL_PROPERTIES = {'index': 'indices',
                      'name': 'names',
                      'type': 'types',
                      'resname': 'resnames',
                      'resid': 'resids',
                      'segid': 'segids',
                      'mass': 'masses',
                      'charge': 'charges',
                      'radius': 'radii',
                      'bfactor': 'bfactors',
                      'resnum': 'resnums',
                      'altLoc': 'altLocs',
                      'serial': 'serials',
                      'value': 'values',
                      'occupancy': 'occupancies'}
# And the return route
_SINGULAR_PROPERTIES = {v: k for k, v in _PLURAL_PROPERTIES.items()}


@functools.total_ordering
class Atom(object):
    """A class representing a single atom.

    :class:`Atom` is the basic building block of all larger data
    structures in MDAnalysis, in particular of the
    :class:`AtomGroup`.

    An :class:`Atom` is typically generated by a :ref:`topology reader
    <Supported topology formats>` from :mod:`MDAnalysis.topology`.

    For performance reasons, only a predefined number of attributes
    are included (and thus it is not possible to add attributes "on
    the fly"; they have to be included in the class definition).

    An :class:`Atom` is bound to a particular :class:`Universe`, but
    via a weak reference only. This means that the :class:`Atom`, and
    any :class:`AtomGroup` it is in, are only relevant while the
    originating :class:`Universe` is in scope.

    .. versionchanged 0.9.0
       Added fragment managed property.
       Changed bonds angles torsions impropers to be a managed property
    .. versionchanged 0.11.0
       Changed references to :class:`Universe` to be weak.
       Renamed atom.number to atom.index
       Renamed atom.torsions to atom.dihedrals
    .. versionchanged:: 0.11.1
       Added occupancy property. Can get and set.
    """

    __slots__ = (
        "index", "id", "name", "type", "resname", "resid", "segid",
        "mass", "charge", "residue", "segment",
        "_universe",
        "radius", "bfactor", "resnum", "serial", "altLoc")

    def __init__(self, index, name, type, resname, resid, segid, mass, charge,
                 residue=None, segment=None, radius=None, bfactor=None,
                 resnum=None, serial=None, altLoc=None, universe=None):
        self.index = index
        self.name = name
        self.altLoc = altLoc
        self.type = str(type)  # always a string (needed for selections)
        self.resname = resname
        self.resid = resid
        self.resnum = resnum
        self.residue = residue  # typically patched in later
        self.segid = segid
        self.segment = segment  # typically patched in later
        self.mass = mass
        self.charge = charge
        self.radius = radius
        self.bfactor = bfactor
        self.serial = serial
        self._universe = universe

    def __repr__(self):
        return ("<Atom {idx}: {name} of type {t} of resname {rname}, "
                "resid {rid} and segid {sid}{altloc}>".format(
                    idx=self.index + 1, name=self.name, t=self.type,
                    rname=self.resname, rid=self.resid, sid=self.segid,
                    altloc="" if not self.altLoc
                    else " and altloc {0}".format(self.altLoc)))

    def __lt__(self, other):
        return self.index < other.index

    def __eq__(self, other):
        return self.index == other.index

    def __hash__(self):
        return hash(self.index)

    def __add__(self, other):
        if not isinstance(other, (Atom, AtomGroup)):
            raise TypeError('Can only add Atoms or AtomGroups (not "{0}")'
                            ' to Atom'.format(other.__class__.__name__))
        if not self.universe is other.universe:
            raise ValueError("Can only add objects from the same Universe")
        if isinstance(other, Atom):
            return AtomGroup([self, other])
        else:
            return AtomGroup([self] + other._atoms)

    @property
    def number(self):
        """The index of this atom"""
        return self.index

    @property
    def pos(self):
        """coordinates of the atom

        Get the current cartesian coordinates of the atom (read-only).

        .. deprecated:: 0.8
           use :attr:`position`
        """
        return self.position

    @property
    def position(self):
        """coordinates of the atom

        Get the current cartesian coordinates of the atom.

        :Returns: a (3,) shape numpy array
        """
        return self.universe.coord.positions[self.index]  # internal numbering starts at 0

    @position.setter
    def position(self, coords):
        """
        Set the current cartesian coordinates of the atom.
        @param coords: a 1x3 numpy array of {x,y,z} coordinates, or optionally
            a single scalar if you should want to set all coordinates to the same value.
        """
        self.universe.coord.positions[self.index, :] = coords  # internal numbering starts at 0

    @property
    def velocity(self):
        """Current velocity of the atom.

        :Returns: a (3,) shape numpy array

        A :exc:`~MDAnalysis.NoDataError` is raised if the trajectory
        does not contain velocities.

        .. versionadded:: 0.7.5
        """
        # TODO: Remove error checking here (and all similar below)
        # and add to Timestep
        try:
            return self.universe.coord.velocities[self.index]
        except (AttributeError, NoDataError):
            raise NoDataError("Timestep does not contain velocities")

    @velocity.setter
    def velocity(self, vals):
        """Set the current velocity of the atom.

        A :exc:`~MDAnalysis.NoDataError` is raised if the trajectory
        does not contain velocities.

        .. versionadded:: 0.9.2
        """
        try:
            self.universe.coord.velocities[self.index] = vals
        except (AttributeError, NoDataError):
            raise NoDataError("Timestep does not contain velocities")

    @property
    def occupancy(self):
        """Access occupancy values.

        If available can have a value between 0 and 1

        Returns
        -------
        float
            occupancy of Atom

        .. versionadded:: 0.11.1
        """
        try:
            return self.universe.coord.data['occupancy'][self.index]
        except KeyError:
            raise NoDataError('Timestep does not contain occupancy')

    @occupancy.setter
    def occupancy(self, _occupancy):
        """Set occupancy for an atom

        If no occupancies are set in the universe of the atom the occupancy
        of the other atoms will be set to 1.
        """
        try:
            self.universe.coord.data['occupancy'][self.index] = _occupancy
        except KeyError:
            n_atoms = self.universe.coord.n_atoms
            self.universe.coord.data['occupancy'] = np.ones(n_atoms)
            self.universe.coord.data['occupancy'][self.index] = _occupancy

    @property
    def force(self):
        """Current force of the atom.

        :Returns: a (3,) shape numpy array

        A :exc:`~MDAnalysis.NoDataError` is raised if the trajectory
        does not contain velocities.

        .. versionadded:: 0.9.2
        """
        try:
            return self.universe.coord.forces[self.index]
        except (AttributeError, NoDataError):
            raise NoDataError("Timestep does not contain forces")

    @force.setter
    def force(self, vals):
        """Set the current force of the atom.

        .. versionadded:: 0.9.2
        """
        try:
            self.universe.coord.forces[self.index] = vals
        except (AttributeError, NoDataError):
            raise NoDataError("Timestep does not contain forces")

    def centroid(self):
        """The centroid of an atom is its position, :attr:`Atom.position`."""
        # centroid exists for compatibility with AtomGroup
        return self.position

    @property
    def universe(self):
        """A pointer back to the Universe of this Atom"""
        if self._universe is None:
            raise NoDataError("This Atom does not belong to a Universe")
        else:
            return self._universe

    @universe.setter
    def universe(self, new):
        self._universe = new

    @property
    def bonded_atoms(self):
        """An AtomGroup of the Atoms that this Atom is bonded to.

        .. versionadded:: 0.9
        """
        return AtomGroup([b.partner(self) for b in self.bonds])

    # The following look up a dictionary stored in the Universe.
    # These dictionaries are lazily built
    @property
    def fragment(self):
        """The fragment that this Atom is part of

        .. versionadded:: 0.9.0
        """
        return self.universe._fragmentDict[self]

    @property
    def bonds(self):
        """A TopologyGroup of the bonds that this Atom is in

        .. versionchanged:: 0.9.0
           Changed to managed property
        """
        return top.TopologyGroup(self.universe._bondDict[self])

    @property
    def angles(self):
        """A TopologyGroup of the angles that this Atom is in

        .. versionchanged:: 0.9.0
           Changed to managed property
        """
        return top.TopologyGroup(self.universe._angleDict[self])

    @property
    def dihedrals(self):
        """A TopologyGroup of the dihedrals that this Atom is in

        .. versionchanged:: 0.9.0
           Changed to managed property
        .. versionchanged:: 0.11.0
           Renamed to dihedrals (was torsions)
        """
        return top.TopologyGroup(self.universe._dihedralDict[self])

    @property
    def impropers(self):
        """A TopologyGroup of the improper dihedrals that this Atom is in

        .. versionchanged:: 0.9.0
           Changed to managed property
        """
        return top.TopologyGroup(self.universe._improperDict[self])


class AtomGroup(object):
    """A group of atoms.

      ag = universe.select_atoms(atom-list)

    The AtomGroup contains a list of atoms; typically, a AtomGroup is generated
    from a selection. It is build from any list-like collection of
    :class:`Atom` instances. It is also possible to create an empty AtomGroup
    from an empty list.

    An AtomGroup can be indexed and sliced like a list::

       ag[0], ag[-1]

    will return the first and the last :class:`Atom` in the group
    whereas the slice ::

       ag[0:6:2]

    returns every second element, corresponding to indices 0, 2, and 4.

    It also supports "advanced slicing" when the argument is a
    :class:`numpy.ndarray` or a :class:`list`::

       aslice = [0, 3, -1, 10, 3]
       ag[aslice]

    will return a new :class:`AtomGroup` containing (ag[0], ag[3], ag[-1],
    ag[10], ag[3]).

    .. Note::

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

    .. Note::

       The name-attribute instant selector access to atoms is mainly
       meant for quick interactive work. Thus it either returns a
       single :class:`Atom` if there is only one matching atom, *or* a
       new :class:`AtomGroup` for multiple matches.  This makes it
       difficult to use the feature consistently in scripts but it is
       much better for interactive work.

    AtomGroup instances are bound to a :class:`Universe`, but only through the
    weak reference :class:`Atom` has. The connection is lost as soon as the
    :class:`Universe` goes out of scope.

    AtomGroups can also be pickled and unpickled. When pickling, the Atom
    indices are serialized alongside the Universe number-of-atoms, and topology
    and trajectory filenames. When unpickling, all built Universes will be
    searched for one with matching number-of-atoms and filenames. Finer control
    over which :class:`Universe` gets chosen to base the unpickling on can
    be exerted using the *is_anchor* and *anchor_name* flags upon Universe creation
    or the methods/attributes :meth:`Universe.make_anchor`, :meth:`Universe.remove_anchor`,
    and :attr:`Universe.anchor_name`.


    .. rubric:: References for analysis methods

    .. [Dima2004] Dima, R. I., & Thirumalai, D. (2004). Asymmetry in the
                  shapes of folded and denatured states of proteins. *J
                  Phys Chem B*, 108(21),
                  6564-6570. doi:`10.1021/jp037128y`_

    .. _10.1021/jp037128y: http://dx.doi.org/10.1021/jp037128y


    .. SeeAlso:: :class:`Universe`
    .. versionchanged:: 0.7.6
       An empty AtomGroup can be created and no longer raises a
       :exc:`NoDataError`.
    .. versionchanged:: 0.9.0
       The size at which cache is used for atom lookup is now stored as variable
       _atomcache_size within the class.
       Added fragments manged property. Is a lazily built, cached entry, similar to residues.
    .. versionchanged:: 0.11.0
       AtomGroups can now be pickled and unpickled provided compatible Universes
       are available.
       The follow methods were changed to properties: ``indices``, ``masses``, ``charges``, ``names``,
       ``types``, ``radii``, ``resids``, ``resnames``, ``resnums``, ``segids``.
       Added ``altLocs`` and ``serials`` properties and setters.
       ``Torsions`` and ``torsion`` renamed to ``dihedral``.
       The ``bond``, ``angle``, ``dihedral`` and ``improper`` methods were removed and replaced
       with properties of the same name which return the corresponding object.
       Deprecated ``selectAtoms`` in favour of ``select_atoms``.
       Setters are now plural to match property names.
       Properties referring to residue (``resids``, ``resnames``, ``resnums``)
       or segment [``segids``] properties now yield arrays of length equal to
       ``n_atoms``
    .. versionchanged:: 0.11.1
       Added occupancies property
    """
    # for generalized __getitem__ __iter__ and __len__
    # (override _containername for ResidueGroup and SegmentGroup)
    _containername = "_atoms"
    _atomcache_size = 10000

    def __init__(self, atoms):
        if len(atoms) > 0:
            # __atoms property is effectively readonly
            # check that atoms is indexable:
            try:
                atoms[0]
                self.__atoms = atoms
            except TypeError:
                self.__atoms = list(atoms)
            # sanity check
            if not isinstance(self.__atoms[0], Atom):
                raise TypeError("atoms must be a Atom or a list of Atoms.")
        else:
            # empty AtomGroup
            self.__atoms = []

        # managed timestep object
        self._ts = None

        # caches:
        # - built on the fly when they are needed
        # - delete entry to invalidate
        self._cache = dict()

        # for generalized __getitem__ (override _containername for ResidueGroup and SegmentGroup)
        self._container = getattr(self, self._containername)
        # Define the Class that gets returned by getitem
        # Override this where return Class differs from Self (ie slicing Residue)
        self._cls = self.__class__

    @property
    def atoms(self):
        """:class:`AtomGroup` of all atoms in this group.

        If this is a :class:`AtomGroup` then it returns itself. Otherwise, it
        will return a new :class:`AtomGroup` based on all :class:`Atom`
        instances contained.

        Apply `:func:`list` to :attr:`atoms` or use :attr:`_atoms` if you
        really only need a list of individual :class:`Atom` instances.
        """
        # Cannot just return self because fails with inheritance from AtomGroup
        if type(self) == AtomGroup:
            return self
        return AtomGroup(self.__atoms)

    @property
    def _atoms(self):
        """a immutable list of references to the atoms in the group"""
        return self.__atoms

    # Universe pointer is important for Selections to work on groups
    @property
    def universe(self):
        """The universe to which the atoms belong (read-only)."""
        try:
            return self._atoms[0].universe
        except IndexError:
            raise NoDataError("Zero length AtomGroup have no Universe")

    def __len__(self):
        """Number of atoms in the group"""
        return len(self._container)

    def __iter__(self):
        return iter(self._container)

    def __getitem__(self, item):
        """Return element (index) or group (slicing).

        .. versionchanged:: 0.8
           :class:`ResidueGroup` and :class:`SegmentGroup`:
           return groups themselves and allow advanced slicing
        .. versionchanged:: 0.9.0
           This method now used by all subclasses.  These subclasses override
           :attr:`_cls` to define the returned class.
        .. versionchanged:: 0.10.0
           Now supports indexing via boolean numpy arrays
        """
        container = self._container
        cls = self._cls

        # consistent with the way list indexing/slicing behaves:
        if isinstance(item, int):
            try:
                return container[item]
            except IndexError:
                raise IndexError(
                    "Index {} is out of bounds for AtomGroup with size {}"
                    "".format(item, len(self)))
        elif isinstance(item, slice):
            return cls(container[item])
        elif isinstance(item, (np.ndarray, list)):
            # advanced slicing, requires array or list
            try:
                if isinstance(item[0], np.bool_):
                    item = np.arange(len(item))[item]
            except IndexError:  # zero length item
                pass
            return cls([container[i] for i in item])
        elif isinstance(item, str):
            return self._get_named_atom(item)
        else:
            raise TypeError("Cannot slice with type: {0}".format(type(item)))

    def __getattr__(self, name):
        try:
            return self._get_named_atom(name)
        except SelectionError:
            raise AttributeError("'{0}' object has no attribute '{1}'".format(
                    self.__class__.__name__, name))

    def _get_named_atom(self, name):
        """Get all atoms with name *name* in the current AtomGroup.

        For more than one atom it returns a list of :class:`Atom`
        instance. A single :class:`Atom` is returned just as such. If
        no atoms are found, a :exc:`SelectionError` is raised.

        .. versionadded:: 0.9.2
        """
        # There can be more than one atom with the same name
        atomlist = [atom for atom in self._atoms if name == atom.name]
        if len(atomlist) == 0:
            raise SelectionError("No atoms with name '{0}'".format(name))
        elif len(atomlist) == 1:
            return atomlist[0]  # XXX: keep this, makes more sense for names
        else:
            return AtomGroup(atomlist)  # XXX: but inconsistent (see residues and Issue 47)

    def __contains__(self, other):
        # If the number of atoms is very large, create a dictionary cache for lookup
        if len(self) > self._atomcache_size and not 'atoms' in self._cache:
            self._cache['atoms'] = dict(((x, None) for x in self.__atoms))
        try:
            return other in self._cache['atoms']
        except KeyError:
            return other in self._atoms

    def __add__(self, other):
        if not isinstance(other, (Atom, AtomGroup)):
            raise TypeError('Can only concatenate Atom or AtomGroup (not "{0}") to'
                            ' AtomGroup'.format(other.__class__.__name__))
        if (self and other) and (not self.universe is other.universe):
            raise ValueError("Can only add objects from the same Universe")
        if isinstance(other, AtomGroup):
            return AtomGroup(self._atoms + other._atoms)
        else:
            return AtomGroup(self._atoms + [other])

    def __repr__(self):
        return "<AtomGroup with {n_atoms} atoms>".format(
            n_atoms=len(self))

    def __getstate__(self):
        try:
            universe = self.universe
        except NoDataError:
            return None, None, None, None, None

        try: # We want to get the ChainReader case, where the trajectory has multiple filenames
            fname = universe.trajectory.filenames
        except AttributeError:
            fname = universe.trajectory.filename
        return (self.indices, universe.anchor_name, len(universe.atoms),
                universe.filename, fname)

    def __setstate__(self, state):
        indices, anchor_name, universe_n_atoms = state[:3]
        if indices is None:
            self.__init__([])
            return
        if np.max(indices) >= universe_n_atoms:
            raise ValueError("Trying to unpickle an inconsistent AtomGroup: "
                      "it has higher indices than the universe it refers to.")
        lookup_set = (MDAnalysis._anchor_universes if anchor_name is None
                      else MDAnalysis._named_anchor_universes)
        # Universes that go out-of-scope but haven't been garbage-collected
        #  yet may be (wrongly) brought back to life by unpickling onto them.
        # gc.collect() makes sure all is clean (Issue 487).
        gc.collect()
        for test_universe in lookup_set:
            if test_universe._matches_unpickling(*state[1:]):
                self.__init__(test_universe.atoms[indices]._atoms)
                return
        raise RuntimeError(("Couldn't find a suitable Universe to unpickle AtomGroup "
                "onto. (needed a universe with {}{} atoms, topology filename: '{}', and "
                "trajectory filename: '{}')").format(
                        "anchor_name: '{}', ".format(anchor_name) if anchor_name is not None else "",
                        *state[2:]))

    @property
    def n_atoms(self):
        """Total number of atoms in the group"""
        return len(self._atoms)

    @property
    def n_residues(self):
        """Total number of residues in the group"""
        return len(self.residues)

    @property
    def n_segments(self):
        """Total number of segments in the group"""
        return len(self.segments)

    @property
    @cached('indices')
    def indices(self):
        """Array of all :attr:`Atom.index` in the group.

        These indices are 0-based and can be used to directly index
        :attr:`Universe.atoms` or the coordinate array
        :attr:`MDAnalysis.coordinates.base.Timestep.positions`.

        .. Note::
           This property is read only

        .. versionchanged:: 0.11.0
           Now a property
        """
        return np.array([atom.index for atom in self._atoms])

    @property
    @cached('masses')
    def masses(self):
        """Array of atomic masses (as defined in the topology)

        .. versionchanged:: 0.11.0
           Now a property
        """
        return np.array([atom.mass for atom in self._atoms])

    @masses.setter
    def masses(self, new):
        self._clear_caches('masses')
        self.set_masses(new)

    def total_mass(self):
        """Total mass of the selection (masses are taken from the topology or guessed)."""
        return np.sum(self.masses, axis=0)

    totalMass = deprecate(total_mass, old_name='totalMass', new_name='total_mass')

    @property
    def occupancies(self):
        """Access occupancies of atoms

        If available can have a value between 0 and 1

        Returns
        -------
        ndarray
            occupancies for all atoms in AtomGroup

        .. versionadded:: 0.11.1
        """
        try:
            return self.universe.coord.data['occupancy'][self.indices]
        except KeyError:
            raise NoDataError('Timestep does not contain occupancy')

    @occupancies.setter
    def occupancies(self, new):
        try:
            self.universe.coord.data['occupancy'][self.indices] = new
        except KeyError:
            n_atoms = self.universe.coord.n_atoms
            self.universe.coord.data['occupancy'] = np.ones(n_atoms)
            self.universe.coord.data['occupancy'][self.indices] = new
    @property
    def charges(self):
        """Array of partial charges of the atoms (as defined in the topology)

        .. versionchanged:: 0.11.0
           Now a property
        """
        return np.array([atom.charge for atom in self._atoms])

    @charges.setter
    def charges(self, new):
        self.set_charges(new)

    def total_charge(self):
        """Sum of all partial charges (must be defined in topology)."""
        return np.sum(self.charges, axis=0)

    totalCharge = deprecate(total_charge, old_name='totalCharge', new_name='total_charge')

    @property
    def names(self):
        """Returns an array of atom names.

        .. versionchanged:: 0.8
           Returns a :class:`numpy.ndarray`
        .. versionchanged:: 0.11.0
           Now a property
        """
        return np.array([a.name for a in self._atoms])

    @names.setter
    def names(self, new):
        self.set_names(new)

    @property
    def types(self):
        """Returns an array of atom types.

        .. versionadded:: 0.9.0
        .. versionchanged:: 0.11.0
           Now a property
        """
        return np.array([a.type for a in self._atoms])

    @types.setter
    def types(self, new):
        self.set_types(new)

    @property
    def radii(self):
        """Array of atomic radii (as defined in the PQR file)

        .. versionchanged:: 0.11.0
           Now a property
        """
        return np.array([atom.radius for atom in self._atoms])

    @radii.setter
    def radii(self, new):
        self.set_radii(new)

    @property
    def bfactors(self):
        """Crystallographic B-factors (from PDB) in A**2.
        """
        return np.array([atom.bfactor for atom in self._atoms])

    @bfactors.setter
    def bfactors(self, new):
        self.set_bfactors(new)

    @property
    def altLocs(self):
        """numpy array of the altLocs for all atoms in this group

        .. versionadded:: 0.11.0
        """
        return np.array([atom.altLoc for atom in self._atoms])

    @altLocs.setter
    def altLocs(self, new):
        self.set_altlocs(new)

    @property
    def serials(self):
        """numpy array of the serials for all atoms in this group

        .. versionadded:: 0.11.0
        """
        return np.array([atom.serial for atom in self._atoms])

    @serials.setter
    def serials(self, new):
        self.set_serials(new)

    @property
    @cached('residues')
    def residues(self):
        """Read-only list of :class:`Residue` objects.

        A :class:`ResidueGroup` of all residues that contain atoms in
        this group.

        .. versionchanged:: 0.9.0
           Now returns strictly a ResidueGroup of the unique Residues that Atoms in this group
           belong to.
        """
        residues = []
        seen_residues = set()
        current_residue = None
        for atom in self._atoms:
            if atom.residue != current_residue and not atom.residue in seen_residues:
                residues.append(atom.residue)
                seen_residues.add(atom.residue)
            current_residue = atom.residue
        return ResidueGroup(residues)

    @property
    @cached('segments')
    def segments(self):
        """Read-only list of :class:`Segment` objects.

        A :class:`SegmentGroup` of all segments that contain atoms in this group.

        .. versionchanged:: 0.9.0
           Now strictly returns a :class:`SegmentGroup` of a set of
           the :class:`Segment` instances from this :class:`AtomGroup`

        """
        segments = []
        seen_segments = set()
        current_segment = None
        for atom in self._atoms:
            if atom.segment != current_segment and not atom.segment in seen_segments:
                segments.append(atom.segment)
                seen_segments.add(atom.segment)
            current_segment = atom.segment
        return SegmentGroup(segments)

    @property
    def resids(self):
        """Returns an array of residue numbers.

        .. versionchanged:: 0.8
           Returns a :class:`numpy.ndarray`
        .. versionchanged:: 0.11.0
           Now a property and returns array of length `len(self)`
        """
        return np.array([a.resid for a in self._atoms])

    @resids.setter
    def resids(self, new):
        self.set_resids(new)

    @property
    def resnames(self):
        """Returns an array of residue names.

        .. versionchanged:: 0.8
           Returns a :class:`numpy.ndarray`
        .. versionchanged:: 0.11.0
           Now a property and returns array of length `len(self)`
        """
        return np.array([a.resname for a in self._atoms])

    @resnames.setter
    def resnames(self, new):
        self.set_resnames(new)

    @property
    def resnums(self):
        """Returns an array of canonical residue numbers.

        .. versionadded:: 0.7.4
        .. versionchanged:: 0.8
           Returns a :class:`numpy.ndarray`
        .. versionchanged:: 0.11.0
           Now a property and returns array of length `len(self)`
        """
        return np.array([a.resnum for a in self._atoms])

    @resnums.setter
    def resnums(self, new):
        self.set_resnums(new)

    @property
    def segids(self):
        """Returns an array of segment names.

        .. versionchanged:: 0.8
           Returns a :class:`numpy.ndarray`
        .. versionchanged:: 0.11.0
           Now a property and returns array of length `len(self)`
        """
        return np.array([a.segid for a in self._atoms])

    @segids.setter
    def segids(self, new):
        self.set_segids(new)

    def sequence(self, **kwargs):
        """Returns the amino acid sequence.

        The format of the sequence is selected with the keyword *format*:

        ============== ============================================
        *format*       description
        ============== ============================================
        'SeqRecord'    :class:`Bio.SeqRecord.SeqRecord` (default)
        'Seq'          :class:`Bio.Seq.Seq`
        'string'       string
        ============== ============================================

        The sequence is returned by default (keyword ``format = 'SeqRecord'``)
        as a :class:`Bio.SeqRecord.SeqRecord` instance, which can then be
        further processed. In this case, all keyword arguments (such as the
        *id* string or the *name* or the *description*) are directly passed to
        :class:`Bio.SeqRecord.SeqRecord`.

        If the keyword *format* is set to ``'Seq'``, all *kwargs* are ignored
        and a :class:`Bio.Seq.Seq` instance is returned. The difference to the
        record is that the record also contains metadata and can be directly
        used as an input for other functions in :mod:`Bio`.

        If the keyword *format* is set to ``'string'``, all *kwargs* are
        ignored and a Python string is returned.

        .. rubric:: Example: Write FASTA file

        Use :func:`Bio.SeqIO.write`, which takes sequence records::

           import Bio.SeqIO

           # get the sequence record of a protein component of a Universe
           protein = u.select_atoms("protein")
           record = protein.sequence(id="myseq1", name="myprotein")

           Bio.SeqIO.write(record, "single.fasta", "fasta")

        A FASTA file with multiple entries can be written with ::

           Bio.SeqIO.write([record1, record2, ...], "multi.fasta", "fasta")

        :Keywords:
            *format*

                - ``"string"``: return sequence as a string of 1-letter codes
                - ``"Seq"``: return a :class:`Bio.Seq.Seq` instance
                - ``"SeqRecord"``: return a :class:`Bio.SeqRecord.SeqRecord`
                  instance

                Default is ``"SeqRecord"``

             *id*
                Sequence ID for SeqRecord (should be different for different
                sequences)
             *name*
                Name of the protein.
             *description*
                Short description of the sequence.
             *kwargs*
                Any other keyword arguments that are understood by
                :class:`Bio.SeqRecord.SeqRecord`.

        :Raises: :exc:`ValueError` if a residue name cannot be converted to a
                 1-letter IUPAC protein amino acid code; make sure to only
                 select protein residues. Raises :exc:`TypeError` if an unknown
                 *format* is selected.

        .. versionadded:: 0.9.0
        """
        import Bio.Seq
        import Bio.SeqRecord
        import Bio.Alphabet
        formats = ('string', 'Seq', 'SeqRecord')

        format = kwargs.pop("format", "SeqRecord")
        if format not in formats:
            raise TypeError("Unknown format='{0}': must be one of: {1}".format(
                    format, ", ".join(formats)))
        try:
            sequence = "".join([util.convert_aa_code(r) for r in self.residues.resnames])
        except KeyError as err:
            raise ValueError("AtomGroup contains a residue name '{0}' that "
                             "does not have a IUPAC protein 1-letter "
                             "character".format(err.message))
        if format == "string":
            return sequence
        seq = Bio.Seq.Seq(sequence, alphabet=Bio.Alphabet.IUPAC.protein)
        if format == "Seq":
            return seq
        return Bio.SeqRecord.SeqRecord(seq, **kwargs)

    @property
    @cached('fragments')
    def fragments(self):
        """Read-only list of fragments.

        Contains all fragments that any Atom in this AtomGroup is part of, the contents of
        the fragments may extend beyond the contents of this AtomGroup.

        .. versionadded 0.9.0
        """
        return tuple(set(a.fragment for a in self._atoms))

    def guess_bonds(self, vdwradii=None):
        """Guess all the bonds that exist within this AtomGroup and add to Universe.

        :Keywords:
          *vdwradii*
            Pass a dict relating atom types to vdwradii.

        .. SeeAlso::
           :func:`MDAnalysis.topology.core.guess_bonds`

        .. versionadded:: 0.10.0
        """
        from ..topology.core import (guess_bonds,
                                     guess_angles,
                                     guess_dihedrals)

        b = guess_bonds(self.atoms, self.atoms.positions, vdwradii=vdwradii)

        # eliminate bonds that already exist
        # have to compare indices not bond objects as same bond which is True and False
        # will hash differently.
        existing_bonds = set(self.universe.bonds.to_indices())
        new_b = set(b).difference(existing_bonds)
        bgroup = top.TopologyGroup.from_indices(new_b, self.universe.atoms,
                                                            bondclass=top.Bond, guessed=True)
        self.universe.bonds += bgroup
        self._clear_caches('bonds')

        a = guess_angles(self.bonds)
        existing_angles = set(self.universe.angles.to_indices())
        new_a = set(a).difference(existing_angles)
        agroup = top.TopologyGroup.from_indices(new_a, self.universe.atoms,
                                                            bondclass=top.Angle, guessed=True)
        self.universe.angles += agroup

        self._clear_caches('angles')

        t = guess_dihedrals(self.angles)
        existing_t = set(self.universe.dihedrals.to_indices())
        new_t = set(t).difference(existing_t)
        tgroup = top.TopologyGroup.from_indices(new_t, self.universe.atoms,
                                                            bondclass=top.Dihedral, guessed=True)
        self.universe.dihedrals += tgroup
        self._clear_caches('dihedrals')

    @property
    @cached('bonds')
    def bonds(self):
        """All the bonds in this AtomGroup

        Note that these bonds might extend out of the AtomGroup, to select
        only bonds which are entirely contained by the AtomGroup use
        u.bonds.atomgroup_intersection(ag, strict=True)

        .. versionadded:: 0.9.0
        .. versionchanged:: 0.10.0
           Now always returns a (possibly empty) TopologyGroup
        """
        mybonds = [b for a in self._atoms for b in a.bonds]

        return top.TopologyGroup(mybonds)

    @property
    @cached('angles')
    def angles(self):
        """All the angles in this AtomGroup

        Note that these angles might extend out of the AtomGroup, to select
        only angles which are entirely contained by the AtomGroup use
        u.angles.atomgroup_intersection(ag, strict=True)

        .. versionadded:: 0.9.0
        .. versionchanged:: 0.10.0
           Now always returns a (possibly empty) TopologyGroup
        """
        mybonds = [b for a in self._atoms for b in a.angles]

        return top.TopologyGroup(mybonds)

    @property
    @cached('dihedrals')
    def dihedrals(self):
        """All the dihedrals in this AtomGroup

        Note that these dihedrals might extend out of the AtomGroup, to select
        only dihedrals which are entirely contained by the AtomGroup use
        u.dihedrals.atomgroup_intersection(ag, strict=True)

        .. versionadded:: 0.9.0
        .. versionchanged:: 0.10.0
           Now always returns a (possibly empty) TopologyGroup
        """
        mybonds = [b for a in self._atoms for b in a.dihedrals]

        return top.TopologyGroup(mybonds)

    @property
    @cached('impropers')
    def impropers(self):
        """All the improper dihedrals in this AtomGroup

        Note that these improper dihedrals might extend out of the AtomGroup,
        to select only dihedrals which are entirely contained by the AtomGroup use
        u.impropers.atomgroup_intersection(ag, strict=True)

        .. versionadded:: 0.9.0
        .. versionchanged:: 0.10.0
           Now always returns a (possibly empty) TopologyGroup
        """
        mybonds = [b for a in self._atoms for b in a.impropers]

        return top.TopologyGroup(mybonds)

    def _set_attribute(self, groupname, name, value, **kwargs):
        """Set attribute *name* to *value* for all elements in *groupname*.

        *groupname* can be *atoms*, *residues*, *segments. ``getattr(self,
        groupname)`` should produce one of the groups in the hierarchy.

        If *value* is a sequence of the same length as the group then each
        element's attribute *name* is set to the corresponding value. If
        *value* is neither of length 1 (or a scalar) nor of the length of the
        group then a :exc:`ValueError` is raised.

        A cache entry ``_cache[groupname]`` is deleted if it exists.

        :Keywords:

          *conversion*
               function such as :func:`str` or :func:`int` that converts the
               argument. ``None`` passes it through unchanged [``None``]

          *cache*
               alternative identifier for the cache, instead of *groupname*

        .. versionadded:: 0.8
        """
        values = util.asiterable(value)
        group = getattr(self, groupname)
        conversion = kwargs.pop('conversion', None)
        cache = kwargs.pop('cache', groupname)
        if not conversion:
            conversion = lambda x: x
        if len(values) == 1:
            for x in group:
                setattr(x, name, conversion(values[0]))
        elif len(group) == len(values):
            for x, value in itertools.izip(group, values):
                setattr(x, name, conversion(value))
        else:
            raise ValueError("set_{0}: can only set all atoms to a single value or each atom to a distinct one "
                             "but len(atoms)={1} whereas len(value)={2}".format(groupname, len(group), len(values)))
        self._clear_caches(cache)

        # big hammer... if we find the time, use this in a more surgical fashion.
        #self.atoms._rebuild_caches()
        #if self.atoms is not self.universe.atoms:
        #    self.universe.atoms._rebuild_caches()

    def _set_atoms(self, name, value, **kwargs):
        """Set attribute *name* to *value* for all atoms in the :class:`AtomGroup`.

        If *value* is a sequence of the same length as the :class:`AtomGroup`
        then each :class:`Atom`'s property *name* is set to the corresponding value. If
        *value* is neither of length 1 (or a scalar) nor of the length of the
        :class:`AtomGroup` then a :exc:`ValueError` is raised.

        .. versionadded:: 0.7.4
        .. versionchanged:: 0.8
           Can set atoms to distinct values by providing a sequence or iterable.
        """
        self._set_attribute("atoms", name, value, **kwargs)

    # override for ResidueGroup, SegmentGroup accordingly
    set = _set_atoms

    def set_occupancies(self, occupancies):
        """Set the occupancy for *all atoms* in the AtomGroup

        If *value* is a sequence of the same length as the :class:`AtomGroup`
        then each :attr:`Atom.name` is set to the corresponding value. If
        *value* is neither of length 1 (or a scalar) nor of the length of the
        :class:`AtomGroup` then a :exc:`ValueError` is raised.

        .. versionadded:: 0.11.1
        """
        self.occupancies = occupancies

    def set_names(self, name):
        """Set the atom names to string for *all atoms* in the AtomGroup.

        If *value* is a sequence of the same length as the :class:`AtomGroup`
        then each :attr:`Atom.name` is set to the corresponding value. If
        *value* is neither of length 1 (or a scalar) nor of the length of the
        :class:`AtomGroup` then a :exc:`ValueError` is raised.

        .. versionadded:: 0.7.4
        .. versionchanged:: 0.8
           Can set atoms to distinct values by providing a sequence or iterable.
        .. versionchanged:: 0.11.0
           Made plural to make consistent with corresponding property
        """
        self.set("name", name, conversion=str)

    set_name = deprecate(set_names, old_name='set_name', new_name='set_names')

    def set_resids(self, resid):
        """Set the resids to integer *resid* for **all atoms** in the :class:`AtomGroup`.

        If *resid* is a sequence of the same length as the :class:`AtomGroup`
        then each :attr:`Atom.resid` is set to the corresponding value together
        with the :attr:`Residue.id` of the residue the atom belongs to. If
        *value* is neither of length 1 (or a scalar) nor of the length of the
        :class:`AtomGroup` then a :exc:`ValueError` is raised.

        .. Note::

           Changing resids can change the topology.

           Assigning the same *resid* to multiple residues will
           **merge** these residues. Assigning different *resid* to
           atoms in the same residue will **split** a residue (and
           potentially merge with another one).

        .. versionadded:: 0.7.4
        .. versionchanged:: 0.7.5
           Also changes the residues.
        .. versionchanged:: 0.8
           Can set atoms and residues to distinct values by providing a
           sequence or iterable and can change the topology via
           :func:`MDAnalysis.topology.core.build_residues`.
        .. versionchanged:: 0.11.0
           Made plural to make consistent with corresponding property
        """
        from MDAnalysis.topology.core import build_residues

        self.set("resid", resid, conversion=int)
        # Note that this also automagically updates THIS AtomGroup;
        # the side effect of build_residues(self.atoms) is to update all Atoms!!!!
        self._fill_cache('residues', ResidueGroup(build_residues(self.atoms)))

        # make sure to update the whole universe: the Atoms are shared but
        # ResidueGroups are not
        if self.atoms is not self.universe.atoms:
            self.universe.atoms._fill_cache(
                    'residues',
                    ResidueGroup(build_residues(self.universe.atoms)))

    set_resid = deprecate(set_resids, old_name='set_resid', new_name='set_resids')

    def set_resnums(self, resnum):
        """Set the resnums to *resnum* for **all atoms** in the :class:`AtomGroup`.

        If *resnum* is a sequence of the same length as the :class:`AtomGroup`
        then each :attr:`Atom.resnum` is set to the corresponding value together
        with the :attr:`Residue.resnum` of the residue the atom belongs to. If
        *value* is neither of length 1 (or a scalar) nor of the length of the
        :class:`AtomGroup` then a :exc:`ValueError` is raised.

        .. Note::

           Changing *resnum* will not affect the topology: you can have
           multiple residues with the same *resnum*.

        .. SeeAlso:: :meth:`set_resid`

        .. versionadded:: 0.7.4
        .. versionchanged:: 0.7.5
           Also changes the residues.
        .. versionchanged:: 0.8
           Can set atoms and residues to distinct values by providing a sequence or iterable.
        .. versionchanged:: 0.11.0
           Made plural to make consistent with corresponding property
        """
        self.set("resnum", resnum)

    set_resnum = deprecate(set_resnums, old_name='set_resnum', new_name='set_resnums')

    def set_resnames(self, resname):
        """Set the resnames to string *resname* for **all atoms** in the :class:`AtomGroup`.

        If *resname* is a sequence of the same length as the :class:`AtomGroup`
        then each :attr:`Atom.resname` is set to the corresponding value together
        with the :attr:`Residue.name` of the residue the atom belongs to. If
        *value* is neither of length 1 (or a scalar) nor of the length of the
        :class:`AtomGroup` then a :exc:`ValueError` is raised.

        .. versionadded:: 0.7.4
        .. versionchanged:: 0.7.5
           Also changes the residues.
        .. versionchanged:: 0.8
           Can set atoms and residues to distinct values by providing a sequence or iterable.
        .. versionchanged:: 0.11.0
           Made plural to make consistent with corresponding property
        """
        self.set("resname", resname, conversion=str)

    set_resname = deprecate(set_resnames, old_name='set_resname', new_name='set_resnames')

    def set_segids(self, segid):
        """Set the segids to *segid* for all atoms in the :class:`AtomGroup`.

        If *segid* is a sequence of the same length as the :class:`AtomGroup`
        then each :attr:`Atom.segid` is set to the corresponding value together
        with the :attr:`Segment.id` of the residue the atom belongs to. If
        *value* is neither of length 1 (or a scalar) nor of the length of the
        :class:`AtomGroup` then a :exc:`ValueError` is raised.

        .. Note::

           :meth:`set_segid` can change the topology.

        .. versionadded:: 0.7.4
        .. versionchanged:: 0.8
           Can set atoms and residues to distinct values by providing a sequence or iterable.
        .. versionchanged:: 0.11.0
           Stale caches are problematic; though it can be expensive, changing segid
           results in Segments being regenerated
        .. versionchanged:: 0.11.0
           Made plural to make consistent with corresponding property
        """
        from MDAnalysis.topology.core import build_segments

        self.set("segid", segid, conversion=str)

        # also updates convenience handles for segments in universe
        segments = self.universe._build_segments()

        # Note that this also automagically updates THIS AtomGroup;
        # the side effect of build_residues(self.atoms) is to update all Atoms!!!!
        self._fill_cache('segments', SegmentGroup(segments))

        # make sure to update the whole universe: the Atoms are shared but
        # ResidueGroups are not
        if self.atoms is not self.universe.atoms:
            self.universe.atoms._fill_cache(
                    'segments',
                    SegmentGroup(segments))

    set_segid = deprecate(set_segids, old_name='set_segid', new_name='set_segids')

    def set_masses(self, mass):
        """Set the atom masses to float *mass* for **all atoms** in the AtomGroup.

        If *mass* is a sequence of the same length as the :class:`AtomGroup`
        then each :attr:`Atom.mass` is set to the corresponding value. If
        *value* is neither of length 1 (or a scalar) nor of the length of the
        :class:`AtomGroup` then a :exc:`ValueError` is raised.

        .. versionadded:: 0.7.4
        .. versionchanged:: 0.8
           Can set atoms and residues to distinct values by providing a sequence or iterable.
        .. versionchanged:: 0.11.0
           Made plural to make consistent with corresponding property
        """
        self.set("mass", mass, conversion=float, cache="masses")

    set_mass = deprecate(set_masses, old_name='set_mass', new_name='set_masses')

    def set_types(self, atype):
        """Set the atom types to *atype* for **all atoms** in the AtomGroup.

        If *atype* is a sequence of the same length as the :class:`AtomGroup`
        then each :attr:`Atom.atype` is set to the corresponding value. If
        *value* is neither of length 1 (or a scalar) nor of the length of the
        :class:`AtomGroup` then a :exc:`ValueError` is raised.

        .. versionadded:: 0.7.4
        .. versionchanged:: 0.8
           Can set atoms and residues to distinct values by providing a sequence or iterable.
        .. versionchanged:: 0.11.0
           Made plural to make consistent with corresponding property
        """
        self.set("type", atype)

    set_type = deprecate(set_types, old_name='set_type', new_name='set_types')

    def set_charges(self, charge):
        """Set the partial charges to float *charge* for **all atoms** in the AtomGroup.

        If *charge* is a sequence of the same length as the :class:`AtomGroup`
        then each :attr:`Atom.charge` is set to the corresponding value. If
        *value* is neither of length 1 (or a scalar) nor of the length of the
        :class:`AtomGroup` then a :exc:`ValueError` is raised.

        .. versionadded:: 0.7.4
        .. versionchanged:: 0.8
           Can set atoms and residues to distinct values by providing a sequence or iterable.
        .. versionchanged:: 0.11.0
           Made plural to make consistent with corresponding property
        """
        self.set("charge", charge, conversion=float)

    set_charge = deprecate(set_charges, old_name='set_charge', new_name='set_charges')

    def set_radii(self, radius):
        """Set the atom radii to float *radius* for **all atoms** in the AtomGroup.

        If *radius* is a sequence of the same length as the :class:`AtomGroup`
        then each :attr:`Atom.radius` is set to the corresponding value. If
        *value* is neither of length 1 (or a scalar) nor of the length of the
        :class:`AtomGroup` then a :exc:`ValueError` is raised.

        .. versionadded:: 0.7.4
        .. versionchanged:: 0.8
           Can set atoms and residues to distinct values by providing a sequence or iterable.
        .. versionchanged:: 0.11.0
           Made plural to make consistent with corresponding property
        """
        self.set("radius", radius, conversion=float)

    set_radius = deprecate(set_radii, old_name='set_radius', new_name='set_radii')

    def set_bfactors(self, bfactor):
        """Set the atom bfactors to float *bfactor* for **all atoms** in the AtomGroup.

        If *bfactor* is a sequence of the same length as the :class:`AtomGroup`
        then each :attr:`Atom.bfactor` is set to the corresponding value. If
        *value* is neither of length 1 (or a scalar) nor of the length of the
        :class:`AtomGroup` then a :exc:`ValueError` is raised.

        .. versionadded:: 0.7.4
        .. versionchanged:: 0.8
           Can set atoms and residues to distinct values by providing a sequence or iterable.
        .. versionchanged:: 0.11.0
           Made plural to make consistent with corresponding property
        """
        self.set("bfactor", bfactor, conversion=float)

    set_bfactor = deprecate(set_bfactors, old_name='set_bfactor', new_name='set_bfactors')

    def set_altLocs(self, altLoc):
        """Set the altLocs to *altLoc for **all atoms** in the AtomGroup.

        If *altLoc* is a sequence of the same length as the :class:`AtomGroup`
        then each :attr:`Atom.altLoc` is set to the corresponding value. If
        *value* is neither of length 1 (or a scalar) nor of the length of the
        :class:`AtomGroup` then a :exc:`ValueError` is raised.

        .. versionadded:: 0.11.0
        """
        self.set("altLoc", altLoc, conversion=str)

    set_altLoc = deprecate(set_altLocs, old_name='set_altLoc', new_name='set_altLocs')

    def set_serials(self, serial):
        """Set the serials to *serial* for **all atoms** in the AtomGroup.

        If *serial* is a sequence of the same length as the :class:`AtomGroup`
        then each :attr:`Atom.serial` is set to the corresponding value. If
        *value* is neither of length 1 (or a scalar) nor of the length of the
        :class:`AtomGroup` then a :exc:`ValueError` is raised.

        .. versionadded:: 0.11.0
        """
        self.set("serial", serial, conversion=int)

    set_serial = deprecate(set_serials, old_name='set_serial', new_name='set_serials')

    def center_of_geometry(self, **kwargs):
        """Center of geometry (also known as centroid) of the selection.

        :Keywords:
          *pbc*
            ``True``: Move all atoms within the primary unit cell before calculation [``False``]

        .. Note::
            The :class:`MDAnalysis.core.flags` flag *use_pbc* when set to ``True`` allows the *pbc*
            flag to be used by default.

        .. versionchanged:: 0.8 Added *pbc* keyword
        """
        pbc = kwargs.pop('pbc', MDAnalysis.core.flags['use_pbc'])
        if pbc:
            return np.sum(self.pack_into_box(inplace=False), axis=0) / self.n_atoms
        else:
            return np.sum(self.positions, axis=0) / self.n_atoms

    centerOfGeometry = deprecate(center_of_geometry, old_name='centerOfGeometry',
                                 new_name='center_of_geometry')

    centroid = center_of_geometry

    def center_of_mass(self, **kwargs):
        """Center of mass of the selection.

        :Keywords:
          *pbc*
            ``True``: Move all atoms within the primary unit cell before calculation [``False``]

        .. Note::
            The :class:`MDAnalysis.core.flags` flag *use_pbc* when set to ``True`` allows the *pbc*
            flag to be used by default.

        .. versionchanged:: 0.8 Added *pbc* keyword
        """
        pbc = kwargs.pop('pbc', MDAnalysis.core.flags['use_pbc'])
        if pbc:
            return np.sum(self.pack_into_box(inplace=False) * self.masses[:, np.newaxis],
                             axis=0) / self.total_mass()
        else:
            return np.sum(self.positions * self.masses[:, np.newaxis], axis=0) / self.total_mass()

    centerOfMass = deprecate(center_of_mass, old_name='centerOfMass', new_name='center_of_mass')

    def radius_of_gyration(self, **kwargs):
        """Radius of gyration.

        :Keywords:
          *pbc*
            ``True``: Move all atoms within the primary unit cell before calculation [``False``]

        .. Note::
            The :class:`MDAnalysis.core.flags` flag *use_pbc* when set to ``True`` allows the *pbc*
            flag to be used by default.

        .. versionchanged:: 0.8 Added *pbc* keyword
        """
        pbc = kwargs.pop('pbc', MDAnalysis.core.flags['use_pbc'])
        masses = self.masses
        if pbc:
            recenteredpos = self.pack_into_box(inplace=False) - self.center_of_mass(pbc=True)
        else:
            recenteredpos = self.positions - self.center_of_mass(pbc=False)
        rog_sq = np.sum(masses * np.sum(np.power(recenteredpos, 2), axis=1)) / self.total_mass()
        return np.sqrt(rog_sq)

    radiusOfGyration = deprecate(radius_of_gyration, old_name='radiusOfGyration', new_name='radius_of_gyration')

    def shape_parameter(self, **kwargs):
        """Shape parameter.

        See [Dima2004]_ for background information.

        :Keywords:
          *pbc*
            ``True``: Move all atoms within the primary unit cell before calculation [``False``]

        .. Note::
            The :class:`MDAnalysis.core.flags` flag *use_pbc* when set to ``True`` allows the *pbc*
            flag to be used by default.

        .. versionadded:: 0.7.7
        .. versionchanged:: 0.8 Added *pbc* keyword
        """
        pbc = kwargs.pop('pbc', MDAnalysis.core.flags['use_pbc'])
        masses = self.masses
        if pbc:
            recenteredpos = self.pack_into_box(inplace=False) - self.center_of_mass(pbc=True)
        else:
            recenteredpos = self.positions - self.center_of_mass(pbc=False)
        tensor = np.zeros((3, 3))
        for x in xrange(recenteredpos.shape[0]):
            tensor += masses[x] * np.outer(recenteredpos[x, :],
                                              recenteredpos[x, :])
        tensor /= self.total_mass()
        eig_vals = np.linalg.eigvalsh(tensor)
        shape = 27.0 * np.prod(eig_vals - np.mean(eig_vals)) / np.power(np.sum(eig_vals), 3)
        return shape

    shapeParameter = deprecate(shape_parameter, old_name='shapeParameter', new_name='shape_parameter')

    def asphericity(self, **kwargs):
        """Asphericity.

        See [Dima2004]_ for background information.

        :Keywords:
          *pbc*
            ``True``: Move all atoms within primary unit cell before calculation [``False``]

        .. Note::
            The :class:`MDAnalysis.core.flags` flag *use_pbc* when set to ``True`` allows the *pbc*
            flag to be used by default.

        .. versionadded:: 0.7.7
        .. versionchanged:: 0.8 Added *pbc* keyword
        """
        pbc = kwargs.pop('pbc', MDAnalysis.core.flags['use_pbc'])
        masses = self.masses
        if pbc:
            recenteredpos = self.pack_into_box(inplace=False) - self.center_of_mass(pbc=True)
        else:
            recenteredpos = self.positions - self.center_of_mass(pbc=False)
        tensor = np.zeros((3, 3))
        for x in xrange(recenteredpos.shape[0]):
            tensor += masses[x] * np.outer(recenteredpos[x, :],
                                              recenteredpos[x, :])
        tensor /= self.total_mass()
        eig_vals = np.linalg.eigvalsh(tensor)
        shape = (3.0 / 2.0) * np.sum(np.power(eig_vals - np.mean(eig_vals), 2)) / np.power(
            np.sum(eig_vals), 2)
        return shape

    def moment_of_inertia(self, **kwargs):
        """Tensor of inertia as 3x3 numpy array.

        :Keywords:
          *pbc*
            ``True``: Move all atoms within the primary unit cell before calculation [``False``]

        .. Note::
            The :class:`MDAnalysis.core.flags` flag *use_pbc* when set to ``True`` allows the *pbc*
            flag to be used by default.

        .. versionchanged:: 0.8 Added *pbc* keyword
        """
        pbc = kwargs.pop('pbc', MDAnalysis.core.flags['use_pbc'])
        # Convert to local coordinates
        if pbc:
            pos = self.pack_into_box(inplace=False) - self.center_of_mass(pbc=True)
        else:
            pos = self.positions - self.center_of_mass(pbc=False)

        masses = self.masses
        # Create the inertia tensor
        # m_i = mass of atom i
        # (x_i, y_i, z_i) = pos of atom i
        # Ixx = sum(m_i*(y_i^2+z_i^2));
        # Iyy = sum(m_i*(x_i^2+z_i^2));
        # Izz = sum(m_i*(x_i^2+y_i^2))
        # Ixy = Iyx = -1*sum(m_i*x_i*y_i)
        # Ixz = Izx = -1*sum(m_i*x_i*z_i)
        # Iyz = Izy = -1*sum(m_i*y_i*z_i)
        tens = np.zeros((3,3), dtype=np.float64)
        # xx
        tens[0][0] = (masses * (pos[:,1] * pos[:,1] + pos[:,2] * pos[:,2])).sum()
        # xy & yx
        tens[0][1] = tens[1][0] = - (masses * pos[:,0] * pos[:,1]).sum()
        # xz & zx
        tens[0][2] = tens[2][0] = - (masses * pos[:,0] * pos[:,2]).sum()
        # yy
        tens[1][1] = (masses * (pos[:,0] * pos[:,0] + pos[:,2] * pos[:,2])).sum()
        # yz + zy
        tens[1][2] = tens[2][1] = - (masses * pos[:,1] * pos[:,2]).sum()
        # zz
        tens[2][2] = (masses * (pos[:,0] * pos[:,0] + pos[:,1] * pos[:,1])).sum()

        return tens

    momentOfInertia = deprecate(moment_of_inertia, old_name='momentOfInertia',
                                new_name='moment_of_inertia')

    def bbox(self, **kwargs):
        """Return the bounding box of the selection.

        The lengths A,B,C of the orthorhombic enclosing box are ::

          L = AtomGroup.bbox()
          A,B,C = L[1] - L[0]

        :Keywords:
          *pbc*
            ``True``: Move all atoms within the primary unit cell before calculation [``False``]

        .. Note::
            The :class:`MDAnalysis.core.flags` flag *use_pbc* when set to ``True`` allows the *pbc*
            flag to be used by default.

        :Returns: [[xmin, ymin, zmin], [xmax, ymax, zmax]]

        .. versionadded:: 0.7.2
        .. versionchanged:: 0.8 Added *pbc* keyword
        """
        pbc = kwargs.pop('pbc', MDAnalysis.core.flags['use_pbc'])
        if pbc:
            x = self.pack_into_box(inplace=False)
        else:
            x = self.coordinates()
        return np.array([x.min(axis=0), x.max(axis=0)])

    def bsphere(self, **kwargs):
        """Return the bounding sphere of the selection.

        The sphere is calculated relative to the centre of geometry.

        :Keywords:
          *pbc*
            ``True``: Move all atoms within primary unit cell before calculation [``False``]

        .. Note::
            The :class:`MDAnalysis.core.flags` flag *use_pbc* when set to ``True`` allows the *pbc*
            flag to be used by default.

        :Returns: `(R, [xcen,ycen,zcen])`

        .. versionadded:: 0.7.3
        .. versionchanged:: 0.8 Added *pbc* keyword
        """
        pbc = kwargs.pop('pbc', MDAnalysis.core.flags['use_pbc'])
        if pbc:
            x = self.pack_into_box(inplace=False)
            centroid = self.center_of_geometry(pbc=True)
        else:
            x = self.coordinates()
            centroid = self.center_of_geometry(pbc=False)
        R = np.sqrt(np.max(np.sum(np.square(x - centroid), axis=1)))
        return R, centroid

    @property
    def bond(self):
        """This AtomGroup represented as a Bond object

        :Returns:
          A :class:`MDAnalysis.core.topologyobjects.Bond` object

        :Raises:
          `ValueError` if the AtomGroup is not length 2

        .. versionadded:: 0.11.0
        """
        if len(self) != 2:
            raise ValueError("bond only makes sense for a group with exactly 2 atoms")
        return top.Bond(self.atoms)

    @property
    def angle(self):
        """This AtomGroup represented as an Angle object

        :Returns:
          A :class:`MDAnalysis.core.topologyobjects.Angle` object

        :Raises:
          `ValueError` if the AtomGroup is not length 3

        .. versionadded:: 0.11.0
        """
        if len(self) != 3:
            raise ValueError("angle only makes sense for a group with exactly 3 atoms")
        return top.Angle(self.atoms)

    @property
    def dihedral(self):
        """This AtomGroup represented as a Dihedral object

        :Returns:
          A :class:`MDAnalysis.core.topologyobjects.Dihedral` object

        :Raises:
          `ValueError` if the AtomGroup is not length 4

        .. versionadded:: 0.11.0
        """
        if len(self) != 4:
            raise ValueError("dihedral only makes sense for a group with exactly 4 atoms")

        return top.Dihedral(self.atoms)

    @property
    def improper(self):
        """This AtomGroup represented as an ImproperDihedral object

        :Returns:
          A :class:`MDAnalysis.core.topologyobjects.ImproperDihedral` object

        :Raises:
          `ValueError` if the AtomGroup is not length 4

        .. versionadded:: 0.11.0
        """
        if len(self) != 4:
            raise ValueError("improper only makes sense for a group with exactly 4 atoms")

        return top.ImproperDihedral(self.atoms)

    def principal_axes(self, **kwargs):
        """Calculate the principal axes from the moment of inertia.

        e1,e2,e3 = AtomGroup.principal_axes()

        The eigenvectors are sorted by eigenvalue, i.e. the first one
        corresponds to the highest eigenvalue and is thus the first principal axes.

        :Keywords:
          *pbc*
            ``True``: Move all atoms within primary unit cell before calculation

        .. Note::
            The :class:`MDAnalysis.core.flags` flag *use_pbc* when set to ``True`` allows the *pbc*
            flag to be used by default.

        :Returns: numpy.array ``v`` with ``v[0]`` as first, ``v[1]`` as second,
                  and ``v[2]`` as third eigenvector.

        .. versionchanged:: 0.8 Added *pbc* keyword
        """
        pbc = kwargs.pop('pbc', MDAnalysis.core.flags['use_pbc'])
        if pbc:
            eigenval, eigenvec = eig(self.moment_of_inertia(pbc=True))
        else:
            eigenval, eigenvec = eig(self.moment_of_inertia(pbc=False))
        # Sort
        indices = np.argsort(eigenval)
        # Return transposed in more logical form. See Issue 33.
        return eigenvec[:, indices].T

    principalAxes = deprecate(principal_axes, old_name='principalAxes', new_name='principal_axes')

    def get_positions(self, ts=None, copy=False, dtype=np.float32):
        """Get a numpy array of the coordinates.

        :Keywords:
           *ts*
               If *ts* is provided then positions are read from that
               :class:`~MDAnalysis.coordinates.base.Timestep` instead of
               the one from the current trajectory belonging to this universe.
               The *ts* is indexed with the indices returned by
               :meth:`~AtomGroup.indices` and it is the user's responsibility
               to provide a time step that has the appropriate dimensions.
           *copy*
               ``True``: always make a copy (slow), ``False``: Try to
               return a array view or reference (faster); note that for
               passing coordinates to C-code it can be necessary to use
               a copy [``False``]
           *dtype*
               numpy Data type of the array; the default is usually
               entirely appropriate. Most C-code actually requires the
               default  [:class:`np.float32`]

        Coordinates can also be directly obtained from the attribute
        :attr:`~AtomGroup.positions`.

        Coordinates can be directly set with :meth:`~AtomGroup.set_positions` or
        by assigning to :attr:`~AtomGroup.positions`.

        This method is identical with :meth:`~AtomGroup.coordinates` but named
        differently for symmetry with with :meth:`~AtomGroup.set_positions`.

        .. versionadded:: 0.7.6
        """
        if ts is None:
            ts = self.universe.trajectory.ts
        return np.array(ts.positions[self.indices], copy=copy, dtype=dtype)

    coordinates = get_positions
    """Np array of the coordinates.

    .. SeeAlso:: :attr:`~AtomGroup.positions` and :meth:`~AtomGroup.get_positions`

    .. deprecated:: 0.7.6
       In new scripts use :meth:`AtomGroup.get_positions` preferrably.
    """
    # coordinates() should NOT be removed as it has been used in many scripts,
    # MDAnalysis itself, and in the paper

    def set_positions(self, coords, ts=None):
        """Set the positions for all atoms in the group.

        :Arguments:
           *coords*
               a Nx3 :class:`numpy.ndarray` where N is the number of
               atoms in this atom group.

        :Keywords:
           *ts*
              :class:`~MDAnalysis.coordinates.base.Timestep`, defaults
              to ``None`` and then the current time step is used.

        .. Note::

           If the group contains N atoms and *coord* is a single point (i.e. an
           array of length 3) then all N atom positions are set to *coord* (due
           to numpy's broadcasting rules), as described for
           :attr:`~AtomGroup.positions`.

        See also :meth:`~AtomGroup.get_positions` and attribute access through
        :attr:`~AtomGroup.positions`.

        .. versionadded:: 0.7.6
        """
        if ts is None:
            ts = self.universe.trajectory.ts
        ts.positions[self.indices, :] = coords

    positions = property(get_positions, set_positions,
                         doc="""
                Coordinates of the atoms in the AtomGroup.

                The positions can be changed by assigning an array of the appropriate
                shape, i.e. either Nx3 to assign individual coordinates or 3, to assign
                the *same* coordinate to all atoms (e.g. ``ag.positions = array([0,0,0])``
                will move all particles to the origin).

                For more control use the :meth:`~AtomGroup.get_positions` and
                :meth:`~AtomGroup.set_positions` methods.

                .. versionadded:: 0.7.6""")

    def get_velocities(self, ts=None, copy=False, dtype=np.float32):
        """numpy array of the velocities.

        Raises a :exc:`NoDataError` if the underlying
        :class:`~MDAnalysis.coordinates.base.Timestep` does not contain
        :attr:`~MDAnalysis.coordinates.base.Timestep.velocities`.

        See also :meth:`AtomGroup.set_velocities` and attribute access through
        :attr:`AtomGroup.velocities`.

        .. versionadded:: 0.7.6
        """
        if ts is None:
            ts = self.universe.trajectory.ts
        try:
            return np.array(ts.velocities[self.indices], copy=copy, dtype=dtype)
        except (AttributeError, NoDataError):
            raise NoDataError("Timestep does not contain velocities")

    def set_velocities(self, v, ts=None):
        """Assign the velocities *v* to the timestep.

        Raises a :exc:`NoDataError` if the underlying
        :class:`~MDAnalysis.coordinates.base.Timestep` does not contain
        :attr:`~MDAnalysis.coordinates.base.Timestep.velocities`.

        See also :meth:`AtomGroup.get_velocities` and :attr:`AtomGroup.velocities` for
        attribute access.

        .. versionadded:: 0.7.6
        """
        if ts is None:
            ts = self.universe.trajectory.ts
        try:
            ts.velocities[self.indices, :] = v
        except AttributeError:
            raise NoDataError("Timestep does not contain velocities")

    velocities = property(get_velocities, set_velocities, doc="""\
        numpy array of the velocities of the atoms in the group.

        If the trajectory does not contain velocity information then a
        :exc:`~MDAnalysis.NoDataError` is raised.

        .. versionadded:: 0.7.5
        .. deprecated:: 0.7.6
           In 0.8 this will become an attribute! You can already use :meth:`get_velocities`
           and :meth:`set_velocities`.
        .. versionchanged:: 0.8
           Became an attribute.
    """)

    def get_forces(self, ts=None, copy=False, dtype=np.float32):
        """
        Get a numpy array of the atomic forces (if available).
        Currently only supported for Gromacs .trr trajectories.

        :Keywords:
           *ts*
               If *ts* is provided then positions are read from that
               :class:`~MDAnalysis.coordinates.base.Timestep` instead of
               the one from the current trajectory belonging to this universe.
               The *ts* is indexed with the indices returned by
               :meth:`~AtomGroup.indices` and it is the user's responsibility
               to provide a time step that has the appropriate dimensions.
           *copy*
               ``True``: always make a copy (slow), ``False``: Try to
               return a array view or reference (faster); note that for
               passing coordinates to C-code it can be necessary to use
               a copy [``False``]
           *dtype*
               numpy data type of the array; the default is usually
               entirely appropriate. Most C-code actually requires the
               default  [:class:`np.float32`]

        Forces can also be directly obtained from the attribute
        :attr:`~AtomGroup.forces`.

        Forces can be directly set with :meth:`~AtomGroup.set_forces` or
        by assigning to :attr:`~AtomGroup.forces`.

        .. versionadded:: 0.7.7
        """
        if ts is None:
            ts = self.universe.trajectory.ts
        try:
            return np.array(ts.forces[self.indices], copy=copy, dtype=dtype)
        except (AttributeError, NoDataError):
            raise NoDataError("Timestep does not contain forces")

    def set_forces(self, forces, ts=None):
        """Set the forces for all atoms in the group.

        :Arguments:
           *forces*
               a Nx3 numpy :class:`numpy.ndarray` where N is the number of
               atoms in this atom group.

        :Keywords:
           *ts*
              :class:`~MDAnalysis.coordinates.base.Timestep`, defaults
              to ``None`` and then the current time step is used.

        .. Note::

           If the group contains N atoms and *force* is a single vector (i.e. an
           array of length 3) then all N atom positions are set to *force* (due
           to numpy's broadcasting rules), as described for
           :attr:`~AtomGroup.forces`.

        See also :meth:`~AtomGroup.get_forces` and attribute access through
        :attr:`~AtomGroup.forces`.

        .. versionadded:: 0.7.7
        """
        if ts is None:
            ts = self.universe.trajectory.ts
        try:
            ts.forces[self.indices, :] = forces
        except AttributeError:
            raise NoDataError("Timestep does not contain forces")

    forces = property(get_forces, set_forces,
                      doc="""
                Forces on the atoms in the AtomGroup.

                The forces can be changed by assigning an array of the appropriate
                shape, i.e. either Nx3 to assign individual force or 3, to assign
                the *same* force to all atoms (e.g. ``ag.forces = array([0,0,0])``
                will set all forces to (0.,0.,0.)).

                For more control use the :meth:`~AtomGroup.get_forces` and
                :meth:`~AtomGroup.set_forces` methods.

                .. versionadded:: 0.7.7""")


    def transform(self, M):
        r"""Apply homogenous transformation matrix *M* to the coordinates.

        The matrix *M* must be a 4x4 matrix, with the rotation in
        `R = `M[:3,:3]`` and the translation in ``t = M[:3,3]``.

        The rotation :math:`\mathsf{R}` is applied before the
        translation :math:`\mathbf{t}`:

        .. math::

           \mathbf{x}' = \mathsf{R}\mathbf{x} + \mathbf{t}

        .. SeeAlso: :mod:`MDAnalysis.core.transformations`
        """
        R = M[:3, :3]
        t = M[:3, 3]
        # changes the coordinates (in place)
        x = self.universe.trajectory.ts.positions
        idx = self.indices
        x[idx] = np.dot(x[idx], R.T)
        x[idx] += t
        return R

    def translate(self, t):
        r"""Apply translation vector *t* to the selection's coordinates.

          >>> AtomGroup.translate(t)
          >>> AtomGroup.translate((A, B))

        The method applies a translation to the AtomGroup from current
        coordinates :math:`\mathbf{x}` to new coordinates :math:`\mathbf{x}'`:

        .. math::

            \mathbf{x}' = \mathbf{x} + \mathbf{t}

        The translation can also be given as a tuple of two MDAnalysis objects
        such as two selections `(selA, selB)`, i.e. two :class:`AtomGroup`, or
        two :class:`Atom`. The translation vector is computed as the
        difference between the centers of geometry (centroid) of B and A::

            t = B.centroid() - A.centroid()
        """
        try:
            sel1, sel2 = t
            x1, x2 = sel1.centroid(), sel2.centroid()
            vector = x2 - x1
        except (ValueError, AttributeError):
            vector = np.asarray(t)
        # changes the coordinates (in place)
        self.universe.trajectory.ts.positions[self.indices] += vector
        return vector

    def rotate(self, R):
        r"""Apply a rotation matrix *R* to the selection's coordinates.

        AtomGroup.rotate(R)

        :math:`\mathsf{R}` is a 3x3 orthogonal matrix that transforms a vector
        :math:`\mathbf{x} \rightarrow \mathbf{x}'`:

        .. math::

            \mathbf{x}' = \mathsf{R}\mathbf{x}
        """
        R = np.matrix(R, copy=False, dtype=np.float32)
        # changes the coordinates (in place)
        x = self.universe.trajectory.ts.positions
        idx = self.indices
        x[idx] = x[idx] * R.T  # R.T acts to the left & is broadcasted N times.
        return R

    def rotateby(self, angle, axis, point=None):
        r"""Apply a rotation to the selection's coordinates.

        AtomGroup.rotateby(angle,axis[,point])

        The transformation from current coordinates :math:`\mathbf{x}`
        to new coordinates :math:`\mathbf{x}'` is

        .. math::

          \mathbf{x}' = \mathsf{R}\,(\mathbf{x}-\mathbf{p}) + \mathbf{p}

        where :math:`\mathsf{R}` is the rotation by *angle* around the
        *axis* going through *point* :math:`\mathbf{p}`.

        :Arguments:
          *angle*
             rotation angle in degrees
          *axis*
             rotation axis vector, a 3-tuple, list, or array, or a 2-tuple of
             two MDAnalysis objects from which the axis is calculated as the
             vector from the first to the second center of geometry.
          *point*
             point on the rotation axis; by default (``None``) the center of
             geometry of the selection is chosen, or, if *axis* is a tuple of
             selections, it defaults to the first point of the axis. *point*
             can be a 3-tuple, list, or array or a MDAnalysis object (in which
             case its :meth:`centroid` is used).

        :Returns: The 4x4 matrix which consists of the rotation matrix ``M[:3,:3]``
                  and the translation vector ``M[:3,3]``.
        """
        alpha = np.radians(angle)
        try:
            sel1, sel2 = axis
            x1, x2 = sel1.centroid(), sel2.centroid()
            v = x2 - x1
            n = v / np.linalg.norm(v)
            if point is None:
                point = x1
        except (ValueError, AttributeError):
            n = np.asarray(axis)
        if point is None:
            p = self.centroid()
        else:
            try:
                p = point.centroid()
            except AttributeError:
                p = np.asarray(point)
        M = transformations.rotation_matrix(alpha, n, point=p)
        self.transform(M)
        return M

    def align_principal_axis(self, axis, vector):
        """Align principal axis with index *axis* with *vector*.

        :Arguments:
          *axis*
            Index of the principal axis (0, 1, or 2), as produced by
            :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.principal_axes`.
          *vector*
            A 3D vector such as the z-axis (``[0,0,1]``); can be
            anything that looks like a list with three entries.

        To align the long axis of a channel (the first principal axis,
        i.e. *axis* = 0) with the z-axis::

          u.atoms.align_principal_axis(0, [0,0,1])
          u.atoms.write("aligned.pdb")
        """
        p = self.principal_axes()[axis]
        angle = np.degrees(mdamath.angle(p, vector))
        ax = transformations.rotaxis(p, vector)
        #print "principal[%d] = %r" % (axis, p)
        #print "axis = %r, angle = %f deg" % (ax, angle)
        return self.rotateby(angle, ax)

    align_principalAxis = deprecate(align_principal_axis,
                                    old_name='align_principalAxis',
                                    new_name='align_principal_axis')

    def pack_into_box(self, box=None, inplace=True):
        r"""Shift all atoms in this group to be within the primary unit cell.

        AtomGroup.packintobox([box, [inplace=True]])

        :Keywords:
          *box*
            Unit cell to move atoms inside of.
          *inplace*
            ``True``: Change coordinates in place and return
            ``False``: Only return the coordinates

        All atoms will be moved so that they lie between 0 and
        boxlength :math:`L_i` in all dimensions, i.e. the lower left
        corner of the simulation box is taken to be at (0,0,0):

        .. math::

           x_i' = x_i - \left\lfloor\frac{x_i}{L_i}\right\rfloor

        The default is to take unit cell information from the
        underlying :class:`~MDAnalysis.coordinates.base.Timestep`
        instance. The optional argument *box* can be used to provide
        alternative unit cell information (in the MDAnalysis standard
        format ``[Lx, Ly, Lz, alpha, beta, gamma]``).

        Works with either orthogonal or triclinic box types.

        By default the coordinates are changed in place and returned

        .. versionadded:: 0.8

        """
        if box is None:  #Try and auto detect box dimensions
            box = self.dimensions  # Can accept any box

        if box.shape == (3, 3):
            # for a vector representation, diagonal cannot be zero
            if (box.diagonal() == 0.0).any():
                raise ValueError("One or more box dimensions is zero."
                                 "  You can specify a boxsize with 'box ='")
        else:
            if (box == 0).any():  #Check that a box dimension isn't zero
                raise ValueError("One or more box dimensions is zero."
                                 "  You can specify a boxsize with 'box='")

        coords = self.universe.coord.positions[self.indices]
        if not inplace:
            return distances.apply_PBC(coords, box)

        self.universe.coord.positions[self.indices] = distances.apply_PBC(coords, box)

        return self.universe.coord.positions[self.indices]

    packIntoBox = deprecate(pack_into_box, old_name='packIntoBox', new_name='pack_into_box')

    def wrap(self, compound="atoms", center="com", box=None):
        """Shift the contents of this AtomGroup back into the unit cell.

        This is a more powerful version of :meth:`pack_into_box`, allowing
        groups of atoms to be kept together through the process.

        :Keywords:
           *compound*
               The group which will be kept together through the shifting
               process. [``atoms``]

               Possible options:

                   - ``atoms``
                   - ``group`` - This AtomGroup
                   - ``residues``
                   - ``segments``
                   - ``fragments``

           *center*
               How to define the center of a given group of atoms [``com``]
           *box*
               Unit cell information.  If not provided, the values from
               Timestep will be used.

        When specifying a *compound*, the translation is calculated based on
        each compound. The same translation is applied to all atoms
        within this compound, meaning it will not be broken by the shift.
        This might however mean that all atoms from the compound are not
        inside the unit cell, but rather the center of the compound is.
        Compounds available for use are *atoms*, *residues*,
        *segments* and *fragments*

        *center* allows the definition of the center of each group to be
        specified.  This can be either 'com' for center of mass, or 'cog'
        for center of geometry.

        *box* allows a unit cell to be given for the transformation.  If not
        specified, an the dimensions information from the current Timestep
        will be used.

        .. Note::
           wrap with all default keywords is identical to :meth:`pack_into_box`

        .. versionadded:: 0.9.2
        """
        if compound.lower() == "atoms":
            return self.pack_into_box(box=box)

        if compound.lower() == 'group':
            objects = [self.atoms]
        elif compound.lower() == 'residues':
            objects = self.residues
        elif compound.lower() == 'segments':
            objects = self.segments
        elif compound.lower() == 'fragments':
            objects = self.fragments
        else:
            raise ValueError("Unrecognised compound definition: {0}"
                             "Please use one of 'group' 'residues' 'segments'"
                             "or 'fragments'".format(compound))

        if center.lower() in ('com', 'centerofmass'):
            centers = np.vstack([o.center_of_mass() for o in objects])
        elif center.lower() in ('cog', 'centroid', 'centerofgeometry'):
            centers = np.vstack([o.center_of_geometry() for o in objects])
        else:
            raise ValueError("Unrecognised center definition: {0}"
                             "Please use one of 'com' or 'cog'".format(center))
        centers = centers.astype(np.float32)

        if box is None:
            box = self.dimensions

        # calculate shift per object center
        dests = distances.apply_PBC(centers, box=box)
        shifts = dests - centers

        for o, s in itertools.izip(objects, shifts):
            # Save some needless shifts
            if not all(s == 0.0):
                o.translate(s)

    def select_atoms(self, sel, *othersel, **selgroups):
        """Selection of atoms using the MDAnalysis selection syntax.

        AtomGroup.select_atoms(selection[,selection[,...]], [groupname=atomgroup[,groupname=atomgroup[,...]]])

        .. SeeAlso:: :meth:`Universe.select_atoms`
        """
        from . import Selection  # can ONLY import in method, otherwise cyclical import!

        atomgrp = Selection.Parser.parse(sel, selgroups).apply(self)
        if len(othersel) == 0:
            return atomgrp
        else:
            # Generate a selection for each selection string
            #atomselections = [atomgrp]
            for sel in othersel:
                atomgrp = atomgrp + Selection.Parser.parse(sel, selgroups).apply(self)
                #atomselections.append(Selection.Parser.parse(sel).apply(self))
            #return tuple(atomselections)
            return atomgrp

    selectAtoms = deprecate(select_atoms, old_name='selectAtoms',
                            new_name='select_atoms')

    def split(self, level):
        """Split atomgroup into a list of atomgroups by *level*.

        *level* can be "atom", "residue", "segment".
        .. versionadded:: 0.9.0
        """
        # CHECK: What happens to duplicate atoms (with advanced slicing)?

        accessors = {'segment': 'segid', 'segid': 'segid',
                     'residue': 'resid', 'resid': 'resid',
                     }

        if level == "atom":
            return [AtomGroup([a]) for a in self]

        # more complicated groupings
        try:
            # use own list comprehension to avoid sorting/compression by eg self.resids
            ids = np.array([getattr(atom, accessors[level]) for atom in self])
        except KeyError:
            raise ValueError("level = '{0}' not supported, must be one of {1}".format(
                    level, accessors.keys()))

        # now sort the resids so that it doesn't matter if the AG contains
        # atoms in random order (i.e. non-sequential resids); groupby needs
        # presorted keys!
        idx = np.argsort(ids)
        sorted_ids = ids[idx]
        # group (index, key) and then pull out the index for each group to form AtomGroups
        # by indexing self (using advanced slicing eg g[[1,2,3]]
        groups = [
            self[[idx_k[0] for idx_k in groupings]]  # one AtomGroup for each residue or segment
            for k, groupings in itertools.groupby(itertools.izip(idx, sorted_ids), lambda v: v[1])
            ]
        return groups

    def write(self, filename=None, format="PDB",
              filenamefmt="%(trjname)s_%(frame)d", **kwargs):
        """Write AtomGroup to a file.

        AtomGroup.write(filename[,format])

        :Keywords:
          *filename*
               ``None``: create TRJNAME_FRAME.FORMAT from filenamefmt [``None``]
          *format*
                PDB, CRD, GRO, VMD (tcl), PyMol (pml), Gromacs (ndx) CHARMM (str)
                Jmol (spt); case-insensitive and can also be supplied as the
                filename extension [PDB]
          *filenamefmt*
                format string for default filename; use substitution tokens
                'trjname' and 'frame' ["%(trjname)s_%(frame)d"]
          *bonds*
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
        frame = trj.ts.frame

        if trj.n_frames == 1: kwargs.setdefault("multiframe", False)

        if filename is None:
            trjname, ext = os.path.splitext(os.path.basename(trj.filename))
            filename = filenamefmt % vars()
        filename = util.filename(filename, ext=format.lower(), keep=True)

        # From the following blocks, one must pass.
        # Both can't pass as the extensions don't overlap.
        try:
            writer = MDAnalysis.coordinates.writer(filename, **kwargs)
        except TypeError:
            coords = False
            pass  # might be selections format
        else:
            coords = True

        try:
            SelectionWriter = MDAnalysis.selections.get_writer(filename, format)
        except (TypeError, NotImplementedError):
            selection = False
            pass
        else:
            writer = SelectionWriter(filename, **kwargs)
            selection = True

        if not (coords or selection):
            raise ValueError("No writer found for format: {0}".format(filename))
        else:
            writer.write(self.atoms)
            if coords:  # only these writers have a close method
                writer.close()

    # TODO: This is _almost_ the same code as write() --- should unify!
    def write_selection(self, filename=None, format="vmd", filenamefmt="%(trjname)s_%(frame)d",
                        **kwargs):
        """Write AtomGroup selection to a file to be used in another programme.

        :Keywords:
          *filename*
                ``None``: create TRJNAME_FRAME.FORMAT from *filenamefmt*
          *format*
                output file format: VMD (tcl), PyMol (pml), Gromacs (ndx), CHARMM (str);
                can also be supplied as the filename extension. Case insensitive. [vmd]
          *filenamefmt*
                format string for default filename; use '%(trjname)s' and '%(frame)s'
                placeholders; the extension is set according to the *format*
                ["%(trjname)s_%(frame)d"]
          *kwargs*
                additional keywords are passed on to the appropriate
                :class:`~MDAnalysis.selections.base.SelectionWriter`

        .. deprecated:: 0.9.0
           Use :meth:`write`
        """
        import MDAnalysis.selections

        SelectionWriter = MDAnalysis.selections.get_writer(filename, format)

        trj = self.universe.trajectory  # unified trajectory API
        frame = trj.ts.frame

        # get actual extension from the static class attribute
        extension = SelectionWriter.ext

        if filename is None:
            trjname, ext = os.path.splitext(os.path.basename(trj.filename))
            filename = filenamefmt % vars()
        filename = util.filename(filename, ext=extension, keep=True)

        writer = SelectionWriter(filename, **kwargs)
        writer.write(self.atoms)  # wants a atomgroup

    @property
    def dimensions(self):
        """Dimensions of the Universe to which the group belongs, at the current time step."""
        if self.universe is not None:
            return self.universe.dimensions
        else:
            raise AttributeError("This AtomGroup does not belong to a Universe with a dimension.")

    @dimensions.setter
    def dimensions(self, box):
        """Pass on to Universe setter

        .. versionadded:: 0.9.0
        """
        self.universe.dimensions = box

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

        if self._ts is None or self._ts.frame != trj_ts.frame:
            # create a timestep of same type as the underlying trajectory
            self._ts = trj_ts.copy_slice(self.indices)
        return self._ts


class Residue(AtomGroup):
    """A group of atoms corresponding to a residue.

    Pythonic access to atoms:
      - Using a atom name as attribute returns the matching atom (a
        :class:`Atom` instance), i.e. ``r.name``. Example::

          >>> from MDAnalysis.tests.datafiles import PSF,DCD
          >>> u = Universe(PSF,DCD)
          >>> print(u.s4AKE.r1.CA)  # C-alpha of M1
          < Atom 5: name 'CA' of type '22' of resname 'MET', resid 1 and segid '4AKE'>

      - ``r['name']`` or ``r[id]`` - returns the atom corresponding to that name

    :Data:
      :attr:`Residue.name`
        Three letter residue name.
      :attr:`Residue.id`
        Numeric (integer) resid, taken from the topology.
      :attr:`Residue.resnum`
        Numeric canonical residue id (e.g. as used in the PDB structure).

    .. Note::

       Creating a :class:`Residue` modifies the underlying :class:`Atom`
       instances. Each :class:`Atom` can only belong to a single
       :class:`Residue`.

    .. versionchanged:: 0.7.4
       Added :attr:`Residue.resnum` attribute and *resnum* keyword argument.
    """
    ## FIXME (see below, Issue 70)
    ##_cache = {}

    def __init__(self, name, id, atoms, resnum=None):
        super(Residue, self).__init__(atoms)
        self.name = name
        self.id = id
        if resnum is not None:
            self.resnum = resnum
        else:
            self.resnum = self.id  # TODO: get resnum from topologies that support it
        self.segment = None
        for i, a in enumerate(atoms):
            a.id = i
            a.resnum = self.resnum
            a.residue = self

        self._cls = AtomGroup
        # Should I cache the positions of atoms within a residue?
        # FIXME: breaks when termini are used to populate the cache; termini typically
        #        have the SAME residue name but different atoms!!! Issue 70
        ##if not Residue._cache.has_key(name):
        ##    Residue._cache[name] = dict([(a.name, i) for i, a in enumerate(self._atoms)])

    def phi_selection(self):
        """AtomGroup corresponding to the phi protein backbone dihedral C'-N-CA-C.

        :Returns: 4-atom selection in the correct order.  If no C'
                  found in the previous residue (by resid) then this
                  method returns ``None``.
        """
        sel = self.universe.select_atoms(
            'segid %s and resid %d and name C' % (self.segment.id, self.id - 1)) + \
              self['N'] + self['CA'] + self['C']
        if len(sel) == 4:  # select_atoms doesnt raise errors if nothing found, so check size
            return sel
        else:
            return None

    def psi_selection(self):
        """AtomGroup corresponding to the psi protein backbone dihedral N-CA-C-N'.

        :Returns: 4-atom selection in the correct order.  If no N'
                  found in the following residue (by resid) then this
                  method returns ``None``.
        """
        sel = self['N'] + self['CA'] + self['C'] + \
              self.universe.select_atoms(
                  'segid %s and resid %d and name N' % (self.segment.id, self.id + 1))
        if len(sel) == 4:
            return sel
        else:
            return None

    def omega_selection(self):
        """AtomGroup corresponding to the omega protein backbone dihedral CA-C-N'-CA'.

        omega describes the -C-N- peptide bond. Typically, it is trans
        (180 degrees) although cis-bonds (0 degrees) are also
        occasionally observed (especially near Proline).

        :Returns: 4-atom selection in the correct order.  If no C'
                  found in the previous residue (by resid) then this
                  method returns ``None``.

        """
        nextres = self.id + 1
        segid = self.segment.id
        sel = self['CA'] + self['C'] + \
              self.universe.select_atoms(
                  'segid %s and resid %d and name N' % (segid, nextres),
                  'segid %s and resid %d and name CA' % (segid, nextres))
        if len(sel) == 4:
            return sel
        else:
            return None

    def chi1_selection(self):
        """AtomGroup corresponding to the chi1 sidechain dihedral N-CA-CB-CG.

        :Returns: 4-atom selection in the correct order.  If no CB and/or CG is
                  found then this method returns ``None``.

        .. versionadded:: 0.7.5
        """
        try:
            return self['N'] + self['CA'] + self['CB'] + self['CG']
        except (SelectionError, NoDataError):
            return None

    def __repr__(self):
        return "<Residue {name}, {id}>".format(
            name=self.name, id=self.id)


class ResidueGroup(AtomGroup):
    """A group of residues.

    Pythonic access to atoms:
      - Using a atom name as attribute returns a list of all atoms (a
        :class:`AtomGroup`) of the same name. Example::

          >>> from MDAnalysis.tests.datafiles import PSF,DCD
          >>> u = Universe(PSF,DCD)
          >>> print(u.s4AKE.MET.CA)  # C-alpha of all Met
          <AtomGroup with 6 atoms>

    :Data: :attr:`ResidueGroup._residues`

    """
    _containername = "_residues"

    def __init__(self, residues):
        """Initialize the ResidueGroup with a list of :class:`Residue` instances."""
        self._residues = residues
        atoms = []
        for res in residues:
            atoms.extend(res.atoms)
        super(ResidueGroup, self).__init__(atoms)
        self._cls = self.__class__

    def _set_residues(self, name, value, **kwargs):
        """Set attribute *name* to *value* for all residues in the :class:`ResidueGroup`.

        If *value* is a sequence of the same length as the
        :class:`ResidueGroup` (:attr:`AtomGroup.residues`) then each
        :class:`Residue`'s property *name* is set to the corresponding
        value. If *value* is neither of length 1 (or a scalar) nor of the
        length of the :class:`ResidueGroup` then a :exc:`ValueError` is raised.

        .. versionadded:: 0.7.5
        .. versionchanged:: 0.8
           Can set residues to distinct values by providing a sequence or iterable.
        """
        values = util.asiterable(value)
        if len(values) == 1:
            self._set_atoms(name, values[0], **kwargs)
        elif len(values) == len(self.residues):
            for r, value in itertools.izip(self.residues, values):
                r._set_atoms(name, value, **kwargs)
        else:
            raise ValueError("set_residues: can only set all atoms to a single value or each atom to a distinct one "
                             "but len(residues)={0} whereas len(value)={1}".format(len(self.residues), len(values)))

        # also fix self --- otherwise users will get confused if the changes are not reflected in the
        # object they are currently using (it works automatically for AtomGroup but not higher order groups)
        #
        # This is a hack to be able to set properties on Atom and Residue
        # instances where they have different names
        attr = {'resname': 'name',
            'resid': 'id'}
        for r, value in itertools.izip(self.residues, itertools.cycle(values)):
            attrname = attr.get(name, name)
            if hasattr(r, attrname):  # should use __slots__ on Residue and try/except here
                setattr(r, attrname, value)

    set = _set_residues

    @property
    def resids(self):
        """Returns an array of residue numbers.

        .. versionchanged:: 0.8
           Returns a :class:`numpy.ndarray`
        .. versionchanged:: 0.11.0
           Now a property and returns array of length `len(self)`
        """
        return np.array([r.id for r in self.residues])

    @property
    def resnames(self):
        """Returns an array of residue names.

        .. versionchanged:: 0.8
           Returns a :class:`numpy.ndarray`
        .. versionchanged:: 0.11.0
           Now a property and returns array of length `len(self)`
        """
        return np.array([r.name for r in self.residues])

    @property
    def resnums(self):
        """Returns an array of canonical residue numbers.

        .. versionadded:: 0.7.4
        .. versionchanged:: 0.8
           Returns a :class:`numpy.ndarray`
        .. versionchanged:: 0.11.0
           Now a property and returns array of length `len(self)`
        """
        return np.array([r.resnum for r in self.residues])

    @property
    def segids(self):
        """Returns an array of segment names.

        .. note:: uses the segment of the first atom in each residue for the
                  segment name returned

        .. versionchanged:: 0.8
           Returns a :class:`numpy.ndarray`
        .. versionchanged:: 0.11.0
           Now a property and returns array of length `len(self)`
        """
        # a bit of a hack to use just
        return np.array([r[0].segid for r in self.residues])

    def set_resids(self, resid):
        """Set the resids to integer *resid* for **all residues** in the
        :class:`ResidueGroup`.

        If *resid* is a sequence of the same length as the :class:`ResidueGroup`
        then each :attr:`Atom.resid` is set to the corresponding value together
        with the :attr:`Residue.id` of the residue the atom belongs to. If
        *value* is neither of length 1 (or a scalar) nor of the length of the
        :class:`AtomGroup` then a :exc:`ValueError` is raised.

        .. Note::

           Changing resids can change the topology.

           Assigning the same *resid* to multiple residues will **merge** these
           residues. The new residue name will be the name of the first old
           residue in the merged residue.

        .. Warning::

           The values of *this* :class:`ResidueGroup` are not being
           changed. You **must create a new** :class:`ResidueGroup` **from the**
           :class:`Universe` --- only :class:`Atom` instances are changed,
           everything else is derived from these atoms.

        .. versionadded:: 0.8
        .. versionchanged:: 0.11.0
           Made plural to make consistent with corresponding property
        """
        super(ResidueGroup, self).set_resids(resid)

    set_resid = deprecate(set_resids, old_name='set_resid', new_name='set_resids')

    def set_resnums(self, resnum):
        """Set the resnums to *resnum* for **all residues** in the :class:`ResidueGroup`.

        If *resnum* is a sequence of the same length as the :class:`ResidueGroup`
        then each :attr:`Atom.resnum` is set to the corresponding value together
        with the :attr:`Residue.resnum` of the residue the atom belongs to. If
        *value* is neither of length 1 (or a scalar) nor of the length of the
        :class:`AtomGroup` then a :exc:`ValueError` is raised.

        .. Note::

           Changing *resnum* will not affect the topology: you can have
           multiple residues with the same *resnum*.

        .. SeeAlso:: :meth:`set_resid`

        .. versionadded:: 0.7.4
        .. versionchanged:: 0.7.5
           Also changes the residues.
        .. versionchanged:: 0.8
           Can set atoms and residues to distinct values by providing a sequence or iterable.
        .. versionchanged:: 0.11.0
           Made plural to make consistent with corresponding property
        """
        super(ResidueGroup, self).set_resnums(resnum)

    set_resnum = deprecate(set_resnums, old_name='set_resnum', new_name='set_resnums')

    def set_resnames(self, resname):
        """Set the resnames to string *resname* for **all residues** in the
        :class:`ResidueGroup`.

        If *resname* is a sequence of the same length as the :class:`ResidueGroup`
        then each :attr:`Atom.resname` is set to the corresponding value together
        with the :attr:`Residue.name` of the residue the atom belongs to. If
        *value* is neither of length 1 (or a scalar) nor of the length of the
        :class:`AtomGroup` then a :exc:`ValueError` is raised.

        .. versionadded:: 0.7.4
        .. versionchanged:: 0.7.5
           Also changes the residues.
        .. versionchanged:: 0.8
           Can set atoms and residues to distinct values by providing a sequence or iterable.
        .. versionchanged:: 0.11.0
           Made plural to make consistent with corresponding property
        """
        super(ResidueGroup, self).set_resnames(resname)

    set_resname = deprecate(set_resnames, old_name='set_resname', new_name='set_resnames')

    # All other AtomGroup.set_xxx() methods should just work as
    # ResidueGroup.set_xxx() because we overrode self.set(); the ones above
    # where kept separate because we can save a call to build_residues()
    # because there is no ambiguity as which residues are changed.

    def __repr__(self):
        return "<ResidueGroup {res}>".format(
            res=repr(list(self.residues)))


class Segment(ResidueGroup):
    """A group of residues corresponding to one segment of the topology.

    Pythonic access to residues:

      - The attribute rN returns the N-th residue :class:`Residue` of the
        segment (numbering starts at N=1). Example::

          >>> from MDAnalysis.tests.datafiles import PSF,DCD
          >>> u = Universe(PSF,DCD)
          >>> print(u.s4AKE.r1)
          <Residue 'MET', 1>

      - Using a residue name as attribute returns a list of all residues (a
        :class:`ResidueGroup`) of the same name. Example::

          >>> from MDAnalysis.tests.datafiles import PSF,DCD
          >>> u = Universe(PSF,DCD)
          >>> print(u.s4AKE.CYS)
          <ResidueGroup [<Residue 'CYS', 77>]>
          >>> print(u.s4AKE.MET)
          <ResidueGroup [<Residue 'MET', 1>, <Residue 'MET', 21>, <Residue 'MET', 34>, <Residue 'MET', 53>,
          <Residue 'MET', 96>, <Residue 'MET', 174>]>

    :Data: :attr:`Segment.name` is the segid from the topology or the
           chain identifier when loaded from a PDB
    """

    def __init__(self, name, residues):
        """Initialize a Segment with segid *name* from a list of :class:`Residue` instances."""
        super(Segment, self).__init__(residues)
        self.name = name
        for res in self.residues:
            res.segment = self
            for atom in res:
                atom.segment = self
        self._cls = ResidueGroup

    @property
    def id(self):
        """Segment id (alias for :attr:`Segment.name`)"""
        return self.name

    @id.setter
    def id(self, x):
        self.name = x

    def __getattr__(self, attr):
        if attr[0] == 'r':
            resnum = int(attr[1:]) - 1  # 1-based for the user, 0-based internally
            return self[resnum]
        else:
            # There can be multiple residues with the same name
            r = []
            for res in self.residues:
                if (res.name == attr):
                    r.append(res)
            if (len(r) == 0):
                return super(Segment, self).__getattr__(attr)
            # elif (len(r) == 1): return r[0]  ## creates unexpected behaviour (Issue 47)
            else:
                return ResidueGroup(r)

    def __repr__(self):
        return "<Segment {name}>".format(
            name=self.name)


class SegmentGroup(ResidueGroup):
    """A group of segments.

    Pythonic access to segments:
      - Using a segid as attribute returns the segment. Because
        of python language rule, any segid starting with a non-letter
        character is prefixed with 's', thus '4AKE' --> 's4AKE'.

        Example::

          >>> from MDAnalysis.tests.datafiles import PSF,DCD
          >>> u = Universe(PSF,DCD)
          >>> print(u.atoms.segments.s4AKE)  # segment 4AKE
          <AtomGroup with 3314 atoms>

      - Indexing the group returns the appropriate segment.

    :Data: :attr:`SegmentGroup._segments`

    """
    _containername = "_segments"

    def __init__(self, segments):
        """Initialize the SegmentGroup with a list of :class:`Segment` instances."""
        self._segments = segments
        residues = []
        for s in segments:
            residues.extend(s.residues)
        super(SegmentGroup, self).__init__(residues)
        self._cls = self.__class__

    def _set_segments(self, name, value, **kwargs):
        """Set attribute *name* to *value* for all :class:`Segment` in this :class:`AtomGroup`.

        If *value* is a sequence of the same length as the
        :class:`SegmentGroup` (:attr:`AtomGroup.residues`) then each
        :class:`Segment`'s property *name* is set to the corresponding
        value. If *value* is neither of length 1 (or a scalar) nor of the
        length of the :class:`SegmentGroup` then a :exc:`ValueError` is raised.

        .. versionadded:: 0.8
        """
        values = util.asiterable(value)
        if len(values) == 1:
            self._set_atoms(name, values[0], **kwargs)
        elif len(values) == len(self.segments):
            for s, value in itertools.izip(self.segments, values):
                s._set_atoms(name, value, **kwargs)
        else:
            raise ValueError("set_segments: can only set all atoms to a single value or each atom to a distinct one "
                             "but len(segments)={0} whereas len(value)={1}".format(len(self.segments), len(values)))

    set = _set_segments

    @property
    def segids(self):
        """Returns an array of segment names.

        .. versionchanged:: 0.8
           Returns a :class:`numpy.ndarray`
        .. versionchanged:: 0.11.0
           Now a property and returns array of length `len(self)`
        """
        return np.array([s.name for s in self.segments])

    def set_segids(self, segid):
        """Set the segids to *segid* for all atoms in the :class:`SegmentGroup`.

        If *segid* is a sequence of the same length as the :class:`SegmentGroup`
        then each :attr:`Atom.segid` is set to the corresponding value together
        with the :attr:`Segment.id` of the segment the atom belongs to. If
        *value* is neither of length 1 (or a scalar) nor of the length of the
        :class:`AtomGroup` then a :exc:`ValueError` is raised.

        .. Note::

           :meth:`set_segid` can change the topology.

           This can be used to join segments or to break groups into multiple
           disjoint segments. Note that each :class:`Atom` can only belong to a
           single :class:`Segment`.

        .. versionadded:: 0.7.4
        .. versionchanged:: 0.8
           Can set atoms and residues to distinct values by providing a sequence or iterable.
        .. versionchanged:: 0.11.0
           Made plural to make consistent with corresponding property
        """
        super(SegmentGroup, self).set_segids(segid)

    set_segid = deprecate(set_segids, old_name='set_segid', new_name='set_segids')

    def __getattr__(self, attr):
        if attr.startswith('s') and attr[1].isdigit():
            attr = attr[1:]  # sNxxx only used for python, the name is stored without s-prefix
        seglist = [segment for segment in self.segments if segment.name == attr]
        if len(seglist) == 0:
            return super(SegmentGroup, self).__getattr__(attr)
        if len(seglist) > 1:
            warnings.warn("SegmentGroup: Multiple segments with the same name {0};"
                          " a combined, NON-CONSECUTIVE "
                          "Segment is returned.".format(attr), category=SelectionWarning)
            #return Segment(sum([s.residues for s in seglist])) ### FIXME: not working yet, need __add__
            return seglist[0]
        return seglist[0]

    def __repr__(self):
        return "<SegmentGroup {segnames}>".format(
            segnames=repr(list(self.segments)))


def as_Universe(*args, **kwargs):
    """Return a universe from the input arguments.

    1. If the first argument is a universe, just return it::

         as_Universe(universe) --> universe

    2. Otherwise try to build a universe from the first or the first
       and second argument::

         as_Universe(PDB, **kwargs) --> Universe(PDB, **kwargs)
         as_Universe(PSF, DCD, **kwargs) --> Universe(PSF, DCD, **kwargs)
         as_Universe(*args, **kwargs) --> Universe(*args, **kwargs)

    :Returns: an instance of :class:`~MDAnalaysis.AtomGroup.Universe`
    """
    if len(args) == 0:
        raise TypeError("as_Universe() takes at least one argument (%d given)" % len(args))
    elif len(args) == 1 and isinstance(args[0], Universe):
        return args[0]
    return Universe(*args, **kwargs)

asUniverse = deprecate(as_Universe, old_name='asUniverse', new_name='as_Universe')

def Merge(*args):
    """Return a :class:`Universe` from two or more :class:`AtomGroup` instances.

    :class:`AtomGroup` instances can come from different Universes, or come
    directly from a :meth:`~Universe.select_atoms` call.

    It can also be used with a single :class:`AtomGroup` if the user wants to,
    for example, re-order the atoms in the Universe.

    :Arguments: One or more :class:`AtomGroup` instances.

    :Returns: an instance of :class:`~MDAnalaysis.AtomGroup.Universe`

    :Raises: :exc:`ValueError` for too few arguments or if an AtomGroup is
             empty and :exc:`TypeError` if arguments are not
             :class:`AtomGroup` instances.

    .. rubric:: Example

    In this example, protein, ligand, and solvent were externally prepared in
    three different PDB files. They are loaded into separate :class:`Universe`
    objects (where they could be further manipulated, e.g. renumbered,
    relabeled, rotated, ...) The :func:`Merge` command is used to combine all
    of them together::

       u1 = Universe("protein.pdb")
       u2 = Universe("ligand.pdb")
       u3 = Universe("solvent.pdb")
       u = Merge(u1.select_atoms("protein"), u2.atoms, u3.atoms)
       u.atoms.write("system.pdb")

    The complete system is then written out to a new PDB file.

    .. Note:: Merging does not create a full trajectory but only a single
              structure even if the input consists of one or more trajectories.

    .. versionchanged 0.9.0::
       Raises exceptions instead of assertion errors.

    """
    import MDAnalysis.topology.core

    if len(args) == 0:
        raise ValueError("Need at least one AtomGroup for merging")

    for a in args:
        if not isinstance(a, AtomGroup):
            raise TypeError(repr(a) + " is not an AtomGroup")
    for a in args:
        if len(a) == 0:
            raise ValueError("cannot merge empty AtomGroup")

    coords = np.vstack([a.coordinates() for a in args])
    trajectory = MDAnalysis.coordinates.base.Reader(None)
    ts = MDAnalysis.coordinates.base.Timestep.from_coordinates(coords)
    setattr(trajectory, "ts", ts)
    trajectory.n_frames = 1

    # create an empty Universe object
    u = Universe()
    u.trajectory = trajectory

    # create a list of Atoms, then convert it to an AtomGroup
    atoms = [copy.copy(a) for gr in args for a in gr]
    for a in atoms:
        a.universe = u

    # adjust the atom numbering
    for i, a in enumerate(atoms):
        a.index = i
        a.serial = i + 1
    u.atoms = AtomGroup(atoms)

    # move over the topology
    offset = 0
    tops = ['bonds', 'angles', 'dihedrals', 'impropers']
    idx_lists = {t:[] for t in tops}
    for ag in args:
        # create a mapping scheme for this atomgroup
        mapping = {a.index:i for i, a in enumerate(ag, start=offset)}
        offset += len(ag)

        for t in tops:
            tg = getattr(ag, t)
            # Create a topology group of only bonds that are within this ag
            # ie we don't want bonds that extend out of the atomgroup
            tg = tg.atomgroup_intersection(ag, strict=True)

            # Map them so they refer to our new indices
            new_idx = [tuple(map(lambda x:mapping[x], entry))
                       for entry in tg.to_indices()]
            idx_lists[t].extend(new_idx)

    for t in tops:
        u._topology[t] = idx_lists[t]

    # adjust the residue and segment numbering (removes any remaining references to the old universe)
    MDAnalysis.topology.core.build_residues(u.atoms)
    MDAnalysis.topology.core.build_segments(u.atoms)

    return u
