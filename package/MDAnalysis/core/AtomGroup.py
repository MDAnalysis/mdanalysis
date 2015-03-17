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

  >>> r = u.selectAtoms("resid 99").residues[0]
  >>> print(r)
  <Residue 'ALA', 99>
  >>> r.name = "UNK"
  >>> print(r)
  <Residue 'UNK', 99>
  >>> rnew = u.selectAtoms("resid 99").residues[0]
  >>> print(rnew)
  <Residue 'UNK', 99>

will typically work as expected. When working with collections such as
:class:`AtomGroup` or :class:`ResidueGroup` it is generally better to use
provided setter methods such as :meth:`AtomGroup.set_resname` or
:meth:`ResidueGroup.set_resname`.

There are two cases when it is very important to use the setters:

* changing *resid*: :meth:`AtomGroup.set_resid` and :meth:`ResidueGroup.set_resid`
* changing *segid*: :meth:`AtomGroup.set_segid` and :meth:`ResidueGroup.set_segid`

Because residues are determined by the :attr:`Atom.resid` and segments by
:attr:`Atom.segid`, the above methods take extra care to rebuild the list of
segments and residues.

.. Note::

   :meth:`AtomGroup.set_resid`, :meth:`ResidueGroup.set_resid`,
   :meth:`AtomGroup.set_segid`, :meth:`ResidueGroup.set_segid` can change the
   topology: they can split or merge residues or segments.

Splitting/merging of residues is probably not very useful because no chemical
rearrangements are carried out. Manipulating segments might be more useful in
order to add additional structure to a :class:`Universe` and provide instant
segment selectors for interactive work::

  u.selectAtoms("protein").set_segid("protein")
  u.selectAtoms("resname POPE or resname POPC").set_segid("lipids")
  u.selectAtoms("resname SOL").set_segid("water")
  u.selectAtoms("resname NA or resname CL").set_segid("ions")

  u.protein.numberOfResidues()
  water_oxygens = u.water.OW

The setter methods have the additional advantage that they can assign
lists. For instance, many MD codes number residues consecutively starting from
1. However, the original structure might be missing a few atoms at the
N-terminus. Let's say that the first residue is really residue 10. In order to
store the canonical residue IDs ("resnum") one could the use ::

  import numpy as np
  protein = u.selectAtoms("protein").residues
  protein.set_resnum(np.array(protein.resnums()) + 9)

.. TODO: correct this example when resnums has become a property that returns a np array
..  protein = u.selectAtoms("protein").residues
..  protein.set_resnum(protein.resnums + 9)

One can then use ``protein.select("resnum 42")`` to select the residue that has
the canonical residue id 42 (instead of ``resid 33``).

One can also read the resids directly from  an original PDB file::

  orig = MDAnalysis.Universe("2jln.pdb")
  protein.set_resnum(orig.selectAtoms("protein").resids())


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
    u = MDAnalysis.Merge(u1.selectAtoms("protein"), u2.atoms, u3.atoms)
    u.atoms.write("system.pdb")

The complete system is then written out to a new PDB file.

Replicating
~~~~~~~~~~~

It is also possible to replicate a molecule to build a system with
multiple copies of the same molecule. In the example, we replicate an
AdK molecule and then translate and rotate the second copy::

    import MDAnalysis; from MDAnalysis.tests.datafiles import *
    u = MDAnalysis.Universe(PSF, DCD)
    p = u.selectAtoms("protein")
    m = MDAnalysis.Merge(p,p)

    # now renumber resids and segids for each copy

    # first copy of the protein (need to use atom indices because currently that's the only reliable property in the
    merged universe)
    p1 = m.selectAtoms("bynum 1:3341")
    # second copy
    p2 = m.selectAtoms("bynum 3342:6682")

    p1.set_segid("A")
    p2.set_segid("B")
    p2.residues.set_resid(p2.residues.resids() + p1.residues.resids()[-1])  # increment resids for p2 with the last
    resid from p1

    # you must regenerate the selections after modifying them (see notes in the docs!)
    # because the changed resids are not reflected in the selection (due to how residues are referenced internally)
    p1 = m.selectAtoms("segid A")       # or as before:  m.selectAtoms("bynum 1:3341")
    p2 = m.selectAtoms("segid B")

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

   .. attribute::     number

      atom number

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

      List of :class:`~MDAnalysis.topogology.core.Bond` instances that
      contains all the bonds that this atom participates in.

      .. versionadded:: 0.8.1

   .. attribute:: angles

      List of :class:`~MDAnalysis.topogology.core.Angle` instances
      that contains all the angles that this atom participates in.

      .. versionadded:: 0.8.1

   .. attribute:: torsions

      List of :class:`~MDAnalysis.topogology.core.Torsion` instances
      that contains all the dihedral torsions that this atom
      participates in.

      .. versionadded:: 0.8.1

.. autoclass:: Residue
   :members:
.. autoclass:: ResidueGroup
   :members:
.. autoclass:: Segment
   :members:
.. autoclass:: SegmentGroup
   :members:

.. autofunction:: asUniverse
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
import numpy
from numpy.linalg import eig
import itertools
from collections import defaultdict
import copy
import logging
import os.path

# Local imports
import MDAnalysis
from .. import SelectionError, NoDataError, SelectionWarning
from . import util
from .util import cached
from .transformations import rotation_matrix, vecangle, rotaxis


logger = logging.getLogger("MDAnalysis.core.AtomGroup")


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

    .. versionchanged 0.9.0
       Added fragment managed property.
       Changed bonds angles torsions impropers to be a managed property
    """

    __slots__ = (
        "number", "id", "name", "type", "resname", "resid", "segid",
        "mass", "charge", "residue", "segment",
        "__universe",
        "radius", "bfactor", "resnum", "serial", "altLoc")

    def __init__(self, number, name, type, resname, resid, segid, mass, charge,
                 residue=None, segment=None, radius=None, bfactor=None,
                 resnum=None, serial=None, altLoc=None):
        self.number = number
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
        self.__universe = None

    def __repr__(self):
        return ("<Atom {idx}: {name} of type {t} of resname {rname}, "
                "resid {rid} and segid {sid}{altloc}>".format(
                    idx=self.number + 1, name=self.name, t=self.type,
                    rname=self.resname, rid=self.resid, sid=self.segid,
                    altloc="" if not self.altLoc
                    else " and altloc {}".format(self.altLoc)))

    def __cmp__(self, other):
        return cmp(self.number, other.number)

    def __eq__(self, other):
        return self.number == other.number

    def __hash__(self):
        return hash(self.number)

    def __add__(self, other):
        if not (isinstance(other, Atom) or isinstance(other, AtomGroup)):
            raise TypeError('Can only concatenate Atoms (not "{}")'
                            ' to AtomGroup'.format(other.__class__.__name__))
        if isinstance(other, Atom):
            return AtomGroup([self, other])
        else:
            return AtomGroup([self] + other._atoms)

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

        :Returns: a numpy 1x3 array
        """
        return self.universe.coord[self.number]  # internal numbering starts at 0

    @position.setter
    def position(self, coords):
        """
        Set the current cartesian coordinates of the atom.
        @param coords: a 1x3 numpy array of {x,y,z} coordinates, or optionally
            a single scalar if you should want to set all coordinates to the same value.
        """
        self.universe.coord._pos[self.number, :] = coords  # internal numbering starts at 0

    @property
    def velocity(self):
        """Current velocity of the atom.

        A :exc:`~MDAnalysis.NoDataError` is raised if the trajectory
        does not contain velocities.

        .. versionadded:: 0.7.5
        """
        try:
            return numpy.array(self.universe.coord._velocities[self.number],
                               dtype=numpy.float32)
        except AttributeError:
            raise NoDataError("Timestep does not contain velocities")

    def centroid(self):
        """The centroid of an atom is its position, :attr:`Atom.pos`."""
        # centroid exists for compatibility with AtomGroup
        return self.pos

    @property
    def universe(self):
        """a pointer back to the Universe"""
        if not self.__universe is None:
            return self.__universe
        else:
            raise AttributeError(
                "Atom {} is not assigned to a Universe".format(self.number))

    @universe.setter
    def universe(self, universe):
        self.__universe = universe

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
        """A list of the bonds that this Atom is in

        .. versionchanged:: 0.9.0
           Changed to managed property
        """
        return self.universe._bondDict[self]

    @property
    def angles(self):
        """A list of the angles that this Atom is in

        .. versionchanged:: 0.9.0
           Changed to managed property
        """
        return self.universe._angleDict[self]

    @property
    def torsions(self):
        """A list of the torsions/dihedrals that this Atom is in

        .. versionchanged:: 0.9.0
           Changed to managed property
        """
        return self.universe._torsionDict[self]

    dihedrals = torsions

    @property
    def impropers(self):
        """A list of the improper torsions that this Atom is in

        .. versionchanged:: 0.9.0
           Changed to managed property
        """
        return self.universe._improperDict[self]


class AtomGroup(object):
    """A group of atoms.

      ag = universe.selectAtoms(atom-list)

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

    .. rubric:: References for analysis methods

    .. [Dima2004] Dima, R. I., & Thirumalai, D. (2004). Asymmetry in the
                  shapes of folded and denatured states of proteins. *J
                  Phys Chem B*, 108(21),
                  6564-6570. doi:`10.1021/jp037128y`_

    .. _10.1021/jp037128y: http://dx.doi.org/10.1021/jp037128y


    .. versionchanged:: 0.7.6
       An empty AtomGroup can be created and no longer raises a
       :exc:`NoDataError`.
    .. versionchanged:: 0.9.0
       The size at which cache is used for atom lookup is now stored as variable
       _atomcache_size within the class.
       Added fragments manged property. Is a lazily built, cached entry, similar to residues.
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
        self.__ts = None

        # caches:
        # - built on the fly when they are needed
        # - delete entry to invalidate
        self._cache = dict()

        # for generalized __getitem__ (override _containername for ResidueGroup and SegmentGroup)
        self._container = getattr(self, self._containername)
        # Define the Class that gets returned by getitem
        # Override this where return Class differs from Self (ie slicing Residue)
        self._cls = self.__class__

    def _rebuild_caches(self):
        """Rebuild all AtomGroup caches.

        A number of lists and attributes are cached. These caches are lazily
        built the first time they are needed. When editing the topology it
        might happen that not all caches were synced properly (even though that
        this is supposed to happen eventually). In this case the user can
        manually force a complete cache rebuild.

        Currently the following caches are used:

        * atoms (for "in" lookup); cache is only built for large systems with
          > 10,000 atoms
        * indices (:meth:`AtomGroup.indices`)
        * masses (:meth:`AtomGroup.masses`)
        * residues (:attr:`AtomGroup.residues`)
        * segments (:attr:`AtomGroup.segments`)
        * bonds (:attr:`AtomGroup.bonds`)
        * angles (:attr:`AtomGroup.angles`)
        * torsions (:attr:`AtomGroup.torsions`)
        * improper torsions (:attr:`AtomGroup.impropers`)

        .. SeeAlso:: :meth:`_clear_caches`

        .. versionadded:: 0.7.5
        .. versionchanged:: 0.9.0
           Added bonds/angles/torsions/impropers to rebuild.
           Reworked how things are rebuilt to avoid code duplication.
        """
        # If the number of atoms is very large, create a dictionary cache for lookup
        if len(self._atoms) > self._atomcache_size:
            self._cache['atoms'] = dict(((x, None) for x in self.__atoms))

        # Delete preexisting cache if exists
        for att in \
                [
                    'indices', 'residues', 'segments', 'masses',
                    'bonds', 'angles', 'torsions', 'impropers'
                ]:
            try:
                del self._cache[att]
            except KeyError:
                pass
        # Call each in turn to force them to build into cache
        # indices
        self._cache['indices'] = self.indices()
        # residue instances
        self._cache['residues'] = self.residues
        # segment instances
        self._cache['segments'] = self.segments
        # masses
        self._cache['masses'] = self.masses()
        # bonds angles torsions impropers
        self._cache['bonds'] = self.bonds
        self._cache['angles'] = self.angles
        self._cache['torsions'] = self.torsions
        self._cache['impropers'] = self.impropers

    def _clear_caches(self, *args):
        """Clear cache for all *args*.

        If no args are provided, all caches are cleared.

        .. SeeAlso:: :meth:`_rebuild_caches`

        .. versionadded:: 0.8
        """
        if len(args) == 0:
            self._cache = {}
        else:
            for name in args:
                try:
                    del self._cache[name]
                except KeyError:
                    pass

    def _fill_cache(self, name, value):
        """Populate _cache[name] with value.

        .. versionadded:: 0.8
        """
        self._cache[name] = value

    # AtomGroup.atoms is guaranteed to be a AtomGroup, too; keeps a consistent API
    # between AtomGroup, Residue, ResidueGroup, Segment; access the list as
    # _atoms (although atoms supports all list-like operations, too).
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
        except AttributeError:
            return None

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
        """
        container = self._container
        cls = self._cls

        # consistent with the way list indexing/slicing behaves:
        if numpy.dtype(type(item)) == numpy.dtype(int):
            return container[item]
        elif type(item) == slice:
            return cls(container[item])
        elif isinstance(item, (numpy.ndarray, list)):
            # advanced slicing, requires array or list
            return cls([container[i] for i in item])
        elif type(item) == str:
            return getattr(self, item)
        else:
            raise TypeError("Cannot slice with type: {0}".format(type(item)))

    def __getattr__(self, name):
        # There can be more than one atom with the same name
        atomlist = [atom for atom in self._atoms if name == atom.name]
        if len(atomlist) == 0:
            raise SelectionError("No atoms or attributes with name " + name)
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
        if not (isinstance(other, Atom) or isinstance(other, AtomGroup)):
            raise TypeError('Can only concatenate AtomGroup (not "{}") to'
                            ' AtomGroup'.format(other.__class__.__name__))
        if isinstance(other, AtomGroup):
            return AtomGroup(self._atoms + other._atoms)
        else:
            return AtomGroup(self._atoms + [other])

    def __repr__(self):
        return "<AtomGroup with {natoms} atoms>".format(
            natoms=len(self))

    def __getstate__(self):
        raise NotImplementedError

    def __setstate__(self, state):
        raise NotImplementedError

    def numberOfAtoms(self):
        """Total number of atoms in the group"""
        return len(self._atoms)

    def numberOfResidues(self):
        """Total number of residues in the group"""
        return len(self.residues)

    def numberOfSegments(self):
        """Total number of segments in the group"""
        return len(self.segments)

    @cached('indices')
    def indices(self):
        """Array of all :attr:`Atom.number` in the group.

        These indices are 0-based and can be used to directly index
        :attr:`Universe.atoms` or the coordinate array
        :attr:`MDAnalysis.coordinates.base.Timestep._pos`.
        """
        return numpy.array([atom.number for atom in self._atoms])

    def names(self):
        """Returns a list of atom names.

        .. versionchanged:: 0.8
           Returns a :class:`numpy.ndarray`
        """
        return numpy.array([a.name for a in self._atoms])

    def types(self):
        """Returns an array of atom types.

        .. versionadded 0.9.0
        """
        return numpy.array([a.type for a in self._atoms])

    @property
    @cached('fragments')
    def fragments(self):
        """Read-only list of fragments.

        Contains all fragments that any Atom in this AtomGroup is part of, the contents of
        the fragments may extend beyond the contents of this AtomGroup.

        .. versionadded 0.9.0
        """
        return tuple(set(a.fragment for a in self._atoms))

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

    def resids(self):
        """Returns a list of residue numbers.

        .. versionchanged:: 0.8
           Returns a :class:`numpy.ndarray`
        """
        return numpy.array([r.id for r in self.residues])

    def resnames(self):
        """Returns a list of residue names.

        .. versionchanged:: 0.8
           Returns a :class:`numpy.ndarray`
        """
        return numpy.array([r.name for r in self.residues])

    def resnums(self):
        """Returns a list of canonical residue numbers.

        .. versionadded:: 0.7.4
        .. versionchanged:: 0.8
           Returns a :class:`numpy.ndarray`
        """
        return numpy.array([r.resnum for r in self.residues])

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

        If the keyword *format* is set to ``'Seq'``, all *kwargs* are ignored and
        a :class:`Bio.Seq.Seq` instance is returned. The difference to the
        record is that the record also contains metadata and can be directly
        used as an input for other functions in :mod:`Bio`.

        If the keyword *format* is set to ``'string'``, all *kwargs* are ignored
        and a Python string is returned.

        .. rubric:: Example: Write FASTA file

        Use :func:`Bio.SeqIO.write`, which takes sequence records::

           import Bio.SeqIO

           # get the sequence record of a protein component of a Universe
           protein = u.selectAtoms("protein")
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
            sequence = "".join([util.convert_aa_code(r) for r in self.resnames()])
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

    def segids(self):
        """Returns a list of segment ids (=segment names).

        .. versionchanged:: 0.8
           Returns a :class:`numpy.ndarray`
        """
        return numpy.array([s.name for s in self.segments])

    @property
    @cached('bonds')
    def bonds(self):
        """All the bonds in this AtomGroup

        Note that these bonds might extend out of the AtomGroup, to select
        only bonds which are entirely contained by the AtomGroup use
        u.bonds.atomgroup_intersection(ag, strict=True)

        .. versionadded 0.9.0
        """
        from MDAnalysis.topology.core import TopologyGroup

        mybonds = [b for a in self._atoms for b in a.bonds]
        if len(mybonds) == 0:
            return None
        else:
            return TopologyGroup(mybonds)

    @property
    @cached('angles')
    def angles(self):
        """All the angles in this AtomGroup

        Note that these angles might extend out of the AtomGroup, to select
        only angles which are entirely contained by the AtomGroup use
        u.angles.atomgroup_intersection(ag, strict=True)

        .. versionadded 0.9.0
        """
        from MDAnalysis.topology.core import TopologyGroup

        mybonds = [b for a in self._atoms for b in a.angles]
        if len(mybonds) == 0:
            return None
        else:
            return TopologyGroup(mybonds)

    @property
    @cached('torsions')
    def torsions(self):
        """All the torsions in this AtomGroup

        Note that these torsions might extend out of the AtomGroup, to select
        only torsions which are entirely contained by the AtomGroup use
        u.torsions.atomgroup_intersection(ag, strict=True)

        .. versionadded 0.9.0
        """
        from MDAnalysis.topology.core import TopologyGroup

        mybonds = [b for a in self._atoms for b in a.torsions]
        if len(mybonds) == 0:
            return None
        else:
            return TopologyGroup(mybonds)

    @property
    @cached('impropers')
    def impropers(self):
        """All the improper torsions in this AtomGroup

        Note that these improper torsions might extend out of the AtomGroup,
        to select only torsions which are entirely contained by the AtomGroup use
        u.impropers.atomgroup_intersection(ag, strict=True)

        .. versionadded 0.9.0
        """
        from MDAnalysis.topology.core import TopologyGroup

        mybonds = [b for a in self._atoms for b in a.impropers]
        if len(mybonds) == 0:
            return None
        else:
            return TopologyGroup(mybonds)

    @cached('masses')
    def masses(self):
        """Array of atomic masses (as defined in the topology)"""
        return numpy.array([atom.mass for atom in self._atoms])

    def totalMass(self):
        """Total mass of the selection (masses are taken from the topology or guessed)."""
        return numpy.sum(self.masses(), axis=0)

    def charges(self):
        """Array of partial charges of the atoms (as defined in the topology)"""
        return numpy.array([atom.charge for atom in self._atoms])

    def totalCharge(self):
        """Sum of all partial charges (must be defined in topology)."""
        return numpy.sum(self.charges(), axis=0)

    def radii(self):
        """Array of atomic radii (as defined in the PQR file)"""
        return numpy.array([atom.radius for atom in self._atoms])

    @property
    def bfactors(self):
        """Crystallographic B-factors (from PDB) in A**2.
        """
        return numpy.array([atom.bfactor for atom in self._atoms])

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

    def set_name(self, name):
        """Set the atom name to string for *all atoms* in the AtomGroup.

        If *value* is a sequence of the same length as the :class:`AtomGroup`
        then each :attr:`Atom.name` is set to the corresponding value. If
        *value* is neither of length 1 (or a scalar) nor of the length of the
        :class:`AtomGroup` then a :exc:`ValueError` is raised.

        .. versionadded:: 0.7.4
        .. versionchanged:: 0.8
           Can set atoms to distinct values by providing a sequence or iterable.
        """
        self.set("name", name, conversion=str)

    def set_resid(self, resid):
        """Set the resid to integer *resid* for **all atoms** in the :class:`AtomGroup`.

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

        """
        from MDAnalysis.topology.core import build_residues

        self.set("resid", resid, conversion=int)
        # Note that this also automagically updates THIS AtomGroup;
        # the side effect of build_residues(self.atoms) is to update all Atoms!!!!
        self._fill_cache('residues', build_residues(self.atoms))
        # make sure to update the whole universe: the Atoms are shared but ResidueGroups are not
        if self.atoms is not self.universe.atoms:
            self.universe.atoms._fill_cache('residues', build_residues(self.universe.atoms))

    def set_resnum(self, resnum):
        """Set the resnum to *resnum* for **all atoms** in the :class:`AtomGroup`.

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
        """
        from MDAnalysis.topology.core import build_residues

        self.set("resnum", resnum)
        # build_residues() is the easiest (though expensive) way to make sure
        # that all residues get their new resnum. There's no easy way to parse
        # a per-atom resnum list (with potential duplicates) into a list of
        # corresponding residues.
        #
        # (This comment applies to the analogous methods below as well.)
        build_residues(self.atoms)

    def set_resname(self, resname):
        """Set the resname to string *resname* for **all atoms** in the :class:`AtomGroup`.

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
        """
        from MDAnalysis.topology.core import build_residues

        self.set("resname", resname, conversion=str)
        build_residues(self.atoms)

    def set_segid(self, segid, buildsegments=True):
        """Set the segid to *segid* for all atoms in the :class:`AtomGroup`.

        If *segid* is a sequence of the same length as the :class:`AtomGroup`
        then each :attr:`Atom.segid` is set to the corresponding value together
        with the :attr:`Segment.id` of the residue the atom belongs to. If
        *value* is neither of length 1 (or a scalar) nor of the length of the
        :class:`AtomGroup` then a :exc:`ValueError` is raised.

        .. Note::

           :meth:`set_segid` can change the topology.

           With the default *buildsegments* = ``True`` it can be used to join
           segments or to break groups into multiple disjoint segments. Note
           that each :class:`Atom` can only belong to a single
           :class:`Segment`.

        For performance reasons, *buildsegments* can be set to ``False``. Then
        one needs to run :meth:`Universe._build_segments` manually later in
        order to update the list of :class:`Segment` instances and regenerate
        the segid instant selectors.

        .. versionadded:: 0.7.4
        .. versionchanged:: 0.8
           Can set atoms and residues to distinct values by providing a sequence or iterable.
        """
        self.set("segid", segid)
        self._clear_caches('segments')
        if buildsegments:
            self.universe._build_segments()
        else:
            # do not even update the local segment
            pass
        self.universe.atoms._clear_caches('segments')

    def set_mass(self, mass):
        """Set the atom mass to float *mass* for **all atoms** in the AtomGroup.

        .. versionadded:: 0.7.4
        .. versionchanged:: 0.8
           Can set atoms and residues to distinct values by providing a sequence or iterable.
        """
        self.set("mass", mass, conversion=float, cache="masses")

    def set_type(self, atype):
        """Set the atom type to *atype* for **all atoms** in the AtomGroup.

        .. versionadded:: 0.7.4
        .. versionchanged:: 0.8
           Can set atoms and residues to distinct values by providing a sequence or iterable.
        """
        self.set("type", atype)

    def set_charge(self, charge):
        """Set the partial charge to float *charge* for **all atoms** in the AtomGroup.

        .. versionadded:: 0.7.4
        .. versionchanged:: 0.8
           Can set atoms and residues to distinct values by providing a sequence or iterable.
        """
        self.set("charge", charge, conversion=float)

    def set_radius(self, radius):
        """Set the atom radius to float *radius* for **all atoms** in the AtomGroup.

        .. versionadded:: 0.7.4
        .. versionchanged:: 0.8
           Can set atoms and residues to distinct values by providing a sequence or iterable.
        """
        self.set("radius", radius, conversion=float)

    def set_bfactor(self, bfactor):
        """Set the atom bfactor to float *bfactor* for **all atoms** in the AtomGroup.

        .. versionadded:: 0.7.4
        .. versionchanged:: 0.8
           Can set atoms and residues to distinct values by providing a sequence or iterable.
        """
        self.set("bfactor", bfactor, conversion=float)

    def centerOfGeometry(self, **kwargs):
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
            return numpy.sum(self.packIntoBox(inplace=False), axis=0) / self.numberOfAtoms()
        else:
            return numpy.sum(self.coordinates(), axis=0) / self.numberOfAtoms()

    centroid = centerOfGeometry

    def centerOfMass(self, **kwargs):
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
            return numpy.sum(self.packIntoBox(inplace=False) * self.masses()[:, numpy.newaxis],
                             axis=0) / self.totalMass()
        else:
            return numpy.sum(self.coordinates() * self.masses()[:, numpy.newaxis], axis=0) / self.totalMass()

    def radiusOfGyration(self, **kwargs):
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
        masses = self.masses()
        if pbc:
            recenteredpos = self.packIntoBox(inplace=False) - self.centerOfMass(pbc=True)
        else:
            recenteredpos = self.coordinates() - self.centerOfMass(pbc=False)
        rog_sq = numpy.sum(masses * numpy.sum(numpy.power(recenteredpos, 2), axis=1)) / self.totalMass()
        return numpy.sqrt(rog_sq)

    def shapeParameter(self, **kwargs):
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
        masses = self.masses()
        if pbc:
            recenteredpos = self.packIntoBox(inplace=False) - self.centerOfMass(pbc=True)
        else:
            recenteredpos = self.coordinates() - self.centerOfMass(pbc=False)
        tensor = numpy.zeros((3, 3))
        for x in xrange(recenteredpos.shape[0]):
            tensor += masses[x] * numpy.outer(recenteredpos[x, :],
                                              recenteredpos[x, :])
        tensor /= self.totalMass()
        eig_vals = numpy.linalg.eigvalsh(tensor)
        shape = 27.0 * numpy.prod(eig_vals - numpy.mean(eig_vals)) / numpy.power(numpy.sum(eig_vals), 3)
        return shape

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
        masses = self.masses()
        if pbc:
            recenteredpos = self.packIntoBox(inplace=False) - self.centerOfMass(pbc=True)
        else:
            recenteredpos = self.coordinates() - self.centerOfMass(pbc=False)
        tensor = numpy.zeros((3, 3))
        for x in xrange(recenteredpos.shape[0]):
            tensor += masses[x] * numpy.outer(recenteredpos[x, :],
                                              recenteredpos[x, :])
        tensor /= self.totalMass()
        eig_vals = numpy.linalg.eigvalsh(tensor)
        shape = (3.0 / 2.0) * numpy.sum(numpy.power(eig_vals - numpy.mean(eig_vals), 2)) / numpy.power(
            numpy.sum(eig_vals), 2)
        return shape

    def momentOfInertia(self, **kwargs):
        """Tensor of inertia as 3x3 NumPy array.

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
            pos = self.packIntoBox(inplace=False) - self.centerOfMass(pbc=True)
        else:
            pos = self.coordinates() - self.centerOfMass(pbc=False)

        masses = self.masses()
        # Create the inertia tensor
        # m_i = mass of atom i
        # (x_i, y_i, z_i) = pos of atom i
        # Ixx = sum(m_i*(y_i^2+z_i^2));
        # Iyy = sum(m_i*(x_i^2+z_i^2));
        # Izz = sum(m_i*(x_i^2+y_i^2))
        # Ixy = Iyx = -1*sum(m_i*x_i*y_i)
        # Ixz = Izx = -1*sum(m_i*x_i*z_i)
        # Iyz = Izy = -1*sum(m_i*y_i*z_i)
        tens = numpy.zeros((3,3), dtype=numpy.float64)
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
            x = self.packIntoBox(inplace=False)
        else:
            x = self.coordinates()
        return numpy.array([x.min(axis=0), x.max(axis=0)])

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
            x = self.packIntoBox(inplace=False)
            centroid = self.centerOfGeometry(pbc=True)
        else:
            x = self.coordinates()
            centroid = self.centerOfGeometry(pbc=False)
        R = numpy.sqrt(numpy.max(numpy.sum(numpy.square(x - centroid), axis=1)))
        return R, centroid

    def bond(self, pbc=False):
        """Returns the distance between atoms in a 2-atom group.

        Distance between atoms 0 and 1::

            0---1

        .. Note::

           Only makes sense for a :class:`AtomGroup` with exactly 2
           :class:`Atom`; anything else will raise a
           :exc:`ValueError`.

        :Keywords:
          *pbc*
            ``True``: Account for minimum image convention when calculating [``False``]

        .. versionadded:: 0.7.3
        .. versionchanged:: 0.8 Added *pbc* keyword
        """
        if len(self) != 2:
            raise ValueError("distance computation only makes sense for a group with exactly 2 atoms")
        if not pbc:
            return numpy.linalg.norm(self[0].pos - self[1].pos)
        else:
            return MDAnalysis.core.distances.self_distance_array(self.coordinates(), box=self.dimensions)[0]

    def angle(self):
        """Returns the angle in degrees between atoms 0, 1, 2.

        Angle between atoms 0 and 2 with apex at 1::

              2
             /
            /
           1------0

        .. Note::

           Only makes sense for a :class:`AtomGroup` with exactly 3
           :class:`Atom`; anything else will raise a :exc:`ValueError`.

        .. versionadded:: 0.7.3
        """
        if len(self) != 3:
            raise ValueError("angle computation only makes sense for a group with exactly 3 atoms")
        a = self[0].pos - self[1].pos
        b = self[2].pos - self[1].pos
        return numpy.rad2deg(numpy.arccos(numpy.dot(a, b) /
                                          (numpy.linalg.norm(a) * numpy.linalg.norm(b))))

    def improper(self):
        """Returns the improper dihedral between 4 atoms.

        The improper dihedral is calculated in the same way as the proper
        :meth:`dihedral`: The angle between the planes formed by atoms (0,1,2)
        and (1,2,3).

        .. Note::

           Only makes sense for a :class:`AtomGroup` with exactly 4
           :class:`Atom`; anything else will raise a :exc:`ValueError`. The
           interpretation of the angle as an "improper" solely depends on the
           selection of atoms and thus the user input!

        .. versionadded:: 0.7.3
        """
        return self.dihedral()

    def dihedral(self):
        """Calculate the dihedral angle in degrees.

        Dihedral angle around axis connecting atoms 1 and 2 (i.e. the angle
        between the planes spanned by atoms (0,1,2) and (1,2,3))::

                  3
                  |
            1-----2
           /
          0

        .. Note::

           Only makes sense for a :class:`AtomGroup` with exactly 4
           :class:`Atom`; anything else will raise a :exc:`ValueError`.

        .. versionadded:: 0.7.0
        """
        if len(self) != 4:
            raise ValueError("dihedral computation only makes sense for a group with exactly 4 atoms")

        A, B, C, D = self.coordinates()[:4]
        ab = A - B
        bc = B - C
        cd = C - D
        return numpy.rad2deg(util.dihedral(ab, bc, cd))

    def principalAxes(self, **kwargs):
        """Calculate the principal axes from the moment of inertia.

        e1,e2,e3 = AtomGroup.principalAxes()

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
            eigenval, eigenvec = eig(self.momentOfInertia(pbc=True))
        else:
            eigenval, eigenvec = eig(self.momentOfInertia(pbc=False))
        # Sort
        indices = numpy.argsort(eigenval)
        # Return transposed in more logical form. See Issue 33.
        return eigenvec[:, indices].T

    def get_positions(self, ts=None, copy=False, dtype=numpy.float32):
        """Get a NumPy array of the coordinates.

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
               NumPy Data type of the array; the default is usually
               entirely appropriate. Most C-code actually requires the
               default  [:class:`numpy.float32`]

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
        return numpy.array(ts[self.indices()], copy=copy, dtype=dtype)

    coordinates = get_positions
    """NumPy array of the coordinates.

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
               a Nx3 NumPy :class:`numpy.ndarray` where N is the number of
               atoms in this atom group.

        :Keywords:
           *ts*
              :class:`~MDAnalysis.coordinates.base.Timestep`, defaults
              to ``None`` and then the current time step is used.

        .. Note::

           If the group contains N atoms and *coord* is a single point (i.e. an
           array of length 3) then all N atom positions are set to *coord* (due
           to NumPy's broadcasting rules), as described for
           :attr:`~AtomGroup.positions`.

        See also :meth:`~AtomGroup.get_positions` and attribute access through
        :attr:`~AtomGroup.positions`.

        .. versionadded:: 0.7.6
        """
        if ts is None:
            ts = self.universe.trajectory.ts
        if hasattr(ts, 'has_x'):  # TRR handling must be told frame now holds valid coord info.
            ts.has_x = True
        ts._pos[self.indices(), :] = coords

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

    def get_velocities(self, ts=None, copy=False, dtype=numpy.float32):
        """NumPy array of the velocities.

        Raises a :exc:`NoDataError` if the underlying
        :class:`~MDAnalysis.coordinates.base.Timestep` does not contain
        :attr:`~MDAnalysis.coordinates.base.Timestep._velocities`.

        See also :meth:`AtomGroup.set_velocities` and attribute access through
        :attr:`AtomGroup.velocities`.

        .. versionadded:: 0.7.6
        """
        if ts is None:
            ts = self.universe.trajectory.ts
        try:
            return numpy.array(ts._velocities[self.indices()], copy=copy, dtype=dtype)
        except AttributeError:
            raise NoDataError("Timestep does not contain velocities")

    def set_velocities(self, v, ts=None):
        """Assign the velocities *v* to the timestep.

        Raises a :exc:`NoDataError` if the underlying
        :class:`~MDAnalysis.coordinates.base.Timestep` does not contain
        :attr:`~MDAnalysis.coordinates.base.Timestep._velocities`.

        See also :meth:`AtomGroup.get_velocities` and :attr:`AtomGroup.velocities` for
        attribute access.

        .. versionadded:: 0.7.6
        """
        if ts is None:
            ts = self.universe.trajectory.ts
        try:
            if hasattr(ts, 'has_v'):  # TRR handling must be told frame now holds valid velocity info.
                ts.has_v = True
            ts._velocities[self.indices(), :] = v
        except AttributeError:
            raise NoDataError("Timestep does not contain velocities")

    velocities = property(get_velocities, set_velocities, doc="""\
        NumPy array of the velocities of the atoms in the group.

        If the trajectory does not contain velocity information then a
        :exc:`~MDAnalysis.NoDataError` is raised.

        .. versionadded:: 0.7.5
        .. deprecated:: 0.7.6
           In 0.8 this will become an attribute! You can already use :meth:`get_velocities`
           and :meth:`set_velocities`.
        .. versionchanged:: 0.8
           Became an attribute.
    """)

    def get_forces(self, ts=None, copy=False, dtype=numpy.float32):
        """
        Get a NumPy array of the atomic forces (if available).
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
               NumPy Data type of the array; the default is usually
               entirely appropriate. Most C-code actually requires the
               default  [:class:`numpy.float32`]

        Forces can also be directly obtained from the attribute
        :attr:`~AtomGroup.forces`.

        Forces can be directly set with :meth:`~AtomGroup.set_forces` or
        by assigning to :attr:`~AtomGroup.forces`.

        .. versionadded:: 0.7.7
        """
        if ts is None:
            ts = self.universe.trajectory.ts
        try:
            return numpy.array(ts._forces[self.indices()], copy=copy, dtype=dtype)
        except AttributeError:
            raise NoDataError("Timestep does not contain forces")

    def set_forces(self, forces, ts=None):
        """Set the forces for all atoms in the group.

        :Arguments:
           *forces*
               a Nx3 NumPy :class:`numpy.ndarray` where N is the number of
               atoms in this atom group.

        :Keywords:
           *ts*
              :class:`~MDAnalysis.coordinates.base.Timestep`, defaults
              to ``None`` and then the current time step is used.

        .. Note::

           If the group contains N atoms and *force* is a single vector (i.e. an
           array of length 3) then all N atom positions are set to *force* (due
           to NumPy's broadcasting rules), as described for
           :attr:`~AtomGroup.forces`.

        See also :meth:`~AtomGroup.get_forces` and attribute access through
        :attr:`~AtomGroup.forces`.

        .. versionadded:: 0.7.7
        """
        if ts is None:
            ts = self.universe.trajectory.ts
        try:
            if hasattr(ts, 'has_f'):  # TRR handling must be told frame now holds valid force info.
                ts.has_f = True
            ts._forces[self.indices(), :] = forces
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
        x = self.universe.trajectory.ts._pos
        idx = self.indices()
        x[idx] = numpy.dot(x[idx], R.T)
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
            vector = numpy.asarray(t)
        # changes the coordinates (in place)
        self.universe.trajectory.ts._pos[self.indices()] += vector
        return vector

    def rotate(self, R):
        r"""Apply a rotation matrix *R* to the selection's coordinates.

        AtomGroup.rotate(R)

        :math:`\mathsf{R}` is a 3x3 orthogonal matrix that transforms a vector
        :math:`\mathbf{x} \rightarrow \mathbf{x}'`:

        .. math::

            \mathbf{x}' = \mathsf{R}\mathbf{x}
        """
        R = numpy.matrix(R, copy=False, dtype=numpy.float32)
        # changes the coordinates (in place)
        x = self.universe.trajectory.ts._pos
        idx = self.indices()
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
        alpha = numpy.radians(angle)
        try:
            sel1, sel2 = axis
            x1, x2 = sel1.centroid(), sel2.centroid()
            v = x2 - x1
            n = v / numpy.linalg.norm(v)
            if point is None:
                point = x1
        except (ValueError, AttributeError):
            n = numpy.asarray(axis)
        if point is None:
            p = self.centroid()
        else:
            try:
                p = point.centroid()
            except AttributeError:
                p = numpy.asarray(point)
        M = rotation_matrix(alpha, n, point=p)
        self.transform(M)
        return M

    def align_principalAxis(self, axis, vector):
        """Align principal axis with index *axis* with *vector*.

        :Arguments:
          *axis*
            Index of the principal axis (0, 1, or 2), as produced by
            :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.principalAxes`.
          *vector*
            A 3D vector such as the z-axis (``[0,0,1]``); can be
            anything that looks like a list with three entries.

        To align the long axis of a channel (the first principal axis,
        i.e. *axis* = 0) with the z-axis::

          u.atoms.align_principalAxis(0, [0,0,1])
          u.atoms.write("aligned.pdb")
        """
        p = self.principalAxes()[axis]
        angle = numpy.degrees(vecangle(p, vector))
        ax = rotaxis(p, vector)
        #print "principal[%d] = %r" % (axis, p)
        #print "axis = %r, angle = %f deg" % (ax, angle)
        return self.rotateby(angle, ax)

    def packIntoBox(self, box=None, inplace=True):
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
            if (box.diagonal() == 0.0).any():  # for a vector representation, diagonal cannot be zero
                raise ValueError("One or more box dimensions is zero.  You can specify a boxsize with 'box ='")
        else:
            if (box == 0).any():  #Check that a box dimension isn't zero
                raise ValueError("One or more box dimensions is zero.  You can specify a boxsize with 'box='")

        coords = self.universe.trajectory.ts._pos[self.indices()]
        if not inplace:
            return MDAnalysis.core.distances.applyPBC(coords, box)

        self.universe.trajectory.ts._pos[self.indices()] = MDAnalysis.core.distances.applyPBC(coords, box)

        return self.universe.trajectory.ts._pos[self.indices()]

    def selectAtoms(self, sel, *othersel, **selgroups):
        """Selection of atoms using the MDAnalysis selection syntax.

        AtomGroup.selectAtoms(selection[,selection[,...]], [groupname=atomgroup[,groupname=atomgroup[,...]]])

        .. SeeAlso:: :meth:`Universe.selectAtoms`
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
            # use own list comprehension to avoid sorting/compression by eg self.resids()
            ids = numpy.array([getattr(atom, accessors[level]) for atom in self])
        except KeyError:
            raise ValueError("level = '{0}' not supported, must be one of {1}".format(
                    level, accessors.keys()))

        # now sort the resids so that it doesn't matter if the AG contains
        # atoms in random order (i.e. non-sequential resids); groupby needs
        # presorted keys!
        idx = numpy.argsort(ids)
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
                PDB, CRD, GRO, VMD (tcl), PyMol (pml), Gromacs (ndx) CHARMM (str);
                case-insensitive and can also be supplied as the filename
                extension [PDB]
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

        trj = self.universe.trajectory  # unified trajectory API
        frame = trj.ts.frame

        if trj.numframes == 1: kwargs.setdefault("multiframe", False)

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
            raise ValueError("No writer found for format: {}".format(filename))
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

    # properties
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
        if self.__ts is None or self.__ts.frame != trj_ts.frame:
            # create a timestep of same type as the underlying trajectory
            self.__ts = trj_ts.copy_slice(self.indices())
        return self.__ts


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
        sel = self.universe.selectAtoms(
            'segid %s and resid %d and name C' % (self.segment.id, self.id - 1)) + \
              self.N + self.CA + self.C
        if len(sel) == 4:  # selectAtoms doesnt raise errors if nothing found, so check size
            return sel
        else:
            return None

    def psi_selection(self):
        """AtomGroup corresponding to the psi protein backbone dihedral N-CA-C-N'.

        :Returns: 4-atom selection in the correct order.  If no N'
                  found in the following residue (by resid) then this
                  method returns ``None``.
        """
        sel = self.N + self.CA + self.C + \
              self.universe.selectAtoms(
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
        sel = self.CA + self.C + \
              self.universe.selectAtoms(
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
            return self.N + self.CA + self.CB + self.CG
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

    def set_resid(self, resid):
        """Set the resid to integer *resid* for **all residues** in the :class:`ResidueGroup`.

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
        """
        from MDAnalysis.topology.core import build_residues

        self.set('resid', resid)
        residues = build_residues(self.atoms)
        # update THIS residue group
        self.residues._residues = residues
        # make sure to update the whole universe: the Atoms are shared but ResidueGroups are not
        if self.atoms is not self.universe.atoms:
            self.universe.atoms._fill_cache('residues', build_residues(self.universe.atoms))

    def set_resnum(self, resnum):
        """Set the resnum to *resnum* for **all residues** in the :class:`ResidueGroup`.

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
        """
        self.set("resnum", resnum)

    def set_resname(self, resname):
        """Set the resname to string *resname* for **all residues** in the :class:`ResidueGroup`.

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
        """
        self.set("resname", resname, conversion=str)

    # All other AtomGroup.set_xxx() methods should just work as
    # ResidueGroup.set_xxx() because we overrode self.set(); the ones above
    # where kept separate because we can save a call to build_residues()
    # because there is no ambiguity as which residues are changed.

    def __repr__(self):
        return "<ResidueGroup {res}>".format(
            res=repr(self._residues))


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
        for res in self._residues:
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
            for res in self._residues:
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

    def set_segid(self, segid, buildsegments=True):
        """Set the segid to *segid* for all atoms in the :class:`SegmentGroup`.

        If *segid* is a sequence of the same length as the :class:`SegmentGroup`
        then each :attr:`Atom.segid` is set to the corresponding value together
        with the :attr:`Segment.id` of the segment the atom belongs to. If
        *value* is neither of length 1 (or a scalar) nor of the length of the
        :class:`AtomGroup` then a :exc:`ValueError` is raised.

        .. Note::

           :meth:`set_segid` can change the topology.

           With the default *buildsegments* = ``True`` it can be used to join
           segments or to break groups into multiple disjoint segments. Note
           that each :class:`Atom` can only belong to a single
           :class:`Segment`.

        For performance reasons, *buildsegments* can be set to ``False``. Then
        one needs to run :meth:`Universe._build_segments` manually later in
        order to update the list of :class:`Segment` instances and regenerate
        the segid instant selectors.

        .. Warning::

           The values of *this* :class:`SegmentGroup` are not being changed
           (i.e. if you assign multiple segids this instance will not be broken
           in multiple segments, rather you will have one :class:`SegmentGroup`
           that groups multiple segments together). You **must create a new**
           :class:`SegmentGroup` **from the** :class:`Universe` --- only
           :class:`Atom` instances are changed, everything else is derived from
           these atoms.


        .. versionadded:: 0.7.4
        .. versionchanged:: 0.8
           Can set atoms and residues to distinct values by providing a sequence or iterable.
        """
        super(SegmentGroup, self).set_segid(segid, buildsegments=buildsegments)

        # Is the following needed? -- orbeckst
        # also fix self --- otherwise users will get confused if the changes are not reflected in the
        # object they are currently using (it works automatically for AtomGroup but not higher order groups)
        #
        # This is a hack to be able to set properties on Segment/Atom
        # instances where they have different names
        #attr = {'segid': 'id'}
        for seg, value in itertools.izip(self.segments, itertools.cycle(util.asiterable(segid))):
            setattr(seg, 'name', value)

    def __getattr__(self, attr):
        if attr.startswith('s') and attr[1].isdigit():
            attr = attr[1:]  # sNxxx only used for python, the name is stored without s-prefix
        seglist = [segment for segment in self._segments if segment.name == attr]
        if len(seglist) == 0:
            return super(SegmentGroup, self).__getattr__(attr)
        if len(seglist) > 1:
            warnings.warn("SegmentGroup: Multiple segments with the same name {}; a combined, NON-CONSECUTIVE "
                          "Segment is returned.".format(attr), category=SelectionWarning)
            #return Segment(sum([s.residues for s in seglist])) ### FIXME: not working yet, need __add__
            return seglist[0]
        return seglist[0]

    def __repr__(self):
        return "<SegmentGroup {segnames}>".format(
            segnames=repr(self._segments))


class Universe(object):
    """The MDAnalysis Universe contains all the information describing the system.

    The system always requires a *topology* file --- in the simplest case just
    a list of atoms. This can be a CHARMM/NAMD PSF file or a simple coordinate
    file with atom informations such as PDB, Gromacs GRO, or CHARMM CRD. See
    :ref:`Supported topology formats` for what kind of topologies can be read.

    A trajectory provides coordinates; the coordinates have to be ordered in
    the same way as the list of atoms in the topology. A trajectory can be a
    single frame such as a PDB, CRD, or GRO file, or it can be a MD trajectory
    (in CHARMM/NAMD/LAMMPS DCD, Gromacs XTC/TRR, or generic XYZ format).  See
    :ref:`Supported coordinate formats` for what can be read as a
    "trajectory".

    As a special case, when the topology is a PDB, GRO or CRD file
    then the coordinates are immediately loaded from the "topology"
    file unless a trajectory is supplied.

    Examples for setting up a universe::

       u = Universe(topology, trajectory)          # read system from file(s)
       u = Universe(pdbfile)                       # read atoms and coordinates from PDB or GRO
       u = Universe(topology, [traj1, traj2, ...]) # read from a list of trajectories
       u = Universe(topology, traj1, traj2, ...)   # read from multiple trajectories

    Load new data into a universe (replaces old trajectory and does *not* append)::

       u.load_new(trajectory)                      # read from a new trajectory file

    Select atoms, with syntax similar to CHARMM (see
    :class:`~Universe.selectAtoms` for details)::

       u.selectAtoms(...)

    *Attributes:*

    - :attr:`Universe.trajectory`: currently loaded trajectory reader;
      :attr:`Universe.trajectory.ts` is the current time step
    - :attr:`Universe.dimensions`: current system dimensions (simulation unit cell, if
      set in the trajectory)
    - bonds, angles, dihedrals, impropers (low level access through :attr:`Universe._psf`)

    .. Note::

       If atom attributes such as element, mass, or charge are not explicitly
       provided in the topology file then MDAnalysis tries to guess them (see
       :mod:`MDAnalysis.topology.tables`). This does not always work and if you
       require correct values (e.g. because you want to calculate the center of
       mass) then you need to make sure that MDAnalysis gets all the
       information needed. Furthermore, the list of bonds is only constructed
       when provided in the topology and never guessed (see `Issue 23`).

    .. _`Issue 23`: http://code.google.com/p/mdanalysis/issues/detail?id=23

    .. versionchanged:: 0.7.5
       Can also read multi-frame PDB files with the
       :class:`~MDAnalysis.coordinates.PDB.PrimitivePDBReader`.

    .. versionchanged:: 0.8
       Parse arbitrary number of arguments as a single topology file and a a sequence
       of trajectories.

    .. versionchanged:: 0.9.0
       Topology information now loaded lazily, but can be forced with
       :meth:`build_topology`
       Changed .bonds attribute to be a :class:`~MDAnalysis.topology.core.TopologyGroup`
       Added .angles and .torsions attribute as :class:`~MDAnalysis.topology.core.TopologyGroup`
       Added fragments to Universe cache
    """

    def __init__(self, *args, **kwargs):
        """Initialize the central MDAnalysis Universe object.

        :Arguments:
          *topologyfile*
             A CHARMM/XPLOR PSF topology file, PDB file or Gromacs GRO file; used to define the
             list of atoms. If the file includes bond information, partial
             charges, atom masses, ... then these data will be available to
             MDAnalysis. A "structure" file (PSF, PDB or GRO, in the sense of a
             topology) is always required.
          *coordinatefile*
             A trajectory (such as CHARMM DCD, Gromacs XTC/TRR/GRO, XYZ, XYZ.bz2) or a PDB that
             will provide coordinates, possibly multiple frames.
             If a **list of filenames** is provided then they are sequentially read and appear
             as one single trajectory to the Universe. The list can contain different file
             formats.

             .. deprecated:: 0.8
                Do not use the *coordinatefile* keyword argument, just provide trajectories as
                positional arguments.

          *permissive*
             currently only relevant for PDB files: Set to ``True`` in order to ignore most errors
             and read typical MD simulation PDB files; set to ``False`` to read with the Bio.PDB reader,
             which can be useful for real Protein Databank PDB files. ``None``  selects the
             MDAnalysis default (which is set in :class:`MDAnalysis.core.flags`) [``None``]
          *topology_format*
             provide the file format of the topology file; ``None`` guesses it from the file
             extension [``None``]
             Can also pass a subclass of :class:`MDAnalysis.topology.base.TopologyReader`
             to define a custom reader to be used on the topology file.
          *format*
             provide the file format of the coordinate or trajectory file;
             ``None`` guesses it from the file extension. Note that this
             keyword has no effect if a list of file names is supplied because
             the "chained" reader has to guess the file format for each
             individual list member. [``None``]
             Can also pass a subclass of :class:`MDAnalysis.coordinates.base.Reader`
             to define a custom reader to be used on the trajectory file.
          *bonds*
             bond handling for PDB files. The default is to read and store the
             CONECT records only. When set to ``True`` it will attempt to guess
             connectivity between all atoms in the Universe.
             Each bond knows if it was guessed or was a CONECT record, so when
             saving out one can specify which ones to write out by ::

               u = Universe("example.pdb")
               u.atoms.write("output.pdb", bonds="conect") # default, only CONECT
               u.atoms.write("output.pdb", bonds="all")
               u.atoms.write("output.pdb", bonds=None)


        This routine tries to do the right thing:

        1. If a pdb/gro file is provided instead of a psf and no *coordinatefile*
           then the coordinates are taken from the first file. Thus you can load
           a functional universe with ::

              u = Universe('1ake.pdb')

           If you want to specify the coordinate file format yourself you can
           do so using the *format* keyword::

              u = Universe('1ake.ent1', format='pdb')

        2. If only a topology file without coordinate information is provided
           one will have to load coordinates manually using
           :meth:`Universe.load_new`. The file format of the topology file
           can be explicitly set with the *topology_format* keyword.

        .. versionchanged:: 0.7.4
           New *topology_format* and *format* parameters to override the file
           format detection.
        """
        from MDAnalysis.topology.core import get_parser_for, guess_format
        from ..topology.base import TopologyReader

        # managed attribute holding Reader
        self.__trajectory = None

        # Cache is used to store objects which are built lazily into Universe
        # Currently cached objects (managed property name and cache key):
        # - bonds
        # - angles
        # - torsions
        # - improper torsions
        # - fragments
        # Cached stuff is handled using util.cached decorator
        self._cache = dict()

        if len(args) == 0:
            # create an empty universe
            self.atoms = AtomGroup([])
            return

        self.filename = args[0]

        # old behaviour (explicit coordfile) overrides new behaviour
        coordinatefile = kwargs.pop('coordinatefile', args[1:])
        topology_format = kwargs.pop('topology_format', None)

        if len(args) == 1 and not coordinatefile:
            # special hacks to treat a coordinate file as a coordinate AND topology file
            # coordinatefile can be None or () (from an empty slice args[1:])
            if kwargs.get('format', None) is None:
                kwargs['format'] = topology_format
            elif topology_format is None:
                topology_format = kwargs.get('format', None)
            if guess_format(self.filename, format=kwargs.get('format', None)) in \
                    MDAnalysis.coordinates._topology_coordinates_readers:
                coordinatefile = self.filename  # hack for pdb/gro/crd - only
            # Fix by SB: make sure coordinatefile is never an empty tuple
            if len(coordinatefile) == 0:
                coordinatefile = None

        # build the topology (or at least a list of atoms)
        try:  # Try and check if the topology format is a TopologyReader
            if issubclass(topology_format, TopologyReader):
                parser = topology_format
        except TypeError:  # But strings/None raise TypeError in issubclass
            perm = kwargs.get('permissive',
                              MDAnalysis.core.flags['permissive_pdb_reader'])
            parser = get_parser_for(self.filename,
                                    permissive=perm,
                                    tformat=topology_format)
        try:
            with parser(self.filename,
                        guess_bonds_mode=kwargs.get('bonds', False)) as p:
                self._psf = p.parse()
        except IOError as err:
            raise IOError("Failed to load from the topology file {}"
                          " with parser {}.\n"
                          "Error: {}".format(self.filename, parser, err))
        except ValueError as err:
            raise ValueError("Failed to construct topology from file {}"
                             " with parser {} \n"
                             "Error: {}".format(self.filename, parser, err))

        # Generate atoms, residues and segments
        self._init_topology()

        # Load coordinates
        self.load_new(coordinatefile, **kwargs)

    def _clear_caches(self, *args):
        """Clear cache for all *args*.

        If not args are provided, all caches are cleared.

        .. versionadded 0.9.0
        """
        if len(args) == 0:
            self._cache = dict()
        else:
            for name in args:
                try:
                    del self._cache[name]
                except KeyError:
                    pass

    def _fill_cache(self, name, value):
        """Populate _cache[name] with value.

        .. versionadded:: 0.9.0
        """
        self._cache[name] = value

    def _init_topology(self):
        """Populate Universe attributes from the structure dictionary *_psf*."""
        self.atoms = AtomGroup(self._psf["_atoms"])

        # XXX: add H-bond information here if available from psf (or other sources)
        # segment instant selectors
        self._build_segments()

        # Let atoms access the universe
        for a in self.atoms:
            a.universe = self

    def _build_segments(self):
        """Parse list of atoms into segments.

        * updates :attr:`Universe.atoms` as a side effect
        * updates :attr:`Universe.segments` and :attr:`Universe.residues`
        * creates the segid instant selectors

        Because of Python's syntax rules, attribute names cannot start with a
        digit and so we prefix any segments starting with a digit with the
        letter 's'. For instance, '4AKE' becomes the segid instant selector
        's4AKE'.
        """
        from MDAnalysis.topology.core import build_segments

        segments = build_segments(self.atoms)
        for seg in segments.keys():
            if seg[0].isdigit():
                newsegname = 's' + seg
                segments[newsegname] = segments.pop(seg)
        self.__dict__.update(segments)
        # convenience access to residues and segments (these are managed attributes
        # (properties) and are built on the fly or read from a cache) -- does this
        # create memory problems?
        self.segments = self.atoms.segments
        self.residues = self.atoms.residues
        self.universe = self  # for Writer.write(universe), see Issue 49

    def _init_bonds(self):
        """Set bond information from u._psf['_bonds']

        .. versionchanged 0.9.0
           Now returns a :class:`~MDAnalysis.topology.core.TopologyGroup`
           Now only accepts list of tuples as input, previously accepted either
           lists of tuples or lists of Bonds.
        """
        from MDAnalysis.topology.core import Bond, TopologyGroup

        def fix_order(bondset):
            """Make sure that all bonds have first index less than second"""
            return set([tuple(sorted(b)) for b in bondset])

        bondlist = set()

        defined_bonds = fix_order(self._psf.get('_bonds', set()))
        guessed_bonds = fix_order(self._psf.get('_guessed_bonds', set()))
        # Some topologies define an order for bonds, this is stored
        # as dict of bondtuple:value
        bondorder = self._psf.get('_bondorder', False)
        # Defined bonds take priority, remove bonds in 'guessed' that are in 'defined'
        guessed_bonds.difference_update(defined_bonds)

        for bondset, guessed in zip([defined_bonds, guessed_bonds],
                                    [False, True]):
            for bond in bondset:
                i, j = bond[0], bond[1]

                if bondorder:
                    order = bondorder.get(bond, None)
                else:
                    order = None

                a1, a2 = self.atoms[i], self.atoms[j]
                bond = Bond([a1, a2], order=order)
                bond.is_guessed = guessed

                bondlist.add(bond)

        if len(bondlist) > 0:
            return TopologyGroup(list(bondlist))
        else:
            return None

    def _init_angles(self):
        """Builds angle information from u._psf['_angles']

        Returns ``None`` if no angle information is present, otherwise
        returns a :class:`~MDAnalysis.topology.core.TopologyGroup`

        .. versionadded 0.9.0
        """
        angle_entries = self._psf.get('_angles', None)
        if angle_entries is None:
            return None
        else:
            from MDAnalysis.topology.core import TopologyGroup, Angle

            angle_list = [Angle([self.atoms[a] for a in entry])
                          for entry in angle_entries]
            if len(angle_list) > 0:
                return TopologyGroup(angle_list)
            else:
                return None

    def _init_torsions(self):
        """Builds torsion information from u._psf['_dihe']

        Returns ``None`` if no torsion information is present, otherwise
        returns a :class:`~MDAnalysis.topology.core.TopologyGroup`

        .. versionadded 0.9.0
        """
        torsion_entries = self._psf.get('_dihe', None)
        if torsion_entries is None:
            return None
        else:
            from MDAnalysis.topology.core import TopologyGroup, Torsion

            torsion_list = [Torsion([self.atoms[a] for a in entry])
                            for entry in torsion_entries]
            if len(torsion_list) > 0:
                return TopologyGroup(torsion_list)
            else:
                return None

    def _init_impropers(self):
        """Build improper torsion information from u._psf['_impr']

        Returns ``None`` if no improper torsion information is present,
        otherwise returns a :class:`~MDAnalysis.topology.core.TopologyGroup`

        .. versionadded 0.9.0
        """
        torsion_entries = self._psf.get('_impr', None)
        if torsion_entries is None:
            return None
        else:
            from MDAnalysis.topology.core import TopologyGroup, Improper_Torsion

            torsion_list = [Improper_Torsion([self.atoms[a] for a in entry])
                            for entry in torsion_entries]
            if len(torsion_list) > 0:
                return TopologyGroup(torsion_list)
            else:
                return None

    def _init_fragments(self):
        """Build all fragments in the Universe

        Generally built on demand by an Atom querying its fragment property.

        .. versionadded 0.9.0
        """
        # Check that bond information is present, else inform
        bonds = self.bonds
        if bonds is None:
            raise NoDataError("Fragments require that the Universe has Bond information")

        # This current finds all fragments from all Atoms
        # Could redo this to only find fragments for a queried atom (ie. only fill out
        # a single fragment).  This would then make it scale better for large systems.
        # eg:
        # try:
        #    return self._fragDict[a]
        # except KeyError:
        #    self._init_fragments(a)  # builds the fragment a belongs to

        class _fragset(object):
            """Normal sets aren't hashable, this is"""

            def __init__(self, ats):
                self.ats = set(ats)

            def __iter__(self):
                return iter(self.ats)

            def add(self, other):
                self.ats.add(other)

            def update(self, other):
                self.ats.update(other.ats)

        f = {a: None for a in self.atoms}  # each atom starts with its own list

        for a1, a2 in bonds:  # Iterate through all bonds
            if not (f[a1] or f[a2]):  # New set made here
                new = _fragset([a1, a2])
                f[a1] = f[a2] = new
            elif f[a1] and not f[a2]:  # If a2 isn't in a fragment, add it to a1's
                f[a1].add(a2)
                f[a2] = f[a1]
            elif not f[a1] and f[a2]:  # If a1 isn't in a fragment, add it to a2's
                f[a2].add(a1)
                f[a1] = f[a2]
            elif f[a1] is f[a2]:  # If they're in the same fragment, do nothing
                continue
            else:  # If they are both in different fragments, combine fragments
                f[a1].update(f[a2])
                f.update({a: f[a1] for a in f[a2]})

                # Lone atoms get their own fragment
        f.update({a: _fragset((a,)) for a, val in f.items() if not val})

        # All the unique values in f are the fragments
        frags = tuple([AtomGroup(list(a.ats)) for a in set(f.values())])

        return frags

    @property
    @cached('fragments')
    def fragments(self):
        """Read only tuple of fragments in the Universe

        .. versionadded 0.9.0
        """
        return self._init_fragments()

    @property
    @cached('bondDict')
    def _bondDict(self):
        """Lazily built dictionary of bonds

        Translates Atom to list of bonds

        .. versionadded:: 0.9.0
        """
        bonds = self.bonds
        bd = defaultdict(list)

        if bonds is None:
            pass
        else:
            for b in bonds:
                for a in b:
                    bd[a].append(b)
        return bd

    @property
    @cached('angleDict')
    def _angleDict(self):
        """Lazily built dictionary of angles

        Translates Atom to list of angles

        .. versionadded:: 0.9.0
        """
        bonds = self.angles
        bd = defaultdict(list)

        if bonds is None:
            pass
        else:
            for b in bonds:
                for a in b:
                    bd[a].append(b)
        return bd

    @property
    @cached('torsionDict')
    def _torsionDict(self):
        """Lazily built dictionary of torsions

        Translates Atom to list of torsions

        .. versionadded:: 0.9.0
        """
        bonds = self.torsions
        bd = defaultdict(list)

        if bonds is None:
            pass
        else:
            for b in bonds:
                for a in b:
                    bd[a].append(b)
        return bd

    @property
    @cached('improperDict')
    def _improperDict(self):
        """Lazily built dictionary of improper torsions

        Translates Atom to list of improper torsions

        .. versionadded:: 0.9.0
        """
        bonds = self.impropers
        bd = defaultdict(list)

        if bonds is None:
            pass
        else:
            for b in bonds:
                for a in b:
                    bd[a].append(b)
        return bd

    @property
    @cached('fragDict')
    def _fragmentDict(self):
        """Lazily built dictionary of fragments.

        Translates :class:`Atom` objects into the fragment they belong to.

        The Atom.fragment managed property queries this dictionary.

        .. versionadded 0.9.0
        """
        frags = self.fragments  # will build if not built
        fd = dict()
        for f in frags:
            for a in f:
                fd[a] = f
        return fd

    def build_topology(self):
        """
        Bond angle and torsion information is lazily constructed into the
        Universe.

        This method forces all this information to be loaded.

        .. versionadded 0.9.0
        """
        if 'bonds' not in self._cache:
            self._cache['bonds'] = self._init_bonds()
        if 'angles' not in self._cache:
            self._cache['angles'] = self._init_angles()
        if 'torsions' not in self._cache:
            self._cache['torsions'] = self._init_torsions()
        if 'impropers' not in self._cache:
            self._cache['impropers'] = self._init_impropers()

    @property
    @cached('bonds')
    def bonds(self):
        """
        Returns a :class:`~MDAnalysis.topology.core.TopologyGroup` of all
        bonds in the Universe.
        If none are found, returns ``None``

        .. versionchaged 0.9.0
           Now a lazily built :class:`~MDAnalysis.topology.core.TopologyGroup`
        """
        return self._init_bonds()

    @bonds.setter
    def bonds(self, bondlist):
        """Can set bonds by supplying an iterable of bond tuples.

        Each bond tuple must contain the zero based indices of the two Atoms in
        the bond

        .. versionadded:: 0.9.0
        """
        del self.bonds
        self._psf['_bonds'] = bondlist
        self._fill_cache('bonds', self._init_bonds())

    @bonds.deleter
    def bonds(self):
        """Delete the bonds from Universe

        This must also remove the per atom record of bonds (bondDict)

        .. versionadded:: 0.9.0
        """
        self._clear_caches('bonds', 'bondDict')

    @property
    @cached('angles')
    def angles(self):
        """
        Returns a :class:`~MDAnalysis.topology.core.TopologyGroup`
        of all angles in the Universe

        If none are found, returns ``None``

        .. versionadded 0.9.0
        """
        return self._init_angles()

    @angles.setter
    def angles(self, bondlist):
        del self.angles
        self._psf['_angles'] = bondlist
        self._fill_cache('angles', self._init_angles())

    @angles.deleter
    def angles(self):
        self._clear_caches('angles', 'angleDict')

    @property
    @cached('torsions')
    def torsions(self):
        """
        Returns a :class:`~MDAnalysis.topology.core.TopologyGroup`
        of all torsions in the Universe

        If none are found, returns ``None``

        .. versionadded 0.9.0
        """
        return self._init_torsions()

    @torsions.setter
    def torsions(self, bondlist):
        del self.torsions
        self._psf['_dihe'] = bondlist
        self._fill_cache('torsions', self._init_torsions())

    @torsions.deleter
    def torsions(self):
        self._clear_caches('torsions', 'torsionDict')

    @property
    @cached('impropers')
    def impropers(self):
        """
        Returns a :class:`~MDAnalysis.topology.core.TopologyGroup`
        of all improper torsions in the Universe

        If none are found, returns ``None``

        .. versionadded 0.9.0
        """
        return self._init_impropers()

    @impropers.setter
    def impropers(self, bondlist):
        del self.impropers
        self._psf['_impr'] = bondlist
        self._fill_cache('impropers', self._init_impropers())

    @impropers.deleter
    def impropers(self):
        self._clear_caches('impropers', 'improperDict')

    def load_new(self, filename, **kwargs):
        """Load coordinates from *filename*, using the suffix to detect file format.

        :Arguments:
             *filename*
                 the coordinate file (single frame or trajectory) *or* a list of
                 filenames, which are read one after another.
             *permissive*
                 currently only relevant for PDB files: Set to ``True`` in order to ignore most errors
                 and read typical MD simulation PDB files; set to ``False`` to read with the Bio.PDB reader,
                 which can be useful for real Protein Databank PDB files. ``None``  selects the
                 MDAnalysis default (which is set in :class:`MDAnalysis.core.flags`) [``None``]
             *format*
                 provide the file format of the coordinate or trajectory file;
                 ``None`` guesses it from the file extension. Note that this
                 keyword has no effect if a list of file names is supplied because
                 the "chained" reader has to guess the file format for each
                 individual list member [``None``]
                 Can also pass a subclass of :class:`MDAnalysis.coordinates.base.Reader`
                 to define a custom reader to be used on the trajectory file.
             *kwargs*
                 Other kwargs are passed to the trajectory reader (only for advanced use)

        :Returns: (filename, trajectory_format) or ``None`` if *filename* == ``None``
        :Raises: :exc:`TypeError` if trajectory format can not be
                  determined or no appropriate trajectory reader found

        .. versionchanged:: 0.8
           If a list or sequence that is provided for *filename*  only contains a single entry
           then it is treated as single coordinate file. This has the consequence that it is
           not read by the :class:`~MDAnalysis.coordinates.base.ChainReader` but directly by
           its specialized file format reader, which typically has more features than the
           :class:`~MDAnalysis.coordinates.base.ChainReader`.
        """
        if filename is None:
            return

        import MDAnalysis.core
        from ..coordinates.core import get_reader_for
        from ..coordinates.base import Reader

        if len(util.asiterable(filename)) == 1:
            # make sure a single filename is not handed to the ChainReader
            filename = util.asiterable(filename)[0]
        logger.debug("Universe.load_new(): loading {0}...".format(filename))

        reader_format = kwargs.pop('format', None)
        perm = kwargs.get('permissive', MDAnalysis.core.flags['permissive_pdb_reader'])
        reader = None
        try:
            if reader_format is not None and issubclass(reader_format, Reader):
                reader = reader_format
        except TypeError:
            pass
        if not reader:
            try:
                reader = get_reader_for(filename,
                                        permissive=perm,
                                        format=reader_format)
            except TypeError as err:
                raise TypeError("Cannot find an appropriate coordinate reader "
                                "for file {}.\n{}".format(filename, err))

        # supply number of atoms for readers that cannot do it for themselves
        kwargs['numatoms'] = self.atoms.numberOfAtoms()
        self.trajectory = reader(filename, **kwargs)    # unified trajectory API
        if self.trajectory.numatoms != self.atoms.numberOfAtoms():
            raise ValueError("The topology and {} trajectory files don't"
                             " have the same number of atoms!\n"
                             "Topology number of atoms {}\n"
                             "Trajectory: {} Number of atoms {}".format(
                                 self.trajectory.format,
                                 len(self.atoms),
                                 filename, self.trajectory.numatoms))

        return filename, self.trajectory.format

    def selectAtoms(self, sel, *othersel, **selgroups):
        """Select atoms using a CHARMM selection string.

        Returns an :class:`AtomGroup` with atoms sorted according to their
        index in the psf (this is to ensure that there aren't any duplicates,
        which can happen with complicated selections).

        Existing :class:`AtomGroup` objects can be passed as named arguments,
        which will then be available to the selection parser.

        Subselections can be grouped with parentheses.

        Example::
           >>> sel = universe.selectAtoms("segid DMPC and not ( name H* or name O* )")
           >>> sel
           <AtomGroup with 3420 atoms>

           >>> universe.selectAtoms("around 10 group notHO", notHO=sel)
           <AtomGroup with 1250 atoms>

        .. Note::

           If exact ordering of atoms is required (for instance, for
           :meth:`~AtomGroup.angle` or :meth:`~AtomGroup.dihedral`
           calculations) then one supplies selections *separately* in the
           required order. Also, when multiple :class:`AtomGroup` instances are
           concatenated with the ``+`` operator then the order of :class:`Atom`
           instances is preserved and duplicates are not removed.

        .. SeeAlso:: :ref:`selection-commands-label` for further details and examples.

        The selection parser understands the following CASE SENSITIVE *keywords*:

        **Simple selections**

            protein, backbone, nucleic, nucleicbackbone
                selects all atoms that belong to a standard set of residues; a protein
                is identfied by a hard-coded set of residue names so it  may not
                work for esoteric residues.
            segid *seg-name*
                select by segid (as given in the topology), e.g. ``segid 4AKE`` or ``segid DMPC``
            resid *residue-number-range*
                resid can take a single residue number or a range of numbers. A range
                consists of two numbers separated by a colon (inclusive) such
                as ``resid 1:5``. A residue number ("resid") is taken directly from the
                topology.
            resnum *resnum-number-range*
                resnum is the canonical residue number; typically it is set to the residue id
                in the original PDB structure.
            resname *residue-name*
                select by residue name, e.g. ``resname LYS``
            name *atom-name*
                select by atom name (as given in the topology). Often, this is force
                field dependent. Example: ``name CA`` (for C&alpha; atoms) or ``name OW`` (for SPC water oxygen)
            type *atom-type*
                select by atom type; this is either a string or a number and depends on
                the force field; it is read from the topology file (e.g. the CHARMM PSF
                file contains numeric atom types). It has non-sensical values when a
                PDB or GRO file is used as a topology.
            atom *seg-name*  *residue-number*  *atom-name*
                a selector for a single atom consisting of segid resid atomname,
                e.g. ``DMPC 1 C2`` selects the C2 carbon of the first residue of the DMPC
                segment
            altloc *alternative-location*
                a selection for atoms where alternative locations are available,
                which is often the case with high-resolution crystal structures
                e.g. `resid 4 and resname ALA and altloc B` selects only the atoms
                of ALA-4 that have an altloc B record.

        **Boolean**

            not
                all atoms not in the selection, e.g. ``not protein`` selects
                all atoms that aren't part of a protein

            and, or
                combine two selections according to the rules of boolean algebra,
                e.g. ``protein and not (resname ALA or resname LYS)`` selects all atoms
                that belong to a protein, but are not in a lysine or alanine residue

        **Geometric**

            around *distance*  *selection*
                selects all atoms a certain cutoff away from another selection,
                e.g. ``around 3.5 protein`` selects all atoms not belonging to protein
                that are within 3.5 Angstroms from the protein
            point *x* *y* *z*  *distance*
                selects all atoms within a cutoff of a point in space, make sure
                coordinate is separated by spaces, e.g. ``point 5.0 5.0 5.0  3.5`` selects
                all atoms within 3.5 Angstroms of the coordinate (5.0, 5.0, 5.0)
            prop [abs] *property*  *operator*  *value*
                selects atoms based on position, using *property*  **x**, **y**, or
                **z** coordinate. Supports the **abs** keyword (for absolute value) and
                the following *operators*: **<, >, <=, >=, ==, !=**. For example, ``prop z >= 5.0``
                selects all atoms with z coordinate greater than 5.0; ``prop abs z <= 5.0``
                selects all atoms within -5.0 <= z <= 5.0.

        **Connectivity**

            byres *selection*
                selects all atoms that are in the same segment and residue as
                selection, e.g. specify the subselection after the byres keyword

        **Index**

            bynum *index-range*
                selects all atoms within a range of (1-based) inclusive indices,
                e.g. ``bynum 1`` selects the first atom in the universe; ``bynum 5:10``
                selects atoms 5 through 10 inclusive. All atoms in the
                :class:`MDAnalysis.Universe` are consecutively numbered, and the index
                runs from 1 up to the total number of atoms.

        **Preexisting selections**

            group *group-name*
                selects the atoms in the :class:`AtomGroup` passed to the function as an
                argument named *group-name*. Only the atoms common to *group-name* and the
                instance :meth:`~selectAtoms` was called from will be considered.
                *group-name* will be included in the parsing just by comparison of atom indices.
                This means that it is up to the user to make sure they were defined in an
                appropriate :class:`Universe`.

            fullgroup *group-name*
                just like the ``group`` keyword with the difference that all the atoms of
                *group-name* are included. The resulting selection may therefore have atoms
                that were initially absent from the instance :meth:`~selectAtoms` was
                called from.


        .. versionchanged:: 0.7.4
           Added *resnum* selection.
        .. versionchanged:: 0.8.1
           Added *group* and *fullgroup* selections.
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

    def __repr__(self):
        return "<Universe with {natoms} atoms{bonds}>".format(
            natoms=len(self.atoms),
            bonds=" and {} bonds".format(len(self.bonds)) if self.bonds else "")

    def __getstate__(self):
        raise NotImplementedError

    def __setstate__(self, state):
        raise NotImplementedError

    # Properties
    @property
    def dimensions(self):
        """Current dimensions of the unitcell"""
        return self.coord.dimensions

    @dimensions.setter
    def dimensions(self, box):
        """Set dimensions if the Timestep allows this

        .. versionadded:: 0.9.0
        """
        # Add fancy error handling here or use Timestep?
        self.coord.dimensions = box

    @property
    def coord(self):
        """Reference to current timestep and coordinates of universe.

        The raw trajectory coordinates are :attr:`Universe.coord._pos`,
        represented as a :class:`numpy.float32` array.

        Because :attr:`coord` is a reference to a
        :class:`~MDAnalysis.coordinates.base.Timestep`, it changes its contents
        while one is stepping through the trajectory.

        .. Note::

           In order to access the coordinates it is probably better to use the
           :meth:`AtomGroup.coordinates` method; for instance, all coordinates
           of the Universe as a numpy array:
           :meth:`Universe.atoms.coordinates`.
        """
        return self.trajectory.ts

    @property
    def trajectory(self):
        """Reference to trajectory reader object containing trajectory data."""
        if not self.__trajectory is None:
            return self.__trajectory
        else:
            raise AttributeError("No trajectory loaded into Universe")

    @trajectory.setter
    def trajectory(self, value):
        del self.__trajectory  # guarantees that files are closed (?)
        self.__trajectory = value

    # NOTE: DO NOT ADD A __del__() method: it somehow keeps the Universe
    #       alive during unit tests and the unit tests run out of memory!
    #### def __del__(self): <------ do not add this! [orbeckst]


def asUniverse(*args, **kwargs):
    """Return a universe from the input arguments.

    1. If the first argument is a universe, just return it::

         as_universe(universe) --> universe

    2. Otherwise try to build a universe from the first or the first
       and second argument::

         asUniverse(PDB, **kwargs) --> Universe(PDB, **kwargs)
         asUniverse(PSF, DCD, **kwargs) --> Universe(PSF, DCD, **kwargs)
         asUniverse(*args, **kwargs) --> Universe(*args, **kwargs)

    :Returns: an instance of :class:`~MDAnalaysis.AtomGroup.Universe`
    """
    if len(args) == 0:
        raise TypeError("asUniverse() takes at least one argument (%d given)" % len(args))
    elif len(args) == 1 and isinstance(args[0], Universe):
        return args[0]
    return Universe(*args, **kwargs)


def Merge(*args):
    """Return a :class:`Universe` from two or more :class:`AtomGroup` instances.

    :class:`AtomGroup` instances can come from different Universes, or come
    directly from a :meth:`~Universe.selectAtoms` call.

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
       u = Merge(u1.selectAtoms("protein"), u2.atoms, u3.atoms)
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

    coords = numpy.vstack([a.coordinates() for a in args])
    trajectory = MDAnalysis.coordinates.base.Reader()
    ts = MDAnalysis.coordinates.base.Timestep(coords)
    setattr(trajectory, "ts", ts)
    trajectory.numframes = 1

    # create an empty Universe object
    u = Universe()
    u.trajectory = trajectory

    # create a list of Atoms, then convert it to an AtomGroup
    atoms = [copy.copy(a) for gr in args for a in gr]
    for a in atoms:
        a.universe = u
    # adjust the atom numbering
    for i, a in enumerate(atoms):
        a.number = i
        a.serial = i + 1
    u.atoms = AtomGroup(atoms)
    # adjust the residue and segment numbering (removes any remaining references to the old universe)
    MDAnalysis.topology.core.build_residues(u.atoms)
    MDAnalysis.topology.core.build_segments(u.atoms)

    return u

