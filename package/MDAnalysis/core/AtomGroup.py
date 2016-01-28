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


