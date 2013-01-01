.. -*- coding: utf-8 -*-
.. _selection-commands-label:

====================
 Selection Commands
====================

Once you have the :meth:`~MDAnalysis.core.Universe` object, you can select
atoms (using a syntax very similar to `CHARMM's atom selection syntax`_)::

  >>> kalp = universe.selectAtoms("segid KALP")

.. _`CHARMM's atom selection syntax`: http://www.charmm.org/html/documentation/c35b1/select.html

:meth:`MDAnalysis.core.Universe.selectAtoms` returns a
:class:`MDAnalysis.core.AtomGroup.AtomGroup`, so you can use all the methods
defined for AtomGroups on them. Selections always return an AtomGroup with
atoms sorted according to their index in the topology (this is to ensure that
there aren't any duplicates, which can happen with complicated selections).

One can group subselections using parentheses::

 >>> universe.selectAtoms("segid DMPC and not (name H* or name O*)")
 <AtomGroup with 3420 atoms>

Almost all the basic CHARMM selections work.

It is also possible to export selections for external software
packages with the help of :ref:`Selection exporters`.


Selection Keywords
==================

The following describes all selection keywords currently understood by the
selection parser. The following applies to all selections:

* Keywords are case sensitive.
* Atoms are automatically sequentially ordered in a resulting selection
  (see notes below on :ref:`ordered-selections-label` for how to circumvent this if
  necessary).
* Selections are parsed left to right and parentheses can be used for
  grouping.
* Currently, only "stemming" is implemented as a primitive form of pattern
  matching: Using the ``*`` character in a string such as ``GL*`` selects
  all strings that start with "GL" such as "GLU", "GLY", "GLX29", "GLN".


Simple selections
-----------------

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

Boolean
-------

    not
        all atoms not in the selection, e.g. ``not protein`` selects all atoms that aren't part of a protein
    and, or
        combine two selections according to the rules of boolean algebra,
        e.g. ``protein and not (resname ALA or resname LYS)`` selects all atoms
        that belong to a protein, but are not in a lysine or alanine residue  

Geometric
---------

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

From version 0.6 onwards, geometric selections can use a k-d tree
based, fast search algorithm (about three times faster than the
previous version). However, it does not take periodicity into
account. The fast algorithm is the default for *around*. Periodicity
is only taken into account with the
:func:`~MDAnalysis.analysis.distances.distance_array` functions via a
minimum image convention (and this only works for rectangular
simulation cells). If periodic boundary conditions should be taken
into account then change the default behaviour of MDAnalysi by setting
these two flags::

  MDAnalysis.core.flags['use_periodic_selections'] = True
  MDAnalysis.core.flags['use_KDTree_routines'] = False


Connectivity
------------

    byres *selection*
        selects all atoms that are in the same segment and residue as
        selection, e.g. specify the subselection after the byres keyword  

Index
-----

    bynum *index-range*
        selects all atoms within a range of (1-based) inclusive indices,
        e.g. ``bynum 1`` selects the first atom in the universe; ``bynum 5:10``
        selects atoms 5 through 10 inclusive. All atoms in the
        :class:`MDAnalysis.Universe` are consecutively numbered, and the index
        runs from 1 up to the total number of atoms.


Instant selectors
=================

For interactive work it becomes rather tedious to type common selection strings
repeatedly. MDAnalysis automatically generates a number of *instant selectors*
as attributes of the :class:`~MDAnalysis.Universe` and number of other levels
of the structural hierarchy, namely for
:class:`~MDAnalysis.AtomGroup.AtomGroup`,
:class:`~MDAnalysis.AtomGroup.Residue`,
:class:`~MDAnalysis.AtomGroup.ResidueGroup`,
:class:`~MDAnalysis.AtomGroup.Segment` and
:class:`~MDAnalysis.AtomGroup.SegmentGroup`.

Segment selector
----------------

- ``universe.<segid>`` or ``universe.s<segid>`` (if *<segid>* starts with a
  number)
- returns a :class:`~MDAnalysis.AtomGroup.Segment`
- works for :class:`~MDAnalysis.Universe` and :class:`~MDAnalysis.AtomGroup.SegmentGroup`
- example
   >>> u.s4AKE
   <Segment '4AKE'>

Resid selector
--------------

- ``seg.r<N>`` selects residue with number ``<N>``
- returns a :class:`~MDAnalysis.AtomGroup.Residue`
- works for :class:`~MDAnalysis.AtomGroup.Segment` and :class:`~MDAnalysis.AtomGroup.SegmentGroup`
- example
    >>>  u.s4AKE.r100
    <Residue 'GLY', 100>
 
Residue name selector
---------------------

- ``seg.<resname>`` selects residues with residue name ``<resname>``
- returns a :class:`~MDAnalysis.AtomGroup.ResidueGroup`
- works for :class:`~MDAnalysis.AtomGroup.Segment` and :class:`~MDAnalysis.AtomGroup.SegmentGroup`
- examples
    >>> u.s4AKE.MET
    <ResidueGroup [<Residue 'MET', 1>, <Residue 'MET', 21>, <Residue 'MET', 34>, <Residue 'MET', 53>, <Residue 'MET', 96>, <Residue 'MET', 174>]>
    >>> u.s4AKE.CYS
    <ResidueGroup [<Residue 'CYS', 77>]>
    >>> u.s4AKE.TRP
    NoDataError: No atoms defined for AtomGroup
- The result is always a :class:`~MDAnalysis.AtomGroup.ResidueGroup`; if no
  residues can be found then a :exc:`MDAnalysis.NoDataError` is raised.

Atom name selector
------------------

- ``g.<atomname>`` selects a single atom or a group of atoms with name
  ``<atomname>``
- returns 
    - a :class:`~MDAnalysis.AtomGroup.Atom` if only a single atom was found,
    - a :class:`~MDAnalysis.AtomGroup.AtomGroup` if more than one atom was
      found, or
    - raises a :exc:`MDAnalysis.SelectionError` if no atom was found.
- works for any group derived from :class:`~MDAnalysis.AtomGroup.AtomGroup`
  (i.e. all the ones mentioned above)
- examples
    >>> u.atoms.CG
    >>> <AtomGroup with 125 atoms>
    >>> u.s4AKE.CG     
    <AtomGroup with 125 atoms>
    >>> u.s4AKE.r100.CA
    < Atom 1516: name 'CA' of type '23' of resname 'GLY', resid 100 and segid '4AKE'>
    >>> u.s4AKE.r100.CB
    SelectionError: No atom in residue GLY with name CB
  	

.. _ordered-selections-label:

Ordered selections
==================

:meth:`~MDAnalysis.Universe.selectAtoms` sorts the atoms in the
:class:`~MDAnalysis.core.AtomGroup.AtomGroup` by atom index before returning them (this is to
eliminate possible duplicates in the selection). If the ordering of atoms is
crucial (for instance when describing angles or dihedrals) or if duplicate
atoms are required then one has to concatenate multiple AtomGroups, which does
not sort them. 

The most straightforward way to concatentate two AtomGroups is by using the
**+** operator::

 >>> ordered = u.selectAtoms("segid DMPC and resid 3 and name P") + u.selectAtoms("segid DMPC and resid 2 and name P")
 >>> print list(ordered)
 [< Atom 570: name 'P' of type '180' of resid 'DMPC', 3 and 'DMPC'>,
 < Atom 452: name 'P' of type '180' of resid 'DMPC', 2 and 'DMPC'>]

A shortcut is to provide *two or more* selections to
:meth:`~MDAnalysis.Universe.selectAtoms`, which then does the concatenation
automatically::

 >>> print list(universe.selectAtoms("segid DMPC and resid 3 and name P", "segid DMPC and resid 2 and name P"))
 [< Atom 570: name 'P' of type '180' of resid 'DMPC', 3 and 'DMPC'>,
 < Atom 452: name 'P' of type '180' of resid 'DMPC', 2 and 'DMPC'>]

Just for comparison to show that a single selection string does not work as one
might expect::

 # WRONG!
 >>> print list(universe.selectAtoms("segid DMPC and ( resid 3 or resid 2 ) and name P"))
 [< Atom 452: name 'P' of type '180' of resid 'DMPC', 2 and 'DMPC'>,
 < Atom 570: name 'P' of type '180' of resid 'DMPC', 3 and 'DMPC'>]
 
