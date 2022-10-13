.. -*- coding: utf-8 -*-
.. _selection-commands-label:

====================
 Selection commands
====================

Once you have the :meth:`~MDAnalysis.core.universe.Universe` object, you can
select atoms (using a syntax very similar to `CHARMM's atom selection
syntax`_)::

  >>> kalp = universe.select_atoms("segid KALP")

.. _`CHARMM's atom selection syntax`:
   https://www.charmm.org/charmm/documentation/by-version/c45b1/select.html

The :meth:`~MDAnalysis.core.groups.AtomGroup.select_atoms` method of a
:class:`~MDAnalysis.core.groups.AtomGroup` or a
:class:`~MDAnalysis.core.universe.Universe` returns a
:class:`~MDAnalysis.core.groups.AtomGroup`, so you can use all the methods
defined for AtomGroups on them. Selections always return an
:class:`~MDAnalysis.core.groups.AtomGroup` with atoms sorted according to their
index in the topology (this is to ensure that there are not any duplicates,
which can happen with complicated selections).

One can group subselections using parentheses::

 >>> universe.select_atoms("segid DMPC and not (name H* or type OW)")
 <AtomGroup with 3420 atoms>

Almost all the basic CHARMM selections work.

It is also possible to export selections for external software
packages with the help of :ref:`Selection exporters`.

.. note::

    By default, atoms are sorted by index in the output AtomGroup.
    For example, the below code will return the first, second, and
    sixth atom in ``ag``::

        >>> ag = u.select_atoms("name N")
        >>> ag2 = ag[[5, 1, 0]]
        >>> ag3 = ag2.select_atoms("name N")
        >>> np.all(ag3.ix == ag2.ix)
        False

    You can turn off sorting behavior with the ``sorted`` keyword::

        >>> ag = u.select_atoms("name N")
        >>> ag2 = ag[[5, 1, 0]]
        >>> ag3 = ag2.select_atoms("name N", sorted=False)
        >>> np.all(ag3.ix == ag2.ix)
        True

    For further details on ordered selections, see :ref:`ordered-selections-label`.
    
    
Selection Keywords
==================

The following describes all selection keywords currently understood by the
selection parser. The following applies to all selections:

* Keywords are case sensitive.
* Atoms are automatically sequentially ordered in a resulting selection (see
  notes below on :ref:`ordered-selections-label` for how to circumvent this if
  necessary).
* Selections are parsed left to right and parentheses can be used for
  grouping.
* You can use the singular name of any topology attribute as a selection
  keyword. `Defined topology attributes`_ are listed in the User Guide.
  Alternatively, you can define a 
  :class:`~MDAnalysis.core.topologyattrs.TopologyAttr` yourself,
  providing that the attribute ``dtype`` is one of ``int``, ``float``, 
  ``str`` (or ``object``), or ``bool``.
  However, the topology must contain this attribute information for
  the selection to work.

    * Selections of attributes that are integers or floats can use the
      syntax "myTopologyAttr 0 - 2", "myTopologyAttr 0:2", or
      "myTopologyAttr 0 to 2", to select a range with
      both ends inclusive. Whitespace and negative numbers are allowed.
    * "myTopologyAttr 0" can be used to select all atoms
      matching the value; however, this can be tricky with floats because of
      precision differences and we recommend using a range like above when
      possible.
    * Boolean selections default to True, so "myTopologyAttr" and
      "myTopologyAttr True" both give all atoms with
      ``myTopologyAttr == True``.

.. seealso::

    Regular expression patterns
    :data:`~MDAnalysis.core.selection.FLOAT_PATTERN` for matching floats;
    :data:`~MDAnalysis.core.selection.INT_PATTERN` for matching integers;
    and :data:`~MDAnalysis.core.selection.RANGE_PATTERN` for matching
    selection ranges.


.. _`Defined topology attributes`: https://userguide.mdanalysis.org/stable/topology_system.html#format-specific-attributes


Simple selections
-----------------

This is a non-exhaustive list of the available selection keywords. As noted
in the dot point above, keywords will be automatically generated for any
suitable :class:`~MDAnalysis.core.topologyattrs.TopologyAttr`. A list of
`Defined topology attributes`_ is available in the User Guide.

protein, backbone, nucleic, nucleicbackbone
    selects all atoms that belong to a standard set of residues; a protein
    is identfied by a hard-coded set of residue names so it  may not
    work for esoteric residues.

segid *seg-name*
    select by segid (as given in the topology), e.g. ``segid 4AKE`` or
    ``segid DMPC``

resid *residue-number-range*
    resid can take a single residue number or a range of numbers. A range
    consists of two numbers separated by a colon (inclusive) such
    as ``resid 1:5``. A residue number ("resid") is taken directly from the
    topology.

resnum *resnum-number-range*
    resnum is the canonical residue number; typically it is set to the
    residue id in the original PDB structure.

resname *residue-name*
    select by residue name, e.g. ``resname LYS``

name *atom-name*
    select by atom name (as given in the topology). Often, this is force
    field dependent. Example: ``name CA`` (for C&alpha; atoms) or ``name
    OW`` (for SPC water oxygen)

type *atom-type*
    select by atom type; this is either a string or a number and depends on
    the force field; it is read from the topology file (e.g. the CHARMM PSF
    file contains numeric atom types). It has non-sensical values when a
    PDB or GRO file is used as a topology.

atom *seg-name*  *residue-number*  *atom-name*
    a selector for a single atom consisting of segid resid atomname,
    e.g. ``DMPC 1 C2`` selects the C2 carbon of the first residue of the
    DMPC segment

altLoc *alternative-location*
    a selection for atoms where alternative locations are available, which is
    often the case with high-resolution crystal structures
    e.g. ``resid 4 and resname ALA and altLoc B`` selects only the atoms of ALA-4
    that have an altLoc B record.

chainID *chain-name*
    a selection for atoms where chainIDs have been defined.

element *element-name*
    a selection for atoms where elements have been defined. e.g. ``element H C``

moltype *molecule-type*
    select by molecule type, e.g. ``moltype Protein_A``. At the moment, only
    the TPR format defines the molecule type.

smarts *SMARTS-query*
    select atoms using Daylight's SMARTS queries, e.g. ``smarts
    [#7;R]`` to find nitrogen atoms in rings. Requires RDKit.
    All matches are combined as a single unique match. The `smarts`
    selection accepts two sets of key word arguments from
    `select_atoms()`: the ``rdkit_kwargs`` are passed internally to
    `RDKitConverter.convert()` and the ``smarts_kwargs`` are passed to
    RDKit's `GetSubstructMatches
    <https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html#rdkit.Chem.rdchem.Mol.GetSubstructMatches>`_.
    By default, the `useChirality` kwarg in ``rdkit_kwargs`` is set to true
    and maxMatches in ``smarts_kwargs`` is ``max(1000, 10 * n_atoms)``, where
    ``n_atoms`` is either ``len(AtomGroup)`` or ``len(Universe.atoms)``,
    whichever is applicable. Note that the number of matches can occasionally
    exceed the default value of maxMatches, causing too few atoms to be
    returned. If this occurs, a warning will be issued. The problem can be
    fixed by increasing the value of maxMatches. This behavior may be updated
    in the future

chiral *R | S*
    select a particular stereocenter. e.g. ``name C and chirality S``
    to select only S-chiral carbon atoms.  Only ``R`` and ``S`` will be
    possible options but other values will not raise an error.

formalcharge *formal-charge*
    select atoms based on their formal charge, e.g.
    ``name O and formalcharge -1`` to select all oxygens with a
    negative 1 formal charge.

Pattern matching
----------------

The pattern matching notation described below is used to specify 
patterns for matching strings (based on :mod:`fnmatch`):

``?`` 
    Is a pattern that will match any single character. For example,
    ``resname T?R`` selects residues named "TYR" and "THR".
``*`` 
    Is a pattern that will match multiple characters.  For example,
    ``GL*`` selects all strings that start with "GL" such as "GLU",
    "GLY", "GLX29", "GLN".
``[seq]``
    Would match any character in seq. For example, "resname GL[NY]" 
    selects all residues named "GLN" or "GLY" but would not select
    "GLU".
``[!seq]``
    Would match any character not in seq. For example, "resname GL[!NY]"
    would match residues named "GLU" but would not match "GLN" or "GLY".

Boolean
-------

not
    all atoms not in the selection, e.g. ``not protein`` selects all atoms
    that aren't part of a protein

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

sphlayer *innerRadius* *externalRadius* *selection*
    selects all atoms within a spherical layer centered in the center of
    geometry (COG) of a given selection, e.g., ``sphlayer 2.4 6.0 ( protein
    and ( resid 130 or resid 80 ) )`` selects the center of geometry of
    protein, resid 130, resid 80 and creates a spherical layer of inner
    radius 2.4 and external radius 6.0 around the COG.

sphzone *externalRadius* *selection*
    selects all atoms within a spherical zone centered in the center of
    geometry (COG) of a given selection, e.g. ``sphzone 6.0 ( protein and (
    resid 130 or resid 80 ) )`` selects the center of geometry of protein,
    resid 130, resid 80 and creates a sphere of radius 6.0 around the COG.

isolayer *inner radius* *outer radius* *selection*
    Similar to sphlayer, but will find layer around all referenced atoms. 
    For example, if the atom types for a polymer backbone were chosen, then
    an isolayer parallel to the backbone will be generated. As another 
    example, if a set of ions were chosen as a reference to isolate the second 
    hydration layer, then they will all be included in the same group.
    However, in the instance that a molecule is in the second hydration layer 
    of one ion and the first hydration layer of another, those atoms will not
    be included.

cylayer *innerRadius* *externalRadius* *zMax* *zMin* *selection*
    selects all atoms within a cylindric layer centered in the center of
    geometry (COG) of a given selection, e.g. ``cylayer 5 10 10 -8
    protein`` selects the center of geometry of protein, and creates a
    cylindrical layer of inner radius 5, external radius 10 centered on the
    COG. In z, the cylinder extends from 10 above the COG to 8
    below. Positive values for *zMin*, or negative ones for *zMax*, are
    allowed.

cyzone *externalRadius* *zMax* *zMin* *selection*
    selects all atoms within a cylindric zone centered in the center of
    geometry (COG) of a given selection, e.g. ``cyzone 15 4 -8 protein and
    resid 42`` selects the center of geometry of protein and resid 42, and
    creates a cylinder of external radius 15 centered on the COG. In z, the
    cylinder extends from 4 above the COG to 8 below. Positive values for
    *zMin*, or negative ones for *zMax*, are allowed.

    .. versionchanged:: 0.10.0
       keywords *cyzone* and *cylayer* now take *zMax* and *zMin* to be
       relative to the COG of *selection*, instead of absolute z-values
       in the box.

point *x* *y* *z*  *distance*
    selects all atoms within a cutoff of a point in space, make sure
    coordinate is separated by spaces, e.g. ``point 5.0 5.0 5.0 3.5``
    selects all atoms within 3.5 Angstroms of the coordinate (5.0, 5.0,
    5.0)

prop [abs] *property*  *operator*  *value*
    selects atoms based on position, using *property* **x**, **y**, or
    **z** coordinate. Supports the **abs** keyword (for absolute value) and
    the following *operators*: **<, >, <=, >=, ==, !=**. For example,
    ``prop z >= 5.0`` selects all atoms with z coordinate greater than 5.0;
    ``prop abs z <= 5.0`` selects all atoms within -5.0 <= z <= 5.0.


.. note::
   By default periodicity **is** taken into account with geometric
   selections, i.e. selections will find atoms that are in different
   periodic images.
   To control this behaviour, use the boolean ``"periodic"`` keyword
   argument of :meth:`~MDAnalysis.core.groups.AtomGroup.select_atoms`.


Similarity and connectivity
---------------------------

same *subkeyword* as *selection*
    selects all atoms that have the same *subkeyword* value as any atom in
    *selection*. Allowed *subkeyword* values are the atom properties: ``name,
    type, resname, resid, segid, mass, charge, radius, bfactor, resnum``, the
    groups an atom belong to: ``residue, segment, fragment``, and the atom
    coordinates ``x, y, z``.

byres *selection*
    selects all atoms that are in the same segment and residue as selection,
    e.g. specify the subselection after the byres keyword.  ``byres`` is a
    shortcut to ``same residue as``

bonded *selection*
    selects all atoms that are bonded to selection
    eg: ``select name H and bonded name O`` selects only hydrogens bonded to
    oxygens

Index
-----

bynum *index-range*
    selects all atoms within a range of (1-based) inclusive indices,
    e.g. ``bynum 1`` selects the first atom in the universe; ``bynum 5:10``
    selects atoms 5 through 10 inclusive. All atoms in the
    :class:`MDAnalysis.Universe` are consecutively numbered, and the index
    runs from 1 up to the total number of atoms.
    
id *index-range*
    selects all atoms in a range of (1-based) inclusive indices, e.g. ``id 1`` selects 
    all the atoms with id 1; ``id 5:7`` selects all atoms with ids 5, all atoms with 
    ids 6 and all atoms with ids 7.
     
index *index-range*
    selects all atoms within a range of (0-based) inclusive indices,
    e.g. ``index 0`` selects the first atom in the universe; ``index 5:10``
    selects atoms 6 through 11 inclusive. All atoms in the
    :class:`MDAnalysis.Universe` are consecutively numbered, and the index
    runs from 0 up to the total number of atoms - 1.
    
    
.. note::
    Conventionally, ``id`` corresponds to the serial number in the PDB format. In contrast 
    to ``bynum``, the ``id`` topology attribute is not necessarily continuous, ordered, or 
    unique, and can be arbitrarily assigned by the user. 
    
     
.. _pre-selections-label:

Preexisting selections and modifiers
------------------------------------

group `group-name`
    selects the atoms in the :class:`AtomGroup` passed to the function as an
    argument named `group-name`. Only the atoms common to `group-name` and the
    instance :meth:`~MDAnalysis.core.groups.AtomGroup.select_atoms` was called
    from will be considered, unless ``group`` is preceded by the ``global``
    keyword. `group-name` will be included in the parsing just by comparison of
    atom indices. This means that it is up to the user to make sure the
    `group-name` group was defined in an appropriate :class:`Universe`.

global *selection*
    by default, when issuing
    :meth:`~MDAnalysis.core.groups.AtomGroup.select_atoms` from an
    :class:`~MDAnalysis.core.groups.AtomGroup`, selections and subselections
    are returned intersected with the atoms of that instance.  Prefixing a
    selection term with ``global`` causes its selection to be returned in its
    entirety.  As an example, the ``global`` keyword allows for
    ``lipids.select_atoms("around 10 global protein")`` --- where ``lipids`` is
    a group that does not contain any proteins. Were ``global`` absent, the
    result would be an empty selection since the ``protein`` subselection would
    itself be empty.  When issuing
    :meth:`~MDAnalysis.core.groups.AtomGroup.select_atoms` from a
    :class:`~MDAnalysis.core.universe.Universe`, ``global`` is ignored.

.. versionchanged:: 1.0.0
   The ``fullgroup`` selection has now been removed. Please use the equivalent
   ``global group`` selection.

Dynamic selections
==================

By default :meth:`~MDAnalysis.core.groups.AtomGroup.select_atoms` returns an
:class:`~MDAnalysis.core.groups.AtomGroup`, in which the list of atoms is
constant across trajectory frame changes. If
:meth:`~MDAnalysis.core.groups.AtomGroup.select_atoms` is invoked with named
argument ``updating`` set to ``True``, an
:class:`~MDAnalysis.core.groups.UpdatingAtomGroup` instance will be returned
instead. It behaves just like an :class:`~MDAnalysis.core.groups.AtomGroup`
object, with the difference that the selection expressions are re-evaluated
every time the trajectory frame changes (this happens lazily, only when the
:class:`~MDAnalysis.core.groups.UpdatingAtomGroup` object is accessed so that
there is no redundant updating going on)::

 # A dynamic selection of corner atoms:
 >>> ag_updating = universe.select_atoms("prop x < 5 and prop y < 5 and prop z < 5", updating=True)
 >>> ag_updating
 <UpdatingAtomGroup with 9 atoms>
 >>> universe.trajectory.next()
 >>> ag_updating
 <UpdatingAtomGroup with 14 atoms>

Using the ``group`` selection keyword for
:ref:`preexisting-selections <pre-selections-label>`, one can
make updating selections depend on
:class:`~MDAnalysis.core.groups.AtomGroup`, or even other
:class:`~MDAnalysis.core.groups.UpdatingAtomGroup`, instances.
Likewise, making an updating selection from an already updating group will
cause later updates to also reflect the updating of the base group::

 >>> chained_ag_updating = ag_updating.select_atoms("resid 1:1000", updating=True)
 >>> chained_ag_updating
 <UpdatingAtomGroup with 3 atoms>
 >>> universe.trajectory.next()
 >>> chained_ag_updating
 <UpdatingAtomGroup with 7 atoms>

Finally, a non-updating selection or a slicing/addition operation made on an
:class:`~MDAnalysis.core.groups.UpdatingAtomGroup` will return a static
:class:`~MDAnalysis.core.groups.AtomGroup`, which will no longer update
across frames::

 >>> static_ag = ag_updating.select_atoms("resid 1:1000")
 >>> static_ag
 <UpdatingAtomGroup with 3 atoms>
 >>> universe.trajectory.next()
 >>> static_ag
 <UpdatingAtomGroup with 3 atoms>

.. _ordered-selections-label:

Ordered selections
==================

:meth:`~MDAnalysis.core.groups.AtomGroup.select_atoms` sorts the atoms
in the :class:`~MDAnalysis.core.groups.AtomGroup` by atom index before
returning them (this is to eliminate possible duplicates in the
selection). If the ordering of atoms is crucial (for instance when
describing angles or dihedrals) or if duplicate atoms are required
then one has to concatenate multiple AtomGroups, which does not sort
them.

The most straightforward way to concatentate two AtomGroups is by using the
``+`` operator::

 >>> ordered = u.select_atoms("segid DMPC and resid 3 and name P") + u.select_atoms("segid DMPC and resid 2 and name P")
 >>> print(list(ordered))
 [< Atom 570: name 'P' of type '180' of resid 'DMPC', 3 and 'DMPC'>,
 < Atom 452: name 'P' of type '180' of resid 'DMPC', 2 and 'DMPC'>]

A shortcut is to provide *two or more* selections to
:meth:`~MDAnalysis.core.universe.Universe.select_atoms`, which then
does the concatenation automatically::

 >>> print(list(universe.select_atoms("segid DMPC and resid 3 and name P", "segid DMPC and resid 2 and name P")))
 [< Atom 570: name 'P' of type '180' of resid 'DMPC', 3 and 'DMPC'>,
 < Atom 452: name 'P' of type '180' of resid 'DMPC', 2 and 'DMPC'>]

Just for comparison to show that a single selection string does not
work as one might expect::

 # WRONG!
 >>> print(list(universe.select_atoms("segid DMPC and (resid 3 or resid 2) and name P")))
 [< Atom 452: name 'P' of type '180' of resid 'DMPC', 2 and 'DMPC'>,
 < Atom 570: name 'P' of type '180' of resid 'DMPC', 3 and 'DMPC'>]
