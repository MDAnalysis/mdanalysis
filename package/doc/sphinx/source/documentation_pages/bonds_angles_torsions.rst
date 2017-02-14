.. -*- coding: utf-8 -*-

.. working with bonds angles and torsions

Working with bonds, angles and torsions
=======================================

Chemical bonds between atoms can be manipulated using MDAnalysis in an
object oriented fashion similar to AtomGroups.  MDAnalysis supports bonds,
angles, dihedrals (also known as torsions) and improper dihedrals.
These are made available as attributes of an AtomGroup
::

     >>> import MDAnalysis as mda
     >>> from MDAnalysisTests.datafiles import PSF, DCD
     >>> u = mda.Universe(PSF, DCD)
     >>> u.atoms.bonds
     <TopologyGroup containing 3365 bonds>
     >>> u.atoms.angles
     <TopologyGroup containing 6123 angles>
     >>> u.atoms.dihedrals
     <TopologyGroup containing 8921 dihedrals>
     >>> u.atoms.impropers
     <TopologyGroup containing 541 impropers>


Guessing bonds
--------------

If bonding information was not included in the topology file, it is also
possible to guess the bonds present in the system based upon the distances
between atoms.
::

     >>> import MDAnalysis as mda
     >>> u = mda.Universe('myfile.gro')
     >>> u.atoms.guess_bonds(vdwradii={'C': 1.61, 'N':1.64}

For large systems this will become very slow as it searches all pairwise
combinations.  To circumvent this you can instead only search within a
smaller group of atoms, for example looping through the segments in a
system.
::

     for segment in u.segments:
	 segment.atoms.guess_bonds()


TopologyGroups
--------------

Regardless of the type of bond that is being dealt with, the main working
object for handling them is the TopologyGroup.  These are designed to
slice and add identically to AtomGroups, so you can use them in much the
same way.
::

     >>> import MDAnalysis as mda
     >>> from MDAnalysisTests.datafiles import PSF, DCD
     >>> u = mda.Universe(PSF, DCD)
     >>> u.atoms.bonds[:10]
     <TopologyGroup containing 10 bonds>
     >>> u.atoms.angles[[5, 8, 9]] + u.atoms.angles[[1, 2, 3]]
     <TopologyGroup containing 6 angles>


The values, bond lengths or angle sizes, can be calculated using the ``.values``
method.  This uses a Cython backend so should be relatively fast.
Note that for angles, results will be returned in units of radians.
::

     >>> import MDAnalysis as mda
     >>> from MDAnalysisTests.datafiles import waterPSF, waterDCD
     >>> import numpy as np
     >>> u = mda.Universe(waterPSF, waterDCD)
     >>> ag = u.select_atoms('type OT')
     >>> np.rad2deg(ag.angles.values())
     array([ 104.52003132,  104.52005773,  104.52001203,  104.52000557,
             104.51993961])

The values method has the ``pbc`` kwarg which will account for periodic
boundaries when performing these calculations.  This is important when
the molecules you are analysing may have been 'chopped' and packed into
the primary unit cell.


Selecting bonds
---------------

Similar to AtomGroups, TopologyGroups have a ``select_bonds`` method which
allows for bonds to be selected based on a tuples of their types.
For example working with a box of ethanol::

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

There are only ``C-C-H`` bonds and no ``H-C-C`` bonds.  Selection however is
aware that sometimes types are reversed::

    >>> u.angles.select_bonds(('H', 'C', 'C'))  # note reversal of type
    <TopologyGroup containing 7365 Angles>

Many different types can be selected by combining different TopologyGroups
::

    >>> tg = u.angles.select_bonds(('C', 'C', 'O')) + u.angles.select_bonds(('C', 'O', 'H'))
    >>> tg.types()
    [('C', 'O', 'H'), ('C', 'C', 'O')]
