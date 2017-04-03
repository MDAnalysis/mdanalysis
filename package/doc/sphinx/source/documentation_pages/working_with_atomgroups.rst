.. -*- coding: utf-8 -*-
.. All about how to use AtomGroups
   very basic as this is one of the first topics in
   the quasi tutorial

Working with AtomGroups
#######################


Manipulating AtomGroups
=======================


Slicing
-------
 * just like numpy
 * fancy indexing
 * boolean indexing too, useful for filtering

Combining AtomGroups
--------------------

 * concatenation
 * subtraction

Set methods
-----------

 * and
 * or
 * xor
 * unique


Using AtomGroups
================

Common attributes
-----------------

 * positions, velocities and forces (of currently loaded frame!)
 * topology information, such as names or masses
 * dimensions
 * residues and segments
 * bonds angles torsions

Common methods
--------------

 * center_of_geometry
 * select_atoms
 * write
