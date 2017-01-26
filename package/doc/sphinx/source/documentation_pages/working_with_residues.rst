.. -*- coding: utf-8 -*-

.. AtomGroups have been introduced, now introduce the idea of a
   hierarchy of different levels

Describe the membership rules, each atom belongs to one and only one
Residue, same with Res to Segment

Residue and Segment attributes

Make it clear that some attributes must be accessed via further
shortcuts, ie Residue.atoms.masses

Make it clear that some attributes 'belong' to a level,
ie Atom.resid belongs only to the Residue and can't fall
out of sync.

Residue and Segment methods

Moving between Residues and Segments
