.. -*- coding: utf-8 -*-

======================================
 Working with residues (and segments)
======================================

In the previous section we introduced AtomGroups and in principle that
is enough to handle almost all situations. However, it is convenient
to refer to collections of atoms in way that relates to their
chemistry. In particular for proteins there is a notion of a *hierarchy
of groupings* or *levels*. 

.. TODO:: Describe the membership rules, each atom belongs to one and only one
          Residue, same with Res to Segment

Residue and Segment attributes
==============================

* TODO: Make it clear that some attributes must be accessed via further
  shortcuts, ie Residue.atoms.masses
* TODO:  Make it clear that some attributes 'belong' to a level,
  ie Atom.resid belongs only to the Residue and can't fall
  out of sync.

Residue and Segment methods
===========================

TODO

Moving between Residues and Segments
====================================

TODO
