.. -*- coding: utf-8 -*-
.. _guessing:

====================
Guessing
====================

When a Universe is created from a Universe, MDAnalysis guesses properties that have not been read from the file. Sometimes these properties are available in the file, but are simply not read by MDAnalysis. For example, :ref:`masses are always guessed <guessing-masses>`.

.. _guessing-masses:

Masses
======

Atom masses are always guessed for every file format. They are guessed from the ``Atom.atom_type``. This attribute represents a number of different values in MDAnalysis, depending on which file format you used to create your Universe. ``Atom.atom_type`` can be force-field specific atom types, from files that provide this information; or it can be an element, guessed from the atom name. `See further discussion here. <https://github.com/MDAnalysis/mdanalysis/issues/2348>`_


.. important:: 

    When an atom mass cannot be guessed from the atom ``atom_type`` or ``name``, the atom is assigned a mass of 0.0. Masses are guessed atom-by-atom, so even if most atoms have been guessed correctly, it is possible that some have been given masses of 0. It is important to check for non-zero masses before using methods that rely on them, such as :meth:`AtomGroup.center_of_mass`.


Types
=====

When atom ``atom_type``\ s are guessed, they represent the atom element. Atom types are always guessed from the atom name. MDAnalysis follows biological naming conventions, where atoms named "CA" are much more likely to represent an alpha-carbon than a calcium atom. This guesser is still relatively fragile for non-traditionally biological atom names.

Bonds, Angles, Dihedrals, Impropers
====================================

MDAnalysis can guess if bonds exist between two atoms, based on the distance between them. A bond is created if the 2 atoms are within 

.. math::

    d < f \cdot (R_1 + R_2)

of each other, where :math:`R_1` and :math:`R_2` are the VdW radii
of the atoms and :math:`f` is an ad-hoc *fudge_factor*. This is
the `same algorithm that VMD uses`_.

Angles can be guessed from the bond connectivity. MDAnalysis assumes that if atoms 1 & 2 are bonded, and 2 & 3 are bonded, then (1,2,3) must be an angle.

::

   1        
    \      
     2 -- 3

Dihedral angles and improper dihedrals can both be guessed from angles. Proper dihedrals are guessed by assuming that if (1,2,3) is an angle, and 3 & 4 are bonded, then (1,2,3,4) must be a dihedral.

::

   1        4
    \      /
     2 -- 3

Likewise, if (1,2,3) is an angle, and 2 & 4 are bonded, then (2, 1, 3, 4) must be an improper dihedral (i.e. the improper dihedral is the angle between the planes formed by (1, 2, 3) and (1, 3, 4))

::

   1        
    \      
     2 -- 3
    /
   4

The method available to users is :meth:`AtomGroup.guess_bonds <MDAnalysis.core.groups.AtomGroup.guess_bonds>`, which allows users to pass in a dictionary of van der Waals' radii for atom types. This guesses bonds, angles, and dihedrals (but not impropers) for the specified AtomGroup and adds it to the underlying Universe.


.. _`same algorithm that VMD uses`:
    http://www.ks.uiuc.edu/Research/vmd/vmd-1.9.1/ug/node26.html
    