.. -*- coding: utf-8 -*-
.. _DATA-format:

=========================================
DATA (LAMMPS)
=========================================

.. include:: classes/DATA.txt

.. important::

    Lennard-Jones units are not implemented. See :ref:`units`
    for other recognized values and the documentation for the LAMMPS
    `units command`_.

.. _units command: http://lammps.sandia.gov/doc/units.html

Reading in
==========

Lammps atoms can have `lots of different formats, and even custom formats`_. By default, MDAnalysis checks:

    * "full" : atoms with 7 fields (reading id, resid, type, and charge)
    * "molecular": atoms with 6 fields (reading id, resid, and type)

.. _`lots of different formats, and even custom formats`: http://lammps.sandia.gov/doc/atom_style.html

Users can pass in their own ``atom_style`` specifications.

    * Required fields: id, type, x, y, z
    * Optional fields: resid, charge

For example::

    u = mda.Universe(LAMMPSDATA, atom_style="id resid type charge element bfactor occupancy x y z")

Only id, resid, charge, type, and coordinate information will be read from the file, even if other topology attributes are specified in the ``atom_style`` argument.

Writing out
===========

MDAnalysis supports writing out the header and applicable sections from Atoms, Masses, Velocities, Bonds, Angles, Dihedrals, and Impropers. The Atoms section is written in the "full" sub-style if charges are available or "molecular" sub-style if they are not. The molecule id is set to 0 for all atoms.


This writer assumes "conventional" or "real" LAMMPS units where length is measured in Angstroms and velocity is measured in Angstroms per femtosecond. To write in different units, specify ``lengthunit`` or ``timeunit``.

For example, to write a certain frame with nanometer units::

    >>> for ts in u.trajectory:
    ...     # analyze frame
    ...     if take_this_frame == True:
    ...         with mda.Writer('frame.data') as W:
    ...             W.write(u.atoms, lengthunit="nm")
    ...         break


If atom types are not already positive integers, the user must set them
to be positive integers, because the writer will not automatically
assign new types.

To preserve numerical atom types when writing a selection, the Masses
section will have entries for each atom type up to the maximum atom type.
If the universe does not contain atoms of some type in
``{1, ... max(atom_types)}``, then the mass for that type will be set to 1.

In order to write bonds, each selected bond type must be explicitly set to
an integer >= 1.
