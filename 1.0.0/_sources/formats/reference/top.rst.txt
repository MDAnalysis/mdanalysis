.. -*- coding: utf-8 -*-
.. _TOP-format:

========================================
TOP, PRMTOP, PARM7 (AMBER topology)
========================================

.. include:: classes/TOP.txt

AMBER specification
===================

.. note::

   The Amber charge is converted to electron charges as used in
   MDAnalysis and other packages. To get back Amber charges, multiply
   by 18.2223.

.. table:: Attributes parsed from AMBER keywords

    +----------------------------+----------------------+
    | AMBER flag                 | MDAnalysis attribute |
    +----------------------------+----------------------+
    | ATOM_NAME                  | names                |
    +----------------------------+----------------------+
    | CHARGE                     | charges              |
    +----------------------------+----------------------+
    | ATOMIC_NUMBER              | elements             |
    +----------------------------+----------------------+
    | MASS                       | masses               |
    +----------------------------+----------------------+
    | BONDS_INC_HYDROGEN         | bonds                |
    | BONDS_WITHOUT_HYDROGEN     |                      |
    +----------------------------+----------------------+
    | ANGLES_INC_HYDROGEN        | angles               |
    | ANGLES_WITHOUT_HYDROGEN    |                      |
    +----------------------------+----------------------+
    | DIHEDRALS_INC_HYDROGEN     | dihedrals / improper |
    | DIHEDRALS_WITHOUT_HYDROGEN |                      |
    +----------------------------+----------------------+
    | ATOM_TYPE_INDEX            | type_indices         |
    +----------------------------+----------------------+
    | AMBER_ATOM_TYPE            | types                |
    +----------------------------+----------------------+
    | RESIDUE_LABEL              | resnames             |
    +----------------------------+----------------------+
    | RESIDUE_POINTER            | residues             |
    +----------------------------+----------------------+


Developer notes
===============

The format is defined in `PARM parameter/topology file
specification`_.  The reader tries to detect if it is a newer
(AMBER 12?) file format by looking for the flag "ATOMIC_NUMBER".

.. _`PARM parameter/topology file specification`: http://ambermd.org/formats.html#topology

