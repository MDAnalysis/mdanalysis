.. -*- coding: utf-8 -*-
.. _TXYZ-format:

========================================
TXYZ, ARC (Tinker)
========================================

.. include:: classes/TXYZ.txt

MDAnalysis can read Tinker_ xyz files .txyz and trajectory .arc files.

Developer notes
===============

Differences between Tinker format_ and normal xyz files:

- there is only one header line containing both the number of atoms and a comment
- column 1 contains atom numbers (starting from 1)
- column 6 contains atoms types
- the following columns indicate connectivity (atoms to which that particular atom is
  bonded, according to numbering in column 1)

.. _format: http://chembytes.wikidot.com/tnk-tut00#toc2
.. _Tinker: https://dasher.wustl.edu/tinker/
