.. -*- coding: utf-8 -*-
.. _DMS-format:

=======================================
DMS (Desmond Molecular Structure files)
=======================================

.. include:: classes/DMS.txt

The DESRES Molecular Structure (DMS) file is an SQLite-format database for storing coordinate and topology information. See the `Desmond Users Guide`_ (chapter 6 and chapter 17) for more information.

.. important:: **Atom ids**

    Unlike most other file formats, Desmond starts atom numbers at 0. This means the first atom in a DMS file will have an :code:`Atom.id` of 0. However, residues are not necessarily numbered from 0. :code:`Residue.resid` numbering can start from 1.

.. _`Desmond Users Guide` : http://www.deshawresearch.com/Desmond_Users_Guide-0.7.pdf