.. -*- coding: utf-8 -*-
.. _CRD-format:

=======================================
CRD (CHARMM CARD files)
=======================================

.. include:: classes/CRD.txt

Reading in
==========

Read a list of atoms from a CHARMM standard or extended CARD coordinate file (CRD_)
to build a basic topology.  Reads atom ids (ATOMNO), atom names (TYPES),
resids (RESID), residue numbers (RESNO), residue names (RESNames), segment ids
(SEGID) and tempfactor (Weighting).  Atom element and mass are guessed based on
the name of the atom.

.. _CRD: https://www.charmmtutorial.org/index.php/CHARMM:The_Basics

Writing out
===========

MDAnalysis automatically writes the CHARMM EXT extended format if there are more than 99,999 atoms.

Writing a CRD file format requires the following attributes to be present:

    - resids
    - resnames
    - names
    - chainIDs
    - tempfactors

If these are not present, then :ref:`default values <topologyattr-defaults>` are provided and a warning is raised. 