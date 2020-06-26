.. -*- coding: utf-8 -*-
.. _standard-selections:

==========================================
Standard residues in MDAnalysis selections
==========================================

.. _protein-selection:

Proteins
========
The residue names listed here are accessible via the "protein" keyword in the :ref:`selections`. 

The below names are drawn from the CHARMM 27, OPLS-AA, GROMOS 53A6, AMBER 03, and AMBER 99sb*-ILDN force fields.

.. include:: generated/selections/protein.txt

----------------
Protein backbone
----------------

Protein backbone atoms in MDAnalysis belong to a recognised protein residue and have the atom names:

.. include:: generated/selections/protein_backbone.txt

.. _nucleic-selection:

Nucleic acids
=============

The residue names listed here are accessible via the "nucleic" keyword in the :ref:`selections`. 

The below names are drawn from largely from the CHARMM force field.

.. include:: generated/selections/nucleic.txt

----------------
Nucleic backbone
----------------

Nucleic backbone atoms in MDAnalysis belong to a recognised nucleic acid residue and have the atom names:

.. include:: generated/selections/nucleic_backbone.txt

.. _nucleobase-selection:

-----------
Nucleobases
-----------

Nucleobase atoms from nucleic acid residues are recognised based on their names in CHARMM.

.. include:: generated/selections/base.txt

--------------
Nucleic sugars
--------------

Nucleic sugar atoms from nucleic acid residues are recognised by MDAnalysis if they have the atom names:

.. include:: generated/selections/nucleic_sugar.txt