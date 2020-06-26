.. -*- coding: utf-8 -*-
.. _PDB-format:

=======================================
PDB, ENT (Standard PDB file)
=======================================

.. include:: classes/PDB.txt

Reading in
==========

MDAnalysis parses the following *PDB records* (see `PDB coordinate section`_ for
details):

    - **CRYST1** for unit cell dimensions A,B,C, alpha,beta,gamma
    - **ATOM** or **HETATM** for serial, name, resName, chainID, resSeq, x, y, z, occupancy, tempFactor, segID
    - **CONECT** records for bonds
    - **HEADER** (:attr:`Universe.trajectory.header`)
    - **TITLE** (:attr:`Universe.trajectory.title`)
    - **COMPND** (:attr:`Universe.trajectory.compound`)
    - **REMARK** (:attr:`Universe.trajectory.remarks`)

All other lines are ignored. Multi-`MODEL`_ PDB files are read as trajectories with a default timestep of 1 ps :ref:`(pass in the dt argument to change this) <universe-kwargs>`. Currently, MDAnalysis `cannot read multi-model PDB files written by VMD`_, as VMD uses the keyword "END" to separate models instead of "MODEL"/"ENDMDL" keywords. 

.. _`cannot read multi-model PDB files written by VMD`: https://github.com/MDAnalysis/mdanalysis/issues/1133

.. important:: 

    MDAnalysis does not read atom elements or charges from a PDB file, even when they are provided. Instead, elements are guessed from atom names.

MDAnalysis attempts to read ``segid`` attributes from the *segID* column. If this column does not contain information, segments are instead created from chainIDs. If chainIDs are also not present, then ``segid``\ s are set to the default ``'SYSTEM'`` value.


.. _PDB-formatted:
    http://www.wwpdb.org/documentation/file-format-content/format32/v3.2.html
.. _PDB coordinate section:
    http://www.wwpdb.org/documentation/file-format-content/format32/sect9.html
.. _MODEL:
    http://www.wwpdb.org/documentation/file-format-content/format32/sect9.html#MODEL


Writing out
===========

MDAnalysis can write both single-frame PDBs and convert trajectories to multi-model PDBs. If the Universe is missing fields that are :ref:`required in a PDB file <pdb-spec>`, MDAnalysis provides default values and raises a warning. There are 2 exceptions to this:

    - ``chainIDs``: if a Universe does not have ``chainIDs``, MDAnalysis uses the first character of the segment ``segid`` instead. 
    - ``elements``: Elements are *always* guessed from the atom name.

These are the default values:

    * names: 'X'
    * altLocs: ''
    * resnames: 'UNK'
    * icodes: ''
    * segids: ''
    * resids: 1
    * occupancies: 1.0
    * tempfactors: 0.0
    

.. _pdb-spec:

PDB specification
====================

.. table:: CRYST1 fields

    =============  ============  ===========  =============================================
    COLUMNS        DATA  TYPE    FIELD        DEFINITION
    =============  ============  ===========  =============================================
    1 -  6         Record name   "CRYST1"
    7 - 15         Real(9.3)     a              a (Angstroms).
    16 - 24        Real(9.3)     b              b (Angstroms).
    25 - 33        Real(9.3)     c              c (Angstroms).
    34 - 40        Real(7.2)     alpha          alpha (degrees).
    41 - 47        Real(7.2)     beta           beta (degrees).
    48 - 54        Real(7.2)     gamma          gamma (degrees).
    =============  ============  ===========  =============================================

.. table:: ATOM/HETATM fields

    =============  ============  ===========  =============================================
    COLUMNS        DATA  TYPE    FIELD        DEFINITION
    =============  ============  ===========  =============================================
    1 -  6         Record name   "ATOM  "
    7 - 11         Integer       serial       Atom  serial number.
    13 - 16        Atom          name         Atom name.
    17             Character     altLoc       Alternate location indicator.
    18 - 21        Residue name  resName      Residue name.
    22             Character     chainID      Chain identifier.
    23 - 26        Integer       resSeq       Residue sequence number.
    27             AChar         iCode        Code for insertion of residues.
    31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
    39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
    47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
    55 - 60        Real(6.2)     occupancy    Occupancy.
    61 - 66        Real(6.2)     tempFactor   Temperature  factor.
    67 - 76        String        segID        (unofficial CHARMM extension ?)
    77 - 78        LString(2)    element      Element symbol, right-justified.
    79 - 80        LString(2)    charge       Charge  on the atom.
    =============  ============  ===========  =============================================
