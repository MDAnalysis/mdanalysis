.. -*- coding: utf-8 -*-
.. _PDBQT-format:

=======================================
PDBQT (Autodock structure)
=======================================

.. include:: classes/PDBQT.txt

Reading in
==========

MDAnalysis reads coordinates from PDBQT_ files and additional optional
data such as B-factors, partial charge and AutoDock_ atom types.  It
is also possible to substitute a PDBQT file for a PSF file in order to
define the list of atoms (but no connectivity information will be
available in this case).

Although PDBQT is a similar file format to PDB, MDAnalysis treats them with several differences:

    * Multi-model PDBQT files are not supported
    * Connectivity is not supported (i.e. bonds are not read)

.. _PDBQT:
   http://autodock.scripps.edu/faqs-help/faq/what-is-the-format-of-a-pdbqt-file
.. _AutoDock:
   http://autodock.scripps.edu/


Writing out
===========

MDAnalysis implements a subset of the PDB_ 3.2 standard and the PDBQT_ spec. Unlike the :ref:`PDB-format` writer, MDAnalysis cannot write multi-frame trajectories to a PDBQT file.

If the Universe is missing fields that are :ref:`required in a PDBQT file <pdbqt-spec>`, MDAnalysis provides default values and raises a warning. There are 2 exceptions to this:

    - ``chainIDs``: if a Universe does not have ``chainIDs``, MDAnalysis uses the first character of the segment ``segid`` instead. 
    - ``elements``: MDAnalysis uses the atom type as the element.

These are the default values:

    * names: 'X'
    * altLocs: ''
    * resnames: 'UNK'
    * icodes: ''
    * segids: ''
    * resids: 1
    * occupancies: 1.0
    * tempfactors: 0.0
    * types (elements): ''
    * charges: 0.0

    .. _PDB: http://www.wwpdb.org/documentation/file-format-content/format32/v3.2.html
    .. _PDBQT: http://autodock.scripps.edu/faqs-help/faq/what-is-the-format-of-a-pdbqt-file


.. _pdbqt-spec:

PDBQT specification
===================

Records read:
    - **CRYST1** for unit cell dimensions A,B,C, alpha,beta,gamma
    - **ATOM** or **HETATM** for serial, name, resName, chainID, resSeq, x, y, z, occupancy, tempFactor, segID

.. _PDB format documentation:
    http://www.wwpdb.org/documentation/file-format-content/format32/v3.2.html
.. _AutoDOCK extensions:
    http://autodock.scripps.edu/faqs-help/faq/what-is-the-format-of-a-pdbqt-file

.. table:: Original `PDB format documentation`_ with `AutoDOCK extensions`_

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

    1 -  6         Record name   "ATOM  "
    7 - 11         Integer       serial       Atom  serial number.
    13 - 16        Atom          name         Atom name.
    17             Character     altLoc       Alternate location indicator. IGNORED
    18 - 21        Residue name  resName      Residue name.
    22             Character     chainID      Chain identifier.
    23 - 26        Integer       resSeq       Residue sequence number.
    27             AChar         iCode        Code for insertion of residues. IGNORED
    31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
    39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
    47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
    55 - 60        Real(6.2)     occupancy    Occupancy.
    61 - 66        Real(6.2)     tempFactor   Temperature  factor.
    67 - 76        Real(10.4)    partialChrg  Gasteiger PEOE partial charge *q*.
    79 - 80        LString(2)    atomType     AutoDOCK atom type *t*.
    =============  ============  ===========  =============================================

We ignore torsion notation and just pull the partial charge and atom type columns:

.. code-block::

    COMPND    NSC7810
    REMARK  3 active torsions:
    REMARK  status: ('A' for Active; 'I' for Inactive)
    REMARK    1  A    between atoms: A7_7  and  C22_23
    REMARK    2  A    between atoms: A9_9  and  A11_11
    REMARK    3  A    between atoms: A17_17  and  C21_21
    ROOT
    123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789. (column reference)
    ATOM      1  A1  INH I           1.054   3.021   1.101  0.00  0.00     0.002 A
    ATOM      2  A2  INH I           1.150   1.704   0.764  0.00  0.00     0.012 A
    ATOM      3  A3  INH I          -0.006   0.975   0.431  0.00  0.00    -0.024 A
    ATOM      4  A4  INH I           0.070  -0.385   0.081  0.00  0.00     0.012 A
    ATOM      5  A5  INH I          -1.062  -1.073  -0.238  0.00  0.00     0.002 A
    ATOM      6  A6  INH I          -2.306  -0.456  -0.226  0.00  0.00     0.019 A
    ATOM      7  A7  INH I          -2.426   0.885   0.114  0.00  0.00     0.052 A
    ATOM      8  A8  INH I          -1.265   1.621   0.449  0.00  0.00     0.002 A
    ATOM      9  A9  INH I          -1.339   2.986   0.801  0.00  0.00    -0.013 A
    ATOM     10  A10 INH I          -0.176   3.667   1.128  0.00  0.00     0.013 A
    ENDROOT
    BRANCH   9  11
    ATOM     11  A11 INH I          -2.644   3.682   0.827  0.00  0.00    -0.013 A
    ATOM     12  A16 INH I          -3.007   4.557  -0.220  0.00  0.00     0.002 A
    ATOM     13  A12 INH I          -3.522   3.485   1.882  0.00  0.00     0.013 A
    ATOM     14  A15 INH I          -4.262   5.209  -0.177  0.00  0.00    -0.024 A
    ATOM     15  A17 INH I          -2.144   4.784  -1.319  0.00  0.00     0.052 A
    ATOM     16  A14 INH I          -5.122   4.981   0.910  0.00  0.00     0.012 A
    ATOM     17  A20 INH I          -4.627   6.077  -1.222  0.00  0.00     0.012 A
    ATOM     18  A13 INH I          -4.749   4.135   1.912  0.00  0.00     0.002 A
    ATOM     19  A19 INH I          -3.777   6.285  -2.267  0.00  0.00     0.002 A
    ATOM     20  A18 INH I          -2.543   5.650  -2.328  0.00  0.00     0.019 A
    BRANCH  15  21
    ATOM     21  C21 INH I          -0.834   4.113  -1.388  0.00  0.00     0.210 C
    ATOM     22  O1  INH I          -0.774   2.915  -1.581  0.00  0.00    -0.644 OA
    ATOM     23  O3  INH I           0.298   4.828  -1.237  0.00  0.00    -0.644 OA
    ENDBRANCH  15  21
    ENDBRANCH   9  11
    BRANCH   7  24
    ATOM     24  C22 INH I          -3.749   1.535   0.125  0.00  0.00     0.210 C
    ATOM     25  O2  INH I          -4.019   2.378  -0.708  0.00  0.00    -0.644 OA
    ATOM     26  O4  INH I          -4.659   1.196   1.059  0.00  0.00    -0.644 OA
    ENDBRANCH   7  24
    TORSDOF 3
    123456789.123456789.123456789.123456789.123456789.123456789.123456789.123456789. (column reference)

