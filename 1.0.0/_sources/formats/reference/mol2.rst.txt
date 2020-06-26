.. -*- coding: utf-8 -*-
.. _MOL2-format:

=======================================
MOL2 (Tripos structure)
=======================================

.. include:: classes/MOL2.txt

The Tripos_ molecule structure format (MOL2_) is a commonly used format. It is used, for instance, by the DOCK_ docking code.

.. _MOL2: http://chemyang.ccnu.edu.cn/ccb/server/AIMMS/mol2.pdf
.. _Tripos: http://www.tripos.com/
.. _DOCK: http://dock.compbio.ucsf.edu/


.. warning::

    :class:`~MDAnalysis.coordinates.MOL2.MOL2Writer` can only be used to write out previously loaded MOL2 files.
    For example, if you're trying to convert a PDB file to MOL2, you should
    use other tools such as rdkit_.

    Here is an example how to use rdkit_ to convert a PDB to MOL::

      from rdkit import Chem
      mol = Chem.MolFromPDBFile("molecule.pdb", removeHs=False)
      Chem.MolToMolFile(mol, "molecule.mol" )

    .. _rdkit: http://www.rdkit.org/docs/GettingStartedInPython.html



MOL2 specification
==================

* Example file::

    #    Name: benzene
    #    Creating user name: tom
    #    Creation time: Wed Dec 28 00:18:30 1988

    #    Modifying user name: tom
    #    Modification time: Wed Dec 28 00:18:30 1988

    @<TRIPOS>MOLECULE
    benzene
    12 12 1  0   0
    SMALL
    NO_CHARGES


    @<TRIPOS>ATOM
    1   C1  1.207   2.091   0.000   C.ar    1   BENZENE 0.000
    2   C2  2.414   1.394   0.000   C.ar    1   BENZENE 0.000
    3   C3  2.414   0.000   0.000   C.ar    1   BENZENE 0.000
    4   C4  1.207   -0.697  0.000   C.ar    1   BENZENE 0.000
    5   C5  0.000   0.000   0.000   C.ar    1   BENZENE 0.000
    6   C6  0.000   1.394   0.000   C.ar    1   BENZENE 0.000
    7   H1  1.207   3.175   0.000   H   1   BENZENE 0.000
    8   H2  3.353   1.936   0.000   H   1   BENZENE 0.000
    9   H3  3.353   -0.542  0.000   H   1   BENZENE 0.000
    10  H4  1.207   -1.781  0.000   H   1   BENZENE 0.000
    11  H5  -0.939  -0.542  0.000   H   1   BENZENE 0.000
    12  H6  -0.939  1.936   0.000   H   1   BENZENE 0.000
    @<TRIPOS>BOND
    1   1   2   ar
    2   1   6   ar
    3   2   3   ar
    4   3   4   ar
    5   4   5   ar
    6   5   6   ar
    7   1   7   1
    8   2   8   1
    9   3   9   1
    10  4   10  1
    11  5   11  1
    12  6   12  1
   @<TRIPOS>SUBSTRUCTURE
    1   BENZENE 1   PERM    0   ****    ****    0   ROOT