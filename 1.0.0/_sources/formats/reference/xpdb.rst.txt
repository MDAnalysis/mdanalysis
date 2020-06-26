.. -*- coding: utf-8 -*-
.. _XPDB-format:

=======================================
XPDB (Extended PDB file)
=======================================

.. include:: classes/XPDB.txt

The extended PDB reader acts virtually the same as the :ref:`PDB-format` reader. The difference is that extended PDB files (MDAnalysis format specifier *XPDB*) may contain residue sequence numbers up to 99,999 by utilizing the insertion character field of the PDB standard. Five-digit residue numbers may take up columns 23 to 27 (inclusive) instead of being confined to 23-26 (with column 27 being reserved for the insertion code in the PDB standard). 

PDB files in this format are written
by popular programs such as VMD_.


As extended PDB files are very similar to PDB files, tell MDAnalysis to use the Extended PDB parser by passing in the ``topology_format`` keyword.

.. ipython:: python

    import MDAnalysis as mda
    from MDAnalysis.tests.datafiles import PDB
    pdb = mda.Universe(PDB)
    pdb.trajectory.format
    xpdb = mda.Universe(PDB, topology_format='XPDB')
    xpdb.trajectory.format

.. _VMD: http://www.ks.uiuc.edu/Research/vmd/