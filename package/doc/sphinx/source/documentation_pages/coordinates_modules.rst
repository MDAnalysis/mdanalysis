.. Contains the formatted docstrings from the coordinates modules located in 'mdanalysis/MDAnalysis/coordinates'
.. _Coordinates:

**************************
Coordinates modules
**************************

The coordinates module contains the classes to read and write
trajectories. Typically, MDAnalysis recognizes :ref:`Supported coordinate
formats` by the file extension and hence most users probably do not need to
concern themselves with classes and functions described here. However,
if MDAnalysis fails to recognize a coordinate file then the user can
provide the format in the keyword argument *format* to
:class:`~MDAnalysis.core.universe.Universe` to force the format.

.. rubric:: Coordinate formats

Trajectory readers and writers
==============================
.. toctree::
   :maxdepth: 1

   coordinates/init

CRD structure files
===================

.. toctree::
   :maxdepth: 1   

   coordinates/CRD

DCD trajectory I/O
==================

.. toctree::
   :maxdepth: 1  

   coordinates/DCD

DL_Poly fromat reader
=====================

.. toctree::
   :maxdepth: 1  

   coordinates/DLPoly

DESRES file format
==================

.. toctree::
   :maxdepth: 1  

   coordinates/DMS

GAMESS trajectory reader
========================

.. toctree::
   :maxdepth: 1  

   coordinates/GMS

GSD trajectory reader
=====================

.. toctree::
   :maxdepth: 1  

   coordinates/GSD

GRO file format
===============

.. toctree::
   :maxdepth: 1  

   coordinates/GRO

H5MD trajectories
=================

.. toctree::
   :maxdepth: 1    

   coordinates/H5MD

INCPRD structure files
======================

.. toctree::
   :maxdepth: 1  

   coordinates/INPCRD

LAMMPS DCD trajectory and DATA I/O
==================================

.. toctree::
   :maxdepth: 1  

   coordinates/LAMMPS

MMTF trajectory reader
======================

.. toctree::
   :maxdepth: 1  

   coordinates/MMTF

MOL2 file format
================

.. toctree::
   :maxdepth: 1 

   coordinates/MOL2

NAMDBIN files format
====================

.. toctree::
   :maxdepth: 1 

   coordinates/NAMDBIN

PDB structure files
===================

.. toctree::
   :maxdepth: 1 

   coordinates/PDB

PDBQT structure files
=====================

.. toctree::
   :maxdepth: 1 

   coordinates/PDBQT

PQR file format
===============

.. toctree::
   :maxdepth: 1 

   coordinates/PQR

TNG trajectory files
====================

.. toctree::
   :maxdepth: 1

   coordinates/TNG

AMBER trajectories
==================

.. toctree::
   :maxdepth: 1 

   coordinates/TRJ

TRR trajectory files
====================

.. toctree::
   :maxdepth: 1 

   coordinates/TRR

TRZ trajectory I/O
==================

.. toctree::
   :maxdepth: 1 

   coordinates/TRZ

TXYZ file format
================

.. toctree::
   :maxdepth: 1 

   coordinates/TXYZ

XTC trajectory files
====================

.. toctree::
   :maxdepth: 1 

   coordinates/XTC

XYZ trajectory reader
=====================

.. toctree::
   :maxdepth: 1 

   coordinates/XYZ

FHI-AIMS file format
====================

.. toctree::
   :maxdepth: 1 

   coordinates/FHIAIMS

Reading trajectories from memory
================================

.. toctree::
   :maxdepth: 1

   coordinates/memory

Reading trajectories with chemfiles
===================================

.. toctree::
   :maxdepth: 1

   coordinates/chemfiles

Null output
===========

.. toctree::
   :maxdepth: 1

   coordinates/null

.. rubric:: Coordinate core modules

The remaining pages are primarily of interest to
developers. Programmers and anyone trying to implement new
functionality should first read the :ref:`Trajectory API`.


Timestep Class
==============

.. toctree::
   :maxdepth: 1 

   coordinates/timestep

Base Classes
============

.. toctree::
   :maxdepth: 1 

   coordinates/base

Common functions for coordinate reading
=======================================

.. toctree::
   :maxdepth: 1 

   coordinates/core

Serialization of Coordinate Readers
===================================

.. toctree::
   :maxdepth: 1

   coordinates/pickle_readers

ChainReader
===========

.. toctree::
   :maxdepth: 1

   coordinates/chain

XDR based trajectory files
==========================

.. toctree::
   :maxdepth: 1
    
   coordinates/XDR

In particular, all trajectory readers have to be 
:ref:`serializable<serialization>` and they should pass all tests
available in the ``MDAnalysisTests.coordinates.base.MultiframeReaderTest`` 
or ``MDAnalysisTests.coordinates.base.BaseReaderTest``.
