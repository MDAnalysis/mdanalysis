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
:class:`~MDAnalysis.core.AtomGroup.Universe` to force the format.

Programmers and anyone trying to implement new functionality should read the
:ref:`Trajectory API` (described in :mod:`MDAnalysis.coordinates`).

.. rubric:: Contents

.. toctree::
   :maxdepth: 1

   coordinates/init
   coordinates/base
   coordinates/core
   coordinates/CRD
   coordinates/DCD
   coordinates/DMS
   coordinates/LAMMPS
   coordinates/GRO
   coordinates/PDB
   coordinates/PDBQT	
   coordinates/PQR
   coordinates/TRJ
   coordinates/TRR
   coordinates/XTC
   coordinates/XYZ
   coordinates/pdbextensions
   coordinates/libxdrfile
