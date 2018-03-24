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

.. toctree::
   :maxdepth: 1

   coordinates/init
   coordinates/CRD
   coordinates/DCD
   coordinates/DLPoly
   coordinates/DMS
   coordinates/GMS
   coordinates/GSD
   coordinates/GRO
   coordinates/INPCRD
   coordinates/LAMMPS
   coordinates/MMTF
   coordinates/MOL2
   coordinates/PDB
   coordinates/PDBQT
   coordinates/PQR
   coordinates/TRJ
   coordinates/TRR
   coordinates/TXYZ
   coordinates/XTC
   coordinates/XYZ
   coordinates/TRZ
   coordinates/memory
   coordinates/null
   coordinates/dummy

.. rubric:: Coordinate core modules

The remaining pages are primarily of interest to
developers. Programmers and anyone trying to implement new
functionality should first read the :ref:`Trajectory API`.


.. toctree::
   :maxdepth: 1

   coordinates/base
   coordinates/core
   coordinates/chain
   coordinates/XDR
