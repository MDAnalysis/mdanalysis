.. _lib:

************************
MDAnalysis library
************************

:mod:`MDAnalysis.lib.distances` contains many high performance
maths functions.

:mod:`MDAnalysis.lib._distances` contains low level access to the
MDAnalysis' Cython functions in `distances`.  These have little to no
error checking done on inputs so should be used with
caution. :mod:`MDAnalysis.lib.parallel` contains parallelized
implementations for some these algorithms.

:mod:`MDAnalysis.lib.transformations` contains a multitude of
matrix operations for manipulating coordinate data.

:mod:`MDAnalysis.lib.qcprot` contains a fast implementation of
superposition by minimizing the RMSD.

:mod:`MDAnalysis.lib.util` contains various file, string
and mathematical utility functions.

:mod:`MDAnalysis.lib.NeighborSearch` contains classes to do neighbor searches
     with MDAnalysis objects.


.. toctree::
   :maxdepth: 1

   ./lib/distances
   ./lib/_distances
   ./lib/NeighborSearch
   ./lib/log
   ./lib/mdamath
   ./lib/parallel
   ./lib/transformations
   ./lib/qcprot
   ./lib/util
