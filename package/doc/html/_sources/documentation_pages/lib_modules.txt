.. _lib:

************************
MDAnalysis library
************************

:mod:`MDAnalysis.lib._distances` contains low level access to
MDAnalysis' math functions.  These have little to
no error checking done on inputs so should be used with caution.

:mod:`MDAnalysis.lib.transformations` contains a multitude of
matrix operations for manipulating coordinate data.

:mod:`MDAnalysis.lib.util` contains various file, string
and mathematical utility functions.

.. toctree::
   :maxdepth: 1

   ./lib/_distances
   ./lib/KDTree_modules
   ./lib/log
   ./lib/transformations
   ./lib/util
