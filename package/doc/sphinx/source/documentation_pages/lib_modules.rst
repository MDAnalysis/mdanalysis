.. _lib:

*******************************************
Library functions --- :mod:`MDAnalysis.lib`
*******************************************

.. module:: MDAnalysis.lib
   :synopsis: ``lib`` collects independent code for re-use in MDAnalysis

:mod:`MDAnalysis.lib` contains code that is independent of the
specific MDAnalysis framework, such as fast calculators for distances
or simple logging support. Modules do not depend on other code inside
MDAnalysis except in :mod:`MDAnalysis.lib` itself (and possibly in
:mod:`MDAnalysis.exceptions`) and thus can be easily imported
elsewhere.

Overview
--------

:mod:`MDAnalysis.lib.distances` contains many high performance maths
functions. Most of them have the keyword *backend* that allows one to
either select serial (single threaded) code (``backend="serial``) or
to use parallelized versions (e.g. ``backend="OpenMP"`` for OpenMP
parallelism).

:mod:`MDAnalysis.lib.transformations` contains a multitude of
matrix operations for manipulating coordinate data.

:mod:`MDAnalysis.lib.qcprot` contains a fast implementation of
superposition by minimizing the RMSD.

:mod:`MDAnalysis.lib.util` contains various file and string utility
functions whereas mathematical functions are to be found in
:mod:`MDAnalysis.lib.mdamath`.

:mod:`MDAnalysis.lib.NeighborSearch` contains classes to do neighbor
searches with MDAnalysis objects.


List of modules
---------------

.. toctree::
   :maxdepth: 1

   ./lib/distances
   ./lib/NeighborSearch
   ./lib/log
   ./lib/mdamath
   ./lib/transformations
   ./lib/qcprot
   ./lib/util
