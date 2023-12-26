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

A number of modules are concerned with finding
neighbors. :mod:`MDAnalysis.lib.NeighborSearch` contains high-level
classes to do neighbor searches with MDAnalysis
objects. :mod:`MDAnalysis.lib.nsgrid` contains a fast implementation
of grid neighbor search whereas :mod:`MDAnalysis.lib.pkdtree` uses
KDTrees (with periodic images) for neighbor searching. Some of the
functions in :mod:`MDAnalysis.lib.distances` user either of these
algorithms to speed up distance calculations.



List of modules
---------------

.. toctree::
   :maxdepth: 1

   ./lib/distances
   ./lib/NeighborSearch
   ./lib/nsgrid
   ./lib/pkdtree	      
   ./lib/log
   ./lib/mdamath
   ./lib/transformations
   ./lib/qcprot
   ./lib/util
   ./lib/correlations
   ./lib/picklable_file_io

Low level file formats
----------------------

The modules in :mod:`MDAnalysis.lib.formats` contain code to access various file
formats in a way that is *independent from other MDAnalysis functionality*
(i.e., they do not use any classes from :mod:`MDAnalysis.core` or
:mod:`MDAnalysis.topology`). This low-level code is used in the
:mod:`MDAnalysis.coordinates` module but can also be re-used by other
Python-based projects.

.. toctree::
   :maxdepth: 1

   ./lib/formats/libmdaxdr
   ./lib/formats/libdcd

Libmdanalysis
-------------

The :file:`__init__.pxd` file in :mod:`MDAnalysis.lib.libmdanalysis` provides a
single place to ``cimport`` MDAnalysis' public Cython headers. This is recommended
for advanced developers only.

For example, imagine we are writing a Cython extension module in
:mod:`MDAnalysis.lib` and we would like to make a function that creates a 
:class:`MDAnalysis.coordinates.timestep.Timestep`

.. code-block:: cython

   from MDAnalysis.lib.libmdanalysis cimport timestep
   # or we could use the relative cimport
   # from .libmdanalysis cimport timestep
   
   cdef timestep.Timestep make_timestep(int natoms):
      cdef timestep.Timestep ts = timestep.Timestep(natoms)
      return ts


Currently modules that are exposed as public Cython headers are:

- :mod:`MDAnalysis.coordinates.timestep`
- :mod:`MDAnalysis.lib.formats.libmdaxdr`
- :mod:`MDAnalysis.lib.formats.libdcd`
- :mod:`MDAnalysis.lib.c_distances`

For more details consult the source :mod:`MDAnalysis.lib.libmdanalysis.__init__.pxd`