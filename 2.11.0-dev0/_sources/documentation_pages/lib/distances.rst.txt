.. automodule:: MDAnalysis.lib.distances



Low-level modules for :mod:`MDAnalysis.lib.distances`
=====================================================

:mod:`MDAnalysis.lib._distances` contains low level access to the
*serial* MDAnalysis Cython functions in `distances`.  These have
little to no error checking done on inputs so should be used with
caution. Similarly, :mod:`MDAnalysis.lib._distances_openmp` contains
the OpenMP-enable functions.

.. toctree::
   :maxdepth: 1

   c_distances
   c_distances_openmp
