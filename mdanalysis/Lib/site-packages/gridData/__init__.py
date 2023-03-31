# gridDataFormats --- python modules to read and write gridded data
# Copyright (c) 2009-2014 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Lesser General Public License, version 3 or later.
# See the files COPYING and COPYING.LESSER for details.

r"""
Handling grids of data --- :mod:`gridData`
==========================================

Overview
--------

This module contains classes that allow importing and exporting of
simple gridded data, A grid is an N-dimensional array that represents
a discrete mesh over a region of space. The array axes are taken to be
parallel to the cartesian axes of this space. Together with this array
we also store the edges, which are are (essentially) the cartesian
coordinates of the intersections of the grid (mesh) lines on the
axes. In this way the grid is anchored in space.

The :class:`~gridData.core.Grid` object can be resampled at arbitrary
resolution (by interpolating the data). Standard algebraic operations
are defined for grids on a point-wise basis (same as for
:class:`numpy.ndarray`).


Description
-----------

The package reads grid data from files, makes them available as a
:class:`~gridData.core.Grid` object, and allows one to write out the data again.

A :class:`~gridData.core.Grid` consists of a rectangular, regular, N-dimensional
array of data. It contains

(1) The position of the array cell edges.
(2) The array data itself.

This is equivalent to knowing

(1) The origin of the coordinate system (i.e. which data cell
    corresponds to (0,0,...,0)
(2) The spacing of the grid in each dimension.
(3) The data on a grid.

:class:`~gridData.core.Grid` objects have some convenient properties:

* The data is represented as a :class:`numpy.ndarray` in
  :attr:`Grid.grid<~gridData.core.Grid.grid>` and thus can be directly
  manipulated with all the tools available in NumPy.

* :class:`Grid` instances can be manipulated arithmetically, e.g. one
  can simply add or subtract two of them and get another one, or
  multiply by a constant. Note that all operations are defined
  point-wise (see the :mod:`numpy` documentation for details) and that
  only grids defined on the same cell edges can be combined.

* A :class:`~gridData.core.Grid` object can also be created from
  within python code e.g. from the output of the
  :func:`numpy.histogramdd` function.

* The representation of the data is abstracted from the format that
  the files are saved in. This makes it straightforward to add
  additional readers for new formats.

* The data can be written out again in formats that are understood by
  other programs such as VMD_, ChimeraX_ or PyMOL_.


Reading grid data files
-----------------------

Some Formats_ can be read directly from a file on disk::

 g = Grid(filename)

*filename* could be, for instance, "density.dx".


Constructing a Grid
-------------------

Data from an n-dimensional array can be packaged as a :class:`~gridData.core.Grid`
for convenient handling (especially export to other formats).  The
:class:`~gridData.core.Grid` class acts as a universal constructor::

 g = Grid(ndarray, edges=edges)                 # from histogramdd
 g = Grid(ndarray, origin=origin, delta=delta)  # from arbitrary data

 g.export(filename, format)   # export to the desire format

See the doc string for :class:`~gridData.core.Grid` for details.


Formats
-------

For the available file formats see :ref:`supported-file-formats`.


.. _VMD: https://www.ks.uiuc.edu/Research/vmd/

.. _PyMOL: https://pymol.org/

.. _ChimeraX: https://www.cgl.ucsf.edu/chimerax/

"""

from .core import Grid
from . import OpenDX
from . import gOpenMol
from . import mrc

__all__ = ['Grid', 'OpenDX', 'gOpenMol', 'mrc']

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

