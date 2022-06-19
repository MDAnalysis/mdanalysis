# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
#


"""
Compiled helpers for group iteration --- :mod:`MDAnalysis.core.group_iterators`
===============================================================================

This module contains Cython extension types that allow various low level
routines to use groups such as an :class:`~MDAnalysis.core.groups.AtomGroup`
or a NumPy :class:`numpy.ndarray` interchangeably. In general this module is
primarily of interest to developers and should be modified with caution.

The Cython extension types in this module are lightweight wrappers around the 
C++ classes defined in `MDAnalysis/lib/include/iterators.h`. These C++ classes
provide a homogenous interface to lower level functions for fast iteration
of groups and arrays and are exposed as the `_iterator` attribute.
See the documentation in `MDAnalysis/lib/include/iterators.h` for more
information.

The classes in this module are primarily responsible for hooking any data
required for iteration onto the `._iterator` C++ class instance which provides
the iteration interface.

Classes
-------

.. autoclass:: AtomGroupIterator
   :members:

.. autoclass:: ArrayIterator
   :members:

"""

import numpy as np
from .groups import AtomGroup
from ..lib cimport iterators
from libcpp.vector cimport vector
from libc.stdint cimport int64_t
cimport cython
cimport numpy as cnp
cnp.import_array()


cdef inline _to_numpy_from_spec(object owner, int ndim, cnp.npy_intp * shape, int npy_type, void * pointer):
    array = cnp.PyArray_SimpleNewFromData(ndim, shape, npy_type, pointer)
    cnp.PyArray_SetBaseObject(array, owner)
    cnp.Py_INCREF(owner)
    return array

cdef class AtomGroupIterator:
    """Iterator for an :class:`~MDAnalysis.core.groups.AtomGroup`

    .. attribute:: _iterator

      :class:`_AtomGroupIterator` C++ class instance from
      `MDAnalysis.lib.include.iterators.h`

    .. versionadded:: 2.3.0
    """

    def __cinit__(self, ag):
        """
        Parameters
        ----------
        ag : :class:`~MDAnalysis.core.groups.AtomGroup`
             atomgroup to be iterated
        """
        self._iterator = iterators._AtomGroupIterator(ag.n_atoms)
        # hook iterator pointer onto base array
        self._iterator.ptr = <float*>cnp.PyArray_DATA(ag.universe.trajectory.ts.positions)
        self._iterator.copy_ix(<int64_t*>cnp.PyArray_DATA(ag.ix_array))

    @property
    def ix(self):
        """
        Read only view of global indices of the
        :class:`~MDAnalysis.core.groups.AtomGroup` being iterated.

        .. versionadded:: 2.3.0
        """
        cdef cnp.npy_intp dims[1]
        dims[0] = self._iterator.n_atoms
        cdef cnp.ndarray arr
        arr = _to_numpy_from_spec(
            self, 1, dims, cnp.NPY_INT64, self._iterator.ix.data())
        return arr

    @property
    def n_atoms(self):
        """
        Read only view of the number of atoms of the
        :class:`~MDAnalysis.core.groups.AtomGroup` being iterated.

        .. versionadded:: 2.3.0
        """
        return self._iterator.n_atoms


cdef class ArrayIterator:
    """Iterator for a :class:`numpy.ndarray`

    .. attribute:: _iterator

      :class:`_ArrayIterator` C++ class instance from
      `MDAnalysis.lib.include.iterators.h`

    .. versionadded:: 2.3.0
    """

    def __cinit__(self, cnp.ndarray arr):
        """
        Parameters
        ----------
        arr : :class:`numpy.ndarry`
             array to be iterated
        """
        self._iterator = iterators._ArrayIterator(arr.shape[0])
        # hook iterator pointer onto base array
        self._iterator.ptr = <float*>cnp.PyArray_DATA(arr)

    @property
    def n_atoms(self):
        """
        Read only view of the number of atoms (size of first dimension) of
        the :class:`numpy.ndarray` being iterated.

        .. versionadded:: 2.3.0
        """
        return self._iterator.n_atoms
