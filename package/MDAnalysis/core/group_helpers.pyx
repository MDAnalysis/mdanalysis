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
Compiled helpers for group iteration --- :mod:`MDAnalysis.core.group_helpers`
=============================================================================

Helpers
"""

from .groups import AtomGroup
from ..libmda cimport iterators
from libcpp.vector cimport vector
from libc.stdint cimport int64_t
cimport cython
cimport numpy as cnp
import numpy as np
cnp.import_array()


cdef inline _to_numpy_from_spec(object owner, int ndim, cnp.npy_intp * shape, int npy_type, void * pointer):
    array = cnp.PyArray_SimpleNewFromData(ndim, shape, npy_type, pointer)
    cnp.PyArray_SetBaseObject(array, owner)
    cnp.Py_INCREF(owner)
    return array

cdef class AtomGroupIterator:

    def __cinit__(self, ag):
        self._iterator = iterators._AtomGroupIterator(ag.n_atoms)
        self._coord_view = ag.universe.trajectory.ts.positions
        self._iterator.ptr = &self._coord_view[0,0]
        self._iterator.copy_ix( < int64_t*>cnp.PyArray_DATA(ag.ix_array))
    
    def print_coords(self):
        print(np.asarray(self._coord_view))

    @property
    def ix(self):
        cdef cnp.npy_intp dims[1]
        dims[0] = self._iterator.n_atoms
        cdef cnp.ndarray arr 
        arr =  _to_numpy_from_spec(self, 1, dims, cnp.NPY_INT64, self._iterator.ix.data())
        return arr

    @property
    def n_atoms(self):
        return self._iterator.n_atoms
    
    def print_ix(self):
        self._iterator.print_ix()


cdef class ArrayIterator:

    def __cinit__(self, cnp.ndarray arr):
        if arr.shape[1] != 3:
            raise ValueError(
                "input array has incorrect second dimension, must be 3")
        self._iterator = iterators._ArrayIterator(arr.shape[0])
        self._coord_view = arr
        self._iterator.ptr = &self._coord_view[0,0]

    def print_coords(self):
        print(np.asarray(self._coord_view))


    @property
    def n_atoms(self):
        return self._iterator.n_atoms
