# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
# Beckstein and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
import numpy as np
cimport numpy as np

from libc.stdlib cimport free
from cpython cimport PyObject, Py_INCREF

np.import_array()


cdef class ArrayWrapper:
    """Helper class to construct an ndarray from a raw C-Pointer. This class
    will take ownership of the past memory and free it up on destruction. So
    make sure that all further access to that memory location happens over
    this object. To deallocate the array `free` will be used, so please make
    sure that the memory was allocated with malloc/realloc.

    See Also
    --------
    ptr_to_ndarray

    Notes
    -----
    Original from Gael Varoquaux,
    https://gist.github.com/GaelVaroquaux/1249305#file-cython_wrapper-pyx
    """
    cdef void* data_ptr
    cdef int* dim
    cdef int ndim
    cdef int data_type

    cdef set_data(self, void* data_ptr, int* dim, int ndim, int data_type):
        """ Set the data of the array
        This cannot be done in the constructor as it must recieve C-level
        arguments.

        Parameters:
        -----------
        data_ptr: void*
            Pointer to the data
        dim: int*
            length in each dimension
        ndim: int
            number of dimensions
        data_type: int
            Numpy DataType enum
        """
        self.data_ptr = data_ptr
        self.dim = dim
        self.data_type = data_type
        self.ndim = ndim

    def __array__(self):
        """ Here we use the __array__ method, that is called when numpy
            tries to get an array from the object."""
        ndarray = np.PyArray_SimpleNewFromData(self.ndim,
                                               <np.npy_intp*> self.dim,
                                               self.data_type,
                                               self.data_ptr)
        return ndarray

    def __dealloc__(self):
        """ Frees the array. This is called by Python when all the
        references to the object are gone. """
        free(<void*>self.data_ptr)
        # free(<int*>self.dim)


cdef np.ndarray ptr_to_ndarray(void* data_ptr, np.int64_t[:] dim, int data_type):
    """convert a pointer to an arbitrary C-pointer to a ndarray. The ndarray is
    constructed so that the array it's holding will be freed when the array is
    destroyed.

    Parameters
    ----------
    data_ptr : void*
        Pointer to the data
    dim : int[:]
        array containing length in each dimension
    data_type : int
        Numpy DataType enum

    Returns
    -------
    ndarray
        Numpy array containing the data

    """
    array_wrapper = ArrayWrapper()
    array_wrapper.set_data(<void*> data_ptr, <int*> &dim[0], dim.size, data_type)

    cdef np.ndarray ndarray = np.array(array_wrapper, copy=False)
    # Assign our object to the 'base' of the ndarray object
    ndarray.base = <PyObject*> array_wrapper
    # Increment the reference count, as the above assignement was done in
    # C, and Python does not know that there is this additional reference
    Py_INCREF(array_wrapper)

    return ndarray
