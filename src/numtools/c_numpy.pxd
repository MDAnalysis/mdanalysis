cdef extern from "Python.h":
                ctypedef int Py_intptr_t

cdef extern from "numpy/arrayobject.h":
                ctypedef class numpy.ndarray [object PyArrayObject]:
                                cdef char *data
                                cdef int nd
                                cdef Py_intptr_t *dimensions
                                cdef Py_intptr_t *strides
                                cdef object base
                                # descr not implemented yet here...
                                cdef int flags
                                cdef int itemsize
                                cdef object weakreflist

                cdef void import_array()
