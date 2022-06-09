from libcpp.vector cimport vector
from libc.stdint cimport uint64_t

cdef extern from "iterators.h":
    cdef cppclass _AtomGroupIterator:
        uint64_t n_atoms
        vector[uint64_t] ix
        uint64_t i
        float *ptr


cdef extern from "iterators.h":
    cdef cppclass _ArrayIterator:
        int64_t n_atoms
        vector[uint64_t] ix
        uint64_t i
        float *ptr
