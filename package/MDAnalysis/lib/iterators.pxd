from libcpp.vector cimport vector
from libc.stdint cimport int64_t

cdef extern from "iterators.h":
    cdef cppclass _AtomGroupIterator:
        int64_t n_atoms
        int64_t *ix
        int64_t i
        float *ptr
        _AtomGroupIterator()
        _AtomGroupIterator(int64_t n_atoms)
        void print_ix()
        void copy_ix(const int64_t* source )
        void load_into_external_buffer(float *buffer, int64_t n_idx)

cdef extern from "iterators.h":
    cdef cppclass _ArrayIterator:
        int64_t n_atoms
        int64_t i
        float *ptr
        _ArrayIterator()
        _ArrayIterator(int64_t n_atoms)
        void load_into_external_buffer(float *buffer, int64_t n_idx)
