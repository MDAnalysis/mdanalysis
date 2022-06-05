from libcpp cimport bool
from libc.stdint cimport uint64_t
cimport numpy as cnp
cnp.import_array()


cdef class Timestep:

    cdef uint64_t _n_atoms
    cdef public int64_t  frame


    # info for numpy C API
    cdef int _typenum
    cdef cnp.npy_intp _particle_dependent_dim[2]

    cdef bool _has_positions
    cdef bool _has_velocities
    cdef bool _has_forces

    # these have to be public for testing
    cdef public cnp.ndarray _unitcell
    cdef public cnp.ndarray _pos
    cdef public cnp.ndarray _velocities
    cdef public cnp.ndarray _forces

    # python objects
    cdef object _dtype
    cdef public dict data
    cdef public object  _reader
    cdef public object aux
