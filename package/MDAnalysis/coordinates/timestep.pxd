from libcpp cimport bool
from libc.stdint cimport uint64_t, int64_t
cimport numpy as cnp
cnp.import_array()


cdef class Timestep:

    # number of atoms in the timestep
    cdef uint64_t _n_atoms
    # the frame number of the timestep
    cdef public int64_t  frame


    # Info for NumPy C API
    # NumPy enumerated type, must be one of cnp.NPY_FLOAT32 or cnp.NPY_FLOAT64
    cdef int _typenum
    # shape of particle dependent arrays, eg positions, set to (_n_atoms, 3)
    # on __cinit__
    cdef cnp.npy_intp _particle_dependent_dim[2]

    # indicates whether timestep has particle associated data
    cdef bool _has_positions
    cdef bool _has_velocities
    cdef bool _has_forces


    # tracks whether particle associated data has been allocated correct shape
    cdef bool _positions_alloc
    cdef bool _velocities_alloc
    cdef bool _forces_alloc

    # unitcell and particle dependent data
    # these have to be public for testing
    cdef public cnp.ndarray _unitcell
    cdef public cnp.ndarray _pos
    cdef public cnp.ndarray _velocities
    cdef public cnp.ndarray _forces

    # Python objects associated with the timestep
    # NumPy dtype object
    cdef object _dtype
    # additional data associated with timestep
    cdef public dict data
    # weakref back to reader that fills this timestep
    cdef public object  _reader
    # auxilliary namespace
    cdef public object aux
