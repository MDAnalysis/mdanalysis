cimport numpy as np

cdef float _dot(float *, float *)
cdef void _cross(float *, float *, float *)
cdef float _norm(float *)
cdef  to_numpy_from_spec(object owner, int ndim, np.npy_intp * shape, int npy_type, void * pointer)