import numpy as np
cimport cython
cimport numpy as np


@cython.boundscheck(False)
@cython.wraparound(False)
def residue_positions(
        np.ndarray[np.float32_t, ndim=2] coords,
        np.ndarray[np.int64_t, ndim=1] indices):
    cdef Py_ssize_t i, j, k

    cdef int n_items = indices.shape[0]
    cdef int n_res = np.unique(indices).shape[0]

    cdef np.ndarray[np.int64_t] counter = np.zeros(n_res, dtype=np.int64)
    cdef np.ndarray[np.float32_t, ndim=2] output = np.zeros((n_res, 3),
                                                            dtype=np.float32)

    for i in range(n_items):
        j = indices[i]
        for k in range(3):
            output[j, k] += coords[i, k]
        counter[j] += 1

    for i in range(n_res):
        for k in range(3):
            output[i, k] /= counter[i]

    return output
