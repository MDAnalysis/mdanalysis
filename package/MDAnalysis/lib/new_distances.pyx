import cython
import numpy as np

from numpy.math import sqrtf as sqrt
from ._cutil cimport norm2
from pbc cimport PBC, minimum_image


@cython.boundscheck(False)
@cython.wraparound(False)
cdef void inner_distance_array(const float[:, :] a, const float[:, :] b,
                               PBC box, float[:, :] result):
    """C level of distance_array"""
    cdef int i, j
    cdef float[3] dx

    for i in range(a.shape[0]):
        for j in range(b.shape[0]):
            dx[0] = a[i, 0] - b[j, 0]
            dx[1] = a[i, 1] - b[j, 1]
            dx[2] = a[i, 2] - b[j, 2]

            minimum_image(dx, box)

            result[i, j] = sqrt(norm2(dx))


def distance_array(coords1, coords2, box):
    cdef PBC cbox
    cdef float[:, :] a, b
    cdef float[:, :] result_view

    a = np.asarray(coords1, dtype=np.float32)
    b = np.asarray(coords2, dtype=np.float32)
    box = np.asarray(box, dtype=np.float32)

    cbox = PBC(box)
    result = np.zeros((a.shape[0], b.shape[0]), dtype=np.float32)
    result_view = result[:]

    inner_distance_array(a, b, cbox, result_view)

    return result
