import cython

from cymem.cymem cimport Pool
from libc.math cimport fabs
from numpy.math cimport deg2radf, sinf as sin, cosf as cos, sqrtf as sqrt

from _cutil cimport norm2

DEF MAX_NTRICVEC=12
DEF EPSILON=1e-5
DEF BOX_MARGIN=1.0010

cdef enum PBC_TYPES:
    ORTHO = 1
    TRICLINIC = 2
    NONE = 3


cdef class PBC:
    """Cython class representing periodic boundaries in system

    Used for minimum_image function
    """
    def __cinit__(self, float[:] box):
        self.mem = None
        self.box = NULL
        self.fbox_diag = NULL
        self.hbox_diag = NULL
        self.mhbox_diag = NULL
        self.tric_vec = NULL

    def __init__(self, float[:] box):
        self.mem = Pool()
        self.box = <float**>self.mem.alloc(3, sizeof(float*))
        for i in range(3):
            self.box[i] = <float*>self.mem.alloc(3, sizeof(float))
        self.fbox_diag = <float*>self.mem.alloc(3, sizeof(float))
        self.hbox_diag = <float*>self.mem.alloc(3, sizeof(float))
        self.mhbox_diag = <float*>self.mem.alloc(3, sizeof(float))
        self.tric_vec = <float**>self.mem.alloc(MAX_NTRICVEC, sizeof(float*))
        for i in range(MAX_NTRICVEC):
            self.tric_vec[i] = <float*>self.mem.alloc(3, sizeof(float))
        self.ntric_vec = 0

        self.define_box(box)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void define_box(self, float[:] box):
        # turn box into triclinic vectors
        cdef bint triclinic

        triclinic = False
        for i in range(3):
            if fabs(box[i] - 90.0) > EPSILON:
                triclinic = True

        if triclinic:
            self.pbc_type = TRICLINIC
            self.triclinic_define(box)
        else:
            self.pbc_type = ORTHO
            self.ortho_define(box)

        self.fbox_diag[0] = box[0]
        self.fbox_diag[1] = box[1]
        self.fbox_diag[2] = box[2]

        self.hbox_diag[0] = 0.5 * self.fbox_diag[0]
        self.hbox_diag[1] = 0.5 * self.fbox_diag[1]
        self.hbox_diag[2] = 0.5 * self.fbox_diag[2]

        self.mhbox_diag[0] = - self.hbox_diag[0]
        self.mhbox_diag[1] = - self.hbox_diag[1]
        self.mhbox_diag[2] = - self.hbox_diag[2]

        self.max_cutoff2 = self.calc_max_cutoff2()

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void ortho_define(self, float[:] box):
        self.box[0][0] = box[0]
        self.box[0][1] = 0.0
        self.box[0][2] = 0.0

        self.box[1][0] = 0.0
        self.box[1][1] = box[1]
        self.box[1][2] = 0.0

        self.box[2][0] = 0.0
        self.box[2][1] = 0.0
        self.box[2][2] = box[2]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void triclinic_define(self, float[:] box):
        cdef float alpha, beta, gamma

        alpha = deg2radf(box[3])
        beta = deg2radf(box[4])
        gamma = deg2radf(box[5])

        self.box[0][0] = box[0]
        self.box[0][1] = 0.0
        self.box[0][2] = 0.0

        self.box[1][0] = box[1] * cos(gamma)
        self.box[1][1] = box[1] * sin(gamma)
        self.box[1][2] = 0.0

        self.box[2][0] = box[2] * cos(beta)
        self.box[2][1] = box[2] * (cos(alpha) - cos(beta) * cos(gamma) / sin(gamma))
        self.box[2][2] = box[2] * sqrt(box[2] * box[2] - self.box[2][0] ** 2 - self.box[2][1] ** 2)

    cdef void find_triclinic_shifts(self):
        cdef int[3] order
        cdef int i, j, k, dd, d, shift
        cdef float[3] trial, pos
        cdef float d2old, d2new, d2new_c
        cdef bint use

        order[0] = 0
        order[1] = -1
        order[2] = 1

        for k in order:
            for j in order:
                for i in order:
                    if j != 0 or k != 0:
                        d2old = 0.0
                        d2new = 0.0

                        for d in range(3):
                            trial[d] = i * self.box[0][d] + j * self.box[1][d] + k * self.box[2][d]
                            if d == 3:
                                trial[d] = 0.0
                                pos[d] = 0.0
                            else:
                                if trial[d] < 0:
                                    pos[d] = min(self.hbox_diag[d], -trial[d])
                                else:
                                    pos[d] = max(self.hbox_diag[d], -trial[d])

                            d2old += pos[d] * pos[d]
                            d2new += (pos[d] + trial[d]) * (pos[d] + trial[d])

                        if BOX_MARGIN * d2new < d2old:
                            use = True

                            for dd in range(3):  # for d in [i, j, k]:
                                if dd == 0:
                                    shift = i
                                elif dd == 1:
                                    shift = j
                                else:
                                    shift = k

                                if shift:
                                    d2new_c = 0
                                    for d in range(3):
                                        d2new_c += (pos[d] + trial[d] - shift * self.box[dd][d])**2

                                    if d2new_c <= BOX_MARGIN * d2new:
                                        use = False

                            if use: # Accept this shift vector.
                                for d in range(3):
                                    self.tric_vec[self.ntric_vec][d] = trial[d]
                                    self.ntric_vec += 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef float calc_max_cutoff2(self):
        cdef float min_hv2, min_ss

        min_hv2 = 0.25 * min(norm2(self.box[0]),
                             norm2(self.box[1]),
                             norm2(self.box[2]))

        min_ss = min(self.box[0][0], fabs(self.box[1][1] - self.box[2][1]), self.box[2][2])

        return min(min_hv2, min_ss * min_ss)


cdef void minimum_image(float* dx, PBC pbc):
    """Apply minimum image convention to a vector *dx*

    Parameters
    ----------
    float* dx
      a single 3d separation
    PBC pbc
      PBC object representing the box

    Modifies dx in place
    """
    cdef int i, j
    cdef float d2min, d2trial
    cdef float[3] dx_start, trial

    if pbc.pbc_type == TRICLINIC:
        for i in range(2, -1, -1):
            while (dx[i] > pbc.hbox_diag[i]):
                for j in range(i, -1, -1):
                    dx[j] -= pbc.box[i][j]
            while (dx[i] <= pbc.mhbox_diag[i]):
                for j in range(i, -1, -1):
                    dx[j] += pbc.box[i][j]
        d2min = norm2(dx)
        if d2min > pbc.max_cutoff2:
            for j in range(3):
                dx_start[j] = dx[j]

            i = 0
            while (d2min > pbc.max_cutoff2) and (i < pbc.ntric_vec):
                for j in range(3):
                    trial[j] = dx_start[j] + pbc.tric_vec[i][j]
                d2trial = norm2(trial)
                if d2trial < d2min:
                    for j in range(3):
                        dx[j] = trial[j]
                    d2min = d2trial
                i += 1

    elif pbc.pbc_type == ORTHO:
        for i in range(3):
            while (dx[i] > pbc.hbox_diag[i]):
                dx[i] -= pbc.fbox_diag[i]
            while (dx[i] <= pbc.mhbox_diag[i]):
                dx[i] += pbc.fbox_diag[i]

    elif pbc.pbc_type == NONE:
        pass
