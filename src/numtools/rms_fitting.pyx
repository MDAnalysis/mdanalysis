
cdef extern from "Numeric/arrayobject.h":
    struct PyArray_Descr:
        int type_num, elsize
        char type

    ctypedef class Numeric.ArrayType [object PyArrayObject]:
        cdef char *data
        cdef int nd
        cdef int *dimensions, *strides
        cdef object base
        cdef PyArray_Descr *descr
        cdef int flags

cdef extern from "mkl_lapack.h":
    int dsyev_(char *jobz, char *uplo, int *n, double *fa, int *lda, double *w, double *work, int *lwork, int *info)

cdef extern from "string.h":
    void* memcpy(void *dst, void *src, int len)

cdef eigen4x4sym(double* mat, double* eval, double* evec):
    cdef double work[264]
    cdef int n, lwork, info

    # Copy mat data to evec since dsyev places the eigenvectors in the original matrix
    memcpy(<void*>evec, <void*>mat, sizeof(double)*16)

    n = 4
    lwork = 264
    info = -1

    dsyev_("V", "U", &n, evec, &n, eval, work, &lwork, &info)
    return info

import Numeric
def rms_rotation_matrix(ArrayType conf, ArrayType ref, ArrayType weights):
    cdef ArrayType ref_cms, pos, cross, k, mat
    cdef ArrayType e, v
    cdef int i, numatoms
    cdef double possq

    if conf.dimensions[0] != ref.dimensions[0]:
        raise Exception("Error: RMS fit - conformation and reference don't have the same number of atoms")
    numatoms = conf.dimensions[0]
    ref_cms = Numeric.zeros((3,), Numeric.Float64)
    for i from 0 <= i < numatoms:
        ref_cms = ref_cms + ref[i]
    pos = Numeric.zeros((3,), Numeric.Float64)
    possq = 0.
    cross = Numeric.zeros((3,3), Numeric.Float64)
    for i from 0 <= i < numatoms:
        r = conf[i]
        r_ref = ref[i]-ref_cms
        w = weights[i]
        pos = pos + w*r
        possq = possq + w*Numeric.add.reduce(r*r) + w*Numeric.add.reduce(r_ref*r_ref)
        cross = cross + w*r[:, Numeric.NewAxis]*r_ref[Numeric.NewAxis, :]
    k = Numeric.zeros((4,4), Numeric.Float64)
    k[0, 0] = -cross[0, 0]-cross[1, 1]-cross[2, 2]
    k[0, 1] = cross[1, 2]-cross[2, 1]
    k[0, 2] = cross[2, 0]-cross[0, 2]
    k[0, 3] = cross[0, 1]-cross[1, 0]
    k[1, 1] = -cross[0, 0]+cross[1, 1]+cross[2, 2]
    k[1, 2] = -cross[0, 1]-cross[1, 0]
    k[1, 3] = -cross[0, 2]-cross[2, 0]
    k[2, 2] = cross[0, 0]-cross[1, 1]+cross[2, 2]
    k[2, 3] = -cross[1, 2]-cross[2, 1]
    k[3, 3] = cross[0, 0]+cross[1, 1]-cross[2, 2]
    for i from 1 <= i < 4:
        for j from 0 <= j < i:
            k[i, j] = k[j, i]
    k = 2.*k
    for i from 0 <= i < 4:
        k[i, i] = k[i, i] + possq - Numeric.add.reduce(pos*pos)
    e = Numeric.zeros((4,), Numeric.Float64)
    v = Numeric.zeros((4,4), Numeric.Float64)
    i = eigen4x4sym(<double*>k.data, <double*>e.data, <double*>v.data)
    v = v[0]
    # Quaternion def: v = qw + iqx + jqy + kqz
    mat = Numeric.zeros((3,3), Numeric.Float64)
    mat[0,0] = 1 - 2*v[2]*v[2] - 2*v[3]*v[3]
    mat[0,1] = 2*v[1]*v[2] - 2*v[0]*v[3]
    mat[0,2] = 2*v[1]*v[3] + 2*v[0]*v[2]
    mat[1,0] = 2*v[1]*v[2] + 2*v[0]*v[3]
    mat[1,1] = 1 - 2*v[1]*v[1] - 2*v[3]*v[3]
    mat[1,2] = 2*v[2]*v[3] - 2*v[0]*v[1]
    mat[2,0] = 2*v[1]*v[3] - 2*v[0]*v[2]
    mat[2,1] = 2*v[2]*v[3] + 2*v[0]*v[1]
    mat[2,2] = 1 - 2*v[1]*v[1] - 2*v[2]*v[2]
    return mat
