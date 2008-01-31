
cdef extern from "Numeric/arrayobject.h":
    cdef struct PyArray_Descr:
        int type_num, elsize
        char type

    ctypedef class Numeric.ArrayType [object PyArrayObject]:
        cdef char *data
        cdef int nd
        cdef int *dimensions, *strides
        cdef object base
        cdef PyArray_Descr *descr
        cdef int flags

    enum ArrayTypes "PyArray_TYPES":
        PyArray_CHAR, PyArray_UBYTE, PyArray_SBYTE,
        PyArray_SHORT, PyArray_USHORT,
        PyArray_INT, PyArray_UINT,
        PyArray_LONG,
        PyArray_FLOAT, PyArray_DOUBLE,
        PyArray_CFLOAT, PyArray_CDOUBLE,
        PyArray_OBJECT,
        PyArray_NTYPES, PyArray_NOTYPE

cdef extern from "string.h":
    void* memcpy(void *dst, void *src, int len)

cdef extern from "distances.h":
    ctypedef float coordinate[3]

    void calc_distance_array(coordinate* ref, int numref, coordinate* conf, int numconf, float* box, double* distances)
    void calc_self_distance_array(coordinate* ref, int numref, float* box, double* distances, int distnum)


import Numeric as N

def distance_array(ArrayType ref, ArrayType conf, ArrayType box, ArrayType result=None):
    cdef ArrayType distances
    cdef int confnum, refnum

    if (conf.nd != 2 and conf.dimensions[1] != 3):
        raise Exception("conf must be a sequence of 3 dimensional coordinates")
    if (ref.nd != 2 and ref.dimensions[1] != 3):
        raise Exception("ref must be a sequence of 3 dimensional coordinates")
    if (conf.descr.type_num != PyArray_FLOAT or ref.descr.type_num != PyArray_FLOAT):
        raise Exception("coordinate data must be of type Float32")
    if (box.nd != 1 and box.dimensions[0] != 3):
        raise Exception("box must be a sequence of 3 dimensional coordinates")
    if (box.descr.type_num != PyArray_FLOAT):
        raise Exception("periodic boundaries must be of type Float32")

    confnum = conf.dimensions[0]
    refnum = ref.dimensions[0]

    if not result is None:
        if (result.nd != 2 and result.dimensions[0] != refnum and results.dimensions[1] != confnum):
            raise Exception("result array has incorrect size or datatype - should be (%dx%d)"%(refnum,confnum))
        distances = N.asarray(result)
    else:
        distances = N.zeros((refnum, confnum), N.Float64)

    calc_distance_array(<coordinate*>ref.data, refnum, <coordinate*>conf.data, confnum, <float*>box.data, <double*>distances.data)
    return distances

def self_distance_array(ArrayType ref, ArrayType box, ArrayType result = None):
    cdef ArrayType distances
    cdef int refnum, distnum

    if (ref.nd != 2 and ref.dimensions[1] != 3):
        raise Exception("ref must be a sequence of 3 dimensional coordinates")
    if (ref.descr.type_num != PyArray_FLOAT):
        raise Exception("coordinate data must be of type Float32")
    if (box.nd != 1 and box.dimensions[0] != 3):
        raise Exception("box must be a sequence of 3 dimensional coordinates")
    if (box.descr.type_num != PyArray_FLOAT):
        raise Exception("periodic boundaries must be of type Float32")

    refnum = ref.dimensions[0]
    distnum = (refnum*(refnum-1))/2

    if not result is None:
        if (result.nd != 1 and result.dimensions[0] != distnum):
            raise Exception("result array has incorrect size or datatype - should be (%d)"%(distnum))
        distances = N.asarray(result)
    else:
        distances = N.zeros((distnum,), N.Float64)
    calc_self_distance_array(<coordinate*>ref.data, refnum, <float*>box.data, <double*>distances.data, distnum)
    return distances
