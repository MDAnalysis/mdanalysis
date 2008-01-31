
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


import Numeric as N

def distance_array(ArrayType ref, ArrayType conf, ArrayType box):
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

    distances = N.zeros((refnum, confnum), N.Float)
    calc_distance_array(<coordinate*>ref.data, refnum, <coordinate*>conf.data, confnum, <float*>box.data, <double*>distances.data)
    return distances
