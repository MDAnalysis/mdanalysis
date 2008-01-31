
cimport c_numpy
c_numpy.import_array()

cdef extern from "string.h":
    void* memcpy(void *dst, void *src, int len)

cdef extern from "calc_distances.h":
    ctypedef float coordinate[3]

    void calc_distance_array(coordinate* ref, int numref, coordinate* conf, int numconf, float* box, double* distances)
    void calc_self_distance_array(coordinate* ref, int numref, float* box, double* distances, int distnum)


import numpy
def distance_array(c_numpy.ndarray ref, c_numpy.ndarray conf, c_numpy.ndarray box, c_numpy.ndarray result=None):
    cdef c_numpy.ndarray distances
    cdef int confnum, refnum

    if (conf.nd != 2 and conf.dimensions[1] != 3):
        raise Exception("conf must be a sequence of 3 dimensional coordinates")
    if (ref.nd != 2 and ref.dimensions[1] != 3):
        raise Exception("ref must be a sequence of 3 dimensional coordinates")
    if (conf.dtype!=numpy.dtype(numpy.float32) and ref.dtype!=numpy.dtype(numpy.float32)):
        raise Exception("coordinate data must be of type float32")
    if (box.nd != 1 and box.dimensions[0] != 3):
        raise Exception("box must be a sequence of 3 dimensional coordinates")
    if (box.dtype!=numpy.dtype(numpy.float32)):
        raise Exception("periodic boundaries must be of type float32")

    confnum = conf.dimensions[0]
    refnum = ref.dimensions[0]

    if not result is None:
        if (result.nd != 2 and result.dimensions[0] != refnum and results.dimensions[1] != confnum):
            raise Exception("result array has incorrect size or datatype - should be (%dx%d)"%(refnum,confnum))
        distances = numpy.asarray(result)
    else:
        distances = numpy.zeros((refnum, confnum), numpy.float64)

    calc_distance_array(<coordinate*>ref.data, refnum, <coordinate*>conf.data, confnum, <float*>box.data, <double*>distances.data)
    return distances

def self_distance_array(c_numpy.ndarray ref, c_numpy.ndarray box, c_numpy.ndarray result = None):
    cdef c_numpy.ndarray distances
    cdef int refnum, distnum

    if (ref.nd != 2 and ref.dimensions[1] != 3):
        raise Exception("ref must be a sequence of 3 dimensional coordinates")
    if (ref.dtype!=numpy.dtype(numpy.float32)):
        raise Exception("coordinate data must be of type float32")
    if (box.nd != 1 and box.dimensions[0] != 3):
        raise Exception("box must be a sequence of 3 dimensional coordinates")
    if (box.dtype!=numpy.dtype(numpy.float32)):
        raise Exception("periodic boundaries must be of type float32")

    refnum = ref.dimensions[0]
    distnum = (refnum*(refnum-1))/2

    if not result is None:
        if (result.nd != 1 and result.dimensions[0] != distnum):
            raise Exception("result array has incorrect size or datatype - should be (%d)"%(distnum))
        distances = numpy.asarray(result)
    else:
        distances = numpy.zeros((distnum,), numpy.float64)
    calc_self_distance_array(<coordinate*>ref.data, refnum, <float*>box.data, <double*>distances.data, distnum)
    return distances
