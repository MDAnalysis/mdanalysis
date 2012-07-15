cimport numpy as np
cimport cython
from cython.parallel import parallel, prange

import numpy as np

# Register a np.float32 as data type 'DTYPE_t'
DTYPE = np.float32
ctypedef np.float32_t DTYPE_t

# Register a C math sqrt function
cdef extern from "math.h":
    float sqrt(double x) nogil


def distance_array_serial(np.ndarray[DTYPE_t, ndim=2] coordA, np.ndarray[DTYPE_t, ndim=2]  coordB):    
    cdef DTYPE_t x, y, z, dist
    cdef Py_ssize_t i, j
    cdef np.ndarray[DTYPE_t, ndim=2] distance
    
    # FIXME assume that coorA and coorB are of tha same length
    rows = coordA.shape[0];
    cols = coordB.shape[0];
    distance = np.empty((rows, cols), dtype=DTYPE)
    
    # The two loops are independent, let's use 
    for i in range(rows) :
        for j in range(cols) :
            
            #if i == j:
            #  distance[i,j] = 0.0;
            #  continue
            
            
            x = coordA[i,0] - coordB[j,0];
            y = coordA[i,1] - coordB[j,1];
            z = coordA[i,2] - coordB[j,2];  
            
            dist = sqrt((x*x)+(y*y)+(z*z));

            distance[i,j] = dist;

    return distance
    
    
@cython.boundscheck(False)
def distance_array(np.ndarray[DTYPE_t, ndim=2] coordA, np.ndarray[DTYPE_t, ndim=2]  coordB):    
    cdef DTYPE_t x, y, z, dist
    cdef Py_ssize_t i, j
    cdef np.ndarray[DTYPE_t, ndim=2] distance
    
    # FIXME assume that coorA and coorB are of tha same length
    rows = coordA.shape[0];
    cols = coordB.shape[0];
    distance = np.empty((rows, cols), dtype=DTYPE)
    with nogil, parallel():
        # The two loops are independent, let's use 
        for i in prange(rows, schedule="dynamic", chunksize=50) :
            for j in range(cols) :
                
                #if i == j:
                #  distance[i,j] = 0.0;
                #  continue       
                
                x = coordA[i,0] - coordB[j,0];
                y = coordA[i,1] - coordB[j,1];
                z = coordA[i,2] - coordB[j,2];  
                
                # FIXME this might not be the optimal thing to do
                dist = sqrt((x*x)+(y*y)+(z*z));
                
                distance[i,j] = dist;
    return distance
    
