# cython: embedsignature=True
# stochasticproxembed.pyx --- Cython wrapper for the stochastic proximity embedding C library
# Copyright (C) 2014 Wouter Boomsma, Matteo Tiberti
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import logging
import numpy
cimport numpy

cimport cstochasticproxembed
cimport cython


@cython.embedsignature(True)

cdef class StochasticProximityEmbedding:
    """
    Stochastic proximity embedding dimensionality reduction algorithm. The algorithm implemented here is described in this paper:

	Dmitrii N. Rassokhin, Dimitris K. Agrafiotis
	A modified update rule for stochastic proximity embedding
	Journal of Molecular Graphics and Modelling 22 (2003) 133–140

    This class is a Cython wrapper for a C implementation (see spe.c)
    """

    def run(self, s, double rco, int dim, double maxlam, double minlam, int ncycle, int nstep, int stressfreq):
        """Run stochastic proximity embedding.

	**Arguments:**
	
	`s` : encore.utils.TriangularMatrix object
                Triangular matrix containing the distance values for each pair of elements in the original space.	

	`rco` : float
		neighborhood distance cut-off

	`dim` : int
		number of dimensions for the embedded space

	`minlam` : float
		final learning parameter

	`maxlam` : float
		starting learning parameter

	`ncycle` : int
		number of cycles. Each cycle is composed of nstep steps. At the end of each cycle, the lerning parameter lambda is updated.

	`nstep` : int
		number of coordinate update steps for each cycle
	
	**Returns:**

	`space` : (float, numpy.array)
		float is the final stress obtained; the array are the coordinates of the elements in the embedded space 
 
	`stressfreq` : int
		calculate and report stress value every stressfreq cycle
	"""

        cdef int nelem = s.size
        cdef double finalstress = 0.0
        
        logging.info("Starting Stochastic Proximity Embedding")

        cdef numpy.ndarray[numpy.float64_t,  ndim=1] matndarray = numpy.ascontiguousarray(s._elements, dtype=numpy.float64)
        cdef numpy.ndarray[numpy.float64_t,   ndim=1] d_coords   = numpy.zeros((nelem*dim),dtype=numpy.float64)
        
        finalstress = cstochasticproxembed.CStochasticProximityEmbedding( <double*>matndarray.data, <double*>d_coords.data, rco, nelem, dim, maxlam, minlam, ncycle, nstep, stressfreq)
        
        logging.info("Stochastic Proximity Embedding finished. Residual stress: %.3f" % finalstress)
          
        return (finalstress, d_coords.reshape((-1,dim)).T)
	
    def __call__(self, *args):
        return self.run(*args)

cdef class kNNStochasticProximityEmbedding:
    """
    k-Nearest Neighbours Stochastic proximity embedding dimensionality reduction algorithm. 
    This is a variation of the SPE algorithm in which neighbourhood is not defined by a distance cut-off; instead, at each step, when a point is randomly chosen to perform coordinate updates, the coordinates of its k nearest neighbours are updated as well.
    This class is a Cython wrapper for a C implementation (see spe.c)
   """ 

    def run(self, s, int kn, int dim, double maxlam, double minlam, int ncycle, int nstep, int stressfreq):
        """Run kNN-SPE.

         **Arguments:**

        `s` : encore.utils.TriangularMatrix object
                Triangular matrix containing the distance values for each pair of elements in the original space.

        `kn` : int
		number of k points to be used as neighbours, in the original space

        `dim` : int
                number of dimensions for the embedded space

        `minlam` : float
                final learning parameter

        `maxlam` : float
                starting learning parameter

        `ncycle` : int
                number of cycles. Each cycle is composed of nstep steps. At the end of each cycle, the lerning parameter lambda is updated.

        `nstep` : int
                number of coordinate update steps for each cycle

        **Returns:**

        `space` : (float, numpy.array)
                float is the final stress obtained; the array are the coordinates of the elements in the embedded space

        `stressfreq` : int
                calculate and report stress value every stressfreq cycle
        """
        
        cdef int nelem = s.size
        cdef double finalstress = 0.0
        
        logging.info("Starting k-Nearest Neighbours Stochastic Proximity Embedding")
        
        cdef numpy.ndarray[numpy.float64_t,  ndim=1] matndarray = numpy.ascontiguousarray(s._elements, dtype=numpy.float64)
        cdef numpy.ndarray[numpy.float64_t,  ndim=1] d_coords   = numpy.zeros((nelem*dim),dtype=numpy.float64)
        
        finalstress = cstochasticproxembed.CkNNStochasticProximityEmbedding(<double*>matndarray.data, <double*>d_coords.data, kn, nelem, dim, maxlam, minlam, ncycle, nstep, stressfreq)
        
        logging.info("Stochastic Proximity Embedding finished. Residual stress: %.3f" % finalstress)
          
        return (finalstress, d_coords.reshape((-1,dim)).T)
