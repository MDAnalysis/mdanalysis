# cutils.pyx --- C-compiled Utility functions for encore package
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

#cython embedsignature=True
import numpy as np
cimport numpy as np
import cython

cdef extern from "math.h":
    double sqrt(double x)
    

@cython.boundscheck(False)
@cython.wraparound(False)

def PureRMSD(np.ndarray[np.float64_t,ndim=2] coordsi,
             np.ndarray[np.float64_t,ndim=2] coordsj,
             int atomsn,
             np.ndarray[np.float64_t,ndim=1] masses,
             double summasses):

    cdef  int k
    cdef double normsum, totmasses

    normsum = 0.0

    for k in xrange(atomsn):
        normsum += masses[k]*((coordsi[k,0]-coordsj[k,0])**2 + (coordsi[k,1]-coordsj[k,1])**2 + (coordsi[k,2]-coordsj[k,2])**2)
    return sqrt(normsum/summasses)

def MinusRMSD(np.ndarray[np.float64_t,ndim=2] coordsi,
             np.ndarray[np.float64_t,ndim=2] coordsj,
             int atomsn,
             np.ndarray[np.float64_t,ndim=1] masses,
             double summasses):

    cdef  int k
    cdef double normsum, totmasses

    normsum = 0.0

    for k in xrange(atomsn):
        normsum += masses[k]*((coordsi[k,0]-coordsj[k,0])**2 + (coordsi[k,1]-coordsj[k,1])**2 + (coordsi[k,2]-coordsj[k,2])**2)
    return -sqrt(normsum/summasses)
    

