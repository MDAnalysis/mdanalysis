# caffinityprop.pxd --- pxd file for affinityprop.pyx
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

cdef extern from *:
    ctypedef char const_char "const char"

cdef extern from "stdio.h":
    #extern int printf (__const char *__restrict __format, ...);
    int printf(const_char *, ...)

cdef extern from "stdlib.h":
    void *calloc (size_t COUNT, size_t ELTSIZE)

cdef extern from "float.h":
    enum:  FLT_MAX

cdef extern from "ap.h":
    int trmIndex(int, int)
    int sqmIndex(int, int, int)
    float pwmax(float, float)
    float pwmin(float, float)
    float min(float*, int)
    float max(float*, int)
    int CAffinityPropagation(double*, int, double, int, int, bint, long*)
