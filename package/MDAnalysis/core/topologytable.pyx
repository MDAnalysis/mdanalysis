# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
#
from libcpp.vector cimport vector
from libcpp.map cimport map as cmap
from libcpp.unordered_set cimport unordered_set
from libcpp.pair cimport pair as cpair

from ..lib._cutil import unique_int_1d
from ..lib._cutil cimport to_numpy_from_spec
import numpy as np
cimport numpy as cnp
cnp.import_array()


cdef class TopologyTable:
    def __cinit__(self,  int nval,  int npair, int[:,:] val, cnp.intp_t[:] flat_vals,  int[:] typ, int[:] guess, int[:] order,  **kwargs):
        """Initialise C++ level parameters of a TopologyTable

        Parameters
        ----------
        order : unsigned int
          The order of the TopologyTable

        .. versionadded:: 2.3.0
           Initialise C++ level parameters
        """
    
        # c++ level objects
        self._nval = nval
        self._npair = npair
        self._types.reserve(self._nval)
        self._guessed.reserve(self._nval)
        self._order.reserve(self._nval)


        self._unique = unique_int_1d(flat_vals) # can we cimport this
        self._nunique = self._unique.shape[0]

        # construct tables
        self._construct_empty_tables_and_maps()
        self._copy_types_guessed_order(typ, guess, order)
        self._parse(val, typ, guess, order)
        self._canonicalise_all()
    
    cdef _construct_empty_tables_and_maps(self):
        cdef int i
        cdef vector[int] tmp
        tmp.reserve(8)

        for i in range(self._nunique):
            self._values.push_back(tmp)
            self.vmap_fwd[self._unique[i]] = i
            self.vmap_rev[i] = self._unique[i]
    
    cdef _copy_types_guessed_order(self, int[:] typ, int[:] guess, int[:] order):
        cdef int i 
        for i in range(self._nval):
            self._types.push_back(typ[i])
            self._guessed.push_back(guess[i])
            self._order.push_back(order[i])
    
    cdef _parse(self, int[:,:] val,  int[:] typ, int[:] guess, int[:] order):
        cdef int i
        for i in range(self._nval):
            self._values[self.vmap_fwd[val[i,0]]].push_back(val[i,1])
            self._values[self.vmap_fwd[val[i,1]]].push_back(val[i,0])
            self._types[self.vmap_fwd[val[i,0]]].push_back(typ[i])
            self._types[self.vmap_fwd[val[i,1]]].push_back(typ[i])  
            self._guessed[self.vmap_fwd[val[i,0]]].push_back(guess[i])
            self._guessed[self.vmap_fwd[val[i,1]]].push_back(guess[i])  
            self._order[self.vmap_fwd[val[i,0]]].push_back(order[i])
            self._order[self.vmap_fwd[val[i,1]]].push_back(order[i])  

    cdef  _canonicalise_all(self):
        cdef int i
        for i in range(self._nunique):
            self.values[i] = self._canonicalise_vec(self._values[i])


        

    cdef vector[int] _canonicalise_vec(self, vector[int] inp):
        cdef unordered_set[int] tmp_set
        cdef int size = inp.size()
        cdef int i
        cdef vector[int] unique_ix
        cdef cpair[vector[int], vector[int]] p
        for i in range(size):
            tmp_set.insert(inp[i])
            unique_ix.push_back(i)
        inp.assign(tmp_set.begin(), tmp_set.end())
        # unique_unique_ix(inp, unique_ix)
        return inp
    
    def query_table(self, int atom):
        cdef int idx
        idx = self.vmap_fwd[atom]
        return self._values[idx]
    
    def add_bonds(self, int[:,:] val,  int[:] typ, int[:] guess, int[:] ord):
        cdef int _nvals
        _nvals = val.shape[0]
        for v in range(_nvals):
            if self.vmap_fwd.find(v) == self.vmap_fwd.end():
                # new entry
                pass
            else:
                #append 
                pass
        pass

    def del_bonds(self, int[:,:] val):
        pass


    def unique_bonds(self):
        a = set()
        # cdef int i, idx
        for i in range(self._nunique):
            for j in range(self._values[i].size()):
                a.add((self.vmap_rev[i], self._values[i][j]))
        return list(a)

    def types(self):
        cdef cnp.npy_intp dims[1]
        dims[0] = self._nval
        cdef cnp.ndarray arr
        arr =  to_numpy_from_spec(self, 1, dims, cnp.NPY_INT32, self._types.data())
        return arr
          
    def guessed(self):
        cdef cnp.npy_intp dims[1]
        dims[0] = self._nval
        cdef cnp.ndarray arr
        arr = to_numpy_from_spec(self, 1, dims, cnp.NPY_INT32, self._guessed.data())
        return arr

    def order(self):
        cdef cnp.npy_intp dims[1]
        dims[0] = self._nval
        cdef cnp.ndarray arr
        arr =  to_numpy_from_spec(self, 1, dims, cnp.NPY_INT32, self._order.data())
        return arr
    
        
    def print_values(self):
        for i in range(self._nunique):
            print(f" {i} {self.vmap_rev[i]} {self._values[i]}")

    def print_types(self):
        for i in range(self._nval):
            print(self._types[i])

    def print_guessed(self):
        for i in range(self._nval):
            print(self._guessed[i])

    def print_order(self):
        for i in range(self._nval):
            print(self._order[i])



    

