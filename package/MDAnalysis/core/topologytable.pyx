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
from libcpp.set cimport set as cset
from libcpp.pair cimport pair as cpair
from libcpp.algorithm cimport sort as csort
from libcpp.algorithm cimport unique as cunique
from libcpp.string cimport string as cstring
from libcpp.iterator cimport iterator
from libcpp cimport bool as cbool
from ..lib._cutil cimport to_numpy_from_spec
from cython.operator cimport dereference as deref
import numpy as np
cimport numpy as cnp
cnp.import_array()


cdef class TopologyTable:
    def __cinit__(self, int[:,:] val, int n_atoms, list typ,  list guess, list order, **kwargs):
        """Initialise C++ level parameters of a TopologyTable

        Parameters
        ----------

        .. versionadded:: 2.3.0
           Initialise C++ level parameters of a TopologyTable
        """
        self.n_atoms = n_atoms
        
        #typ_ = np.asarray(["-1" if i is None else i for i in typ], dtype=string)
        guess_ = np.asarray([int(i) for i in guess], dtype=np.int32)
        order_ = np.asarray([-1 if i is None else i for i in order], dtype=np.int32)
        typ_ = np.arange(order_.shape[0], dtype=np.int32)
        self._generate_bix(val, typ_,  guess_, order_)



    cdef void _pairsort(self, vector[cpair[int,int]] &a, vector[int] &b):
        cdef vector[cpair[cpair[int,int], int]] pair_arr
        cdef cpair[cpair[int,int], int] pair
        cdef cpair[int,int] ptmp
        cdef int itmp
        for i in range(a.size()):
            ptmp = a[i]
            itmp = b[i]
            pair = cpair[cpair[int,int], int](ptmp, itmp)
            pair_arr.push_back(pair)

        csort(pair_arr.begin(), pair_arr.end())

        for i in range(a.size()):
            a[i] = pair_arr[i].first
            b[i] = pair_arr[i].second
   
    cdef void _generate_bix(self, int[:,:] val, int[:] typ, int[:] guess, int[:] order):
        """Generate bond indicies, spans and value arrays for the TopologyTable
        
        """
        cdef vector[cpair[int, int]] _values
        cdef vector[cpair[int, int]] _values_copy

        cdef vector[int]  _type
        cdef vector[int]  _order
        cdef vector[int]  _guessed

        
        # the bond, forward and reverse
        cdef cpair[int, int] bond
        cdef cpair[int, int] rev_bond

        cdef int i 
        cdef cset[cpair[int, int]] _tmp_set
        for i in range(val.shape[0]):
            bond = cpair[int, int](val[i,0], val[i,1])
            _tmp_set.insert(bond)



        # make sure each bond is in the table both forwards and backwards IF IT
        # ISNT ALREADY
        for i in range(val.shape[0]):
            bond = cpair[int, int](val[i,0], val[i,1])
            rev_bond = cpair[int, int](val[i,1], val[i,0])
            _values.push_back(bond)
            _type.push_back(typ[i])
            _guessed.push_back(guess[i])
            _order.push_back(order[i])

            if _tmp_set.count(rev_bond):
                pass # the reverse bond is already present, do not add
            else:
                _values.push_back(rev_bond)
                _type.push_back(typ[i])
                _order.push_back(order[i])
                _guessed.push_back(guess[i])


        # pairwise sort each array can these be done together all at once ?
        _values_copy = _values
        self._pairsort(_values, _type)
        _values = _values_copy
        self._pairsort(_values, _order)
        _values = _values_copy
        self._pairsort(_values, _guessed)
  
        # deduplicate??

        # initialise spans
        cdef int lead_val
        cdef int prev_val = val[0,0]
        self._span_map[prev_val] = 0

        # unique value counter
        cdef int bix_counter = 0

        for i in range(_values.size()):
            bond = _values[i]
            rev_bond = cpair[int, int](bond.second, bond.first)
            if self._mapping.count(bond):
                # the value is already in the map, grab forward value
                # and that we will read second element
                self._bix.push_back(self._mapping[bond])
            
            elif self._mapping.count(rev_bond):
                # the reversed value is already in the map, grab reverse value
                # and that we will read first element
                self._bix.push_back(self._mapping[rev_bond])

            else:
                # new value
                self._mapping.insert(cpair[cpair[int,int], int](bond, bix_counter))
                
                self._bix.push_back(bix_counter)

                # increment unique values counter
                bix_counter += 1
                # save new value to ix_pair array
                self._ix_pair_array.push_back(bond)
                self._type.push_back(_type[i])
                self._order.push_back(_order[i])
                self._guessed.push_back(_order[i])
            
            lead_val = bond.first
            
            if lead_val != prev_val:
                self._span_map[lead_val] = i
            
            prev_val = lead_val

        self._span_map[prev_val] = i



       # need to find the maximum key and value to correctly cap span array
        cdef int max_val, max_key

        max_key = deref(self._span_map.rbegin()).first
        max_val = deref(self._span_map.rbegin()).second

        # pop fake element into map to correctly cap the span
        self._span_map[max_key +1] = max_val +1  

        # sort out the spans so that each atoms has a span
        prev_val = 0
        for i in range(self.n_atoms):
            if self._span_map.count(i):
                    self._spans.push_back(self._span_map[i])
                    prev_val = self._span_map[i]
            else:
                self._spans.push_back(prev_val) 

        # remove fake element 
        self._span_map.erase(max_key +1)


    def get_bonds(self, int target):
        cdef vector[int] bonds
        bonds = self._get_bond(target)
        return bonds
    

    
    def get_bonds_slice(self, cnp.int64_t[:] targets):
        cdef vector[vector[int]] bonds
        cdef vector[int] row 
        cdef int i, j
        for i in range(targets.shape[0]):
            row = self._get_bond(targets[i])
            for j in range(row.size()):
                
                bonds[i].push_back(row[j])
        return bonds

    @property
    def bonds(self):
        return self._ix_pair_array


    @property
    def types(self):
        cdef cnp.npy_intp size[1]
        size[0] = self._type.size()
        cdef cnp.ndarray arr
        arr =  to_numpy_from_spec(self, 1, size, cnp.NPY_INT32, &self._type[0])
        return arr


    @property
    def orders(self):
        cdef cnp.npy_intp size[1]
        size[0] = self._order.size()
        cdef cnp.ndarray arr
        arr =  to_numpy_from_spec(self, 1, size, cnp.NPY_INT32, &self._order[0])
        return arr

    @property
    def guessed(self):
        cdef cnp.npy_intp size[1]
        size[0] = self._guessed.size()
        cdef cnp.ndarray arr
        arr =  to_numpy_from_spec(self, 1, size, cnp.NPY_INT32, &self._guessed[0])
        return arr.astype(bool)



    cdef vector[int] _get_bond(self, int target):
        cdef int start = self._spans[target]
        cdef int end = self._spans[target+1]
        cdef int i, b_ix, first, second
        cdef vector[int] bonds
        if start == end:
            return bonds
        
        else:
            for i in range(start, end, 1):
                b_ix = self._bix[i]
                first = self._ix_pair_array[b_ix].first
                second = self._ix_pair_array[b_ix].second
                if first != target:
                    bonds.push_back(first)
                else:
                    bonds.push_back(second)
            return bonds

    @property
    def bix(self):
        return self._bix

    @property
    def ix_pair_array(self):
        return self._ix_pair_array
  
    @property
    def spans(self):
        return self._spans

    @property
    def span_map(self):
        return self._span_map    
