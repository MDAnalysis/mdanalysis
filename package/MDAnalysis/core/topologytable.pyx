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
from libcpp.pair cimport pair as cpair
from libcpp.algorithm cimport sort as csort
from libcpp.algorithm cimport unique as cunique
from ..lib._cutil cimport to_numpy_from_spec
from cython.operator cimport dereference as deref
import numpy as np
cimport numpy as cnp
cnp.import_array()


cdef class TopologyTable:
    def __cinit__(self, int[:,:] val, int n_atoms, **kwargs):
        """Initialise C++ level parameters of a TopologyTable

        Parameters
        ----------

        .. versionadded:: 2.3.0
           Initialise C++ level parameters
        """
        self.n_atoms = n_atoms
        self._generate_bix(val)
        
    
    cdef void _generate_bix(self, int[:,:] val):
        # track whether we have seen this bond before
        
        # the bond, foraward and reverse
        cdef cpair[int, int] bond
        cdef cpair[int, int] rev_bond

        # make sure each bond is in the table both forwards and backwards
        cdef int i 
        for i in range(val.shape[0]):
            bond = cpair[int, int](val[i,0], val[i,1])
            rev_bond = cpair[int, int](val[i,1], val[i,0])
            self._values.push_back(bond)
            self._values.push_back(rev_bond)

        # remove duplicates with a sort, can do these together for speed
        csort(self._values.begin(), self._values.end())
        self._values.erase( cunique( self._values.begin(), self._values.end() ), self._values.end() )


        # initialise spans
        cdef int lead_val
        cdef int prev_val = val[0,0]
        self._span_map[prev_val] = 0

        # unique value counter
        cdef int bix_counter = 0


        for i in range(self._values.size()):
            bond = self._values[i]
            rev_bond = cpair[int, int](bond.second, bond.first)
            if self._mapping.count(bond):
                # the value is already in the map, grab forward value
                # and that we will read second element
                self._bix.push_back(self._mapping[bond])
                self._access.push_back(1)
            
            elif self._mapping.count(rev_bond):
                # the reversed value is already in the map, grab reverse value
                # and that we will read first element
                self._bix.push_back(self._mapping[rev_bond])
                self._access.push_back(0)

            else:
                # new value
                self._mapping.insert(cpair[cpair[int,int], int](bond, bix_counter))
                self._bix.push_back(bix_counter)
                self._access.push_back(1)

                # increment unique values counter
                bix_counter += 1
                # save new value to ix_pair array
                self._ix_pair_array.push_back(bond)
            
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
        bonds = self._get_bonds(target)
        return bonds

    cdef vector[int] _get_bonds(self, int target):
        cdef int start = self._spans[target]
        cdef int end = self._spans[target+1]
        cdef int i, b_ix, atom
        cdef vector[int] bonds
        if start == end:
            return bonds
        
        else:
            for i in range(start, end, 1):
                b_ix = self._bix[i]
                #NOTE: Can we get rid of this ugly access thing?
                if self._access[i] == 0:
                    bond = self._ix_pair_array[b_ix].first
                else:
                    bond = self._ix_pair_array[b_ix].second
                bonds.push_back(bond)
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

    @property
    def values(self):
        return self._values









        

    

        
            
