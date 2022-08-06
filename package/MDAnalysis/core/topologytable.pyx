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
from ..lib._cutil cimport to_numpy_from_spec
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

        # initialise spans
        self.spans.push_back(0)
        cdef int lead_val
        cdef int prev_val = val[0,0]

        self.span_map[prev_val] = 0
        # unique value counter
        cdef int bix_counter = 0

        cdef int i 
        for i in range(val.shape[0]):
            bond = cpair[int, int](val[i,0], val[i,1])
            rev_bond = cpair[int, int](val[i,1], val[i,0])
            self.atoms_set.insert(val[i,0])
            self.atoms_set.insert(val[i,1])
            if self.mapping.count(bond):
                # the value is already in the map, grab forward value
                # and that we will read second element
                self.bix.push_back(self.mapping[bond])
                self.access.push_back(1)
            
            elif self.mapping.count(rev_bond):
                # the reversed value is already in the map, grab reverse value
                # and that we will read first element
                self.bix.push_back(self.mapping[rev_bond])
                self.access.push_back(0)

            else:
                # new value
                self.mapping.insert(cpair[cpair[int,int], int](bond, bix_counter))
                self.bix.push_back(bix_counter)
                self.access.push_back(1)

                # increment unique values counter
                bix_counter += 1
                # save new value to ix_pair array
                self.ix_pair_array.push_back(bond)
            
            # sort out spans
            lead_val = bond.first
            
            if lead_val != prev_val:
                self.spans.push_back(i)
                self.span_map[lead_val] = i
            
            prev_val = lead_val

            self.aidx[]
        
        self.spans.push_back(val.shape[0])


            

    def get_bonds(self, int target):
        cdef vector[int] bonds
        bonds = self._get_bonds(target)
        return bonds

    cdef vector[int] _get_bonds(self, int target):
        cdef int start = self.spans[target]
        cdef int end = self.spans[target+1]
        cdef int i, b_ix, atom
        cdef vector[int] bonds
        for i in range(start, end, 1):
            b_ix = self.bix[i]
            #NOTE: Can we git rid of this ugly access thing?
            if self.access[i] == 0:
                bond = self.ix_pair_array[b_ix].first
            else:
                bond = self.ix_pair_array[b_ix].second
            bonds.push_back(bond)
        return bonds
    
    @property
    def span(self):
        return self.spans

    @property
    def idx(self):
        return self.aidx
    
    @property
    def bix(self):
        return self.bix

    @property
    def ix_pair_array(self):
        return self.ix_pair_array


    @property
    def span_mp(self):
        return self.span_map









        

    

        
            
