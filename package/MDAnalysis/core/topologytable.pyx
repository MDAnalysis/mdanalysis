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
from libcpp.map cimport map
from ..lib._cutil import unique_int_1d
import numpy as np
cimport numpy as cnp
cnp.import_array()


cdef class TopologyTable:
    def __cinit__(self, unsigned int nval, unsigned int npair, int[:,:] val, int[:] flat_vals,  int[:] typ, int[:] guess, int[:] order,  **kwargs):
        """Initialise C++ level parameters of a TopologyTable

        Parameters
        ----------
        order : unsigned int
          The order of the TopologyTable

        .. versionadded:: 2.3.0
           Initialise C++ level parameters
        """
    
        # c++ level objects
        print("cinit")
        self._nval = nval
        self._npair = npair
        self._unique = unique_int_1d(flat_vals)
        self._nunique = unique.shape[0]
        self._construct_empty_tables()
        self._parse(val, typ, guess, order)
        print("parsed")
    
    cdef _construct_empty_tables(self):
        cdef size_t i
        cdef vector[int] tmp
        tmp.reserve(8)
        # cdef int t
        # for _ in range(10):
        #     tmp.push_back(0)

        for i in range(self._nunique):
            self.values.push_back(tmp)
            self.types.push_back(tmp)
            self.guessed.push_back(tmp)
    
    def _parse(self, int[:,:] val,  int[:] typ, int[:] guess, int[:] ord):
        cdef int i
        for i in range(self._npair):
            print(f" {i}, {val[i,0]}, {val[i,1]}")
            self.values[val[i,0]].push_back(val[i,1])
            self.values[val[i,1]].push_back(val[i,0])  


    def print_values(self):
        for i in range(self._nval):
            print(f" {i} {self.values[i]}")

    def print_types(self):
        for i in range(self._nval):
            print(self.types[i])

    def print_guessed(self):
        for i in range(self._nval):
            print(self.guessed[i])

    def print_order(self):
        for i in range(self._nval):
            print(self.order[i])



    

