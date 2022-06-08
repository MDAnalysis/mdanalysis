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


"""
Compiled helpers for group iteration --- :mod:`MDAnalysis.core.group_helpers`
=============================================================================

Helpers
"""

import numpy
from libcpp.vector cimport vector
from libc.stdint cimport uint64_t, UINT64_MAX
cimport cython
cimport numpy as cnp
cnp.import_array()


cdef class AtomGroupIterHelper:

    def __cinit__(self, uint64_t n_atoms ** kwargs):
        self.n_atoms = n_atoms
        self._i = 0
        # self._coord_view = 

    # buffer passed in must be 3x as large as n_idx
    cdef void _load_to_external_buffer(self, float * buffer, uint64_t n_idx):
        cdef uint64_t i
        for i in range(n_idx):
            buffer[3*self._i] = self._coord_view[3*self._i,0]
            buffer[3*self._i + 1] = self._coord_view[3*self._i + 1,1]
            buffer[3*self._i + 2] = self._coord_view[3*self._i + 2,2]
            self._i += 1

    cdef void _reset_iteration(self):
        self._i = 0

    
        
