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

from ..libmda cimport iterators
import numpy
from libcpp.vector cimport vector
from libc.stdint cimport uint64_t, UINT64_MAX
cimport cython
cimport numpy as cnp
cnp.import_array()



cdef class AtomGroupIterator:

    def __cinit__(self, uint64_t n_atoms ** kwargs):
        self.thisptr = new iterators._AtomGroupIterator(n_atoms)
    
    cdef iterators._AtomGroupIterator get_iterator(self):
        return   self.thisptr[0] 

    def __dealloc__(self):
        del self.thisptr

cdef class ArrayIterator:

    def __cinit__(self, uint64_t n_atoms ** kwargs):
        self.thisptr = new iterators._ArrayIterator(n_atoms)

    cdef iterators._ArrayIterator get_iterator(self):
        return   self.thisptr[0] 

    def __dealloc__(self):
        del self.thisptr
    

