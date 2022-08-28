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


from libcpp cimport bool as cbool
from libcpp.map cimport map as cmap
from libcpp.pair cimport pair as cpair
from libcpp.set cimport set as cset
from libcpp.vector cimport vector

import numpy as np

cimport numpy as cnp

cnp.import_array()


cdef class BondTable:
    # attributes are public for testing but are not intended to be accessed
    
    # maximum index in the table 
    cdef public int max_index
    # whether the table has any values in it 
    cdef public cbool _is_empty
    # which bond index entry in _ix_pair_array is value in input
    cdef public vector[int] _bix
    # span of each atom in _bix
    cdef public vector[int] _spans_start
    cdef public vector[int] _spans_end

    # unique list of bonds
    cdef public vector[cpair[int, int]] _ix_pair_array
    # unique list of types
    cdef public list _type
    # unique list of orders
    cdef public list _order
    # unique list of guesses
    cdef public vector[int] _guessed

    # convert inputs to table representation 
    cdef void _generate_bix(self, int[:, :] val, list typ, int[:] guess,
                            list order)

    # utility method to sort vector a and the vector b by the same index.
    cdef void _pairsort(self, vector[cpair[int, int]] & a, vector[int] & b)
    # utility method to sort vector a and the list b by the same index.
    cdef list _pairsort_list(self, vector[cpair[int, int]] a, list b)

