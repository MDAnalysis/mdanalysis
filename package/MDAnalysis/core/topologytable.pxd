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


from libcpp.vector cimport vector
from libcpp.map cimport map as cmap
from libcpp.pair cimport pair as cpair
from libcpp.set cimport set as cset
from libcpp.string cimport string as cstring
import numpy as np
cimport numpy as cnp
cnp.import_array()


cdef class TopologyTable:
    # number of atoms in the table
    cdef int max_index
    # which bond index entry in _ix_pair_array is value in input
    cdef vector[int] _bix
    # span of each atom in _bix
    cdef vector[int] _spans
    # which element of each bond do we need first or second
    cdef cmap[cpair[int, int], int] _mapping
    cdef vector[int] _type
    cdef vector[int] _order
    cdef vector[int] _guessed
    cdef vector[cpair[int, int]] _ix_pair_array
    cdef cmap[int, int] _span_map
    cdef void _generate_bix(self, int[:, :] val, int[:] typ, int[:] guess,
                            int[:] order)

    cdef  vector[int] _get_bond(self, int target)

    cdef void _pairsort(self, vector[cpair[int, int]] & a, vector[int] & b)
