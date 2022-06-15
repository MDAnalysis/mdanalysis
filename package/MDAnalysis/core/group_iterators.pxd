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

from ..libmda cimport iterators
import numpy
from libc.stdint cimport int64_t 
cimport cython
cimport numpy as cnp
cnp.import_array()

"""
Cython header for iterators of groups, used in low level functions in `lib`.
Light wrapper around C++ iterator classes in `lib/include/iterators.h`
For full documentation, see `group_iterators.pyx`.
"""

cdef class AtomGroupIterator:
    # C++ class for iteration of an AtomGroup 
    cdef iterators._AtomGroupIterator _iterator
    # number of atoms in the AtomGroup
    cdef int64_t n_atoms
    # view of coordinates, requires target array C contiguous
    cdef float [:, ::1] _coord_view



cdef class ArrayIterator:
    # C++ class for iteration of a NumPy Array
    cdef iterators._ArrayIterator _iterator
    # number of atoms in the array
    cdef int64_t n_atoms
    # view of coordinates, requires target array C contiguous
    cdef float[:, ::1] _coord_view

