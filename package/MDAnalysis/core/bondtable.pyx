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

"""\
BondTable --- :mod:`MDAnalysis.core.bondtable`
========================================================

.. versionadded:: 2.3.0

:class:`BondTable` is a table-like representation of bonds and bond
attributes that allows fast lookup.  This is done by mapping the bonds to a
linear representation that can be looked up using a span. This gives better
scaling than using a key:value paired datastructure. 

"""

from cython.operator cimport dereference as deref
from cpython.ref cimport PyObject
from libcpp cimport bool as cbool
from libcpp.algorithm cimport sort as csort
from libcpp.map cimport map as cmap
from libcpp.pair cimport pair as cpair
from libcpp.set cimport set as cset
from libcpp.vector cimport vector

from ..lib._cutil cimport to_numpy_from_spec

import numpy as np
import cython 

cimport numpy as cnp

cnp.import_array()

# index sort used for sorting lists in _pairsort_list
cdef extern from *:
    """
    #include <algorithm>
    #include <vector>
    void sort_via_score(std::vector<size_t>& indices, const std::vector<std::pair<int,int>>& scores){
        std::sort(indices.begin(), indices.end(),
                  [&scores](int i, int j){return scores.at(i)<scores.at(j);}
                 );
    }
    """
    void sort_via_score(vector[size_t] & indices, vector[cpair[int, int]] & scores)


cdef class BondTable:
    """
    A table-like structure used to look up bonds with better scaling than a
    dictionary or map. This is acheived by mapping the input array of indices
    to a unique array of indices which are then accessed using two spans,
    bypassing the need to look up the bonds in a key:value paired
    datastructure. 

    Notes
    -----

    Using `__cinit__` is avoided here to enable pickling of the table.

    .. versionadded:: 2.3.0
    """

    def __init__(self, val, typ,  guess, order, **kwargs):
        """
        Initialise a BondTable

        Parameters
        ----------
        val: np.ndarray
            2D array of integers of size (nbonds, 2) that define the input
            bonds for the BondTable
        typ: list
            A list of corresponding bond type objects of length (nbonds)
        guess: list[bool]
            A list of bools of whether the bonds were guessed or not of length
            (nbonds).
        order: list
            A list of corresponding bond order objects of length (nbonds).

        .. versionadded:: 2.3.0
        """

        self._is_empty = False
        if len(val) == 0:
            self._is_empty = True
        if not self._is_empty:
            # make sure the inputs are sane
            guess_ = np.asarray([int(i) for i in guess], dtype=np.int32)
            order_ = list(order)
            typ_ = list(typ)

            # make table
            self._type = []
            self._order = []
            self._generate_bix(val, typ_,  guess_, order_)

    cdef void _generate_bix(self, int[:, :] val, list typ, int[:] guess,
                            list order):
        """
        Generate bond indicies, spans and value arrays for the BondTable
        """
        cdef vector[cpair[int, int]] _values

        # types must be a list as we allow flexibility in what it can be
        cdef list  _type = []
        cdef list  _order = []
        cdef vector[int]  _guessed

        # a bond, forward and reverse
        cdef cpair[int, int] bond
        cdef cpair[int, int] rev_bond

        # set used to check if reverse bonds already present
        cdef cset[cpair[int, int]] _tmp_set

        # map of unique spans
        cdef cmap[int, int] _span_map
        cdef cmap[int, int] _span_map_end

        # mapping of bonds:bix
        cdef cmap[cpair[int, int], int] _mapping

        cdef int i
        for i in range(val.shape[0]):
            bond = cpair[int, int](val[i, 0], val[i, 1])
            _tmp_set.insert(bond)

        # make sure each bond is in the table both forwards and backwards if
        # it isn't already
        for i in range(val.shape[0]):
            bond = cpair[int, int](val[i, 0], val[i, 1])
            rev_bond = cpair[int, int](val[i, 1], val[i, 0])
            _values.push_back(bond)
            _type.append(typ[i])
            _guessed.push_back(guess[i])
            _order.append(order[i])

            if not _tmp_set.count(rev_bond):
                # reverse bond not present
                _values.push_back(rev_bond)
                _type.append(typ[i])
                _order.append(order[i])
                _guessed.push_back(guess[i])

        # pairwise sort each array by the array of values
        _type = self._pairsort_list(_values, _type)
        _order = self._pairsort_list(_values, _order)
        self._pairsort(_values, _guessed)

        # initialise spans
        cdef int lead_val
        cdef int prev_val = val[0, 0]
        _span_map[prev_val] = 0

        # unique value counter
        cdef int bix_counter = 0

        for i in range(_values.size()):
            bond = _values[i]
            rev_bond = cpair[int, int](bond.second, bond.first)
            if _mapping.count(bond):
                # the value is already in the map, grab forward value
                self._bix.push_back(_mapping[bond])

            elif _mapping.count(rev_bond):
                # the reversed value is already in the map, grab reverse value
                self._bix.push_back(_mapping[rev_bond])

            else:
                # new value
                _mapping.insert(
                    cpair[cpair[int, int], int](bond, bix_counter))

                self._bix.push_back(bix_counter)

                # increment unique values counter
                bix_counter += 1
                # save new value to ix_pair array
                self._ix_pair_array.push_back(bond)
                self._type.append(_type[i])
                self._order.append(_order[i])
                self._guessed.push_back(_guessed[i])

            lead_val = bond.first

            if lead_val != prev_val:
                _span_map[lead_val] = i
                _span_map_end[prev_val] = i

            prev_val = lead_val

        _span_map_end[lead_val] = i + 1

       # need to find the maximum key and value to make sure no overrun
        cdef int max_val, max_key

        max_key = deref(_span_map.rbegin()).first
        max_val = deref(_span_map.rbegin()).second

        self.max_index = max_key
        # sort out the spans so that each atoms has a span
        prev_val = -1
        for i in range(self.max_index + 1):
            if _span_map.count(i):
                self._spans_start.push_back(_span_map[i])
                self._spans_end.push_back(_span_map_end[i])
            else:
                self._spans_start.push_back(-1)
                self._spans_end.push_back(-1)

    @cython.boundscheck(False)
    @cython.wraparound(False) 
    def get_b_t_g_o_slice(self, targets):
        """
        Get a slice of all four properties (bonds, types, guesses, orders)
        for an array of indices

        Parameters
        ----------
        targets: np.ndarray 
            the atom indices to get properties for

        Returns
        -------
        pairs: np.ndarray
            An array of bonds that contain the target index 
        types: np.ndarray
            An array of bond type objects that contain the target index
        guesses: np.ndarray
            An array booleans for bonds that contain the target index
        orders: np.ndarray
            An array of bond order objects that contain the target index
        """
        # if empty return empty arrays
        if self._is_empty:
            return np.empty((0), dtype=np.int32), np.empty((0), dtype=object),\
                np.empty((0), dtype=bool), np.empty((0), dtype=object)
        # if scalar, make it an array
        if np.isscalar(targets):
            targets = np.asarray([targets])
        # objects for pairs
        cdef vector[cpair[int, int]] bonds
        # objects for types
        cdef list types = []
        # objects for guesses
        cdef vector[int] guesses
        # objects for orders
        cdef list orders = []
        # iteration
        cdef int i, j, idx, start, end, b_ix, guess
        for i in range(targets.shape[0]):
            idx = targets[i]
            if idx <= self.max_index:
                start = self._spans_start[idx]
                end = self._spans_end[idx]
                if start != end:
                    for j in range(start, end, 1):
                        b_ix = self._bix[j]
                        bonds.push_back(self._ix_pair_array[b_ix])
                        types.append(self._type[b_ix])
                        orders.append(self._order[b_ix])
                        guesses.push_back(self._guessed[b_ix]) 

        return np.asarray(bonds, dtype=np.int32), \
            np.asarray(types, dtype=object), \
            np.asarray(guesses, dtype=bool), \
            np.asarray(orders, dtype=object)

    cdef void _pairsort(self, vector[cpair[int, int]] & a, vector[int] & b):
        """
        Sort a vector of integer pairs **a** and sort a vector of integers
        **b** by the corresponding index in **a**. Uses sort from the
        <algorithm> header and operates inplace on both input arguments.
        """
        cdef vector[cpair[cpair[int, int], int]] pair_arr
        cdef cpair[cpair[int, int], int] pair
        for i in range(a.size()):
            pair = cpair[cpair[int, int], int](a[i], b[i])
            pair_arr.push_back(pair)

        csort(pair_arr.begin(), pair_arr.end())

        for i in range(a.size()):
            a[i] = pair_arr[i].first
            b[i] = pair_arr[i].second

    cdef list _pairsort_list(self, vector[cpair[int, int]] a, list b):
        """
        Sort a vector of integer pairs **a** and a list of Python objects
        (PyObject*) by the corresponding index in **a**. Uses an index sort and
        only returns the sorted list input of **b**, leaving the input
        arguments unchanged.
        """
        cdef int size = a.size()
        cdef size_t index
        cdef vector[size_t] indices

        for i in range(size):
            indices.push_back(i)

        sort_via_score(indices, a)

        cdef list b_ = []
        for i in range(size):
            index = indices[i]
            b_.append(b[index])

        return b_
