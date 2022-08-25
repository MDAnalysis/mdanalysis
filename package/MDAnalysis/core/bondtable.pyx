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
# cython: linetrace=True

"""\
BondTable --- :mod:`MDAnalysis.core.bondtable`
========================================================

.. versionadded:: 2.3.0

:class:`BondTable` is a table-like representation of bonds and bond
attributes that allows fast lookup. 

TODO: Add in-depth discussion.
"""

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
    to a unique array of indices which are then accessed using a span, bypassing
    the need to look up the bonds in a key:value paired datastructure. 

    Notes
    -----

    Using `__cinit__` is avoided here to enable pickling of the table.
    """
    def __init__(self, val, typ,  guess, order, **kwargs):
        """Initialise a BondTable

        Parameters
        ----------
        val: 2D memoryview of integers
            2D array of integers of size (nbonds, 2) that define the input
            bonds for the topologytable
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
            guess_ = np.asarray([int(i) for i in guess], dtype=np.int32)
            order_ = list(order)
            typ_ = list(typ)
            self._type = []
            self._order = []
            self._generate_bix(val, typ_,  guess_, order_)

    cdef void _pairsort(self, vector[cpair[int, int]] & a, vector[int] & b):
        """
        Sort a vector of integer pairs ** a ** and sort a vector of integers
        ** b ** by the corresponding index in **a**. Uses sort from the
        <algorithm> header and operates inplace on both input arguments.
        """
        cdef vector[cpair[cpair[int, int], int]] pair_arr
        cdef cpair[cpair[int, int], int] pair
        cdef cpair[int, int] ptmp
        cdef int itmp
        for i in range(a.size()):
            ptmp = a[i]
            itmp = b[i]
            pair = cpair[cpair[int, int], int](ptmp, itmp)
            pair_arr.push_back(pair)

        csort(pair_arr.begin(), pair_arr.end())

        for i in range(a.size()):
            a[i] = pair_arr[i].first
            b[i] = pair_arr[i].second

    def _pairsort_list(self, vector[cpair[int, int]] a, list b):
        """
        Sort a vector of integer pairs ** a ** and a list of Python objects
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

    cdef void _generate_bix(self, int[:, :] val, list typ, int[:] guess, list order):
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

        cdef int i
        cdef cset[cpair[int, int]] _tmp_set
        for i in range(val.shape[0]):
            bond = cpair[int, int](val[i, 0], val[i, 1])
            _tmp_set.insert(bond)

        # make sure each bond is in the table both forwards and backwards IF IT
        # ISNT ALREADY
        for i in range(val.shape[0]):
            bond = cpair[int, int](val[i, 0], val[i, 1])
            rev_bond = cpair[int, int](val[i, 1], val[i, 0])
            _values.push_back(bond)
            _type.append(typ[i])
            _guessed.push_back(guess[i])
            _order.append(order[i])

            if _tmp_set.count(rev_bond):
                pass  # the reverse bond is already present, do not add
            else:
                _values.push_back(rev_bond)
                _type.append(typ[i])
                _order.append(order[i])
                _guessed.push_back(guess[i])

        # pairwise sort each array
        _type = self._pairsort_list(_values, _type)
        _order = self._pairsort_list(_values, _order)
        self._pairsort(_values, _guessed)

        # print(_values)
        # initialise spans
        cdef int lead_val
        cdef int prev_val = val[0, 0]
        self._span_map[prev_val] = 0

        # unique value counter
        cdef int bix_counter = 0

        for i in range(_values.size()):
            bond = _values[i]
            rev_bond = cpair[int, int](bond.second, bond.first)
            # print(bond)
            if self._mapping.count(bond):
                # the value is already in the map, grab forward value
                # and that we will read second element
                self._bix.push_back(self._mapping[bond])

            elif self._mapping.count(rev_bond):
                # the reversed value is already in the map, grab reverse value
                self._bix.push_back(self._mapping[rev_bond])

            else:
                # new value
                self._mapping.insert(
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
                self._span_map[lead_val] = i
                self._span_map_end[prev_val] = i

            prev_val = lead_val
        
        print(self._span_map)
        print(prev_val)

        #self._span_map[prev_val] = i
        # WORKS?
        self._span_map_end[lead_val] = i +1


       # need to find the maximum key and value to correctly cap span array
        cdef int max_val, max_key

        max_key = deref(self._span_map.rbegin()).first
        max_val = deref(self._span_map.rbegin()).second
        # print("span_map")
        # print(self._span_map)
        self.max_index = max_key #+ 2  # +2 so that the fake key made below is used
        # pop fake element into map to correctly cap the span
        # self._span_map[max_key + 1] = max_val + 1

        # print("span_map")
        # print(self._span_map)
        # print("max_index")
        # print(self.max_index)
        # sort out the spans so that each atoms has a span
        prev_val = -1
        for i in range(self.max_index +1 ):
            if self._span_map.count(i):
                self._spans.push_back(self._span_map[i])
                self._spans_end.push_back(self._span_map_end[i])
            else:
                self._spans.push_back(-1)
                self._spans_end.push_back(-1)


        # remove fake element and decrement maximum index
        # self._span_map.erase(max_key + 1)
        # self.max_index -= 2
        # print("span_map")
        # print(self._span_map)
        # print("max_index")
        # print(self.max_index)

    def get_bonds(self, int target):
        """
        Get the bonded atoms for an atom in the BondTable

        Parameters
        ----------
        target: int
            The index of the atom to get bonds for.
        
        Returns
        -------
        bonds: list[int]
            A list of indices of atoms bonded to the target index.
        """
        cdef vector[int] bonds
        if self._is_empty:
            return bonds
        if target <= self.max_index:
            bonds = self._get_bond(target)
        else:
            pass
        return bonds

    def get_pairs(self, int target):
        """
        Get the bond pairs for an atom in the BondTable

        Parameters
        ----------
        target: int
            The index of the atom to get bonds for.
        
        Returns
        -------
        pairs: list
            A list of bonds that contain the target index 
        """
        cdef vector[cpair[int, int]] pairs
        cdef cpair[vector[cpair[int, int]], cbool] ret_struct
        if self._is_empty:
            return []
        if target <= self.max_index:
            ret_struct = self._get_pair(target)
            if ret_struct.second:
                pairs = ret_struct.first
                return pairs
            else:
                return []
        else:
            return []

    def get_bonds_slice(self, targets):
        """
        Get bonded atoms for an array of indices.

        Parameters
        ----------
        targets: np.ndarray or scalar index
            the atoms to get bonded atoms for 

        Returns
        -------
        bonds: list
            The bonded atoms for the indices in targets
        """
        if self._is_empty:
            return np.empty((0), dtype=np.int32)
        if np.isscalar(targets):
            targets = np.asarray([targets])
        cdef vector[vector[int]] bonds
        cdef vector[int] row, tmp
        cdef int i, j, size_j
        for i in range(targets.shape[0]):
            index = targets[i]
            if index <= self.max_index:
                row = self._get_bond(index)
                for j in range(row.size()):
                    size_j += 1
                    tmp = vector[int](2)
                    tmp[0] = targets[i]
                    tmp[1] = row[j]
                    bonds.push_back(tmp)
            else:
                pass
        return np.asarray(bonds, dtype=np.int32)

    def get_pairs_slice(self,  targets):
        """
        Get the bond pairs for an array of indices
        
        Parameters
        ----------
        targets: np.ndarray or scalar index
            the atoms to get bond pairs for

        Returns
        -------
        pairs: np.ndarray
            An array of bonds that contain the target index 
        """
        if self._is_empty:
            return np.empty((0), dtype=np.int32)
        if np.isscalar(targets):
            targets = np.asarray([targets])
        cdef vector[cpair[int, int]] bonds
        cdef vector[cpair[int, int]] pair_arr
        cdef cpair[vector[cpair[int, int]], cbool] ret_struct
        cdef cbool has_bonds
        cdef int i, j
        for i in range(targets.shape[0]):
            if targets[i] <= self.max_index:
                ret_struct = self._get_pair(targets[i])
                pair_arr = ret_struct.first
                has_bonds = ret_struct.second
                if has_bonds:
                    for j in range(pair_arr.size()):
                        bonds.push_back(pair_arr[j])
            else:
                pass
        return np.asarray(bonds, dtype=np.int32)

    def get_types_slice(self, targets):
        """
        Get bond types for an array of indices
        
        Parameters
        ----------
        targets: np.ndarray or scalar index
            the atoms to get types for

        Returns
        -------
        types: np.ndarray
            An array of bond type objects that contain the target index
        """
        if self._is_empty:
            return np.empty((0), dtype=object)
        if np.isscalar(targets):
            targets = np.asarray([targets])
        types = []
        typ_arr = []
        cdef cpair[vector[int], cbool] ret_struct
        cdef cbool has_typ
        cdef int i, j
        cdef cnp.ndarray arr
        for i in range(targets.shape[0]):
            if targets[i] <= self.max_index:
                typ_arr, has_typ = self._get_typ(targets[i])
                if has_typ:
                    for j in range(len(typ_arr)):
                        types.append(typ_arr[j])
            else:
                pass
        arr = np.asarray(types, dtype=object)
        return arr

    def get_guess_slice(self, targets):
        """
        Get whether bonds have been guessed for an array of indices
        
        Parameters
        ----------
        targets: np.ndarray or scalar index
            the atoms to get guesses for

        Returns
        -------
        guesses: np.ndarray
            An array booleans for bonds that contain the target index
        """
        if self._is_empty:
            return np.empty((0), dtype=bool)
        if np.isscalar(targets):
            targets = np.asarray([targets])
        cdef vector[int] guesses
        cdef vector[int] guess_arr
        cdef cpair[vector[int], cbool] ret_struct
        cdef cbool has_guess
        cdef int i, j
        for i in range(targets.shape[0]):
            if targets[i] <= self.max_index:
                ret_struct = self._get_guess(targets[i])
                guess_arr = ret_struct.first
                has_guess = ret_struct.second
                if has_guess:
                    for j in range(guess_arr.size()):
                        guesses.push_back(guess_arr[j])
            else:
                pass

        return np.asarray(guesses, dtype=bool)

    def get_order_slice(self, targets):
        """
        Get bond orders for an array of indices
        
        Parameters
        ----------
        targets: np.ndarray or scalar index
            the atoms to get orders for

        Returns
        -------
        orders: np.ndarray
            An array of bond order objects that contain the target index
        """
        if self._is_empty:
            return np.empty((0), dtype=object)
        if np.isscalar(targets):
            targets = np.asarray([targets])
        orders = []
        orders_arr = []
        cdef cpair[vector[int], cbool] ret_struct
        cdef cbool has_order
        cdef int i, j
        cdef cnp.ndarray arr
        for i in range(targets.shape[0]):
            if targets[i] <= self.max_index:
                orders_arr, has_order = self._get_ord(targets[i])
                if has_order:
                    for j in range(len(orders_arr)):
                        orders.append(orders_arr[j])
            else:
                pass

        arr = np.asarray(orders).astype(object)
        return arr

    def get_b_t_g_o_slice(self, targets):
        """
        Get a slice of all four properties (bonds, types, guesses, orders)
        for an array of indices

        Parameters
        ----------
        targets: np.ndarray or scalar index
            the atoms to get properties for

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
        if self._is_empty:
            return np.empty((0), dtype=np.int32), np.empty((0), dtype=object), np.empty((0),dtype=bool), np.empty((0), dtype=object)
        if np.isscalar(targets):
            targets = np.asarray([targets])
        cdef vector[cpair[int, int]] bonds
        cdef vector[cpair[int, int]] pair_arr
        cdef cpair[vector[cpair[int, int]], cbool] ret_struct_bond
        cdef int i, j
        cdef list types = []
        cdef list typ_arr = []
        cdef vector[int] guesses
        cdef vector[int] guess_arr
        cdef list orders = []
        cdef list orders_arr = []
        cdef cpair[vector[int], cbool] ret_struct_guess
        cdef cbool has_bonds, has_typ, has_guess, has_order
        for i in range(targets.shape[0]):
            if targets[i] <= self.max_index:
                
                ret_struct_bond = self._get_pair(targets[i])
                pair_arr = ret_struct_bond.first
                has_bonds = ret_struct_bond.second
                if has_bonds:
                    for j in range(pair_arr.size()):
                        bonds.push_back(pair_arr[j])
                typ_arr, has_typ = self._get_typ(targets[i])
                if has_typ:
                    for j in range(len(typ_arr)):
                        types.append(typ_arr[j])
                ret_struct_guess = self._get_guess(targets[i])
                guess_arr = ret_struct_guess.first
                has_guess = ret_struct_guess.second
                if has_guess:
                    for j in range(guess_arr.size()):
                        guesses.push_back(guess_arr[j])
                orders_arr, has_order = self._get_ord(targets[i])
                if has_order:
                    for j in range(len(orders_arr)):
                        orders.append(orders_arr[j])
            else:
                pass
        return np.asarray(bonds, dtype=np.int32), np.asarray(types, dtype=object), np.asarray(guesses, dtype=bool), np.asarray(orders, dtype=object)

    @property
    def bonds(self):
        """
        The unique bonds in the BondTable
        
        Returns
        -------
        bonds: np.ndarray
            Unique bond pairs in the BondTable
        """
        return np.asarray(self._ix_pair_array, dtype=np.int32)

    @property
    def types(self):
        """
        The types of the unique bonds in the BondTable 

        Returns
        -------
        types: np.ndarray
            The types of the unique bonds
        """
        return np.asarray(self._type, dtype=object)

    @property
    def orders(self):
        """
        The orders of the unique bonds in the BondTable

        Returns
        -------
        orders: np.ndarray
            The orders of the unique bonds
        """
        return np.asarray(self._order, dtype=object)

    @property
    def guessed(self):
        """
        Whether the bond has been guessed for the unique bonds in the BondTable

        Returns
        -------
        guesses: np.ndarray
            Whether each unique bond has been guessed
        """
        return np.asarray(self._guessed, dtype=bool)

    cdef vector[int] _get_bond(self, int target):
        """
        Low level utility to get the bonds for a single index

        Parameters
        ----------
        target: int
            The target atom to get bonds for
        
        Returns
        -------
        bonds: vector[int]
            The bonded atoms for the target index
        """
        # private, does not check for target < self.max_index
        cdef int start = self._spans[target]
        cdef int end = self._spans_end[target]
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

    cdef cpair[vector[cpair[int, int]], cbool] _get_pair(self, int target):
        """
        Low level utility to get the bond pairs for a single index

        Parameters
        ----------
        target: int
            The target atom to get bonds for
        
        Returns
        -------
        struct: pair[vector[pair[int,int]], bool]
            Pair containing the vector of pairs, and whether the target index
            had any bonded atoms. 
        """
        # private, does not check for target < self.max_index
        cdef int start = self._spans[target]
        cdef int end = self._spans_end[target]
        cdef int i, b_ix, first, second
        cdef vector[cpair[int, int]] bonds
        if start == end:
            return cpair[vector[cpair[int, int]], cbool](bonds, False)
        else:
            for i in range(start, end, 1):
                b_ix = self._bix[i]
                pair = self._ix_pair_array[b_ix]
                bonds.push_back(pair)

            return cpair[vector[cpair[int, int]], cbool](bonds, True)

    cdef _get_typ(self, int target):
        """
        Low level utility to get the types of the bonds for a single index

        Parameters
        ----------
        target: int
            The target atom to get bonds for
        
        Returns
        -------
        types: list
            The type objects of the bonds
        has_bonds: bool
            If there are any bonded atoms for the target index
        """
        # private, does not check for target < self.max_index
        cdef int start = self._spans[target]
        cdef int end = self._spans_end[target]
        cdef int i, b_ix
        cdef list types = []
        if start == end:
            return types, False
        else:
            for i in range(start, end, 1):
                b_ix = self._bix[i]
                typ = self._type[b_ix]
                types.append(typ)

            return types, True

    cdef _get_ord(self, int target):
        """
        Low level utility to get the orders of the bonds for a single index

        Parameters
        ----------
        target: int
            The target atom to get bonds for
        
        Returns
        -------
        orders: list
            The order objects of the bonds
        has_bonds: bool
            If there are any bonded atoms for the target index
        """
        # private, does not check for target < self.max_index
        cdef int start = self._spans[target]
        cdef int end = self._spans_end[target]
        cdef int i, b_ix
        cdef list orders = []
        if start == end:
            return orders, False
        else:
            for i in range(start, end, 1):
                b_ix = self._bix[i]
                ord = self._order[b_ix]
                orders.append(ord)

            return orders, True

    cdef cpair[vector[int], cbool] _get_guess(self, int target):
        """
        Low level utility to get the guesses of the bonds for a single index

        Parameters
        ----------
        target: int
            The target atom to get bonds for
        
        Returns
        -------
        guesses: list[bool]
            If the bond for this index was guessed
        struct: pair[vector[int], bool]
            Pair containing the vector of guesses (cast to int), and whether the
            target index had any bonded atoms. 
        """
        # private, does not check for target < self.max_index
        cdef int start = self._spans[target]
        cdef int end = self._spans_end[target]
        cdef int i, b_ix, guess
        cdef vector[int] guesses
        if start == end:
            return cpair[vector[int], cbool](guesses, False)
        else:
            for i in range(start, end, 1):
                b_ix = self._bix[i]
                guess = self._guessed[b_ix]
                guesses.push_back(guess)

            return cpair[vector[int], cbool](guesses, True)