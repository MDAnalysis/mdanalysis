# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
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
"""Tests for MDAnalysis.core.bondtable.

"""
import numpy as np

from numpy.testing import (
    assert_equal,
    assert_almost_equal,
)

from MDAnalysis.core.topologyattrs import Bonds
from MDAnalysis.core.bondtable import BondTable


import MDAnalysis as mda
import pytest


class TestBondTableInputs:

    class PyObjectMock:
        pass

    @pytest.mark.parametrize('inp', [(), [], np.empty(0)])
    def test_empty_table(self, inp):
        bonds = BondTable(inp, inp, inp, inp)
        assert bonds._is_empty

    @pytest.mark.parametrize('inp', [[[[1, 2], [3, 4], [5, 6]], [1, 2], [1, 3], [1]],
                                     [[[1, 2], [0, 3]], [1, 2], [1, 3], [1]]])
    def test_mismatch_inputs(self, inp):
        with pytest.raises(ValueError, match="BondTable must be made from inputs"):
            bonds = BondTable(*inp)

    @pytest.mark.parametrize('val', [np.zeros((1, 1, 1)), np.empty(1), [[[1]]], [1]])
    def test_val_not_2D(self, val):
        with pytest.raises(ValueError, match="values argument for a BondTable must be two"):
            bonds = BondTable(val, [1], [2], [3])

    @pytest.mark.parametrize('val', [((0, 1, 3),), [[0, 1, 2, 3, 4]]])
    def test_val_second_dim_not_2(self, val):
        with pytest.raises(ValueError, match="values argument for a BondTable must have a second dimension"):
            bonds = BondTable(val, [1], [2], [3])

    @pytest.mark.parametrize('val', [np.asarray([[0.1234, 1.1234]]), np.asarray([[1+2j, 3+4j]]), np.asarray([['abc', 'def']])])
    def test_val_not_castable(self, val):
        with pytest.raises(TypeError, match="values argument cannot be cast to np.int32"):
            bonds = BondTable(val, [1], [2], [3])

    def test_guess_not_castable(self):
        with pytest.raises(TypeError, match="guesses argument cannot be cast to np.int32"):
            bonds = BondTable([[1, 2]], [1], [['abc', 'def']], [2])

    @pytest.mark.parametrize('types', [['blah'], [1], (2,), [PyObjectMock()]])
    def test_types_allows_pyobject(self, types):
        bonds = BondTable([[1, 2]], types, [2], [3])

    @pytest.mark.parametrize('orders', [['blah'], [1], (2,), [PyObjectMock()]])
    def test_orders_allows_pyobject(self, orders):
        bonds = BondTable([[1, 2]], [1], [2], orders)


class TestSimpleBonds(object):
    # all values unique so is a direct copy
    _bonds = [[0,1], [0,2], [1,2], [2,3]]
    _types = ['0', '1', '2', '3']
    _guesses = [True, False, True, False]
    _orders = ['zero', 'parsnip', 'two', 'three']

    @pytest.fixture(scope='class')
    def bonds(self):
        return Bonds(self._bonds, self._types, self._guesses, self._orders)

    def test_bonds(self, bonds):
        assert_equal(np.asarray(bonds._bondtable._ix_pair_array), np.asarray(self._bonds))
    
    def test_type(self, bonds):
        assert(bonds._bondtable._type == self._types)

    def test_order(self, bonds):
        assert(bonds._bondtable._order == self._orders)
    
    def test_max_index(self, bonds):
        assert(bonds._bondtable.max_index == 3)

    @pytest.mark.parametrize('input, expected', [(0, [np.asarray([[0,1],[0,2]]), ['0', '1'], [True, False], ['zero', 'parsnip']]),
                             (1, [np.asarray([[0,1],[1,2]]), ['0', '2'], [True, True], ['zero', 'two']]),
                             (2, [np.asarray([[0,2],[1,2],[2,3]]), ['1', '2', '3'], [False, True, False], ['parsnip', 'two', 'three']])])
    def test_b_t_g_o_scalar(self, bonds, input, expected):
        b, t, g, o = bonds._bondtable.get_b_t_g_o_slice(input)
        assert_equal(b, expected[0])
        assert_equal(t, expected[1])
        assert_equal(g, expected[2])
        assert_equal(o, expected[3])

