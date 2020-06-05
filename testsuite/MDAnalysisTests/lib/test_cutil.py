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
from __future__ import absolute_import

import pytest
import numpy as np
from numpy.testing import assert_equal

from MDAnalysis.lib._cutil import unique_int_1d, find_fragments


@pytest.mark.parametrize('values', (
    [],  # empty array
    [1, 1, 1, 1, ],  # all identical
    [2, 3, 5, 7, ],  # all different, monotonic
    [5, 2, 7, 3, ],  # all different, non-monotonic
    [1, 2, 2, 4, 4, 6, ],  # duplicates, monotonic
    [1, 2, 2, 6, 4, 4, ],  # duplicates, non-monotonic
))
def test_unique_int_1d(values):
    array = np.array(values, dtype=np.intp)
    ref = np.unique(array)
    res = unique_int_1d(array)
    assert_equal(res, ref)
    assert type(res) == type(ref)
    assert res.dtype == ref.dtype


@pytest.mark.parametrize('edges,ref', [
    ([[0, 1], [1, 2], [2, 3], [3, 4]],
     [[0, 1, 2, 3, 4]]),  # linear chain
    ([[0, 1], [1, 2], [2, 3], [3, 4], [4, 10]],
     [[0, 1, 2, 3, 4]]),  # unused edge (4, 10)
    ([[0, 1], [1, 2], [2, 3]],
     [[0, 1, 2, 3], [4]]),  # lone atom
    ([[0, 1], [1, 2], [2, 0], [3, 4], [4, 3]],
     [[0, 1, 2], [3, 4]]),  # circular
])
def test_find_fragments(edges, ref):
    atoms = np.arange(5)

    fragments = find_fragments(atoms, edges)

    assert len(fragments) == len(ref)
    for frag, r in zip(fragments, ref):
        assert_equal(frag, r)
