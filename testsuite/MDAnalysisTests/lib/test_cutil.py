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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import absolute_import

import pytest
import numpy as np
from numpy.testing import assert_equal

from MDAnalysis.lib._cutil import unique_int_1d


@pytest.mark.parametrize('values', (
    [],  # empty array
    [1, 1, 1, 1, ],  # all identical
    [2, 3, 5, 7, ],  # all different, monotonic
    [5, 2, 7, 3, ],  # all different, non-monotonic
    [1, 2, 2, 4, 4, 6, ],  # duplicates, monotonic
    [1, 2, 2, 6, 4, 4, ],  # duplicates, non-monotonic
))
def test_unique_int_1d(values):
    array = np.array(values, dtype=np.int64)
    assert_equal(unique_int_1d(array), np.unique(array))
