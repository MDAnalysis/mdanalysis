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

import pytest
import numpy as np
from numpy.testing import assert_equal

import MDAnalysis as mda

def test_transform_StoR_pass():
    box = np.array([10, 7, 3, 45, 60, 90], dtype=np.float32)
    s = np.array([[0.5, -0.1, 0.5]], dtype=np.float32)

    original_r = np.array([[ 5.75,  0.36066014, 0.75000012]], dtype=np.float32)

    test_r = mda.lib.distances.transform_StoR(s, box)

    assert_equal(original_r, test_r)

def test_transform_StoR_fail():
    box = np.array([10, 7, 3, 45, 60, 90], dtype=np.float32)
    s = np.array([[0.5, -0.1, 0.5]])

    with pytest.raises(TypeError, match='S must be of type float32'):
        r = mda.lib.distances.transform_StoR(s, box)
