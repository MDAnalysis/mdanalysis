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

from __future__ import print_function, absolute_import
from six.moves import zip

import numpy as np
from numpy.testing import assert_almost_equal

from MDAnalysis.lib import distances


def test_apply_PBC():
    #
    # test triclinic cell
    #
    box = np.array([10, 7, 3, 45, 60, 90], dtype=np.float32)
    s = np.array([0.5, -0.1, 0.5], dtype=np.float32)  # fractional coordinates
    r = distances.transform_StoR(s.reshape((1,3)), box)
    r_in_cell = distances.apply_PBC(r,box)

    # Now move vector s to the central cell
    s_in_cell = np.array([0.5, 0.9, 0.5], dtype=np.float32)
    v_in_cell = distances.transform_StoR(s_in_cell.reshape((1,3)), box)

    assert_almost_equal(r_in_cell, v_in_cell, decimal=5)

