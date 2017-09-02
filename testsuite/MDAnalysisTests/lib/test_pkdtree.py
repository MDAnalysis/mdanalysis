# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.mdanalysis.org
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

import numpy as np
from numpy.testing import TestCase, assert_array_equal

from MDAnalysis.lib.pkdtree import PeriodicKDTree


class TestPeriodicKDTree(TestCase):
    def setUp(self):
        self.box = np.array([10, 10, 10], dtype=np.float32)
        self.coords = np.array([[2, 2, 2],
                                [5, 5, 5],
                                [1.1, 1.1, 1.1],
                                [11, 11, 11],  # wrapped to [1, 1, 1]
                                [21, 21, 3]],  # wrapped to [1, 1, 3]
                                dtype=np.float32)

    def tearDown(self):
        pass

    def test_init(self):
        box = np.array([10, np.inf, -3], dtype=np.float32)
        tree = PeriodicKDTree(box)
        assert_array_equal(tree.box,
                           np.array([10, 0, 0], dtype=np.float32))

    def test_set_coords(self):
        tree = PeriodicKDTree(self.box)
        tree.set_coords(self.coords)
        print('hello')

    def test_find_centers(self):
        pass

    def test_search(self):
        pass

if __name__ == "__main__" :
    import sys
    np.testing.run_module_suite(argv=sys.argv)