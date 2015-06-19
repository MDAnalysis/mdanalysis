# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

import operator
import numpy as np
from numpy.testing import assert_allclose, assert_equal, TestCase
import numpy.testing.utils
import MDAnalysis.lib.KDTree.NeighborSearch as KDNS

def assert_array_lessequal(x, y, err_msg='', verbose=True):
    """Raise an assertion if two array_like objects are not ordered by less-equal than."""
    np.testing.utils.assert_array_compare(
        operator.__le__, x, y, err_msg=err_msg,
        verbose=True,
        header='Arrays are not less-equal-ordered')

def assert_array_greater(x, y, err_msg='', verbose=True):
    """Raise an assertion if two array_like objects are not ordered by greater than."""
    np.testing.utils.assert_array_compare(
        operator.__le__, x, y, err_msg=err_msg,
        verbose=True,
        header='Arrays are not greater-ordered')

class TestKDTree(TestCase):
    def setUp(self):
        self.coords = np.array([
                [0, 0, 0],
                [1, 1, -1], [-1, 1, -1], [-1, -1, -1], [1, -1, -1],
                [1, 1, 0], [-1, 1, 0], [-1, -1, 0], [1, -1, 0],
                [1, 1, 1], [-1, 1, 1], [-1, -1, 1], [1, -1, 1],
                ], dtype=np.float32)

    def _test_point_neigborsearch(self, point, R):
        """Find points of a 2x2x2 cube (+origin) within R of point."""
        CNS = KDNS.CoordinateNeighborSearch(self.coords)
        center = np.asarray(point)

        def _assert_distances(assertfunc, indices, msg, center=center, R=R):
            diff = self.coords[found_indices] - center[np.newaxis, :]
            distances = np.sqrt(np.sum(diff * diff, axis=1))
            ref = R * np.ones_like(distances)
            assertfunc(distances, ref, err_msg=msg.format(point, R))

        found_indices = CNS.search(center, R)
        _assert_distances(assert_array_lessequal, found_indices,
                          "point {0} not found within distance {1} of the test coordinates.")

        all_indices = set(range(len(self.coords)))
        notfound_indices = all_indices.difference(found_indices)
        _assert_distances(assert_array_greater, notfound_indices,
                          "point {0} wrongly found within distance {1} of the test coordinates.")

    def test_CoordinateNeighborSearch(self, repeats=20):
        radii = np.random.uniform(0, 1.5, 10)
        points = np.random.uniform(-3,3, 3*repeats).reshape(repeats, 3)
        for radius in radii:
            for point in points:
                self._test_point_neigborsearch(point, radius)

    # add more/better tests

def test_KDTree_old_module():
    """test that MDAnalysis.KDTree is still importable (deprecated for 1.0)"""
    try:
        import MDAnalysis.KDTree.NeighborSearch as KDNS
    except (ImportError, NameError):
        raise AssertionError("MDAnalysis.KDTree not importable. Only remove for 1.0")

    # NOTE: removed this test with release 1.0 when we remove the stub


