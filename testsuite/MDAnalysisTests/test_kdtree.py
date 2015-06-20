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
import numpy
from numpy.testing import assert_allclose, assert_equal, TestCase
import numpy.testing.utils
import MDAnalysis.lib.KDTree.NeighborSearch as KDNS

import warnings

from MDAnalysis.lib.KDTree import KDTree, CKDTree

def assert_array_lessequal(x, y, err_msg='', verbose=True):
    """Raise an assertion if two array_like objects are not ordered by less-equal than."""
    numpy.testing.utils.assert_array_compare(
        operator.__le__, x, y, err_msg=err_msg,
        verbose=True,
        header='Arrays are not less-equal-ordered')

def assert_array_greater(x, y, err_msg='', verbose=True):
    """Raise an assertion if two array_like objects are not ordered by greater than."""
    numpy.testing.utils.assert_array_compare(
        operator.__le__, x, y, err_msg=err_msg,
        verbose=True,
        header='Arrays are not greater-ordered')

class TestKDTree(TestCase):
    def _test_point_neighborsearch(self, point, R, nr_points=100, dim=3):
        """Find points of a 2x2x2 cube (+origin) within R of point."""

        coords = numpy.random.uniform(-3, 3, nr_points*dim).reshape(nr_points, dim).astype(numpy.float32)
        CNS = KDNS.CoordinateNeighborSearch(coords)
        center = numpy.asarray(point)

        def _assert_distances(assertfunc, indices, msg, center=center, R=R):
            distances = self._dist(coords[found_indices], center[numpy.newaxis, :])
            ref = R * numpy.ones_like(distances)
            assertfunc(distances, ref, err_msg=msg.format(point, R))

        found_indices = CNS.search(center, R)
        _assert_distances(assert_array_lessequal, found_indices,
                          "point {0} not found within distance {1} of the test coordinates.")

        all_indices = set(range(len(coords)))
        notfound_indices = all_indices.difference(found_indices)
        _assert_distances(assert_array_greater, notfound_indices,
                          "point {0} wrongly found within distance {1} of the test coordinates.")

    def test_CoordinateNeighborSearch(self, repeats=20):
        radii = numpy.random.uniform(0, 1.5, 10)
        points = numpy.random.uniform(-3,3, 3*repeats).reshape(repeats, 3)
        for radius in radii:
            for point in points:
                self._test_point_neighborsearch(point, radius)

    # from original KDTree.py

    def _dist(self, p, q):
        diff = p - q
        axis = 1 if len(diff.shape) == 2 else None
        return numpy.sqrt(numpy.sum(diff * diff, axis=axis))

    def test_neighbor_search(self, nr_points=100, dim=3, bucket_size=10, radius=0.05):
        """ Test all fixed radius neighbor search.

        Test all fixed radius neighbor search using the
        KD tree C module.

        o nr_points - number of points used in test
        o dim - dimension of coords
        o bucket_size - nr of points per tree node
        o radius - radius of search (typically 0.05 or so)
        """
        # KD tree search
        kdt = CKDTree.KDTree(dim, bucket_size)
        coords = numpy.random.random((nr_points, dim)).astype("f")
        kdt.set_data(coords, nr_points)
        kdt.neighbor_search(radius)
        r = kdt.neighbor_get_radii()
        if r is None:
            l1 = 0
        else:
            l1 = len(r)
        # now do a slow search to compare results
        kdt.neighbor_simple_search(radius)
        r = kdt.neighbor_get_radii()
        if r is None:
            l2 = 0
        else:
            l2 = len(r)
        assert_equal(l1, l2, err_msg="CKDTree neighbor_search() does not find"
                     "same number of points as manual search "
                     "{0} <> {1}.".format(l1, l2))


    def test_search_center_radius(self, nr_points=100, dim=3, bucket_size=10, radius=0.05):
        """Test neighbor search.

        Test neighbor search using the KD tree C module.

        o nr_points - number of points used in test
        o dim - dimension of coords
        o bucket_size - nr of points per tree node
        o radius - radius of search (typically 0.05 or so)
        """
        # kd tree search
        kdt = CKDTree.KDTree(dim, bucket_size)
        coords = numpy.random.random((nr_points, dim)).astype(numpy.float32)
        center = coords[0]
        kdt.set_data(coords, nr_points)
        kdt.search_center_radius(center, radius)
        r = kdt.get_indices()
        if r is None:
            l1 = 0
        else:
            l1 = len(r)
        l2 = 0
        # now do a manual search to compare results
        for i in range(0, nr_points):
            p = coords[i]
            if self._dist(p, center) <= radius:
                l2 = l2 + 1
        assert_equal(l1, l2, err_msg="CKDTree search_center_radius() does not find"
                     "same number of points as manual search "
                     "{0} <> {1}.".format(l1, l2))

    def test_all_search(self, nr_points=100, dim=3, bucket_size=10, query_radius=10):
        coords = 200 * numpy.random.random((nr_points, dim)).astype(numpy.float32)
        kdtree = KDTree(dim, bucket_size)
        # enter coords
        kdtree.set_coords(coords)

        # Find all point pairs within radius
        kdtree.all_search(query_radius)

        # get indices & radii of points

        # indices is a list of tuples. Each tuple contains the
        # two indices of a point pair within query_radius of
        # each other.
        indices = kdtree.all_get_indices()
        radii = kdtree.all_get_radii()

        #print "Found %i point pairs within radius %f." % (len(indices), query_radius)
        # Do 10 individual queries
        for _ in xrange(0, 10):
            # pick a random center
            center = numpy.random.random(dim).astype(numpy.float32)
            # search neighbors
            kdtree.search(center, query_radius)
            # get indices & radii of points
            indices = kdtree.get_indices()
            radii = kdtree.get_radii()
            # test that the points are really neighbors of center
            points = coords[indices]
            distances = self._dist(points, center[numpy.newaxis, :])
            ref_distances = query_radius * numpy.ones_like(distances)
            assert_array_lessequal(distances, ref_distances,
                                    err_msg="KDTree all_search() gives "
                                    "incorrect results for query_radius {0}".format(query_radius))

def test_KDTree_old_module():
    """test that MDAnalysis.KDTree is still importable (deprecated for 1.0)"""
    try:
        #import MDAnalysis.KDTree.NeighborSearch   # <--- does not work anymore
        import MDAnalysis.KDTree
        from MDAnalysis.KDTree import NeighborSearch
        from MDAnalysis.KDTree import KDTree
    except (ImportError, NameError):
        raise AssertionError("MDAnalysis.KDTree not importable. Only remove for 1.0")

    # NOTE: removed this test with release 1.0 when we remove the stub




