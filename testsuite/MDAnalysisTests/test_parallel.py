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
from MDAnalysis.lib.parallel.distances import distance_array, distance_array_serial
from MDAnalysis.lib.distances import distance_array as distance_array_reference

import numpy as np
from numpy.testing import *


try:
    from numpy.testing import assert_allclose
except ImportError:
    # not implemented in numpy 1.3 ... 1.4.1
    def assert_allclose(actual, desired, rtol=1e-7, atol=0, err_msg='', verbose=True):
        # from numpy/testing/utils.py (1.6.1)
        import numpy.testing.utils

        def compare(x, y):
            return np.allclose(x, y, rtol=rtol, atol=atol)

        actual, desired = np.asanyarray(actual), np.asanyarray(desired)
        header = 'Not equal to tolerance rtol=%g, atol=%g' % (rtol, atol)
        numpy.testing.utils.assert_array_compare(compare, actual, desired, err_msg=str(err_msg), verbose=verbose,
                                                 header=header)


class TestDistances(TestCase):
    """
    Test if cython (serial and parallel) distance matrices are the same as
    MDAnalysis.lib.distances.distance_array
    """

    def setUp(self):
        # Create two random sets of 3N coordinates, in the same format
        # as delivered by MDAnalysis coordinate readers
        self.coord = [
            np.random.random((100, 3)).astype(np.float32),
            np.random.random((100, 3)).astype(np.float32),
        ]
        self.ref = distance_array_reference(self.coord[0], self.coord[1])
        self.box = (np.random.random((3)) * np.random.random()).astype(np.float32)
        self.ref_pbc = distance_array_reference(self.coord[0], self.coord[1], box=self.box)

    def test_PBC2(self):
        a = np.array([7.90146923, -13.72858524, 3.75326586], dtype=np.float32)
        b = np.array([-1.36250901, 13.45423985, -0.36317623], dtype=np.float32)
        box = np.array([5.5457325, 5.5457325, 5.5457325], dtype=np.float32)

        def mindist(a, b, box):
            x = a - b
            return np.linalg.norm(x - np.rint(x / box) * box)

        ref = mindist(a, b, box)
        val = distance_array(np.array([a]), np.array([b]), box)[0, 0]
        assert_allclose(val, ref, rtol=1e-6, atol=1e-6,
                        err_msg="Issue 151 not correct (PBC in distance array)")

    def test_distance_array_parallel(self):
        cython_parallel = distance_array(self.coord[0], self.coord[1])
        assert_allclose(cython_parallel, self.ref, rtol=1e-6, atol=1e-6,
                        err_msg="Cython parallel distance matrix does not match C")

    def test_distance_array_parallel_results(self):
        result = np.empty((self.coord[0].shape[0], self.coord[0].shape[0])).astype(np.float32)
        cython_parallel = distance_array(self.coord[0], self.coord[1], result=result)
        assert_allclose(cython_parallel, self.ref, rtol=1e-6, atol=1e-6,
                        err_msg="Cython parallel distance matrix does not match C")

    def test_distance_array_parallel_pbc(self):
        cython_parallel = distance_array(self.coord[0], self.coord[1], box=self.box)
        assert_allclose(cython_parallel, self.ref_pbc, rtol=1e-6, atol=1e-6,
                        err_msg="Cython parallel distance matrix does not match C")

    def test_distance_array_serial(self):
        coord = self.coord
        cython_serial = distance_array_serial(self.coord[0], self.coord[1])
        assert_allclose(cython_serial, self.ref, rtol=1e-6, atol=1e-6,
                        err_msg="Cython serial distance matrix does not match C")

    def test_distance_array_serial_results(self):
        result = np.empty((self.coord[0].shape[0], self.coord[0].shape[0])).astype(np.float32)
        cython_parallel = distance_array_serial(self.coord[0], self.coord[1], result=result)
        assert_allclose(cython_parallel, self.ref, rtol=1e-6, atol=1e-6,
                        err_msg="Cython parallel distance matrix does not match C")

    def test_distance_array_serial_pbc(self):
        cython_parallel = distance_array_serial(self.coord[0], self.coord[1], box=self.box)
        assert_allclose(cython_parallel, self.ref_pbc, rtol=1e-6, atol=1e-6,
                        err_msg="Cython parallel distance matrix does not match C")
