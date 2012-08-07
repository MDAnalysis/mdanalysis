# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2011 Naveen Michaud-Agrawal,
#               Elizabeth J. Denning, Oliver Beckstein,
#               and contributors (see website for details)
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
#     N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
#     O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
#     Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
#     doi:10.1002/jcc.21787
#
from MDAnalysis.core.parallel.distances import distance_array, distance_array_serial
from MDAnalysis.core.distances import distance_array as distance_array_reference

import numpy as np
from numpy.testing import *
from numpy.ma.testutils import assert_allclose

class TestDistances(TestCase):
    """
    Test if cython (serial and parallel) distance matrices are the same as
    MDAnalysis.core.distances.distance_array
    """
    def setUp(self):
        # Create two random sets of 3N coordinates, in the same format
        # as delivered by MDAnalysis coordinate readers
        self.coord = [
            np.random.random((100,3)).astype(np.float32),
            np.random.random((100,3)).astype(np.float32),
            ]
        self.ref = distance_array_reference(self.coord[0], self.coord[1])
        self.box = (np.random.random((3)) * np.random.random()).astype(np.float32)
        self.ref_pbc = distance_array_reference(self.coord[0], self.coord[1], box=self.box)
        

    def test_distance_array_parallel(self):
        cython_parallel = distance_array(self.coord[0], self.coord[1])
        assert_allclose(cython_parallel, self.ref, rtol=1e-6, atol=0,
                        err_msg="Cython parallel distance matrix does not match C")

    def test_distance_array_parallel_results(self):
        result = np.empty((self.coord[0].shape[0], self.coord[0].shape[0])).astype(np.float32)
        cython_parallel = distance_array(self.coord[0], self.coord[1], result=result)
        assert_allclose(cython_parallel, self.ref, rtol=1e-6, atol=0,
                        err_msg="Cython parallel distance matrix does not match C")

    def test_distance_array_parallel_pbc(self):
        cython_parallel = distance_array(self.coord[0], self.coord[1], box=self.box)
        assert_allclose(cython_parallel, self.ref_pbc, rtol=1e-6, atol=0,
                        err_msg="Cython parallel distance matrix does not match C")

    def test_distance_array_serial(self):
        coord = self.coord
        cython_serial = distance_array_serial(self.coord[0], self.coord[1])
        assert_allclose(cython_serial, self.ref, rtol=1e-6, atol=0,
                        err_msg="Cython serial distance matrix does not match C")

    def test_distance_array_serial_results(self):
        result = np.empty((self.coord[0].shape[0], self.coord[0].shape[0])).astype(np.float32)
        cython_parallel = distance_array_serial(self.coord[0], self.coord[1], result=result)
        assert_allclose(cython_parallel, self.ref, rtol=1e-6, atol=0,
                        err_msg="Cython parallel distance matrix does not match C")

    def test_distance_array_serial_pbc(self):
        cython_parallel = distance_array_serial(self.coord[0], self.coord[1], box=self.box)
        assert_allclose(cython_parallel, self.ref_pbc, rtol=1e-6, atol=0,
                        err_msg="Cython parallel distance matrix does not match C")
