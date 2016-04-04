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
from __future__ import print_function

import MDAnalysis
from MDAnalysisTests import module_not_found

from numpy.testing import TestCase, assert_equal, dec
import numpy as np


class TestContactMatrix(TestCase):

    @dec.skipif(module_not_found('scipy'),
                "Test skipped because scipy is not available.")
    def setUp(self):
        import MDAnalysis.analysis.distances
        self.coord = np.array([[1, 1, 1],
                                  [5, 5, 5],
                                  [1.1, 1.1, 1.1],
                                  [11, 11, 11],  # neighboring image with pbc
                                  [21, 21, 21]],  # non neighboring image with pbc
                                 dtype=np.float32)
        self.box = np.array([10, 10, 10], dtype=np.float32)
        self.shape = (5, 5)
        self.res_no_pbc = np.array([[1, 0, 1, 0, 0],
                                       [0, 1, 0, 0, 0],
                                       [1, 0, 1, 0, 0],
                                       [0, 0, 0, 1, 0],
                                       [0, 0, 0, 0, 1]], dtype=np.bool)
        self.res_pbc = np.array([[1, 0, 1, 1, 1],
                                    [0, 1, 0, 0, 0],
                                    [1, 0, 1, 1, 1],
                                    [1, 0, 1, 1, 1],
                                    [1, 0, 1, 1, 1]], dtype=np.bool)

    def test_np(self):
        contacts = MDAnalysis.analysis.distances.contact_matrix(
            self.coord, cutoff=1, returntype="numpy")
        assert_equal(contacts.shape, self.shape,
                     "wrong shape (should be {0})".format(self.shape))
        assert_equal(contacts, self.res_no_pbc)

    def test_sparse(self):
        contacts = MDAnalysis.analysis.distances.contact_matrix(
            self.coord, cutoff=1.5, returntype="sparse")
        assert_equal(contacts.shape, self.shape,
                     "wrong shape (should be {0})".format(self.shape))
        assert_equal(contacts.toarray(), self.res_no_pbc)

    def test_box_numpy(self):
        contacts = MDAnalysis.analysis.distances.contact_matrix(
            self.coord, box=self.box, cutoff=1)
        assert_equal(contacts.shape, self.shape,
                     "wrong shape (should be {0})".format(self.shape))
        assert_equal(contacts, self.res_pbc)

    def test_box_sparse(self):
        contacts = MDAnalysis.analysis.distances.contact_matrix(
            self.coord, box=self.box, cutoff=1, returntype='sparse')
        assert_equal(contacts.shape, self.shape,
                     "wrong shape (should be {0})".format(self.shape))
        assert_equal(contacts.toarray(), self.res_pbc)

