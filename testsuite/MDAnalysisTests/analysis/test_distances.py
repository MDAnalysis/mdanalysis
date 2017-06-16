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

import scipy
import scipy.spatial

import MDAnalysis
from MDAnalysisTests import module_not_found
from MDAnalysisTests.datafiles import GRO
from MDAnalysisTests.util import block_import

import MDAnalysis.analysis.distances

from numpy.testing import TestCase, assert_equal, dec
import numpy as np

import warnings
import sys


class TestContactMatrix(TestCase):

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

class TestDist(TestCase):
    def setUp(self):
        self.u = MDAnalysis.Universe(GRO)
        self.ag = self.u.atoms[:20]
        self.u2 = MDAnalysis.Universe(GRO)
        self.ag2 = self.u2.atoms[:20]
        self.ag2.positions = np.random.shuffle(self.ag2.positions)
        self.expected = np.diag(scipy.spatial.distance.cdist(
                                                self.ag.positions,
                                                self.ag2.positions))

    def tearDown(self):
        del self.u
        del self.ag
        del self.u2
        del self.ag2
        del self.expected

    def test_pairwise_dist(self):
        '''Ensure that pairwise distances between atoms are
        correctly calculated.'''
        actual = MDAnalysis.analysis.distances.dist(self.ag, self.ag2)[2]
        assert_equal(actual, self.expected)

    def test_pairwise_dist_offset_effect(self):
        '''Test that feeding in offsets to dist() doesn't alter
        pairwise distance matrix.'''
        actual = MDAnalysis.analysis.distances.dist(self.ag, self.ag2,
                                                    offset=229)[2]
        assert_equal(actual, self.expected)


    def test_offset_calculation(self):
        '''Test that offsets fed to dist() are correctly calculated.'''
        actual = MDAnalysis.analysis.distances.dist(self.ag, self.ag2,
                offset=33)[:2]
        assert_equal(actual, np.array([self.ag.atoms.resids + 33,
                                       self.ag2.atoms.resids + 33]))

    def test_mismatch_exception(self):
        '''A ValueError should be raised if the two atomgroups
        don't have the same number of atoms.'''
        with self.assertRaises(ValueError):
            MDAnalysis.analysis.distances.dist(self.ag[:19], self.ag2)

class TestBetween(TestCase):
    def setUp(self):
        self.u = MDAnalysis.Universe(GRO)
        self.ag = self.u.atoms[:10]
        self.ag2 = self.u.atoms[12:33]
        self.group = self.u.atoms[40:]
        self.distance = 5.9
        self.distance_matrix_1 = scipy.spatial.distance.cdist(self.group.positions,
                                                              self.ag.positions)
        self.mask_1 = np.unique(np.where(self.distance_matrix_1 <= self.distance)[0])
        self.group_filtered = self.group[self.mask_1]
        self.distance_matrix_2 = scipy.spatial.distance.cdist(self.group_filtered.positions,
                                                              self.ag2.positions)
        self.mask_2 = np.unique(np.where(self.distance_matrix_2 <= self.distance)[0])
        self.expected = self.group_filtered[self.mask_2].indices

    def tearDown(self):
        del self.u
        del self.ag
        del self.ag2
        del self.group
        del self.distance
        del self.distance_matrix_1
        del self.distance_matrix_2
        del self.mask_1
        del self.mask_2
        del self.group_filtered
        del self.expected

    def test_between_simple_case_indices_only(self):
        '''Test MDAnalysis.analysis.distances.between() for
        a simple input case. Checks the sorted atom indices
        of returned AtomGroup against sorted expected index
        values.'''
        actual = sorted(MDAnalysis.analysis.distances.between(self.group,
                                                              self.ag,
                                                              self.ag2,
                                                              self.distance).indices)
        assert_equal(actual, self.expected)
