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
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
import pytest
import scipy
import scipy.spatial

import MDAnalysis
from MDAnalysisTests.datafiles import GRO

import MDAnalysis.analysis.distances

from numpy.testing import (assert_equal ,assert_allclose)
import numpy as np


class TestContactMatrix(object):
    @staticmethod
    @pytest.fixture()
    def coord():
        return np.array([[1, 1, 1],
                      [5, 5, 5],
                      [1.1, 1.1, 1.1],
                      [11, 11, 11],  # neighboring image with pbc
                      [21, 21, 21]],  # non neighboring image with pbc
                     dtype=np.float32)
    
    @staticmethod
    @pytest.fixture()
    def box():
        return np.array([10, 10, 10, 90, 90, 90], dtype=np.float32)
    
    @staticmethod
    @pytest.fixture()
    def shape():
        return 5, 5
    
    @staticmethod
    @pytest.fixture()
    def res_no_pbc():
        return np.array([[1, 0, 1, 0, 0],
                           [0, 1, 0, 0, 0],
                           [1, 0, 1, 0, 0],
                           [0, 0, 0, 1, 0],
                           [0, 0, 0, 0, 1]], dtype=bool)
    
    @staticmethod
    @pytest.fixture()
    def res_pbc():
        return np.array([[1, 0, 1, 1, 1],
                        [0, 1, 0, 0, 0],
                        [1, 0, 1, 1, 1],
                        [1, 0, 1, 1, 1],
                        [1, 0, 1, 1, 1]], dtype=bool)

    def test_np(self, coord, shape, res_no_pbc):
        contacts = MDAnalysis.analysis.distances.contact_matrix(
            coord, cutoff=1, returntype="numpy"
        )
        assert contacts.shape == shape, \
            "wrong shape (should be {0})".format(shape)
        assert_equal(contacts, res_no_pbc)

    def test_sparse(self, coord, shape, res_no_pbc):
        contacts = MDAnalysis.analysis.distances.contact_matrix(
            coord, cutoff=1.5, returntype="sparse"
        )
        assert contacts.shape == shape, \
            "wrong shape (should be {0})".format(shape)
        assert_equal(contacts.toarray(), res_no_pbc)

    def test_box_numpy(self, coord, box, shape, res_pbc):
        contacts = MDAnalysis.analysis.distances.contact_matrix(
            coord, box=box, cutoff=1
        )
        assert contacts.shape == shape, \
            "wrong shape (should be {0})".format(shape)
        assert_equal(contacts, res_pbc)

    def test_box_sparse(self, coord, box, shape, res_pbc):
        contacts = MDAnalysis.analysis.distances.contact_matrix(
            coord, box=box, cutoff=1, returntype='sparse'
        )
        assert contacts.shape == shape, \
            "wrong shape (should be {0})".format(shape)
        assert_equal(contacts.toarray(), res_pbc)



class TestDist(object):

    @staticmethod
    @pytest.fixture()
    def ag():
        u = MDAnalysis.Universe(GRO)
        return u.atoms[1:10]

    # TODO: How are ag and ag2 different?!
    @staticmethod
    @pytest.fixture()
    def ag2():
        u2 = MDAnalysis.Universe(GRO)
        return u2.atoms[11:20]

    @staticmethod
    @pytest.fixture()
    def box():
        return np.array([8, 8, 8, 90, 90, 90], dtype=np.float32)
    
    @staticmethod
    @pytest.fixture()
    def expected(ag, ag2):
        
        return np.diag(scipy.spatial.distance.cdist(
            ag.positions, ag2.positions)
        )

    @staticmethod
    @pytest.fixture()
    def expected_box(ag, ag2, box):

        rp = np.abs(ag.positions - ag2.positions)
        box_2d = box[np.newaxis, 0:3]
        rp = np.where(rp  > box_2d / 2, box_2d - rp, rp)
        return np.sqrt(np.square(rp).sum(axis=1))

    def test_pairwise_dist(self, ag, ag2, expected):
        '''Ensure that pairwise distances between atoms are
        correctly calculated.'''
        actual = MDAnalysis.analysis.distances.dist(ag, ag2)[2]
        assert_allclose(actual, expected)

    def test_pairwise_dist_box(self, ag, ag2, expected_box, box):
        '''Ensure that pairwise distances between atoms are
        correctly calculated.'''
        actual = MDAnalysis.analysis.distances.dist(ag, ag2, 0, box)[2]
        assert_allclose(actual, expected_box, rtol=1e-05, atol=10)

    def test_pairwise_dist_offset_effect(self, ag, ag2, expected):
        '''Test that feeding in offsets to dist() doesn't alter
        pairwise distance matrix.'''
        actual = MDAnalysis.analysis.distances.dist(
            ag, ag2, offset=229)[2]
        assert_allclose(actual, expected)

    def test_offset_calculation(self, ag, ag2):
        '''Test that offsets fed to dist() are correctly calculated.'''
        actual = MDAnalysis.analysis.distances.dist(ag, ag2,
                                                    offset=33)[:2]
        assert_equal(actual, np.array([ag.atoms.resids + 33,
                                       ag2.atoms.resids + 33]))

    def test_mismatch_exception(self, ag, ag2, expected):
        '''A ValueError should be raised if the two atomgroups
        don't have the same number of atoms.'''
        with pytest.raises(ValueError):
            MDAnalysis.analysis.distances.dist(ag[:8], ag2)


class TestBetween(object):
    @staticmethod
    @pytest.fixture()
    def u():
        return MDAnalysis.Universe(GRO)

    @staticmethod
    @pytest.fixture()
    def ag(u):
        return u.atoms[:10]

    @staticmethod
    @pytest.fixture()
    def ag2(u):
        return u.atoms[12:33]

    @staticmethod
    @pytest.fixture()
    def group(u):
        return u.atoms[40:]

    distance = 5.9

    @pytest.fixture()
    def expected(self, group, ag, ag2):
        distance_matrix_1 = scipy.spatial.distance.cdist(group.positions,
                                                         ag.positions)
        mask_1 = np.unique(np.where(distance_matrix_1 <= self.distance)[0])
        group_filtered = group[mask_1]
        distance_matrix_2 = scipy.spatial.distance.cdist(group_filtered.positions,
                                                         ag2.positions)
        mask_2 = np.unique(np.where(distance_matrix_2 <= self.distance)[0])
        return group_filtered[mask_2].indices

    def test_between_simple_case_indices_only(self, group, ag, ag2, expected):
        '''Test MDAnalysis.analysis.distances.between() for
        a simple input case. Checks atom indices
        of returned AtomGroup against sorted expected index
        values.'''
        actual = MDAnalysis.analysis.distances.between(
            group,
            ag,
            ag2,
            self.distance
        ).indices
        assert_equal(actual, expected)

    @pytest.mark.parametrize('dists', [5.9, 0.0])
    def test_between_return_type(self, dists, group, ag, ag2):
        '''Test that MDAnalysis.analysis.distances.between() 
        returns an AtomGroup even when the returned group is empty.'''
        actual = MDAnalysis.analysis.distances.between(
            group,
            ag,
            ag2,
            dists
        )
        assert isinstance(actual, MDAnalysis.core.groups.AtomGroup)
