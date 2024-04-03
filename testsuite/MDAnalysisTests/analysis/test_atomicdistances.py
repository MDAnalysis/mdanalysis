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

import MDAnalysis as mda

import MDAnalysis.analysis.atomicdistances as ad
from MDAnalysis.lib.distances import calc_bonds
import MDAnalysis.transformations.boxdimensions as bd

from numpy.testing import assert_allclose
import numpy as np


class TestAtomicDistances(object):
    @staticmethod
    @pytest.fixture()
    def ad_u():
        # Universe from scratch with 100 atoms
        natoms = 100
        u = mda.Universe.empty(natoms, trajectory=True)

        # randomly generate 10 frames of x, y, z coordinates of 0 to 100
        rng = np.random.default_rng()
        coords = rng.integers(101, size=(10, natoms, 3))

        # load frames and coordinates into Universe
        u.load_new(coords)

        # set box dimensions to dim of [lx, ly, lz, alpha, beta, gamma]
        dim = np.array([80, 80, 80, 60, 60, 90])
        transform = bd.set_dimensions(dim)
        u.trajectory.add_transformations(transform)
        return u

    @staticmethod
    @pytest.fixture()
    def ad_ag1(ad_u):
        return ad_u.atoms[10:20]

    @staticmethod
    @pytest.fixture()
    def ad_ag2(ad_u):
        return ad_u.atoms[70:80]

    @staticmethod
    @pytest.fixture()
    def ad_ag3(ad_u):
        # more atoms than expected (exception)
        return ad_u.atoms[:7]

    @staticmethod
    @pytest.fixture()
    def ad_ag4():
        u_other = mda.Universe.empty(20, trajectory=True)
        # different trajectory (exception)
        return u_other.atoms[-10:]

    @staticmethod
    @pytest.fixture()
    def expected_dist(ad_ag1, ad_ag2):
        expected = np.zeros((len(ad_ag1.universe.trajectory),
                             ad_ag1.atoms.n_atoms))

        # calculate distances without PBCs using dist()
        for i, ts in enumerate(ad_ag1.universe.trajectory):
            expected[i] = calc_bonds(ad_ag1, ad_ag2)
        return expected

    @staticmethod
    @pytest.fixture()
    def expected_pbc_dist(ad_ag1, ad_ag2):
        expected = np.zeros((len(ad_ag1.universe.trajectory),
                             ad_ag1.atoms.n_atoms))

        # calculate distances with PBCs using dist()
        for i, ts in enumerate(ad_ag1.universe.trajectory):
            expected[i] = calc_bonds(ad_ag1, ad_ag2, box=ad_ag1.dimensions)
        return expected

    def test_ad_exceptions(self, ad_ag1, ad_ag3, ad_ag4):
        '''Test exceptions raised when number of atoms do not
        match or AtomGroups come from different trajectories.'''

        # number of atoms do not match
        match_exp_atoms = "AtomGroups do not"
        with pytest.raises(ValueError, match=match_exp_atoms):
            ad_atoms = ad.AtomicDistances(ad_ag1, ad_ag3)

        # AtomGroups come from different trajectories
        match_exp_traj = "AtomGroups are not"
        with pytest.raises(ValueError, match=match_exp_traj):
            ad_traj = ad.AtomicDistances(ad_ag1, ad_ag4)

    # only need to test that this class correctly applies distance calcs
    # calc_bonds() is tested elsewhere
    def test_ad_pairwise_dist(self, ad_ag1, ad_ag2,
                              expected_dist):
        '''Ensure that pairwise distances between atoms are
        correctly calculated without PBCs.'''
        pairwise_no_pbc = (ad.AtomicDistances(ad_ag1, ad_ag2,
                                              pbc=False).run())
        actual = pairwise_no_pbc.results

        # compare with expected values from dist()
        assert_allclose(actual, expected_dist)

    def test_ad_pairwise_dist_pbc(self, ad_ag1, ad_ag2,
                                  expected_pbc_dist):
        '''Ensure that pairwise distances between atoms are
        correctly calculated with PBCs.'''
        pairwise_pbc = (ad.AtomicDistances(ad_ag1, ad_ag2).run())
        actual = pairwise_pbc.results

        # compare with expected values from dist()
        assert_allclose(actual, expected_pbc_dist)
