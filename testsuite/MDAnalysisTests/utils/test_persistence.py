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
import pickle

import MDAnalysis as mda
from numpy.testing import assert_equal



from MDAnalysisTests.datafiles import PDB_small


class TestAtomGroupPickle(object):
    """Set up hopefully unique universes."""

    # _n marks named universes/atomgroups/pickled strings

    @staticmethod
    @pytest.fixture()
    def universe():
        return mda.Universe(PDB_small, PDB_small, PDB_small)

    @staticmethod
    @pytest.fixture()
    def ag(universe):
        return universe.atoms[:20]

    @staticmethod
    @pytest.fixture()
    def ag_2(universe):
        return universe.atoms[10:20]

    @staticmethod
    @pytest.fixture()
    def pickle_str(ag):
        return pickle.dumps(ag, protocol=pickle.HIGHEST_PROTOCOL)

    @staticmethod
    @pytest.fixture()
    def pickle_str_two_ag(ag, ag_2):
        return pickle.dumps((ag, ag_2), protocol=pickle.HIGHEST_PROTOCOL)

    @staticmethod
    @pytest.fixture()
    def pickle_str_ag_with_universe_f(ag, universe):
        return pickle.dumps((universe, ag), protocol=pickle.HIGHEST_PROTOCOL)

    @staticmethod
    @pytest.fixture()
    def pickle_str_ag_with_universe(ag, universe):
        return pickle.dumps((ag, universe), protocol=pickle.HIGHEST_PROTOCOL)

    def test_unpickle(self, pickle_str, ag, universe):
        """Test that an AtomGroup can be unpickled (Issue 293)"""
        newag = pickle.loads(pickle_str)
        # Can unpickle
        assert_equal(ag.indices, newag.indices)

    def test_pickle_unpickle_empty(self, universe):
        """Test that an empty AtomGroup can be pickled/unpickled (Issue 293)"""
        ag = universe.atoms[[]]
        pickle_str = pickle.dumps(ag, protocol=pickle.HIGHEST_PROTOCOL)
        newag = pickle.loads(pickle_str)
        assert len(newag) == 0

    def test_unpickle_two_ag(self, pickle_str_two_ag):
        newag, newag2 = pickle.loads(pickle_str_two_ag)
        assert newag.universe is newag2.universe, (
            "Two AtomGroups are unpickled to two different Universes"
        )

    def test_unpickle_ag_with_universe_f(self,
                                         pickle_str_ag_with_universe_f):
        newu, newag = pickle.loads(pickle_str_ag_with_universe_f)
        assert newag.universe is newu, (
            "AtomGroup is not unpickled to the bound Universe"
            "when Universe is pickled first"
        )

    def test_unpickle_ag_with_universe(self,
                                       pickle_str_ag_with_universe):
        newag, newu = pickle.loads(pickle_str_ag_with_universe)
        assert newag.universe is newu, (
                "AtomGroup is not unpickled to the bound Universe"
                "when AtomGroup is pickled first"
        )


class TestPicklingUpdatingAtomGroups(object):

    @staticmethod
    @pytest.fixture()
    def u():
        return mda.Universe(PDB_small)

    def test_pickling_uag(self, u):
        ag = u.atoms[:100]
        uag = ag.select_atoms('name C', updating=True)
        pickle_str = pickle.dumps(uag, protocol=pickle.HIGHEST_PROTOCOL)
        new_uag = pickle.loads(pickle_str)

        assert_equal(uag.indices, new_uag.indices)

    def test_pickling_uag_of_uag(self, u):
        uag1 = u.select_atoms('name C or name H', updating=True)
        uag2 = uag1.select_atoms('name C', updating=True)
        pickle_str = pickle.dumps(uag2, protocol=pickle.HIGHEST_PROTOCOL)
        new_uag2 = pickle.loads(pickle_str)

        assert_equal(uag2.indices, new_uag2.indices)
