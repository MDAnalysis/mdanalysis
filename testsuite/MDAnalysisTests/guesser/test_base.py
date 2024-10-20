# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2 or any higher version
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
import numpy as np
import MDAnalysis as mda
from MDAnalysis.guesser.base import GuesserBase, get_guesser
from MDAnalysis.core.topology import Topology
from MDAnalysis.core.topologyattrs import Masses, Atomnames, Atomtypes
import MDAnalysis.tests.datafiles as datafiles
from numpy.testing import assert_allclose, assert_equal


class TesttBaseGuesser():

    def test_get_guesser(self):
        class TestGuesser1(GuesserBase):
            context = 'test1'

        class TestGuesser2(GuesserBase):
            context = 'test2'

        assert get_guesser(TestGuesser1).context == 'test1'
        assert get_guesser('test1').context == 'test1'
        assert get_guesser(TestGuesser2()).context == 'test2'

    def test_get_guesser_with_universe(self):
        class TestGuesser1(GuesserBase):
            context = 'test1'

        u = mda.Universe.empty(n_atoms=5)
        guesser = get_guesser(TestGuesser1(), u, foo=1)

        assert len(guesser._universe.atoms) == 5
        assert 'foo' in guesser._kwargs

    def test_guess_invalid_attribute(self):
        with pytest.raises(ValueError,
                           match='default guesser can not guess '
                                 'the following attribute: foo'):
            mda.Universe(datafiles.PDB, to_guess=['foo'])

    def test_guess_attribute_with_missing_parent_attr(self):
        names = Atomnames(np.array(['C', 'HB', 'HA', 'O'], dtype=object))
        masses = Masses(
            np.array([np.nan, np.nan, np.nan, np.nan], dtype=np.float64))
        top = Topology(4, 1, 1, attrs=[names, masses, ])
        u = mda.Universe(top, to_guess=['masses'])
        assert_allclose(u.atoms.masses, np.array(
            [12.01100, 1.00800, 1.00800, 15.99900]), atol=0)

    def test_force_guessing(self):
        names = Atomnames(np.array(['C', 'H', 'H', 'O'], dtype=object))
        types = Atomtypes(np.array(['1', '2', '3', '4'], dtype=object))
        top = Topology(4, 1, 1, attrs=[names, types, ])
        u = mda.Universe(top, force_guess=['types'])
        assert_equal(u.atoms.types, ['C', 'H', 'H', 'O'])

    def test_partial_guessing(self):
        types = Atomtypes(np.array(['C', 'H', 'H', 'O'], dtype=object))
        masses = Masses(np.array([0, np.nan, np.nan, 0], dtype=np.float64))
        top = Topology(4, 1, 1, attrs=[types, masses, ])
        u = mda.Universe(top, to_guess=['masses'])
        assert_allclose(u.atoms.masses, np.array(
            [0, 1.00800, 1.00800, 0]), atol=0)

    def test_force_guess_priority(self):
        "check that passing the attribute to force_guess have higher power"
        types = Atomtypes(np.array(['C', 'H', 'H', 'O'], dtype=object))
        masses = Masses(np.array([0, np.nan, np.nan, 0], dtype=np.float64))
        top = Topology(4, 1, 1, attrs=[types, masses, ])
        u = mda.Universe(top, to_guess=['masses'], force_guess=['masses'])
        assert_allclose(u.atoms.masses, np.array(
            [12.01100, 1.00800, 1.00800, 15.99900]), atol=0)

    def test_partial_guess_attr_with_unknown_no_value_label(self):
        "trying to partially guess attribute tha doesn't have declared"
        "no_value_label should gives no effect"
        names = Atomnames(np.array(['C', 'H', 'H', 'O'], dtype=object))
        types = Atomtypes(np.array(['', '', '', ''], dtype=object))
        top = Topology(4, 1, 1, attrs=[names, types, ])
        u = mda.Universe(top, to_guess=['types'])
        assert_equal(u.atoms.types, ['', '', '', ''])


def test_Universe_guess_bonds_deprecated():
    with pytest.warns(
        DeprecationWarning,
        match='`guess_bonds` keyword is deprecated'
    ):
        u = mda.Universe(datafiles.PDB_full, guess_bonds=True)
