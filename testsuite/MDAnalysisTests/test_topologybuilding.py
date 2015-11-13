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
"""Tests for building systems from scratch

Tests:
 - DummyParser
 - DummyReader
 - Universe.create_dummy()
"""
from numpy.testing import (
    assert_,
    assert_raises,
)

import MDAnalysis as mda


class TestDummyParser(object):
    def test_dummy(self):
        with mda.topology.dummy.DummyParser(10) as p:
            s = p.parse()
        assert_(type(s) == type(dict()))
        assert_(len(s['atoms']) == 10)

    def test_create_dummy(self):
        u = mda.Universe(10, topology_format='dummy')

        assert_(len(u.atoms) == 10)
        # Check unique indices
        assert_(len(set(u.atoms.indices)) == 10)


class TestDummyReader(object):
    def test_dummyreader(self):
        r = mda.coordinates.dummy.DummyReader(10)

        assert_(r.filename == 10)
        assert_(r.n_atoms == 10)

    def test_dummy_ts(self):
        r = mda.coordinates.dummy.DummyReader(10)

        assert_(len(r.ts) == 10)
        assert_(r.ts.has_positions == False)
        assert_raises(mda.NoDataError, getattr, r.ts, 'positions')
        assert_(r.ts.has_velocities == False)
        assert_raises(mda.NoDataError, getattr, r.ts, 'velocities')
        assert_(r.ts.has_forces == False)
        assert_raises(mda.NoDataError, getattr, r.ts, 'forces')


class TestDummyUniverse(object):
    def test_create_dummy(self):
        u = mda.Universe.create_dummy(10)

        assert_(len(u.atoms) == 10)
        assert_(u.atoms.universe is u)
        assert_raises(mda.NoDataError, getattr, u.atoms, 'positions')
