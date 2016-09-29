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
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from numpy.testing import (
    assert_,
)

import MDAnalysis as mda
from MDAnalysis.core.topology import Topology


class ParserBase(object):
    """Base class for testing Topology parsers.

    All Parsers must subclass this class!
    """
    def setUp(self):
        with self.parser(self.filename) as p:
            self.top = p.parse()

    def tearDown(self):
        del self.top

    def test_output(self):
        """Testing the call signature"""
        with self.parser(self.filename) as p:
            top = p.parse()

        assert_(isinstance(top, Topology))

    def test_expected_attributes(self):
        for attr in self.expected_attrs:
            assert_(hasattr(self.top, attr),
                    'Missing attribute: {}'.format(attr))

    def test_size(self):
        """Check that the Topology is correctly sized"""
        assert_(self.top.n_atoms == self.expected_n_atoms,
                '{} atoms read, {} expected in {}'
                .format(self.top.n_atoms, self.expected_n_atoms,
                        self.__class__.__name__))
        assert_(self.top.n_residues == self.expected_n_residues,
                '{} residues read, {} expected in {}'
                .format(self.top.n_residues, self.expected_n_residues,
                        self.__class__.__name__))
        assert_(self.top.n_segments == self.expected_n_segments,
                '{} segment read, {} expected in {}'
                .format(self.top.n_segments, self.expected_n_segments,
                        self.__class__.__name__))

    def test_tt_size(self):
        """Check that the transtable is appropriately sized"""
        assert_(self.top.tt.size == (self.expected_n_atoms,
                                     self.expected_n_residues,
                                     self.expected_n_segments))

    def test_creates_universe(self):
        """Check that Universe works with this Parser"""
        u = mda.Universe(self.filename)
        assert_(isinstance(u, mda.Universe))
