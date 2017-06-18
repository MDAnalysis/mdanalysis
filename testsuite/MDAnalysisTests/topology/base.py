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
from __future__ import absolute_import
from numpy.testing import (
    assert_,
)

import MDAnalysis as mda
from MDAnalysis.core.topology import Topology


class ParserBase(object):
    """Base class for testing Topology parsers.

    All Parsers must subclass this class!
    """
    expected_attrs = []
    guessed_attrs = []

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

    def test_mandatory_attributes(self):
        # attributes required as part of the API
        # ALL parsers must provide these
        mandatory_attrs = ['ids', 'masses', 'types',
                           'resids', 'resnums', 'segids']

        for attr in mandatory_attrs:
            assert_(hasattr(self.top, attr),
                    'Missing required attribute: {}'.format(attr))

    def test_expected_attributes(self):
        # Extra attributes as declared in specific implementations
        for attr in self.expected_attrs:
            assert_(hasattr(self.top, attr),
                    'Missing expected attribute: {}'.format(attr))

    def test_guessed_attributes(self):
        # guessed attributes must be declared as guessed
        for attr in self.top.attrs:
            val = attr.is_guessed
            if not val in (True, False):  # only for simple yes/no cases
                continue
            assert_(val == (attr.attrname in self.guessed_attrs),
                    'Attr "{}" guessed= {}'.format(attr, val))

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
