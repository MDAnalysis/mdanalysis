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

import pytest

import MDAnalysis as mda
from MDAnalysis.core.topology import Topology
from MDAnalysisTests.datafiles import CRD


@pytest.mark.parametrize('parser, filename, expected_attrs, guessed_attrs, expected_n_atoms, expected_n_residues, expected_n_segments',  [
    (mda.topology.CRDParser.CRDParser,
        CRD,
        ['ids', 'names', 'tempfactors', 'resids', 'resnames', 'resnums', 'segids'],
        ['masses', 'types'],
        3341,
        214,
        1),
])
class TestParserBase(object):
    """Base class for testing Topology parsers.

    All Parsers must subclass this class!
    """

    # expected_attrs = []
    # guessed_attrs = []

    @pytest.fixture()
    def top(self, parser, filename):
        with parser(filename) as p:
            yield p.parse()

    def test_output(self, parser, filename, expected_attrs, guessed_attrs, expected_n_atoms, expected_n_residues, expected_n_segments, top):
        """Testing the call signature"""
        with parser(filename) as p:
            top = p.parse()

        assert isinstance(top, Topology)

    def test_mandatory_attributes(self, parser, filename, expected_attrs, guessed_attrs, expected_n_atoms, expected_n_residues, expected_n_segments, top):
        # attributes required as part of the API
        # ALL parsers must provide these
        mandatory_attrs = ['ids', 'masses', 'types',
                           'resids', 'resnums', 'segids']

        for attr in mandatory_attrs:
            assert hasattr(top, attr), 'Missing required attribute: {}'.format(attr)

    def test_expected_attributes(self, parser, filename, expected_attrs, guessed_attrs, expected_n_atoms, expected_n_residues, expected_n_segments, top):
        # Extra attributes as declared in specific implementations
        for attr in expected_attrs:
            assert hasattr(top, attr), 'Missing expected attribute: {}'.format(attr)

    def test_guessed_attributes(self, parser, filename, expected_attrs, guessed_attrs, expected_n_atoms, expected_n_residues, expected_n_segments, top):
        # guessed attributes must be declared as guessed
        for attr in top.attrs:
            val = attr.is_guessed
            if not val in (True, False):  # only for simple yes/no cases
                continue
            assert val == (attr.attrname in guessed_attrs), 'Attr "{}" guessed= {}'.format(attr, val)

    def test_size(self, parser, filename, expected_attrs, guessed_attrs, expected_n_atoms, expected_n_residues, expected_n_segments, top):
        """Check that the Topology is correctly sized"""
        assert top.n_atoms == expected_n_atoms, '{} atoms read, {} expected in {}'.format(
            top.n_atoms, expected_n_atoms, self.__class__.__name__)

        assert top.n_residues == expected_n_residues, '{} residues read, {} expected in {}'.format(
            top.n_residues, expected_n_residues, self.__class__.__name__)

        assert top.n_segments == expected_n_segments, '{} segment read, {} expected in {}'.format(
            top.n_segments, expected_n_segments, self.__class__.__name__)

    def test_tt_size(self, parser, filename, expected_attrs, guessed_attrs, expected_n_atoms, expected_n_residues, expected_n_segments, top):
        """Check that the transtable is appropriately sized"""
        assert top.tt.size == (expected_n_atoms, expected_n_residues, expected_n_segments)

    def test_creates_universe(self, parser, filename, expected_attrs, guessed_attrs, expected_n_atoms, expected_n_residues, expected_n_segments):
        """Check that Universe works with this Parser"""
        u = mda.Universe(filename)
        assert isinstance(u, mda.Universe)


# TODO: Get rid of this class once all formats have been ported
class ParserBase(object):
    """Base class for testing Topology parsers.

    All Parsers must subclass this class!
    """

    expected_attrs = []
    guessed_attrs = []

    @pytest.fixture()
    def top(self):
        with self.parser(self.filename) as p:
            yield p.parse()

    def test_output(self):
        """Testing the call signature"""
        with self.parser(self.filename) as p:
            top = p.parse()

        assert isinstance(top, Topology)

    def test_mandatory_attributes(self, top):
        # attributes required as part of the API
        # ALL parsers must provide these
        mandatory_attrs = ['ids', 'masses', 'types',
                           'resids', 'resnums', 'segids']

        for attr in mandatory_attrs:
            assert hasattr(top, attr), 'Missing required attribute: {}'.format(attr)

    def test_expected_attributes(self, top):
        # Extra attributes as declared in specific implementations
        for attr in self.expected_attrs:
            assert hasattr(top, attr), 'Missing expected attribute: {}'.format(attr)

    def test_guessed_attributes(self, top):
        # guessed attributes must be declared as guessed
        for attr in top.attrs:
            val = attr.is_guessed
            if not val in (True, False):  # only for simple yes/no cases
                continue
            assert val == (attr.attrname in self.guessed_attrs), 'Attr "{}" guessed= {}'.format(attr, val)

    def test_size(self, top):
        """Check that the Topology is correctly sized"""
        assert top.n_atoms == self.expected_n_atoms, '{} atoms read, {} expected in {}'.format(
            top.n_atoms, self.expected_n_atoms, self.__class__.__name__)

        assert top.n_residues == self.expected_n_residues, '{} residues read, {} expected in {}'.format(
            top.n_residues, self.expected_n_residues, self.__class__.__name__)

        assert top.n_segments == self.expected_n_segments, '{} segment read, {} expected in {}'.format(
            top.n_segments, self.expected_n_segments, self.__class__.__name__)

    def test_tt_size(self, top):
        """Check that the transtable is appropriately sized"""
        assert top.tt.size == (self.expected_n_atoms, self.expected_n_residues, self.expected_n_segments)

    def test_creates_universe(self):
        """Check that Universe works with this Parser"""
        u = mda.Universe(self.filename)
        assert isinstance(u, mda.Universe)