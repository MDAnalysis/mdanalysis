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
from numpy.testing import assert_equal,  assert_allclose

import MDAnalysis as mda
from MDAnalysis.core.topology import Topology
from MDAnalysis.guesser import DefaultGuesser

mandatory_attrs = ['ids', 'resids', 'resnums', 'segids']


class ParserBase(object):
    """Base class for testing Topology parsers.

    All Parsers must subclass this class!
    """

    expected_attrs = []
    guessed_attrs = []

    @pytest.fixture
    def filename(self):
        return self.ref_filename

    @pytest.fixture()
    def top(self, filename):
        with self.parser(filename) as p:
            yield p.parse()

    @pytest.fixture
    def guessed_types(self, top):
        return DefaultGuesser(None).guess_types(atoms=top.names.values)

    @pytest.fixture
    def guessed_masses(self, top):
        return DefaultGuesser(None).guess_masses(atoms=top.types.values)

    def test_output(self, filename):
        """Testing the call signature"""
        with self.parser(filename) as p:
            top = p.parse()

        assert isinstance(top, Topology)

    def test_mandatory_attributes(self, top):
        # attributes required as part of the API
        # ALL parsers must provide these
        for attr in mandatory_attrs:
            assert hasattr(top, attr), 'Missing required attribute: {}'.format(attr)

    def test_expected_attributes(self, top):
        # Extra attributes as declared in specific implementations
        for attr in self.expected_attrs:
            assert hasattr(top, attr), 'Missing expected attribute: {}'.format(attr)

    def test_no_unexpected_attributes(self, top):
        attrs = set(self.expected_attrs
                    + mandatory_attrs
                    + ['indices', 'resindices', 'segindices'] + self.guessed_attrs)
        for attr in top.attrs:
            assert attr.attrname in attrs, 'Unexpected attribute: {}'.format(attr.attrname)

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

    def test_creates_universe(self, filename):
        """Check that Universe works with this Parser"""
        u = mda.Universe(filename)
        assert isinstance(u, mda.Universe)

    def test_guessed_attributes(self, filename):
        """check that the universe created with certain parser have the same
        guessed attributes as  when it was guessed inside the parser"""
        u = mda.Universe(filename)
        for attr in self.guessed_attrs:
            assert hasattr(u.atoms, attr)

    @pytest.mark.skipif('names' not in expected_attrs,
                        reason="topology doesn't have names attribute")
    def test_guessed_types(self, filename, guessed_types):
        """check that type values from universe creation have the same
        expected values after removing mass and type guessing from parsers"""
        u = mda.Universe(filename)
        assert_equal(u.atoms.types, guessed_types)

    def test_guessed_masses(self, filename, guessed_masses):
        """check that mass values from universe creation have the same expected
        values after removing mass and type guessing from parsers"""
        u = mda.Universe(filename)
        assert_allclose(u.atoms.masses, guessed_masses, rtol=1e-3, atol=0)
