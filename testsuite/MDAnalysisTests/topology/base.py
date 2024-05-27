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
from MDAnalysis.core.topology import Topology

mandatory_attrs = ['ids', 'masses', 'types', 
                   'resids', 'resnums', 'segids']


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
        for attr in self.expected_attrs+self.guessed_attrs:
            assert hasattr(top, attr), 'Missing expected attribute: {}'.format(attr)
    
    def test_no_unexpected_attributes(self, top):
        attrs = set(self.expected_attrs
                    + self.guessed_attrs
                    + mandatory_attrs
                    + ['indices', 'resindices', 'segindices'])
        for attr in top.attrs:
            assert attr.attrname in attrs, 'Unexpected attribute: {}'.format(attr.attrname)

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

    def test_creates_universe(self, filename):
        """Check that Universe works with this Parser"""
        u = mda.Universe(filename)
        assert isinstance(u, mda.Universe)

    def test_pathlib_input(self, filename):
        """Check that pathlib.Path objects are accepted"""
        import pathlib
        path = pathlib.Path(filename)
        u_str = mda.Universe(filename)
        u_path = mda.Universe(path)
        assert u_str.atoms.n_atoms == u_path.atoms.n_atoms
