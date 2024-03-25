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
from numpy.testing import assert_almost_equal
import MDAnalysis as mda

from MDAnalysisTests.topology.base import ParserBase
from MDAnalysisTests.datafiles import HoomdXMLdata


class TestHoomdXMLParser(ParserBase):
    parser = mda.topology.HoomdXMLParser.HoomdXMLParser
    ref_filename = HoomdXMLdata
    expected_attrs = [
        'types', 'masses', 'charges', 'radii', 'bonds', 'angles', 'dihedrals', 'impropers'
    ]

    expected_n_atoms = 769
    expected_n_residues = 1
    expected_n_segments = 1

    def test_attr_size(self, top):
        assert len(top.types) == top.n_atoms
        assert len(top.charges) == top.n_atoms
        assert len(top.masses) == top.n_atoms

    def test_bonds(self, top):
        assert len(top.bonds.values) == 704
        assert isinstance(top.bonds.values[0], tuple)

    def test_angles(self, top):
        assert len(top.angles.values) == 640
        assert isinstance(top.angles.values[0], tuple)

    def test_dihedrals(self, top):
        assert len(top.dihedrals.values) == 576
        assert isinstance(top.dihedrals.values[0], tuple)

    def test_impropers(self, top):
        assert len(top.impropers.values) == 0

    def test_bonds_identity(self, top):
        vals = top.bonds.values
        for b in ((0, 1), (1, 2), (2, 3), (3, 4)):
            assert (b in vals) or (b[::-1] in vals)
        assert ((0, 450) not in vals)

    def test_angles_identity(self, top):
        vals = top.angles.values
        for b in ((0, 1, 2), (1, 2, 3), (2, 3, 4), (3, 4, 5)):
            assert (b in vals) or (b[::-1] in vals)
        assert ((0, 350, 450) not in vals)

    def test_dihedrals_identity(self, top):
        vals = top.dihedrals.values
        for b in ((0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5), (3, 4, 5, 6)):
            assert (b in vals) or (b[::-1] in vals)
        assert ((0, 250, 350, 450) not in vals)

    def test_read_masses(self, top):
        assert_almost_equal(top.masses.values, 1.0)

    def test_read_charges(self, top):
        # note: the example topology file contains 0 for all charges which
        # is the same as the default so this test does not fully test
        # reading of charges from the file (#2888)
        assert_almost_equal(top.charges.values, 0.0)
