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
from __future__ import absolute_import

import MDAnalysis as mda

from MDAnalysisTests.topology.base import ParserBase
from MDAnalysisTests.datafiles import HoomdXMLdata


class TestHoomdXMLParser(ParserBase):
    parser = mda.topology.HoomdXMLParser.HoomdXMLParser
    ref_filename = HoomdXMLdata
    expected_attrs = [
        'types', 'masses', 'charges', 'radii', 'bonds', 'angles', 'dihedrals'
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
