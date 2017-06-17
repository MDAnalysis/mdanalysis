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

from MDAnalysisTests.topology.base import ParserBase
from MDAnalysisTests.datafiles import (
    HoomdXMLdata,
)


class TestHoomdXMLParser(ParserBase):
    parser = mda.topology.HoomdXMLParser.HoomdXMLParser
    filename = HoomdXMLdata
    expected_attrs = ['types', 'masses', 'charges', 'radii',
                      'bonds', 'angles', 'dihedrals']
    expected_n_atoms = 769
    expected_n_residues = 1
    expected_n_segments = 1

    def test_attr_size(self):
        assert_(len(self.top.types) == self.top.n_atoms)
        assert_(len(self.top.charges) == self.top.n_atoms)
        assert_(len(self.top.masses) == self.top.n_atoms)

    def test_bonds(self):
        assert_(len(self.top.bonds.values) == 704)

    def test_angles(self):
        assert_(len(self.top.angles.values) == 640)

    def test_dihedrals(self):
        assert_(len(self.top.dihedrals.values) == 576)
