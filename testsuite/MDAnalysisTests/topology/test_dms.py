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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import absolute_import

import MDAnalysis as mda

from MDAnalysisTests.topology.base import ParserBase
from MDAnalysisTests.datafiles import DMS


class TestDMSParser(ParserBase):
    parser = mda.topology.DMSParser.DMSParser
    ref_filename = DMS
    expected_attrs = ['ids', 'names', 'bonds', 'charges',
                      'masses', 'resids', 'resnames', 'segids',
                      'chainIDs', 'atomnums']
    guessed_attrs = ['types']
    expected_n_atoms = 3341
    expected_n_residues = 214
    expected_n_segments = 1

    def test_number_of_bonds(self, top):
        assert len(top.bonds.values) == 3365

    def test_atomsels(self, filename):
        # Desired value taken from VMD atomsel
        u = mda.Universe(filename)

        s0 = u.select_atoms("name CA")
        assert len(s0) == 214

        s1 = u.select_atoms("resid 33")
        assert len(s1) == 12

        s2 = u.select_atoms("segid 4AKE")
        assert len(s2) == 3341

        s3 = u.select_atoms("resname ALA")
        assert len(s3) == 190
