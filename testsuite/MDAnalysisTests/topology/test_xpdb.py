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
import MDAnalysis as mda

from MDAnalysisTests.topology.base import ParserBase
from MDAnalysisTests.datafiles import (
    XPDB_small,
)

from numpy.testing import assert_equal, assert_allclose


class TestXPDBParser(ParserBase):
    parser = mda.topology.ExtendedPDBParser.ExtendedPDBParser
    ref_filename = XPDB_small
    expected_attrs = ['ids', 'names', 'record_types', 'resids',
                      'resnames', 'altLocs', 'icodes', 'occupancies',
                      'tempfactors', 'chainIDs']
    guessed_attrs = ['masses', 'types']
    expected_n_atoms = 5
    expected_n_residues = 5
    expected_n_segments = 1

    def test_guessed_masses(self, filename):
        u = mda.Universe(filename)
        expected = [15.999, 15.999, 15.999, 15.999, 15.999]
        assert_allclose(u.atoms.masses, expected)

    def test_guessed_types(self, filename):
        u = mda.Universe(filename)
        expected = ['O', 'O', 'O', 'O', 'O']
        assert_equal(u.atoms.types, expected)
