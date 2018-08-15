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
from MDAnalysisTests.datafiles import (
    PRM,  # ache.prmtop
    PRM12,  # anti.top
    PRM7,  # tz2.truncoct.parm7.bz2
    PRMpbc,
)


class TOPBase(ParserBase):
    parser = mda.topology.TOPParser.TOPParser
    expected_attrs = [
        "names", "types", "type_indices", "charges", "masses", "resnames",
        "bonds", "angles", "dihedrals", "impropers"
    ]
    expected_n_segments = 1

    def test_attr_size(self, top):
        assert len(top.names) == self.expected_n_atoms
        assert len(top.types) == self.expected_n_atoms
        assert len(top.type_indices) == self.expected_n_atoms
        assert len(top.charges) == self.expected_n_atoms
        assert len(top.masses) == self.expected_n_atoms
        assert len(top.resnames) == self.expected_n_residues
        assert len(top.bonds.values) == self.expected_n_bonds
        assert len(top.angles.values) == self.expected_n_angles
        assert len(top.dihedrals.values) == self.expected_n_dihedrals
        assert len(top.impropers.values) == self.expected_n_impropers


class TestPRMParser(TOPBase):
    ref_filename = PRM
    expected_n_atoms = 252
    expected_n_residues = 14
    expected_n_bonds = 259
    expected_n_angles = 456
    expected_n_dihedrals = 673
    expected_n_impropers = 66
    guessed_attrs = ['elements']


class TestPRM12Parser(TOPBase):
    ref_filename = PRM12
    expected_attrs = [
        "names", "types", "type_indices", "charges", "masses", "resnames",
        "bonds", "angles", "dihedrals", "impropers"
    ]
    expected_n_atoms = 8923
    expected_n_residues = 2861
    expected_n_bonds = 8947
    expected_n_angles = 756
    expected_n_dihedrals = 1128
    expected_n_impropers = 72
    ref_proteinatoms = 0


class TestParm7Parser(TOPBase):
    ref_filename = PRM7
    expected_n_atoms = 5827
    expected_n_residues = 1882
    expected_n_bonds = 5834
    expected_n_angles = 402
    expected_n_dihedrals = 602
    expected_n_impropers = 55
    guessed_attrs = ['elements']


class TestPRM2(TOPBase):
    ref_filename = PRMpbc
    expected_n_atoms = 5071
    expected_n_residues = 1686
    ref_proteinatoms = 22
    expected_n_bonds = 5070
    expected_n_angles = 36
    expected_n_dihedrals = 41
    expected_n_impropers = 4
    guessed_attrs = ['elements']
