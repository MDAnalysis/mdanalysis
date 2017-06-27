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
    PRM,  # ache.prmtop
    PRM12,  # anti.top
    PRM7,  # tz2.truncoct.parm7.bz2
    PRMpbc,
)


class TOPBase(ParserBase):
    parser = mda.topology.TOPParser.TOPParser
    expected_attrs = ["names", "types", "type_indices", "charges", "masses",
                      "resnames"]
    expected_n_segments = 1

    def test_attr_size(self):
        assert_(len(self.top.names) == self.expected_n_atoms)
        assert_(len(self.top.types) == self.expected_n_atoms)
        assert_(len(self.top.type_indices) == self.expected_n_atoms)
        assert_(len(self.top.charges) == self.expected_n_atoms)
        assert_(len(self.top.masses) == self.expected_n_atoms)
        assert_(len(self.top.resnames) == self.expected_n_residues)


class TestPRMParser(TOPBase):
    filename = PRM
    expected_n_atoms = 252
    expected_n_residues = 14
    guessed_attrs = ['elements']


class TestPRM12Parser(TOPBase):
    expected_attrs = ["names", "types", "type_indices", "charges", "masses",
                      "resnames"]
    filename = PRM12
    expected_n_atoms = 8923
    expected_n_residues = 2861
    ref_proteinatoms = 0


class TestParm7Parser(TOPBase):
    filename = PRM7
    expected_n_atoms = 5827
    expected_n_residues = 1882
    guessed_attrs = ['elements']


class TestPRM2(TOPBase):
    filename = PRMpbc
    expected_n_atoms = 5071
    expected_n_residues = 1686
    ref_proteinatoms = 22
    guessed_attrs = ['elements']
