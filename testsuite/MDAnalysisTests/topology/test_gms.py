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
    assert_array_equal,
)

import MDAnalysis as mda

from MDAnalysisTests.topology.base import ParserBase
from MDAnalysisTests.datafiles import (
    GMS_ASYMOPT,  # c1opt.gms.gz
    GMS_SYMOPT,  # symopt.gms
    GMS_ASYMSURF,  # surf2wat.gms
)


class GMSBase(ParserBase):
    parser = mda.topology.GMSParser.GMSParser
    expected_attrs = ['names', 'atomiccharges']
    guessed_attrs = ['masses', 'types']
    expected_n_residues = 1
    expected_n_segments = 1


class TestGMSASYMOPT(GMSBase):
    filename = GMS_ASYMOPT
    expected_n_atoms = 6

    def test_names(self):
        assert_array_equal(self.top.names.values,
                           ['O', 'H', 'H', 'O', 'H', 'H'])

    def test_types(self):
        assert_array_equal(self.top.atomiccharges.values,
                           [8, 1, 1, 8, 1, 1])


class TestGMSSYMOPT(GMSBase):
    filename = GMS_SYMOPT
    expected_n_atoms = 4

    def test_names(self):
        assert_array_equal(self.top.names.values,
                           ['CARBON', 'CARBON', 'HYDROGEN', 'HYDROGEN'])

    def test_types(self):
        assert_array_equal(self.top.atomiccharges.values,
                           [6, 6, 1, 1])

class TestGMSASYMSURF(TestGMSASYMOPT):
    filename = GMS_ASYMSURF
