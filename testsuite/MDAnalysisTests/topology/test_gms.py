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
from numpy.testing import assert_equal, assert_allclose

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
    expected_n_atoms = 6
    ref_filename = GMS_ASYMOPT

    def test_names(self, top):
        assert_equal(top.names.values,
                     ['O', 'H', 'H', 'O', 'H', 'H'])

    def test_types(self, top):
        assert_equal(top.atomiccharges.values,
                     [8, 1, 1, 8, 1, 1])

    def test_guessed_masses(self, filename):
        u = mda.Universe(filename)
        expected = [15.999, 1.008, 1.008, 15.999, 1.008, 1.008]
        assert_allclose(u.atoms.masses, expected)

    def test_guessed_types(self, filename):
        u = mda.Universe(filename)
        expected = ['O', 'H', 'H', 'O', 'H', 'H']
        assert (u.atoms.types == expected).all()


class TestGMSSYMOPT(GMSBase):
    expected_n_atoms = 4
    ref_filename = GMS_SYMOPT

    def test_names(self, top):
        assert_equal(top.names.values,
                     ['CARBON', 'CARBON', 'HYDROGEN', 'HYDROGEN'])

    def test_types(self, top):
        assert_equal(top.atomiccharges.values,
                     [6, 6, 1, 1])


class TestGMSASYMSURF(TestGMSASYMOPT):
    ref_filename = GMS_ASYMSURF
