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
import pytest

from MDAnalysisTests.topology.base import ParserBase
from MDAnalysisTests.datafiles import (
    XYZ,
    XYZ_bz2,
    XYZ_mini,
)

from numpy.testing import assert_equal, assert_allclose


class XYZBase(ParserBase):
    parser = mda.topology.XYZParser.XYZParser
    expected_n_residues = 1
    expected_n_segments = 1
    expected_attrs = ['names', 'elements']
    guessed_attrs = ['masses', 'types']


class TestXYZMini(XYZBase):
    ref_filename = XYZ_mini
    expected_n_atoms = 3


class TestXYZParser(XYZBase):
    expected_n_atoms = 1284

    @pytest.fixture(params=[XYZ, XYZ_bz2])
    def filename(self, request):
        return request.param

    def test_guessed_masses(self, filename):
        u = mda.Universe(filename)
        expected = [1.008, 1.008, 1.008, 1.008, 1.008, 1.008, 1.008]
        assert_allclose(u.atoms.masses[:7], expected)

    def test_guessed_types(self, filename):
        u = mda.Universe(filename)
        expected = ['H', 'H', 'H', 'H', 'H', 'H', 'H']
        assert_equal(u.atoms.types[:7], expected)
