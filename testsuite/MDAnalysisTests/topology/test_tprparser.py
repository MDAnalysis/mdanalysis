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

import pytest
from numpy.testing import assert_equal
import functools

import numpy as np

from MDAnalysis.tests.datafiles import (
    TPR,
    TPR400, TPR402, TPR403, TPR404, TPR405, TPR406, TPR407,
    TPR450, TPR451, TPR452, TPR453, TPR454, TPR455, TPR455Double,
    TPR460, TPR461, TPR502, TPR504, TPR505, TPR510, TPR510_bonded,
    TPR2016, TPR2016_bonded,
)
from MDAnalysisTests.topology.base import ParserBase
import MDAnalysis.topology.TPRParser


class TPRAttrs(ParserBase):
    parser = MDAnalysis.topology.TPRParser.TPRParser
    expected_attrs = ['ids', 'names', 'resids', 'resnames', 'moltypes']
    guessed_attrs = ['elements']

    def test_moltypes(self, top):
        moltypes = top.moltypes.values
        assert_equal(moltypes, self.ref_moltypes)


class TestTPR(TPRAttrs):
    """
    this test the data/adk_oplsaa.tpr which is of tpx version 58
    """
    expected_n_atoms = 47681
    expected_n_residues = 11302
    expected_n_segments = 3
    ref_moltypes = np.array(['AKeco'] * 214 + ['SOL'] * 11084 + ['NA+'] * 4,
                            dtype=object)

    @pytest.fixture()
    def filename(self):
        return TPR


# The follow test the same system grompped by different version of gromacs
# FORMAT: TPRABC, where numbers ABC indicates the version of gromacs that
# generates the corresponding tpr file

class TPRBase(TPRAttrs):
    expected_n_atoms = 2263
    expected_n_residues = 230
    expected_n_segments = 2
    ref_moltypes = np.array(['Protein_A'] * 129 + ['SOL'] * 101, dtype=object)


# All these classes should be generated in a loop.
class TestTPR400(TPRBase):
    @pytest.fixture()
    def filename(self):
        return TPR400


class TestTPR402(TPRBase):
    @pytest.fixture()
    def filename(self):
        return TPR400


class TestTPR403(TPRBase):
    @pytest.fixture()
    def filename(self):
        return TPR403


class TestTPR404(TPRBase):
    @pytest.fixture()
    def filename(self):
        return TPR404


class TestTPR405(TPRBase):
    @pytest.fixture()
    def filename(self):
        return TPR405


class TestTPR406(TPRBase):
    @pytest.fixture()
    def filename(self):
        return TPR406


class TestTPR407(TPRBase):
    @pytest.fixture()
    def filename(self):
        return TPR407


class TestTPR450(TPRBase):
    @pytest.fixture()
    def filename(self):
        return TPR450


class TestTPR451(TPRBase):
    @pytest.fixture()
    def filename(self):
        return TPR451


class TestTPR452(TPRBase):
    @pytest.fixture()
    def filename(self):
        return TPR452


class TestTPR453(TPRBase):
    @pytest.fixture()
    def filename(self):
        return TPR453


class TestTPR454(TPRBase):
    @pytest.fixture()
    def filename(self):
        return TPR454


class TestTPR455(TPRBase):
    @pytest.fixture()
    def filename(self):
        return TPR455


class TPRDouble(TPRAttrs):
    expected_n_atoms = 21692
    expected_n_residues = 4352
    expected_n_segments = 7
    ref_moltypes = np.array(['DOPC'] * 21 + ['DPPC'] * 10 + ['CHOL'] * 3
                            + ['DOPC'] * 21 + ['DPPC'] * 10 + ['CHOL'] * 3
                            + ['SOL'] * 4284,
                            dtype=object)


class TestTPR455Double(TPRDouble):
    @pytest.fixture()
    def filename(self):
        return TPR455Double


class TPR46xBase(TPRAttrs):
    expected_n_atoms = 44052
    expected_n_residues = 10712
    expected_n_segments = 8
    ref_moltypes = np.array(['Protein_A'] * 27 + ['Protein_B'] * 27
                            + ['Protein_C'] * 27 + ['Protein_D'] * 27
                            + ['Protein_E'] * 27
                            + ['SOL'] * 10530 + ['NA+'] * 26 + ['CL-'] * 21,
                            dtype=object)


class TestTPR460(TPR46xBase):
    @pytest.fixture()
    def filename(self):
        return TPR460


class TestTPR461(TPR46xBase):
    @pytest.fixture()
    def filename(self):
        return TPR461


class TestTPR502(TPRBase):
    @pytest.fixture()
    def filename(self):
        return TPR502


class TestTPR504(TPRBase):
    @pytest.fixture()
    def filename(self):
        return TPR504


class TestTPR505(TPRBase):
    @pytest.fixture()
    def filename(self):
        return TPR505


class TestTPR510(TPRBase):
    @pytest.fixture()
    def filename(self):
        return TPR510


class TestTPR2016(TPRBase):
    @pytest.fixture()
    def filename(self):
        return TPR2016


def _test_is_in_topology(name, elements, topology_path, topology_section):
    """
    Test if an interaction appears as expected in the topology
    """
    universe = MDAnalysis.Universe(topology_path)
    parser = MDAnalysis.topology.TPRParser.TPRParser(topology_path)
    top = parser.parse()
    for element in elements:
        assert element in getattr(top, topology_section).values, \
            'Interaction type "{}" not found'.format(name)


@pytest.mark.parametrize('topology', (
        TPR510_bonded,
        TPR2016_bonded
))
@pytest.mark.parametrize('bond', (
        ('BONDS', [(0, 1)]),
        ('G96BONDS', [(1, 2)]),
        ('MORSE', [(2, 3)]),
        ('CUBICBONDS', [(3, 4)]),
        ('CONNBONDS', [(4, 5)]),
        ('HARMONIC', [(5, 6)]),
        ('FENEBONDS', [(6, 7)]),
        ('RESTRAINTPOT', [(7, 8)]),
        ('TABBONDS', [(8, 9)]),
        ('TABBONDSNC', [(9, 10)]),
        ('CONSTR', [(10, 11)]),
        ('CONSTRNC', [(11, 12)]),
))
def test_all_bonds(topology, bond):
    """Test that all bond types are parsed as expected"""
    bond_type_in_topology = functools.partial(_test_is_in_topology,
                                              topology_section='bonds')
    bond_type, elements = bond
    bond_type_in_topology(bond_type, elements, topology)


@pytest.mark.parametrize('topology', (
    TPR510_bonded,
    TPR2016_bonded
))
@pytest.mark.parametrize('angle', (
    ('ANGLES',[(0, 1, 2)]),
    ('G96ANGLES',[(1, 2, 3)]),
    ('CROSS_BOND_BOND',[(2, 3, 4)]),
    ('CROSS_BOND_ANGLE',[(3, 4, 5)]),
    ('UREY_BRADLEY',[(4, 5, 6)]),
    ('QANGLES',[(5, 6, 7)]),
    ('RESTRANGLES',[(6, 7, 8)]),
    ('TABANGLES',[(7, 8, 9)]),
))
def test_all_angles(topology, angle):
    angle_type_in_topology = functools.partial(_test_is_in_topology,
                                               topology_section='angles')
    angle_type, elements = angle
    angle_type_in_topology(angle_type, elements, topology)


@pytest.mark.parametrize('topology', (
    TPR510_bonded,
    TPR2016_bonded
))
@pytest.mark.parametrize('dih', (
        ('PDIHS',[(0, 1, 2, 3), (1, 2, 3, 4), (7, 8, 9, 10)]),
        ('RBDIHS',[(4, 5, 6, 7)]),
        ('RESTRDIHS',[(8, 9, 10, 11)]),
        ('CBTDIHS',[(9, 10, 11, 12)]),
        ('FOURDIHS',[(6, 7, 8, 9)]),
        ('TABDIHS',[(10, 11, 12, 13)]),
))
def test_all_dihedrals(topology, dih):
    dih_type_in_topology = functools.partial(_test_is_in_topology,
                                             topology_section='dihedrals')
    dih_type, elements = dih
    dih_type_in_topology(dih_type, elements, topology)


@pytest.mark.parametrize('topology', (
    TPR510_bonded,
    TPR2016_bonded
))
@pytest.mark.parametrize('impr', (
    ('IDIHS', [(2, 3, 4, 5), (3, 4, 5, 6)]),
    ('PIDIHS', [(5, 6, 7, 8)])
))
def test_all_impropers(topology, impr):
    impr_type_in_topology = functools.partial(_test_is_in_topology,
                                              topology_section='impropers')

    impr_type, elements = impr
    impr_type_in_topology(impr_type, elements, topology)
