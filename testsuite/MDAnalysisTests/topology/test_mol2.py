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

import pytest
from numpy.testing import assert_equal

import MDAnalysis as mda

from MDAnalysisTests.topology.base import ParserBase
from MDAnalysisTests.datafiles import (
    mol2_molecule,
    mol2_molecules,
)
from numpy.testing import assert_equal


class TestMOL2Base(ParserBase):
    parser = mda.topology.MOL2Parser.MOL2Parser
    expected_attrs = [
        'ids', 'names', 'types', 'charges', 'resids', 'resnames', 'bonds'
    ]
    guessed_attrs = ['masses']
    expected_n_atoms = 49
    expected_n_residues = 1
    expected_n_segments = 1

    def test_attr_size(self, top):
        assert len(top.ids) == top.n_atoms
        assert len(top.names) == top.n_atoms
        assert len(top.types) == top.n_atoms
        assert len(top.charges) == top.n_atoms
        assert len(top.resids) == top.n_residues
        assert len(top.resnames) == top.n_residues

    def test_bonds(self, top):
        assert len(top.bonds) == 49  # bonds for 49 atoms
        assert len(top.bonds.values) == 51  # this many bonds

    @pytest.fixture(params=[mol2_molecule, mol2_molecules])
    def filename(self, request):
        return request.param


def test_bond_orders():
    ref_orders = ('am 1 1 2 1 2 1 1 am 1 1 am 2 2 '
                  '1 1 1 1 1 1 1 1 1 1 1 1 1 1 '
                  'ar ar ar 1 ar 1 ar 1 ar 1 1 1 '
                  '2 1 1 1 1 2 1 1 2 1 1').split()
    u = mda.Universe(mol2_molecule)
    orders = [bond.order for bond in u.atoms.bonds]
    assert_equal(orders, ref_orders)
