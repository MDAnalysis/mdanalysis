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
import pytest
from numpy.testing import assert_equal

import MDAnalysis as mda

from MDAnalysisTests.topology.base import ParserBase
from MDAnalysisTests.datafiles import (
    mol2_molecule,
    mol2_molecules,
)

import numpy as np
from io import StringIO

class TestMOL2Base(ParserBase):
    parser = mda.topology.MOL2Parser.MOL2Parser
    expected_attrs = [
        'ids', 'names', 'types', 'charges', 'resids', 'resnames', 'bonds',
        'elements',
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
        assert len(top.elements) == top.n_atoms

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


def test_elements():
    u = mda.Universe(mol2_molecule)

    assert_equal(
        u.atoms.elements[:5],
        np.array(["N", "S", "N", "N", "O"], dtype="U3")
    )


"""
# Test for #2927
def test_elements_selection():
    u = mda.Universe(mol2_molecule)
    ag = u.select_atoms("element S")

    assert_equal(
        ag.elements,
        np.array(["S", "S"], dtype="U3")
    )
"""

# Bond information is needed
# See #3057
mol2_wrong_element = """\
@<TRIPOS>MOLECULE
FXA101_1
49 51 1 0 0
SMALL
USER_CHARGES


@<TRIPOS>ATOM
  1 N1       6.8420     9.9900    22.7430 N.am  1 Q101  -0.8960
  2 S1       8.1400     9.2310    23.3330 X.o2  1 Q101   1.3220
  3 N2       4.4000     9.1300    20.4710 X.am  1 Q101  -0.3970
@<TRIPOS>BOND
  1   1   2  am
"""


def test_wrong_elements_warnings():
    with pytest.warns(UserWarning, match='Unknown element X found') as record:
        u = mda.Universe(StringIO(mol2_wrong_element), format='MOL2')

    # One warning from invalid elements, one from invalid masses
    assert len(record) == 2

    expected = np.array(['N', '', ''], dtype=object)
    assert_equal(u.atoms.elements, expected)
