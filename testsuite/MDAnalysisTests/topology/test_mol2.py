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
from io import StringIO

import numpy as np
from numpy.testing import assert_equal, assert_allclose
import pytest

import MDAnalysis as mda
from MDAnalysisTests.topology.base import ParserBase

from MDAnalysisTests.datafiles import (
    mol2_molecule,
    mol2_molecules,
)


mol2_wo_opt_col = """\
@<TRIPOS>MOLECULE
mol2_wo_opt_col
2
SMALL
NO_CHARGES
@<TRIPOS>ATOM
  1 N1       6.8420     9.9900    22.7430 N.am
  2 N2       4.4000     9.1300    20.4710 N.am
"""

mol2_partial_opt_col = """\
@<TRIPOS>MOLECULE
mol2_partial_opt_col
2
SMALL
NO_CHARGES
@<TRIPOS>ATOM
  1 N1       6.8420     9.9900    22.7430 N.am  1
  2 N2       4.4000     9.1300    20.4710 N.am  2
"""

mol2_resname_unformat = """\
@<TRIPOS>MOLECULE
mol2_resname_unformat
2 1 0 0 0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
  1 O        3.0000    0.0000    0.0000 O.3     1
  2 C        2.0000    0.0000    0.0000 C.2     1  UNL1
"""

mol2_wo_required_col = """\
@<TRIPOS>MOLECULE
mol2_wo_required_col
2
SMALL
NO_CHARGES
@<TRIPOS>ATOM
  1 N1       6.8420     9.9900    22.7430
  2 N2       4.4000     9.1300    20.4710
"""

mol2_no_charge_error1 = """\
@<TRIPOS>MOLECULE
mol2_no_charge_error1
2 1 0 0 0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
  1 O        3.0000    0.0000    0.0000 O.3     1  UNL1       -0.0520
  2 C        2.0000    0.0000    0.0000 C.2     1  UNL1        0.0520
"""

mol2_no_charge_error2 = """\
@<TRIPOS>MOLECULE
mol2_no_charge_error2
2 1 0 0 0
SMALL
MMFF94_CHARGES

@<TRIPOS>ATOM
  1 O        3.0000     0.0000    0.0000 O.3    1  UNL1
  2 C        2.0000     0.0000    0.0000 C.2    1  UNL1
"""

mol2_wrong_element = """\
@<TRIPOS>MOLECULE
mol2_wrong_element
3 0 0 0 0
SMALL
USER_CHARGES


@<TRIPOS>ATOM
  1 N1       6.8420     9.9900    22.7430 N.am  1 Q101  -0.8960
  2 S1       8.1400     9.2310    23.3330 X.o2  1 Q101   1.3220
  3 N2       4.4000     9.1300    20.4710 XX.am  1 Q101  -0.3970
"""

mol2_all_wrong_elements = """\
@<TRIPOS>MOLECULE
mol2_all_wrong_elements
3 0 0 0 0
SMALL
USER_CHARGES


@<TRIPOS>ATOM
  1 N1       6.8420     9.9900    22.7430 X.am  1 Q101  -0.8960
  2 S1       8.1400     9.2310    23.3330 X.o2  1 Q101   1.3220
  3 N2       4.4000     9.1300    20.4710 XX.am  1 Q101  -0.3970
"""


mol2_fake = """\
@<TRIPOS>MOLECULE
mol2_fake
26 0 0 0 0
SMALL
USER_CHARGES


@<TRIPOS>ATOM
  1 H1       0.0000     1.0000    10.0000 H.spc  1 XXXX   5.0000
  2 H2       0.0000     1.0000    10.0000 H.t3p  1 XXXX   5.0000
  3 H3       0.0000     1.0000    10.0000 H.xyz  1 XXXX   5.0000
  4 C1       0.0000     1.0000    10.0000 C.1  1 XXXX   5.0000
  5 C2       0.0000     1.0000    10.0000 C.2  1 XXXX   5.0000
  5 C3       0.0000     1.0000    10.0000 C.3  1 XXXX   5.0000
  6 C4       0.0000     1.0000    10.0000 C.ar  1 XXXX   5.0000
  7 C5       0.0000     1.0000    10.0000 C.cat  1 XXXX   5.0000
  8 C6       0.0000     1.0000    10.0000 C.xyz  1 XXXX   5.0000
  9 N1       0.0000     1.0000    10.0000 N.1  1 XXXX   5.0000
  10 N2       0.0000     1.0000    10.0000 N.2  1 XXXX   5.0000
  11 N3       0.0000     1.0000    10.0000 N.3  1 XXXX   5.0000
  12 N4       0.0000     1.0000    10.0000 N.ar  1 XXXX   5.0000
  13 O1       0.0000     1.0000    10.0000 O.2  1 XXXX   5.0000
  14 O2       0.0000     1.0000    10.0000 O.3  1 XXXX   5.0000
  15 O3       0.0000     1.0000    10.0000 O.co2  1 XXXX   5.0000
  16 O4       0.0000     1.0000    10.0000 O.spc  1 XXXX   5.0000
  16 O5       0.0000     1.0000    10.0000 O.t3p  1 XXXX   5.0000
  17 S1       0.0000     1.0000    10.0000 S.3  1 XXXX   5.0000
  18 S2       0.0000     1.0000    10.0000 S.2  1 XXXX   5.0000
  19 S3       0.0000     1.0000    10.0000 S.O  1 XXXX   5.0000
  20 S4       0.0000     1.0000    10.0000 S.O2  1 XXXX   5.0000
  21 S5       0.0000     1.0000    10.0000 S.o  1 XXXX   5.0000
  22 S6       0.0000     1.0000    10.0000 S.o2  1 XXXX   5.0000
  23 P1       0.0000     1.0000    10.0000 P.3  1 XXXX   5.0000
  24 Cr1       0.0000     1.0000    10.0000 Cr.th  1 XXXX   5.0000
  25 Cr2       0.0000     1.0000    10.0000 Cr.oh  1 XXXX   5.0000
  26 Co1       0.0000     1.0000    10.0000 Co.oh  1 XXXX   5.0000
"""


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


# Test for #2927
def test_elements_selection():
    u = mda.Universe(mol2_molecule)
    ag = u.select_atoms("element S")

    assert_equal(
        ag.elements,
        np.array(["S", "S"], dtype="U3")
    )


def test_wrong_elements_warnings():
    with pytest.warns(UserWarning, match='Unknown elements found') as record:
        u = mda.Universe(StringIO(mol2_wrong_element), format='MOL2')

    # One warning from invalid elements, one from masses PendingDeprecationWarning
    assert len(record) == 3

    expected_elements = np.array(['N', '', ''], dtype=object)
    guseed_masses = np.array([14.007, 0.0, 0.0], dtype=float)
    gussed_types = np.array(['N.am', 'X.o2', 'XX.am'])

    assert_equal(u.atoms.elements, expected_elements)
    assert_equal(u.atoms.types, gussed_types)
    assert_allclose(u.atoms.masses, guseed_masses)


def test_all_wrong_elements_warnings():
    with pytest.warns(UserWarning, match='Unknown elements found'):
        u = mda.Universe(StringIO(mol2_all_wrong_elements), format='MOL2')

    with pytest.raises(mda.exceptions.NoDataError,
                       match='This Universe does not contain element '
                       'information'):

        u.atoms.elements


def test_all_elements():
    with pytest.warns(UserWarning, match='Unknown elements found'):
        u = mda.Universe(StringIO(mol2_fake), format='MOL2')

    expected = ["H"] * 2 + [""] + ["C"] * 5 + [""] + ["N"] * 4 + ["O"] * 5 + \
        ["S"] * 6 + ["P"] + ["Cr"] * 2 + ["Co"]
    expected = np.array(expected, dtype=object)
    assert_equal(u.atoms.elements, expected)


# Test for Issue #3385 / PR #3598
def test_wo_optional_columns():
    u = mda.Universe(StringIO(mol2_wo_opt_col), format='MOL2')
    assert_equal(u.atoms.resids, np.array([1, 1]))
    with pytest.raises(mda.exceptions.NoDataError):
        u.atoms.resnames
    with pytest.raises(mda.exceptions.NoDataError):
        u.atoms.charges


def test_partial_optional_columns():
    u = mda.Universe(StringIO(mol2_partial_opt_col), format='MOL2')
    assert_equal(u.atoms.resids, np.array([1, 2]))
    with pytest.raises(mda.exceptions.NoDataError):
        u.atoms.resnames
    with pytest.raises(mda.exceptions.NoDataError):
        u.atoms.charges


def test_mol2_wo_required_columns():
    with pytest.raises(ValueError,
                       match='The @<TRIPOS>ATOM block in mol2 file'):
        u = mda.Universe(StringIO(mol2_wo_required_col), format='MOL2')


def test_mol2_no_charges():
    with pytest.raises(ValueError,
                       match='indicates no charges'):
        u = mda.Universe(StringIO(mol2_no_charge_error1), format='MOL2')
    with pytest.raises(ValueError,
                       match='indicates a charge model'):
        u = mda.Universe(StringIO(mol2_no_charge_error2), format='MOL2')


def test_unformat():
    with pytest.raises(ValueError,
                       match='Some atoms in the mol2 file'):
        u = mda.Universe(StringIO(mol2_resname_unformat), format='MOL2')


def test_guessed_masses():
    u = mda.Universe(mol2_molecules)
    assert_allclose(u.atoms.masses[:7], [14.007, 32.06,
                    14.007, 14.007, 15.999, 15.999, 12.011])
