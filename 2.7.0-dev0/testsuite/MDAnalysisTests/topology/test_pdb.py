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

import pytest
import numpy as np
from numpy.testing import assert_equal
import MDAnalysis as mda

from MDAnalysisTests.topology.base import ParserBase
from MDAnalysisTests.datafiles import (
    PDB,
    PDB_HOLE,
    PDB_small,
    PDB_conect,
    PDB_conect2TER,
    PDB_singleconect,
    PDB_chainidnewres,
    PDB_sameresid_diffresname,
    PDB_helix,
    PDB_elements,
    PDB_charges,
)
from MDAnalysis.topology.PDBParser import PDBParser
from MDAnalysis import NoDataError

_PDBPARSER = mda.topology.PDBParser.PDBParser

hybrid36 = [
    ("A0000", 100000),
    ("MEGAN", 20929695),
    ("J0NNY", 15247214),
    ("DREW6", 6417862),
    ("ST3V3", 31691119),
    ("ADA8M", 719798),
    ("a0000", 43770016),
    ("megan", 64599711),
    ("j0nny", 58917230),
    ("drew6", 50087878),
    ("st3v3", 75361135),
    ("ada8m", 44389814),
    ("    6", 6),
    ("   24", 24),
    ("  645", 645),
    (" 4951", 4951),
    ("10267", 10267)
]


@pytest.mark.parametrize('hybrid, integer', hybrid36)
def test_hy36decode(hybrid, integer):
    assert mda.topology.PDBParser.hy36decode(5, hybrid) == integer


class PDBBase(ParserBase):
    expected_attrs = ['ids', 'names', 'record_types', 'resids',
                      'resnames', 'altLocs', 'icodes', 'occupancies',
                      'tempfactors', 'chainIDs']
    guessed_attrs = ['types', 'masses']


class TestPDBParser(PDBBase):
    """This one has neither chainids or segids"""
    parser = mda.topology.PDBParser.PDBParser
    ref_filename = PDB
    expected_n_atoms = 47681
    expected_n_residues = 11302
    expected_n_segments = 1


class TestPDBParserSegids(PDBBase):
    """Has segids"""
    parser = mda.topology.PDBParser.PDBParser
    ref_filename = PDB_small
    expected_n_atoms = 3341
    expected_n_residues = 214
    expected_n_segments = 1


class TestPDBConect(object):
    """Testing PDB topology parsing (PDB)"""

    def test_conect_parser(self):
        lines = ("CONECT1233212331",
                 "CONECT123331233112334",
                 "CONECT123341233312335",
                 "CONECT123351233412336",
                 "CONECT12336123271233012335",
                 "CONECT12337 7718 84081234012344",
                 "CONECT1233812339123401234112345")
        results = ((12332, [12331]),
                   (12333, [12331, 12334]),
                   (12334, [12333, 12335]),
                   (12335, [12334, 12336]),
                   (12336, [12327, 12330, 12335]),
                   (12337, [7718, 8408, 12340, 12344]),
                   (12338, [12339, 12340, 12341, 12345]))
        for line, res in zip(lines, results):
            bonds = mda.topology.PDBParser._parse_conect(line)
            assert_equal(bonds[0], res[0])
            for bond, res_bond in zip(bonds[1], res[1]):
                assert_equal(bond, res_bond)

    def test_conect_parser_runtime(self):
        with pytest.raises(RuntimeError):
            mda.topology.PDBParser._parse_conect('CONECT12337 7718 '
                                                 '84081234012344123')

    def test_conect_topo_parser(self):
        """Check that the parser works as intended,
        and that the returned value is a dictionary
        """
        with _PDBPARSER(PDB_conect) as p:
            top = p.parse()
            assert isinstance(top, mda.core.topology.Topology)


def test_conect2ter():
    def parse():
        with PDBParser(PDB_conect2TER) as p:
            struc = p.parse()
        return struc

    with pytest.warns(UserWarning):
        struc = parse()

    assert hasattr(struc, 'bonds')
    assert len(struc.bonds.values) == 4


def test_single_conect():
    def parse():
        with PDBParser(PDB_singleconect) as p:
            struc = p.parse()
        return struc

    with pytest.warns(UserWarning):
        struc = parse()
    assert hasattr(struc, 'bonds')
    assert len(struc.bonds.values) == 2


def test_new_chainid_new_res():
    # parser must start new residue when chainid starts

    u = mda.Universe(PDB_chainidnewres)

    assert len(u.residues) == 4
    assert_equal(u.residues.resids, [1, 2, 3, 3])
    assert len(u.segments) == 4
    assert_equal(u.segments.segids, ['A', 'B', 'C', 'D'])
    assert len(u.segments[0].atoms) == 5
    assert len(u.segments[1].atoms) == 5
    assert len(u.segments[2].atoms) == 5
    assert len(u.segments[3].atoms) == 7


def test_sameresid_diffresname():
    with _PDBPARSER(PDB_sameresid_diffresname) as p:
        top = p.parse()
    resids = [9, 9]
    resnames = ['GLN', 'POPC']
    for i, (resid, resname) in enumerate(zip(resids, resnames)):
        assert top.resids.values[i] == resid
        assert top.resnames.values[i] == resname


def test_PDB_record_types():
    u = mda.Universe(PDB_HOLE)

    assert u.atoms[0].record_type == 'ATOM'
    assert u.atoms[132].record_type == 'HETATM'

    assert_equal(u.atoms[10:20].record_types, 'ATOM')
    assert_equal(u.atoms[271:].record_types, 'HETATM')


PDB_noresid = """\
REMARK For testing reading of CRYST
REMARK This has MODELs then CRYST entries
CRYST1   80.000   80.017   80.017  90.00  90.00  90.00 P 1           1
MODEL        1
ATOM      1  H2  TIP3           10.000  44.891  14.267  1.00  0.00      TIP3
ATOM      2  OH2 TIP3           67.275  48.893  23.568  1.00  0.00      TIP3
ATOM      3  H1  TIP3           66.641  48.181  23.485  1.00  0.00      TIP3
ATOM      4  H2  TIP3           66.986  49.547  22.931  1.00  0.00      TIP3
ENDMDL
"""


def test_PDB_no_resid():
    u = mda.Universe(StringIO(PDB_noresid), format='PDB')

    assert len(u.atoms) == 4
    assert len(u.residues) == 1
    # should have default resid of 1
    assert u.residues[0].resid == 1


PDB_hex = """\
REMARK For testing reading of hex atom numbers
REMARK This has MODELs then hex atom numbers entries
CRYST1   80.000   80.017   80.017  90.00  90.00  90.00 P 1           1
MODEL        1
HETATM    1  H     2 L 400      20.168  00.034  40.428
HETATMA0000  H     2 L 400      40.168  50.034  40.428
HETATMA0001  H     2 L 400      30.453  60.495  50.132
HETATMA0002  H     2 L 400      20.576  40.354  60.483
HETATMA0003  H     2 L 400      10.208  30.067  70.045
END
"""

# this causes an error on Win64/Python 3.8 on Azure when loaded
# in as a file instead
PDB_metals = """\
HETATM    1 CU    CU A   1      00.000  00.000  00.000  1.00 00.00          Cu
HETATM    2 FE    FE A   2      03.000  03.000  03.000  1.00 00.00          Fe
HETATM    3 Ca    Ca A   3      03.000  03.000  03.000  1.00 00.00          Ca
HETATM    3 Mg    Mg A   3      03.000  03.000  03.000  1.00 00.00          Mg
"""


def test_PDB_hex():
    u = mda.Universe(StringIO(PDB_hex), format='PDB')
    assert len(u.atoms) == 5
    assert u.atoms[0].id == 1
    assert u.atoms[1].id == 100000
    assert u.atoms[2].id == 100001
    assert u.atoms[3].id == 100002
    assert u.atoms[4].id == 100003


@pytest.mark.filterwarnings("error:Failed to guess the mass")
def test_PDB_metals():
    from MDAnalysis.topology import tables

    u = mda.Universe(StringIO(PDB_metals), format='PDB')

    assert len(u.atoms) == 4
    assert u.atoms[0].mass == pytest.approx(tables.masses["CU"])
    assert u.atoms[1].mass == pytest.approx(tables.masses["FE"])
    assert u.atoms[2].mass == pytest.approx(tables.masses["CA"])
    assert u.atoms[3].mass == pytest.approx(tables.masses["MG"])


def test_PDB_elements():
    """The test checks whether elements attribute are assigned
    properly given a PDB file with valid elements record.
    """
    u = mda.Universe(PDB_elements, format='PDB')
    element_list = np.array(['N', 'C', 'C', 'O', 'C', 'C', 'O', 'N', 'H',
                             'H', 'H', 'H', 'H', 'H', 'H', 'H', 'Cu', 'Fe',
                             'Mg', 'Ca', 'S', 'O', 'C', 'C', 'S', 'O', 'C',
                             'C'], dtype=object)
    assert_equal(u.atoms.elements, element_list)


def test_missing_elements_noattribute():
    """Check that:

    1) a warning is raised if elements are missing
    2) the elements attribute is not set
    """
    wmsg = ("Element information is missing, elements attribute will not be "
            "populated")
    with pytest.warns(UserWarning, match=wmsg):
        u = mda.Universe(PDB_small)
    with pytest.raises(AttributeError):
        _ = u.atoms.elements


PDB_wrong_ele = """\
REMARK For testing warnings of wrong elements
REMARK This file represent invalid elements in the elements column
ATOM      1  N   ASN A   1      -8.901   4.127  -0.555  1.00  0.00           N
ATOM      2  CA  ASN A   1      -8.608   3.135  -1.618  1.00  0.00
ATOM      3  C   ASN A   1      -7.117   2.964  -1.897  1.00  0.00           C
ATOM      4  O   ASN A   1      -6.634   1.849  -1.758  1.00  0.00           O
ATOM      5  X   ASN A   1      -9.437   3.396  -2.889  1.00  0.00          XX
TER       6
HETATM    7 CU    CU A   2      03.000  00.000  00.000  1.00 00.00          CU
HETATM    8 FE    FE A   3      00.000  03.000  00.000  1.00 00.00          Fe
HETATM    9 Mg    Mg A   4      03.000  03.000  03.000  1.00 00.00          MG
TER       10
"""


def test_wrong_elements_warnings():
    """The test checks whether there are invalid elements in the elements
    column which have been parsed and returns an appropriate warning.
    """
    with pytest.warns(UserWarning, match='Unknown element XX found'):
        u = mda.Universe(StringIO(PDB_wrong_ele), format='PDB')

    expected = np.array(['N', '', 'C', 'O', '', 'Cu', 'Fe', 'Mg'],
                        dtype=object)
    assert_equal(u.atoms.elements, expected)


def test_nobonds_error():
    """Issue #2832: PDB without CONECT record should not have a bonds
    attribute and raises NoDataError on access"""
    u = mda.Universe(PDB_helix)

    errmsg = "This Universe does not contain bonds information"

    with pytest.raises(NoDataError, match=errmsg):
        u.atoms.bonds


def test_PDB_charges():
    """The test checks whether formalcharges attribute are assigned
    properly given a PDB file with a valid formal charges record.
    """
    u = mda.Universe(PDB_charges)
    formal_charges = np.array([0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0], dtype=int)
    assert_equal(u.atoms.formalcharges, formal_charges)


PDB_charges_nosign = """\
REMARK Invalid charge assignment - no sign for MG2+
HETATM    1 CU    CU A   1      03.000  00.000  00.000  1.00 00.00          CU2+
HETATM    2 FE    FE A   2      00.000  03.000  00.000  1.00 00.00          Fe2+
HETATM    3 Mg    Mg A   3      03.000  03.000  03.000  1.00 00.00          MG2
END
"""


PDB_charges_invertsign = """\
REMARK Invalid charge format for MG2+
HETATM    1 CU    CU A   1      03.000  00.000  00.000  1.00 00.00          CU2+
HETATM    2 FE    FE A   2      00.000  03.000  00.000  1.00 00.00          Fe2+
HETATM    3 Mg    Mg A   3      03.000  03.000  03.000  1.00 00.00          MG+2
END
"""


@pytest.mark.parametrize('infile,entry', [
        [PDB_charges_nosign, r'2'],
        [PDB_charges_invertsign, r'\+2']
])
def test_PDB_bad_charges(infile, entry):
    """
    Test that checks that a warning is raised and formal charges are not set:
        * If there are missing signs for a given formal charge entry
        * If there is an unrecognised formal charge entry
    """
    wmsg = f"Unknown entry {entry} encountered in formal charge field."
    with pytest.warns(UserWarning, match=wmsg):
        u = mda.Universe(StringIO(infile), format='PDB')
        assert not hasattr(u, 'formalcharges')
