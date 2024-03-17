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
import numpy as np
from numpy.testing import assert_equal
from MDAnalysisTests.topology.base import ParserBase
from MDAnalysisTests.datafiles import (
    PRM,  # ache.prmtop
    PRM_chainid_bz2,  # multi_anche.prmtop.bz2
    PRM12,  # anti.top
    PRM7,  # tz2.truncoct.parm7.bz2
    PRMpbc,
    PRMNCRST,
    PRMNEGATIVE,
    PRMErr1,
    PRMErr2,
    PRMErr3,
    PRMErr4,
    PRMErr5,
    PRM_UreyBradley,
    PRM19SBOPC,
)


class TOPBase(ParserBase):
    parser = mda.topology.TOPParser.TOPParser
    expected_attrs = [
        "names", "types", "type_indices", "charges", "masses", "resnames",
        "bonds", "angles", "dihedrals", "impropers", "elements"
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

    def test_bonds_atom_counts(self, filename):
        u = mda.Universe(filename)
        assert len(u.atoms[[0]].bonds) == self.expected_n_zero_bonds
        assert len(u.atoms[[self.atom_i]].bonds) == self.expected_n_i_bonds

    def test_angles_atom_counts(self, filename):
        u = mda.Universe(filename)
        assert len(u.atoms[[0]].angles) == self.expected_n_zero_angles
        assert len(u.atoms[[self.atom_i]].angles) == self.expected_n_i_angles

    def test_dihedrals_atom_counts(self, filename):
        u = mda.Universe(filename)
        assert len(u.atoms[[0]].dihedrals) == self.expected_n_zero_dihedrals
        assert len(u.atoms[[self.atom_i]].dihedrals) == \
            self.expected_n_i_dihedrals

    def test_impropers_atom_counts(self, filename):
        u = mda.Universe(filename)
        assert len(u.atoms[[0]].impropers) == self.expected_n_zero_impropers
        assert len(u.atoms[[self.atom_i]].impropers) == \
            self.expected_n_i_impropers

    def test_bonds_identity(self, top):
        vals = top.bonds.values
        for bond in self.atom_zero_bond_values:
            assert (bond in vals) or (bond[::-1] in vals)
        for bond in self.atom_i_bond_values:
            assert (bond in vals) or (bond[::-1] in vals)

    def test_angles_identity(self, top):
        vals = top.angles.values
        for ang in self.atom_zero_angle_values:
            assert (ang in vals) or (ang[::-1] in vals)
        for ang in self.atom_i_angle_values:
            assert (ang in vals) or (ang[::-1] in vals)

    def test_dihedrals_identity(self, top):
        vals = top.dihedrals.values
        for dih in self.atom_zero_dihedral_values:
            assert (dih in vals) or (dih[::-1] in vals)
        for dih in self.atom_i_dihedral_values:
            assert (dih in vals) or (dih[::-1] in vals)

    def test_impropers_identity(self, top):
        vals = top.impropers.values
        for imp in self.atom_zero_improper_values:
            assert (imp in vals) or (imp[::-1] in vals)
        for imp in self.atom_i_improper_values:
            assert (imp in vals) or (imp[::-1] in vals)

    def test_angle_atoms_bonded(self, top):
        vals = top.bonds.values
        for ang in top.angles.values:
            for b in ((ang[0], ang[1]), (ang[1], ang[2])):
                assert (b in vals) or (b[::-1] in vals)

    def test_dihedral_atoms_bonded(self, top):
        vals = top.bonds.values
        for dih in top.dihedrals.values:
            for b in ((dih[0], dih[1]), (dih[1], dih[2]), (dih[2], dih[3])):
                assert (b in vals) or (b[::-1] in vals)

    def test_improper_atoms_bonded(self, top):
        vals = top.bonds.values
        for imp in top.impropers.values:
            forward = ((imp[0], imp[2]), (imp[1], imp[2]), (imp[2], imp[3]))
            backward = ((imp[0], imp[1]), (imp[1], imp[2]), (imp[1], imp[3]))
            for a, b in zip(forward, backward):
                assert ((b in vals) or (b[::-1] in vals) or
                        (a in vals) or (a[::-1] in vals))

    def test_elements(self, top):
        """Tests elements attribute.

        If elements present, loops over ranges of the topology elements list
        and compared against a provided list of expected values.
        Otherwise, checks that elements are not in the topology attributes.
        """

        if self.expected_elems:
            for erange, evals in zip(self.elems_ranges, self.expected_elems):
                assert_equal(top.elements.values[erange[0]:erange[1]], evals,
                             "unexpected element match")
        else:
            assert not hasattr(top, 'elements'), 'Unexpected elements attr'


class TestPRMParser(TOPBase):
    ref_filename = PRM
    # Does not contain an ATOMIC_NUMBER record, so no elements
    expected_attrs = [
        "names", "types", "type_indices", "charges", "masses", "resnames",
        "bonds", "angles", "dihedrals", "impropers"
    ]
    expected_n_atoms = 252
    expected_n_residues = 14
    expected_n_bonds = 259
    expected_n_angles = 456
    expected_n_dihedrals = 673
    expected_n_impropers = 66
    atom_i = 79
    expected_n_zero_bonds = 4
    expected_n_i_bonds = 3
    expected_n_zero_angles = 9
    expected_n_i_angles = 9
    expected_n_zero_dihedrals = 14
    expected_n_i_dihedrals = 15
    expected_n_zero_impropers = 0
    expected_n_i_impropers = 4
    atom_zero_bond_values = ((0, 4), (0, 1), (0, 2), (0, 3))
    atom_i_bond_values = ((79, 80), (79, 83), (77, 79))
    atom_zero_angle_values = ((0, 4, 6), (0, 4, 10), (3, 0, 4),
                              (2, 0, 3), (2, 0, 4), (1, 0, 2),
                              (1, 0, 3), (1, 0, 4), (0, 4, 5))
    atom_i_angle_values = ((80, 79, 83), (77, 79, 80), (77, 79, 83),
                           (74, 77, 79), (79, 80, 81), (79, 80, 82),
                           (79, 83, 84), (79, 83, 85), (78, 77, 79))
    atom_zero_dihedral_values = ((0, 4, 10, 11), (0, 4, 10, 12),
                                 (3, 0, 4, 5), (3, 0, 4, 6),
                                 (3, 0, 4, 10), (2, 0, 4, 5),
                                 (2, 0, 4, 6), (2, 0, 4, 10),
                                 (1, 0, 4, 5), (1, 0, 4, 6),
                                 (1, 0, 4, 10), (0, 4, 6, 7),
                                 (0, 4, 6, 8), (0, 4, 6, 9))
    atom_i_dihedral_values = ((71, 74, 77, 79), (74, 77, 79, 80),
                              (74, 77, 79, 83), (75, 74, 77, 79),
                              (76, 74, 77, 79), (77, 79, 80, 81),
                              (77, 79, 80, 82), (77, 79, 83, 84),
                              (77, 79, 83, 85), (78, 77, 79, 80),
                              (78, 77, 79, 83), (80, 79, 83, 84),
                              (80, 79, 83, 85), (81, 80, 79, 83),
                              (82, 80, 79, 83))
    atom_zero_improper_values = ()
    atom_i_improper_values = ((74, 79, 77, 78), (77, 80, 79, 83),
                              (79, 81, 80, 82), (79, 84, 83, 85))
    expected_elems = None


class TestPRMChainidParser(TOPBase):
    ref_filename = PRM_chainid_bz2
    # Checks the reading of %FLAG RESIDUE_CHAINID. See PR #4007
    expected_attrs = [
        "names",
        "types",
        "type_indices",
        "charges",
        "masses",
        "resnames",
        "bonds",
        "angles",
        "dihedrals",
        "impropers",
        "elements",
        "chainIDs",
    ]
    expected_n_atoms = 677
    expected_n_residues = 38
    expected_n_segments = 3
    expected_n_bonds = 695
    expected_n_angles = 1220
    expected_n_dihedrals = 1797
    expected_n_impropers = 189
    atom_i = 79
    expected_n_zero_bonds = 4
    expected_n_i_bonds = 3
    expected_n_zero_angles = 9
    expected_n_i_angles = 9
    expected_n_zero_dihedrals = 14
    expected_n_i_dihedrals = 15
    expected_n_zero_impropers = 0
    expected_n_i_impropers = 4
    atom_zero_bond_values = ((0, 4), (0, 1), (0, 2), (0, 3))
    atom_i_bond_values = ((79, 80), (79, 83), (77, 79))
    atom_zero_angle_values = (
        (0, 4, 6),
        (0, 4, 10),
        (3, 0, 4),
        (2, 0, 3),
        (2, 0, 4),
        (1, 0, 2),
        (1, 0, 3),
        (1, 0, 4),
        (0, 4, 5),
    )
    atom_i_angle_values = (
        (80, 79, 83),
        (77, 79, 80),
        (77, 79, 83),
        (74, 77, 79),
        (79, 80, 81),
        (79, 80, 82),
        (79, 83, 84),
        (79, 83, 85),
        (78, 77, 79),
    )
    atom_zero_dihedral_values = (
        (0, 4, 10, 11),
        (0, 4, 10, 12),
        (3, 0, 4, 5),
        (3, 0, 4, 6),
        (3, 0, 4, 10),
        (2, 0, 4, 5),
        (2, 0, 4, 6),
        (2, 0, 4, 10),
        (1, 0, 4, 5),
        (1, 0, 4, 6),
        (1, 0, 4, 10),
        (0, 4, 6, 7),
        (0, 4, 6, 8),
        (0, 4, 6, 9),
    )
    atom_i_dihedral_values = (
        (71, 74, 77, 79),
        (74, 77, 79, 80),
        (74, 77, 79, 83),
        (75, 74, 77, 79),
        (76, 74, 77, 79),
        (77, 79, 80, 81),
        (77, 79, 80, 82),
        (77, 79, 83, 84),
        (77, 79, 83, 85),
        (78, 77, 79, 80),
        (78, 77, 79, 83),
        (80, 79, 83, 84),
        (80, 79, 83, 85),
        (81, 80, 79, 83),
        (82, 80, 79, 83),
    )
    atom_zero_improper_values = ()
    atom_i_improper_values = (
        (74, 79, 77, 78),
        (77, 80, 79, 83),
        (79, 81, 80, 82),
        (79, 84, 83, 85),
    )
    elems_ranges = [[0, 9], [250, 257], [500, 508]]

    expected_elems = [
        np.array(
            [
                "N",
                "H",
                "H",
                "H",
                "C",
                "H",
                "C",
                "H",
                "H",
            ],
            dtype=object,
        ),
        np.array(
            [
                "O",
                "O",
                "N",
                "H",
                "H",
                "H",
                "C",
            ],
            dtype=object,
        ),
        np.array(["H", "C", "O", "O", "N", "H", "H", "H"], dtype=object),
    ]

    expected_chainIDs = np.array(
        [
            "A",
            "A",
            "A",
            "A",
            "A",
            "A",
            "A",
            "A",
            "A",
            "A",
            "A",
            "A",
            "A",
            "A",
            "B",
            "B",
            "B",
            "B",
            "B",
            "B",
            "B",
            "B",
            "B",
            "B",
            "B",
            "B",
            "B",
            "B",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
            "C",
        ]
    )

    def test_chainIDs(self, filename):
        """Tests chainIDs attribute.

        If RESIDUE_CHAINID present, residue chainIDs are compared against a
        provided list of expected values.
        Otherwise, checks that elements are not in the topology attributes.
        """

        u = mda.Universe(filename)
        if hasattr(self, "expected_chainIDs"):
            reschainIDs = [atomchainIDs[0] for atomchainIDs in u.residues.chainIDs]
            assert_equal(
                reschainIDs, self.expected_chainIDs, "unexpected element match"
            )
        else:
            assert not hasattr(u.atoms, "chainIDs"), "Unexpected chainIDs attr"


class TestPRM12Parser(TOPBase):
    ref_filename = PRM12
    expected_n_atoms = 8923
    expected_n_residues = 2861
    expected_n_bonds = 8947
    expected_n_angles = 756
    expected_n_dihedrals = 1128
    expected_n_impropers = 72
    expected_n_zero_bonds = 1
    expected_n_i_bonds = 4
    expected_n_zero_angles = 1
    expected_n_i_angles = 12
    expected_n_zero_dihedrals = 3
    expected_n_i_dihedrals = 28
    expected_n_zero_impropers = 0
    expected_n_i_impropers = 1
    atom_i = 335
    ref_proteinatoms = 0
    atom_zero_bond_values = ((0, 1),)
    atom_i_bond_values = ((335, 337), (335, 354),
                          (334, 335), (335, 336))
    atom_zero_angle_values = ((0, 1, 2),)
    atom_i_angle_values = ((337, 335, 354), (335, 337, 338),
                           (335, 337, 351), (335, 354, 352),
                           (334, 335, 337), (334, 335, 354),
                           (332, 334, 335), (336, 335, 337),
                           (336, 335, 354), (335, 354, 355),
                           (335, 354, 356), (334, 335, 336))
    atom_zero_dihedral_values = ((0, 1, 2, 3), (0, 1, 2, 4),
                                 (0, 1, 2, 5))
    atom_i_dihedral_values = ((329, 332, 334, 335), (332, 334, 335, 336),
                              (332, 334, 335, 337), (332, 334, 335, 354),
                              (332, 352, 354, 335), (333, 332, 334, 335),
                              (334, 335, 337, 338), (334, 335, 337, 351),
                              (334, 335, 354, 352), (334, 335, 354, 355),
                              (334, 335, 354, 356), (335, 334, 332, 352),
                              (335, 337, 338, 339), (335, 337, 338, 340),
                              (335, 337, 351, 341), (335, 337, 351, 350),
                              (335, 354, 352, 353), (335, 354, 352, 357),
                              (336, 335, 337, 338), (336, 335, 337, 351),
                              (336, 335, 354, 352), (336, 335, 354, 355),
                              (336, 335, 354, 356), (337, 335, 354, 352),
                              (337, 335, 354, 355), (337, 335, 354, 356),
                              (338, 337, 335, 354), (351, 337, 335, 354))
    atom_zero_improper_values = ()
    atom_i_improper_values = ((335, 337, 338, 351),)
    elems_ranges = [[0, 36], [351, 403]]
    expected_elems = [np.array(["H", "O", "C", "H", "H", "C", "H", "O", "C",
                                "H", "N", "C", "H", "N", "C", "C", "O", "N",
                                "H", "C", "N", "H", "H", "N", "C", "C", "H",
                                "C", "H", "H", "O", "P", "O", "O", "O", "C"],
                      dtype=object),
                      np.array(["C", "C", "H", "C", "H", "H", "O", "P", "O",
                                "O", "O", "C", "H", "H", "C", "H", "O", "C",
                                "H", "N", "C", "H", "N", "C", "C", "O", "N",
                                "H", "C", "N", "H", "H", "N", "C", "C", "H",
                                "C", "H", "H", "O", "H", "Na", "Na", "Na",
                                "Na", "Na", "Na", "Na", "Na", "O", "H", "H"],
                      dtype=object)]


class TestParm7Parser(TOPBase):
    ref_filename = PRM7
    # Does not contain an ATOMIC_NUMBER record, so no elements
    expected_attrs = [
        "names", "types", "type_indices", "charges", "masses", "resnames",
        "bonds", "angles", "dihedrals", "impropers"
    ]
    expected_n_atoms = 5827
    expected_n_residues = 1882
    expected_n_bonds = 5834
    expected_n_angles = 402
    expected_n_dihedrals = 602
    expected_n_impropers = 55
    atom_i = 135
    expected_n_zero_bonds = 4
    expected_n_i_bonds = 4
    expected_n_zero_angles = 9
    expected_n_i_angles = 13
    expected_n_zero_dihedrals = 14
    expected_n_i_dihedrals = 27
    expected_n_zero_impropers = 0
    expected_n_i_impropers = 2
    atom_zero_bond_values = ((0, 4), (0, 1), (0, 2), (0, 3))
    atom_i_bond_values = ((135, 137), (135, 155), (133, 135),
                          (135, 136))
    atom_zero_angle_values = ((0, 4, 6), (0, 4, 11), (3, 0, 4),
                              (2, 0, 3), (2, 0, 4), (1, 0, 2),
                              (1, 0, 3), (1, 0, 4), (0, 4, 5))
    atom_i_angle_values = ((131, 133, 135), (137, 135, 155),
                           (135, 137, 140), (135, 155, 156),
                           (135, 155, 157), (133, 135, 137),
                           (133, 135, 155), (136, 135, 137),
                           (136, 135, 155), (135, 137, 138),
                           (135, 137, 139), (134, 133, 135),
                           (133, 135, 136))
    atom_zero_dihedral_values = ((0, 4, 6, 7), (0, 4, 6, 8),
                                 (0, 4, 6, 9), (0, 4, 11, 12),
                                 (0, 4, 11, 13), (1, 0, 4, 5),
                                 (1, 0, 4, 6), (1, 0, 4, 11),
                                 (2, 0, 4, 5), (2, 0, 4, 6),
                                 (2, 0, 4, 11), (3, 0, 4, 5),
                                 (3, 0, 4, 6), (3, 0, 4, 11))
    atom_i_dihedral_values = ((113, 131, 133, 135), (131, 133, 135, 136),
                              (131, 133, 135, 137), (131, 133, 135, 155),
                              (132, 131, 133, 135), (133, 135, 137, 138),
                              (133, 135, 137, 139), (133, 135, 137, 140),
                              (133, 135, 155, 156), (133, 135, 155, 157),
                              (134, 133, 135, 136), (134, 133, 135, 137),
                              (134, 133, 135, 155), (135, 137, 140, 141),
                              (135, 137, 140, 154), (135, 155, 157, 158),
                              (135, 155, 157, 159), (136, 135, 137, 138),
                              (136, 135, 137, 139), (136, 135, 137, 140),
                              (136, 135, 155, 156), (136, 135, 155, 157),
                              (137, 135, 155, 156), (137, 135, 155, 157),
                              (138, 137, 135, 155), (139, 137, 135, 155),
                              (140, 137, 135, 155))
    atom_zero_improper_values = ()
    atom_i_improper_values = ((131, 135, 133, 134), (135, 157, 155, 156))
    expected_elems = None


class TestPRM2(TOPBase):
    ref_filename = PRMpbc
    # Does not contain an ATOMIC_NUMBER record, so no elements
    expected_attrs = [
        "names", "types", "type_indices", "charges", "masses", "resnames",
        "bonds", "angles", "dihedrals", "impropers"
    ]
    expected_n_atoms = 5071
    expected_n_residues = 1686
    ref_proteinatoms = 22
    expected_n_bonds = 5070
    expected_n_angles = 36
    expected_n_dihedrals = 41
    expected_n_impropers = 4
    atom_i = 14
    expected_n_zero_bonds = 1
    expected_n_i_bonds = 3
    expected_n_zero_angles = 3
    expected_n_i_angles = 8
    expected_n_zero_dihedrals = 2
    expected_n_i_dihedrals = 18
    expected_n_zero_impropers = 0
    expected_n_i_impropers = 2
    atom_zero_bond_values = ((0, 1),)
    atom_i_bond_values = ((14, 15), (14, 16), (8, 14))
    atom_zero_angle_values = ((0, 1, 2), (0, 1, 3), (0, 1, 4))
    atom_i_angle_values = ((15, 14, 16), (14, 16, 18), (10,  8, 14),
                           (8, 14, 15), (8, 14, 16), (6,  8, 14),
                           (14, 16, 17), (9,  8, 14))
    atom_zero_dihedral_values = ((0, 1, 4, 5), (0, 1, 4, 6))
    atom_i_dihedral_values = ((4, 6, 8, 14), (6, 8, 14, 15),
                              (6, 8, 14, 16), (7, 6, 8, 14),
                              (8, 14, 16, 17), (8, 14, 16, 18),
                              (9, 8, 14, 15), (9, 8, 14, 16),
                              (10, 8, 14, 15), (10, 8, 14, 16),
                              (11, 10, 8, 14), (12, 10, 8, 14),
                              (13, 10, 8, 14), (14, 16, 18, 19),
                              (14, 16, 18, 20), (14, 16, 18, 21),
                              (15, 14, 16, 17), (15, 14, 16, 18))
    atom_zero_improper_values = ()
    atom_i_improper_values = ((8, 16, 14, 15), (14, 18, 16, 17))
    expected_elems = None


class TestPRMNCRST(TOPBase):
    # Test case of PARM7 with no non-hydrogen including dihedrals
    ref_filename = PRMNCRST
    expected_n_atoms = 6
    expected_n_residues = 1
    ref_proteinatoms = 6
    expected_n_bonds = 5
    expected_n_angles = 7
    expected_n_dihedrals = 3
    expected_n_impropers = 0
    atom_i = 4
    expected_n_zero_bonds = 1
    expected_n_i_bonds = 2
    expected_n_zero_angles = 3
    expected_n_i_angles = 4
    expected_n_zero_dihedrals = 1
    expected_n_i_dihedrals = 3
    expected_n_zero_impropers = 0
    expected_n_i_impropers = 0
    atom_zero_bond_values = ((0, 1),)
    atom_i_bond_values = ((1, 4), (4, 5))
    atom_zero_angle_values = ((0, 1, 2), (0, 1, 3), (0, 1, 4))
    atom_i_angle_values = ((0, 1, 4), (1, 4, 5), (2, 1, 4), (3, 1, 4))
    atom_zero_dihedral_values = ((0, 1, 4, 5),)
    atom_i_dihedral_values = ((0, 1, 4, 5), (2, 1, 4, 5), (3, 1, 4, 5))
    atom_zero_improper_values = ()
    atom_i_improper_values = ()
    elems_ranges = [[0, 6], ]
    expected_elems = [np.array(["H", "C", "H", "H", "C", "O"], dtype=object), ]


class TestPRMNCRST_negative(TOPBase):
    # Same as above but with negative ATOMIC_NUMBER values (Issue 2306)
    ref_filename = PRMNEGATIVE
    expected_n_atoms = 6
    expected_n_residues = 1
    ref_proteinatoms = 6
    expected_n_bonds = 5
    expected_n_angles = 7
    expected_n_dihedrals = 3
    expected_n_impropers = 0
    atom_i = 4
    expected_n_zero_bonds = 1
    expected_n_i_bonds = 2
    expected_n_zero_angles = 3
    expected_n_i_angles = 4
    expected_n_zero_dihedrals = 1
    expected_n_i_dihedrals = 3
    expected_n_zero_impropers = 0
    expected_n_i_impropers = 0
    atom_zero_bond_values = ((0, 1),)
    atom_i_bond_values = ((1, 4), (4, 5))
    atom_zero_angle_values = ((0, 1, 2), (0, 1, 3), (0, 1, 4))
    atom_i_angle_values = ((0, 1, 4), (1, 4, 5), (2, 1, 4), (3, 1, 4))
    atom_zero_dihedral_values = ((0, 1, 4, 5),)
    atom_i_dihedral_values = ((0, 1, 4, 5), (2, 1, 4, 5), (3, 1, 4, 5))
    atom_zero_improper_values = ()
    atom_i_improper_values = ()
    elems_ranges = [[0, 6], ]
    expected_elems = [np.array(["H", "", "H", "H", "C", ""], dtype=object), ]


class TestPRMEP(TOPBase):
    # Issue #2449
    ref_filename = PRM19SBOPC
    expected_n_atoms = 46
    expected_n_residues = 9
    ref_proteinatoms = 23
    expected_n_bonds = 45
    expected_n_angles = 36
    expected_n_dihedrals = 41
    expected_n_impropers = 4
    atom_i = 25  # testing EPW atom
    expected_n_zero_bonds = 1
    expected_n_i_bonds = 1
    expected_n_zero_angles = 3
    expected_n_i_angles = 0
    expected_n_zero_dihedrals = 2
    expected_n_i_dihedrals = 0
    expected_n_zero_impropers = 0
    expected_n_i_impropers = 0
    atom_zero_bond_values = ((0, 1),)
    atom_i_bond_values = ((25, 22),)
    atom_zero_angle_values = ((0, 1, 2), (0, 1, 3), (0, 1, 4))
    atom_i_angle_values = ()
    atom_zero_dihedral_values = ((0, 1, 4, 5), (0, 1, 4, 6))
    atom_i_dihedral_values = ()
    atom_zero_improper_values = ()
    atom_i_improper_values = ()
    elems_ranges = [[0, 8], [20, 28]]
    expected_elems = [np.array(["H", "C", "H", "H", "C", "O", "N", "H"],
                      dtype=object),
                      np.array(["H", "H", "O", "H", "H", "", "O", "H"],
                      dtype=object)]


class TestErrorsAndWarnings(object):

    ATOMIC_NUMBER_MSG = (
        "ATOMIC_NUMBER record not found, elements attribute will not be populated"
    )
    MISSING_ELEM_MSG = (
        "Unknown ATOMIC_NUMBER value found for some atoms, "
        "these have been given an empty element record"
    )
    COORDINATE_READER_MSG = "No coordinate reader found"
    RESIDUE_CHAINID_MSG = (
        "Number of residues (38) does not match number of "
        "%RESIDUE_CHAINID (37). Skipping section."
    )

    @pytest.mark.parametrize(
        "parm,errmatch",
        (
            [PRMErr1, "%VE Missing in header"],
            [PRMErr2, "'TITLE' missing in header"],
            [PRMErr4, "Invalid header line."],
            [PRM_UreyBradley, "Chamber-style TOP file"],
        ),
    )
    def test_value_errors(self, parm, errmatch):
        with pytest.raises(ValueError, match=errmatch):
            u = mda.Universe(parm)

    def test_flag_index_error(self):
        with pytest.raises(IndexError, match="%FLAG section not found"):
            u = mda.Universe(PRMErr3)

    @pytest.mark.parametrize(
        "parm, errmsgs",
        (
            [PRM, [ATOMIC_NUMBER_MSG, COORDINATE_READER_MSG]],
            [PRM7, [ATOMIC_NUMBER_MSG, COORDINATE_READER_MSG]],
            [PRMpbc, [ATOMIC_NUMBER_MSG, COORDINATE_READER_MSG]],
            [PRMNEGATIVE, [MISSING_ELEM_MSG, COORDINATE_READER_MSG]],
            [PRM19SBOPC, [MISSING_ELEM_MSG, COORDINATE_READER_MSG]],
            [PRMErr5, [RESIDUE_CHAINID_MSG, COORDINATE_READER_MSG]],
        ),
    )
    def test_warning(self, parm, errmsgs):
        with pytest.warns(UserWarning) as record:
            u = mda.Universe(parm)

        messages = [rec.message.args[0] for rec in record]

        for msg in errmsgs:
            assert any(msg in recmsg for recmsg in messages)
