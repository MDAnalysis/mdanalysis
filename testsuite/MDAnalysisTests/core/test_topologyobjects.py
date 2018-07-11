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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from __future__ import absolute_import

import numpy as np
from numpy.testing import (
    assert_almost_equal,
    assert_equal,
)
import pytest

import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_bonds, calc_angles, calc_dihedrals
from MDAnalysis.core.topologyobjects import (
    TopologyGroup, TopologyObject, TopologyDict,
    # TODO: the following items are not used
    Bond, Angle, Dihedral, ImproperDihedral,
)


from MDAnalysisTests.datafiles import PSF, DCD, TRZ_psf, TRZ


@pytest.fixture(scope='module')
def PSFDCD():
    return mda.Universe(PSF, DCD)

class TestTopologyObjects(object):
    """Test the base TopologyObject funtionality

    init
    repr
    eq
    ne
    iter
    len
    """
    precision = 3  # see Issue #271 and #1556

    @staticmethod
    @pytest.fixture
    def a1(PSFDCD):
        return PSFDCD.atoms[1:3]

    @staticmethod
    @pytest.fixture
    def a2(PSFDCD):
        return PSFDCD.atoms[3:5]

    @staticmethod
    @pytest.fixture
    def TO1(a1):
        return TopologyObject(a1.indices, a1.universe)

    @staticmethod
    @pytest.fixture
    def TO2(a2):
        return TopologyObject(a2.indices, a2.universe)

    @staticmethod
    @pytest.fixture
    def b(PSFDCD):
        return PSFDCD.atoms[12].bonds[0]

    def test_repr(self, TO1):
        assert_equal(repr(TO1), '<TopologyObject between: Atom 1, Atom 2>')

    def test_eq(self, a1 ,TO1, TO2, PSFDCD):
        TO1_b = TopologyObject(a1.indices, PSFDCD)

        assert_equal(TO1 == TO1_b, True)
        assert_equal(TO1 == TO2, False)

    def test_ne(self, TO1, TO2, a1, PSFDCD):
        TO1_b = TopologyObject(a1.indices, PSFDCD)

        assert_equal(TO1 != TO1_b, False)
        assert_equal(TO1 != TO2, True)

    def test_gt(self, TO1, TO2):
        assert_equal(TO1 > TO2, False)

    def test_lt(self, TO1, TO2):
        assert_equal(TO1 < TO2, True)

    def test_iter(self, a1, TO1):
        assert_equal(list(a1), list(TO1))

    def test_len(self, a1):
        assert_equal(len(a1), 2)

    def test_hash(self, TO1, TO2, a1, PSFDCD):
        assert hash(TO1) != hash(TO2)

        # Different universe should yield different hash
        u = mda.Universe(PSF, DCD)
        ag = u.atoms[1:3]
        TO3 = TopologyObject(ag.indices, u)
        assert hash(TO1) != hash(TO3)

        # Different order should yield different hash
        u = mda.Universe(PSF, DCD)
        ag = u.atoms[[2, 1]]
        TO3 = TopologyObject(ag.indices, u)
        assert hash(TO1) != hash(TO3)

        # should work with TO from the same atomgroup
        TO3 = TopologyObject(a1.indices, PSFDCD)
        assert_equal(hash(TO1), hash(TO3))

    def test_indices(self, b):
        assert_equal(b.indices, tuple([a.index for a in b.atoms]))

    # Bond class checks
    def test_partner(self, b):
        a1, a2 = b
        assert_equal(b.partner(a1), a2)
        assert_equal(b.partner(a2), a1)

    def test_partner_VE(self, PSFDCD, b):
        a3 = PSFDCD.atoms[0]
        with pytest.raises(ValueError):
            b.partner(a3)

    def test_bondlength(self, b):
        assert_almost_equal(b.length(), 1.7661301556941993, self.precision)

    def test_bondrepr(self, b):
        assert_equal(repr(b), '<Bond between: Atom 9, Atom 12>')

    # Angle class checks
    def test_angle(self, PSFDCD):
        angle = PSFDCD.atoms[210].angles[0]

        assert_almost_equal(angle.angle(), 107.20893, self.precision)
        assert_almost_equal(angle.value(), 107.20893, self.precision)

    def test_angle_repr(self, PSFDCD):
        angle = PSFDCD.atoms[[30, 10, 20]].angle

        assert_equal(repr(angle), '<Angle between: Atom 20, Atom 10, Atom 30>')

    def test_angle_180(self):
        # we edit the coordinates, so make our own universe
        u = mda.Universe(PSF, DCD)
        angle = u.atoms[210].angles[0]
        coords = np.array([[1, 1, 1],
                           [2, 1, 1],
                           [3, 1, 1]],
                          dtype=np.float32)

        angle.atoms.positions = coords

        assert_almost_equal(angle.value(), 180.0, self.precision)

    # Dihedral class check
    def test_dihedral(self, PSFDCD):
        dihedral = PSFDCD.atoms[14].dihedrals[0]

        assert_almost_equal(dihedral.dihedral(), 18.317778, self.precision)
        assert_almost_equal(dihedral.value(), 18.317778, self.precision)

    def test_dihedral_repr(self, PSFDCD):
        dihedral = PSFDCD.atoms[[4, 7, 8, 1]].dihedral

        assert_equal(repr(dihedral),
                     '<Dihedral between: Atom 1, Atom 8, Atom 7, Atom 4>')

    # Improper_Dihedral class check
    def test_improper(self, PSFDCD):
        imp = PSFDCD.atoms[4].impropers[0]

        assert_almost_equal(imp.improper(), -3.8370631, self.precision)
        assert_almost_equal(imp.value(), -3.8370631, self.precision)

    def test_improper_repr(self, PSFDCD):
        imp = PSFDCD.atoms[[4, 7, 8, 1]].improper

        assert_equal(
            repr(imp),
            '<ImproperDihedral between: Atom 1, Atom 8, Atom 7, Atom 4>')

class TestTopologyGroup(object):
    """Tests TopologyDict and TopologyGroup classes with psf input"""

    @staticmethod
    @pytest.fixture
    def res1(PSFDCD):
        return PSFDCD.residues[0]

    @staticmethod
    @pytest.fixture
    def res2(PSFDCD):
        return PSFDCD.residues[1]

    @staticmethod
    @pytest.fixture
    def b_td(PSFDCD):
        return PSFDCD.atoms.bonds.topDict

    @staticmethod
    @pytest.fixture
    def a_td(PSFDCD):
        return PSFDCD.atoms.angles.topDict

    @staticmethod
    @pytest.fixture
    def t_td(PSFDCD):
        return PSFDCD.atoms.dihedrals.topDict

    # Checking TopologyDict functionality
    # * check that enough types are made
    # * check the identity of one of the keys
    # * check uniqueness of keys (ie reversed doesnt exist)
    # * then check same key reversed is accepted
    # * then select based on key
    # * then select based on reversed key and check output is the same
    # all for Bonds Angles and Dihedrals
    def test_td_len(self, b_td):
        assert len(b_td) == 57

    def test_td_iter(self, b_td):
        assert list(b_td) == list(b_td.dict.keys())

    def test_td_keyerror(self, b_td):
        with pytest.raises(KeyError):
            b_td[('something', 'stupid')]

    def test_td_universe(self, b_td, PSFDCD):
        assert b_td.universe is PSFDCD

    def test_bonds_types(self, PSFDCD, res1):
        """Tests TopologyDict for bonds"""
        assert len(PSFDCD.atoms.bonds.types()) == 57
        assert len(res1.atoms.bonds.types()) == 12

    def test_bonds_contains(self, b_td):
        assert ('57', '2') in b_td

    def test_bond_uniqueness(self, PSFDCD):
        bondtypes = PSFDCD.atoms.bonds.types()
        # check that a key doesn't appear in reversed format in keylist
        # have to exclude case of b[::-1] == b as this is false positive
        assert not any([b[::-1] in bondtypes
                        for b in bondtypes if b[::-1] != b])


    def test_bond_reversal(self, PSFDCD, b_td):
        bondtypes = PSFDCD.atoms.bonds.types()
        b = bondtypes[1]
        assert all([b in b_td, b[::-1] in b_td])

        tg1 = b_td[b]
        tg2 = b_td[b[::-1]]
        assert tg1 == tg2

    def test_angles_types(self, PSFDCD):
        """TopologyDict for angles"""
        assert len(PSFDCD.atoms.angles.types()) == 130

    def test_angles_contains(self, a_td):
        assert ('23', '73', '1') in a_td

    def test_angles_uniqueness(self, a_td):
        bondtypes = a_td.keys()
        assert not any(b[::-1] in bondtypes
                       for b in bondtypes if b[::-1] != b)

    def test_angles_reversal(self, a_td):
        bondtypes = list(a_td.keys())
        b = bondtypes[1]
        assert all([b in a_td, b[::-1] in a_td])

        tg1 = a_td[b]
        tg2 = a_td[b[::-1]]
        assert tg1 == tg2

    def test_dihedrals_types(self, PSFDCD):
        """TopologyDict for dihedrals"""
        assert len(PSFDCD.atoms.dihedrals.types()) == 220

    def test_dihedrals_contains(self, t_td):
        assert ('30', '29', '20', '70') in t_td

    def test_dihedrals_uniqueness(self, t_td):
        bondtypes = t_td.keys()
        assert not any(b[::-1] in bondtypes for b in bondtypes if b[::-1] != b)

    def test_dihedrals_reversal(self, t_td):
        bondtypes = list(t_td.keys())
        b = bondtypes[1]
        assert all([b in t_td, b[::-1] in t_td])

        tg1 = t_td[b]
        tg2 = t_td[b[::-1]]
        assert tg1 == tg2

    def test_bad_creation(self):
        """Test making a TopologyDict out of nonsense"""
        inputlist = ['a', 'b', 'c']
        with pytest.raises(TypeError):
            TopologyDict(inputlist)

    def test_bad_creation_TG(self):
        """Test making a TopologyGroup out of nonsense"""
        inputlist = ['a', 'b', 'c']
        with pytest.raises(TypeError):
            TopologyGroup(inputlist)

    def test_tg_creation_bad_btype(self, PSFDCD):
        vals = np.array([[0, 10], [5, 15]])
        with pytest.raises(ValueError):
            TopologyGroup(vals, PSFDCD, btype='apple')

    def test_bond_tg_creation_notype(self, PSFDCD):
        vals = np.array([[0, 10], [5, 15]])
        tg = TopologyGroup(vals, PSFDCD)

        assert tg.btype == 'bond'
        assert_equal(tg[0].indices, (0, 10))
        assert_equal(tg[1].indices, (5, 15))

    def test_angle_tg_creation_notype(self, PSFDCD):
        vals = np.array([[0, 5, 10], [5, 10, 15]])
        tg = TopologyGroup(vals, PSFDCD)

        assert tg.btype == 'angle'
        assert_equal(tg[0].indices, (0, 5, 10))
        assert_equal(tg[1].indices, (5, 10, 15))

    def test_dihedral_tg_creation_notype(self, PSFDCD):
        vals = np.array([[0, 2, 4, 6], [5, 7, 9, 11]])

        tg = TopologyGroup(vals, PSFDCD)

        assert tg.btype == 'dihedral'
        assert_equal(tg[0].indices, (0, 2, 4, 6))
        assert_equal(tg[1].indices, (5, 7, 9, 11))

    def test_create_guessed_tg(self, PSFDCD):
        vals = np.array([[0, 10], [5, 15]])
        tg = TopologyGroup(vals, PSFDCD, guessed=True)

        assert_equal(tg._guessed, np.array([[True], [True]]))

    def test_create_guessed_tg_2(self, PSFDCD):
        vals = np.array([[0, 10], [5, 15]])
        tg = TopologyGroup(vals, PSFDCD, guessed=False)

        assert_equal(tg._guessed, np.array([[False], [False]]))

    def test_TG_equality(self, PSFDCD):
        """Make two identical TGs,
        * check they're equal
        * change one very slightly and see if they notice
        """
        tg = PSFDCD.atoms.bonds.selectBonds(('23', '3'))
        tg2 = PSFDCD.atoms.bonds.selectBonds(('23', '3'))

        assert tg == tg2

        tg3 = PSFDCD.atoms.bonds.selectBonds(('81', '10'))

        assert not (tg == tg3)
        assert tg != tg3

    def test_create_TopologyGroup(self, res1, PSFDCD):
        res1_tg = res1.atoms.bonds.select_bonds(('23', '3'))  # make a tg
        assert len(res1_tg) == 4  # check size of tg
        testbond = PSFDCD.atoms[7].bonds[0]
        assert testbond in res1_tg  # check a known bond is present

        res1_tg2 = res1.atoms.bonds.select_bonds(('23', '3'))
        assert res1_tg == res1_tg2

    @pytest.mark.parametrize('attr',
                             ['bonds', 'angles', 'dihedrals', 'impropers'])
    def test_TG_loose_intersection(self, PSFDCD, attr):
        """Pull bonds from a TG which are at least partially in an AG"""
        ag = PSFDCD.atoms[10:60]

        TG = getattr(PSFDCD.atoms, attr)
        ref = getattr(ag, attr)
        newTG = TG.atomgroup_intersection(ag)

        assert newTG == ref

    def test_TG_strict_intersection(self, PSFDCD):
        """Pull bonds from TG which are fully in an AG"""

        def check_strict_intersection(topg, atomg):
            new_topg = topg.atomgroup_intersection(atomg, strict=True)

            return all([all([a in atomg for a in b.atoms]) for b in new_topg])

        def manual(topg, atomg):
            """Gives a set of the Bonds that should be found"""
            if len(atomg) == 1:  # hack for Atom input
                atomg = [atomg]
            man = []
            for b in topg:
                if all([a in atomg for a in b.atoms]):
                    man.append(b)

            return set(man)

        testinput = PSFDCD.atoms[10:60]

        # bonds
        assert check_strict_intersection(PSFDCD.atoms.bonds, testinput)
        assert (manual(PSFDCD.atoms.bonds, testinput) ==
                set(PSFDCD.atoms.bonds.atomgroup_intersection(
                    testinput, strict=True)))
        # angles
        assert check_strict_intersection(PSFDCD.atoms.angles, testinput)
        assert (manual(PSFDCD.atoms.angles, testinput) ==
                set(PSFDCD.atoms.angles.atomgroup_intersection(
                    testinput, strict=True)))
        # dihedrals
        assert check_strict_intersection(PSFDCD.atoms.dihedrals, testinput)
        assert (manual(PSFDCD.atoms.dihedrals, testinput) ==
                set(PSFDCD.atoms.dihedrals.atomgroup_intersection(
                    testinput, strict=True)))

    def test_add_TopologyGroups(self, res1, res2, PSFDCD):
        res1_tg = res1.atoms.bonds.selectBonds(('23', '3'))
        res2_tg = res2.atoms.bonds.selectBonds(('23', '3'))

        combined_tg = res1_tg + res2_tg  # add tgs together
        assert len(combined_tg) == 10

        big_tg = PSFDCD.atoms.bonds.selectBonds(('23', '3'))
        big_tg += combined_tg  # try and add some already included bonds
        assert len(big_tg) == 494  # check len doesn't change

    def test_add_empty_to_TG(self, PSFDCD):
        tg1 = PSFDCD.bonds[10:15]
        tg2 = PSFDCD.bonds[:0]  # empty
        tg3 = tg1 + tg2

        assert tg1 == tg3

    def test_add_TO_to_empty_TG(self, PSFDCD):
        tg1 = PSFDCD.bonds[:0]  # empty
        to = PSFDCD.bonds[5]
        tg3 = tg1 + to

        assert_equal(tg3.indices,  to.indices[None, :])

    def test_add_TG_to_empty_TG(self, PSFDCD):
        tg1 = PSFDCD.bonds[:0]  # empty
        tg2 = PSFDCD.bonds[5:7]
        tg3 = tg1 + tg2

        assert tg2 == tg3

    def test_add_singleitem(self, PSFDCD):
        tg = PSFDCD.atoms.bonds[:10]
        to = PSFDCD.atoms.bonds[55]

        assert len(tg + to) == 11

    def test_add_wrongtype_TopologyGroup(self, PSFDCD):
        tg = PSFDCD.atoms.bonds[:10]  # TG of bonds
        ang = PSFDCD.atoms.angles[10]  # single angle
        angg = PSFDCD.atoms.angles[:10]  # TG of angles

        for other in [ang, angg]:
            with pytest.raises(TypeError):
                this = tg + other

    def test_bad_add_TopologyGroup(self, PSFDCD):
        tg = PSFDCD.atoms.bonds[:10]  # TopologyGroup
        ag = PSFDCD.atoms[:10]  # AtomGroup
        with pytest.raises(TypeError):
            this = tg + ag

    def test_TG_getitem_single(self, PSFDCD):
        tg = PSFDCD.atoms.bonds[:10]
        bondlist = list(tg)
        bond = tg[0]

        assert bond == bondlist[0]

    def test_TG_getitem_slice(self, PSFDCD):
        tg = PSFDCD.atoms.bonds[:10]
        tg2 = tg[1:4]

        assert_equal(tg2.indices, tg.indices[1:4])

    def test_TG_getitem_fancy(self, PSFDCD):
        tg = PSFDCD.atoms.bonds[:10]
        tg2 = tg[[1, 4, 5]]

        manual = TopologyGroup(tg.indices[[1, 4, 5]],
                               tg.universe, tg.btype)

        assert list(tg2) == list(manual)

    def test_TG_getitem_bool(self, PSFDCD):
        # Issue #282
        sel = np.array([True, False, True])
        tg = PSFDCD.atoms.bonds[:3]
        tg2 = tg[sel]

        assert len(tg2) == 2
        for b in [tg[0], tg[2]]:
            assert b in tg2

    def test_TG_getitem_bool_IE(self, PSFDCD):
        sel = []
        tg = PSFDCD.atoms.bonds[10:13]
        tg2 = tg[sel]

        assert len(tg2) == 0

    # atomX access
    def test_atom1(self, PSFDCD):
        tg = PSFDCD.bonds[:5]
        a1 = tg.atom1

        assert len(tg) == len(a1)
        for (atom, bond) in zip(a1, tg):
            assert atom == bond[0]

    def test_atom2(self, PSFDCD):
        tg = PSFDCD.bonds[:5]
        a2 = tg.atom2

        assert len(tg) == len(a2)
        for (atom, bond) in zip(a2, tg):
            assert atom == bond[1]

    def test_atom3_IE(self, PSFDCD):
        tg = PSFDCD.bonds[:5]
        with pytest.raises(IndexError):
            tg.atom3

    def test_atom3(self, PSFDCD):
        tg = PSFDCD.angles[:5]
        a3 = tg.atom3
        assert len(tg) == len(a3)
        for (atom, bond) in zip(a3, tg):
            assert atom == bond[2]

    def test_atom4_IE(self, PSFDCD):
        tg = PSFDCD.bonds[:5]
        with pytest.raises(IndexError):
            tg.atom4

    def test_atom4(self, PSFDCD):
        tg = PSFDCD.dihedrals[:5]

        a4 = tg.atom4
        assert len(tg) == len(a4)
        for (atom, bond) in zip(a4, tg):
            assert atom == bond[3]


class TestTopologyGroup_Cython(object):
    """
    Check that the shortcut to all cython functions:
     - work (return proper values)
     - catch errors
    """
    @staticmethod
    @pytest.fixture
    def bgroup(PSFDCD):
        return PSFDCD.atoms[:5].bonds

    @staticmethod
    @pytest.fixture
    def agroup(PSFDCD):
        return PSFDCD.atoms[:5].angles

    @staticmethod
    @pytest.fixture
    def dgroup(PSFDCD):
        return PSFDCD.atoms[:5].dihedrals

    @staticmethod
    @pytest.fixture
    def igroup(PSFDCD):
        return PSFDCD.atoms[:5].impropers

    # bonds
    def test_wrong_type_bonds(self, agroup, dgroup, igroup):
        for tg in [agroup, dgroup, igroup]:
            with pytest.raises(TypeError):
                tg.bonds()

    def test_right_type_bonds(self, bgroup, PSFDCD):
        assert_equal(bgroup.bonds(),
                     calc_bonds(bgroup.atom1.positions,
                                bgroup.atom2.positions))
        assert_equal(bgroup.bonds(pbc=True),
                     calc_bonds(bgroup.atom1.positions,
                                bgroup.atom2.positions,
                                box=PSFDCD.dimensions))
        assert_equal(bgroup.values(),
                     calc_bonds(bgroup.atom1.positions,
                                bgroup.atom2.positions))
        assert_equal(bgroup.values(pbc=True),
                     calc_bonds(bgroup.atom1.positions,
                                bgroup.atom2.positions,
                                box=PSFDCD.dimensions))

    # angles
    def test_wrong_type_angles(self, bgroup, dgroup, igroup):
        for tg in [bgroup, dgroup, igroup]:
            with pytest.raises(TypeError):
                tg.angles()

    def test_right_type_angles(self, agroup, PSFDCD):
        assert_equal(agroup.angles(),
                     calc_angles(agroup.atom1.positions,
                                 agroup.atom2.positions,
                                 agroup.atom3.positions))
        assert_equal(agroup.angles(pbc=True),
                     calc_angles(agroup.atom1.positions,
                                 agroup.atom2.positions,
                                 agroup.atom3.positions,
                                 box=PSFDCD.dimensions))
        assert_equal(agroup.values(),
                     calc_angles(agroup.atom1.positions,
                                 agroup.atom2.positions,
                                 agroup.atom3.positions))
        assert_equal(agroup.values(pbc=True),
                     calc_angles(agroup.atom1.positions,
                                 agroup.atom2.positions,
                                 agroup.atom3.positions,
                                 box=PSFDCD.dimensions))

    # dihedrals & impropers
    def test_wrong_type_dihedrals(self, bgroup, agroup):
        for tg in [bgroup, agroup]:
            with pytest.raises(TypeError):
                tg.dihedrals()

    def test_right_type_dihedrals(self, dgroup, PSFDCD):
        assert_equal(dgroup.dihedrals(),
                     calc_dihedrals(dgroup.atom1.positions,
                                   dgroup.atom2.positions,
                                   dgroup.atom3.positions,
                                   dgroup.atom4.positions))
        assert_equal(dgroup.dihedrals(pbc=True),
                     calc_dihedrals(dgroup.atom1.positions,
                                   dgroup.atom2.positions,
                                   dgroup.atom3.positions,
                                   dgroup.atom4.positions,
                                   box=PSFDCD.dimensions))
        assert_equal(dgroup.values(),
                     calc_dihedrals(dgroup.atom1.positions,
                                   dgroup.atom2.positions,
                                   dgroup.atom3.positions,
                                   dgroup.atom4.positions))
        assert_equal(dgroup.values(pbc=True),
                     calc_dihedrals(dgroup.atom1.positions,
                                   dgroup.atom2.positions,
                                   dgroup.atom3.positions,
                                   dgroup.atom4.positions,
                                   box=PSFDCD.dimensions))

    def test_right_type_impropers(self, igroup, PSFDCD):
        assert_equal(igroup.dihedrals(),
                     calc_dihedrals(igroup.atom1.positions,
                                   igroup.atom2.positions,
                                   igroup.atom3.positions,
                                   igroup.atom4.positions))
        assert_equal(igroup.dihedrals(pbc=True),
                     calc_dihedrals(igroup.atom1.positions,
                                   igroup.atom2.positions,
                                   igroup.atom3.positions,
                                   igroup.atom4.positions,
                                   box=PSFDCD.dimensions))
        assert_equal(igroup.values(),
                     calc_dihedrals(igroup.atom1.positions,
                                   igroup.atom2.positions,
                                   igroup.atom3.positions,
                                   igroup.atom4.positions))
        assert_equal(igroup.values(pbc=True),
                     calc_dihedrals(igroup.atom1.positions,
                                   igroup.atom2.positions,
                                   igroup.atom3.positions,
                                   igroup.atom4.positions,
                                   box=PSFDCD.dimensions))


def test_bond_length_pbc():
    u = mda.Universe(TRZ_psf, TRZ)

    ref = u.bonds[0].length()

    # move an atom a box width in all dimensions
    u.atoms[0].position += u.dimensions[:3]

    assert_almost_equal(ref, u.bonds[0].length(pbc=True), decimal=6)

def test_cross_universe_eq():
    u1 = mda.Universe(PSF)
    u2 = mda.Universe(PSF)

    assert not (u1.bonds[0] == u2.bonds[0])

def test_zero_size_TG_indices_bonds():
    u = mda.Universe.empty(10)

    u.add_TopologyAttr('bonds', values=[(1, 2), (2, 3)])

    ag = u.atoms[[0]]

    idx = ag.bonds.to_indices()

    assert idx.shape == (0, 2)
    assert idx.dtype == np.int32

def test_zero_size_TG_indices_angles():
    u = mda.Universe.empty(10)

    u.add_TopologyAttr('angles', values=[(1, 2, 3), (2, 3, 4)])

    ag = u.atoms[[0]]

    idx = ag.angles.to_indices()

    assert idx.shape == (0, 3)
    assert idx.dtype == np.int32
