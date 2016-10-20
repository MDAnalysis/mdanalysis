# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

import numpy as np
from numpy.testing import (
    assert_array_equal,
    assert_almost_equal,
    assert_equal,
    assert_,
    assert_raises,
)

import MDAnalysis as mda
from MDAnalysis.core.topologyobjects import (
    TopologyGroup, TopologyObject, TopologyDict,
    Bond, Angle, Dihedral, ImproperDihedral,
)


from MDAnalysisTests.datafiles import PSF, DCD


class TestTopologyObjects(object):
    """Test the base TopologyObject funtionality

    init
    repr
    eq
    ne
    iter
    len
    """

    def setUp(self):
        self.precision = 3  # rather lenient but see #271
        self.u = mda.Universe(PSF, DCD)
        self.a1 = self.u.atoms[1:3]
        self.a2 = self.u.atoms[3:5]

        self.TO1 = TopologyObject(self.a1.indices, self.u)
        self.TO2 = TopologyObject(self.a2.indices, self.u)
        # this atom only has one bond, so the order bonds are done (random
        # because of sets) won't come back to bite us
        self.b = self.u.atoms[12].bonds[0]

    def tearDown(self):
        del self.u
        del self.a1
        del self.a2
        del self.TO1
        del self.TO2
        del self.b

    def test_repr(self):
        assert_equal(repr(self.TO1),
                     '<TopologyObject between: Atom 1, Atom 2>')

    def test_eq(self):
        TO1_b = TopologyObject(self.a1.indices, self.u)

        assert_equal(self.TO1 == TO1_b, True)
        assert_equal(self.TO1 == self.TO2, False)

    def test_ne(self):
        TO1_b = TopologyObject(self.a1.indices, self.u)

        assert_equal(self.TO1 != TO1_b, False)
        assert_equal(self.TO1 != self.TO2, True)

    def test_gt(self):
        assert_equal(self.TO1 > self.TO2, False)

    def test_lt(self):
        assert_equal(self.TO1 < self.TO2, True)

    def test_iter(self):
        assert_equal(list(self.a1), list(self.TO1))

    def test_len(self):
        assert_equal(len(self.a1), 2)

    def test_hash(self):
        assert_(hash(self.TO1) != hash(self.TO2))

        # Different universe should yield different hash
        u = mda.Universe(PSF, DCD)
        ag = u.atoms[1:3]
        TO3 = TopologyObject(ag.indices, u)
        assert_(hash(self.TO1) != hash(TO3))

        # Different order should yield different hash
        u = mda.Universe(PSF, DCD)
        ag = u.atoms[[2, 1]]
        TO3 = TopologyObject(ag.indices, u)
        assert_(hash(self.TO1) != hash(TO3))

        # should work with TO from the same atomgroup
        TO3 = TopologyObject(self.a1.indices, self.u)
        assert_equal(hash(self.TO1), hash(TO3))

    def test_indices(self):
        assert_equal(self.b.indices, tuple([b.index for b in self.b.atoms]))

    # Bond class checks
    def test_partner(self):
        a1, a2 = self.b
        assert_equal(self.b.partner(a1), a2)
        assert_equal(self.b.partner(a2), a1)

    def test_partner_VE(self):
        a3 = self.u.atoms[0]
        assert_raises(ValueError, self.b.partner, a3)

    def test_bondlength(self):
        assert_almost_equal(self.b.length(), 1.7661301556941993, self.precision)

    def test_bondrepr(self):
        assert_equal(repr(self.b),
                     '<Bond between: Atom 9, Atom 12>')

    # Angle class checks
    def test_angle(self):
        angle = self.u.atoms[210].angles[0]

        assert_almost_equal(angle.angle(), 107.20893, self.precision)
        assert_almost_equal(angle.value(), 107.20893, self.precision)

    # Dihedral class check
    def test_dihedral(self):
        dihedral = self.u.atoms[14].dihedrals[0]

        assert_almost_equal(dihedral.dihedral(), 18.317778, self.precision)
        assert_almost_equal(dihedral.value(), 18.317778, self.precision)

    # Improper_Dihedral class check
    def test_improper(self):
        imp = self.u.atoms[4].impropers[0]

        assert_almost_equal(imp.improper(), -3.8370631, self.precision)
        assert_almost_equal(imp.value(), -3.8370631, self.precision)


class TestTopologyGroup(object):
    """Tests TopologyDict and TopologyGroup classes with psf input"""

    def setUp(self):
        topology = PSF
        self.universe = mda.Universe(topology)
        # force the loading of topology
        self.res1 = self.universe.residues[0]
        self.res2 = self.universe.residues[1]
        # topologydicts for testing
        self.b_td = self.universe.atoms.bonds.topDict
        self.a_td = self.universe.atoms.angles.topDict
        self.t_td = self.universe.atoms.dihedrals.topDict

    def tearDown(self):
        del self.universe
        del self.res1
        del self.res2
        del self.b_td
        del self.a_td
        del self.t_td

    # Checking TopologyDict functionality
    # * check that enough types are made
    # * check the identity of one of the keys
    # * check uniqueness of keys (ie reversed doesnt exist)
    # * then check same key reversed is accepted
    # * then select based on key
    # * then select based on reversed key and check output is the same
    # all for Bonds Angles and Dihedrals
    def test_td_len(self):
        assert_equal(len(self.b_td), 57)

    def test_td_iter(self):
        assert_equal(list(self.b_td), list(self.b_td.dict.keys()))

    def test_td_keyerror(self):
        assert_raises(KeyError, self.b_td.__getitem__, ('something', 'stupid'))

    def test_bonds_types(self):
        """Tests TopologyDict for bonds"""
        assert_equal(len(self.universe.atoms.bonds.types()), 57)
        assert_equal(len(self.res1.atoms.bonds.types()), 12)

    def test_bonds_contains(self):
        assert_equal(('57', '2') in self.b_td, True)

    def test_bond_uniqueness(self):
        bondtypes = self.universe.atoms.bonds.types()
        # check that a key doesn't appear in reversed format in keylist
        # have to exclude case of b[::-1] == b as this is false positive
        assert_equal(any([b[::-1] in bondtypes for b in bondtypes if b[::-1] != b]),
                     False)

    def test_bond_reversal(self):
        bondtypes = self.universe.atoms.bonds.types()
        b = bondtypes[1]
        assert_equal(all([b in self.b_td, b[::-1] in self.b_td]), True)

        tg1 = self.b_td[b]
        tg2 = self.b_td[b[::-1]]
        assert_equal(tg1 == tg2, True)

    def test_angles_types(self):
        """TopologyDict for angles"""
        assert_equal(len(self.universe.atoms.angles.types()), 130)

    def test_angles_contains(self):
        assert_equal(('23', '73', '1') in self.a_td, True)

    def test_angles_uniqueness(self):
        bondtypes = self.a_td.keys()
        assert_equal(any([b[::-1] in bondtypes for b in bondtypes if b[::-1] != b]),
                     False)

    def test_angles_reversal(self):
        bondtypes = self.a_td.keys()
        b = bondtypes[1]
        assert_equal(all([b in self.a_td, b[::-1] in self.a_td]), True)

        tg1 = self.a_td[b]
        tg2 = self.a_td[b[::-1]]
        assert_equal(tg1 == tg2, True)

    def test_dihedrals_types(self):
        """TopologyDict for dihedrals"""
        assert_equal(len(self.universe.atoms.dihedrals.types()), 220)

    def test_dihedrals_contains(self):
        assert_equal(('30', '29', '20', '70') in self.t_td, True)

    def test_dihedrals_uniqueness(self):
        bondtypes = self.t_td.keys()
        assert_equal(any([b[::-1] in bondtypes for b in bondtypes if b[::-1] != b]),
                     False)

    def test_dihedrals_reversal(self):
        bondtypes = self.t_td.keys()
        b = bondtypes[1]
        assert_equal(all([b in self.t_td, b[::-1] in self.t_td]), True)

        tg1 = self.t_td[b]
        tg2 = self.t_td[b[::-1]]
        assert_equal(tg1 == tg2, True)

    def test_bad_creation(self):
        """Test making a TopologyDict out of nonsense"""
        inputlist = ['a', 'b', 'c']
        assert_raises(TypeError, TopologyDict, inputlist)

    def test_bad_creation_TG(self):
        """Test making a TopologyGroup out of nonsense"""
        inputlist = ['a', 'b', 'c']
        assert_raises(TypeError, TopologyGroup, inputlist)

    def test_TG_equality(self):
        """Make two identical TGs,
        * check they're equal
        * change one very slightly and see if they notice
        """
        tg = self.universe.atoms.bonds.selectBonds(('23', '3'))
        tg2 = self.universe.atoms.bonds.selectBonds(('23', '3'))

        assert_equal(tg == tg2, True)

        tg3 = self.universe.atoms.bonds.selectBonds(('81', '10'))
        assert_equal(tg == tg3, False)
        assert_equal(tg != tg3, True)

    def test_create_TopologyGroup(self):
        res1_tg = self.res1.atoms.bonds.select_bonds(('23', '3'))  # make a tg
        assert_equal(len(res1_tg), 4)  # check size of tg
        testbond = self.universe.atoms[7].bonds[0]
        assert_equal(testbond in res1_tg, True)  # check a known bond is present

        res1_tg2 = self.res1.atoms.bonds.select_bonds(('23', '3'))
        assert_equal(res1_tg == res1_tg2, True)

    def test_TG_loose_intersection(self):
        """Pull bonds from a TG which are at least partially in an AG"""
        u = self.universe
        ag = self.universe.atoms[10:60]

        # Check that every bond has at least one atom in the atomgroup
        for TG, ref in [(u.atoms.bonds, ag.bonds),
                        (u.atoms.angles, ag.angles),
                        (u.atoms.dihedrals, ag.dihedrals),
                        (u.atoms.impropers, ag.impropers)]:
            newTG = TG.atomgroup_intersection(ag)

            assert_(newTG == ref,
                    "Loose intersection failed with: " + TG.btype)

    def test_TG_strict_intersection(self):
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

        u = self.universe
        testinput = self.universe.atoms[10:60]

        # bonds
        assert_(check_strict_intersection(u.atoms.bonds, testinput))
        assert_equal(manual(u.atoms.bonds, testinput),
                     set(u.atoms.bonds.atomgroup_intersection(
                         testinput, strict=True)))
        # angles
        assert_(check_strict_intersection(u.atoms.angles, testinput))
        assert_equal(manual(u.atoms.angles, testinput),
                     set(u.atoms.angles.atomgroup_intersection(
                         testinput, strict=True)))
        # dihedrals
        assert_(check_strict_intersection(u.atoms.dihedrals, testinput))
        assert_equal(manual(u.atoms.dihedrals, testinput),
                     set(u.atoms.dihedrals.atomgroup_intersection(
                         testinput, strict=True)))

    def test_add_TopologyGroups(self):
        res1_tg = self.res1.atoms.bonds.selectBonds(('23', '3'))
        res2_tg = self.res2.atoms.bonds.selectBonds(('23', '3'))

        combined_tg = res1_tg + res2_tg  # add tgs together
        assert_equal(len(combined_tg), 10)

        big_tg = self.universe.atoms.bonds.selectBonds(('23', '3'))

        big_tg += combined_tg  # try and add some already included bonds
        assert_equal(len(big_tg), 494)  # check len doesn't change

    def test_add_singleitem(self):
        tg = self.universe.atoms.bonds[:10]
        to = self.universe.atoms.bonds[55]

        assert_equal(len(tg + to), 11)

    def test_add_wrongtype_TopologyGroup(self):
        def adder(a, b):
            return a + b

        tg = self.universe.atoms.bonds[:10]  # TG of bonds

        ang = self.universe.atoms.angles[10]  # single angle
        angg = self.universe.atoms.angles[:10]  # TG of angles

        for other in [ang, angg]:
            assert_raises(TypeError, adder, tg, other)

    def test_bad_add_TopologyGroup(self):
        def adder(a, b):
            return a + b

        tg = self.universe.atoms.bonds[:10]  # TopologyGroup

        ag = self.universe.atoms[:10]  # AtomGroup

        assert_raises(TypeError, adder, tg, ag)

    def test_TG_getitem_single(self):
        tg = self.universe.atoms.bonds[:10]

        bondlist = list(tg)
        bond = tg[0]

        assert_equal(bond, bondlist[0])

    def test_TG_getitem_slice(self):
        tg = self.universe.atoms.bonds[:10]

        tg2 = tg[1:4]

        assert_array_equal(tg2.indices, tg.indices[1:4])

    def test_TG_getitem_fancy(self):
        tg = self.universe.atoms.bonds[:10]

        tg2 = tg[[1, 4, 5]]

        manual = TopologyGroup(tg.indices[[1, 4, 5]],
                               tg.universe, tg.btype)

        assert_equal(list(tg2), list(manual))

    def test_TG_getitem_bool(self):
        # Issue #282
        sel = np.array([True, False, True])
        tg = self.universe.atoms.bonds[:3]
        tg2 = tg[sel]
        assert_equal(len(tg2), 2)
        for b in [tg[0], tg[2]]:
            assert_equal(b in tg2, True)

    def test_TG_getitem_bool_IE(self):
        sel = []
        tg = self.universe.atoms.bonds[10:13]
        tg2 = tg[sel]
        assert_equal(len(tg2), 0)
