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

from numpy.testing import (
    assert_array_equal,
    assert_almost_equal,
    assert_equal,
    assert_,
    assert_raises,
)

import MDAnalysis
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
        self.u = MDAnalysis.Universe(PSF, DCD)
        self.a1 = self.u.atoms[1:3]
        self.a2 = self.u.atoms[3:5]

        self.TO1 = TopologyObject(self.a1.indices, self.u)
        self.TO2 = TopologyObject(self.a2.indices, self.u)
        # this atom only has one bond, so the order bonds are done (random because of sets)
        # won't come back to bite us
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
                     '<Bond between: Atom 12, Atom 9>')

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


