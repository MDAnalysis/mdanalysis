# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
# Beckstein and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

from numpy.testing import (
    assert_,
    assert_raises,
)


import MDAnalysis as mda
from MDAnalysis.core.topologyobjects import (
    Bond,
    Angle,
    Dihedral,
    ImproperDihedral,
)

from MDAnalysis.tests.datafiles import (PSF, DCD)


class TestAtomGroupToTopology(object):
    """Test the conversion of AtomGroup to TopologyObjects"""
    def setUp(self):
        self.u = mda.Universe(PSF, DCD)

    def tearDown(self):
        del self.u

    def test_bond(self):
        ag = self.u.atoms[:2]
        bond = ag.bond
        assert_(isinstance(bond, Bond))

    def test_angle(self):
        ag = self.u.atoms[:3]
        angle = ag.angle
        assert_(isinstance(angle, Angle))

    def test_dihedral(self):
        ag = self.u.atoms[:4]
        dih = ag.dihedral
        assert_(isinstance(dih, Dihedral))

    def test_improper(self):
        ag = self.u.atoms[:4]
        imp = ag.improper
        assert_(isinstance(imp, ImproperDihedral))

    def _check_VE(self, btype):
        ag = self.u.atoms[:10]

        assert_raises(ValueError, getattr, ag, btype)

    def test_VEs(self):
        for btype in ('bond', 'angle', 'dihedral', 'improper'):
            yield self._check_VE, btype

