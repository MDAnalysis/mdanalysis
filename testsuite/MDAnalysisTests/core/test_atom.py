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
import numpy as np
from numpy.testing import (
    dec,
    raises,
    assert_almost_equal,
    assert_equal,
    assert_raises,
)

import MDAnalysis as mda
from MDAnalysis import NoDataError

from MDAnalysisTests.datafiles import (
    PSF, DCD,
    XYZ_mini,
)
from MDAnalysisTests import parser_not_found


class TestAtom(object):
    # Legacy tests from before 363
    """Tests of Atom."""

    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def setUp(self):
        """Set up the standard AdK system in implicit solvent."""
        self.universe = mda.Universe(PSF, DCD)
        self.atom = self.universe.atoms[1000]  # Leu67:CG
        self.known_pos = np.array([3.94543672, -12.4060812, -7.26820087],
                                  dtype=np.float32)

    def tearDown(self):
        del self.universe
        del self.atom
        del self.known_pos

    def test_attributes_names(self):
        a = self.atom
        assert_equal(a.name, 'CG')
        assert_equal(a.resname, 'LEU')

    def test_setting_attribute_name(self):
        self.atom.name = 'AA'
        assert_equal(self.atom.name, 'AA')

    def test_setting_attribute_type(self):
        self.atom.type = 'Z'
        assert_equal(self.atom.type, 'Z')

    def test_setting_attribute_mass(self):
        self.atom.mass = 13
        assert_equal(self.atom.mass, 13)

    def test_setting_attributes_charge(self):
        self.atom.charge = 6
        assert_equal(self.atom.charge, 6)

    def test_attributes_positions(self):
        a = self.atom
        # new position property (mutable)
        assert_almost_equal(a.position, self.known_pos)
        pos = a.position + 3.14
        a.position = pos
        assert_almost_equal(a.position, pos)

    def test_atom_selection(self):
        asel = self.universe.select_atoms('atom 4AKE 67 CG').atoms[0]
        assert_equal(self.atom, asel)

    def test_hierarchy(self):
        u = self.universe
        a = self.atom
        assert_equal(a.segment, u.s4AKE)
        assert_equal(a.residue, u.residues[66])

    def test_bad_add(self):
        def bad_add():
            return self.atom + 1

        assert_raises(TypeError, bad_add)

    def test_add_AG(self):
        ag = self.universe.atoms[:2]

        ag2 = self.atom + ag

        for at in [self.atom, ag[0], ag[1]]:
            assert_equal(at in ag2, True)

    def test_no_velo(self):
        def lookup_velo():
            return self.atom.velocity

        assert_raises(NoDataError, lookup_velo)

    def test_bonded_atoms(self):
        at = self.universe.atoms[0]
        ref = [b.partner(at) for b in at.bonds]
        assert_equal(ref, list(at.bonded_atoms))

    @raises(AttributeError)
    def test_undefined_occupancy(self):
        self.universe.atoms[0].occupancy


class TestAtomNoForceNoVel(object):
    def setUp(self):
        self.u = mda.Universe(XYZ_mini)
        self.a = self.u.atoms[0]

    def tearDown(self):
        del self.u

    def test_velocity_fail(self):
        assert_raises(NoDataError, getattr, self.a, 'velocity')

    def test_force_fail(self):
        assert_raises(NoDataError, getattr, self.a, 'force')

    def test_velocity_set_fail(self):
        assert_raises(NoDataError, setattr, self.a, 'velocity',
                      [1.0, 1.0, 1.0])

    def test_force_set_fail(self):
        assert_raises(NoDataError, setattr, self.a, 'force', [1.0, 1.0, 1.0])


