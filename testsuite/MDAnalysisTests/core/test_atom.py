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
import pickle
import pytest
import numpy as np

import MDAnalysis as mda
from MDAnalysis import NoDataError
from MDAnalysisTests.datafiles import (
    PSF, DCD,
    XYZ_mini,
)
from numpy.testing import assert_almost_equal


class TestAtom(object):
    # Legacy tests from before 363
    """Tests of Atom."""

    """Set up the standard AdK system in implicit solvent."""

    @staticmethod
    @pytest.fixture()
    def universe():
        return mda.Universe(PSF, DCD)

    @staticmethod
    @pytest.fixture()
    def atom(universe):
        # Leu67:CG
        return universe.atoms[1000]

    def test_attributes_names(self, atom):
        a = atom
        assert a.name == 'CG'
        assert a.resname == 'LEU'

    def test_setting_attribute_name(self, atom):
        atom.name = 'AA'
        assert atom.name == 'AA'

    def test_setting_attribute_type(self, atom):
        atom.type = 'Z'
        assert atom.type == 'Z'

    def test_setting_attribute_mass(self, atom):
        atom.mass = 13
        assert atom.mass == 13

    def test_setting_attributes_charge(self, atom):
        atom.charge = 6
        assert atom.charge == 6

    def test_attributes_positions(self, atom):
        known_pos = np.array([3.94543672, -12.4060812, -7.26820087], dtype=np.float32)
        a = atom
        # new position property (mutable)
        assert_almost_equal(a.position, known_pos)
        pos = a.position + 3.14
        a.position = pos
        assert_almost_equal(a.position, pos)

    def test_atom_selection(self, universe, atom):
        asel = universe.select_atoms('atom 4AKE 67 CG').atoms[0]
        assert atom == asel

    def test_hierarchy(self, universe, atom):
        u = universe
        a = atom
        assert a.segment == u.select_atoms('segid 4AKE').segments[0]
        assert a.residue == u.residues[66]

    def test_bad_add(self, atom):
        with pytest.raises(TypeError):
            atom + 1

    def test_add_AG(self, universe, atom):
        ag = universe.atoms[:2]

        ag2 = atom + ag

        for at in [atom, ag[0], ag[1]]:
            assert at in ag2

    def test_no_velo(self, atom):
        with pytest.raises(NoDataError):
            atom.velocity

    def test_bonded_atoms(self, universe):
        at = universe.atoms[0]
        ref = [b.partner(at) for b in at.bonds]
        assert ref == list(at.bonded_atoms)

    def test_undefined_occupancy(self, universe):
        with pytest.raises(AttributeError):
            universe.atoms[0].occupancy

    @pytest.mark.parametrize("ix", (1, -1))
    def test_atom_pickle(self, universe, ix):
        atm_out = universe.atoms[ix]
        atm_in = pickle.loads(pickle.dumps(atm_out))
        assert atm_in == atm_out


class TestAtomNoForceNoVel(object):
    @staticmethod
    @pytest.fixture()
    def a():
        u = mda.Universe(XYZ_mini)
        return u.atoms[0]

    def test_velocity_fail(self, a):
        with pytest.raises(NoDataError):
            getattr(a, 'velocity')

    def test_force_fail(self, a):
        with pytest.raises(NoDataError):
            getattr(a, 'force')

    def test_velocity_set_fail(self, a):
        with pytest.raises(NoDataError):
            setattr(a, 'velocity', [1.0, 1.0, 1.0])

    def test_force_set_fail(self, a):
        with pytest.raises(NoDataError):
            setattr(a, 'force', [1.0, 1.0, 1.0])
