# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2024 The MDAnalysis Development Team and contributors
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
import numpy as np

import MDAnalysis
from MDAnalysis.transformations import unwrap
from MDAnalysis.guesser.tables import vdwradii
from MDAnalysis.exceptions import NoDataError
from MDAnalysisTests.datafiles import TPR, PDB, XTC, PDB_sasa, PDB_rsasa

from MDAnalysis.analysis.sasa import _sasa, SASA, RSASA


@pytest.fixture
def protein():
    """Returns AdK protein AtomGroup"""
    u = MDAnalysis.Universe(TPR, PDB)
    u.trajectory.add_transformations(unwrap(u.atoms))
    if not hasattr(u, "radii"):
        u.add_TopologyAttr("radii")
        u.atoms.radii = [vdwradii.get(e, 2.0) for e in u.atoms.elements]
    protein = u.select_atoms("protein")
    return protein


@pytest.fixture
def residue(protein):
    """Returns single residue of AdK protein AtomGroup"""
    return protein.select_atoms("resid 1")


@pytest.fixture
def trajectory():
    """Returns AdK protein AtomGroup"""
    u = MDAnalysis.Universe(TPR, XTC)
    u.trajectory.add_transformations(unwrap(u.atoms))
    if not hasattr(u, "radii"):
        u.add_TopologyAttr("radii")
        u.atoms.radii = [vdwradii.get(e, 2.0) for e in u.atoms.elements]
    protein = u.select_atoms("protein")
    return protein


@pytest.fixture
def co():
    """Return Carbon Monoxide (CO) AtomGroup"""
    u = MDAnalysis.Universe.empty(2, trajectory=True)
    u.add_TopologyAttr("names", ["O", "C"])
    u.add_TopologyAttr("radii", [1.52, 1.70])
    u.atoms.positions = [[ 0.5285, 0., 0.],
                         [-0.5285, 0., 0.]]
    return u.atoms


@pytest.fixture
def co_area_exact(co):
    """Analytical solution for two particle overlap"""
    # See: https://en.wikipedia.org/wiki/Spherical_cap
    r1, r2 = co.atoms[0].radius, co.atoms[1].radius
    d = np.linalg.norm(co.atoms[0].position - co.atoms[1].position)
    h1 = 0.5 * (r2 - r1 + d) * (r2 + r1 - d) / d
    h2 = 0.5 * (r1 - r2 + d) * (r1 + r2 - d) / d
    a1 = (4 * np.pi * r1 ** 2) - (2 * np.pi * r1 * h1)
    a2 = (4 * np.pi * r2 ** 2) - (2 * np.pi * r2 * h2)
    return [a1, a2]


def test_analytical(co, co_area_exact):
    """Test sasa function against analytical result"""
    co_area = _sasa(co.positions, co.radii, probe_radius=0, n_dots=1024)
    assert np.allclose(co_area, co_area_exact, 0.01)


class TestSASA(object):
    # Parameters used for testing
    PROBE_RADIUS = 1.4
    N_DOTS = 256
    RADII_DICT = {"H": 1.1, "C": 1.7, "N": 1.55, "O": 1.52, "S": 1.8}

    def test_sasa_co(self, co, co_area_exact):
        """Test for CO against analytical result"""
        r = SASA(co, probe_radius=0., n_dots=1_024).run()
        assert np.allclose(r.results.area.mean(axis=0),
                           co_area_exact, 0.01)

    def test_sasa(self, protein):
        """Test against reference values"""
        r = SASA(protein, self.PROBE_RADIUS, self.N_DOTS).run()
        assert r.results.area.shape[0] == 1
        assert r.results.area.shape[1] == protein.n_atoms
        actual = r.results.area.mean(axis=0)
        desired = np.load(PDB_sasa)
        # Difference should be less than 1% or 0.1 Å²
        np.testing.assert_allclose(actual, desired, 0.01, 0.1)

    def test_sasa_traj(self, trajectory):
        """Test for trajectory"""
        # Three residues is stressful enough
        ag = trajectory.select_atoms("resid 1-3")
        r = SASA(ag).run()
        assert r.results.area.shape[0] == len(ag.universe.trajectory)
        assert r.results.area.shape[1] == ag.n_atoms

    @pytest.mark.parametrize("n_dots", [16, 128])
    def test_sasa_dots(self, residue, n_dots):
        """Test for different number of dots"""
        r = SASA(residue, n_dots=n_dots).run()
        assert r.results.area.sum() != 0
        assert r.n_dots == n_dots

    def test_sasa_n_dots_error(self, residue):
        """Should raise error if number of dots less than 1"""
        with pytest.raises(ValueError):
            SASA(residue, n_dots=0)

    @pytest.mark.parametrize("probe_radius", [0., 1.4])
    def test_sasa_probe_radius(self, residue, probe_radius):
        """Test for different probe radius"""
        r = SASA(residue, probe_radius=probe_radius).run()
        assert r.results.area.sum() != 0
        assert r.probe_radius == probe_radius

    def test_sasa_probe_radius_error(self, residue):
        """Should raise error if probe radius less than 0"""
        with pytest.raises(ValueError):
            SASA(residue, probe_radius=-1.0)

    def test_sasa_updating_group_error(self, residue):
        """Should raise error if UpdatingAtomGroup"""
        updating = residue.select_atoms("prop x < 1", updating=True)
        with pytest.raises(TypeError):
            SASA(updating)

    def test_sasa_no_radii_error(self):
        """Should raise error if no radii attribute"""
        u = MDAnalysis.Universe.empty(2, trajectory=True)
        assert not hasattr(u, "radii")
        with pytest.raises(NoDataError):
            SASA(u.atoms)
        # Should also raise error for radii is less than zero
        u.add_TopologyAttr("radii")
        u.atoms.radii = 0
        with pytest.raises(ValueError):
            SASA(u.atoms)


class TestRSASA(object):
    # Parameters used for testing
    PROBE_RADIUS = 1.4
    N_DOTS = 256
    RADII_DICT = {"H": 1.1, "C": 1.7, "N": 1.55, "O": 1.52, "S": 1.8}

    def test_rsasa_residue(self, residue):
        """Test for a single residue"""
        r = RSASA(residue).run()
        assert r.results.relative_area == 1.0

    def test_rsasa(self, protein):
        """Test against reference values"""
        r = RSASA(protein, None, self.PROBE_RADIUS, self.N_DOTS).run()
        assert r.results.relative_area.shape[0] == 1
        assert r.results.relative_area.shape[1] == protein.n_residues
        actual = r.results.relative_area.mean(axis=0)
        desired = np.load(PDB_rsasa)
        # Difference should be less than 5% rel and 1% abs
        np.testing.assert_allclose(actual, desired, 0.05, 0.01)

    def test_rsasa_traj(self, trajectory):
        """Test for trajectory"""
        # Three residues is stressful enough
        ag = trajectory.select_atoms("resid 1-3")
        r = RSASA(ag).run()
        assert r.results.relative_area.shape[0] == len(ag.universe.trajectory)
        assert r.results.relative_area.shape[1] == ag.n_residues

    @pytest.mark.parametrize("n_dots", [16, 128])
    def test_rsasa_dots(self, residue, n_dots):
        """Test for different number of dots"""
        r = RSASA(residue, n_dots=n_dots).run()
        assert r.results.relative_area.sum() == 1
        assert r.n_dots == n_dots

    def test_rsasa_n_dots_error(self, residue):
        """Should raise error if number of dots less than 1"""
        with pytest.raises(ValueError):
            RSASA(residue, n_dots=0)

    @pytest.mark.parametrize("probe_radius", [0., 1.4])
    def test_rsasa_probe_radius(self, residue, probe_radius):
        """Test for different probe radius"""
        r = RSASA(residue, probe_radius=probe_radius).run()
        assert r.results.relative_area.sum() == 1
        assert r.probe_radius == probe_radius

    def test_rsasa_probe_radius_error(self, residue):
        """Should raise error if probe radius less than 0"""
        with pytest.raises(ValueError):
            RSASA(residue, probe_radius=-1.0)

    def test_rsasa_updating_group_error(self, residue):
        """Should raise error if UpdatingAtomGroup"""
        updating = residue.select_atoms("prop x < 1", updating=True)
        with pytest.raises(TypeError):
            RSASA(updating)

    @pytest.mark.parametrize("subsele", [None, "all", "not backbone"])
    def test_rsasa_subsele(self, residue, subsele):
        """Should raise error if UpdatingAtomGroup"""
        r = RSASA(residue, subsele=subsele).run()

    def test_rsasa_no_radii_error(self):
        """Should raise error if no radii attribute"""
        u = MDAnalysis.Universe.empty(2, trajectory=True)
        u.add_TopologyAttr("names", ["O", "C"])
        u.atoms.positions = [[ 0.5285, 0., 0.],
                             [-0.5285, 0., 0.]]
        u.guess_TopologyAttrs(to_guess=["bonds"])
        assert not hasattr(u, "radii")
        with pytest.raises(NoDataError):
            RSASA(u.atoms)
        # Should also raise error for radii less than zero
        u.add_TopologyAttr("radii")
        u.atoms.radii = 0
        with pytest.raises(ValueError):
            RSASA(u.atoms)

    def test_rsasa_no_bonds_error(self):
        """Should raise error if no bonds attribute"""
        u = MDAnalysis.Universe.empty(2, trajectory=True)
        u.add_TopologyAttr("names", ["O", "C"])
        u.add_TopologyAttr("radii", [1.52, 1.70])
        u.atoms.positions = [[ 0.5285, 0., 0.],
                             [-0.5285, 0., 0.]]
        assert not hasattr(u, "bonds")
        with pytest.raises(NoDataError):
            RSASA(u.atoms)
