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

import MDAnalysis
from MDAnalysis.analysis.sasa import SASA, RSASA
from MDAnalysis.transformations import unwrap
from MDAnalysisTests.datafiles import TPR, XTC, GRO

import logging
import numpy as np
from numpy.testing import assert_almost_equal


@pytest.fixture
def u():
    """Returns AdK protein AtomGroup"""
    u = MDAnalysis.Universe(TPR, GRO)
    u.trajectory.add_transformations(unwrap(u.atoms))
    return u


@pytest.fixture
def protein(u):
    """Returns AdK protein AtomGroup"""
    return u.select_atoms("protein")


@pytest.fixture
def residue(u):
    """Returns AdK 1st residue AtomGroup"""
    return u.select_atoms("resid 1")


@pytest.fixture
def co():
    """Return Carbon Monoxide (CO) AtomGroup"""
    u = MDAnalysis.Universe.empty(2, trajectory=True)
    u.add_TopologyAttr("elements", ["O", "C"])
    u.atoms.positions = [[0.5285, 0., 0.], [-0.5285, 0., 0.]]
    return u.atoms


@pytest.fixture
def traj():
    """Returns AdK protein trajectory AtomGroup"""
    u = MDAnalysis.Universe(TPR, XTC)
    u.trajectory.add_transformations(unwrap(u.atoms))
    return u.select_atoms("protein")


@pytest.fixture
def radii_dict():
    """Returns custom radii dictionary"""
    # This is to ensure same result if the built-in vdw table ever updates
    return {"H": 1.1, "C": 1.7, "N": 1.55, "O": 1.52}


@pytest.fixture
def analytical_co(co, radii_dict):
    """Analytical solution for two particle overlap"""
    r1 = radii_dict[co.atoms[0].element]
    r2 = radii_dict[co.atoms[1].element]
    d = np.linalg.norm(co.atoms[0].position - co.atoms[1].position)
    h1 = 0.5 * (r2 - r1 + d) * (r2 + r1 - d) / d
    h2 = 0.5 * (r1 - r2 + d) * (r1 + r2 - d) / d
    a1 = (4 * np.pi * r1 ** 2) - (2 * np.pi * r1 * h1)
    a2 = (4 * np.pi * r2 ** 2) - (2 * np.pi * r2 * h2)
    return [a1, a2]


class TestSASA(object):

    def test_sasa_analytical(self, co, analytical_co, radii_dict):
        """Compare to analytical result"""
        R = SASA(co, n_dots=1_024, probe_radius=0, radii_dict=radii_dict).run()
        assert_almost_equal(R.results.area, analytical_co, decimal=1)

    def test_sasa_residue(self, residue, radii_dict):
        """Test for a single residue"""
        R = SASA(residue, radii_dict=radii_dict).run()

    @pytest.mark.parametrize('n_dots', [1, 1000])
    def test_sasa_dots(self, residue, n_dots):
        """Test for different number of dots"""
        R = SASA(residue, n_dots=n_dots).run()

    @pytest.mark.parametrize('probe_radius', [0., 2.])
    def test_sasa_probe_radius(self, residue, probe_radius):
        """Test for different probe radius"""
        R = SASA(residue, probe_radius=probe_radius).run()

    def test_sasa_single_frame(self, protein, radii_dict):
        """Test calculation on a single frame"""
        R = SASA(protein, radii_dict=radii_dict).run()

    def test_sasa(self, traj, radii_dict):
        """Test calculation on a trajectory"""
        # 2 frames are stressful enough
        R = SASA(traj, radii_dict=radii_dict).run(step=5)

    def test_sasa_updating_atom_group(self, protein):
        """Error if UpdatingAtomGroup"""
        updating = protein.select_atoms("prop x < 1", updating=True)
        with pytest.raises(TypeError):
            R = SASA(updating)

    def test_sasa_no_element_warning(self, residue, caplog):
        """Check if warning is logged if unknown element"""
        unknown = residue.copy()
        unknown.atoms[0].element = "X"
        R = SASA(unknown)
        assert any(
            record.levelno == logging.WARNING for record in caplog.records)

    def test_sasa_wrong_probe_radius(self, protein):
        """Error if probe radius less than 0"""
        with pytest.raises(ValueError):
            R = SASA(protein, probe_radius=-1.0)

    def test_sasa_wrong_n_dots(self, protein):
        """Error if number of dots less than 1"""
        with pytest.raises(ValueError):
            R = SASA(protein, n_dots=0)

    def test_sasa_no_elements_attribute(self, protein):
        """Error if AtomGroup has no elements attribute"""
        no_elements = protein.copy()
        no_elements.universe.del_TopologyAttr("elements")
        with pytest.raises(ValueError):
            R = SASA(no_elements)


class TestSASA_r(object):

    def test_rsasa_analytical(self, protein, co, analytical_co, radii_dict):
        """Compare to analytical result"""
        # HACK: Initiate class but call private method on CO data instead
        RS = RSASA(protein, n_dots=1_024,
                   probe_radius=0, radii_dict=radii_dict)
        area = RS._get_sasa(co)
        assert_almost_equal(area, analytical_co, decimal=1)

    def test_rsasa_residue(self, residue, radii_dict):
        """Test for a single residue"""
        RS = RSASA(residue.atoms, radii_dict=radii_dict).run()
        assert RS.results.relative_area == 1.0

    @pytest.mark.parametrize('n_dots', [1, 1000])
    def test_rsasa_dots(self, residue, n_dots):
        """Test for different number of dots"""
        RS = RSASA(residue, n_dots=n_dots).run()
        assert all(0 <= a <= 1 for a in RS.results.relative_area)

    @pytest.mark.parametrize('probe_radius', [0., 2.])
    def test_rsasa_probe_radius(self, residue, probe_radius):
        """Test for different probe radius"""
        RS = RSASA(residue, probe_radius=probe_radius).run()
        assert all(0 <= a <= 1 for a in RS.results.relative_area)

    def test_rsasa_subsele(self, protein):
        """Test for a sub-selection"""
        RS = RSASA(protein, subsele="not backbone").run()
        assert all(0 <= a <= 1 for a in RS.results.relative_area)

    def test_rsasa_single_frame(self, protein, radii_dict):
        """Test calculation on a single frame"""
        RS = RSASA(protein, radii_dict=radii_dict).run()
        assert all(0 <= a <= 1 for a in RS.results.relative_area)

    def test_rsasa(self, traj, radii_dict):
        """Test calculation on a trajectory"""
        # 2 frames are stressful enough
        RS = RSASA(traj, radii_dict=radii_dict).run(step=5)
        assert all(0 <= a <= 1 for a in RS.results.relative_area)

    def test_rsasa_updating_atom_group(self, protein):
        """Error if UpdatingAtomGroup"""
        updating = protein.select_atoms("prop x < 1", updating=True)
        with pytest.raises(TypeError):
            RS = RSASA(updating).run()

    def test_rsasa_no_element_warning(self, residue, caplog):
        """Check if warning is logged if unknown element"""
        unknown = residue.copy()
        unknown.atoms[0].element = "X"
        RS = RSASA(unknown)
        assert any(
            record.levelno == logging.WARNING for record in caplog.records)

    def test_rsasa_wrong_probe_radius(self, protein):
        """Error if probe radius less than 0"""
        with pytest.raises(ValueError):
            RS = RSASA(protein, probe_radius=-1.0).run()

    def test_rsasa_wrong_n_dots(self, protein):
        """Error if number of dots less than 1"""
        with pytest.raises(ValueError):
            RS = RSASA(protein, n_dots=0).run()

    def test_rsasa_no_elements_attribute(self, protein):
        """Error if AtomGroup has no elements attribute"""
        no_elements = protein.copy()
        no_elements.universe.del_TopologyAttr("elements")
        with pytest.raises(ValueError):
            RS = RSASA(no_elements).run()

    def test_rsasa_no_bonds_attribute(self, co):
        """Error if AtomGroup has no bonds attribute"""
        with pytest.raises(ValueError):
            RS = RSASA(co).run()
