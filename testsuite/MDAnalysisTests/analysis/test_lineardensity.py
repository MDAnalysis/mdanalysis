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
import numpy as np

from MDAnalysisTests.datafiles import waterPSF, waterDCD
from MDAnalysis.analysis.lineardensity import LinearDensity
from numpy.testing import assert_almost_equal



def test_grouping_atoms():
    """For testing the case of grouping='atoms'"""
    universe = mda.Universe(waterPSF, waterDCD)
    sel_string = 'all'
    selection = universe.select_atoms(sel_string)

    ld = LinearDensity(selection, grouping="atoms", binsize=5).run()

    expected_masses = np.array([15.9994, 1.008, 1.008, 15.9994, 1.008, 1.008,
                                15.9994, 1.008, 1.008, 15.9994, 1.008, 1.008,
                                15.9994, 1.008, 1.008])

    expected_charges = np.array([-0.834, 0.417, 0.417, -0.834, 0.417,
                                 0.417, -0.834, 0.417, 0.417, -0.834,
                                 0.417, 0.417, -0.834, 0.417, 0.417])

    xpos = np.array([0., 0., 0., 0.0072334, 0.00473299, 0.,
                          0., 0., 0., 0.])

    assert_almost_equal(ld.masses, expected_masses)
    assert_almost_equal(ld.charges, expected_charges)
    assert_almost_equal(ld.results['x']['pos'], xpos)

def test_grouping_residues():
    """For testing the case of grouping='residues'"""
    universe = mda.Universe(waterPSF, waterDCD)
    sel_string = 'all'
    selection = universe.select_atoms(sel_string)
    ld = LinearDensity(selection, grouping="residues", binsize=5).run()

    expected_masses = np.array([18.0154, 18.0154, 18.0154, 18.0154, 18.0154])
    expected_charges = np.array([0, 0, 0, 0, 0])

    assert_almost_equal(ld.masses, expected_masses)
    assert_almost_equal(ld.charges, expected_charges)


def test_grouping_segments():
    """For testing the case of grouping='segments'"""
    universe = mda.Universe(waterPSF, waterDCD)
    sel_string = 'all'
    selection = universe.select_atoms(sel_string)
    ld = LinearDensity(selection, grouping="segments", binsize=5).run()

    expected_masses = np.array([90.0770])
    expected_charges = np.array([0])
    assert_almost_equal(ld.masses, expected_masses)
    assert_almost_equal(ld.charges, expected_charges)

def test_grouping_fragments():
    """For testing the case of grouping='fragments'"""
    universe = mda.Universe(waterPSF, waterDCD)
    sel_string = 'all'
    selection = universe.select_atoms(sel_string)
    ld = LinearDensity(selection, grouping="fragments", binsize=5).run()

    expected_masses = np.array([18.0154, 18.0154, 18.0154, 18.0154, 18.0154])
    expected_charges = np.array([0, 0, 0, 0, 0])

    assert_almost_equal(ld.masses, expected_masses)
    assert_almost_equal(ld.charges, expected_charges)

