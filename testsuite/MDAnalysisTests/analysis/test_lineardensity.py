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
import pytest

from MDAnalysisTests.datafiles import waterPSF, waterDCD
from MDAnalysis.analysis.lineardensity import LinearDensity
from numpy.testing import assert_allclose


def test_invalid_grouping():
    """Invalid groupings raise AttributeError"""
    universe = mda.Universe(waterPSF, waterDCD)
    sel_string = 'all'
    selection = universe.select_atoms(sel_string)
    with pytest.raises(AttributeError):
        # centroid is attribute of AtomGroup, but not valid here
        ld = LinearDensity(selection, grouping="centroid", binsize=5)
        ld.run()


# test data for grouping='atoms'
expected_masses_atoms = np.array([15.9994, 1.008, 1.008, 15.9994, 1.008, 1.008,
                                  15.9994, 1.008, 1.008, 15.9994, 1.008, 1.008,
                                  15.9994, 1.008, 1.008])
expected_charges_atoms = np.array([-0.834, 0.417, 0.417, -0.834, 0.417,
                                   0.417, -0.834, 0.417, 0.417, -0.834,
                                   0.417, 0.417, -0.834, 0.417, 0.417])
expected_xpos_atoms = np.array([0., 0., 0., 0.0072334, 0.00473299, 0.,
                                0., 0., 0., 0.])
expected_xchar_atoms = np.array([0., 0., 0., 2.2158751e-05, -2.2158751e-05,
                                 0., 0., 0., 0., 0.])

# test data for grouping='residues'
expected_masses_residues = np.array([18.0154, 18.0154, 18.0154, 18.0154,
                                     18.0154])
expected_charges_residues = np.array([0, 0, 0, 0, 0])
expected_xpos_residues = np.array([0., 0., 0., 0.00717983, 0.00478656,
                                   0., 0., 0., 0., 0.])
expected_xchar_residues = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])

# test data for grouping='segments'
expected_masses_segments = np.array([90.0770])
expected_charges_segments = np.array([0])
expected_xpos_segments = np.array([0., 0., 0., 0.01196639, 0.,
                                   0., 0., 0., 0., 0.])
expected_xchar_segments = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])

# test data for grouping='fragments'
expected_masses_fragments = np.array([18.0154, 18.0154, 18.0154, 18.0154,
                                      18.0154])
expected_charges_fragments = np.array([0, 0, 0, 0, 0])
expected_xpos_fragments = np.array([0., 0., 0., 0.00717983, 0.00478656,
                                   0., 0., 0., 0., 0.])
expected_xchar_fragments = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])


@pytest.mark.parametrize("grouping, expected_masses, expected_charges,\
                         expected_xpos, expected_xchar", [
                         ("atoms",
                          expected_masses_atoms,
                          expected_charges_atoms,
                          expected_xpos_atoms,
                          expected_xchar_atoms),
                         ("residues",
                          expected_masses_residues,
                          expected_charges_residues,
                          expected_xpos_residues,
                          expected_xchar_residues),
                         ("segments",
                          expected_masses_segments,
                          expected_charges_segments,
                          expected_xpos_segments,
                          expected_xchar_segments),
                         ("fragments",
                          expected_masses_fragments,
                          expected_charges_fragments,
                          expected_xpos_fragments,
                          expected_xchar_fragments)
                         ])
def test_lineardensity(grouping, expected_masses, expected_charges,
                       expected_xpos, expected_xchar):
    universe = mda.Universe(waterPSF, waterDCD)
    sel_string = 'all'
    selection = universe.select_atoms(sel_string)
    ld = LinearDensity(selection, grouping, binsize=5).run()
    assert_allclose(ld.masses, expected_masses)
    assert_allclose(ld.charges, expected_charges)
    # rtol changed here due to floating point imprecision
    assert_allclose(ld.results['x']['pos'], expected_xpos, rtol=1e-06)
    assert_allclose(ld.results['x']['char'], expected_xchar)
