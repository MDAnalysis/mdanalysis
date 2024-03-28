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
from MDAnalysis.core._get_readers import get_reader_for
from MDAnalysisTests.util import no_deprecated_call


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
expected_xmass_atoms = np.array([0., 0., 0., 0.00723323, 0.00473288, 0.,
                                0., 0., 0., 0.])
expected_xcharge_atoms = np.array([0., 0., 0., 2.21582311e-05,
                                   -2.21582311e-05, 0., 0., 0., 0., 0.])

# test data for grouping='residues'
expected_masses_residues = np.array([18.0154, 18.0154, 18.0154, 18.0154,
                                     18.0154])
expected_charges_residues = np.array([0, 0, 0, 0, 0])
expected_xmass_residues = np.array([0., 0., 0., 0.00717967, 0.00478644,
                                   0., 0., 0., 0., 0.])
expected_xcharge_residues = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])

# test data for grouping='segments'
expected_masses_segments = np.array([90.0770])
expected_charges_segments = np.array([0])
expected_xmass_segments = np.array([0., 0., 0., 0.01196611, 0.,
                                   0., 0., 0., 0., 0.])
expected_xcharge_segments = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])

# test data for grouping='fragments'
expected_masses_fragments = np.array([18.0154, 18.0154, 18.0154, 18.0154,
                                      18.0154])
expected_charges_fragments = np.array([0, 0, 0, 0, 0])
expected_xmass_fragments = np.array([0., 0., 0., 0.00717967, 0.00478644,
                                     0., 0., 0., 0., 0.])
expected_xcharge_fragments = np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])


@pytest.mark.parametrize("grouping, expected_masses, expected_charges, expected_xmass, expected_xcharge", [
    ("atoms", expected_masses_atoms, expected_charges_atoms,
     expected_xmass_atoms, expected_xcharge_atoms),
    ("residues", expected_masses_residues, expected_charges_residues,
     expected_xmass_residues, expected_xcharge_residues),
    ("segments", expected_masses_segments, expected_charges_segments,
     expected_xmass_segments, expected_xcharge_segments),
    ("fragments", expected_masses_fragments, expected_charges_fragments,
     expected_xmass_fragments, expected_xcharge_fragments)
])
def test_lineardensity(grouping, expected_masses, expected_charges,
                       expected_xmass, expected_xcharge):
    universe = mda.Universe(waterPSF, waterDCD)
    sel_string = 'all'
    selection = universe.select_atoms(sel_string)
    ld = LinearDensity(selection, grouping, binsize=5).run()
    assert_allclose(ld.masses, expected_masses)
    assert_allclose(ld.charges, expected_charges)
    # rtol changed here due to floating point imprecision
    assert_allclose(ld.results.x.mass_density, expected_xmass, rtol=1e-06)
    assert_allclose(ld.results.x.charge_density, expected_xcharge)


@pytest.fixture(scope="module")
def testing_Universe():
    """Generate a universe for testing whether LinearDensity works with
       updating atom groups. Also used for parallel analysis test."""
    n_atoms = 3
    u = mda.Universe.empty(n_atoms=n_atoms,
                           n_residues=n_atoms,
                           n_segments=n_atoms,
                           atom_resindex=np.arange(n_atoms),
                           residue_segindex=np.arange(n_atoms))

    for attr in ["charges", "masses"]:
        u.add_TopologyAttr(attr, values=np.ones(n_atoms))

    coords = np.array([
                      [[1., 1., 1.], [1., 2., 1.], [2., 1., 1.]],
                      [[1., 1., 2.], [1., 2., 1.], [2., 1., 1.]],
                      [[1., 1., 3.], [1., 2., 1.], [2., 1., 1.]],
                      [[1., 1., 4.], [1., 2., 1.], [2., 1., 1.]],
                      [[1., 1., 5.], [1., 2., 1.], [2., 1., 1.]]
                      ])

    u.trajectory = get_reader_for(coords)(coords, order='fac', n_atoms=n_atoms)

    for ts in u.trajectory:
        ts.dimensions = np.array([2, 2, 6, 90, 90, 90])

    return u


def test_updating_atomgroup(testing_Universe):
    expected_z_pos = np.array([0., 0.91329641, 0.08302695, 0., 0., 0.])
    u = testing_Universe
    selection = u.select_atoms("prop z < 3", updating=True)
    ld = LinearDensity(selection, binsize=1).run()
    assert_allclose(ld.results.z.mass_density, expected_z_pos)
    # Test whether histogram bins are saved correctly.
    expected_bin_edges = np.arange(0, 7)
    assert_allclose(ld.results.x.hist_bin_edges, expected_bin_edges)


testdict = {"pos": "mass_density",
            "pos_std": "mass_density_stddev",
            "char": "charge_density",
            "char_std": "charge_density_stddev"}


# TODO: Remove in 3.0.0
def test_old_name_deprecations():
    universe = mda.Universe(waterPSF, waterDCD)
    sel_string = 'all'
    selection = universe.select_atoms(sel_string)
    ld = LinearDensity(selection, binsize=5).run()
    with pytest.warns(DeprecationWarning):
        assert_allclose(ld.results.x.pos, ld.results.x.mass_density)
        assert_allclose(ld.results.x.pos_std, ld.results.x.mass_density_stddev)
        assert_allclose(ld.results.x.char, ld.results.x.charge_density)
        assert_allclose(ld.results.x.char_std,
                        ld.results.x.charge_density_stddev)
        for key in testdict.keys():
            assert_allclose(ld.results["x"][key],
                            ld.results["x"][testdict[key]])

    # Check that no DeprecationWarning is raised with new attributes
    with no_deprecated_call():
        ld.results.x.mass_density
        ld.results.x.mass_density_stddev
        ld.results.x.charge_density
        ld.results.x.charge_density_stddev


# TODO: deprecated, remove in 3.0.0
def test_parallel_analysis(testing_Universe):
    """tests _add_other_result() method. Runs LinearDensity for all atoms of
    a universe and for two subsets, then adds the results of the two subsets
    and checks the results are the same."""
    u = testing_Universe
    selection1 = u.select_atoms("prop x < 1.1")
    selection2 = u.select_atoms("prop x >= 1.1")
    selection_whole = u.select_atoms("all")
    ld1 = LinearDensity(selection1, binsize=1).run()
    ld2 = LinearDensity(selection2, binsize=1).run()
    ld_whole = LinearDensity(selection_whole, binsize=1).run()
    with pytest.warns(DeprecationWarning,
                      match="`_add_other_results` is deprecated!"):
        ld1._add_other_results(ld2)
    assert_allclose(ld1.results.z.mass_density, ld_whole.results.z.mass_density)
    assert_allclose(ld1.results.x.mass_density, ld_whole.results.x.mass_density)
