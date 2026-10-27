# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
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
import warnings

from MDAnalysisTests.datafiles import waterPSF, waterDCD
from MDAnalysis.analysis.lineardensity import LinearDensity
from numpy.testing import assert_allclose
from MDAnalysis.core._get_readers import get_reader_for
from MDAnalysisTests.util import no_deprecated_call


def test_invalid_grouping(client_LinearDensity):
    """Invalid groupings raise AttributeError"""
    universe = mda.Universe(waterPSF, waterDCD)
    sel_string = "all"
    selection = universe.select_atoms(sel_string)
    with pytest.raises(AttributeError):
        # centroid is attribute of AtomGroup, but not valid here
        ld = LinearDensity(selection, grouping="centroid", binsize=5)
        ld.run(**client_LinearDensity)


# test data for grouping='atoms'
expected_masses_atoms = np.array(
    [
        15.9994,
        1.008,
        1.008,
        15.9994,
        1.008,
        1.008,
        15.9994,
        1.008,
        1.008,
        15.9994,
        1.008,
        1.008,
        15.9994,
        1.008,
        1.008,
    ]
)
expected_charges_atoms = np.array(
    [
        -0.834,
        0.417,
        0.417,
        -0.834,
        0.417,
        0.417,
        -0.834,
        0.417,
        0.417,
        -0.834,
        0.417,
        0.417,
        -0.834,
        0.417,
        0.417,
    ]
)
expected_xmass_atoms = np.array(
    [0.0, 0.0, 0.0, 0.00723323, 0.00473288, 0.0, 0.0, 0.0, 0.0, 0.0]
)
expected_xcharge_atoms = np.array(
    [0.0, 0.0, 0.0, 2.21582311e-05, -2.21582311e-05, 0.0, 0.0, 0.0, 0.0, 0.0]
)

# test data for grouping='residues'
expected_masses_residues = np.array(
    [18.0154, 18.0154, 18.0154, 18.0154, 18.0154]
)
expected_charges_residues = np.array([0, 0, 0, 0, 0])
expected_xmass_residues = np.array(
    [0.0, 0.0, 0.0, 0.00717967, 0.00478644, 0.0, 0.0, 0.0, 0.0, 0.0]
)
expected_xcharge_residues = np.array(
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
)

# test data for grouping='segments'
expected_masses_segments = np.array([90.0770])
expected_charges_segments = np.array([0])
expected_xmass_segments = np.array(
    [0.0, 0.0, 0.0, 0.01196611, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
)
expected_xcharge_segments = np.array(
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
)

# test data for grouping='fragments'
expected_masses_fragments = np.array(
    [18.0154, 18.0154, 18.0154, 18.0154, 18.0154]
)
expected_charges_fragments = np.array([0, 0, 0, 0, 0])
expected_xmass_fragments = np.array(
    [0.0, 0.0, 0.0, 0.00717967, 0.00478644, 0.0, 0.0, 0.0, 0.0, 0.0]
)
expected_xcharge_fragments = np.array(
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
)


@pytest.mark.parametrize(
    "grouping, expected_masses, expected_charges, expected_xmass, expected_xcharge",
    [
        (
            "atoms",
            expected_masses_atoms,
            expected_charges_atoms,
            expected_xmass_atoms,
            expected_xcharge_atoms,
        ),
        (
            "residues",
            expected_masses_residues,
            expected_charges_residues,
            expected_xmass_residues,
            expected_xcharge_residues,
        ),
        (
            "segments",
            expected_masses_segments,
            expected_charges_segments,
            expected_xmass_segments,
            expected_xcharge_segments,
        ),
        (
            "fragments",
            expected_masses_fragments,
            expected_charges_fragments,
            expected_xmass_fragments,
            expected_xcharge_fragments,
        ),
    ],
)
def test_lineardensity(
    grouping,
    expected_masses,
    expected_charges,
    expected_xmass,
    expected_xcharge,
    client_LinearDensity,
):
    universe = mda.Universe(waterPSF, waterDCD)
    sel_string = "all"
    selection = universe.select_atoms(sel_string)
    ld = LinearDensity(selection, grouping, binsize=5).run(
        **client_LinearDensity
    )
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
    u = mda.Universe.empty(
        n_atoms=n_atoms,
        n_residues=n_atoms,
        n_segments=n_atoms,
        atom_resindex=np.arange(n_atoms),
        residue_segindex=np.arange(n_atoms),
    )

    for attr in ["charges", "masses"]:
        u.add_TopologyAttr(attr, values=np.ones(n_atoms))

    coords = np.array(
        [
            [[1.0, 1.0, 1.0], [1.0, 2.0, 1.0], [2.0, 1.0, 1.0]],
            [[1.0, 1.0, 2.0], [1.0, 2.0, 1.0], [2.0, 1.0, 1.0]],
            [[1.0, 1.0, 3.0], [1.0, 2.0, 1.0], [2.0, 1.0, 1.0]],
            [[1.0, 1.0, 4.0], [1.0, 2.0, 1.0], [2.0, 1.0, 1.0]],
            [[1.0, 1.0, 5.0], [1.0, 2.0, 1.0], [2.0, 1.0, 1.0]],
        ]
    )

    u.trajectory = get_reader_for(coords)(coords, order="fac", n_atoms=n_atoms)

    for ts in u.trajectory:
        ts.dimensions = np.array([2, 2, 6, 90, 90, 90])

    return u


def test_updating_atomgroup(testing_Universe, client_LinearDensity):
    expected_z_pos = np.array([0.0, 0.91329641, 0.08302695, 0.0, 0.0, 0.0])
    u = testing_Universe
    selection = u.select_atoms("prop z < 3", updating=True)
    ld = LinearDensity(selection, binsize=1).run(**client_LinearDensity)
    assert_allclose(ld.results.z.mass_density, expected_z_pos)
    # Test whether histogram bins are saved correctly.
    expected_bin_edges = np.arange(0, 7)
    assert_allclose(ld.results.x.hist_bin_edges, expected_bin_edges)


testdict = {
    "pos": "mass_density",
    "pos_std": "mass_density_stddev",
    "char": "charge_density",
    "char_std": "charge_density_stddev",
}


# TODO: Remove in 3.0.0
def test_old_name_deprecations():
    universe = mda.Universe(waterPSF, waterDCD)
    sel_string = "all"
    selection = universe.select_atoms(sel_string)
    ld = LinearDensity(selection, binsize=5).run()
    with pytest.warns(DeprecationWarning):
        assert_allclose(ld.results.x.pos, ld.results.x.mass_density)
        assert_allclose(ld.results.x.pos_std, ld.results.x.mass_density_stddev)
        assert_allclose(ld.results.x.char, ld.results.x.charge_density)
        assert_allclose(
            ld.results.x.char_std, ld.results.x.charge_density_stddev
        )
        for key in testdict.keys():
            assert_allclose(
                ld.results["x"][key], ld.results["x"][testdict[key]]
            )

    # Check that no DeprecationWarning is raised with new attributes
    with no_deprecated_call():
        ld.results.x.mass_density
        ld.results.x.mass_density_stddev
        ld.results.x.charge_density
        ld.results.x.charge_density_stddev


# TODO: deprecated, remove in 3.0.0
# the parallelization here is not related to the parallelization through
# the AnalysisBase, so it is tested only in serial
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
    with pytest.warns(
        DeprecationWarning, match="`_add_other_results` is deprecated!"
    ):
        ld1._add_other_results(ld2)
    assert_allclose(
        ld1.results.z.mass_density, ld_whole.results.z.mass_density
    )
    assert_allclose(
        ld1.results.x.mass_density, ld_whole.results.x.mass_density
    )


def test_class_is_parallelizable():
    assert (
        mda.analysis.lineardensity.LinearDensity._analysis_algorithm_is_parallelizable
        == True
    )


def test_supported_backends():
    assert (
        mda.analysis.lineardensity.LinearDensity.get_supported_backends()
        == (
            "serial",
            "multiprocessing",
            "dask",
        )
    )


@pytest.mark.parametrize(
    "binsize, expected_nbins",
    [
        (0.1, 200),
        (0.2, 100),
        (0.05, 400),
        (0.04, 500),
        (0.02, 1000),
    ],
)
def test_lineardensity_floating_point_bins(binsize, expected_nbins):
    """Test that decimal bin sizes do not lose bins due to float truncation (Issue #4476)."""
    n_atoms = 4
    u = mda.Universe.empty(
        n_atoms=n_atoms,
        n_residues=1,
        atom_resindex=[0] * n_atoms,
    )
    u.add_TopologyAttr("masses", np.ones(n_atoms))
    u.add_TopologyAttr("charges", np.zeros(n_atoms))
    coords = np.array(
        [
            [
                [1.0, 1.0, 1.0],
                [5.0, 5.0, 5.0],
                [10.0, 10.0, 10.0],
                [15.0, 15.0, 15.0],
            ]
        ]
    )
    u.trajectory = get_reader_for(coords)(coords, order="fac", n_atoms=n_atoms)
    u.trajectory[0].dimensions = np.array([20.0, 20.0, 20.0, 90.0, 90.0, 90.0])

    ld = LinearDensity(u.atoms, binsize=binsize).run()
    assert ld.nbins == expected_nbins
    assert len(ld.results.x.hist_bin_edges) == expected_nbins + 1
    assert_allclose(ld.results.x.slice_volume, 8000.0 / expected_nbins)


@pytest.mark.parametrize(
    "binsize, expected_nbins",
    [
        (3.0, 6),
        (7.0, 2),
        (6.0, 3),
    ],
)
def test_lineardensity_non_clean_bins(binsize, expected_nbins):
    """Test that non-divisible bin sizes preserve floor division behavior (Issue #4476)."""
    n_atoms = 4
    u = mda.Universe.empty(
        n_atoms=n_atoms,
        n_residues=1,
        atom_resindex=[0] * n_atoms,
    )
    u.add_TopologyAttr("masses", np.ones(n_atoms))
    u.add_TopologyAttr("charges", np.zeros(n_atoms))
    coords = np.array(
        [
            [
                [1.0, 1.0, 1.0],
                [5.0, 5.0, 5.0],
                [10.0, 10.0, 10.0],
                [15.0, 15.0, 15.0],
            ]
        ]
    )
    u.trajectory = get_reader_for(coords)(coords, order="fac", n_atoms=n_atoms)
    u.trajectory[0].dimensions = np.array([20.0, 20.0, 20.0, 90.0, 90.0, 90.0])

    ld = LinearDensity(u.atoms, binsize=binsize).run()
    assert ld.nbins == expected_nbins
    assert len(ld.results.x.hist_bin_edges) == expected_nbins + 1
    assert_allclose(ld.results.x.slice_volume, 8000.0 / expected_nbins)


@pytest.mark.parametrize("binsize", [0.1, 0.25, 0.05])
def test_lineardensity_no_grid_aliasing_spikes(binsize):
    """Histogram bin edges must not produce spike/dip artifacts
    when atom coordinates coincide with the discrete float grid.

    Without the internal edge offset, histogram boundaries that
    coincide with discrete coordinate lattice points can cause
    atoms on boundaries to be systematically assigned to the
    adjacent bin due to float precision, creating artificial
    spike/dip dipole artifacts (Issue #4476).
    """
    box_len = 20.0
    nbins = int(np.floor(box_len / binsize + 1e-12))
    n_per_bin = int(binsize / 0.01)

    # Place one atom at every 0.01 A position (XTC precision)
    xs = np.arange(0.0, box_len, 0.01).astype(np.float32)
    n_atoms = len(xs)

    u = mda.Universe.empty(
        n_atoms=n_atoms,
        n_residues=1,
        atom_resindex=[0] * n_atoms,
    )
    u.add_TopologyAttr("masses", np.ones(n_atoms))
    u.add_TopologyAttr("charges", np.zeros(n_atoms))

    coords = np.zeros((1, n_atoms, 3), dtype=np.float32)
    coords[0, :, 0] = xs
    coords[0, :, 1] = 1.0
    coords[0, :, 2] = 1.0

    u.trajectory = get_reader_for(coords)(coords, order="fac", n_atoms=n_atoms)
    u.trajectory[0].dimensions = np.array(
        [box_len, box_len, box_len, 90.0, 90.0, 90.0]
    )

    ld = LinearDensity(u.atoms, binsize=binsize).run()

    # With offset internal bin edges, every bin should get exactly
    # n_per_bin atoms of unit mass, so mass_density should
    # be uniform across all occupied bins (the x dimension).
    x_mass = ld.results.x.mass_density[:nbins]
    expected_val = x_mass.mean()
    assert_allclose(
        x_mass,
        expected_val,
        atol=expected_val * 0.01,
        err_msg=(
            f"Grid-aliasing spike detected at binsize={binsize}: "
            f"min={x_mass.min()}, max={x_mass.max()}, "
            f"expected uniform ~{expected_val}"
        ),
    )
