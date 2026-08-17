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
import pytest

import MDAnalysis as mda
from MDAnalysis.analysis import polymer
from MDAnalysis.exceptions import NoDataError
from MDAnalysis.core.topologyattrs import Bonds

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from numpy.testing import assert_allclose, assert_almost_equal, assert_equal

from MDAnalysisTests.datafiles import Plength, TRZ_psf, TRZ


def test_class_is_parallelizable():
    assert polymer.PersistenceLength._analysis_algorithm_is_parallelizable


def test_supported_backends():
    assert polymer.PersistenceLength.get_supported_backends() == (
        "serial",
        "multiprocessing",
        "dask",
    )


class TestPersistenceLength(object):
    @staticmethod
    @pytest.fixture()
    def u():
        return mda.Universe(Plength)

    @pytest.fixture()
    def p(self, u):
        ags = [r.atoms.select_atoms("name C* N*") for r in u.residues]
        p = polymer.PersistenceLength(ags)
        return p

    @pytest.fixture()
    def p_run(self, p, client_PersistenceLength):
        return p.run(**client_PersistenceLength)

    def test_ag_ValueError(self, u):
        ags = [u.atoms[:10], u.atoms[10:110]]
        with pytest.raises(ValueError):
            polymer.PersistenceLength(ags)

    def test_run(self, p_run):
        assert len(p_run.results.bond_autocorrelation) == 280

    def test_lb(self, p_run):
        assert_almost_equal(p_run.results.lb, 1.485, 3)

    def test_fit(self, p_run):
        assert_almost_equal(p_run.results.lp, 6.504, 3)
        assert len(p_run.results.fit) == len(
            p_run.results.bond_autocorrelation
        )

    def test_raise_NoDataError(self, p):
        # Ensure that a NoDataError is raised if perform_fit()
        # is called before the run() method of AnalysisBase
        with pytest.raises(NoDataError):
            p._perform_fit()

    def test_plot_ax_return(self, p_run):
        """Ensure that a matplotlib axis object is
        returned when plot() is called."""
        assert isinstance(p_run.plot(), matplotlib.axes.Axes)

    def test_plot_with_ax(self, p_run):
        fig, ax = plt.subplots()
        ax2 = p_run.plot(ax=ax)
        assert ax2 is ax

    def test_current_axes(self, p_run):
        fig, ax = plt.subplots()
        ax2 = p_run.plot(ax=None)
        assert ax2 is not ax

    @pytest.mark.parametrize("attr", ("lb", "lp", "fit"))
    def test(self, p_run, attr):
        wmsg = f"The `{attr}` attribute was deprecated in MDAnalysis 2.0.0"
        with pytest.warns(DeprecationWarning, match=wmsg):
            getattr(p_run, attr) is p_run.results[attr]


class TestFitExponential(object):
    x = np.linspace(0, 250, 251)
    a_ref = 20.0
    y = np.exp(-x / a_ref)

    def test_fit_simple(self):
        a = polymer.fit_exponential_decay(self.x, self.y)
        assert a == self.a_ref

    def test_fit_noisy(self):
        noise = np.sin(self.x) * 0.01
        y2 = noise + self.y

        a = polymer.fit_exponential_decay(self.x, y2)

        assert_almost_equal(a, self.a_ref, decimal=3)
        # assert np.rint(a) == self.a_ref


class TestSortBackbone(object):
    @staticmethod
    @pytest.fixture(scope="class")
    def u():
        return mda.Universe(TRZ_psf, TRZ)

    def test_missing_bonds(self):
        u = mda.Universe(Plength)

        with pytest.raises(NoDataError):
            polymer.sort_backbone(u.atoms[:10])

    def test_sortbb(self, u):
        # grab backbone atoms out of order
        # 0 1 4 6 8 - correct
        ag = u.atoms[[4, 1, 0, 8, 6]]

        s_ag = polymer.sort_backbone(ag)

        assert_equal(s_ag.ids, [0, 1, 4, 6, 8])

    def test_branches(self, u):
        # includes side branches, can't sort
        bad_ag = u.atoms[:10]  # include -H etc

        with pytest.raises(ValueError, match="branches or isolated atoms"):
            polymer.sort_backbone(bad_ag)

    def test_isolated(self, u):
        u = mda.Universe.empty(4, trajectory=True)
        bondlist = [(0, 1), (1, 2)]
        u.add_TopologyAttr(Bonds(bondlist))
        with pytest.raises(ValueError, match="branches or isolated atoms"):
            polymer.sort_backbone(u.atoms)

    def test_missing_internal(self, u):
        u = mda.Universe.empty(4, trajectory=True)
        bondlist = [(0, 1), (2, 3)]
        u.add_TopologyAttr(Bonds(bondlist))
        with pytest.raises(ValueError, match="Backbone connectivity invalid"):
            polymer.sort_backbone(u.atoms)

    def test_circular(self):
        u = mda.Universe.empty(6, trajectory=True)
        # circular structure
        bondlist = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]
        u.add_TopologyAttr(Bonds(bondlist))
        with pytest.raises(ValueError, match="Cyclical"):
            polymer.sort_backbone(u.atoms)


def _reference_bond_autocorrelation(chains, frame_indices):
    """Independently compute the expected normalized bond autocorrelation
    for exactly the given sequence of trajectory frame indices.

    Mirrors the mathematical definition directly from atom positions,
    without calling into :class:`~MDAnalysis.analysis.polymer.
    PersistenceLength` or reusing any of its results, so this stays a
    meaningful check even if the class under test is itself broken.
    Frame indices are visited explicitly (in the order given, duplicates
    allowed) rather than via a slice, so this can independently mirror
    any of `run`'s `start`/`stop`/`step`/`frames` selections. See Issue
    #5453.
    """
    chainlength = len(chains[0])
    raw = np.zeros(chainlength - 1, dtype=np.float64)

    universe = chains[0].universe
    n_frames = 0
    for frame in frame_indices:
        universe.trajectory[frame]
        n_frames += 1
        for chain in chains:
            vecs = chain.positions[1:] - chain.positions[:-1]
            vecs = vecs / np.sqrt((vecs * vecs).sum(axis=1))[:, None]
            inner_pr = np.inner(vecs, vecs)
            for i in range(chainlength - 1):
                raw[: (chainlength - 1) - i] += inner_pr[i, i:]

    norm = np.linspace(chainlength - 1, 1, chainlength - 1)
    norm *= len(chains) * n_frames
    return raw / norm, n_frames


class TestPersistenceLengthSlicedNormalization(object):
    # Regression test for Issue #5453.
    @staticmethod
    @pytest.fixture(scope="class")
    def u():
        return mda.Universe(TRZ_psf, TRZ)

    @staticmethod
    @pytest.fixture(scope="class")
    def chains(u):
        backbones = [
            chain.select_atoms("not name O* H*") for chain in u.atoms.fragments
        ]
        return [polymer.sort_backbone(bb) for bb in backbones]

    def test_full_trajectory_unchanged(self, chains, u):
        p = polymer.PersistenceLength(chains).run()

        assert p.n_frames == u.trajectory.n_frames

        expected, n_frames = _reference_bond_autocorrelation(
            chains, range(u.trajectory.n_frames)
        )
        assert n_frames == u.trajectory.n_frames
        # atol/rtol (looser than a plain decimal=6 comparison) account for
        # PersistenceLength accumulating raw_bond_autocorr in float32 while
        # this reference accumulates in float64; their rounding drift is
        # a couple of float32 ULPs at this trajectory's frame/chain count,
        # not a sign of disagreement.
        assert_allclose(
            p.results.bond_autocorrelation, expected, rtol=1e-4, atol=5e-6
        )

    def test_sliced_run_normalized_by_frames_analyzed(self, chains, u):
        n_sliced_frames = 3
        assert n_sliced_frames < u.trajectory.n_frames

        p = polymer.PersistenceLength(chains).run(stop=n_sliced_frames)

        assert p.n_frames == n_sliced_frames
        assert p.n_frames != u.trajectory.n_frames

        expected, n_frames = _reference_bond_autocorrelation(
            chains, range(n_sliced_frames)
        )
        assert n_frames == n_sliced_frames
        # atol/rtol (looser than a plain decimal=6 comparison) account for
        # PersistenceLength accumulating raw_bond_autocorr in float32 while
        # this reference accumulates in float64; their rounding drift is
        # a couple of float32 ULPs at this trajectory's frame/chain count,
        # not a sign of disagreement.
        assert_allclose(
            p.results.bond_autocorrelation, expected, rtol=1e-4, atol=5e-6
        )

        # Sanity check: reproduce the old (buggy) behavior by normalizing
        # the *sliced* raw sum with the *total* trajectory frame count
        # instead of n_sliced_frames, and confirm that no longer matches
        # the class's output -- i.e. this test is sensitive to the bug
        # described in Issue #5453.
        wrong = expected * n_sliced_frames / u.trajectory.n_frames
        assert not np.allclose(p.results.bond_autocorrelation, wrong)

    @pytest.mark.parametrize(
        "run_kwargs, frame_indices",
        [
            pytest.param({"step": 2}, [0, 2, 4], id="step_only"),
            pytest.param(
                {"start": 1, "stop": 6, "step": 2},
                [1, 3, 5],
                id="start_stop_step",
            ),
            pytest.param({"frames": [0, 2, 5]}, [0, 2, 5], id="frames_list"),
            pytest.param(
                {"frames": [4, 1, 3]},
                [4, 1, 3],
                id="frames_list_unordered",
            ),
        ],
    )
    def test_various_frame_selections_normalized_correctly(
        self, chains, u, run_kwargs, frame_indices
    ):
        # Regression test for Issue #5453: the normalization fix relies
        # on self.n_frames, which AnalysisBase derives the same way
        # regardless of whether the run was sliced via start/stop/step
        # or via an explicit frames= list. Exercise more than just a
        # plain stop= slice to guard against that assumption breaking.
        p = polymer.PersistenceLength(chains).run(**run_kwargs)

        assert p.n_frames == len(frame_indices)
        assert p.n_frames != u.trajectory.n_frames

        expected, n_frames = _reference_bond_autocorrelation(
            chains, frame_indices
        )
        assert n_frames == len(frame_indices)
        # atol/rtol (looser than a plain decimal=6 comparison) account for
        # PersistenceLength accumulating raw_bond_autocorr in float32 while
        # this reference accumulates in float64; their rounding drift is
        # a couple of float32 ULPs at this trajectory's frame/chain count,
        # not a sign of disagreement.
        assert_allclose(
            p.results.bond_autocorrelation, expected, rtol=1e-4, atol=5e-6
        )

    def test_sliced_run_parallel_matches_serial(self, chains, u):
        # Regression test for Issue #5453: the bug was introduced by the
        # multi-worker refactor in PR #5074, but neither existing test
        # actually exercised backend="multiprocessing" together with a
        # sliced run. _conclude() only ever runs once, on the main
        # process, after worker results are merged -- so self.n_frames
        # at normalization time should be the total sliced frame count
        # regardless of how work was chunked across workers. Confirm
        # that holds, rather than only checking it by hand.
        n_sliced_frames = 5
        assert n_sliced_frames < u.trajectory.n_frames

        serial = polymer.PersistenceLength(chains).run(stop=n_sliced_frames)
        parallel = polymer.PersistenceLength(chains).run(
            stop=n_sliced_frames, backend="multiprocessing", n_workers=2
        )

        assert serial.n_frames == parallel.n_frames == n_sliced_frames
        assert_allclose(
            serial.results.bond_autocorrelation,
            parallel.results.bond_autocorrelation,
            rtol=1e-4,
            atol=5e-6,
        )

    def test_empty_frame_selection_raises_clear_error(self, chains):
        # Edge case found while verifying the fix for Issue #5453: an
        # empty frame selection drives self.n_frames to 0, which used
        # to silently divide by the (nonzero) total trajectory frame
        # count and return a fabricated all-zero curve with no warning.
        # After the fix, the zero-frame divisor would instead produce
        # NaNs that surface as an opaque scipy error deep inside
        # curve_fit. Raise a clear, immediate error instead.
        with pytest.raises(ValueError, match="analyzed zero frames"):
            polymer.PersistenceLength(chains).run(start=2, stop=2)

    def test_single_atom_chain_unaffected_by_normalization_fix(self, u):
        # Not a regression from Issue #5453's fix: a chain with only one
        # atom has no bonds to correlate (chainlength - 1 == 0), so
        # raw_bond_autocorr and norm are both length-0 arrays regardless
        # of the run() call's start/stop/step. This has always raised
        # from inside _perform_fit()'s curve_fit call, independent of
        # which frame count _conclude() normalizes by. Documented here
        # so it doesn't get mistaken for something this PR should fix.
        one_atom_chain = [u.atoms.fragments[0][:1]]
        with pytest.raises(ValueError):
            polymer.PersistenceLength(one_atom_chain).run()

    @pytest.mark.xfail(
        reason=(
            "Pre-existing AnalysisBase bug, independent of Issue #5453: "
            "_setup_computation_groups checks `isinstance(obj, bool)`, "
            "which is False for numpy bool scalars, so a numpy boolean "
            "frames= mask is treated as an integer index array (True/"
            "False read as 1/0) instead of a mask. A plain Python list "
            "of bools is unaffected. Left xfail to document the "
            "framework-level issue without attempting to fix "
            "AnalysisBase in this PR."
        ),
        strict=True,
    )
    def test_numpy_bool_mask_frames_selection(self, chains, u):
        mask = np.zeros(u.trajectory.n_frames, dtype=bool)
        mask[[0, 2, 4]] = True

        p = polymer.PersistenceLength(chains).run(frames=mask)

        assert p.n_frames == int(mask.sum())
        expected, n_frames = _reference_bond_autocorrelation(
            chains, np.flatnonzero(mask)
        )
        assert_allclose(
            p.results.bond_autocorrelation, expected, rtol=1e-4, atol=5e-6
        )
