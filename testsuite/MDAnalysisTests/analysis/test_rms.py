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
import os

import MDAnalysis
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align

from numpy.testing import assert_equal, assert_almost_equal

import numpy as np
import pytest

from MDAnalysis.exceptions import SelectionError, NoDataError
from MDAnalysisTests.datafiles import GRO, XTC, rmsfArray, PSF, DCD


class Testrmsd(object):
    shape = (5, 3)
    # vectors with length one
    ones = np.ones(shape) / np.sqrt(3)

    @pytest.fixture()
    def a(self):
        return self.ones * np.arange(1, 6)[:, np.newaxis]

    @pytest.fixture()
    def b(self, a):
        return a + self.ones

    @pytest.fixture()
    def u(self):
        u = mda.Universe(PSF, DCD)
        return u

    @pytest.fixture()
    def u2(self):
        u = mda.Universe(PSF, DCD)
        return u

    @pytest.fixture()
    def p_first(self, u):
        u.trajectory[0]
        return u.select_atoms('protein')

    @pytest.fixture()
    def p_last(self, u2):
        u2.trajectory[-1]
        return u2.select_atoms('protein')

    # internal test
    def test_p_frames(self, p_first, p_last):
        # check that these fixtures are really different
        assert p_first.universe.trajectory.ts.frame != p_last.universe.trajectory.ts.frame
        assert not np.allclose(p_first.positions, p_last.positions)

    def test_no_center(self, a, b):
        rmsd = rms.rmsd(a, b, center=False)
        assert_almost_equal(rmsd, 1.0)

    def test_center(self, a, b):
        rmsd = rms.rmsd(a, b, center=True)
        assert_almost_equal(rmsd, 0.0)

    def test_list(self, a, b):
        rmsd = rms.rmsd(a.tolist(),
                        b.tolist(),
                        center=False)
        assert_almost_equal(rmsd, 1.0)

    @pytest.mark.parametrize('dtype', [np.float32, np.float64])
    def test_superposition(self, a, b, u, dtype):
        bb = u.atoms.select_atoms('backbone')
        a = bb.positions.copy()
        u.trajectory[-1]
        b = bb.positions.copy()

        a, b = a.astype(dtype), b.astype(dtype)

        rmsd = rms.rmsd(a, b, superposition=True)
        assert_almost_equal(rmsd, 6.820321761927005)

    def test_weights(self, a, b):
        weights = np.zeros(len(a))
        weights[0] = 1
        weights[1] = 1
        weighted = rms.rmsd(a, b, weights=weights)
        firstCoords = rms.rmsd(a[:2], b[:2])
        assert_almost_equal(weighted, firstCoords)

    def test_weights_and_superposition_1(self, u):
        weights = np.ones(len(u.trajectory[0]))
        weighted = rms.rmsd(u.trajectory[0], u.trajectory[1],
                            weights=weights, superposition=True)
        firstCoords = rms.rmsd(u.trajectory[0], u.trajectory[1],
                               superposition=True)
        assert_almost_equal(weighted, firstCoords, decimal=5)

    def test_weights_and_superposition_2(self, u):
        weights = np.zeros(len(u.trajectory[0]))
        weights[:100] = 1
        weighted = rms.rmsd(u.trajectory[0], u.trajectory[-1],
                            weights=weights, superposition=True)
        firstCoords = rms.rmsd(u.trajectory[0][:100],
                               u.trajectory[-1][:100],
                               superposition=True)
        # very close to zero, change significant decimal places to 5
        assert_almost_equal(weighted, firstCoords, decimal=5)

    def test_unequal_shape(self):
        a = np.ones((4, 3))
        b = np.ones((5, 3))
        with pytest.raises(ValueError):
            rms.rmsd(a, b)

    def test_wrong_weights(self, a, b):
        w = np.ones(2)
        with pytest.raises(ValueError):
            rms.rmsd(a, b, w)

    def test_with_superposition_smaller(self, p_first, p_last):
        A = p_first.positions
        B = p_last.positions
        rmsd = rms.rmsd(A, B)
        rmsd_superposition = rms.rmsd(A, B, center=True, superposition=True)
        # by design the super positioned rmsd is smaller
        assert rmsd > rmsd_superposition

    def test_with_superposition_equal(self, p_first, p_last):
        align.alignto(p_first, p_last)
        A = p_first.positions
        B = p_last.positions
        rmsd = rms.rmsd(A, B)
        rmsd_superposition = rms.rmsd(A, B, center=True, superposition=True)
        assert_almost_equal(rmsd, rmsd_superposition, decimal=6)


class TestRMSD(object):
    @pytest.fixture()
    def universe(self):
        return MDAnalysis.Universe(PSF, DCD)

    @pytest.fixture()
    def outfile(self, tmpdir):
        return os.path.join(str(tmpdir), 'rmsd.txt')

    @pytest.fixture()
    def correct_values(self):
        return [[0, 1, 0], [49, 50, 4.68953]]

    @pytest.fixture()
    def correct_values_mass(self):
        return [[0, 1, 0], [49, 50, 4.74920]]

    @pytest.fixture()
    def correct_values_mass_add_ten(self):
        return [[0, 1, 0.0632], [49, 50, 4.7710]]

    @pytest.fixture()
    def correct_values_group(self):
        return [[0, 1, 0, 0, 0],
                [49, 50, 4.7857, 4.7048, 4.6924]]

    @pytest.fixture()
    def correct_values_backbone_group(self):
        return [[0, 1, 0, 0, 0],
                [49, 50, 4.6997, 1.9154, 2.7139]]

    def test_rmsd(self, universe, correct_values, client_RMSD):
        # client_RMSD is defined in testsuite/analysis/conftest.py
        # among with other testing fixtures. During testing, it will
        # collect all possible backends and reasonable number of workers
        # for a given AnalysisBase subclass, and extend the tests
        # to run with all of them.
        RMSD = MDAnalysis.analysis.rms.RMSD(universe, select='name CA')
        RMSD.run(step=49, **client_RMSD)
        assert_almost_equal(RMSD.results.rmsd, correct_values, 4,
                            err_msg="error: rmsd profile should match" +
                            "test values")

    def test_rmsd_frames(self, universe, correct_values, client_RMSD):
        RMSD = MDAnalysis.analysis.rms.RMSD(universe, select='name CA')
        RMSD.run(frames=[0, 49], **client_RMSD)
        assert_almost_equal(RMSD.results.rmsd, correct_values, 4,
                            err_msg="error: rmsd profile should match" +
                            "test values")

    def test_rmsd_unicode_selection(self, universe, correct_values, client_RMSD):
        RMSD = MDAnalysis.analysis.rms.RMSD(universe, select=u'name CA')
        RMSD.run(step=49, **client_RMSD)
        assert_almost_equal(RMSD.results.rmsd, correct_values, 4,
                            err_msg="error: rmsd profile should match" +
                            "test values")

    def test_rmsd_atomgroup_selections(self, universe, client_RMSD):
        # see Issue #1684
        R1 = MDAnalysis.analysis.rms.RMSD(universe.atoms,
                                          select="resid 1-30").run(**client_RMSD)
        R2 = MDAnalysis.analysis.rms.RMSD(universe.atoms.select_atoms("name CA"),
                                          select="resid 1-30").run(**client_RMSD)
        assert not np.allclose(R1.results.rmsd[:, 2], R2.results.rmsd[:, 2])

    def test_rmsd_single_frame(self, universe, client_RMSD):
        RMSD = MDAnalysis.analysis.rms.RMSD(universe, select='name CA',
                                            ).run(start=5, stop=6, **client_RMSD)
        single_frame = [[5, 6, 0.91544906]]
        assert_almost_equal(RMSD.results.rmsd, single_frame, 4,
                            err_msg="error: rmsd profile should match" +
                            "test values")

    def test_mass_weighted(self, universe, correct_values, client_RMSD):
        # mass weighting the CA should give the same answer as weighing
        # equally because all CA have the same mass
        RMSD = MDAnalysis.analysis.rms.RMSD(universe, select='name CA',
                                            weights='mass').run(step=49, **client_RMSD)

        assert_almost_equal(RMSD.results.rmsd, correct_values, 4,
                            err_msg="error: rmsd profile should match"
                            "test values")

    def test_custom_weighted(self, universe, correct_values_mass, client_RMSD):
        RMSD = MDAnalysis.analysis.rms.RMSD(universe, weights="mass").run(step=49, **client_RMSD)

        assert_almost_equal(RMSD.results.rmsd, correct_values_mass, 4,
                            err_msg="error: rmsd profile should match"
                            "test values")

    def test_weights_mass_is_mass_weighted(self, universe, client_RMSD):
        RMSD_mass = MDAnalysis.analysis.rms.RMSD(universe,
                                                 weights="mass").run(step=49, **client_RMSD)
        RMSD_cust = MDAnalysis.analysis.rms.RMSD(universe,
                                                 weights=universe.atoms.masses).run(step=49, **client_RMSD)
        assert_almost_equal(RMSD_mass.results.rmsd, RMSD_cust.results.rmsd, 4,
                            err_msg="error: rmsd profiles should match for 'mass' "
                            "and universe.atoms.masses")

    def test_custom_weighted_list(self, universe, correct_values_mass, client_RMSD):
        weights = universe.atoms.masses
        RMSD = MDAnalysis.analysis.rms.RMSD(universe,
                                            weights=list(weights)).run(step=49, **client_RMSD)
        assert_almost_equal(RMSD.results.rmsd, correct_values_mass, 4,
                            err_msg="error: rmsd profile should match" +
                            "test values")

    def test_custom_groupselection_weights_applied_1D_array(self, universe, client_RMSD):
        RMSD = MDAnalysis.analysis.rms.RMSD(universe,
                                            select='backbone',
                                            groupselections=['name CA and resid 1-5', 'name CA and resid 1'],
                                            weights=None,
                                            weights_groupselections=[[1, 0, 0, 0, 0], None]).run(step=49, 
                                                                                                 **client_RMSD
                                                                                                )

        assert_almost_equal(RMSD.results.rmsd.T[3], RMSD.results.rmsd.T[4], 4,
                            err_msg="error: rmsd profile should match "
                            "for applied weight array and selected resid")

    def test_custom_groupselection_weights_applied_mass(self, universe, correct_values_mass, client_RMSD):
        RMSD = MDAnalysis.analysis.rms.RMSD(universe,
                                            select='backbone',
                                            groupselections=['all', 'all'],
                                            weights=None,
                                            weights_groupselections=['mass',
                                                                     universe.atoms.masses]).run(step=49, 
                                                                                                **client_RMSD
                                                                                                )

        assert_almost_equal(RMSD.results.rmsd.T[3], RMSD.results.rmsd.T[4], 4,
                            err_msg="error: rmsd profile should match "
                            "between applied mass and universe.atoms.masses")

    def test_rmsd_scalar_weights_raises_ValueError(self, universe):
        with pytest.raises(ValueError):
            RMSD = MDAnalysis.analysis.rms.RMSD(
                universe, weights=42)

    def test_rmsd_string_weights_raises_ValueError(self, universe):
        with pytest.raises(ValueError):
            RMSD = MDAnalysis.analysis.rms.RMSD(
                universe, weights="Jabberwock")

    def test_rmsd_mismatched_weights_raises_ValueError(self, universe):
        with pytest.raises(ValueError):
            RMSD = MDAnalysis.analysis.rms.RMSD(
                universe, weights=universe.atoms.masses[:-1])

    def test_rmsd_misuse_weights_for_groupselection_raises_TypeError(self, universe):
        with pytest.raises(TypeError):
            RMSD = MDAnalysis.analysis.rms.RMSD(
                universe, groupselections=['all'],
                weights=[universe.atoms.masses, universe.atoms.masses[:-1]])

    def test_rmsd_mismatched_weights_in_groupselection_raises_ValueError(self, universe):
        with pytest.raises(ValueError):
            RMSD = MDAnalysis.analysis.rms.RMSD(
                universe, groupselections=['all'],
                weights=universe.atoms.masses,
                weights_groupselections = [universe.atoms.masses[:-1]])

    def test_rmsd_list_of_weights_wrong_length(self, universe):
        with pytest.raises(ValueError):
            RMSD = MDAnalysis.analysis.rms.RMSD(
                universe, groupselections=['backbone', 'name CA'],
                weights='mass',
                weights_groupselections=[None])

    def test_rmsd_group_selections(self, universe, correct_values_group, client_RMSD):
        RMSD = MDAnalysis.analysis.rms.RMSD(universe,
                                            groupselections=['backbone', 'name CA']
                                            ).run(step=49, **client_RMSD)
        assert_almost_equal(RMSD.results.rmsd, correct_values_group, 4,
                            err_msg="error: rmsd profile should match"
                            "test values")

    def test_rmsd_backbone_and_group_selection(self, universe,
                                               correct_values_backbone_group,
                                               client_RMSD):
        RMSD = MDAnalysis.analysis.rms.RMSD(
            universe,
            reference=universe,
            select="backbone",
            groupselections=['backbone and resid 1:10',
                             'backbone and resid 10:20']).run(step=49, **client_RMSD)
        assert_almost_equal(
            RMSD.results.rmsd, correct_values_backbone_group, 4,
            err_msg="error: rmsd profile should match test values")

    def test_ref_length_unequal_len(self, universe):
        reference = MDAnalysis.Universe(PSF, DCD)
        reference.atoms = reference.atoms[:-1]
        with pytest.raises(SelectionError):
            RMSD = MDAnalysis.analysis.rms.RMSD(universe,
                                                reference=reference)

    def test_mass_mismatches(self, universe):
        reference = MDAnalysis.Universe(PSF, DCD)
        reference.atoms.masses = 10
        with pytest.raises(SelectionError):
            RMSD = MDAnalysis.analysis.rms.RMSD(universe,
                                                reference=reference)

    def test_ref_mobile_mass_mismapped(self, universe,correct_values_mass_add_ten, client_RMSD):
        reference = MDAnalysis.Universe(PSF, DCD)
        universe.atoms.masses = universe.atoms.masses + 10
        RMSD = MDAnalysis.analysis.rms.RMSD(universe,
                                                reference=reference,
                                                select='all',
                                                weights='mass',
                                                tol_mass=100)
        RMSD.run(step=49, **client_RMSD)
        assert_almost_equal(RMSD.results.rmsd, correct_values_mass_add_ten, 4,
                            err_msg="error: rmsd profile should match "
                            "between true values and calculated values")

    def test_group_selections_unequal_len(self, universe):
        reference = MDAnalysis.Universe(PSF, DCD)
        reference.atoms[0].residue.resname = 'NOTMET'
        with pytest.raises(SelectionError):
            RMSD = MDAnalysis.analysis.rms.RMSD(universe,
                                                reference=reference,
                                                groupselections=['resname MET', 'type NH3'])

    def test_rmsd_attr_warning(self, universe, client_RMSD):
        RMSD = MDAnalysis.analysis.rms.RMSD(
                universe, select='name CA').run(stop=2, **client_RMSD)

        wmsg = "The `rmsd` attribute was deprecated in MDAnalysis 2.0.0"
        with pytest.warns(DeprecationWarning, match=wmsg):
            assert_equal(RMSD.rmsd, RMSD.results.rmsd)


class TestRMSF(object):
    @pytest.fixture()
    def universe(self):
        return mda.Universe(GRO, XTC)

    def test_rmsf(self, universe, client_RMSF):
        rmsfs = rms.RMSF(universe.select_atoms('name CA'))
        rmsfs.run(**client_RMSF)
        test_rmsfs = np.load(rmsfArray)

        assert_almost_equal(rmsfs.results.rmsf, test_rmsfs, 5,
                            err_msg="error: rmsf profile should match test "
                            "values")

    def test_rmsf_single_frame(self, universe, client_RMSF):
        rmsfs = rms.RMSF(universe.select_atoms('name CA')).run(start=5, stop=6, **client_RMSF)

        assert_almost_equal(rmsfs.results.rmsf, 0, 5,
                            err_msg="error: rmsfs should all be zero")

    def test_rmsf_identical_frames(self, universe, tmpdir, client_RMSF):

        outfile = os.path.join(str(tmpdir), 'rmsf.xtc')

        # write a dummy trajectory of all the same frame
        with mda.Writer(outfile, universe.atoms.n_atoms) as W:
            for _ in range(universe.trajectory.n_frames):
                W.write(universe)

        universe = mda.Universe(GRO, outfile)
        rmsfs = rms.RMSF(universe.select_atoms('name CA'))
        rmsfs.run(**client_RMSF)
        assert_almost_equal(rmsfs.results.rmsf, 0, 5,
                            err_msg="error: rmsfs should all be 0")

    def test_rmsf_attr_warning(self, universe, client_RMSF):
        rmsfs = rms.RMSF(universe.select_atoms('name CA')).run(stop=2, **client_RMSF)

        wmsg = "The `rmsf` attribute was deprecated in MDAnalysis 2.0.0"
        with pytest.warns(DeprecationWarning, match=wmsg):
            assert_equal(rmsfs.rmsf, rmsfs.results.rmsf)


@pytest.mark.parametrize(
    "classname,is_parallelizable",
    [
        (MDAnalysis.analysis.rms.RMSD, True),
        (MDAnalysis.analysis.rms.RMSF, False),
    ]
)
def test_not_parallelizable(classname, is_parallelizable):
    assert classname._analysis_algorithm_is_parallelizable == is_parallelizable


@pytest.mark.parametrize(
    "classname,backends",
    [
        (MDAnalysis.analysis.rms.RMSD,  ('serial', 'multiprocessing', 'dask',)),
        (MDAnalysis.analysis.rms.RMSF, ('serial',)),
    ]
)
def test_supported_backends(classname, backends):
    assert classname.get_supported_backends() == backends
