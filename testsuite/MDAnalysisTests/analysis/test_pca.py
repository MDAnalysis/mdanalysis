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
from __future__ import print_function, absolute_import

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.pca import (PCA, cosine_content,
                                     rmsip, cumulative_overlap)

from numpy.testing import (assert_almost_equal, assert_equal,
                           assert_array_almost_equal)

from MDAnalysisTests.datafiles import (PSF, DCD, RANDOM_WALK, RANDOM_WALK_TOPO,
                                       waterPSF, waterDCD)
import pytest

SELECTION = 'backbone and name CA and resid 1-10'


@pytest.fixture(scope='module')
def u():
    return mda.Universe(PSF, DCD)


@pytest.fixture(scope='module')
def u_aligned():
    u = mda.Universe(PSF, DCD, in_memory=True)
    align.AlignTraj(u, u, select=SELECTION).run()
    return u


@pytest.fixture(scope='module')
def pca(u):
    u.transfer_to_memory()
    return PCA(u, select=SELECTION).run()


@pytest.fixture(scope='module')
def pca_aligned(u):
    return PCA(u, select=SELECTION, align=True).run()


def test_cov(pca, u):
    atoms = u.select_atoms(SELECTION)
    xyz = np.zeros((pca.n_frames, atoms.n_atoms * 3))
    for i, ts in enumerate(u.trajectory):
        xyz[i] = atoms.positions.ravel()
    cov = np.cov(xyz, rowvar=0)
    assert_array_almost_equal(pca.cov, cov, 4)


def test_cum_var(pca):
    assert_almost_equal(pca.cumulated_variance[-1], 1)
    l = pca.cumulated_variance
    l = np.sort(l)
    assert_almost_equal(pca.cumulated_variance, l, 5)


def test_pcs(pca):
    assert_equal(pca.p_components.shape, (pca._n_atoms * 3, pca._n_atoms * 3))


def test_pcs_n_components(u):
    pca = PCA(u, select=SELECTION).run()
    assert_equal(pca.n_components, pca._n_atoms*3)
    assert_equal(pca.p_components.shape, (pca._n_atoms * 3, pca._n_atoms * 3))
    pca.n_components = 10
    assert_equal(pca.n_components, 10)
    assert_equal(pca.p_components.shape, (pca._n_atoms * 3, 10))


def test_different_steps(pca, u):
    atoms = u.select_atoms(SELECTION)
    dot = pca.transform(atoms, start=5, stop=7, step=1)
    assert_equal(dot.shape, (2, atoms.n_atoms*3))


def test_transform_different_atoms(pca, u):
    atoms = u.select_atoms('backbone and name N and resid 1-10')
    with pytest.warns(UserWarning):
        pca.transform(atoms, start=5, stop=7, step=1)


def test_transform_rerun(u):
    atoms = u.select_atoms('bynum 1-10')
    u.transfer_to_memory()
    pca = PCA(u, select='bynum 1-10').run(stop=5)
    dot = pca.transform(atoms)
    assert_equal(dot.shape, (98, atoms.n_atoms * 3))


def test_pca_not_run(u):
    atoms = u.select_atoms('bynum 1-10')
    u.transfer_to_memory()
    pca = PCA(u, select='bynum 1-10')
    with pytest.raises(ValueError):
        dot = pca.transform(atoms, stop=5)


def test_no_frames(u):
    atoms = u.select_atoms(SELECTION)
    u.transfer_to_memory()
    with pytest.raises(ValueError):
        PCA(u, select=SELECTION).run(stop=1)


def test_transform(pca, u):
    ag = u.select_atoms(SELECTION)
    pca_space = pca.transform(ag, n_components=1)
    assert_equal(pca_space.shape, (u.trajectory.n_frames, 1))


def test_transform_mismatch(pca, u):
    with pytest.raises(ValueError):
        pca.transform(u, n_components=1)


def test_transform_universe():
    u1 = mda.Universe(waterPSF, waterDCD)
    u2 = mda.Universe(waterPSF, waterDCD)
    pca_test = PCA(u1).run()
    pca_test.transform(u2)


def test_cosine_content():
    rand = mda.Universe(RANDOM_WALK_TOPO, RANDOM_WALK)
    pca_random = PCA(rand).run()
    dot = pca_random.transform(rand.atoms)
    content = cosine_content(dot, 0)
    assert_almost_equal(content, .99, 1)


def test_mean_shape(pca_aligned, u):
    atoms = u.select_atoms(SELECTION)
    assert_equal(pca_aligned.mean.shape[0], atoms.n_atoms * 3)


def test_calculate_mean(pca_aligned, u, u_aligned):
    ag = u_aligned.select_atoms(SELECTION)
    coords = u_aligned.trajectory.coordinate_array[:, ag.ix]
    assert_almost_equal(pca_aligned.mean, coords.mean(
        axis=0).ravel(), decimal=5)


def test_given_mean(pca_aligned, u):
    pca = PCA(u, select=SELECTION, align=False,
              mean=pca_aligned._mean).run()
    assert_almost_equal(pca_aligned.cov, pca.cov, decimal=5)


def test_alignment(pca_aligned, u, u_aligned):
    pca_pre_align = PCA(u_aligned, select=SELECTION, align=False).run()
    assert_almost_equal(pca_aligned.mean, pca_pre_align.mean)
    assert_almost_equal(pca_aligned.cov, pca_pre_align.cov)


def test_covariance_norm(pca_aligned, u):
    assert_almost_equal(np.linalg.norm(pca_aligned.cov), 0.96799758, decimal=5)


def test_pca_rmsip_self(pca):
    assert_almost_equal(pca.rmsip(pca), 1.0)


def test_rmsip_ortho(pca):
    value = rmsip(pca.p_components[:, :10].T,
                  pca.p_components[:, 10:20].T)
    assert_almost_equal(value, 0.0)


def test_pytest_too_many_components(pca):
    with pytest.raises(ValueError) as exc:
        pca.rmsip(pca, n_components=(1, 2, 3))
    assert 'Too many values' in str(exc.value)


def test_asymmetric_rmsip(pca):
    a = pca.rmsip(pca, n_components=(10, 4))
    b = pca.rmsip(pca, n_components=(4, 10))

    assert abs(a-b) > 0.1, 'RMSIP should be asymmetric'
    assert_almost_equal(b, 1.0)


def test_pca_cumulative_overlap_self(pca):
    value = pca.cumulative_overlap(pca, i=1)
    assert_almost_equal(value, 1.0)


def test_cumulative_overlap_ortho(pca):
    pcs = pca.p_components
    value = cumulative_overlap(pcs[:, 11].T, pcs.T, n_components=10)
    assert_almost_equal(value, 0.0)


@pytest.mark.parametrize(
    'method', ['rmsip',
               'cumulative_overlap'])
def test_compare_not_run_other(u, pca, method):
    pca2 = PCA(u)
    func = getattr(pca, method)
    with pytest.raises(ValueError) as exc:
        func(pca2)
    assert 'Call run()' in str(exc.value)


@pytest.mark.parametrize(
    'method', ['rmsip',
               'cumulative_overlap'])
def test_compare_not_run_self(u, pca, method):
    pca2 = PCA(u)
    func = getattr(pca2, method)
    with pytest.raises(ValueError) as exc:
        func(pca)
    assert 'Call run()' in str(exc.value)


@pytest.mark.parametrize(
    'method', ['rmsip',
               'cumulative_overlap'])
def test_compare_wrong_class(u, pca, method):
    func = getattr(pca, method)
    with pytest.raises(ValueError) as exc:
        func(3)
    assert 'must be another PCA class' in str(exc.value)
