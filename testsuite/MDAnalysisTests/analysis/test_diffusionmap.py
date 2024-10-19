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
import MDAnalysis
import MDAnalysis.analysis.diffusionmap as diffusionmap
import numpy as np
import pytest
from MDAnalysisTests.datafiles import PDB, XTC
from numpy.testing import assert_array_almost_equal, assert_allclose


@pytest.fixture(scope='module')
def u():
    return MDAnalysis.Universe(PDB, XTC)


@pytest.fixture(scope='module')
def dist(u):
    return diffusionmap.DistanceMatrix(u, select='backbone')


@pytest.fixture(scope='module')
def dmap(dist):
    d_map = diffusionmap.DiffusionMap(dist)
    d_map.run()
    return d_map


def test_eg(dist, dmap):
    eigvals = dmap.eigenvalues
    # number of frames is trajectory is now 10 vs. 98
    assert eigvals.shape == (dist.n_frames,)
    # makes no sense to test values here, no physical meaning


def test_dist_weights(u, client_DistanceMatrix):
    backbone = u.select_atoms('backbone')
    weights_atoms = np.ones(len(backbone.atoms))
    dist = diffusionmap.DistanceMatrix(u,
                                       select='backbone',
                                       weights=weights_atoms)
    dist.run(**client_DistanceMatrix, step=3)
    dmap = diffusionmap.DiffusionMap(dist)
    dmap.run()
    assert_array_almost_equal(dmap.eigenvalues, [1, 1, 1, 1], 4)
    assert_array_almost_equal(dmap._eigenvectors,
                              ([[0, 0, 1, 0],
                                [0, 0, 0, 1],
                                [-.707, -.707, 0, 0],
                                [.707, -.707, 0, 0]]), 2)


def test_dist_weights_frames(u, client_DistanceMatrix):
    backbone = u.select_atoms('backbone')
    weights_atoms = np.ones(len(backbone.atoms))
    dist = diffusionmap.DistanceMatrix(u,
                                       select='backbone',
                                       weights=weights_atoms)
    frames = np.arange(len(u.trajectory))
    dist.run(**client_DistanceMatrix, frames=frames[::3])
    dmap = diffusionmap.DiffusionMap(dist)
    dmap.run()
    assert_array_almost_equal(dmap.eigenvalues, [1, 1, 1, 1], 4)
    assert_array_almost_equal(dmap._eigenvectors,
                              ([[0, 0, 1, 0],
                                [0, 0, 0, 1],
                                [-.707, -.707, 0, 0],
                                [.707, -.707, 0, 0]]), 2)


def test_distvalues_ag_universe(u, client_DistanceMatrix):
    dist_universe = diffusionmap.DistanceMatrix(u, select='backbone').run(
        **client_DistanceMatrix
    )
    ag = u.select_atoms('backbone')
    dist_ag = diffusionmap.DistanceMatrix(ag).run(**client_DistanceMatrix)
    assert_allclose(dist_universe.results.dist_matrix,
                    dist_ag.results.dist_matrix)


def test_distvalues_ag_select(u, client_DistanceMatrix):
    dist_universe = diffusionmap.DistanceMatrix(u, select='backbone').run(
        **client_DistanceMatrix
    )
    ag = u.select_atoms('protein')
    dist_ag = diffusionmap.DistanceMatrix(ag, select='backbone').run(
        **client_DistanceMatrix
    )
    assert_allclose(dist_universe.results.dist_matrix,
                    dist_ag.results.dist_matrix)
                    

def test_different_steps(u):
    dmap = diffusionmap.DiffusionMap(u, select='backbone')
    dmap.run(step=3)
    assert dmap._eigenvectors.shape == (4, 4)


def test_transform(u, dmap):
    eigvects = dmap._eigenvectors
    n_eigenvectors = 4
    dmap = diffusionmap.DiffusionMap(u)
    dmap.run()
    diffusion_space = dmap.transform(n_eigenvectors, 1)
    assert diffusion_space.shape == (eigvects.shape[0], n_eigenvectors)


def test_long_traj(u):
    with pytest.warns(UserWarning, match='The distance matrix is very large'):
        dmap = diffusionmap.DiffusionMap(u)
        dmap._dist_matrix.run(stop=1)
        dmap._dist_matrix.n_frames = 5001
        dmap.run()


def test_updating_atomgroup(u):
    with pytest.warns(UserWarning, match='U must be a static AtomGroup'):
        resid_select = 'around 5 resname ALA'
        ag = u.select_atoms(resid_select, updating=True)
        dmap = diffusionmap.DiffusionMap(ag)
        dmap.run()

def test_not_universe_atomgroup_error(u):
    trj_only = u.trajectory
    with pytest.raises(ValueError, match='U is not a Universe or AtomGroup'):
        diffusionmap.DiffusionMap(trj_only)


def test_DistanceMatrix_attr_warning(u):
    dist = diffusionmap.DistanceMatrix(u, select='backbone').run(step=3)
    wmsg = f"The `dist_matrix` attribute was deprecated in MDAnalysis 2.0.0"
    with pytest.warns(DeprecationWarning, match=wmsg):
        assert getattr(dist, "dist_matrix") is dist.results.dist_matrix
