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

import MDAnalysis
import MDAnalysis.analysis.diffusionmap as diffusionmap
import numpy as np
import pytest
from MDAnalysisTests.datafiles import PDB, XTC
from numpy.testing import assert_array_almost_equal


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

weighted_evals = [1, 1, 1, 1]

weighted_evecs = ([[0, 0, 1, 0],
                   [0, 0, 0, 1],
                   [-.707, -.707, 0, 0],
                   [.707, -.707, 0, 0]])

def test_dist_weights(u):
    backbone = u.select_atoms('backbone')
    weights_atoms = np.ones(len(backbone.atoms))
    dist = diffusionmap.DistanceMatrix(u, select='backbone', weights=weights_atoms)
    dist.run(step=3)
    dmap = diffusionmap.DiffusionMap(dist)
    dmap.run()
    assert_array_almost_equal(dmap.eigenvalues, weighted_evals, 4)
    assert_array_almost_equal(dmap._eigenvectors,
                              weighted_evecs, 2)


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

def test_distance_array(u):
    backbone = u.select_atoms('backbone')
    weights_atoms = np.ones(len(backbone.atoms))
    dist = diffusionmap.DistanceMatrix(u, select='backbone', weights=weights_atoms)
    dist.run(step=3)
    dmap = diffusionmap.DiffusionMap(dist.dist_matrix)
    dmap.run()
    assert_array_almost_equal(dmap.eigenvalues, weighted_evals, 4)
    assert_array_almost_equal(dmap._eigenvectors,
                              weighted_evecs, 2)

def test_subset_frames(u):
    backbone = u.select_atoms('backbone')
    weights_atoms = np.ones(len(backbone.atoms))
    dist = diffusionmap.DistanceMatrix(u, select='backbone', weights=weights_atoms)
    dist.run(step=3)
    dmap = diffusionmap.DiffusionMap(dist.dist_matrix)
    dmap.run(step=2)
    assert_array_almost_equal(dmap.eigenvalues, weighted_evals[::2], 4)
    assert dmap._eigenvectors.shape == (2, 2)


@pytest.mark.parametrize(
    'dmat_run, dmap_run, frames', [
        ((5, 30, None), (4, 29, None), "4"),
        ((None, None, 2), (None, None, 1), "1 3 5")
    ])
def test_wrong_frames_error(dist, dmat_run, dmap_run, frames):
    dist.run(*dmat_run)
    dmap = diffusionmap.DiffusionMap(dist)
    with pytest.raises(ValueError) as exc:
        dmap.run(*dmap_run)
    assert "not in the distance matrix" in str(exc.value)
    assert frames in str(exc.value)

def test_too_many_frames_error(dist):
    dist.run(start=4, stop=10)
    dmap = diffusionmap.DiffusionMap(dist.dist_matrix)
    with pytest.raises(ValueError) as exc:
        dmap.run(stop=10)
    assert "not in the distance matrix" in str(exc.value)
    # run sets start=0 by default
    assert "6 7 8 9" in str(exc.value)

def not_square_array_error():
    with pytest.raises(ValueError) as exc:
        dmap = diffusionmap.DiffusionMap(np.zeros((3, 4)))
    assert "square" in str(exc.value)

def test_get_plotly_graphs(dmap):
    plotly = pytest.importorskip("plotly")
    n_frames = 5
    fig = dmap.plot_animated_transform(n_frames=n_frames, x=3, y=4)
    assert isinstance(fig, plotly.graph_objects.Figure)
    assert len(fig.layout['sliders'][0]['steps']) == n_frames
    assert fig.layout['xaxis']['title']['text'] == 'DC 3'
    assert fig.layout['yaxis']['title']['text'] == 'DC 4'