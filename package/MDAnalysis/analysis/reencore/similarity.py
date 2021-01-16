# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
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

import itertools

import numpy as np
from scipy.stats import gaussian_kde

from ...core.ensemble import Ensemble
from ..clustering import Clusters, methods
from ..align import AlignTraj
from ..rms import rmsd
from ..diffusionmap import DistanceMatrix
from ..pca import PCA
from . import utils

def clusters_to_indices(clusters, outlier_label=-1, outlier_index=-1):
    # convert labels to indices for ease
    cl_ids, cl_inv = np.unique(clusters, return_inverse=True)
    n_cl = len(cl_ids)
    i_ids = np.arange(n_cl)
    n_outliers = 0
    try:
        o = np.where(cl_ids == outlier_label)[0][0]
    except IndexError:
        pass
    else:
        i_ids[o] = outlier_index
        n_outliers = len(o)
        n_cl -= 1
    i_labels = i_ids[cl_inv]
    return i_labels, n_cl, n_outliers


def prepare_frames(ensemble, frames=None, start=None, stop=None, step=None,
                   estimate_error=False, n_bootstrap_samples=10,):
    if frames is None:
        start, stop, step = ensemble.trajectory.check_slice_indices(start, stop, step)
        frames = np.arange(start, stop, step)

    frames = np.array([frames])

    if estimate_error:
        edges = ensemble._frame_edges
        indices = [np.arange(edges[i], edges[i+1]) for i in range(len(edges)-1)]

        bs = [utils.get_bootstrap_frames(x, n_samples=n_bootstrap_samples) for x in indices]
        bs = [np.concatenate(x) for x in zip(*bs)]
        frames = np.concatenate([frames, bs])

    return frames


def prepare_ces(ensemble, clusters, start=None, stop=None, step=None,
                similarity_matrix=None, metric=rmsd, frames=None,
                estimate_error=False, n_bootstrap_samples=10,
                **kwargs):
    if callable(clusters):
        clusters = Clusters(clusters, **kwargs)

    frames = prepare_frames(ensemble, frames=frames, start=start, stop=stop,
                            step=step, estimate_error=estimate_error,
                            n_bootstrap_samples=n_bootstrap_samples)

    data = np.zeros((len(frames), frames[0].shape[0]))

    if isinstance(clusters, Clusters):
        if similarity_matrix is None:
            if clusters._data is not None:
                similarity_matrix = clusters._data
            else:
                dist_mat = DistanceMatrix(ensemble, metric=metric, **kwargs).run()
                similarity_matrix = -dist_mat.dist_matrix

        for i, fr in enumerate(frames):
            clusters.run(similarity_matrix[frames[0]][:, frames[0]])
            data[i] = clusters.data_labels

    else:
        for i, fr in enumerate(frames):
            data[i] = clusters[fr]
        
    return data
    

def ces(ensemble, clusters=methods.AffinityPropagation,
        select=None,
        outlier_label=-1, outlier_index=-1,
        estimate_error=False, **kwargs):
    if select is not None:
        ensemble = ensemble.select_atoms(select)
    clusters = prepare_ces(ensemble, clusters, estimate_error=estimate_error,
                           **kwargs)

    results = np.zeros((len(clusters), ensemble.n_universes, ensemble.n_universes))
    for i, group in enumerate(clusters):
        results[i] = _ces(ensemble, group, outlier_label=outlier_label,
                          outlier_index=outlier_index)
    
    if estimate_error:
        return results.mean(axis=0), results.std(axis=0)

    return results[0]


def _ces(ensemble, data, outlier_label=-1, outlier_index=-1):
    i_labels, n_cl, n_outliers = clusters_to_indices(data, outlier_label)
    
    # count membership in each cluster
    frames_per_cl = np.zeros((ensemble.n_universes, n_cl + n_outliers),
                             dtype=float)
    outlier_i = n_cl  # marker for new outlier columns
    for row, data in zip(frames_per_cl, ensemble.split_array(i_labels)):
        labels_, counts_ = np.unique(data, return_counts=True)
        # treat outliers as individual clusters
        if labels_[0] == outlier_index:
            n_outlier_ = counts_[0]
            row[outlier_i:outlier_i + n_outlier_] = 1
            outlier_i += n_outlier_
            labels_ = labels_[1:]
            counts_ = counts_[1:]
        # normalise over number of frames
        row[labels_] = counts_
        row /= len(data)            
    return utils.discrete_js_matrix(frames_per_cl)


def prepare_dres(ensemble, subspace, start=None, stop=None, step=None,
                 frames=None,
                 n_components=3, estimate_error=False,
                 n_bootstrap_samples=10,
                 similarity_matrix=None, metric=rmsd, **kwargs):
    if isinstance(subspace, type):
        try:
            # scikit-learn class
            subspace = subspace(n_components=n_components, **kwargs)
        except TypeError:
            # MDAnalysis analysis class
            subspace = subspace(ensemble, n_components=n_components, **kwargs)

    frames = prepare_frames(ensemble, frames=frames, start=start, stop=stop,
                            step=step, estimate_error=estimate_error,
                            n_bootstrap_samples=n_bootstrap_samples)

    if n_components is None:
        _, b, c = frames[0].shape[0]
        n_components = b*c

    data = np.zeros((len(frames), frames[0].shape[0], n_components))

    if isinstance(subspace, PCA):
        for i, fr in enumerate(frames):
            subspace.run(frames=fr)
            data[i] = subspace.transform(subspace._atoms,
                                         n_components=n_components,
                                         frames=fr)
    # in case it's scikit learn
    elif hasattr(subspace, "fit_transform"):
        coordinates = ensemble.timeseries(frames=frames[0], order="fac")
        coordinates = coordinates.reshape((coordinates.shape[0], -1))
        for i, fr in enumerate(frames):
            data[i] = subspace.fit_transform(coordinates)
    else:
        for i, fr in enumerate(frames):
            data[i] = subspace[fr]

    if n_components is not None:
        data = data[:, :, :n_components]

    return data


def dres(ensemble, subspace=PCA, select=None,
         n_resample=1000, seed=None, estimate_error=False, **kwargs):
    if select is not None:
        ensemble = ensemble.select_atoms(select)

    subspaces = prepare_dres(ensemble, subspace,
                             estimate_error=estimate_error, **kwargs)

    results = np.zeros((len(subspaces), ensemble.n_universes, ensemble.n_universes))
    for i, space in enumerate(subspaces):
        data = [x.T for x in ensemble.split_array(space)]
        results[i] = _dres(data, n_resample=n_resample, seed=seed)
    
    if estimate_error:
        return results.mean(axis=0), results.std(axis=0)

    return results[0]

def _dres(data, n_resample=1000, seed=None):
    n_u = len(data)
    kdes = [gaussian_kde(x) for x in data]
    resamples = [k.resample(n_resample, seed=seed) for k in kdes]
    pdfs = np.zeros((n_u, n_u, n_resample))
    for i, k in enumerate(kdes):
        for j, x in enumerate(resamples):
            pdfs[i, j] = k.evaluate(x)

    logpdfs = np.log(pdfs).mean(axis=-1)
    logmean = np.broadcast_to(np.diag(logpdfs), (n_u, n_u))
    dix = np.diag_indices(n_u)

    sum_pdfs = pdfs + pdfs[dix]
    ln_pq_exp_pq = np.log(0.5*sum_pdfs).mean(axis=-1)
    return 0.5 * (logmean + logmean.T - ln_pq_exp_pq - ln_pq_exp_pq.T)


def hes(ensemble, select=None, weights="mass", estimator="shrinkage",
        align=False, estimate_error=False, n_bootstrap_samples=10):
    if select is not None:
        ensemble = ensemble.select_atoms(select)
    # check given arguments
    if estimator == "shrinkage":
        estimator = utils.shrinkage_covariance
    elif estimator == "maximum_likelihood":
        estimator = utils.max_likelihood_covariance
    else:
        if not callable(estimator):
            raise ValueError("estimator must be 'shrinkage', "
                             "'maximum_likelihood' or be callable.")
    
    n_u = ensemble.n_universes
    n_a = ensemble.n_atoms
    n_a3 = n_a * 3

    if weights is None:
        weights = np.ones((ensemble.n_universes, ensemble.n_atoms))
    elif weights == "mass":
        weights = [ag.masses for ag in ensemble._ags]
    else:
        n_w = len(weights)
        if not len(weights) == n_u:
            if len(weights) == n_a:
                weights = [weights] * n_u
            else:
                raise ValueError("Weights must be provided for every "
                                 f"Universe. Given {n_w} weights for "
                                 f"{n_u} universes.")
    
    # make weight matrix
    weights_ = np.repeat(np.array(weights), 3, axis=1)
    n_f, n_w = weights_.shape
    weights3 = np.zeros((n_f, n_w, n_w))
    dix = np.diag_indices(n_w)
    weights3[:, dix, dix] = weights_ ** 0.5

    ensemble.transfer_to_memory()
    if align:
        AlignTraj(ensemble, ensemble, weights=weights[0],
                  in_memory=True).run()

    frames = [u.trajectory.timeseries(ag, order="fac")
              for ag, u in zip(ensemble._ags, ensemble.universes)]

    if estimate_error:
        bs = [utils.get_bootstrap_frames(f, n_samples=n_bootstrap_samples)
              for f in frames]
        frames = [np.array(list(x)) for x in zip(*bs)]
    
    else:
        frames = [frames]

    n_s = len(frames)
    avgs = np.zeros((n_s, n_u, n_u, n_a3), dtype=np.float64)
    covs = np.zeros((n_s, n_u, n_u, n_a3, n_a3), dtype=np.float64)
    
    for s in range(n_s):
        for i, (coords, w) in enumerate(zip(frames[s], weights3)):
            avgs[s, i] = coords.mean(axis=0).flatten()
            cov = estimator(coords.reshape(len(coords), -1))
            try:
                cov = np.dot(w, np.dot(cov, w))
            except ValueError:
                raise ValueError("weights dimensions don't match selected atoms")
            covs[s, i] = cov
    
    inv_covs = np.zeros((n_s, n_u, n_u, n_a3, n_a3))
    for i, sub in enumerate(covs):
        for j, arr in enumerate(sub):
            inv_covs[i, j] = np.linalg.pinv(arr[0])

    diff = avgs - avgs.transpose((0, 2, 1, 3))
    cov_prod = covs @ inv_covs.transpose((0, 2, 1, 3, 4))
    cov_prod += cov_prod.transpose((0, 2, 1, 3, 4))
    trace = np.trace(cov_prod, axis1=-1, axis2=-2) - (2 * n_a3)
    inv_cov_ = inv_covs + inv_covs.transpose((0, 2, 1, 3, 4))
    prod = np.einsum('ijklm,ijkl->ijkl', inv_cov_, diff)
    similarity = np.einsum('ijkl,ijkl->ijk', diff, prod)
    similarity = 0.25 * (similarity + trace)
    if estimate_error:
        return similarity.mean(axis=0), similarity.std(axis=0)
    return similarity[0]

def gen_window_indices(universe, window_size=10):
    n_frames = len(universe.trajectory)
    frame_ends = np.arange(window_size, n_frames+1, window_size)
    frame_ends[-1] = n_frames
    return [np.arange(i) for i in frame_ends]


def convergence(universe, data, func, window_size=10, **kwargs):
    windows = gen_window_indices(universe, window_size)
    data = data[np.concatenate(windows)]
    n_frames = windows[-1][-1] + 1
    traj_frames = [i*n_frames + x for i, x in enumerate(windows)]
    ensemble = Ensemble([universe] * len(windows),
                        frames=np.concatenate(traj_frames))
    return func(ensemble, data, **kwargs)


def ces_convergence(universe, clusters=methods.AffinityPropagation,
                    window_size=10, **kwargs):
    clusters = prepare_ces(universe, clusters, estimate_error=False)[0]
    return convergence(universe, clusters, ces, window_size=window_size,
                       **kwargs)[-1]


def dres_convergence(universe, subspace=PCA, window_size=10, seed=None,
                     **kwargs):
    subspace = prepare_dres(universe, subspace, estimate_error=False,
                            **kwargs)[0]
    return convergence(universe, subspace, dres, window_size=window_size,
                       seed=seed, **kwargs)[-1]
