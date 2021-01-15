import itertools

import numpy as np
from scipy.stats import gaussian_kde

from ...core.ensemble import Ensemble
from ..align import AlignTraj
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

def ces(ensemble, clusters, outlier_label=-1, outlier_index=-1):
    try:
        clusters = clusters.data_labels
    except AttributeError:
        pass

    i_labels, n_cl, n_outliers = clusters_to_indices(clusters, outlier_label)
    
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


def dres(ensemble, subspace=PCA, n_components=3, n_resample=1000,
         start=None, stop=None, step=None, **kwargs):
    if isinstance(subspace, type):
        subspace = subspace(ensemble, **kwargs)
    if isinstance(subspace, PCA):
        if not subspace._calculated:
            subspace.run(start=start, stop=stop, step=step)
        subspace = subspace.transform(subspace._atoms,
                                      n_components=n_components,
                                      start=start, stop=stop,
                                      step=step)
    if n_components is not None:
        subspace = subspace[:, :n_components]

    n_u = ensemble.n_universes
    kdes = [gaussian_kde(x.T) for x in ensemble.split_array(subspace)]
    resamples = [k.resample(size=n_resample) for k in kdes]
    # resamples = np.concatenate([k.resample(size=n_resample) for k in kdes],
    #                            axis=1)  # (n_dim, n_u*n_samples)
    # pdfs = np.array([k.evaluate(resamples) for k in kdes])
    # pdfs = np.array(np.split(pdfs, n_u, axis=1)) #pdfs.reshape((n_u, n_u, n_resample))
    pdfs = np.zeros((n_u, n_u, n_resample))
    for i, k in enumerate(kdes):
        pdfs[i] = [k.evaluate(x) for x in resamples]
    logpdfs = np.log(pdfs).mean(axis=-1)
    pdfsT = pdfs.transpose((1, 0, 2))

    sum_pdfs = pdfs + pdfsT
    ln_pq_exp_pq = np.log(0.5 * (pdfs + pdfsT)).mean(axis=-1)
    print(pdfs.shape)
    print(logpdfs)

    return 0.5 * (logpdfs + logpdfs.T - ln_pq_exp_pq - ln_pq_exp_pq.T)


def hes(ensemble, weights="mass", estimator="shrinkage", align=False,
        estimate_error=False, n_samples=100):
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
    # print(np.array(weights).shape)

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
        bs = [utils.get_bootstrap_frames(f, n_samples=n_samples)
              for f in frames]
        frames = zip(bs)
    else:
        frames = [frames]

    n_s = len(frames)
    avgs = np.zeros((n_s, n_u, n_u, n_a3), dtype=np.float64)
    covs = np.zeros((n_s, n_u, n_u, n_a3, n_a3), dtype=np.float64)
    
    for s in range(n_s):
        for i, (coords, w) in enumerate(zip(frames[s], weights3)):
            avgs[s, i] = coords.mean(axis=0).flatten()
            cov = estimator(coords.reshape(len(coords), -1))
            # print(cov[:5, :5])
            try:
                cov = np.dot(w, np.dot(cov, w))
            except ValueError:
                raise ValueError("weights dimensions don't match selected atoms")
            # print(cov[:5, :5])
            covs[s, i] = cov
    
    inv_covs = np.zeros((n_s, n_u, n_u, n_a3, n_a3))
    for i, sub in enumerate(covs):
        for j, arr in enumerate(sub):
            inv_covs[i, j] = np.linalg.pinv(arr[0])
            # if j == 0:
            #     print(arr[0])
            #     print(inv_covs[i, j])

    # print(avgs.shape)
    diff = avgs - avgs.transpose((0, 2, 1, 3))

    # print(avgs[0, 0, 0][:10])
    # print(avgs[0, 1, 0][:10])
    # print((avgs[0, 0, 0] - avgs[0, 1, 0])[:10])

    # print(diff[0][0, 1][:20])

    cov_prod = covs @ inv_covs.transpose((0, 2, 1, 3, 4))
    cov_prod += cov_prod.transpose((0, 2, 1, 3, 4))
    # print(cov_prod[0, 0, 1, : 10])
    # trace = np.trace(cov_prod - 2 * np.identity(n_a3), axis1=-1, axis2=-2)
    trace = np.trace(cov_prod, axis1=-1, axis2=-2) - (2 * n_a3)
    # print(trace)

    inv_cov_ = inv_covs + inv_covs.transpose((0, 2, 1, 3, 4))
    # print(inv_cov_.dtype)
    prod = np.einsum('ijklm,ijkl->ijkl', inv_cov_, diff)
    similarity = np.einsum('ijkl,ijkl->ijk', diff, prod)
    similarity = 0.25 * (similarity + trace)
    if estimate_error:
        return similarity.mean(axis=0), similarity.std(axis=0)
    return similarity[0]

def gen_window_indices(universe, window_size=10):
    n_frames = len(universe.trajectory)
    frame_ends = np.arange(window_size, n_frames, window_size)
    frame_ends[-1] = n_frames
    return [np.arange(i) for i in frame_ends]


def convergence(universe, data, func, window_size=10, **kwargs):
    windows = gen_window_indices(universe, window_size)
    data = data[np.ravel(windows)]
    n_frames = windows[-1][-1] + 1
    traj_frames = [i*n_frames + x for i, x in enumerate(windows)]
    ensemble = Ensemble([universe] * len(windows), frames=traj_frames)
    return func(ensemble, data, **kwargs)


def ces_convergence(universe, clusters, window_size=10, **kwargs):
    try:
        clusters = clusters.cluster_indices
    except AttributeError:
        pass

    return convergence(universe, clusters, ces, window_size=window_size,
                       **kwargs)

    # i_labels, n_cl, n_outliers = clusters_to_indices(clusters, outlier_label)

    
    # ag = universe.select_atoms(select)
    # n_frames = len(ag.universe.trajectory)
    # frame_ends = np.arange(window_size, n_frames, window_size)
    # frame_ends[-1] = n_frames
    # n_windows = len(frame_ends)
    # frames_per_cl = np.zeros((n_windows, n_cl + n_outliers),
    #                          dtype=float)
    
    # for i, end in enumerate(frame_ends[::-1], 1):
    #     row = frames_per_cl[i]
    #     labels_, counts_ = np.unique(i_labels[:end], return_counts=True)
    #     if labels_[0] == outlier_index:
    #         n_outlier_ = counts_[0]
    #         row[n_cl:n_cl + n_outlier_] = 1
    #         labels_ = labels_[1:]
    #         counts_ = counts_[1:]
        
    #     # normalise over number of frames
    #     row[labels_] = counts_
    #     row /= end
    
    # return utils.discrete_js_matrix(frames_per_cl)[0][::-1]


def dres_convergence(universe, subspace, window_size=10, **kwargs):
    return convergence(universe, subspace, dres, window_size=window_size,
                       **kwargs)

    
    # if n_components is not None:
    #     subspace = subspace[:, :n_components]

    # n_frames = len(universe.trajectory)
    # frame_ends = np.arange(window_size, n_frames, window_size)
    # frame_ends[-1] = n_frames
    # n_w = len(frame_ends)


    # kdes = [gaussian_kde(subspace[:end].T) for end in frame_ends[::-1]]
    # resamples = np.concatenate([k.resample(size=n_resample) for k in kdes],
    #                            axis=1)  # (n_dim, n_u*n_samples)
    # pdfs = np.array([k.evaluate(resamples) for k in kdes])
    # pdfs = pdfs.reshape((n_w, n_resample))
    # logpdfs = np.broadcast_to(np.log(pdfs).mean(axis=1), (n_w, n_w))

    # pdfs_ = np.broadcast_to(pdfs, (n_w, n_w, n_resample))
    # sum_pdfs = pdfs_ + np.transpose(pdfs_, axes=(1, 0, 2))
    # ln_pq_exp_pq = np.log(0.5 * sum_pdfs).mean(axis=-1)

    # return 0.5 * (logpdfs + logpdfs.T - ln_pq_exp_pq - ln_pq_exp_pq.T)[0][::-1]