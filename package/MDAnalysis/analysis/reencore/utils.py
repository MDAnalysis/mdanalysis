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

import numpy as np
import scipy

EPS = 1e-15
xlogy = np.vectorize(
    lambda x, y: 0.0 if (x <= EPS and y <= EPS) else x * np.log(y))


def discrete_kl_div(a, b):
    return np.sum(xlogy(a, a / b))

def discrete_js_div(a, b):
    ab = (a+b) * 0.5
    mask = ab > 0
    a = a[mask]
    b = b[mask]
    ab = ab[mask]
    return 0.5 * (discrete_kl_div(a, ab) + discrete_kl_div(b, ab))


def discrete_kl_matrix(matrix):
    log = np.log(matrix)
    cross_entropy = -(matrix @ log.T)
    entropy = -np.einsum('ij,ij->i', matrix, log)
    kl = cross_entropy - entropy[..., None]
    return kl


def discrete_js_matrix(matrix):
    x = matrix.shape[0]
    matrix = np.array([matrix] * x)
    matrixT = matrix.transpose((1, 0, 2))
    half = (matrix + matrixT) * 0.5
    pB = matrix / half
    loghalf = np.log(pB)
    kl = (matrix * loghalf)
    mask = (half < EPS) | (matrix < EPS) | (pB < EPS)
    kl[mask] = 0
    klsum = kl.sum(axis=-1)
    return 0.5 * (klsum + klsum.T)


def max_likelihood_covariance(coordinates, center=True):
    if center:
        coordinates -= coordinates.mean(axis=0)
    cov = np.einsum('ij,ik->jk', coordinates, coordinates)
    cov /= coordinates.shape[0]
    return cov


def shrinkage_covariance(coordinates, shrinkage_parameter=None):
    offset = coordinates - coordinates.mean(axis=0)
    t, n = coordinates.shape
    xmkt = offset.mean(axis=1)[..., None]
    coords = np.hstack([offset, xmkt])
    sample = max_likelihood_covariance(coords, False)
    sample *= (t-1)/float(t)

    # Split covariance matrix into components
    covmkt = sample[:n, n]
    varmkt = sample[n, n]
    sample = sample[:n, :n]

    covouter = np.outer(covmkt, covmkt)
    prior = covouter / varmkt
    dix = np.diag_indices(n)
    prior[dix] = sample[dix]

    if shrinkage_parameter is None:
        c = np.linalg.norm(sample - prior, ord="fro") ** 2
        inv_t = 1/t
        y = offset ** 2
        p = inv_t * (y.T @ y).sum() - (sample ** 2).sum()
        rdiag = inv_t * (y ** 2).sum() - np.trace(sample ** 2)
        z = offset * np.repeat(xmkt, n, axis=1)
        covmkt_ = np.repeat(covmkt[..., None], n, axis=1)
        v1 = inv_t * (y.T @ z) - covmkt_ * sample
        roff1 = ((v1 * covmkt_.T).sum() - np.trace(v1*covmkt)) / varmkt
        v3 = inv_t * (z.T @ z) - (varmkt * sample)
        roff3 = ((v3 * covouter).sum() - np.trace(v3 * covmkt ** 2))
        roff3 /= varmkt ** 2
        roff = 2 * roff1 - roff3
        r = rdiag + roff
        k = (p - r) / c
        shrinkage_parameter = max(0, min(1, k * inv_t))
        
    cov = shrinkage_parameter * prior + (1 - shrinkage_parameter) * sample
    return cov


def get_bootstrap_frames(frames, n_samples=100):
    n = len(frames)
    indices = [np.random.randint(0, high=n, size=n) for i in range(n_samples)]
    arr = np.array([frames[x] for x in indices])
    return arr