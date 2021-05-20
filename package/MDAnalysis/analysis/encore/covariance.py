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
"""
Covariance calculation --- :mod:`encore.covariance`
=====================================================================

The module contains functions to estimate the covariance matrix of
an ensemble of structures.

:Author: Matteo Tiberti, Wouter Boomsma, Tone Bengtsen

.. versionadded:: 0.16.0
"""
import numpy as np

def ml_covariance_estimator(coordinates, reference_coordinates=None):
    """
    Standard maximum likelihood estimator of the covariance matrix.

    Parameters
    ----------

    coordinates : numpy.array
        Flattened array of coordiantes

    reference_coordinates : numpy.array
        Optional reference to use instead of mean

    Returns
    -------

    cov_mat : numpy.array
        Estimate of  covariance matrix

    """

    if reference_coordinates is not None:

        # Offset from reference
        coordinates_offset = coordinates - reference_coordinates

    else:
        # Normal covariance calculation: distance to the average
        coordinates_offset = coordinates - np.average(coordinates, axis=0)

    # Calculate covariance manually
    coordinates_cov = np.zeros((coordinates.shape[1],
                                coordinates.shape[1]))
    for frame in coordinates_offset:
        coordinates_cov += np.outer(frame, frame)
    coordinates_cov /= coordinates.shape[0]

    return coordinates_cov

def shrinkage_covariance_estimator( coordinates,
                                    reference_coordinates=None,
                                    shrinkage_parameter=None):
    """
    Shrinkage estimator of the covariance matrix using the method described in

    Improved Estimation of the Covariance Matrix of Stock Returns With an
    Application to Portfolio Selection. Ledoit, O.; Wolf, M., Journal of
    Empirical Finance, 10, 5, 2003

    This implementation is based on the matlab code made available by Olivier
    Ledoit on his website:
    http://www.ledoit.net/ole2_abstract.htm

    Parameters
    ----------

    coordinates : numpy.array
        Flattened array of coordinates

    reference_coordinates: numpy.array
        Optional reference to use instead of mean

    shrinkage_parameter: None or float
        Optional shrinkage parameter

    Returns
    --------

    cov_mat : nump.array
        Covariance matrix
    """

    x = coordinates
    t = x.shape[0]
    n = x.shape[1]

    mean_x = np.average(x, axis=0)

    # Use provided coordinates as "mean" if provided
    if reference_coordinates is not None:
        mean_x = reference_coordinates

    x = x - mean_x
    xmkt = np.average(x, axis=1)

    # Call maximum likelihood estimator (note the additional column)
    sample = ml_covariance_estimator(np.hstack([x, xmkt[:, np.newaxis]]), 0)\
        * (t-1)/float(t)

    # Split covariance matrix into components
    covmkt = sample[0:n, n]
    varmkt = sample[n, n]
    sample = sample[:n, :n]

    # Prior
    prior = np.outer(covmkt, covmkt)/varmkt
    prior[np.ma.make_mask(np.eye(n))] = np.diag(sample)

    # If shrinkage parameter is not set, estimate it
    if shrinkage_parameter is None:

        # Frobenius norm
        c = np.linalg.norm(sample - prior, ord='fro')**2

        y = x**2
        p = 1/float(t)*np.sum(np.dot(np.transpose(y), y))\
            - np.sum(np.sum(sample**2))
        rdiag = 1/float(t)*np.sum(np.sum(y**2))\
            - np.sum(np.diag(sample)**2)
        z = x * np.repeat(xmkt[:, np.newaxis], n, axis=1)
        v1 = 1/float(t) * np.dot(np.transpose(y), z) \
            - np.repeat(covmkt[:, np.newaxis], n, axis=1)*sample
        roff1 = (np.sum(
            v1*np.transpose(
                np.repeat(
                    covmkt[:, np.newaxis], n, axis=1)
                )
            )/varmkt -
                 np.sum(np.diag(v1)*covmkt)/varmkt)
        v3 = 1/float(t)*np.dot(np.transpose(z), z) - varmkt*sample
        roff3 = (np.sum(v3*np.outer(covmkt, covmkt))/varmkt**2 -
                 np.sum(np.diag(v3)*covmkt**2)/varmkt**2)
        roff = 2*roff1-roff3
        r = rdiag+roff

        # Shrinkage constant
        k = (p-r)/c
        shrinkage_parameter = max(0, min(1, k/float(t)))

    # calculate covariance matrix
    sigma = shrinkage_parameter*prior+(1-shrinkage_parameter)*sample

    return sigma




def covariance_matrix(ensemble,
                      select="name CA",
                      estimator=shrinkage_covariance_estimator,
                      weights='mass',
                      reference=None):
    """
    Calculates (optionally mass weighted) covariance matrix

    Parameters
    ----------
    ensemble : Universe object
        The structural ensemble
    select : str (optional)
        Atom selection string in the MDAnalysis format.
    estimator : function (optional)
        Function that estimates the covariance matrix. It requires at least
        a "coordinates" numpy array (of shape (N,M,3), where N is the number
        of frames and M the number of atoms). See ml_covariance_estimator and
        shrinkage_covariance_estimator for reference.
    weights : str/array_like (optional)
        specify weights. If ``'mass'`` then chose masses of ensemble atoms, if ``None`` chose uniform weights
    reference : MDAnalysis.Universe object (optional)
        Use the distances to a specific reference structure rather than the
        distance to the mean.

    Returns
    -------
    cov_mat : numpy.array
        Covariance matrix

    """
    # Extract coordinates from ensemble
    coordinates = ensemble.trajectory.timeseries(
        ensemble.select_atoms(select),
        order='fac')

    # Flatten coordinate matrix into n_frame x n_coordinates
    coordinates = np.reshape(coordinates, (coordinates.shape[0], -1))

    # Extract coordinates from reference structure, if specified
    reference_coordinates = None
    if reference:
        # Select the same atoms in reference structure
        reference_atom_selection = reference.select_atoms(
            ensemble.get_atom_selection_string())
        reference_coordinates = reference_atom_selection.atoms.coordinates()

        # Flatten reference coordinates
        reference_coordinates = reference_coordinates.flatten()

    sigma = estimator(coordinates, reference_coordinates)

    # Optionally correct with weights
    if weights is not None:
        # Calculate mass-weighted covariance matrix
        if not isinstance(weights, (list, tuple, np.ndarray)) and weights == 'mass':
            if select:
                weights = ensemble.select_atoms(select).masses
            else:
                weights = ensemble.atoms.masses
        else:
            if select:
                req_len = ensemble.select_atoms(select).n_atoms
            else:
                req_len = ensemble.atoms.n_atoms
            if req_len != len(weights):
                raise ValueError("number of weights is unequal to number of "
                                 "atoms in ensemble")

        # broadcast to a (len(weights), 3) array
        weights = np.repeat(weights, 3)

        weight_matrix = np.sqrt(np.identity(len(weights))*weights)
        sigma = np.dot(weight_matrix, np.dot(sigma, weight_matrix))

    return sigma
