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
=================================================================================
Ensemble Similarity Calculations --- :mod:`MDAnalysis.analysis.encore.similarity`
=================================================================================

:Author: Matteo Tiberti, Wouter Boomsma, Tone Bengtsen

.. versionadded:: 0.16.0

.. deprecated:: 2.8.0
   This module is deprecated in favour of the 
   MDAKit `mdaencore <https://mdanalysis.org/mdaencore/>`_ and will be removed
   in MDAnalysis 3.0.0.

The module contains implementations of similarity measures between protein
ensembles described in :footcite:p:`LindorffLarsen2009`. The implementation and
examples are described in :footcite:p:`Tiberti2015`.

The module includes facilities for handling ensembles and trajectories through
the :class:`Universe` class, performing clustering or dimensionality reduction
of the ensemble space, estimating multivariate probability distributions from
the input data, and more. ENCORE can be used to compare experimental and
simulation-derived ensembles, as well as estimate the convergence of
trajectories from time-dependent simulations.

ENCORE includes three different methods for calculations of similarity measures
between ensembles implemented in individual functions:

+ **Harmonic Ensemble Similarity** : :func:`hes`
+ **Clustering Ensemble Similarity** : :func:`ces`
+ **Dimensional Reduction Ensemble Similarity** : :func:`dres`

as well as two methods to evaluate the convergence of trajectories:

+ **Clustering based convergence evaluation** : :func:`ces_convergence`
+ **Dimensionality-reduction based convergence evaluation** : :func:`dres_convergence`

When using this module in published work please cite :footcite:p:`Tiberti2015`.

.. rubric:: References

.. footbibliography::

.. _Examples:
Examples
========

The examples show how to use ENCORE to calculate a similarity measurement
of two simple ensembles. The ensembles are obtained from the MDAnalysis
test suite for two different simulations of the protein AdK.

To calculate the Harmonic Ensemble Similarity (:func:`hes`)
two ensemble objects are first created and then used for calculation:

    >>> from MDAnalysis import Universe
    >>> import MDAnalysis.analysis.encore as encore
    >>> from MDAnalysis.tests.datafiles import PSF, DCD, DCD2
    >>> ens1 = Universe(PSF, DCD)
    >>> ens2 = Universe(PSF, DCD2)
    >>> HES, details = encore.hes([ens1, ens2])
    >>> print(HES)
    [[       0.         38279540.04524205]
     [38279540.04524205        0.        ]]

HES can assume any non-negative value, i.e. no upper bound exists and the
measurement can therefore be used as an absolute scale.

The calculation of the Clustering Ensemble Similarity (:func:`ces`)
is computationally more expensive. It is based on clustering algorithms that in
turn require a similarity matrix between the frames the ensembles are made
of. The similarity matrix is derived from a distance matrix (By default a RMSD
matrix; a full RMSD matrix between each pairs of elements needs to be computed).
The RMSD matrix is automatically calculated:

    >>> from MDAnalysis import Universe
    >>> import MDAnalysis.analysis.encore as encore
    >>> from MDAnalysis.tests.datafiles import PSF, DCD, DCD2
    >>> ens1 = Universe(PSF, DCD)
    >>> ens2 = Universe(PSF, DCD2)
    >>> CES, details = encore.ces([ens1, ens2])
    >>> print(CES)
    [[0.         0.68070702]
     [0.68070702 0.        ]]

The RMSD matrix can also be separately calculated to reuse it, e.g. for running
CES with different parameters or running the
Dimensional Reduction Ensemble Similarity (:func:`dres`) method.
DRES is based on the estimation of the probability density in
a dimensionally-reduced conformational space of the ensembles, obtained from
the original space using either the Stochastic Proximity Embedding algorithm or
the Principal Component Analysis.
In the following example the dimensions are reduced to 3 using the
RMSD matrix and the default SPE dimensional reduction method:

    >>> from MDAnalysis import Universe
    >>> import MDAnalysis.analysis.encore as encore
    >>> from MDAnalysis.tests.datafiles import PSF, DCD, DCD2
    >>> ens1 = Universe(PSF, DCD)
    >>> ens2 = Universe(PSF, DCD2)
    >>> rmsd_matrix = encore.get_distance_matrix(
    ...                             encore.utils.merge_universes([ens1, ens2]))
    >>> DRES,details = encore.dres([ens1, ens2],
    ...                             distance_matrix = rmsd_matrix)

The RMSD matrix can also be saved on disk with the option ``save_matrix``: ::

    rmsd_matrix = encore.get_distance_matrix(
                                    encore.utils.merge_universes([ens1, ens2]),
                                    save_matrix="rmsd.npz")

It can then be loaded and reused at a later time instead of being recalculated: ::

    rmsd_matrix = encore.get_distance_matrix(
                                    encore.utils.merge_universes([ens1, ens2]),
                                    load_matrix="rmsd.npz")

In addition to the quantitative similarity estimate, the dimensional reduction
can easily be visualized, see the ``Example`` section in
:mod:`MDAnalysis.analysis.encore.dimensionality_reduction.reduce_dimensionality`.
Due to the stochastic nature of SPE, two identical ensembles will not
necessarily result in an exactly 0 estimate of the similarity, but will be very
close. For the same reason, calculating the similarity with the :func:`dres`
twice will not result in necessarily identical values but rather two very close
values.

It should be noted that both in :func:`ces` and :func:`dres` the similarity is
evaluated using the Jensen-Shannon divergence resulting in an upper bound of
ln(2), which indicates no similarity between the ensembles and a lower bound
of 0.0 signifying two identical ensembles. In contrast, the :func:`hes` function uses
a symmetrized version of the Kullback-Leibler divergence, which is unbounded.


Functions for ensemble comparisons
==================================

.. autofunction:: hes
   :noindex:

.. autofunction:: ces
   :noindex:

.. autofunction:: dres
   :noindex:

Function reference
==================

.. All functions are included via automodule :members:.

"""
import warnings
import logging

import numpy as np
import scipy.stats

import MDAnalysis as mda

from ...coordinates.memory import MemoryReader
from .confdistmatrix import get_distance_matrix
from .bootstrap import (get_distance_matrix_bootstrap_samples,
                        get_ensemble_bootstrap_samples)
from .clustering.cluster import cluster
from .clustering.ClusteringMethod import AffinityPropagationNative
from .dimensionality_reduction.DimensionalityReductionMethod import (
    StochasticProximityEmbeddingNative)
from .dimensionality_reduction.reduce_dimensionality import (
    reduce_dimensionality)
from .covariance import (
    covariance_matrix, ml_covariance_estimator, shrinkage_covariance_estimator)
from .utils import merge_universes
from .utils import trm_indices_diag, trm_indices_nodiag

# Low boundary value for log() argument - ensure no nans
EPSILON = 1E-15

xlogy = np.vectorize(
    lambda x, y: 0.0 if (x <= EPSILON and y <= EPSILON) else x * np.log(y))


def discrete_kullback_leibler_divergence(pA, pB):
    """Kullback-Leibler divergence between discrete probability distribution.
    Notice that since this measure is not symmetric ::
    :math:`d_{KL}(p_A,p_B) != d_{KL}(p_B,p_A)`

    Parameters
    ----------

    pA : iterable of floats
        First discrete probability density function

    pB : iterable of floats
        Second discrete probability density function

    Returns
    -------

    dkl : float
        Discrete Kullback-Liebler divergence
    """

    return np.sum(xlogy(pA, pA / pB))


# discrete dJS
def discrete_jensen_shannon_divergence(pA, pB):
    """Jensen-Shannon divergence between discrete probability distributions.

    Parameters
    ----------

    pA : iterable of floats
        First discrete probability density function

    pB : iterable of floats
        Second discrete probability density function

    Returns
    -------

    djs : float
        Discrete Jensen-Shannon divergence
"""
    return 0.5 * (discrete_kullback_leibler_divergence(pA, (pA + pB) * 0.5) +
                  discrete_kullback_leibler_divergence(pB, (pA + pB) * 0.5))


# calculate harmonic similarity
def harmonic_ensemble_similarity(sigma1,
                                 sigma2,
                                 x1,
                                 x2):
    """
    Calculate the harmonic ensemble similarity measure
    as defined in :footcite:p:`Tiberti2015`.

    Parameters
    ----------

    sigma1 : numpy.array
        Covariance matrix for the first ensemble.

    sigma2 : numpy.array
        Covariance matrix for the second ensemble.

    x1: numpy.array
        Mean for the estimated normal multivariate distribution of the first
        ensemble.

    x2: numpy.array
        Mean for the estimated normal multivariate distribution of the second
        ensemble.

    Returns
    -------

        dhes : float
            harmonic similarity measure
    """

    # Inverse covariance matrices
    sigma1_inv = np.linalg.pinv(sigma1)
    sigma2_inv = np.linalg.pinv(sigma2)

    # Difference between average vectors
    d_avg = x1 - x2

    # Distance measure
    trace = np.trace(np.dot(sigma1, sigma2_inv) +
                        np.dot(sigma2, sigma1_inv)
                        - 2 * np.identity(sigma1.shape[0]))

    d_hes = 0.25 * (np.dot(np.transpose(d_avg),
                              np.dot(sigma1_inv + sigma2_inv,
                                        d_avg)) + trace)
    return d_hes


def clustering_ensemble_similarity(cc, ens1, ens1_id, ens2, ens2_id,
                                   select="name CA"):
    """Clustering ensemble similarity: calculate the probability densities from
     the clusters and calculate discrete Jensen-Shannon divergence.

    Parameters
    ----------

    cc : encore.clustering.ClustersCollection
        Collection from cluster calculated by a clustering algorithm
        (e.g. Affinity propagation)

    ens1 : :class:`~MDAnalysis.core.universe.Universe`
        First ensemble to be used in comparison

    ens1_id : int
        First ensemble id as detailed in the ClustersCollection metadata

    ens2 : :class:`~MDAnalysis.core.universe.Universe`
        Second ensemble to be used in comparison

    ens2_id : int
        Second ensemble id as detailed in the ClustersCollection metadata

    select : str
        Atom selection string in the MDAnalysis format. Default is "name CA".

    Returns
    -------

    djs : float
        Jensen-Shannon divergence between the two ensembles, as calculated by
        the clustering ensemble similarity method
    """
    ens1_coordinates = ens1.trajectory.timeseries(ens1.select_atoms(select),
                                                  order='fac')
    ens2_coordinates = ens2.trajectory.timeseries(ens2.select_atoms(select),
                                                  order='fac')
    tmpA = np.array([np.where(c.metadata['ensemble_membership'] == ens1_id)[
                            0].shape[0] / float(ens1_coordinates.shape[0]) for
                        c in cc])
    tmpB = np.array([np.where(c.metadata['ensemble_membership'] == ens2_id)[
                            0].shape[0] / float(ens2_coordinates.shape[0]) for
                        c in cc])

    # Exclude clusters which have 0 elements in both ensembles
    pA = tmpA[tmpA + tmpB > EPSILON]
    pB = tmpB[tmpA + tmpB > EPSILON]

    return discrete_jensen_shannon_divergence(pA, pB)


def cumulative_clustering_ensemble_similarity(cc, ens1_id, ens2_id,
                                              ens1_id_min=1, ens2_id_min=1):
    """
    Calculate clustering ensemble similarity between joined ensembles.
    This means that, after clustering has been performed, some ensembles are
    merged and the dJS is calculated between the probability distributions of
    the two clusters groups. In particular, the two ensemble groups are defined
    by their ensembles id: one of the two joined ensembles will comprise all
    the ensembles with id [ens1_id_min, ens1_id], and the other ensembles will
    comprise all the ensembles with id [ens2_id_min, ens2_id].

    Parameters
    ----------

    cc : encore.ClustersCollection
            Collection from cluster calculated by a clustering algorithm
            (e.g. Affinity propagation)

    ens1_id : int
            First ensemble id as detailed in the ClustersCollection
            metadata

    ens2_id : int
            Second ensemble id as detailed in the ClustersCollection
            metadata

    Returns
    -------

    djs : float
            Jensen-Shannon divergence between the two ensembles, as
            calculated by the clustering ensemble similarity method

"""

    ensA = [np.where(np.logical_and(
        c.metadata['ensemble_membership'] <= ens1_id,
        c.metadata['ensemble_membership'])
                     >= ens1_id_min)[0].shape[0] for c in cc]
    ensB = [np.where(np.logical_and(
        c.metadata['ensemble_membership'] <= ens2_id,
        c.metadata['ensemble_membership'])
                     >= ens2_id_min)[0].shape[0] for c in cc]
    sizeA = float(np.sum(ensA))
    sizeB = float(np.sum(ensB))

    tmpA = np.array(ensA) / sizeA
    tmpB = np.array(ensB) / sizeB

    # Exclude clusters which have 0 elements in both ensembles
    pA = tmpA[tmpA + tmpB > EPSILON]
    pB = tmpB[tmpA + tmpB > EPSILON]

    return discrete_jensen_shannon_divergence(pA, pB)


def gen_kde_pdfs(embedded_space, ensemble_assignment, nensembles,
                 nsamples):
    """
    Generate Kernel Density Estimates (KDE) from embedded spaces and
    elaborate the coordinates for later use.

    Parameters
    ----------

    embedded_space : numpy.array
        Array containing the coordinates of the embedded space

    ensemble_assignment : numpy.array
        Array containing one int per ensemble conformation. These allow to
        distinguish, in the complete embedded space, which conformations
        belong to each ensemble. For instance if ensemble_assignment
        is [1,1,1,1,2,2], it means that the first four conformations belong
        to ensemble 1 and the last two to ensemble 2

    nensembles : int
        Number of ensembles

    nsamples : int
        samples to be drawn from the ensembles. Will be required in
        a later stage in order to calculate dJS.

    Returns
    -------

    kdes : scipy.stats.gaussian_kde
        KDEs calculated from ensembles

    resamples : list of numpy.array
        For each KDE, draw samples according to the probability distribution
        of the KDE mixture model

    embedded_ensembles : list of numpy.array
        List of numpy.array containing, each one, the elements of the
        embedded space belonging to a certain ensemble
    """
    kdes = []
    embedded_ensembles = []
    resamples = []

    for i in range(1, nensembles + 1):
        this_embedded = embedded_space.transpose()[
            np.where(np.array(ensemble_assignment) == i)].transpose()
        embedded_ensembles.append(this_embedded)
        kdes.append(scipy.stats.gaussian_kde(this_embedded))

    # # Set number of samples
    # if not nsamples:
    #     nsamples = this_embedded.shape[1] * 10

    # Resample according to probability distributions
    for this_kde in kdes:
        resamples.append(this_kde.resample(nsamples))

    return (kdes, resamples, embedded_ensembles)


def dimred_ensemble_similarity(kde1, resamples1, kde2, resamples2,
                               ln_P1_exp_P1=None, ln_P2_exp_P2=None,
                               ln_P1P2_exp_P1=None, ln_P1P2_exp_P2=None):
    r"""Calculate the Jensen-Shannon divergence according the Dimensionality
    reduction method.

    In this case, we have continuous probability densities, this we need to
    integrate over the measurable space. The aim is to first calculate the
    Kullback-Liebler divergence, which is defined as:

    .. math::

       D_{KL}(P(x) || Q(x)) =
           \int_{-\infty}^{\infty}P(x_i) ln(P(x_i)/Q(x_i)) =
           \langle{}ln(P(x))\rangle{}_P - \langle{}ln(Q(x))\rangle{}_P

    where the :math:`\langle{}.\rangle{}_P` denotes an expectation calculated
    under the distribution P. We can, thus, just estimate the expectation
    values of the components to get an estimate of dKL.  Since the
    Jensen-Shannon distance is actually more complex, we need to estimate four
    expectation values:

    .. math::
         \langle{}log(P(x))\rangle{}_P

         \langle{}log(Q(x))\rangle{}_Q

         \langle{}log(0.5*(P(x)+Q(x)))\rangle{}_P

         \langle{}log(0.5*(P(x)+Q(x)))\rangle{}_Q

    Parameters
    ----------

    kde1 : scipy.stats.gaussian_kde
        Kernel density estimation for ensemble 1

    resamples1 : numpy.array
        Samples drawn according do kde1. Will be used as samples to
        calculate the expected values according to 'P' as detailed before.

    kde2 : scipy.stats.gaussian_kde
            Kernel density estimation for ensemble 2

    resamples2 : numpy.array
        Samples drawn according do kde2. Will be used as sample to
        calculate the expected values according to 'Q' as detailed before.

    ln_P1_exp_P1 : float or None
        Use this value for :math:`\langle{}log(P(x))\rangle{}_P`; if ``None``,
        calculate it instead

    ln_P2_exp_P2 : float or None
        Use this value for :math:`\langle{}log(Q(x))\rangle{}_Q`; if
        ``None``, calculate it instead

    ln_P1P2_exp_P1 : float or None
        Use this value for
        :math:`\langle{}log(0.5*(P(x)+Q(x)))\rangle{}_P`;
        if ``None``, calculate it instead

    ln_P1P2_exp_P2 : float or None
        Use this value for
        :math:`\langle{}log(0.5*(P(x)+Q(x)))\rangle{}_Q`;
        if ``None``, calculate it instead

    Returns
    -------
    djs : float
        Jensen-Shannon divergence calculated according to the dimensionality
        reduction method

    """

    if not ln_P1_exp_P1 and not ln_P2_exp_P2 and not ln_P1P2_exp_P1 and not \
            ln_P1P2_exp_P2:
        ln_P1_exp_P1 = np.average(np.log(kde1.evaluate(resamples1)))
        ln_P2_exp_P2 = np.average(np.log(kde2.evaluate(resamples2)))
        ln_P1P2_exp_P1 = np.average(np.log(
            0.5 * (kde1.evaluate(resamples1) + kde2.evaluate(resamples1))))
        ln_P1P2_exp_P2 = np.average(np.log(
            0.5 * (kde1.evaluate(resamples2) + kde2.evaluate(resamples2))))

    return 0.5 * (
        ln_P1_exp_P1 - ln_P1P2_exp_P1 + ln_P2_exp_P2 - ln_P1P2_exp_P2)


def cumulative_gen_kde_pdfs(embedded_space, ensemble_assignment, nensembles,
                            nsamples, ens_id_min=1, ens_id_max=None):
    """
    Generate Kernel Density Estimates (KDE) from embedded spaces and
    elaborate the coordinates for later use. However, consider more than
    one ensemble as the space on which the KDE will be generated. In
    particular, will use ensembles with ID [ens_id_min, ens_id_max].


    Parameters
    ----------

    embedded_space : numpy.array
            Array containing the coordinates of the embedded space

    ensemble_assignment : numpy.array
            array containing one int per ensemble conformation. These allow
            to distinguish, in the complete embedded space, which
            conformations belong to each ensemble. For instance if
            ensemble_assignment is [1,1,1,1,2,2], it means that the first
            four conformations belong to ensemble 1 and the last two
            to ensemble 2

    nensembles : int
            Number of ensembles

    nsamples : int
        Samples to be drawn from the ensembles. Will be required in a later
        stage in order to calculate dJS.

    ens_id_min : int
        Minimum ID of the ensemble to be considered; see description

    ens_id_max : int
        Maximum ID of the ensemble to be considered; see description. If None,
        it will be set to the maximum possible value given the number of
        ensembles.

    Returns
    -------

    kdes : scipy.stats.gaussian_kde
            KDEs calculated from ensembles

    resamples : list of numpy.array
            For each KDE, draw samples according to the probability
            distribution of the kde mixture model

    embedded_ensembles : list of numpy.array
            List of numpy.array containing, each one, the elements of the
            embedded space belonging to a certain ensemble


    """

    kdes = []
    embedded_ensembles = []
    resamples = []
    if not ens_id_max:
        ens_id_max = nensembles + 1
    for i in range(ens_id_min, ens_id_max):
        this_embedded = embedded_space.transpose()[np.where(
            np.logical_and(ensemble_assignment >= ens_id_min,
                              ensemble_assignment <= i))].transpose()
        embedded_ensembles.append(this_embedded)
        kdes.append(scipy.stats.gaussian_kde(this_embedded))

    # Resample according to probability distributions
    for this_kde in kdes:
        resamples.append(this_kde.resample(nsamples))

    return (kdes, resamples, embedded_ensembles)


def write_output(matrix, base_fname=None, header="", suffix="",
                 extension="dat"):
    """
    Write output matrix with a nice format, to stdout and optionally a file.

    Parameters
    ----------

    matrix : encore.utils.TriangularMatrix
        Matrix containing the values to be printed

    base_fname : str
        Basic filename for output. If None, no files will be written, and
        the matrix will be just printed on standard output

    header : str
        Text to be written just before the matrix

    suffix : str
        String to be concatenated to basename, in order to get the final
        file name

    extension : str
        Extension for the output file

    """

    if base_fname is not None:
        fname = base_fname + "-" + suffix + "." + extension
    else:
        fname = None
    matrix.square_print(header=header, fname=fname)


def prepare_ensembles_for_convergence_increasing_window(ensemble,
                                                        window_size,
                                                        select="name CA"):
    """
    Generate ensembles to be fed to ces_convergence or dres_convergence
    from a single ensemble. Basically, the different slices the algorithm
    needs are generated here.

    Parameters
    ----------

    ensemble : :class:`~MDAnalysis.core.universe.Universe` object
        Input ensemble

    window_size : int
        size of the window (in number of frames) to be used

    select : str
        Atom selection string in the MDAnalysis format. Default is "name CA"

    Returns
    -------

    tmp_ensembles :
        The original ensemble is divided into different ensembles, each being
        a window_size-long slice of the original ensemble. The last
        ensemble will be bigger if the length of the input ensemble
        is not exactly divisible by window_size.

    """

    ens_size = ensemble.trajectory.timeseries(ensemble.select_atoms(select),
                                              order='fac').shape[0]

    rest_slices = ens_size // window_size
    residuals = ens_size % window_size
    slices_n = [0]

    tmp_ensembles = []

    for rs in range(rest_slices - 1):
        slices_n.append(slices_n[-1] + window_size)
    slices_n.append(slices_n[-1] + residuals + window_size)

    for s,sl in enumerate(slices_n[:-1]):
        tmp_ensembles.append(mda.Universe(
            ensemble.filename,
            ensemble.trajectory.timeseries(order='fac')
            [slices_n[s]:slices_n[s + 1], :, :],
            format=MemoryReader))

    return tmp_ensembles


def hes(ensembles,
        select="name CA",
        cov_estimator="shrinkage",
        weights='mass',
        align=False,
        estimate_error=False,
        bootstrapping_samples=100,
        calc_diagonal=False):
    r"""Calculates the Harmonic Ensemble Similarity (HES) between ensembles.

    The HES is calculated with the symmetrized version of Kullback-Leibler
    divergence as described in :footcite:p:`Tiberti2015`.

    Parameters
    ----------
    ensembles : list
        List of Universe objects for similarity measurements.
    select : str, optional
        Atom selection string in the MDAnalysis format. Default is "name CA"
    cov_estimator : str, optional
        Covariance matrix estimator method, either shrinkage, `shrinkage`,
        or Maximum Likelyhood, `ml`. Default is shrinkage.
    weights : str/array_like, optional
        specify optional weights. If ``mass`` then chose masses of ensemble atoms
    align : bool, optional
        Whether to align the ensembles before calculating their similarity.
        Note: this changes the ensembles in-place, and will thus leave your
        ensembles in an altered state.
        (default is False)
    estimate_error : bool, optional
        Whether to perform error estimation (default is False).
    bootstrapping_samples : int, optional
        Number of times the similarity matrix will be bootstrapped (default
        is 100), only if estimate_error is True.
    calc_diagonal : bool, optional
        Whether to calculate the diagonal of the similarity scores
        (i.e. the similarities of every ensemble against itself).
        If this is False (default), 0.0 will be used instead.

    Returns
    -------
    hes, details : numpy.array, dictionary
        Harmonic similarity measurements between each pair of ensembles,
        and dict containing mean and covariance matrix for each ensemble

    Notes
    -----
    The method assumes that each ensemble is derived from a multivariate normal
    distribution. The mean and covariance matrix are, thus, estimatated from
    the distribution of each ensemble and used for comparision by the
    symmetrized version of Kullback-Leibler divergence defined as:

    .. math::
       D_{KL}(P(x) || Q(x)) =
           \int_{-\infty}^{\infty}P(x_i) ln(P(x_i)/Q(x_i)) =
           \langle{}ln(P(x))\rangle{}_P - \langle{}ln(Q(x))\rangle{}_P


    where the :math:`\langle{}.\rangle{}_P` denotes an expectation
    calculated under the distribution :math:`P`.

    For each ensemble, the  mean conformation is estimated as the average over
    the ensemble, and the covariance matrix is calculated by default using a
    shrinkage estimation method (or by a maximum-likelihood method,
    optionally).

    Note that the symmetrized version of the Kullback-Leibler divergence has no
    upper bound (unlike the Jensen-Shannon divergence used by for instance CES and DRES).

    When using this similarity measure, consider whether you want to align
    the ensembles first (see example below).

    Example
    -------

    To calculate the Harmonic Ensemble similarity, two ensembles are created
    as Universe objects from a topology file and two trajectories. The
    topology- and trajectory files used are obtained from the MDAnalysis
    test suite for two different simulations of the protein AdK.
    You can use the ``align=True`` option to align the ensembles first. This will
    align everything to the current timestep in the first ensemble. Note that
    this changes the ``ens1`` and ``ens2`` objects:

        >>> from MDAnalysis import Universe
        >>> import MDAnalysis.analysis.encore as encore
        >>> from MDAnalysis.tests.datafiles import PSF, DCD, DCD2
        >>> ens1 = Universe(PSF, DCD)
        >>> ens2 = Universe(PSF, DCD2)
        >>> HES, details = encore.hes([ens1, ens2])
        >>> print(HES)
        [[       0.         38279540.04524205]
         [38279540.04524205        0.        ]]
        >>> print(encore.hes([ens1, ens2], align=True)[0])
        [[   0.         6889.89729056]
         [6889.89729056    0.        ]]

    Alternatively, for greater flexibility in how the alignment should be done
    you can call use an :class:`~MDAnalysis.analysis.align.AlignTraj` object
    manually:

        >>> from MDAnalysis import Universe
        >>> import MDAnalysis.analysis.encore as encore
        >>> from MDAnalysis.tests.datafiles import PSF, DCD, DCD2
        >>> from MDAnalysis.analysis import align
        >>> ens1 = Universe(PSF, DCD)
        >>> ens2 = Universe(PSF, DCD2)
        >>> _ = align.AlignTraj(ens1, ens1, select="name CA", in_memory=True).run()
        >>> _ = align.AlignTraj(ens2, ens1, select="name CA", in_memory=True).run()
        >>> print(encore.hes([ens1, ens2])[0])
        [[   0.         6889.89729056]
         [6889.89729056    0.        ]]


    .. versionchanged:: 1.0.0
       ``hes`` doesn't accept the `details` argument anymore, it always returns
       the details of the calculation instead, in the form of a dictionary

    """

    if not isinstance(weights, (list, tuple, np.ndarray)) and weights == 'mass':
        weights = ['mass' for _ in range(len(ensembles))]
    elif weights is not None:
        if len(weights) != len(ensembles):
            raise ValueError("need weights for every ensemble")
    else:
        weights = [None for _ in range(len(ensembles))]

    # Ensure in-memory trajectories either by calling align
    # with in_memory=True or by directly calling transfer_to_memory
    # on the universe.
    if align:
        for e, w in zip(ensembles, weights):
            mda.analysis.align.AlignTraj(e, ensembles[0],
                                         select=select,
                                         weights=w,
                                         in_memory=True).run()
    else:
        for ensemble in ensembles:
            ensemble.transfer_to_memory()

    if calc_diagonal:
        pairs_indices = list(trm_indices_diag(len(ensembles)))
    else:
        pairs_indices = list(trm_indices_nodiag(len(ensembles)))

    logging.info("Chosen metric: Harmonic similarity")
    if cov_estimator == "shrinkage":
        covariance_estimator = shrinkage_covariance_estimator
        logging.info("    Covariance matrix estimator: Shrinkage")
    elif cov_estimator == "ml":
        covariance_estimator = ml_covariance_estimator
        logging.info("    Covariance matrix estimator: Maximum Likelihood")
    else:
        logging.error(
            "Covariance estimator {0} is not supported. "
            "Choose between 'shrinkage' and 'ml'.".format(cov_estimator))
        return None

    out_matrix_eln = len(ensembles)

    xs = []
    sigmas = []

    if estimate_error:
        data = []
        ensembles_list = []
        for i, ensemble in enumerate(ensembles):
            ensembles_list.append(
                get_ensemble_bootstrap_samples(
                    ensemble,
                    samples=bootstrapping_samples))
        for t in range(bootstrapping_samples):
            logging.info("The coordinates will be bootstrapped.")

            xs = []
            sigmas = []
            values = np.zeros((out_matrix_eln, out_matrix_eln))
            for i, e_orig in enumerate(ensembles):
                xs.append(np.average(
                    ensembles_list[i][t].trajectory.timeseries(
                        e_orig.select_atoms(select),
                        order=('fac')),
                    axis=0).flatten())
                sigmas.append(covariance_matrix(ensembles_list[i][t],
                                                weights=weights[i],
                                                estimator=covariance_estimator,
                                                select=select))

            for pair in pairs_indices:
                value = harmonic_ensemble_similarity(x1=xs[pair[0]],
                                                     x2=xs[pair[1]],
                                                     sigma1=sigmas[pair[0]],
                                                     sigma2=sigmas[pair[1]])
                values[pair[0], pair[1]] = value
                values[pair[1], pair[0]] = value
            data.append(values)
        avgs = np.average(data, axis=0)
        stds = np.std(data, axis=0)

        return (avgs, stds)

    # Calculate the parameters for the multivariate normal distribution
    # of each ensemble
    values = np.zeros((out_matrix_eln, out_matrix_eln))

    for e, w in zip(ensembles, weights):
        # Extract coordinates from each ensemble
        coordinates_system = e.trajectory.timeseries(e.select_atoms(select),
                                                     order='fac')

        # Average coordinates in each system
        xs.append(np.average(coordinates_system, axis=0).flatten())

        # Covariance matrices in each system
        sigmas.append(covariance_matrix(e,
                                        weights=w,
                                        estimator=covariance_estimator,
                                        select=select))

    for i, j in pairs_indices:
        value = harmonic_ensemble_similarity(x1=xs[i],
                                             x2=xs[j],
                                             sigma1=sigmas[i],
                                             sigma2=sigmas[j])
        values[i, j] = value
        values[j, i] = value

    # Save details as required
    details = {}
    for i in range(out_matrix_eln):
        details['ensemble{0:d}_mean'.format(i + 1)] = xs[i]
        details['ensemble{0:d}_covariance_matrix'.format(i + 1)] = sigmas[i]

    return values, details


def ces(ensembles,
        select="name CA",
        clustering_method=AffinityPropagationNative(
            preference=-1.0,
            max_iter=500,
            convergence_iter=50,
            damping=0.9,
            add_noise=True),
        distance_matrix=None,
            estimate_error=False,
        bootstrapping_samples=10,
        ncores=1,
        calc_diagonal=False,
        allow_collapsed_result=True):
    """

    Calculates the Clustering Ensemble Similarity (CES) between ensembles
    using the Jensen-Shannon divergence as described in
    :footcite:p:`Tiberti2015`.

    Parameters
    ----------

    ensembles : list
        List of ensemble objects for similarity measurements

    select : str, optional
        Atom selection string in the MDAnalysis format. Default is "name CA"

    clustering_method :
        A single or a list of instances of the
        :class:`MDAnalysis.analysis.encore.clustering.ClusteringMethod` classes
        from the clustering module. Different parameters for the same clustering
        method can be explored by adding different instances of the same
        clustering class. Clustering methods options are the
        Affinity Propagation (default), the DBSCAN and the KMeans. The latter
        two methods need the sklearn python module installed.

    distance_matrix : encore.utils.TriangularMatrix
        Distance matrix clustering methods. If this parameter
        is not supplied the matrix will be calculated on the fly.

    estimate_error :  bool, optional
        Whether to perform error estimation (default is False).
        Only bootstrapping mode is supported.

    bootstrapping_samples : int, optional
        number of samples to be used for estimating error.

    ncores : int, optional
        Maximum number of cores to be used (default is 1).

    calc_diagonal : bool, optional
        Whether to calculate the diagonal of the similarity scores
        (i.e. the similarities of every ensemble against itself).
        If this is False (default), 0.0 will be used instead.

    allow_collapsed_result: bool, optional
        Whether a return value of a list of one value should be collapsed
        into just the value.



    Returns
    -------

    ces, details : numpy.array, numpy.array

        ces contains the similarity values, arranged in a numpy.array.
        If only one clustering_method is provided the output will be a
        2-dimensional square symmetrical numpy.array. The order of the matrix
        elements depends on the order of the input ensembles: for instance, if

            ensemble = [ens1, ens2, ens3]

        the matrix elements [0,2] and [2,0] will both contain the similarity
        value between ensembles ens1 and ens3.
        Elaborating on the previous example, if *n* ensembles are given and *m*
        clustering_methods are provided the output will be a list of *m* arrays
        ordered by the input sequence of methods, each with a *n*x*n*
        symmetrical similarity matrix.

        details contains information on the clustering: the individual size of
        each cluster, the centroids and the frames associated with each cluster.


    Notes
    -----

    In the Jensen-Shannon divergence the upper bound of ln(2) signifies
    no similarity between the two ensembles, the lower bound, 0.0,
    signifies identical ensembles.

    To calculate the CES, the affinity propagation method (or others, if
    specified) is used to partition the whole space of conformations. The
    population of each ensemble in each cluster is then taken as a probability
    density function. Different probability density functions from each
    ensemble are finally compared using the Jensen-Shannon divergence measure.

    Examples
    --------
    To calculate the Clustering Ensemble similarity, two ensembles are
    created as Universe object using a topology file and two trajectories. The
    topology- and trajectory files used are obtained from the MDAnalysis
    test suite for two different simulations of the protein AdK.
    To use a different clustering method, set the parameter clustering_method
    (Note that the sklearn module must be installed). Likewise, different parameters
    for the same clustering method can be explored by adding different
    instances of the same clustering class.
    Here the simplest case of just two instances of :class:`Universe` is illustrated:

        >>> from MDAnalysis import Universe
        >>> import MDAnalysis.analysis.encore as encore
        >>> from MDAnalysis.tests.datafiles import PSF, DCD, DCD2
        >>> ens1 = Universe(PSF, DCD)
        >>> ens2 = Universe(PSF, DCD2)
        >>> CES, details = encore.ces([ens1,ens2])
        >>> print(CES)
        [[0.         0.68070702]
         [0.68070702 0.        ]]
        >>> CES, details = encore.ces([ens1,ens2],
        ...                           clustering_method = [encore.DBSCAN(eps=0.45),
        ...                                                encore.DBSCAN(eps=0.50)])
        >>> print("eps=0.45: ", CES[0])
        eps=0.45:  [[0.         0.20447236]
         [0.20447236 0.        ]]

        >>> print("eps=0.5: ", CES[1])
        eps=0.5:  [[0.         0.25331629]
         [0.25331629 0.        ]]

    """

    for ensemble in ensembles:
        ensemble.transfer_to_memory()

    if calc_diagonal:
        pairs_indices = list(trm_indices_diag(len(ensembles)))
    else:
        pairs_indices = list(trm_indices_nodiag(len(ensembles)))

    clustering_methods = clustering_method
    if not hasattr(clustering_method, '__iter__'):
        clustering_methods = [clustering_method]

    any_method_accept_distance_matrix = \
        np.any([method.accepts_distance_matrix for method in clustering_methods])
    all_methods_accept_distance_matrix = \
        np.all([method.accepts_distance_matrix for method in clustering_methods])

    # Register which ensembles the samples belong to
    ensemble_assignment = []
    for i, ensemble in enumerate(ensembles):
        ensemble_assignment += [i+1]*len(ensemble.trajectory)

    # Calculate distance matrix if not provided
    if any_method_accept_distance_matrix and not distance_matrix:
        distance_matrix = get_distance_matrix(merge_universes(ensembles),
                                              select=select,
                                              ncores=ncores)
    if estimate_error:
        if any_method_accept_distance_matrix:
            distance_matrix = \
                get_distance_matrix_bootstrap_samples(
                    distance_matrix,
                    ensemble_assignment,
                    samples=bootstrapping_samples,
                    ncores=ncores)
        if not all_methods_accept_distance_matrix:
            ensembles_list = []
            for i, ensemble in enumerate(ensembles):
                ensembles_list.append(
                    get_ensemble_bootstrap_samples(
                        ensemble,
                        samples=bootstrapping_samples))
            ensembles = []
            for j in range(bootstrapping_samples):
                ensembles.append([])
                for i, e in enumerate(ensembles_list):
                    ensembles[-1].append(e[j])
        else:
            # if all methods accept distances matrices, duplicate
            # ensemble so that it matches size of distance matrices
            # (no need to resample them since they will not be used)
            ensembles = [ensembles]*bootstrapping_samples


    # Call clustering procedure
    ccs = cluster(ensembles,
                  method= clustering_methods,
                  select=select,
                  distance_matrix = distance_matrix,
                  ncores = ncores,
                  allow_collapsed_result=False)

    # Do error analysis
    if estimate_error:
        k = 0
        values = {}
        avgs = []
        stds = []
        for i, p in enumerate(clustering_methods):
            failed_runs = 0
            values[i] = []
            for j in range(bootstrapping_samples):
                if ccs[k].clusters is None:
                    failed_runs += 1
                    k += 1
                    continue
                values[i].append(np.zeros((len(ensembles[j]),
                                           len(ensembles[j]))))

                for pair in pairs_indices:
                    # Calculate dJS
                    this_djs = \
                        clustering_ensemble_similarity(ccs[k],
                                                       ensembles[j][
                                                           pair[0]],
                                                       pair[0] + 1,
                                                       ensembles[j][
                                                           pair[1]],
                                                       pair[1] + 1,
                                                       select=select)
                    values[i][-1][pair[0], pair[1]] = this_djs
                    values[i][-1][pair[1], pair[0]] = this_djs
                k += 1
            outs = np.array(values[i])
            avgs.append(np.average(outs, axis=0))
            stds.append(np.std(outs, axis=0))

        if hasattr(clustering_method, '__iter__'):
            pass
        else:
            avgs = avgs[0]
            stds = stds[0]

        return avgs, stds

    values = []
    details = {}
    for i, p in enumerate(clustering_methods):
        if ccs[i].clusters is None:
            continue
        else:
            values.append(np.zeros((len(ensembles), len(ensembles))))

            for pair in pairs_indices:
                # Calculate dJS
                this_val = \
                    clustering_ensemble_similarity(ccs[i],
                                                   ensembles[pair[0]],
                                                   pair[0] + 1,
                                                   ensembles[pair[1]],
                                                   pair[1] + 1,
                                                   select=select)
                values[-1][pair[0], pair[1]] = this_val
                values[-1][pair[1], pair[0]] = this_val

    details['clustering'] = ccs

    if allow_collapsed_result and not hasattr(clustering_method, '__iter__'):
        values = values[0]

    return values, details


def dres(ensembles,
         select="name CA",
         dimensionality_reduction_method = StochasticProximityEmbeddingNative(
             dimension=3,
             distance_cutoff = 1.5,
             min_lam=0.1,
             max_lam=2.0,
             ncycle=100,
             nstep=10000),
         distance_matrix=None,
         nsamples=1000,
         estimate_error=False,
         bootstrapping_samples=100,
         ncores=1,
         calc_diagonal=False,
         allow_collapsed_result=True):
    """

    Calculates the Dimensional Reduction Ensemble Similarity (DRES) between
    ensembles using the Jensen-Shannon divergence as described in
    :footcite:p:`Tiberti2015`.


    Parameters
    ----------

    ensembles : list
        List of ensemble objects for similarity measurements

    select : str, optional
        Atom selection string in the MDAnalysis format. Default is "name CA"

    dimensionality_reduction_method :
        A single or a list of instances of the DimensionalityReductionMethod
        classes from the dimensionality_reduction module. Different parameters
        for the same method can be explored by adding different instances of
        the same dimensionality reduction class. Provided methods are the
        Stochastic Proximity Embedding (default) and the Principal Component
        Analysis.

    distance_matrix : encore.utils.TriangularMatrix
        conformational distance matrix, It will be calculated on the fly
        from the ensemble data if it is not provided.

    nsamples : int, optional
        Number of samples to be drawn from the ensembles (default is 1000).
        This is used to resample the density estimates and calculate the
        Jensen-Shannon divergence between ensembles.

    estimate_error : bool, optional
        Whether to perform error estimation (default is False)

    bootstrapping_samples : int, optional
        number of samples to be used for estimating error.

    ncores : int, optional
        Maximum number of cores to be used (default is 1).

    calc_diagonal : bool, optional
        Whether to calculate the diagonal of the similarity scores
        (i.e. the simlarities of every ensemble against itself).
        If this is False (default), 0.0 will be used instead.

    allow_collapsed_result: bool, optional
        Whether a return value of a list of one value should be collapsed
        into just the value.


    Returns
    -------

    dres, details : numpy.array, numpy.array
        dres contains the similarity values, arranged in numpy.array.
        If one number of dimensions is provided as an integer,
        the output will be a 2-dimensional square symmetrical numpy.array.
        The order of the matrix elements depends on the order of the
        input ensemble: for instance, if

            ensemble = [ens1, ens2, ens3]

        then the matrix elements [0,2] and [2,0] will both contain the
        similarity value between ensembles ens1 and ens3.
        Elaborating on the previous example, if *n* ensembles are given and *m*
        methods are provided the output will be a list of *m* arrays
        ordered by the input sequence of methods, each with a *n*x*n*
        symmetrical similarity matrix.

        details provide an array of the reduced_coordinates.

    Notes
    -----

    To calculate the similarity, the method first projects the ensembles into
    lower dimensions by using the Stochastic Proximity Embedding (or others)
    algorithm. A gaussian kernel-based density estimation method is then used
    to estimate the probability density for each ensemble which is then used
    to compute the Jensen-Shannon divergence between each pair of ensembles.

    In the Jensen-Shannon divergence the upper bound of ln(2) signifies
    no similarity between the two ensembles, the lower bound, 0.0,
    signifies identical ensembles. However, due to the stochastic nature of
    the dimensional reduction in :func:`dres`, two identical ensembles will
    not necessarily result in an exact 0.0 estimate of the similarity but
    will be very close. For the same reason, calculating the similarity with
    the :func:`dres` twice will not result in two identical numbers; small
    differences have to be expected.

    Examples
    --------

    To calculate the Dimensional Reduction Ensemble similarity, two ensembles
    are created as Universe objects from a topology file and two trajectories.
    The topology- and trajectory files used are obtained from the MDAnalysis
    test suite for two different simulations of the protein AdK.
    To use a different dimensional reduction methods, simply set the
    parameter dimensionality_reduction_method. Likewise, different parameters
    for the same clustering method can be explored by adding different
    instances of the same method  class.
    Here the simplest case of comparing just two instances of :class:`Universe` is
    illustrated:

        >>> from MDAnalysis import Universe
        >>> import MDAnalysis.analysis.encore as encore
        >>> from MDAnalysis.tests.datafiles import PSF, DCD, DCD2
        >>> ens1 = Universe(PSF,DCD)
        >>> ens2 = Universe(PSF,DCD2)
        >>> DRES, details = encore.dres([ens1,ens2])
        >>> PCA_method = encore.PrincipalComponentAnalysis(dimension=2)
        >>> DRES, details = encore.dres([ens1,ens2],
        ...                             dimensionality_reduction_method=PCA_method)

    In addition to the quantitative similarity estimate, the dimensional
    reduction can easily be visualized, see the ``Example`` section in
    :mod:`MDAnalysis.analysis.encore.dimensionality_reduction.reduce_dimensionality``

    """

    for ensemble in ensembles:
        ensemble.transfer_to_memory()

    if calc_diagonal:
        pairs_indices = list(trm_indices_diag(len(ensembles)))
    else:
        pairs_indices = list(trm_indices_nodiag(len(ensembles)))

    dimensionality_reduction_methods = dimensionality_reduction_method
    if not hasattr(dimensionality_reduction_method, '__iter__'):
        dimensionality_reduction_methods = [dimensionality_reduction_method]

    any_method_accept_distance_matrix = \
        np.any([method.accepts_distance_matrix for method in dimensionality_reduction_methods])
    all_methods_accept_distance_matrix = \
        np.all([method.accepts_distance_matrix for method in dimensionality_reduction_methods])

    # Register which ensembles the samples belong to
    ensemble_assignment = []
    for i, ensemble in enumerate(ensembles):
        ensemble_assignment += [i+1]*len(ensemble.trajectory)

    # Calculate distance matrix if not provided
    if any_method_accept_distance_matrix and not distance_matrix:
        distance_matrix = get_distance_matrix(merge_universes(ensembles),
                                              select=select,
                                              ncores=ncores)
    if estimate_error:
        if any_method_accept_distance_matrix:
            distance_matrix = \
                get_distance_matrix_bootstrap_samples(
                    distance_matrix,
                    ensemble_assignment,
                    samples=bootstrapping_samples,
                    ncores=ncores)
        if not all_methods_accept_distance_matrix:
            ensembles_list = []
            for i, ensemble in enumerate(ensembles):
                ensembles_list.append(
                    get_ensemble_bootstrap_samples(
                        ensemble,
                        samples=bootstrapping_samples))
            ensembles = []
            for j in range(bootstrapping_samples):
                ensembles.append(ensembles_list[i, j] for i
                                 in range(ensembles_list.shape[0]))
        else:
            # if all methods accept distances matrices, duplicate
            # ensemble so that it matches size of distance matrices
            # (no need to resample them since they will not be used)
            ensembles = [ensembles] * bootstrapping_samples

    # Call dimensionality reduction procedure
    coordinates, dim_red_details = reduce_dimensionality(
        ensembles,
        method=dimensionality_reduction_methods,
        select=select,
        distance_matrix = distance_matrix,
        ncores = ncores,
        allow_collapsed_result = False)

    details = {}
    details["reduced_coordinates"] = coordinates
    details["dimensionality_reduction_details"] = dim_red_details

    if estimate_error:
        k = 0
        values = {}
        avgs = []
        stds = []
        for i,method in enumerate(dimensionality_reduction_methods):
            values[i] = []
            for j in range(bootstrapping_samples):

                values[i].append(np.zeros((len(ensembles[j]),
                                           len(ensembles[j]))))

                kdes, resamples, embedded_ensembles = gen_kde_pdfs(
                    coordinates[k],
                    ensemble_assignment,
                    len(ensembles[j]),
                    nsamples=nsamples)

                for pair in pairs_indices:
                    this_value = dimred_ensemble_similarity(kdes[pair[0]],
                                                            resamples[pair[0]],
                                                            kdes[pair[1]],
                                                            resamples[pair[1]])
                    values[i][-1][pair[0], pair[1]] = this_value
                    values[i][-1][pair[1], pair[0]] = this_value

                k += 1
            outs = np.array(values[i])
            avgs.append(np.average(outs, axis=0))
            stds.append(np.std(outs, axis=0))

        if hasattr(dimensionality_reduction_method, '__iter__'):
            pass
        else:
            avgs = avgs[0]
            stds = stds[0]

        return avgs, stds

    values = []

    for i,method in enumerate(dimensionality_reduction_methods):

        values.append(np.zeros((len(ensembles), len(ensembles))))
        kdes, resamples, embedded_ensembles = gen_kde_pdfs(coordinates[i],
                                                           ensemble_assignment,
                                                           len(ensembles),
                                                           nsamples=nsamples)

        for pair in pairs_indices:
            this_value = dimred_ensemble_similarity(kdes[pair[0]],
                                                    resamples[pair[0]],
                                                    kdes[pair[1]],
                                                    resamples[pair[1]])
            values[-1][pair[0], pair[1]] = this_value
            values[-1][pair[1], pair[0]] = this_value

    if allow_collapsed_result and not hasattr(dimensionality_reduction_method,
                                              '__iter__'):
        values = values[0]

    return values, details


def ces_convergence(original_ensemble,
                    window_size,
                    select="name CA",
                    clustering_method=AffinityPropagationNative(
                        preference=-1.0,
                        max_iter=500,
                        convergence_iter=50,
                        damping=0.9,
                        add_noise=True),
                    ncores=1):
    """
    Use the CES to evaluate the convergence of the ensemble/trajectory.
    CES will be calculated between the whole trajectory contained in an
    ensemble and windows of such trajectory of increasing sizes, so that
    the similarity values should gradually drop to zero. The rate at which
    the value reach zero will be indicative of how much the trajectory
    keeps on resampling the same regions of the conformational space, and
    therefore of convergence.


    Parameters
    ----------

    original_ensemble : :class:`~MDAnalysis.core.universe.Universe` object
        ensemble containing the trajectory whose convergence has to estimated

    window_size : int
        Size of window to be used, in number of frames

    select : str, optional
        Atom selection string in the MDAnalysis format. Default is "name CA"

    clustering_method : MDAnalysis.analysis.encore.clustering.ClusteringMethod
        A single or a list of instances of the ClusteringMethod classes from
        the clustering module. Different parameters for the same clustering
        method can be explored by adding different instances of the same
        clustering class.

    ncores : int, optional
        Maximum number of cores to be used (default is 1).


    Returns
    -------

    out : np.array
        array of shape (number_of_frames / window_size, preference_values).


    Example
    --------
    To calculate the convergence of a trajectory using the clustering ensemble
    similarity method a Universe object is created from a topology file and the
    trajectory. The topology- and trajectory files used are obtained from the
    MDAnalysis test suite for two different simulations of the protein AdK.
    Here the simplest case of evaluating the convergence is illustrated by
    splitting the trajectory into a window_size of 10 frames:

        >>> from MDAnalysis import Universe
        >>> import MDAnalysis.analysis.encore as encore
        >>> from MDAnalysis.tests.datafiles import PSF, DCD, DCD2
        >>> ens1 = Universe(PSF,DCD)
        >>> ces_conv = encore.ces_convergence(ens1, 10)
        >>> print(ces_conv)
        [[0.48194205]
         [0.40284672]
         [0.31699026]
         [0.25220447]
         [0.19829817]
         [0.14642725]
         [0.09911411]
         [0.05667391]
         [0.        ]]

    """

    ensembles = prepare_ensembles_for_convergence_increasing_window(
        original_ensemble, window_size, select=select)

    ccs = cluster(ensembles,
                  select=select,
                  method=clustering_method,
                  allow_collapsed_result=False,
                  ncores=ncores)

    out = []
    for cc in ccs:
        if cc.clusters is None:
            continue
        out.append(np.zeros(len(ensembles)))
        for j, ensemble in enumerate(ensembles):
            out[-1][j] = cumulative_clustering_ensemble_similarity(
                cc,
                len(ensembles),
                j + 1)

    out = np.array(out).T
    return out


def dres_convergence(original_ensemble,
                     window_size,
                     select="name CA",
                     dimensionality_reduction_method = \
                            StochasticProximityEmbeddingNative(
                                dimension=3,
                                distance_cutoff=1.5,
                                min_lam=0.1,
                                max_lam=2.0,
                                ncycle=100,
                                nstep=10000
                            ),
                     nsamples=1000,
                     ncores=1):
    """
    Use the DRES to evaluate the convergence of the ensemble/trajectory.
    DRES will be calculated between the whole trajectory contained in an
    ensemble and windows of such trajectory of increasing sizes, so that
    the similarity values should gradually drop to zero. The rate at which
    the value reach zero will be indicative of how much the trajectory
    keeps on resampling the same ares of the conformational space, and
    therefore of convergence.

    Parameters
    ----------

    original_ensemble : :class:`~MDAnalysis.core.universe.Universe` object
        ensemble containing the trajectory whose convergence has to estimated

    window_size : int
        Size of window to be used, in number of frames

    select : str, optional
        Atom selection string in the MDAnalysis format. Default is "name CA"

    dimensionality_reduction_method :
        A single or a list of instances of the DimensionalityReductionMethod
        classes from the dimensionality_reduction module. Different parameters
        for the same method can be explored by adding different instances of
        the same dimensionality reduction class.

    nsamples : int, optional
        Number of samples to be drawn from the ensembles (default is 1000).
        This is akin to the nsamples parameter of dres().

    ncores  : int, optional
        Maximum number of cores to be used (default is 1).


    Returns
    -------

    out : np.array
        array of shape (number_of_frames / window_size, preference_values).



    Example
    --------
    To calculate the convergence of a trajectory using the DRES
    method, a Universe object is created from a topology file and the
    trajectory. The topology- and trajectory files used are obtained from the
    MDAnalysis test suite for two different simulations of the protein AdK.
    Here the simplest case of evaluating the convergence is illustrated by
    splitting the trajectory into a window_size of 10 frames:

        >>> from MDAnalysis import Universe
        >>> import MDAnalysis.analysis.encore as encore
        >>> from MDAnalysis.tests.datafiles import PSF, DCD, DCD2
        >>> ens1 = Universe(PSF,DCD)
        >>> dres_conv = encore.dres_convergence(ens1, 10)

    Here, the rate at which the values reach zero will be indicative of how
    much the trajectory keeps on resampling the same ares of the conformational
    space, and therefore of convergence.
    """

    ensembles = prepare_ensembles_for_convergence_increasing_window(
        original_ensemble, window_size, select=select)

    coordinates, dimred_details = \
        reduce_dimensionality(
            ensembles,
            select=select,
            method=dimensionality_reduction_method,
            allow_collapsed_result=False,
            ncores=ncores)

    ensemble_assignment = []
    for i, ensemble in enumerate(ensembles):
        ensemble_assignment += [i+1]*len(ensemble.trajectory)
    ensemble_assignment = np.array(ensemble_assignment)

    out = []
    for i, _ in enumerate(coordinates):

        out.append(np.zeros(len(ensembles)))

        kdes, resamples, embedded_ensembles = \
            cumulative_gen_kde_pdfs(
                coordinates[i],
                ensemble_assignment=ensemble_assignment,
                nensembles=len(ensembles),
                nsamples=nsamples)

        for j, ensemble in enumerate(ensembles):
            out[-1][j] = dimred_ensemble_similarity(kdes[-1],
                                                    resamples[-1],
                                                    kdes[j],
                                                    resamples[j])

    out = np.array(out).T
    return out
