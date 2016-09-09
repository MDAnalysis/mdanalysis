# similarity.py --- Similarity measures between protein ensembles
# Copyright (C) 2014 Wouter Boomsma, Matteo Tiberti
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Ensemble Similarity Calculations --- :mod:`MDAnalysis.analysis.encore.similarity`
=================================================================================

:Author: Matteo Tiberti, Wouter Boomsma, Tone Bengtsen
:Year: 2015-2016
:Copyright: GNU Public License v3
:Maintainer: Matteo Tiberti <matteo.tiberti@gmail.com>, mtiberti on github

.. versionadded:: 0.16.0

The module contains implementations of similarity measures between protein
ensembles described in [Lindorff-Larsen2009]_. The implementation and examples
are described in [Tiberti2015]_.

The module includes facilities for handling ensembles and trajectories through
the :class:`Ensemble` class, performing clustering or dimensionality reduction
of the ensemble space, estimating multivariate probability distributions from
the input data, and more. ENCORE can be used to compare experimental and
simulation-derived ensembles, as well as estimate the convergence of
trajectories from time-dependent simulations.

ENCORE includes three different methods for calculations of similarity measures
between ensembles implemented in individual functions as well as a class to
handle the ensembles:


    + **Harmonic Ensemble Similarity** : :func:`hes`
    + **Clustering Ensemble Similarity** : :func:`ces`
    + **Dimensional Reduction Ensemble Similarity** : :func:`dres`
    + **Ensemble class** : :class:`Ensemble`

When using this module in published work please cite [Tiberti2015]_.

References
----------

    .. [Lindorff-Larsen2009] Similarity Measures for Protein Ensembles. Lindorff-Larsen, K. Ferkinghoff-Borg, J. PLoS ONE 2008, 4, e4203.

    .. [Tiberti2015] ENCORE: Software for Quantitative Ensemble Comparison. Matteo Tiberti, Elena Papaleo, Tone Bengtsen, Wouter Boomsma, Kresten Lindorff- Larsen. PLoS Comput Biol. 2015, 11

.. _Examples:
Examples
--------

The examples show how to use ENCORE to calculate a similarity measurement
of two simple ensembles. The ensembles are obtained from the MDAnalysis
test suite for two different simulations of the protein AdK. To run the
examples first execute: ::

    >>> import MDAnalysis.analysis.encore as encore
    >>> from MDAnalysis.tests.datafiles import PDB_small, DCD, DCD2

To calculate the Harmonic Ensemble Similarity (:func:`hes`)
two ensemble objects are first created and then used for calculation: ::

    >>> ens1 = encore.Ensemble(topology=PDB_small, trajectory=DCD)
    >>> ens2 = encore.Ensemble(topology=PDB_small, trajectory=DCD2)
    >>> print encore.hes([ens1, ens2])
    (array([[        0.        ,  13946090.57640726],
            [ 13946090.57640726,         0.        ]]), None)

Here None is returned in the array as the default details parameter is False. 
HES can assume any non-negative value, i.e. no upper bound exists and the
measurement can therefore be used as an absolute scale.

The calculation of the Clustering Ensemble Similarity (:func:`ces`)
is computationally more expensive. It is based on the Affinity Propagation
clustering algorithm that in turns requires a similarity matrix between
the frames the ensembles are made of (By default we use -RMSD; therefore
a full RMSD matrix between each pairs of elements needs to be computed.)
To decrease the computational load the :class:`Ensemble` object can be 
initialized by only loading every nth frame from the trajectory using the 
parameter `frame_interval`. Additionally, by saving the calculated 
 matrix using the `save_matrix` parameter, the computational cost
can be reduced for future calculations using e.g. different parameters
for the clustering algorithm, or can be reused for DRES: ::

    >>> ens1 = encore.Ensemble( topology = PDB_small, 
                                trajectory = DCD, 
                                frame_interval=3 )
    >>> ens2 = encore.Ensemble( topology = PDB_small, 
                                trajectory = DCD2, 
                                frame_interval=3)
    >>> print encore.ces([ens1, ens2], save_matrix = "minusrmsd.npz")
    (array([[ 0.        ,  0.08093055],
           [ 0.08093055,  0.        ]]), None)


In the above example the negative RMSD-matrix was saved as minusrmsd.npz and
can now be used as an input in further calculations of the
Dimensional Reduction Ensemble Similarity (:func:`dres`).
DRES is based on the estimation of the probability density in 
a dimensionally-reduced conformational space of the ensembles, obtained from 
the original space using the Stochastic proximity embedding algorithm.
As SPE requires the distance matrix calculated on the original space, we 
can reuse the previously-calculated -RMSD matrix with sign changed.
In the following example the dimensions are reduced to 3: ::

    >>> print encore.dres(  [ens1, ens2], 
                            dimensions = 3, 
                            load_matrix = "minusrmsd.npz", 
                            change_sign = True )
    (array([[ 0.        ,  0.68108127],
           [ 0.68108127,  0.        ]]), None)

Due to the stocastic nature of SPE, two
identical ensembles will not necessarily result in excatly 0 estimate of 
the similarity, but will be very close. For the same reason, calculating the 
similarity with the :func:`dres` twice will not result in 
necessarily identical values. 

It should be noted that both in :func:`ces` and :func:`dres` 
the similarity is evaluated using the Jensen-Shannon 
divergence resulting in an upper bound of ln(2), which indicates no similarity
between the ensembles and a lower bound of 0.0 signifying two identical 
ensembles. Therefore using CES and DRES ensembles can be compared in a more
relative sense respect to HES, i.e. they can be used to understand whether
ensemble A is closer to ensemble B respect to C, but absolute 
values are less meaningful as they also depend on the chosen parameters.


Functions
---------

.. autofunction:: hes

.. autofunction:: ces

.. autofunction:: dres



"""
from __future__ import print_function
import MDAnalysis as mda
import numpy as np
import warnings
import logging
from time import sleep
try:
    from scipy.stats import gaussian_kde
except ImportError:
    gaussian_kde = None
    msg = "scipy.stats.gaussian_kde could not be imported. " \
          "Dimensionality reduction ensemble comparisons will not " \
          "be available."
    warnings.warn(msg,
                  category=ImportWarning)
    logging.warn(msg)
    del msg

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
from .covariance import covariance_matrix, ml_covariance_estimator, shrinkage_covariance_estimator
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
def harmonic_ensemble_similarity(sigma1=None,
                                 sigma2=None,
                                 x1=None,
                                 x2=None):
    """
    Calculate the harmonic ensemble similarity measure
    as defined in 

        Similarity Measures for Protein Ensembles. Lindorff-Larsen, K.;
        Ferkinghoff-Borg, J. PLoS ONE 2009, 4, e4203.

    Parameters
    ----------

    sigma1 : numpy.array
        Covariance matrix for the first ensemble. If this None, calculate
        it from ensemble1 using covariance_estimator

    sigma2 : numpy.array
        Covariance matrix for the second ensemble. If this None, calculate
        it from ensemble1 using covariance_estimator

    x1: numpy.array
        Mean for the estimated normal multivariate distribution of the first
        ensemble. If this is None, calculate it from ensemble1

    x2: numpy.array
        Mean for the estimated normal multivariate distribution of the first
        ensemble.. If this is None, calculate it from ensemble2

    mass_weighted : bool
        Whether to perform mass-weighted covariance matrix estimation

    covariance_estimator : function
        Covariance estimator to be used

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
                                   selection="name CA"):
    """Clustering ensemble similarity: calculate the probability densities from
     the clusters and calculate discrete Jensen-Shannon divergence.

    Parameters
    ----------

    cc : encore.ClustersCollection
        Collection from cluster calculated by a clustering algorithm
        (e.g. Affinity propagation)

    ens1 : :class:`~MDAnalysis.core.AtomGroup.Universe`
        First ensemble to be used in comparison

    ens2 : :class:`~MDAnalysis.core.AtomGroup.Universe`
        Second ensemble to be used in comparison

    ens1_id : int
        First ensemble id as detailed in the ClustersCollection metadata

    ens2_id : int
        Second ensemble id as detailed in the ClustersCollection metadata

    selection : str
        Atom selection string in the MDAnalysis format. Default is "name CA"

    Returns
    -------

    djs : float
        Jensen-Shannon divergence between the two ensembles, as calculated by
        the clustering ensemble similarity method
    """
    ens1_coordinates = ens1.trajectory.timeseries(ens1.select_atoms(selection),
                                                  format='fac')
    ens2_coordinates = ens2.trajectory.timeseries(ens2.select_atoms(selection),
                                                  format='fac')
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
    """ Calculate clustering ensemble similarity between joined ensembles.
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

    nesensembles : int
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

    if gaussian_kde is None:
        # hack: if we are running with minimal dependencies then scipy was
        # not imported and we have to bail here (see scipy import at top)
        raise ImportError("For Kernel Density Estimation functionality you"
                          "need to import scipy")

    for i in range(1, nensembles + 1):
        this_embedded = embedded_space.transpose()[
            np.where(np.array(ensemble_assignment) == i)].transpose()
        embedded_ensembles.append(this_embedded)
        kdes.append(gaussian_kde(
            this_embedded))

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
    """
    Calculate the Jensen-Shannon divergence according the the
    Dimensionality reduction method. In this case, we have continuous
    probability densities we have to integrate over the measureable space.
    Our target is calculating Kullback-Liebler, which is defined as:

    .. math::
        D_{KL}(P(x) || Q(x)) = \\int_{-\\infty}^{\\infty}P(x_i) ln(P(x_i)/Q(x_i)) = \\langle{}ln(P(x))\\rangle{}_P - \\langle{}ln(Q(x))\\rangle{}_P

    where the :math:`\\langle{}.\\rangle{}_P` denotes an expectation calculated
    under the distribution P. We can, thus, just estimate the expectation
    values of the components to get an estimate of dKL.
    Since the Jensen-Shannon distance is actually  more complex, we need to
    estimate four expectation values:

    .. math::
         \\langle{}log(P(x))\\rangle{}_P

         \\langle{}log(Q(x))\\rangle{}_Q

         \\langle{}log(0.5*(P(x)+Q(x)))\\rangle{}_P

         \\langle{}log(0.5*(P(x)+Q(x)))\\rangle{}_Q

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
        Use this value for :math:`\\langle{}log(P(x))\\rangle{}_P; if None,
        calculate it instead

    ln_P2_exp_P2 : float or None
        Use this value for :math:`\\langle{}log(Q(x))\\rangle{}_Q`; if
        None, calculate it instead

    ln_P1P2_exp_P1 : float or None
        Use this value for
        :math:`\\langle{}log(0.5*(P(x)+Q(x)))\\rangle{}_P`;
         if None, calculate it instead

    ln_P1P2_exp_P1 : float or None
        Use this value for
        :math:`\\langle{}log(0.5*(P(x)+Q(x)))\\rangle{}_Q`;
        if None, calculate it instead

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
                            nsamples=None, ens_id_min=1, ens_id_max=None):
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
        Maximum ID of the ensemble to be considered; see description

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

    if gaussian_kde is None:
        # hack: if we are running with minimal dependencies then scipy was
        # not imported and we have to bail here (see scipy import at top)
        raise ImportError("For Kernel Density Estimation functionality you"
                          "need to import scipy")

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
        kdes.append(
            gaussian_kde(this_embedded))

    # Set number of samples
    if not nsamples:
        nsamples = this_embedded.shape[1] * 10

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
        the matrix will be just printed on screen

    header : str
            Line to be written just before the matrix

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


def bootstrap_coordinates(coords, times):
    """
    Bootstrap conformations in a :class:`~MDAnalysis.core.AtomGroup.Universe`.
    This means drawing from the encore.Ensemble.coordinates numpy array with
    replacement "times" times and returning the outcome.

    Parameters
    ----------

    coords : numpy.array
         3-dimensional coordinates array

    times : int
        Number of times the coordinates will be bootstrapped

    Returns
    -------

    out : list
        Bootstrapped coordinates list. len(out) = times.
        """
    out = []
    for t in range(times):
        this_coords = np.zeros(coords.shape)
        for c in range(this_coords.shape[0]):
            this_coords[c, :, :] = \
                coords[np.random.randint(low=0,
                                            high=this_coords.shape[0]),
                       :,
                       :]
        out.append(this_coords)
    return out



def prepare_ensembles_for_convergence_increasing_window(ensemble,
                                                        window_size,
                                                        selection="name CA"):
    """
    Generate ensembles to be fed to ces_convergence or dres_convergence
    from a single ensemble. Basically, the different slices the algorithm
    needs are generated here.

    Parameters
    ----------

    ensemble : :class:`~MDAnalysis.core.AtomGroup.Universe` object
        Input ensemble

    window_size : int
        size of the window (in number of frames) to be used

    selection : str
        Atom selection string in the MDAnalysis format. Default is "name CA"

    Returns
    -------

    tmp_ensembles : 
        the original ensemble is divided into ensembles, each being
        a window_size-long slice of the original ensemble. The last
        ensemble will be bigger if the length of the input ensemble
        is not exactly divisible by window_size.

    """

    ens_size = ensemble.trajectory.timeseries(ensemble.select_atoms(selection),
                                              format='fac').shape[0]

    rest_slices = ens_size / window_size
    residuals = ens_size % window_size
    slices_n = [0]

    tmp_ensembles = []

    for rs in range(rest_slices - 1):
        slices_n.append(slices_n[-1] + window_size)
    slices_n.append(slices_n[-1] + residuals + window_size)

    for s,sl in enumerate(slices_n[:-1]):
        tmp_ensembles.append(mda.Universe(
            ensemble.filename,
            ensemble.trajectory.timeseries()
            [:, slices_n[s]:slices_n[s + 1], :],
            format=MemoryReader))

    return tmp_ensembles


def hes(ensembles,
        selection="name CA",
        cov_estimator="shrinkage",
        mass_weighted=True,
        align=False,
        details=False,
        estimate_error=False,
        bootstrapping_samples=100,
        calc_diagonal=False):
    """

    Calculates the Harmonic Ensemble Similarity (HES) between ensembles using
    the symmetrized version of Kullback-Leibler divergence as described
    in [Lindorff-Larsen2009]_.

    Parameters
    ----------

    ensembles : list
        List of universe objects for similarity measurements.

    selection : str
        Atom selection string in the MDAnalysis format. Default is "name CA"

    cov_estimator : str, optional
        Covariance matrix estimator method, either shrinkage, `shrinkage`,
        or Maximum Likelyhood, `ml`. Default is shrinkage.

    mass_weighted : bool, optional
        Whether to perform mass-weighted covariance matrix estimation
        (default is True).

    align : bool, optional
        Whether to align the ensembles before calculating their similarity.
        Note: this changes the ensembles in-place, and will thus leave your
        ensembles in an altered state.
        (default is False)

    details : bool, optional
        Save the mean and covariance matrix for each
        ensemble in a numpy array (default is False).

    estimate_error : bool, optional
        Whether to perform error estimation (default is False).

    bootstrapping_samples : int, optional
        Number of times the similarity matrix will be bootstrapped (default
        is 100).


    Returns
    -------

    numpy.array (bidimensional)
        Harmonic similarity measurements between each pair of ensembles.

    Notes
    -----

    The method assumes that each ensemble is derived from a multivariate normal
    distribution. The mean and covariance matrix are, thus, estimatated from
    the distribution of each ensemble and used for comparision by the
    symmetrized version of Kullback-Leibler divergence defined as:

    .. math::
        D_{KL}(P(x) || Q(x)) = \\int_{-\\infty}^{\\infty}P(x_i)
        ln(P(x_i)/Q(x_i)) = \\langle{}ln(P(x))\\rangle{}_P -
        \\langle{}ln(Q(x))\\rangle{}_P


    where the :math:`\\langle{}.\\rangle{}_P` denotes an expectation
    calculated under the distribution P.

    For each ensemble, the  mean conformation is estimated as the average over
    the ensemble, and the covariance matrix is calculated by default using a
    shrinkage estimate method (or by a maximum-likelihood method, optionally).

    In the Harmonic Ensemble Similarity measurement no upper bound exists and
    the measurement can therefore best be used for relative comparison between
    multiple ensembles.

    When using this similarity measure, consider whether you want to align
    the ensembles first (see example below)

    Example
    -------

    To calculate the Harmonic Ensemble similarity, two ensembles are created
    as Universe object from a topology file and two trajectories. The
    topology- and trajectory files used are obtained from the MDAnalysis
    test suite for two different simulations of the protein AdK. To run the
    examples see the module `Examples`_ for how to import the files: ::

    >>> ens1 = Universe(PDB_small, DCD)
    >>> ens2 = Universe(PDB_small, DCD2)
    >>> print encore.hes([ens1, ens2])
    (array([[        0.        ,  13946090.57640726],
           [ 13946090.57640726,         0.        ]]), None)


    Here None is returned in the array as no details has been requested.

    You can use the align=True option to align the ensembles first. This will
    align everything to the current timestep in the first ensemble. Note that
    this changes the ens1 and ens2 objects:

    >>> print encore.hes([ens1, ens2], align=True)
    (array([[    0.        ,  6868.27953491],
           [ 6868.27953491,     0.        ]]), None)
    Alternatively, for greater flexibility in how the alignment should be done
    you can call the rms_fit_trj function manually:

    >>> align.rms_fit_trj(ens1, ens1, select="name CA", in_memory=True)
    >>> align.rms_fit_trj(ens2, ens1, select="name CA", in_memory=True)
    >>> print encore.hes([ens1, ens2])
    (array([[    0.        ,  6935.99303895],
           [ 6935.99303895,     0.        ]]), None)
    """

    # Ensure in-memory trajectories either by calling align
    # with in_memory=True or by directly calling transfer_to_memory
    # on the universe.
    if align:
        for ensemble in ensembles:
            mda.analysis.align.rms_fit_trj(ensemble, ensembles[0],
                                           select=selection,
                                           mass_weighted=True,
                                           in_memory=True)
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

    if calc_diagonal:
        pairs_indices = list(trm_indices_diag(out_matrix_eln))
    else:
        pairs_indices = list(trm_indices_nodiag(out_matrix_eln))
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
                        e_orig.select_atoms(selection),
                        format=('fac')),
                    axis=0).flatten())
                sigmas.append(covariance_matrix(ensembles_list[i][t],
                                                mass_weighted=True,
                                                estimator=covariance_estimator,
                                                selection=selection))

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

    for e in ensembles:

        # Extract coordinates from each ensemble
        coordinates_system = e.trajectory.timeseries(e.select_atoms(selection),
                                                     format='fac')

        # Average coordinates in each system
        xs.append(np.average(coordinates_system, axis=0).flatten())

        # Covariance matrices in each system
        sigmas.append(covariance_matrix(e,
                                        mass_weighted=mass_weighted,
                                        estimator=covariance_estimator,
                                        selection=selection))

    for i, j in pairs_indices:
        value = harmonic_ensemble_similarity(x1=xs[i],
                                             x2=xs[j],
                                             sigma1=sigmas[i],
                                             sigma2=sigmas[j])
        values[i, j] = value
        values[j, i] = value

    # Save details as required
    if details:
        kwds = {}
        for i in range(out_matrix_eln):
            kwds['ensemble{0:d}_mean'.format(i + 1)] = xs[i]
            kwds['ensemble{0:d}_covariance_matrix'.format(i + 1)] = sigmas[i]
        details = np.array(kwds)

    else:
        details = None

    return values, details


def ces(ensembles,
        selection="name CA",
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
    [Lindorff-Larsen2009]_.


    Parameters
    ----------

    ensembles : list
        List of ensemble objects for similarity measurements

    selection : str
        Atom selection string in the MDAnalysis format. Default is "name CA"

    clustering_method :
        A single or a list of instances of the ClusteringMethod classes from
        the clustering module. Different parameters for the same clustering
        method can be explored by adding different instances of the same
        clustering class.

    distance_matrix : encore.utils.TriangularMatrix
        distance matrix for affinity propagation. If this parameter
        is not supplied the matrix will be calculated on the fly.

    estimate_error :  bool, optional
        Whether to perform error estimation (default is False).
        Only bootstrapping mode is supported.

    bootstrapping_samples : int
        number of samples to be used for estimating error.

    ncores : int, optional
        Maximum number of cores to be used (default is 1).

    calc_diagonal : bool
        Whether to calculate the diagonal of the similarity scores
        (i.e. the simlarities of every ensemble against itself).
        If this is False (default), 0.0 will be used instead.

    allow_collapsed_result: bool
        Whether a return value of a list of one value should be collapsed
        into just the value.



    Returns
    -------

    ces, details : numpy.array, numpy.array
        ces contains the similarity values, arranged in a numpy.array.
        if one similarity value is provided as a floating point number,
        the output will be a 2-dimensional square symmetrical numpy.array.
        the order of the matrix elements depends on the order of the input
        ensemble: for instance, if

            ensemble = [ens1, ens2, ens3]

        the matrix elements [0,2] and [2,0] will contain the similarity values
        between ensembles ens1 and ens3.
        If similarity values are supplied as a list, the array will be 3-d
        with the first two dimensions running over the ensembles and
        the third dimension running over the values of the preferences
        parameter.
        Elaborating on the previous example, if preference_values are provided
        as [-1.0, -2.0] the output will be a (3,3,2) array, with element [0,2]
        corresponding to the similarity values between ens1 and ens2, and
        consisting of a 1-d array with similarity values ordered according to
        the preference_values parameters. This means that [0,2,0] will
        correspond to the similarity score between ens1 and ens3, using -1.0
        as the preference value.


    Notes
    -----

    In the Jensen-Shannon divergence the upper bound of ln(2) signifies
    no similarity between the two ensembles, the lower bound, 0.0,
    signifies identical ensembles.

    To calculate the CES, the affinity propagation method are used
    for clustering to partition the whole space of conformations in to clusters
    of structures. After the structures are clustered, the population of each
    ensemble in each cluster as a probability distribution of conformations are
    calculated. The obtained probability distribution are then compared using
    the Jensen-Shannon divergence measure between probability distributions.


    Example
    -------
    To calculate the Clustering Ensemble similarity, two ensembles are
    created as Universe object using a topology file and two trajectories. The
    topology- and trajectory files used are obtained from the MDAnalysis
    test suite for two different simulations of the protein AdK. To run the
    examples see the module `Examples`_ for how to import the files.
    Here the simplest case of just two :class:`Ensemble`s used for comparison
    are illustrated: ::

        >>> ens1 = Universe(PDB_small, DCD)
        >>> ens2 = Universe(PDB_small, DCD2)
        >>> CES = encore.ces([ens1,ens2])
        >>> print CES
            (array([[[ 0.          0.55392484]
                     [ 0.55392484  0.        ]]],None)

    Here None is returned in the array as no details has been requested.

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
                                              selection=selection,
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
                for i in range(len(ensembles_list)):
                    ensembles[-1].append(ensembles_list[i][j])
        else:
            # if all methods accept distances matrices, duplicate
            # ensemble so that it matches size of distance matrices
            # (no need to resample them since they will not be used)
            ensembles = [ensembles]*bootstrapping_samples


    # Call clustering procedure
    ccs = cluster(ensembles,
                  method= clustering_methods,
                  selection=selection,
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
                                                       selection=selection)
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
                                                   selection=selection)
                values[-1][pair[0], pair[1]] = this_val
                values[-1][pair[1], pair[0]] = this_val

    details['clustering'] = ccs

    if allow_collapsed_result and not hasattr(clustering_method, '__iter__'):
        values = values[0]

    return values, details


def dres(ensembles,
         selection="name CA",
         dimensionality_reduction_method = StochasticProximityEmbeddingNative(
             dimension=3,
             distance_cutoff = 1.5,
             min_lam=0.1,
             max_lam=2.0,
             ncycle=100,
             nstep=10000
         ),
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
    [Lindorff-Larsen2009]_.


    Parameters
    ----------

    ensembles : list
        List of ensemble objects for similarity measurements

    selection : str
        Atom selection string in the MDAnalysis format. Default is "name CA"

    dimensionality_reduction_method :
        A single or a list of instances of the DimensionalityReductionMethod
        classes from the dimensionality_reduction module. Different parameters
        for the same method can be explored by adding different instances of
        the same dimensionality reduction class.

    distance_matrix : encore.utils.TriangularMatrix
        conformational distance matrix

    nsamples : int, optional
        Number of samples to be drawn from the ensembles (default is 1000).
        Parameter used in Kernel Density Estimates (KDE) from embedded
        spaces.

    estimate_error : bool, optional
        Whether to perform error estimation (default is False)

    bootstrapping_samples : int
        number of samples to be used for estimating error.

    ncores : int, optional
        Maximum number of cores to be used (default is 1).

    calc_diagonal : bool
        Whether to calculate the diagonal of the similarity scores
        (i.e. the simlarities of every ensemble against itself).
        If this is False (default), 0.0 will be used instead.

    allow_collapsed_result: bool
        Whether a return value of a list of one value should be collapsed
        into just the value.


    Returns
    -------

    dres, details : numpy.array, numpy.array
        dres contains the similarity values, arranged in numpy.array.
        if one number of dimensions is provided as an integer,
        the output will be a 2-dimensional square symmetrical numpy.array.
        the order of the matrix elements depends on the order of the 
        input ensemble: for instance, if

            ensemble = [ens1, ens2, ens3]

        then the matrix elements [0,2] and [2,0] will contain the similarity
        values between ensembles ens1 and ens3.
        If numbers of dimensions are supplied as a list, the array will be
        3-dimensional with the first two dimensions running over the ensembles
        and the third dimension running over the number of dimensions.
        Elaborating on the previous example, if dimensions are provided
        as [2, 3] the output will be a (3,3,2) array, with element [0,2]
        corresponding to the similarity values between ens1 and ens2, and
        consisting of a 1-d array with similarity values ordered according to
        the dimensions parameters. This means that [0,2,0] will correspond to
        the similarity score between ens1 and ens3, using 2 as the number
        of dimensions.

    Notes
    -----

    To calculate the similarity the method first projects the ensembles into
    lower dimensions by using the Stochastic Proximity Embedding algorithm. A
    gaussian kernel-based density estimation method is then used to estimate
    the probability density for each ensemble which is then used to estimate
    the Jensen-shannon divergence between each pair of ensembles.

    In the Jensen-Shannon divergence the upper bound of ln(2) signifies
    no similarity between the two ensembles, the lower bound, 0.0,
    signifies identical ensembles. However, due to the stocastic nature of 
    the dimensional reduction in :func:`dres`, two identical ensembles will 
    not necessarily result in an exact 0.0 estimate of the similarity but 
    will be very close. For the same reason, calculating the similarity with
    the :func:`dres` twice will not result in two identical numbers but 
    instead small differences.  

    Example
    -------

    To calculate the Dimensional Reduction Ensemble similarity, two ensembles
    are created as Universe objects from a topology file and two trajectories. The
    topology- and trajectory files used are obtained from the MDAnalysis
    test suite for two different simulations of the protein AdK. To run the
    examples see the module `Examples`_ for how to import the files.
    Here the simplest case of comparing just two :class:`Ensemble`s are
    illustrated: ::


        >>> ens1 = Universe(PDB_small,DCD)
        >>> ens2 = Universe(PDB_small,DCD2)
        >>> DRES = encore.dres([ens1,ens2])
        >>> print DRES
           (array( [[[ 0.          0.67383396]
                 [ 0.67383396  0.        ]], None]

    Here None is returned in the array as no details has been requested.

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
                                              selection=selection,
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
        selection=selection,
        distance_matrix = distance_matrix,
        ncores = ncores,
        allow_collapsed_result = False)

    details = {}
    details["reduced_coordinates"] = coordinates
    details["dimensionality_reduction_details"] = details

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
                    selection="name CA",
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
    keeps on resampling the same ares of the conformational space, and
    therefore of convergence.


    Parameters
    ----------

    original_ensemble : :class:`~MDAnalysis.core.AtomGroup.Universe` object
        ensemble containing the trajectory whose convergence has to estimated

    window_size : int
        Size of window to be used, in number of frames

    selection : str
        Atom selection string in the MDAnalysis format. Default is "name CA"

    clustering_method :
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
    """

    ensembles = prepare_ensembles_for_convergence_increasing_window(
        original_ensemble, window_size, selection=selection)

    ccs = cluster(ensembles,
                  selection=selection,
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
                     selection="name CA",
                     dimensionality_reduction_method=StochasticProximityEmbeddingNative(
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

    original_ensemble : :class:`~MDAnalysis.core.AtomGroup.Universe` object
        ensemble containing the trajectory whose convergence has to estimated

    window_size : int
        Size of window to be used, in number of frames

    selection : str
        Atom selection string in the MDAnalysis format. Default is "name CA"

    dimensionality_reduction_method :
        A single or a list of instances of the DimensionalityReductionMethod
        classes from the dimensionality_reduction module. Different parameters
        for the same method can be explored by adding different instances of
        the same dimensionality reduction class.

    nsamples : int, optional
        Number of samples to be drawn from the ensembles (default is 1000).
        Parameter used in Kernel Density Estimates (KDE) from embedded
        spaces.

    ncores  : int, optional
        Maximum number of cores to be used (default is 1).


    Returns
    -------

    out : np.array
        array of shape (number_of_frames / window_size, preference_values).

    """
    
    ensembles = prepare_ensembles_for_convergence_increasing_window(
        original_ensemble, window_size, selection=selection)

    coordinates, dimred_details = \
        reduce_dimensionality(
            ensembles,
            selection=selection,
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
