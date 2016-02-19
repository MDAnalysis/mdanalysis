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

    >>> from MDAnalysis import *
    >>> from MDAnalysis.analysis.encore.similarity import *
    >>> from MDAnalysis.tests.datafiles import PDB_small, DCD, DCD2


To calculate the Harmonic Ensemble Similarity (:func:`hes`)
two ensemble objects are first created and then used for calculation: ::

    >>> ens1 = Ensemble(topology=PDB_small, trajectory=DCD)
    >>> ens2 = Ensemble(topology=PDB_small, trajectory=DCD2)
    >>> print hes([ens1, ens2])
    13946090.5764

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

    >>> ens1 = Ensemble(topology = PDB_small, trajectory = DCD, frame_interval=3)
    >>> ens2 = Ensemble(topology = PDB_small, trajectory = DCD2, frame_interval=3)
    >>> print ces([ens1, ens2], save_matrix = "minusrmsd.npz")
    0.55392484    

In the above example the negative RMSD-matrix was saved as minusrmsd.npz and
can now be used as an input in further calculations of the
Dimensional Reduction Ensemble Similarity (:func:`dres`), thereby reducing the 
computational cost. DRES is based on the estimation of the probability density in 
a dimensionally-reduced conformational space of the ensembles, obtained from 
the original space using the Stochastic proximity embedding algorithm.
As SPE requires the distance matrix calculated on the original space, we 
can reuse the previously-calculated -RMSD matrix with sign changed.
In the following example the dimensions are reduced to 3: ::

    >>> print dres([ens1, ens2], dimensions = 3, load_matrix = "minusrmsd.npz", change_sign = True)
    0.648772821

Due to the stocastic nature of SPE, two
identical ensembles will not necessarily result in an exact 0.0 estimate of 
the similarity but will be very close. For the same reason, calculating the 
similarity with the :func:`dres` twice will not result in 
necessarily identical values. 

It should be noted that both in :func:`ces` and :func:`dres` 
the similarity is evaluated using the Jensen-Shannon 
divergence resulting in an upper bound of ln(2), which indicates no similarity
between the ensembles and a lower bound of 0.0 signifying two identical 
ensembles. Therefore using CES and DRES ensembles can be compared in a more relative sense
respect to HES, i.e. they can be used to understand whether 
ensemble A is closer to ensemble B respect to C, but absolute 
values are less meaningful as they also depend on the chosen parameters.


Functions
---------

.. autofunction:: hes

.. autofunction:: ces

.. autofunction:: dres



"""
import optparse
import numpy
import warnings
import logging
from time import sleep
from MDAnalysis import Universe
from .Ensemble import Ensemble
from .clustering.Cluster import ClustersCollection
from .clustering.affinityprop import AffinityPropagation
from .dimensionality_reduction.stochasticproxembed import \
    StochasticProximityEmbedding, kNNStochasticProximityEmbedding
from .confdistmatrix import MinusRMSDMatrixGenerator, RMSDMatrixGenerator
from .covariance import covariance_matrix, EstimatorShrinkage, EstimatorML
from multiprocessing import cpu_count
from .utils import *
from scipy.stats import gaussian_kde
from random import randint
import sys

# Silence deprecation warnings - scipy problem
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

# Low boundary value for log() argument - ensure no nans 
EPSILON = 1E-15

# x*log(y) with the assumption that 0*(log(0)) = 0
xlogy = numpy.vectorize(
    lambda x, y: 0.0 if (x <= EPSILON and y <= EPSILON) else x * numpy.log(y))


def is_int(n):
    try:
        int(n)
        return True
    except:
        return False


# discrete dKL
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

    return numpy.sum(xlogy(pA, pA / pB))


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
def harmonic_ensemble_similarity(ensemble1=None,
                                 ensemble2=None,
                                 sigma1=None,
                                 sigma2=None,
                                 x1=None,
                                 x2=None,
                                 mass_weighted=True,
                                 covariance_estimator=EstimatorShrinkage()):
    """
    Calculate the harmonic ensemble similarity measure
    as defined in 

        Similarity Measures for Protein Ensembles. Lindorff-Larsen, K.;
        Ferkinghoff-Borg, J. PLoS ONE 2009, 4, e4203.

    Parameters
    ----------

        ensemble1 : encore.Ensemble or None
            First ensemble to be compared. If this is None, sigma1 and x1 must be provided.

        ensemble2 : encore.Ensemble or None
            Second ensemble to be compared. If this is None, sigma2 and x2 must be provided.

        sigma1 : numpy.array
            Covariance matrix for the first ensemble. If this None, calculate it from ensemble1 using covariance_estimator

        sigma2 : numpy.array
            Covariance matrix for the second ensemble. If this None, calculate it from ensemble1 using covariance_estimator

        x1: numpy.array
            Mean for the estimated normal multivariate distribution of the first ensemble. If this is None, calculate it from ensemble1

        x2: numpy.array
                    Mean for the estimated normal multivariate distribution of the first ensemble.. If this is None, calculate it from ensemble2

        mass_weighted : bool
            Whether to perform mass-weighted covariance matrix estimation

        covariance_estimator : either EstimatorShrinkage or EstimatorML objects
            Which covariance estimator to use

    Returns
    -------

        dhes : float
            harmonic similarity measure
    """

    # If matrices and means are specified, use them
    if x1 == None or x2 == None or sigma1 == None or sigma2 == None:
        if ensemble1 == None or ensemble2 == None:
            raise RuntimeError

        # Extract coordinates from ensembles
        coordinates_system1 = ensemble1.coordinates
        coordinates_system2 = ensemble2.coordinates

        # Average coordinates in the two systems
        x1 = numpy.average(coordinates_system1, axis=0).flatten()
        x2 = numpy.average(coordinates_system2, axis=0).flatten()

        # Covariance matrices in the two systems
        sigma1 = covariance_matrix(ensemble1,
                                   mass_weighted=mass_weighted,
                                   estimator=covariance_estimator)
        sigma2 = covariance_matrix(ensemble2,
                                   mass_weighted=mass_weighted,
                                   estimator=covariance_estimator)

    # Inverse covariance matrices
    sigma1_inv = numpy.linalg.pinv(sigma1)
    sigma2_inv = numpy.linalg.pinv(sigma2)

    # Difference between average vectors
    d_avg = x1 - x2

    # Sigma
    sigma = sigma1_inv + sigma2_inv

    # Distance measure
    trace = numpy.trace(numpy.dot(sigma1, sigma2_inv) +
                        numpy.dot(sigma2, sigma1_inv)
                        - 2 * numpy.identity(sigma1.shape[0]))

    d_hes = 0.25 * (numpy.dot(numpy.transpose(d_avg),
                              numpy.dot(sigma1_inv + sigma2_inv,
                                        d_avg)) + trace)
    return d_hes


def clustering_ensemble_similarity(cc, ens1, ens1_id, ens2, ens2_id):
    """Clustering ensemble similarity: calculate the probability densities from
     the clusters and calculate discrete Jensen-Shannon divergence.

    Parameters
    ----------

        cc : encore.ClustersCollection
            Collection from cluster calculated by a clustering algorithm
            (e.g. Affinity propagation)

        ens1 : encore.Ensemble
            First ensemble to be used in comparison

            ens2 : encore.Ensemble
                    Second ensemble to be used in comparison

        ens1_id : int
            First ensemble id as detailed in the ClustersCollection metadata

        ens2_id : int
            Second ensemble id as detailed in the ClustersCollection metadata

    Returns
    -------

        djs : float
            Jensen-Shannon divergence between the two ensembles, as calculated by
            the clustering ensemble similarity method
    """
    tmpA = numpy.array([numpy.where(c.metadata['ensemble'] == ens1_id)[
                            0].shape[0] / float(ens1.coordinates.shape[0]) for
                        c in cc])
    tmpB = numpy.array([numpy.where(c.metadata['ensemble'] == ens2_id)[
                            0].shape[0] / float(ens2.coordinates.shape[0]) for
                        c in cc])

    # Exclude clusters which have 0 elements in both ensembles    
    pA = tmpA[tmpA + tmpB > EPSILON]
    pB = tmpB[tmpA + tmpB > EPSILON]

    return discrete_jensen_shannon_divergence(pA, pB)


def cumulative_clustering_ensemble_similarity(cc, ens1, ens1_id, ens2, ens2_id,
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

        ens1 : encore.Ensemble
                First ensemble to be used in comparison

        ens2 : encore.Ensemble
                Second ensemble to be used in comparison

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

    ensA = [numpy.where(numpy.logical_and(c.metadata['ensemble'] <= ens1_id,
                                          c.metadata[
                                              'ensemble']) >= ens1_id_min)[
                0].shape[0] for c in cc]
    ensB = [numpy.where(numpy.logical_and(c.metadata['ensemble'] <= ens2_id,
                                          c.metadata[
                                              'ensemble']) >= ens2_id_min)[
                0].shape[0] for c in cc]
    sizeA = float(numpy.sum(ensA))
    sizeB = float(numpy.sum(ensB))
    # sizeA = float( numpy.sum(
    # [numpy.where( numpy.logical_and(c.metadata['ensemble'] <=
    # ens1_id, c.metadata['ensemble']) >= ens1_id_min)[0].shape[0] for c in cc])
    # sizeB = float(numpy.sum(
    # [numpy.where( numpy.logical_and(c.metadata['ensemble']
    # <= ens2_id, c.metadata['ensemble']) >= ens2_id_min)[0].shape[0] for c in cc])

    tmpA = numpy.array(ensA) / sizeA
    tmpB = numpy.array(ensB) / sizeB

    # Exclude clusters which have 0 elements in both ensembles
    pA = tmpA[tmpA + tmpB > EPSILON]
    pB = tmpB[tmpA + tmpB > EPSILON]

    return discrete_jensen_shannon_divergence(pA, pB)


def gen_kde_pdfs(embedded_space, ensemble_assignment, nensembles,
                 nsamples=None, **kwargs):
    """ 
    Generate Kernel Density Estimates (KDE) from embedded spaces and
    elaborate the coordinates for later use.

    Parameters
    ----------

        embedded_space : numpy.array
            Array containing the coordinates of the embedded space

        ensemble_assignment : numpy.array
            Array containing one int per ensemble conformation. These allow to
            distinguish, in the complete embedded space, which conformations belong
            to each ensemble. For instance if ensemble_assignment is [1,1,1,1,2,2],
            it means that the first four conformations belong to ensemble 1
            and the last two to ensemble 2

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
            List of numpy.array containing, each one, the elements of the embedded
             space belonging to a certain ensemble
    """
    kdes = []
    embedded_ensembles = []
    resamples = []

    for i in range(1, nensembles + 1):
        this_embedded = embedded_space.transpose()[
            numpy.where(ensemble_assignment == i)].transpose()
        embedded_ensembles.append(this_embedded)
        kdes.append(gaussian_kde(
            this_embedded))  # XXX support different bandwidth values

    # Set number of samples
    if not nsamples:
        nsamples = this_embedded.shape[1] * 10

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
    under the distribution P. We can, thus, just estimate the expectation values
    of the 	components to get an estimate of dKL.
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
            Samples drawn according do kde1. Will be used as samples to calculate
            the expected values according to 'P' as detailed before.

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
        ln_P1_exp_P1 = numpy.average(numpy.log(kde1.evaluate(resamples1)))
        ln_P2_exp_P2 = numpy.average(numpy.log(kde2.evaluate(resamples2)))
        ln_P1P2_exp_P1 = numpy.average(numpy.log(
            0.5 * (kde1.evaluate(resamples1) + kde2.evaluate(resamples1))))
        ln_P1P2_exp_P2 = numpy.average(numpy.log(
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
                array containing one int per ensemble conformation. These allow to
                distinguish, in the complete embedded space, which conformations
                belong to each ensemble. For instance if ensemble_assignment is
                [1,1,1,1,2,2], it means that the first four conformations belong
                to ensemble 1 and the last two to ensemble 2

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


    kdes = []
    embedded_ensembles = []
    resamples = []
    if not ens_id_max:
        ens_id_max = nensembles + 1
    for i in range(ens_id_min, ens_id_max + 1):
        this_embedded = embedded_space.transpose()[numpy.where(
            numpy.logical_and(ensemble_assignment >= ens_id_min,
                              ensemble_assignment <= i))].transpose()
        embedded_ensembles.append(this_embedded)
        kdes.append(
            gaussian_kde(this_embedded))  # XXX support different bandwidth values

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
            Basic filename for output. If None, no files will be written, and the
            matrix will be just printed on screen

        header : str
                Line to be written just before the matrix

        suffix : str
            String to be concatenated to basename, in order to get the final file
            name

        extension : str
            Extension for the output file

    """

    if base_fname != None:
        fname = base_fname + "-" + suffix + "." + extension
    else:
        fname = None
    matrix.square_print(header=header, fname=fname)


def write_output_line(value, fhandler=None, suffix="", label="win.", number=0,
                      rawline=None):
    """
    Write a line of data with a fixed format to standard output and optionally
    file. The line will be appended or written to a file object.
    The format is (in the Python str.format specification language):
     '{:s}{:d}\t{:.3f}', with the first element being the label, the second
     being
    a number that identifies the data point, and the third being the number
    itself. For instance:

    win.3	0.278

    Parameters
    ----------

        value : float
            Value to be printed.

        fhandler : file object
            File object in which the line will be written. if None, nothing will
            be written to file, and the value will be just printed on screen

        label : str
            Label to be written before the data

        number : int
            Number that identifies the data being written in this line.

        rawline : str
            If rawline is not None, write rawline to fhandler instead of the
            formatted number line. rawline can be any arbitrary string.
            """

    if fhandler == None:
        fh = Tee(sys.stdout)
    else:
        fh = Tee(sys.stdout, fhandler)

    if rawline != None:
        print >> fh, rawline
        return

    print >> fh, "{:s}{:d}\t{:.3f}".format(label, number, value)


def bootstrap_coordinates(coords, times):
    """
    Bootstrap conformations in a encore.Ensemble. This means drawing from the
    encore.Ensemble.coordinates numpy array with replacement "times" times
    and returning the outcome.

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
        this_coords = numpy.zeros(coords.shape)
        for c in range(this_coords.shape[0]):
            this_coords[c, :, :] = coords[numpy.random.randint(low=0, high=
            this_coords.shape[0]), :, :]
        out.append(this_coords)
    return out


def bootstrapped_matrix(matrix, ensemble_assignment):
    """
    Bootstrap an input square matrix. The resulting matrix will have the same
    shape as the original one, but the order of its elements will be drawn
    (with repetition). Separately bootstraps each ensemble.

    Parameters
    ----------

        matrix : encore.utils.TriangularMatrix
            similarity/dissimilarity matrix

    Returns
    -------

        this_m : encore.utils.TriangularMatrix
            bootstrapped similarity/dissimilarity matrix
    """
    ensemble_identifiers = numpy.unique(ensemble_assignment)
    this_m = TriangularMatrix(size=matrix.size)
    indexes = []
    for ens in ensemble_identifiers:
        old_indexes = numpy.where(ensemble_assignment == ens)[0]
        indexes.append(numpy.random.randint(low=numpy.min(old_indexes),
                                            high=numpy.max(old_indexes) + 1,
                                            size=old_indexes.shape[0]))

    indexes = numpy.hstack(indexes)
    for j in range(this_m.size):
        for k in range(j):
            this_m[j, k] = matrix[indexes[j], indexes[k]]

    logging.info("Matrix bootstrapped.")
    return this_m


def get_similarity_matrix(ensembles,
                          similarity_mode="minusrmsd",
                          load_matrix=None,
                          change_sign=False,
                          save_matrix=None,
                          superimpose=True,
                          superimposition_subset="name CA",
                          mass_weighted=True,
                          bootstrap_matrix=False,
                          bootstrapping_samples=100,
                          np=1):
    """
    Retrieves or calculates the similarity or conformational distance (RMSD) matrix.
    The similarity matrix is calculated between all the frames of all the 
    encore.Ensemble objects given as input. The order of the matrix elements depends on
    the order of the coordinates of the ensembles AND on the order of the 
    input ensembles themselves, therefore the ordering of the input list is significant.

    The similarity matrix can either be calculated from input Ensembles or
    loaded from an input numpy binary file. The signs of the elements of 
    the loaded matrix elements can be inverted using by the option `change_sign`.

    Please notice that the .npz file does not contain a bidimensional array,
    but a flattened representation that is meant to represent the elements of 
    an encore.utils.TriangularMatrix object.


    Parameters
    ----------
        ensembles : list
            List of ensembles

        similarity_mode : str, optional
            whether input matrix is smilarity matrix (minus RMSD) or
            a conformational distance matrix (RMSD). Accepted values 
            are "minusrmsd" and "rmsd".
        
        load_matrix : str, optional
            Load similarity/dissimilarity matrix from numpy binary file instead
            of calculating it (default is None). A filename is required.

        change_sign : bool, optional
            Change the sign of the elements of loaded matrix (default is False).
            Useful to switch between similarity/distance matrix.

        save_matrix : bool, optional
            Save calculated matrix as numpy binary file (default is None). A
            filename is required.

        superimpose : bool, optional
            Whether to superimpose structures before calculating distance
            (default is True).

        superimposition_subset : str, optional
            Group for superimposition using MDAnalysis selection syntax
            (default is Calpha atoms "name CA")

        mass_weighted : bool, optional
            calculate a mass-weighted RMSD (default is True). If set to False
            the superimposition will also not be mass-weighted.

        bootstrap_matrix : bool, optional
            Whether to bootstrap the similarity matrix (default is False).

        bootstrapping_samples : int, optional
            Number of times to bootstrap the similarity matrix (default is
            100).

        np : int, optional
            Maximum number of cores to be used (default is 1)

    Returns
    -------
        confdistmatrix : encore.utils.TriangularMatrix or list of encore.utils.TriangularMatrix
            Conformational distance or similarity matrix. If bootstrap_matrix
            is true, bootstrapping_samples matrixes are bootstrapped from the
            original one and they are returned as a list.
            
    


    """

    trajlist = []
    ensemble_assignment = []

    nensembles = len(ensembles)

    # Define ensemble assignments as required on the joined ensemble
    for i in range(1, nensembles + 1):
        ensemble_assignment += [i for j in ensembles[i - 1].coordinates]
    ensemble_assignment = numpy.array(ensemble_assignment)

    # Joined ensemble
    joined_ensemble = Ensemble(topology=ensembles[0].topology_filename,
                               trajectory=[ensembles[0].topology_filename],
                               atom_selection_string="all",
                               superimposition_selection_string=ensembles[
                                   0].superimposition_selection_string)

    # Joined ensemble coordinates as a concatenation of single ensembles
    # -  faster this way
    joined_ensemble.coordinates = numpy.concatenate(
        tuple([e.coordinates for e in ensembles]))
    joined_ensemble.superimposition_coordinates = numpy.concatenate(
        tuple([e.superimposition_coordinates for e in ensembles]))

    # Define metadata dictionary
    metadata = {'ensemble': ensemble_assignment}

    # Choose distance metric
    if similarity_mode == "minusrmsd":
        logging.info("    Similarity matrix: -RMSD matrix")
        matrix_builder = MinusRMSDMatrixGenerator()
    elif similarity_mode == "rmsd":
        logging.info("    Similarity matrix: RMSD matrix")
        matrix_builder = RMSDMatrixGenerator()
    else:
        logging.error(
            "Supported conformational distance measures are rmsd and minusrmsd")
        return None

    # Load the matrix if required
    if load:
        logging.info("        Loading similarity matrix from: %s" % load)
        confdistmatrix = TriangularMatrix(
            size=joined_ensemble.coordinates.shape[0], loadfile=load)
        logging.info("        Done!")
        for key in confdistmatrix.metadata.dtype.names:
            logging.info("        %s : %s" % (
                key, str(confdistmatrix.metadata[key][0])))

        # Change matrix sign if required. Useful to switch between
        # similarity/distance matrix.
        if change_sign:
            logging.info("        The matrix sign will be changed.")
            confdistmatrix.change_sign()

        # Check matrix size for consistency
        if not confdistmatrix.size == joined_ensemble.coordinates.shape[0]:
            logging.error(
                "ERROR: The size of the loaded matrix and of the ensemble"
                " do not match")
            return None

    # Calculate the matrix  
    else:
        logging.info(
            "        Perform pairwise alignment: %s" % str(superimpose))
        logging.info("        Mass-weighted alignment and RMSD: %s" % str(
            mass_weighted))
        if superimpose:
            logging.info(
                "        Atoms subset for alignment: %s" % superimposition_subset)
        logging.info("    Calculating similarity matrix . . .")

        # Use superimposition subset, if necessary. If the pairwise alignment is not required, it will not be performed anyway.
        if superimposition_subset:
            confdistmatrix = matrix_builder(
                joined_ensemble,
                pairwise_align=superimpose,
                align_subset_coordinates=
                                joined_ensemble.superimposition_coordinates,
                mass_weighted=mass_weighted,
                ncores=np)

        else:
            confdistmatrix = matrix_builder(joined_ensemble,
                                            pairwise_align=superimpose,
                                            mass_weighted=mass_weighted,
                                            ncores=np)

        logging.info("    Done!")

        if save:
            confdistmatrix.savez(save)

    if bootstrap_matrix:
        bs_args = [tuple([confdistmatrix, ensemble_assignment]) for i in
                   range(bootstrapping_samples)]

        pc = ParallelCalculation(np, bootstrapped_matrix, bs_args)

        pc_results = pc.run()

        bootstrap_matrices = zip(*pc_results)[1]

        return bootstrap_matrices

    return confdistmatrix


def prepare_ensembles_for_convergence_increasing_window(ensembles,
                                                        window_size):
    """
    Generate ensembles to be fed to ces_convergence or dres_convergence
    from a single ensemble. Basically, the different slices the algorithm
    needs are generated here.

    Parameters
    ----------

        ensembles : encore.Ensemble object
            Input ensemble

        window_size : int
            size of the window (in number of frames) to be used

    Returns
    -------

        tmp_ensembles : 
            the original ensemble is divided into ensembles, each being
            a window_size-long slice of the original ensemble. The last
            ensemble will be bigger if the length of the input ensemble
            is not exactly divisible by window_size.
    
    """

    ens_size = ensembles.coordinates.shape[0]

    rest_slices = ens_size / window_size
    residuals = ens_size % window_size
    slices_n = [0]

    tmp_ensembles = []

    for rs in range(rest_slices - 1):
        slices_n.append(slices_n[-1] + window_size)
    # if residuals != 0:
    #    slices_n.append(slices_n[-1] + residuals + window_size)
    # else:
    #    slices_n.append(slices_n[-1] + window_size)
    slices_n.append(slices_n[-1] + residuals + window_size)
    for s in range(len(slices_n) - 1):
        tmp_ensembles.append(Ensemble(
            topology=ensembles.topology_filename,
            trajectory=[ensembles.topology_filename],
            atom_selection_string=ensembles.atom_selection_string,
            superimposition_selection_string=ensembles.superimposition_selection_string))
        # print slices_n
        tmp_ensembles[-1].coordinates = ensembles.coordinates[
                                        slices_n[s]:slices_n[s + 1], :, :]

    return tmp_ensembles


def hes(ensembles,
        cov_estimator="shrinkage",
        mass_weighted=True,
        details=False,
        estimate_error=False,
        bootstrapping_runs=100,
        calc_diagonal=False):
    """

    Calculates the Harmonic Ensemble Similarity (HES) between ensembles using
    the symmetrized version of Kullback-Leibler divergence as described
    in [Lindorff-Larsen2009]_.


    Parameters
    ----------
        ensembles : list
            List of ensemble objects for similarity measurements.

        cov_estimator : str, optional
            Covariance matrix estimator method, either shrinkage, `shrinkage`,
            or Maximum Likelyhood, `ml`. Default is shrinkage.

        mass_weighted : bool, optional
            Whether to perform mass-weighted covariance matrix estimation
            (default is True).

        details : bool, optional 
            Save the mean and covariance matrix for each
            ensemble in a numpy array (default is False).

        estimate_error : bool, optional
            Whether to perform error estimation (default is False).

        bootstrapping_runs : int, optional
            Number of times the similarity matrix will be bootstrapped (default
            is 100).


    Returns
    -------
        hes : numpy.array (bidimensional)
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


    Example
    -------

    To calculate the Harmonic Ensemble similarity, two Ensemble objects are
    created from a topology file and two trajectories. The
    topology- and trajectory files used are obtained from the MDAnalysis
    test suite for two different simulations of the protein AdK. To run the
    examples see the module `Examples`_ for how to import the files: ::

    >>> ens1 = Ensemble(topology=PDB_small, trajectory=DCD)
    >>> ens2 = Ensemble(topology=PDB_small, trajectory=DCD2)
    >>> print hes([ens1, ens2])
    13946090.5764


    """

    logging.info("Chosen metric: Harmonic similarity")
    if cov_estimator == "shrinkage":
        covariance_estimator = EstimatorShrinkage()
        logging.info("    Covariance matrix estimator: Shrinkage")
    elif cov_estimator == "ml":
        covariance_estimator = EstimatorML()
        logging.info("    Covariance matrix estimator: Maximum Likelihood")
    else:
        logging.error(
            "Covariance estimator %s is not supported. "
            "Choose between 'shrinkage' and 'ml'." % cov_estimator)
        return None

    out_matrix_eln = len(ensembles)

    if calc_diagonal:
        pairs_indeces = list(trm_indeces_diag(out_matrix_eln))
    else:
        pairs_indeces = list(trm_indeces_nodiag(out_matrix_eln))
    xs = []
    sigmas = []

    if estimate_error:
        data = []
        for t in range(bootstrapping_runs):
            logging.info("The coordinates will be bootstrapped.")
            xs = []
            sigmas = []
            values = numpy.zeros((out_matrix_eln, out_matrix_eln))
            for e in ensembles:
                this_coords = bootstrap_coordinates(e.coordinates, 1)[0]
                xs.append(numpy.average(this_coords, axis=0).flatten())
                sigmas.append(covariance_matrix(e,
                                                mass_weighted=True,
                                                estimator=covariance_estimator))
            for i, j in pairs_indeces:
                value = harmonic_ensemble_similarity(x1=xs[i],
                                                     x2=xs[j],
                                                     sigma1=sigmas[i],
                                                     sigma2=sigmas[j])
                values[i, j] = value
                values[j, i] = value
            data.append(values)
        outs = numpy.array(data)
        avgs = numpy.average(data, axis=0)
        stds = numpy.std(data, axis=0)

        return (avgs, stds)

    # Calculate the parameters for the multivariate normal distribution
    # of each ensemble
    values = numpy.zeros((out_matrix_eln, out_matrix_eln))

    for e in ensembles:
        print e
        # Extract coordinates from each ensemble
        coordinates_system = e.coordinates

        # Average coordinates in each system
        xs.append(numpy.average(coordinates_system, axis=0).flatten())

        # Covariance matrices in each system
        sigmas.append(covariance_matrix(e,
                                        mass_weighted=mass_weighted,
                                        estimator=covariance_estimator))

    for i, j in pairs_indeces:
        value = harmonic_ensemble_similarity(x1=xs[i],
                                             x2=xs[j],
                                             sigma1=sigmas[i],
                                             sigma2=sigmas[j])
        values[i, j] = value
        values[j, i] = value

    if values.shape[0] == 2:
        values = values[0,1]
        
    # Save details as required
    if details:
        kwds = {}
        for i in range(out_matrix_eln):
            kwds['ensemble%d_mean' % (i + 1)] = xs[i]
            kwds['ensemble%d_covariance_matrix' % (i + 1)] = sigmas[i]
        details = numpy.array(kwds)

    else:
        details = None


        
    return values, details


def ces(ensembles,
        preference_values=-1.0,
        max_iterations=500,
        convergence=50,
        damping=0.9,
        noise=True,
        clustering_mode="ap",
        similarity_mode="minusrmsd",
        similarity_matrix=None,
        cluster_collections=None,
        estimate_error=False,
        bootstrapping_samples=100,
        details=False,
        np=1,
        calc_diagonal=False,
        **kwargs):
    """

    Calculates the Clustering Ensemble Similarity (CES) between ensembles
    using the Jensen-Shannon divergence as described in
    [Lindorff-Larsen2009]_.


    Parameters
    ----------

        ensembles : list
            List of ensemble objects for similarity measurements

        preference_values : float or iterable of floats, optional
            Preference parameter used in the Affinity Propagation algorithm for
            clustering  (default -1.0). A high preference value results in
            many clusters, a low preference will result in fewer numbers of
            clusters. Inputting a list of different preference values results
            in multiple calculations of the CES, one for each preference
            clustering.

        max_iterations : int, optional
            Parameter in the Affinity Propagation for
            clustering (default is 500).

        convergence : int, optional
            Minimum number of unchanging iterations to achieve convergence
            (default is 50). Parameter in the Affinity Propagation for
            clustering.

        damping : float, optional
            Damping factor (default is 0.9). Parameter in the Affinity
            Propagation for clustering.

        noise : bool, optional
            Apply noise to similarity matrix (default is True).

        clustering_mode : str, optional
            Choice of clustering algorithm. Only Affinity Propagation,`ap`,
            is implemented so far (default).

        similarity_matrix : By default, 

        estimate_error :  bool, optional
            Whether to perform error estimation (default is False).
            Only bootstrapping mode is supported so far.

        boostrapped_matrices : XXX

        details : XXX

        np : int, optional
            Maximum number of cores to be used (default is 1).

        **kwargs : XXX


    Returns
    -------
        ces : dict
            dictionary with the input preferences as keys and numpy.array of
            the clustering similarity as values.
            The clustering ensemble similarity are between each pair of
            ensembles measured by the Jensen-Shannon divergence.


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
    To calculate the Clustering Ensemble similarity, two Ensemble objects are
    created from a topology file and two trajectories. The
    topology- and trajectory files used are obtained from the MDAnalysis
    test suite for two different simulations of the protein AdK. To run the
    examples see the module `Examples`_ for how to import the files.
    Here the simplest case of just two :class:`Ensemble`s used for comparison
    are illustrated: ::

        >>> ens1 = Ensemble( topology = PDB_small, trajectory = DCD)
        >>> ens2 = Ensemble(topology = PDB_small, trajectory = DCD2)
        >>> CES = ces([ens1,ens2])
        >>> print CES
            (array([[[ 0.          0.55392484]
                     [ 0.55392484  0.        ]]])

    """


    if not hasattr(preference_values, '__iter__'):
        preference_values = [preference_values]
        full_output = False
    else:
        full_output = True
    try:
        preference_values = numpy.array(preference_values, dtype=numpy.float)
    except:
        raise TypeError("preferences expects a float or an iterable of numbers, such as a list of floats or a numpy.array")
        

    ensemble_assignment = []
    for i in range(1, len(ensembles) + 1):
        ensemble_assignment += [i for j in ensembles[i - 1].coordinates]
    ensemble_assignment = numpy.array(ensemble_assignment)

    metadata = {'ensemble': ensemble_assignment}

    out_matrix_eln = len(ensembles)

    if calc_diagonal:
        pairs_indeces = list(trm_indeces_diag(out_matrix_eln))
    else:
        pairs_indeces = list(trm_indeces_nodiag(out_matrix_eln))

    if similarity_matrix:
        confdistmatrix = similarity_matrix
    else:
        kwargs['similarity_mode'] = similarity_mode
        if not estimate_error:
            confdistmatrix = get_similarity_matrix(ensembles, **kwargs)
        else:
            confdistmatrix = get_similarity_matrix(
                ensembles,
                bootstrapping_samples=bootstrapping_samples,
                bootstrap_matrix=True,
                **kwargs)

    if clustering_mode == "ap":

        preferences = map(float, preference_values)

        logging.info("    Clustering algorithm: Affinity Propagation")
        logging.info("        Preference values: %s" % ", ".join(
            map(lambda x: "%3.2f" % x, preferences)))
        logging.info("        Maximum iterations: %d" % max_iterations)
        logging.info("        Convergence: %d" % convergence)
        logging.info("        Damping: %1.2f" % damping)
        logging.info("        Apply noise to similarity matrix: %s" % str(noise))

        # Choose clustering algorithm
        clustalgo = AffinityPropagation()

        # Prepare input for parallel calculation
        if estimate_error:
            bootstrap_matrices = confdistmatrix
            confdistmatrixs = []
            lams = []
            max_iterationss = []
            convergences = []
            noises = []
            real_prefs = []
            nmat = len(bootstrap_matrices)
            for p in preferences:
                confdistmatrixs.extend(bootstrap_matrices)
                lams.extend([damping] * nmat)
                max_iterationss.extend([max_iterations] * nmat)
                noises.extend([noise] * nmat)
                convergences.extend([convergence] * nmat)
                real_prefs.extend([p] * nmat)
            old_prefs = preferences
            preferences = real_prefs
        else:
            confdistmatrixs = [confdistmatrix for i in preferences]
            lams = [damping for i in preferences]
            max_iterationss = [max_iterations for i in preferences]
            convergences = [convergence for i in preferences]
            noises = [int(noise) for i in preferences]

        args = zip(confdistmatrixs, preferences, lams, max_iterationss,
                   convergences, noises)
        logging.info("    Starting affinity propagation runs . . .")

        # Do it
        pc = ParallelCalculation(np, clustalgo, args)

        results = pc.run()

        # Create clusters collections from clustering results, one for each cluster.
        #  None if clustering didn't work.
        ccs = [ClustersCollection(clusters[1], metadata=metadata) for clusters in
               results]

        if estimate_error:
            preferences = old_prefs
            k = 0
            values = {}
            avgs = {}
            stds = {}
            for i, p in enumerate(preferences):
                failed_runs = 0
                values[p] = []
                for j in range(len(bootstrap_matrices)):
                    if ccs[k].clusters == None:
                        failed_runs += 1
                        k += 1
                        continue
                    values[p].append(numpy.zeros((out_matrix_eln, out_matrix_eln)))

                    for pair in pairs_indeces:
                        # Calculate dJS
                        this_djs = clustering_ensemble_similarity(ccs[k],
                                                                  ensembles[
                                                                      pair[0]],
                                                                  pair[0] + 1,
                                                                  ensembles[
                                                                      pair[1]],
                                                                  pair[1] + 1)
                        values[p][-1][pair[0], pair[1]] = this_djs
                        values[p][-1][pair[1], pair[0]] = this_djs
                    k += 1
                outs = numpy.array(values[p])
                avgs[p] = numpy.average(outs, axis=0)
                stds[p] = numpy.std(outs, axis=0)

            return (avgs, stds)

        values = []
        kwds = {}
        for i, p in enumerate(preferences):
            if ccs[i].clusters == None:
                continue
            else:
                values.append(numpy.zeros((out_matrix_eln, out_matrix_eln)))

                for pair in pairs_indeces:
                    # Calculate dJS
                    this_val = clustering_ensemble_similarity(ccs[i],
                                                              ensembles[pair[0]],
                                                              pair[0] + 1,
                                                              ensembles[pair[1]],
                                                              pair[1] + 1)
                    values[-1][pair[0], pair[1]] = this_val
                    values[-1][pair[1], pair[0]] = this_val

            if details:
                kwds['centroids_pref%.3f' % p] = numpy.array(
                    [c.centroid for c in ccs[i]])
                kwds['ensemble_sizes'] = numpy.array(
                    [e.coordinates.shape[0] for e in ensembles])
                for cln, cluster in enumerate(ccs[i]):
                    kwds["cluster%d_pref%.3f" % (cln + 1, p)] = numpy.array(
                        cluster.elements)

    if full_output:
        values = numpy.array(values).swapaxes(0, 2)
    else:
        if len(ensembles) == 2:
            values = values[0][0, 1]
        else:
            values = values[0]
    
    if details:
        details = numpy.array(kwds)
    else:
        details = None

    return values, details


def dres(ensembles,
         conf_dist_mode="rmsd",
         conf_dist_matrix=None,
         mode='vanilla',
         dimensions=[3],
         maxlam=2.0,
         minlam=0.1,
         ncycle=100,
         nstep=10000,
         neighborhood_cutoff=1.5,
         kn=100,
         nsamples=1000,
         estimate_error=False,
         bootstrapping_samples=100,
         details=False,
         np=1,
         calc_diagonal = False,
         **kwargs):
    """

    Calculates the Dimensional Reduction Ensemble Similarity (DRES) between
    ensembles using the Jensen-Shannon divergence as described in
    [Lindorff-Larsen2009]_.


    Parameters
    ----------
        ensembles : list
            List of ensemble objects for similarity measurements

        conf_dist_matrix :

        mode : str, opt
            Which algorithm to use for dimensional reduction. Three options:
                - Stochastic Proximity Embedding (`vanilla`)  (default)
                - Random Neighborhood Stochastic Proximity Embedding (`rn`)
                - k-Nearest Neighbor Stochastic Proximity Embedding (`knn`)



        dimensions : int, optional
            Number of dimensions for reduction (default is 3)

        maxlam : float, optional
            Starting lambda learning rate parameter (default is 2.0). Parameter
            for Stochastic Proximity Embedding calculations.

        minlam : float, optional
            Final lambda learning rate (default is 0.1). Parameter
            for Stochastic Proximity Embedding calculations.

        ncycle : int, optional
            Number of cycles per run (default is 100). At the end of every
            cycle, lambda is changed.

        nstep : int, optional
            Number of steps per cycle (default is 10000)

        neighborhood_cutoff : float, optional
            Neighborhood cutoff (default is 1.5).

        kn : int, optional
            Number of neighbours to be considered (default is 100)

        estimate_error : bool, optional
            Whether to perform error estimation (default is False)

        boostrapped_matrices :
            XXX

        nsamples : int, optional
            Number of samples to be drawn from the ensembles (default is 1000).
            Parameter used in Kernel Density Estimates (KDE) from embedded
            spaces.

        details : bool, optional
            XXX

        np : int, optional
            Maximum number of cores to be used (default is 1).

        **kwargs :

    Returns
    -------

        dres : dict
            dictionary with the input dimensions as keys and numpy.array of
            the dres similarity as values.
            The similiarity is calculated betweem each pair of
            ensembles measured by the Jensen-Shannon divergence.

    Notes
    -----
    To calculate the similarity the method first projects the ensembles into lower
    dimensions by using the Stochastic Proximity Embedding algorithm. A
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
    To calculate the Dimensional Reduction Ensemble similarity, two Ensemble
    objects are created from a topology file and two trajectories. The
    topology- and trajectory files used are obtained from the MDAnalysis
    test suite for two different simulations of the protein AdK. To run the
    examples see the module `Examples`_ for how to import the files.
    Here the simplest case of comparing just two :class:`Ensemble`s are
    illustrated: ::


        >>> ens1 = Ensemble(topology=PDB_small,trajectory=DCD)
        >>> ens2 = Ensemble(topology=PDB_small,trajectory=DCD2)
        >>> DRES = dres([ens1,ens2])
        >>> print DRES
           (array( [[[ 0.          0.67383396]
                 [ 0.67383396  0.        ]]]



    """

    if not hasattr(dimensions, '__iter__'):
        dimensions = [dimensions]
        full_output = False
    else:
        full_output = True
    try:
        dimensions = numpy.array(dimensions, dtype=numpy.int)
    except:
        raise TypeError("preferences expects a float or an iterable of numbers, such as a list of floats or a numpy.array")

    stressfreq = -1

    out_matrix_eln = len(ensembles)

    if calc_diagonal:
        pairs_indeces = list(trm_indeces_diag(out_matrix_eln))
    else:
        pairs_indeces = list(trm_indeces_nodiag(out_matrix_eln))

    ensemble_assignment = []
    for i in range(1, len(ensembles) + 1):
        ensemble_assignment += [i for j in ensembles[i - 1].coordinates]
    ensemble_assignment = numpy.array(ensemble_assignment)

    metadata = {'ensemble': ensemble_assignment}

    if conf_dist_matrix:
        confdistmatrix = conf_dist_matrix
    else:
        kwargs['similarity_mode'] = conf_dist_mode
        if not estimate_error:
            confdistmatrix = get_similarity_matrix(ensembles, **kwargs)
        else:
            confdistmatrix = get_similarity_matrix(
                ensembles,
                bootstrapping_samples=bootstrapping_samples,
                bootstrap_matrix=True,
                **kwargs)

    dimensions = map(int, dimensions)

    # prepare runs. (e.g.: runs = [1,2,3,1,2,3,1,2,3, ...])
    if estimate_error:
        runs = []
        bootstrapped_matrices = confdistmatrix
        for d in dimensions:
            runs.extend([d] * len(bootstrapped_matrices))
        matrices = bootstrapped_matrices * len(bootstrapped_matrices)
    else:
        runs = dimensions
        matrices = [confdistmatrix for i in runs]

    # Choose algorithm and prepare options
    embedding_options = []
    if mode == 'vanilla':
        embedder = StochasticProximityEmbedding()
        for r in range(len(runs)):
            embedding_options += [(matrices[r],
                                   neighborhood_cutoff,
                                   runs[r],
                                   maxlam,
                                   minlam,
                                   ncycle,
                                   nstep,
                                   stressfreq)]

    if mode == 'rn':
        embedder = RandomNeighborhoodStochasticProximityEmbedding()
        for r in range(len(runs)):
            embedding_options += [(matrices[r],
                                   neighborhood_cutoff,
                                   kn,
                                   runs[r],
                                   maxlam,
                                   minlam,
                                   ncycle,
                                   stressfreq)]

    if mode == 'knn':
        embedder = kNNStochasticProximityEmbedding()
        for r in range(len(runs)):
            embedding_options += [(matrices[r],
                                   kn,
                                   runs[r],
                                   maxlam,
                                   minlam,
                                   ncycle,
                                   nstep,
                                   stressfreq)]

    pc = ParallelCalculation(np, embedder, embedding_options)

    # Run parallel calculation
    results = pc.run()
    sleep(1)

    embedded_spaces_perdim = {}
    stresses_perdim = {}

    # Sort out obtained spaces and their residual stress values

    if estimate_error:  # if bootstrap
        avgs = {}
        stds = {}
        values = {}
        k = 0
        for ndim in dimensions:
            values[ndim] = []
            for i in range(len(bootstrapped_matrices)):

                values[ndim].append(numpy.zeros((out_matrix_eln, out_matrix_eln)))

                embedded_stress = results[k][1][0]
                embedded_space = results[k][1][1]

                kdes, resamples, embedded_ensembles = gen_kde_pdfs(
                    embedded_space,
                    ensemble_assignment,
                    out_matrix_eln,
                    nsamples=nsamples)

                for pair in pairs_indeces:
                    this_value = dimred_ensemble_similarity(kdes[pair[0]],
                                                            resamples[pair[0]],
                                                            kdes[pair[1]],
                                                            resamples[pair[1]])
                    values[ndim][-1][pair[0], pair[1]] = this_value
                    values[ndim][-1][pair[1], pair[0]] = this_value

                k += 1
                outs = numpy.array(values[ndim])
                avgs[ndim] = numpy.average(outs, axis=0)
                stds[ndim] = numpy.std(outs, axis=0)

        return (avgs, stds)

    values = []

    for i in range(len(dimensions)):
        stresses_perdim[dimensions[i]] = []
        embedded_spaces_perdim[dimensions[i]] = []
        for j in range(1):
            stresses_perdim[dimensions[i]].append(
                results[j * len(dimensions) + i][1][0])
            embedded_spaces_perdim[dimensions[i]].append(
                results[j * len(dimensions) + i][1][1])

    kwds = {}

    for ndim in dimensions:

        values.append(numpy.zeros((len(ensembles), len(ensembles))))

        embedded_spaces = embedded_spaces_perdim[ndim]
        embedded_stresses = stresses_perdim[ndim]

        embedded_stress = embedded_stresses[numpy.argmin(embedded_stresses)]
        embedded_space = embedded_spaces[numpy.argmin(embedded_stresses)]

        kdes, resamples, embedded_ensembles = gen_kde_pdfs(embedded_space,
                                                           ensemble_assignment,
                                                           len(ensembles),
                                                           nsamples=nsamples)

        for pair in pairs_indeces:
            this_value = dimred_ensemble_similarity(kdes[pair[0]],
                                                    resamples[pair[0]],
                                                    kdes[pair[1]],
                                                    resamples[pair[1]])
            values[-1][pair[0], pair[1]] = this_value
            values[-1][pair[1], pair[0]] = this_value

        if details:
            kwds["stress_%ddims" % ndim] = numpy.array([embedded_stress])
            for en, e in enumerate(embedded_ensembles):
                kwds["ensemble%d_%ddims" % (en, ndim)] = e

    if full_output:
        values = numpy.array(values).swapaxes(0, 2)
    else:
        if len(ensembles) == 2:
            values = values[0][0, 1]
        else:
            values = values[0]

    if details:
        details = numpy.array(kwds)
    else:
        details = None

    return values, details


def ces_convergence(original_ensemble,
                    window_size,
                    similarity_mode="minusrmsd",
                    preference_values=[1.0],
                    max_iterations=500,
                    convergence=50,
                    damping=0.9,
                    noise=True,
                    save_matrix=None,
                    load_matrix=None,
                    np=1,
                    **kwargs):

    """ 
    Use the Clustering comparison measure to evaluate the convergence of the ensemble/trajectory


    Parameters
    ----------

        window_size : XXX
            Size of window  XXX

        preference_values : list , optional
            Preference parameter used in the Affinity Propagation algorithm for
            clustering  (default [-1.0]). A high preference value results in
            many clusters, a low preference will result in fewer numbers of
            clusters. Inputting a list of different preference values results
            in multiple calculations of the CES, one for each preference
            clustering.

        max_iterations : int, optional
            Parameter in the Affinity Propagation for
            clustering (default is 500).

        convergence : int, optional
            Minimum number of unchanging iterations to achieve convergence
            (default is 50). Parameter in the Affinity Propagation for
            clustering.

        damping : float, optional
            Damping factor (default is 0.9). Parameter in the Affinity
            Propagation for clustering.

        noise : bool, optional
            Apply noise to similarity matrix (default is True).

        save_matrix : bool, optional
            Save calculated matrix as numpy binary file (default None). A
            filename is required.

        load_matrix : str, optional
            Load similarity/dissimilarity matrix from numpy binary file instead
            of calculating it (default is None). A filename is required.

        np : int, optional
            Maximum number of cores to be used (default is 1).

        kwargs : XXX

    

    Returns
    -------

        out : XXX


   """

    if not hasattr(preference_values, '__iter__'):
        preferences = [preference_values]
        full_output = False
    else:
        full_output = True
    try:
        preferences = numpy.array(preference_values, dtype=numpy.float)
    except:
        raise TypeError("preferences expects a float or an iterable of numbers, such as a list of floats or a numpy.array")

    ensembles = prepare_ensembles_for_convergence_increasing_window(
        original_ensemble, window_size)

    kwargs['similarity_mode'] = similarity_mode
    confdistmatrix = get_similarity_matrix([original_ensemble], **kwargs)

    ensemble_assignment = []
    for i in range(1, len(ensembles) + 1):
        ensemble_assignment += [i for j in ensembles[i - 1].coordinates]
    ensemble_assignment = numpy.array(ensemble_assignment)

    metadata = {'ensemble': ensemble_assignment}

    preferences = preference_values

    logging.info("    Clustering algorithm: Affinity Propagation")
    logging.info("        Preference values: %s" % ", ".join(
        map(lambda x: "%3.2f" % x, preferences)))
    logging.info("        Maximum iterations: %d" % max_iterations)
    logging.info("        Convergence: %d" % convergence)
    logging.info("        Damping: %1.2f" % damping)
    logging.info("        Apply noise to similarity matrix: %s" % str(noise))

    confdistmatrixs = [confdistmatrix for i in preferences]
    lams = [damping for i in preferences]
    max_iterationss = [max_iterations for i in preferences]
    convergences = [convergence for i in preferences]
    noises = [int(noise) for i in preferences]

    clustalgo = AffinityPropagation()

    args = zip(confdistmatrixs, preferences, lams, max_iterationss,
               convergences, noises)

    logging.info("    Starting affinity propagation runs . . .")

    pc = ParallelCalculation(np, clustalgo, args)

    results = pc.run()

    logging.info("\n    Done!")
    ccs = [ClustersCollection(clusters[1], metadata=metadata) for clusters in
           results]

    
    out = []

    for i, p in enumerate(preferences):
        if ccs[i].clusters == None:
            continue
        out.append(numpy.zeros(len(ensembles)))
        for j in range(0, len(ensembles)):
            out[-1][j] = cumulative_clustering_ensemble_similarity(
                ccs[i],
                ensembles[ -1],
                len(ensembles) + 1,
                ensembles[j],
                j + 1)

    out = numpy.array(out).T
    return out



def dres_convergence(original_ensemble,
                     window_size,
                     conf_dist_mode='rmsd',
                     mode='vanilla',
                     dimensions=[3],
                     maxlam=2.0,
                     minlam=0.1,
                     ncycle=100,
                     nstep=10000,
                     neighborhood_cutoff=1.5,
                     kn=100,
                     nsamples=1000,
                     estimate_error=False,
                     bootstrapping_samples=100,
                     details=False,
                     np=1,
                     **kwargs):

    """
    Use the Dimensional Reduction comparison measure to evaluate the convergence of the ensemble/trajectory.

    Parameters
    ----------

        original_ensemble : XXX

        window_size : XXX

        mode : str, opt
            Which algorithm to use for dimensional reduction. Three options:
                - Stochastic Proximity Embedding (`vanilla`)  (default)
                - Random Neighborhood Stochastic Proximity Embedding (`rn`)
                - k-Nearest Neighbor Stochastic Proximity Embedding (`knn`)

        dimensions  : int, optional
            Number of dimensions for reduction (default is 3)

        maxlam : float, optional
            Starting lambda learning rate parameter (default is 2.0). Parameter
            for Stochastic Proximity Embedding calculations.

        minlam : float, optional
            Final lambda learning rate (default is 0.1). Parameter
            for Stochastic Proximity Embedding calculations.

        ncycle  :  int, optional
            Number of cycles per run (default is 100). At the end of every
            cycle, lambda is changed.

        nstep  : int, optional
            Number of steps per cycle (default is 10000)

        neighborhood_cutoff  : float, optional
            Neighborhood cutoff (default is 1.5).

        kn  : int, optional
            Number of neighbours to be considered (default is 100)

        nsamples : int, optional
            Number of samples to be drawn from the ensembles (default is 1000).
            Parameter used in Kernel Density Estimates (KDE) from embedded
            spaces.

        estimate_error  :  bool, optional
            Whether to perform error estimation (default is False)

        bootstrapping_samples  :  int, optional
            Number of bootstrapping runs XXX (default is 100).

        details  : bool, optional
            XXX  (default is False)

        np  : int, optional
            Maximum number of cores to be used (default is 1).

        kwargs : XXX


    Returns:
    --------

        out : XXX
            XXX




    """

    ensembles = prepare_ensembles_for_convergence_increasing_window(
        original_ensemble, window_size)

    kwargs['similarity_mode'] = conf_dist_mode
    confdistmatrix = get_similarity_matrix([original_ensemble], **kwargs)

    ensemble_assignment = []
    for i in range(1, len(ensembles) + 1):
        ensemble_assignment += [i for j in ensembles[i - 1].coordinates]
    ensemble_assignment = numpy.array(ensemble_assignment)

    out_matrix_eln = len(ensembles)

    runs = dimensions
    matrices = [confdistmatrix for i in runs]

    stressfreq = -1

    embedding_options = []
    if mode == 'vanilla':
        embedder = StochasticProximityEmbedding()
        for r in range(len(runs)):
            embedding_options += [(matrices[r],
                                   neighborhood_cutoff,
                                   runs[r],
                                   maxlam,
                                   minlam,
                                   ncycle,
                                   nstep,
                                   stressfreq)]

    if mode == 'rn':
        embedder = RandomNeighborhoodStochasticProximityEmbedding()
        for r in range(len(runs)):
            embedding_options += [(matrices[r],
                                   neighborhood_cutoff,
                                   kn,
                                   runs[r],
                                   maxlam,
                                   minlam,
                                   ncycle,
                                   stressfreq)]

    if mode == 'knn':
        embedder = kNNStochasticProximityEmbedding()
        for r in range(len(runs)):
            embedding_options += [(matrices[r],
                                   kn,
                                   runs[r],
                                   maxlam,
                                   minlam,
                                   ncycle,
                                   nstep,
                                   stressfreq)]

    pc = ParallelCalculation(np, embedder, embedding_options)

    results = pc.run()
    sleep(1)

    embedded_spaces_perdim = {}
    stresses_perdim = {}
    out = []

    for i in range(len(dimensions)):
        stresses_perdim[dimensions[i]] = []
        embedded_spaces_perdim[dimensions[i]] = []
        for j in range(1):
            stresses_perdim[dimensions[i]].append(
                results[j * len(dimensions) + i][1][0])
            embedded_spaces_perdim[dimensions[i]].append(
                results[j * len(dimensions) + i][1][1])

    # Run parallel calculation

    for ndim in dimensions:

        out.append(numpy.zeros(out_matrix_eln))

        embedded_spaces = embedded_spaces_perdim[ndim]
        embedded_stresses = stresses_perdim[ndim]

        embedded_stress = embedded_stresses[numpy.argmin(embedded_stresses)]
        embedded_space = embedded_spaces[numpy.argmin(embedded_stresses)]

        # For every chosen dimension value:

        kdes, resamples, embedded_ensembles = cumulative_gen_kde_pdfs(
            embedded_space, ensemble_assignment, out_matrix_eln - 1,
            nsamples=nsamples)

        for j in range(0, out_matrix_eln):
            out[-1][j] = dimred_ensemble_similarity(kdes[-1],
                                                      resamples[-1],
                                                      kdes[j],
                                                      resamples[j])

    out = numpy.array(out).T
    return out
