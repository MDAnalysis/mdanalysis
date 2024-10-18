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
clustering frontend --- :mod:`MDAnalysis.analysis.encore.clustering.ClusteringMethod`
=====================================================================================

The module defines classes for interfacing to various clustering algorithms.
One has been implemented natively, and will always be available, while
others are available only if scikit-learn is installed

:Author: Matteo Tiberti, Wouter Boomsma, Tone Bengtsen

.. versionadded:: 0.16.0

.. deprecated:: 2.8.0
   This module is deprecated in favour of the 
   MDAKit `mdaencore <https://mdanalysis.org/mdaencore/>`_ and will be removed
   in MDAnalysis 3.0.0.

"""
import numpy as np
import warnings
import logging

# Import native affinity propagation implementation
from . import affinityprop

# Attempt to import scikit-learn clustering algorithms
try:
    import sklearn.cluster
except ImportError:
    sklearn = None
    msg = "sklearn.cluster could not be imported: some functionality will " \
          "not be available in encore.fit_clusters()"
    warnings.warn(msg, category=ImportWarning)
    logging.warning(msg)
    del msg


def encode_centroid_info(clusters, cluster_centers_indices):
    """
    Adjust cluster indices to include centroid information
    as described in documentation for ClusterCollection
    """
    values, indices = np.unique(clusters, return_inverse=True)
    for c_center in cluster_centers_indices:
        if clusters[c_center] != c_center:
            values[indices[c_center]] = c_center
    return values[indices]


class ClusteringMethod (object):
    """
    Base class for any Clustering Method
    """

    # Whether the method accepts a distance matrix
    accepts_distance_matrix=True

    def __call__(self, x):
        """
        Parameters
        ----------

        x
            either trajectory coordinate data (np.array) or an
            encore.utils.TriangularMatrix, encoding the conformational
            distance matrix

        Raises
        ------
        NotImplementedError
           Method or behavior needs to be defined by a subclass    
        
        """
        raise NotImplementedError("Class {0} doesn't implement __call__()"
                                  .format(self.__class__.__name__))


class AffinityPropagationNative(ClusteringMethod):
    """
    Interface to the natively implemented Affinity propagation procedure.
    """
    def __init__(self,
                 damping=0.9, preference=-1.0,
                 max_iter=500, convergence_iter=50,
                 add_noise=True):
        """
        Parameters
        ----------

        damping : float, optional
            Damping factor (default is 0.9). Parameter for the Affinity
            Propagation for clustering.

        preference : float, optional
            Preference parameter used in the Affinity Propagation algorithm for
            clustering  (default -1.0). A high preference value results in
            many clusters, a low preference will result in fewer numbers of
            clusters.

        max_iter : int, optional
            Maximum number of iterations for affinity propagation (default is
            500).

        convergence_iter : int, optional
            Minimum number of unchanging iterations to achieve convergence
            (default is 50). Parameter in the Affinity Propagation for
            clustering.

        add_noise : bool, optional
            Apply noise to similarity matrix before running clustering
            (default is True)

        """
        self.damping = damping
        self.preference = preference
        self.max_iter = max_iter
        self.convergence_iter = convergence_iter
        self.add_noise = add_noise

    def __call__(self, distance_matrix):
        """
        Parameters
        ----------

        distance_matrix : encore.utils.TriangularMatrix
            conformational distance matrix


        Returns
        -------
        numpy.array : array, shape(n_elements) 
            centroid frames of the clusters for all of the elements

        .. versionchanged:: 1.0.0
           This method no longer returns ``details``
        """
        clusters = affinityprop.AffinityPropagation(
            s=distance_matrix * -1.,   # invert sign
            preference=self.preference,
            lam=self.damping,
            max_iterations = self.max_iter,
            convergence = self.convergence_iter,
            noise=int(self.add_noise))
        
        return clusters
if sklearn:

    class AffinityPropagation(ClusteringMethod):
        """
        Interface to the Affinity propagation clustering procedure implemented
        in sklearn.
        """

        def __init__(self,
                     damping=0.9, preference=-1.0,
                     max_iter=500, convergence_iter=50,
                     **kwargs):
            """
            Parameters
            ----------

            damping : float, optional
                Damping factor (default is 0.9). Parameter for the Affinity
                Propagation for clustering.

            preference : float, optional
                Preference parameter used in the Affinity Propagation algorithm
                for clustering  (default -1.0). A high preference value results
                in many clusters, a low preference will result in fewer numbers
                of clusters.

            max_iter : int, optional
                Maximum number of iterations for affinity propagation (default
                is 500).

            convergence_iter : int, optional
                Minimum number of unchanging iterations to achieve convergence
                (default is 50). Parameter in the Affinity Propagation for
                clustering.

            **kwargs : optional
                Other keyword arguments are passed to :class:`sklearn.cluster.AffinityPropagation`.

            """
            self.ap = \
                sklearn.cluster.AffinityPropagation(
                    damping=damping,
                    preference=preference,
                    max_iter=max_iter,
                    convergence_iter=convergence_iter,
                    affinity="precomputed",
                    **kwargs)

        def __call__(self, distance_matrix):
            """
            Parameters
            ----------

            distance_matrix : encore.utils.TriangularMatrix
                conformational distance matrix

            Returns
            -------
            numpy.array : array, shape(n_elements) 
                centroid frames of the clusters for all of the elements

            .. versionchanged:: 1.0.0
               This method no longer returns ``details``
            """
            logging.info("Starting Affinity Propagation: {0}".format
                         (self.ap.get_params()))

            # Convert from distance matrix to similarity matrix
            similarity_matrix = distance_matrix.as_array() * -1
            clusters = self.ap.fit_predict(similarity_matrix)
            clusters = encode_centroid_info(clusters,
                                            self.ap.cluster_centers_indices_)
            
            return clusters



    class DBSCAN(ClusteringMethod):
        """
        Interface to the DBSCAN clustering procedure implemented in sklearn.
        """
        def __init__(self,
                     eps=0.5,
                     min_samples=5,
                     algorithm="auto",
                     leaf_size=30,
                     **kwargs):
            """
            Parameters
            ----------

            eps : float, optional (default = 0.5)
                The maximum distance between two samples for them to be
                considered as in the same neighborhood.

            min_samples : int, optional (default = 5)
                The number of samples (or total weight) in a neighborhood for
                a point to be considered as a core point. This includes the
                point itself.

            algorithm : {'auto', 'ball_tree', 'kd_tree', 'brute'}, optional
                The algorithm to be used by the NearestNeighbors module
                to compute pointwise distances and find nearest neighbors.
                See NearestNeighbors module documentation for details.

            leaf_size : int, optional (default = 30)
                Leaf size passed to BallTree or cKDTree. This can affect the
                speed of the construction and query, as well as the memory
                required to store the tree. The optimal value depends
                on the nature of the problem.

            sample_weight : array, shape (n_samples,), optional
                Weight of each sample, such that a sample with a weight of at
                least ``min_samples`` is by itself a core sample; a sample with
                negative weight may inhibit its eps-neighbor from being core.
                Note that weights are absolute, and default to 1.

            """

            self.dbscan = sklearn.cluster.DBSCAN(eps=eps,
                                                 min_samples = min_samples,
                                                 algorithm=algorithm,
                                                 leaf_size = leaf_size,
                                                 metric="precomputed",
                                                 **kwargs)

        def __call__(self, distance_matrix):
            """
            Parameters
            ----------

            distance_matrix : encore.utils.TriangularMatrix
                conformational distance matrix


            Returns
            -------
            numpy.array : array, shape(n_elements) 
                centroid frames of the clusters for all of the elements

            .. versionchanged:: 1.0.0
               This method no longer returns ``details``
            """
            logging.info("Starting DBSCAN: {0}".format(
                self.dbscan.get_params()))
            clusters = self.dbscan.fit_predict(distance_matrix.as_array())
            if np.min(clusters == -1):
                clusters += 1
            # No centroid information is provided by DBSCAN, so we just
            # pick random members
            cluster_representatives = np.unique(clusters, return_index=True)[1]
            clusters = encode_centroid_info(clusters,
                                            cluster_representatives)
          
            return clusters

    class KMeans(ClusteringMethod):

        # Whether the method accepts a distance matrix
        accepts_distance_matrix = False

        """
        Interface to the KMeans clustering procedure implemented in sklearn.
        """
        def __init__(self,
                     n_clusters,
                     max_iter=300,
                     n_init=10,
                     init='k-means++',
                     algorithm="auto",
                     tol=1e-4,
                     verbose=False,
                     random_state=None,
                     copy_x=True,
                     **kwargs):
            """
            Parameters
            ----------
            n_clusters : int
                The number of clusters to form as well as the number of
                centroids to generate.

            max_iter : int, optional (default 300)
                Maximum number of iterations of the k-means algorithm to run.

            n_init : int, optional (default 10)
                Number of time the k-means algorithm will be run with different
                centroid seeds. The final results will be the best output of
                n_init consecutive runs in terms of inertia.

            init : {'k-means++', 'random', or ndarray, or a callable}, optional
                Method for initialization, default to 'k-means++':
                'k-means++' : selects initial cluster centers for k-mean
                clustering in a smart way to speed up convergence. See section
                Notes in k_init for more details.
                'random': generate k centroids from a Gaussian with mean and
                variance estimated from the data.
                If an ndarray is passed, it should be of shape
                (n_clusters, n_features) and gives the initial centers.
                If a callable is passed, it should take arguments X, k and
                and a random state and return an initialization.

            tol : float, optional (default 1e-4)
                The relative increment in the results before declaring
                convergence.

            verbose : boolean, optional (default False)
                Verbosity mode.

            random_state : integer or numpy.RandomState, optional
                The generator used to initialize the centers. If an integer is
                given, it fixes the seed. Defaults to the global numpy random
                number generator.

            copy_x : boolean, optional
                When pre-computing distances it is more numerically accurate to
                center the data first.  If copy_x is True, then the original
                data is not modified.  If False, the original data is modified,
                and put back before the function returns, but small numerical
                differences may be introduced by subtracting and then adding
                the data mean.

            """
            self.kmeans = sklearn.cluster.KMeans(n_clusters=n_clusters,
                                                 max_iter=max_iter,
                                                 n_init=n_init,
                                                 init=init,
                                                 tol=tol,
                                                 verbose=verbose,
                                                 random_state=random_state,
                                                 copy_x=copy_x,
                                                 **kwargs)

        def __call__(self, coordinates):
            """
            Parameters
            ----------

            coordinates : np.array
                trajectory atom coordinates


            Returns
            -------
            numpy.array : array, shape(n_elements) 
                centroid frames of the clusters for all of the elements

            .. versionchanged:: 1.0.0
               This method no longer returns ``details``
            """
            logging.info("Starting Kmeans: {0}".format(
                         (self.kmeans.get_params())))
            clusters = self.kmeans.fit_predict(coordinates)
            distances = self.kmeans.transform(coordinates)
            cluster_center_indices = np.argmin(distances, axis=0)
            clusters = encode_centroid_info(clusters,
                                             cluster_center_indices)
            
            return clusters
