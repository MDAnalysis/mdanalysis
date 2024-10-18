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
clustering frontend --- :mod:`MDAnalysis.analysis.encore.clustering.cluster`
============================================================================

The module defines a function serving as front-end for various clustering
algorithms, wrapping them to allow them to be used interchangably.

:Author: Matteo Tiberti, Wouter Boomsma, Tone Bengtsen

.. versionadded:: 0.16.0

.. deprecated:: 2.8.0
   This module is deprecated in favour of the 
   MDAKit `mdaencore <https://mdanalysis.org/mdaencore/>`_ and will be removed
   in MDAnalysis 3.0.0.

"""
import numpy as np
from ..utils import ParallelCalculation, merge_universes
from .ClusterCollection import ClusterCollection
from ..confdistmatrix import get_distance_matrix
from . import ClusteringMethod


def cluster(ensembles,
            method = ClusteringMethod.AffinityPropagationNative(),
            select="name CA",
            distance_matrix=None,
            allow_collapsed_result=True,
            ncores=1,
            **kwargs):
    """Cluster frames from one or more ensembles, using one or more
    clustering methods. The function optionally takes pre-calculated distances
    matrices as an argument. Note that not all clustering procedure can work
    directly on distance matrices, so the distance matrices might be ignored
    for particular choices of method.


    Parameters
    ----------

    ensembles : MDAnalysis.Universe, or list, or list of list thereof
        The function takes either a single Universe object, a list of Universe
        objects or a list of lists of Universe objects. If given a single
        universe, it simply clusters the conformations in the trajectory. If
        given a list of ensembles, it will merge them and cluster them together,
        keeping track of the ensemble to which each of the conformations belong.
        Finally, if passed a list of list of ensembles, the function will just
        repeat the functionality just described - merging ensembles for each
        ensemble in the outer loop.

    method: encore.ClusteringMethod or list thereof, optional
        A single or a list of instances of the Clustering classes from
        the clustering module. A separate analysis will be run for each
        method. Note that different parameters for the same clustering method
        can be explored by adding different instances of the same clustering
        class.

    select : str, optional
        Atom selection string in the MDAnalysis format. Default is "name CA"

    distance_matrix : encore.utils.TriangularMatrix or list thereof, optional
        Distance matrix used for clustering. If this parameter
        is not supplied the matrix will be calculated on the fly.
        If several distance matrices are supplied, an analysis will be done
        for each of them. The number of provided distance matrices should
        match the number of provided ensembles.

    allow_collapsed_result: bool, optional
        Whether a return value of a list of one value should be collapsed
        into just the value.

    ncores : int, optional
        Maximum number of cores to be used (default is 1).


    Returns
    -------

    list of ClustersCollection objects (or potentially a single
    ClusteringCollection object if allow_collapsed_result is set to True)


    Example
    -------
    Two ensembles are created as Universe object using a topology file and
    two trajectories. The topology- and trajectory files used are obtained
    from the MDAnalysis test suite for two different simulations of the protein
    AdK.
    Here, we cluster two ensembles ::

        >>> from MDAnalysis import Universe
        >>> import MDAnalysis.analysis.encore as encore
        >>> from MDAnalysis.tests.datafiles import PSF, DCD, DCD2
        >>> ens1 = Universe(PSF, DCD)
        >>> ens2 = Universe(PSF, DCD2)
        >>> cluster_collection = encore.cluster([ens1,ens2])
        >>> print cluster_collection

    You can change the parameters of the clustering method by explicitly
    specifying the method ::

        >>> cluster_collection =
                encore.cluster(
                     [ens1,ens2],
                     method=encore.AffinityPropagationNative(preference=-2.))

    Here is an illustration using DBSCAN algorithm, instead
    of the default clustering method ::

        >>> cluster_collection =
                encore.cluster(
                     [ens1,ens2],
                     method=encore.DBSCAN())

    You can also combine multiple methods in one call ::

        >>> cluster_collection =
                encore.cluster(
                     [ens1,ens2],
                     method=[encore.AffinityPropagationNative(preference=-1.),
                             encore.AffinityPropagationNative(preference=-2.)])

    In addition to standard cluster membership information, the
    `cluster_collection` output keep track of the origin of each
    conformation, so you check how the different trajectories are
    represented in each cluster. Here, for brevity, we print just the
    members of the two first clusters ::

        >>> print [cluster.metadata["ensemble_membership"]
                     for cluster in cluster_collection][:2]
        [array([1, 1, 1, 1, 2]), array([1, 1, 1, 1, 1])]

    """

    # Internally, ensembles are always transformed to a list of lists
    if ensembles is not None:
        if not hasattr(ensembles, '__iter__'):
            ensembles = [ensembles]

        ensembles_list = ensembles
        if not hasattr(ensembles[0], '__iter__'):
            ensembles_list = [ensembles]

        # Calculate merged ensembles and transfer to memory
        merged_ensembles = []
        for ensembles in ensembles_list:
            # Transfer ensembles to memory
            for ensemble in ensembles:
                ensemble.transfer_to_memory()
            merged_ensembles.append(merge_universes(ensembles))

    methods = method
    if not hasattr(method, '__iter__'):
        methods = [method]

    # Check whether any of the clustering methods can make use of a distance
    # matrix
    any_method_accept_distance_matrix = \
        np.any([_method.accepts_distance_matrix for _method in methods])

    # If distance matrices are provided, check that it matches the number
    # of ensembles
    if distance_matrix:
        if not hasattr(distance_matrix, '__iter__'):
            distance_matrix = [distance_matrix]
        if ensembles is not None and \
                        len(distance_matrix) != len(merged_ensembles):
            raise ValueError("Dimensions of provided list of distance matrices "
                             "does not match that of provided list of "
                             "ensembles: {0} vs {1}"
                             .format(len(distance_matrix),
                                     len(merged_ensembles)))

    else:
        # Calculate distance matrices for all merged ensembles - if not provided
        if any_method_accept_distance_matrix:
            distance_matrix = []
            for merged_ensemble in merged_ensembles:
                distance_matrix.append(get_distance_matrix(merged_ensemble,
                                                           select=select,
                                                           **kwargs))

    args = []
    for method in methods:
        if method.accepts_distance_matrix:
            args += [(d,) for d in distance_matrix]
        else:
            for merged_ensemble in merged_ensembles:
                coordinates = merged_ensemble.trajectory.timeseries(order="fac")

                # Flatten coordinate matrix into n_frame x n_coordinates
                coordinates = np.reshape(coordinates,
                                         (coordinates.shape[0], -1))

                args.append((coordinates,))

    # Execute clustering procedure
    pc = ParallelCalculation(ncores, methods, args)

    # Run parallel calculation
    results = pc.run()

    # Keep track of which sample belongs to which ensembles
    metadata = None
    if ensembles is not None:
        ensemble_assignment = []
        for i, ensemble in enumerate(ensembles):
            ensemble_assignment += [i+1]*len(ensemble.trajectory)
        ensemble_assignment = np.array(ensemble_assignment)
        metadata = {'ensemble_membership': ensemble_assignment}

    # Create clusters collections from clustering results,
    # one for each cluster. None if clustering didn't work.
    ccs = [ClusterCollection(clusters[1],
                             metadata=metadata) for clusters in results]

    if allow_collapsed_result and len(ccs) == 1:
        ccs = ccs[0]

    return ccs
