# cluster.py --- Common function for calling clustering algorithms
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
clustering frontend --- :mod:`MDAnalysis.analysis.encore.clustering.cluster`
=====================================================================

The module defines a function serving as front-end for various clustering
algorithms, wrapping them to allow them to be used interchangably.

:Author: Matteo Tiberti, Wouter Boomsma, Tone Bengtsen
:Year: 2015--2016
:Copyright: GNU Public License v3
:Mantainer: Matteo Tiberti <matteo.tiberti@gmail.com>, mtiberti on github

.. versionadded:: 0.16.0

"""

import numpy as np
from ..utils import ParallelCalculation, merge_universes
from .ClusterCollection import ClusterCollection
from ..confdistmatrix import get_distance_matrix
from . import ClusteringMethod


def cluster(ensembles,
            method = ClusteringMethod.AffinityPropagationNative(),
            selection="name CA",
            distance_matrix=None,
            allow_collapsed_result=True,
            ncores=1,
            **kwargs):
    """
    Cluster frames from one or more ensembles, using one or more
    clustering methods. The function optionally takes pre-calculated distances
    matrices as an argument. Note that not all clustering procedure can work
    directly on distance matrices, so the distance matrices might be ignored
    for particular choices of method.


    Parameters
    ----------

    ensembles : MDAnalysis.Universe, or list or list of list thereof
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

    selection : str, optional
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
    AdK. To run the examples see the module `Examples`_ for how to import the
    files.
    Here, we reduce cluster two ensembles ::
        >>> ens1 = Universe(PSF, DCD)
        >>> ens2 = Universe(PSF, DCD2)
        >>> cluster_collection = encore.cluster([ens1,ens2])
        >>> print cluster_collection

    You can change the parameters of the clustering method by explicitly
    specifying the method ::

        >>> cluster_collection = \
                encore.cluster( \
                     [ens1,ens2], \
                     method=encore.AffinityPropagationNative(preference=-2.))

    Here is an illustration using DBSCAN algorithm, instead
    of the default clustering method ::

        >>> cluster_collection = \
                encore.cluster( \
                     [ens1,ens2], \
                     method=encore.DBSCAN())

    You can also combine multiple methods in one call ::

        >>> cluster_collection = \
                encore.cluster( \
                     [ens1,ens2], \
                     method=[encore.AffinityPropagationNative(preference=-1.), \
                             encore.AffinityPropagationNative(preference=-2.)])

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
        np.any([method.accepts_distance_matrix for method in methods])

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
                                                           selection=selection,
                                                           **kwargs))

    args = []
    for method in methods:
        if method.accepts_distance_matrix:
            args += [(d,) for d in distance_matrix]
        else:
            for merged_ensemble in merged_ensembles:
                coordinates = merged_ensemble.trajectory.timeseries(format="fac")

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
        for i in range(len(ensembles)):
            ensemble_assignment += [i+1]*len(ensembles[i].trajectory)
        ensemble_assignment = np.array(ensemble_assignment)
        metadata = {'ensemble_membership': ensemble_assignment}

    # Create clusters collections from clustering results,
    # one for each cluster. None if clustering didn't work.
    ccs = [ClusterCollection(clusters[1][0],
                             metadata=metadata) for clusters in results]

    if allow_collapsed_result and len(ccs) == 1:
        ccs = ccs[0]

    return ccs