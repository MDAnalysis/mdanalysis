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
Cluster representation --- :mod:`MDAnalysis.analysis.encore.clustering.ClusterCollection`
=========================================================================================

The module contains the Cluster and ClusterCollection classes which are
designed to store results from clustering algorithms.

:Author: Matteo Tiberti, Wouter Boomsma, Tone Bengtsen

.. versionadded:: 0.16.0

"""
import numpy as np


class Cluster(object):
    """
    Generic Cluster class for clusters with centroids.

    Attributes
    ----------

    id : int
        Cluster ID number. Useful for the ClustersCollection class

    metadata : iterable
        dict of lists or numpy.array, containing metadata for the cluster
        elements. The iterable must return the same number of elements as
        those that belong to the cluster.

    size : int
        number of elements.

    centroid : element object
        cluster centroid.

    elements : numpy.array
        array containing the cluster elements.
   """

    def __init__(self, elem_list=None, centroid=None, idn=None, metadata=None):
        """Class constructor. If elem_list is None, an empty cluster is created
            and the remaining arguments ignored.

        Parameters
        ----------

        elem_list : numpy.array or None
            numpy array of cluster elements

        centroid : None or element object
            centroid

        idn : int
            cluster ID

        metadata : iterable
            metadata, one value for each cluster element. The iterable
            must have the same length as the elements array.

    """

        self.id = idn

        if elem_list is None:
            self.size = 0
            self.elements = np.array([])
            self.centroid = None
            self.metadata = {}
            return

        self.metadata = {}
        self.elements = elem_list
        if centroid not in self.elements:
            raise LookupError("Centroid of cluster not found in the element list")

        self.centroid = centroid
        self.size = self.elements.shape[0]
        if metadata:
            for name, data in metadata.items():
                if len(data) != self.size:
                    raise TypeError('Size of metadata having label "{0}" '
                                    'is not equal to the number of cluster '
                                    'elements'.format(name))
            self.add_metadata(name, data)

    def __iter__(self):
        """
        Iterate over elements in cluster
        """
        return iter(self.elements)

    def __len__(self):
        """
        Size of cluster
        """
        return len(self.elements)

    def add_metadata(self, name, data):
        if len(data) != self.size:
            raise TypeError("Size of metadata is not equal to the number of "
                            "cluster elements")
        self.metadata[name] = np.array(data)

    def __repr__(self):
        """
        Textual representation
        """
        if self.size == 0:
            return "<Cluster with no elements>"
        else:
            return "<Cluster with {0} elements, centroid={1}, id={2}>".format(
                                                               self.size,
                                                               self.centroid,
                                                               self.id)

class ClusterCollection(object):
    """Clusters collection class; this class represents the results of a full
    clustering run. It stores a group of clusters defined as
    encore.clustering.Cluster objects.

    Attributes
    ----------

    clusters : list
        list of of Cluster objects which are part of the Cluster collection

"""

    def __init__(self, elements=None, metadata=None):
        """Class constructor. If elements is None, an empty cluster collection
            will be created. Otherwise, the constructor takes as input an
            iterable of ints, for instance:

            [ a, a, a, a, b, b, b, c, c, ... , z, z ]

            the variables a,b,c,...,z are cluster centroids, here as cluster
            element numbers (i.e. 3 means the 4th element of the ordered input
            for clustering). The array maps a correspondence between
            cluster elements (which are implicitly associated with the
            position in the array) with centroids, i. e. defines clusters.
            For instance:

            [ 1, 1, 1, 4, 4, 5 ]

            means that elements 0, 1, 2 form a cluster which has 1 as centroid,
            elements 3 and 4 form a cluster which has 4 as centroid, and
            element 5 has its own cluster.


            Parameters
            ----------

            elements : iterable of ints or None
                clustering results. See the previous description for details

            metadata : {str:list, str:list,...} or None
                metadata for the data elements. The list must be of the same
                size as the elements array, with one value per element.

        """
        idn = 0
        if elements is None:
            self.clusters = None
            return

        if not len(set((type(el) for el in elements))) == 1:
            raise TypeError("all the elements must have the same type")
        self.clusters = []
        elements_array = np.array(elements)
        centroids = np.unique(elements_array)
        for i in centroids:
            if elements[i] != i:
                raise ValueError("element {0}, which is a centroid, doesn't "
                                 "belong to its own cluster".format(
                                     elements[i]))
        for c in centroids:
            this_metadata = {}
            this_array = np.where(elements_array == c)
            if metadata:
                for k, v in metadata.items():
                    this_metadata[k] = np.asarray(v)[this_array]
            self.clusters.append(
                Cluster(elem_list=this_array[0], idn=idn, centroid=c,
                        metadata=this_metadata))

            idn += 1

    def get_ids(self):
        """
        Get the ID numbers of the clusters

        Returns
        -------

        ids : list of int
        list of cluster ids
        """
        return [v.id for v in self.clusters]

    def get_centroids(self):
        """
        Get the centroids of the clusters

        Returns
        -------

        centroids : list of cluster element objects
        list of cluster centroids
        """

        return [v.centroid for v in self.clusters]

    def __iter__(self):
        """
        Iterate over clusters

        """
        return iter(self.clusters)

    def __len__(self):
        """
        Length of clustering collection
        """
        return len(self.clusters)

    def __repr__(self):
        """
        Textual representation
        """
        if self.clusters is None:
            return "<ClusterCollection with no clusters>"
        else:
            return "<ClusterCollection with {0} clusters>".format(
                                             len(self.clusters))
