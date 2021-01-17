# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
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

from ...lib.util import iterable


class Clusters:
    """
    Class used to contain groups of clusters. Used in
    LeafletFinder algorithms. Can also be used in clustering.

    Parameters
    ----------
    predictor: function, object or class (optional)
        The associated predictor with the data. Must be provided
        if you want to use the ``run`` method


    Attributes
    ----------
    predictor: function, class or object
        Provided predictor
    cluster_indices: list of iterables of integers
        Indices of the provided data assigned to each cluster
    clusters_by_size: list of iterables of integers
        Indices of the provided data assigned to each cluster.
        The list is sorted by the size of the cluster (largest first)
    outlier_indices: iterable of indices
        Indices of the provided data considered to be outliers
    data_labels: iterable
        List of cluster labels in the order of the provided data
    """
    def __init__(self, predictor=None, **kwargs):
        if isinstance(predictor, type):
            predictor = predictor(**kwargs)
        self.predictor = predictor
        self.cluster_indices = []
        self.clusters_by_size = []
        self.outlier_indices = []
        self.data_labels = []
        self._data = None

    def __len__(self):
        return len(self.cluster_indices)

    def run(self, data, **kwargs):
        """
        Run the given predictor on the given data,
        and assign clusters.
        """
        try:
            self.data_labels = self.predictor.fit_predict(data)
        except AttributeError:
            try:
                self.data_labels = self.predictor(data, **kwargs)
            except TypeError:
                raise ValueError("predictor must have a `fit_predict` or be "
                                "a method that returns an iterable of "
                                f"labels. Given {self.predictor}")
        if not iterable(self.data_labels):
            raise ValueError("predictor must return an iterable "
                             "of labels")
        self._data = data
        ix = np.argsort(self.data_labels)
        indices = np.arange(len(data))
        splix = np.where(np.ediff1d(self.data_labels[ix]))[0] + 1
        cluster_indices = np.split(indices[ix], splix)
        self.cluster_indices = [np.sort(x) for x in cluster_indices]
        if self.data_labels[ix[0]] == -1:
            self.outlier_indices = self.cluster_indices.pop(0)
        self.clusters_by_size = sorted(self.cluster_indices,
                                       key=lambda x: len(x), reverse=True)

    def set_clusters(self, cluster_indices):
        """
        Set cluster data from values. No predictor is necessary here.
        """
        self.cluster_indices = cluster_indices
        self.clusters_by_size = sorted(self.cluster_indices,
                                       key=lambda x: len(x), reverse=True)
        labels = np.zeros(sum(map(len, cluster_indices)))
        for i, cl in enumerate(cluster_indices):
            labels[cl] = i
        self.data_labels = labels
