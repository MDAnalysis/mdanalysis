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
Cython wrapper for the C implementation of the Affinity Perturbation clustering algorithm.

:Author: Matteo Tiberti, Wouter Boomsma, Tone Bengtsen

"""
from ..utils import TriangularMatrix
import logging
import numpy
cimport numpy
cimport cython

numpy.import_array()

cdef extern from "ap.h":
    int CAffinityPropagation(float*, int, float, int, int, bint, long*)

@cython.boundscheck(False)
@cython.wraparound(False)
def AffinityPropagation(s, preference, float lam, int max_iterations, int convergence, int noise=1):
    """Affinity propagation clustering algorithm. This class is a Cython wrapper around the Affinity propagation algorithm, which is implement as a C library (see ap.c). The implemented algorithm is described in the paper:

    Clustering by Passing Messages Between Data Points.
    Brendan J. Frey and Delbert Dueck, University of Toronto
    Science 315, 972–976, February 2007

    Parameters
    ----------
    s : encore.utils.TriangularMatrix object
        Triangular matrix containing the similarity values for each pair of
        clustering elements. Notice that the current implementation does not
        allow for asymmetric values (i.e. similarity(a,b) is assumed to be
        equal to similarity(b,a))
    preference : numpy.array of floats or float
        Preference values, which the determine the number of clusters. If a
        single value is given, all the preference values are set to that.
        Otherwise, the list is used to set the preference values (one value per
        element, so the list must be of the same size as the number of
        elements)
    lam : float
        Floating point value that defines how much damping is applied to the
        solution at each iteration. Must be ]0,1]
    max_iterations : int
        Maximum number of iterations
    convergence : int
        Number of iterations in which the cluster centers must remain the same
        in order to reach convergence
    noise : int
        Whether to apply noise to the input s matrix, such there are no equal
        values. 1 is for yes, 0 is for no.

    Returns
    -------
    elements : list of int or None
        List of cluster-assigned elements, which can be used by
        encore.utils.ClustersCollection to generate Cluster objects. See these
        classes for more details.

    """
    cdef int cn = s.size
    cdef float cpreference = preference

    # Assign preference values to diagonal
    try:
        for i in xrange(s.size):
            s[i,i] = <float>preference[i]
    except:
        pass

    if type(preference) == float:
        for i in xrange(s.size):
            s[i,i] = <float>preference
    else:
        raise TypeError ("Preference should be of type float")

    logging.info("Preference %3.2f: starting Affinity Propagation" % (preference))

    # Prepare input and ouput arrays
    cdef numpy.ndarray[numpy.float32_t,  ndim=1] matndarray = numpy.ascontiguousarray(s._elements, dtype=numpy.float32)
    cdef numpy.ndarray[long,   ndim=1] clusters   = numpy.zeros((s.size),dtype=numpy.dtype("long"))

    # run C module Affinity Propagation
    iterations = CAffinityPropagation(<float*>matndarray.data, cn, lam, max_iterations, convergence, noise, <long*>clusters.data)

    # Provide warning in case of lack of convergence
    if iterations == 0:
        logging.info("Preference %3.2f: could not converge in %d iterations" % (preference, -iterations))
        import warnings
        warnings.warn("Clustering with preference {0:3.2f} did not fully converge in {1:d} iterations".format(preference, -iterations))

    # Find centroids
    centroids = numpy.unique(clusters)
    for k in numpy.arange(centroids.shape[0]):
        ii = numpy.where(clusters == centroids[k])[0]
        small_mat = numpy.zeros((ii.shape[0], ii.shape[0]))
        for ii1 in numpy.arange(ii.shape[0]):
            for ii2 in numpy.arange(ii.shape[0]):
                small_mat[ii1,ii2] = s[ ii[ii1], ii[ii2] ]
        j = numpy.argmax(numpy.sum(small_mat, axis=0))

        centroids[k] = ii[j]

    # Similarity to centroids
    S_centroids = numpy.zeros((s.size, centroids.shape[0]))
    for line in numpy.arange(s.size):
        for c in numpy.arange(centroids.shape[0]):
            S_centroids[line,c] = s[line, centroids[c]]

    # Center values for each observation
    c = numpy.argmax(S_centroids, axis=1)

    # Centroids should point to themselves
    c[centroids] = numpy.arange(centroids.shape[0])

    # Assign centroid indices to all observables
    clusters = centroids[c]

    logging.info("Preference %3.2f: converged in %d iterations" % (preference, iterations))

    return clusters
