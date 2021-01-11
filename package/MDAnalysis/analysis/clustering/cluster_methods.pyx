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
import warnings

import numpy as np
cimport numpy as np
cimport cython

cdef extern from "ap.h":
    int CAffinityPropagation(float*, int, float, int, int, int, long*)

@cython.boundscheck(False)
@cython.wraparound(False)
def affinity_propagation(similarity, preference, float lam, int max_iter, int conv_threshold, bint noise=True):
    """Affinity propagation clustering algorithm. This class is a Cython wrapper around the Affinity propagation algorithm, which is implement as a C library (see ap.c). The implemented algorithm is described in the paper:

    Clustering by Passing Messages Between Data Points.
    Brendan J. Frey and Delbert Dueck, University of Toronto
    Science 315, 972–976, February 2007

    Parameters
    ----------
    similarity : ndarray
        Square array
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
    
    cdef int cn = len(similarity)
    cdef int n_iter

    similarity[np.diag_indices(cn)] = preference
    indices = np.tril_indices(cn)

    # print(similarity)
    # print(np.tril(similarity))
    
    cdef np.ndarray[np.float32_t, ndim=1] sim = np.ravel(similarity[indices]).astype(np.float32)
    # sim = np.ascontiguousarray(sim)

    print(similarity[-2:])

    cdef np.ndarray[long, ndim=1] clusters = np.zeros((cn), dtype=long)

    # run C module Affinity Propagation
    n_iter = CAffinityPropagation(<float*>sim.data, cn, lam, max_iter, conv_threshold,
                                  noise, <long*>clusters.data)

    # Provide warning in case of lack of convergence
    if n_iter == 0:
        msg = ("Clustering with preference {:3.2f} "
               "did not fully converge in {:d} iterations")
        warnings.warn(msg.format(preference, -n_iter))

    return clusters
