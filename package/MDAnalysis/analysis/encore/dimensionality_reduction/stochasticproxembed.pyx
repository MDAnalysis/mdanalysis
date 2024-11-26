# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
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
Cython wrapper for the C implementation of the Stochastic Proximity Embedding
dimensionality reduction algorithm.

:Author: Matteo Tiberti, Wouter Boomsma
"""


import logging
import numpy
cimport numpy

numpy.import_array()
cimport cython

cdef extern from "spe.h":
    double CStochasticProximityEmbedding(double*, double*, double, int, int, double, double, int, int, int)


@cython.embedsignature(True)
def StochasticProximityEmbedding(s, double rco, int dim, double maxlam, double minlam, int ncycle, int nstep, int stressfreq):
    """
    Stochastic proximity embedding dimensionality reduction algorithm. The
    algorithm implemented here is described in this paper:

    Dmitrii N. Rassokhin, Dimitris K. Agrafiotis
    A modified update rule for stochastic proximity embedding
    Journal of Molecular Graphics and Modelling 22 (2003) 133–140

    This class is a Cython wrapper for a C implementation (see spe.c)

    Parameters
    ----------
    s : encore.utils.TriangularMatrix object
        Triangular matrix containing the distance values for each pair of
        elements in the original space.
    rco : float
        neighborhood distance cut-off
    dim : int
        number of dimensions for the embedded space
    minlam  : float
        final learning parameter
    maxlam  : float
        starting learning parameter
    ncycle : int
        number of cycles. Each cycle is composed of nstep steps. At the end
        of each cycle, the lerning parameter lambda is updated.
    nstep : int
        number of coordinate update steps for each cycle


    Returns
    -------
    space : (float, numpy.array)
        float is the final stress obtained; the array are the coordinates of
        the elements in the embedded space
    stressfreq : int
        calculate and report stress value every stressfreq cycle


    """

    cdef int nelem = s.size
    cdef double finalstress = 0.0

    logging.info("Starting Stochastic Proximity Embedding")

    cdef numpy.ndarray[numpy.float64_t,  ndim=1] matndarray = numpy.ascontiguousarray(s._elements, dtype=numpy.float64)
    cdef numpy.ndarray[numpy.float64_t,  ndim=1] d_coords   = numpy.zeros((nelem*dim),dtype=numpy.float64)

    finalstress = CStochasticProximityEmbedding(<double*>matndarray.data, <double*>d_coords.data, rco, nelem, dim, maxlam, minlam, ncycle, nstep, stressfreq)

    logging.info("Stochastic Proximity Embedding finished. Residual stress: %.3f" % finalstress)

    return (finalstress, d_coords.reshape((-1,dim)).T)
