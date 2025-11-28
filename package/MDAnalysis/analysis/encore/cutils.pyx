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
# cython embedsignature=True

"""
Mixed Cython utils for ENCORE

:Author: Matteo Tiberti, Wouter Boomsma, Tone Bengtsen

.. versionadded:: 0.16.0
"""


import numpy as np
cimport numpy as cnp
import cython
from libc.math cimport sqrt

cnp.import_array()


@cython.boundscheck(False)
@cython.wraparound(False)
def PureRMSD(cnp.ndarray[cnp.float64_t, ndim=2] coordsi,
             cnp.ndarray[cnp.float64_t, ndim=2] coordsj,
             int atomsn,
             cnp.ndarray[cnp.float64_t, ndim=1] masses,
             double summasses):

    cdef  int k
    cdef double normsum

    normsum = 0.0

    for k in xrange(atomsn):
        normsum += masses[k] * (
            (coordsi[k, 0] - coordsj[k, 0])**2 +
            (coordsi[k, 1] - coordsj[k, 1])**2 +
            (coordsi[k, 2] - coordsj[k, 2])**2
        )

    return sqrt(normsum/summasses)
