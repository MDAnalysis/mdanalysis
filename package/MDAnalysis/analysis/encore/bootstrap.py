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
bootstrap procedures --- :mod:`MDAnalysis.analysis.ensemble.bootstrap`
======================================================================


The module contains functions for bootstrapping either ensembles (Universe
objects) or distance matrices, by resampling with replacement.

:Author: Matteo Tiberti, Wouter Boomsma, Tone Bengtsen

.. versionadded:: 0.16.0

.. deprecated:: 2.8.0
   This module is deprecated in favour of the 
   MDAKit `mdaencore <https://mdanalysis.org/mdaencore/>`_ and will be removed
   in MDAnalysis 3.0.0.

"""
import numpy as np
import logging
import MDAnalysis as mda
from .utils import TriangularMatrix, ParallelCalculation


def bootstrapped_matrix(matrix, ensemble_assignment):
    """
    Bootstrap an input square matrix. The resulting matrix will have the same
    shape as the original one, but the order of its elements will be drawn
    (with repetition). Separately bootstraps each ensemble.

    Parameters
    ----------

    matrix : encore.utils.TriangularMatrix
        similarity/dissimilarity matrix

    ensemble_assignment: numpy.array
        array of ensemble assignments. This array must be matrix.size long.

    Returns
    -------

    this_m : encore.utils.TriangularMatrix
        bootstrapped similarity/dissimilarity matrix
    """
    ensemble_identifiers = np.unique(ensemble_assignment)
    this_m = TriangularMatrix(size=matrix.size)
    indexes = []
    for ens in ensemble_identifiers:
        old_indexes = np.where(ensemble_assignment == ens)[0]
        indexes.append(
            np.random.randint(
                low=np.min(old_indexes),
                high=np.max(old_indexes) + 1,
                size=old_indexes.shape[0],
            )
        )

    indexes = np.hstack(indexes)
    for j in range(this_m.size):
        for k in range(j):
            this_m[j, k] = matrix[indexes[j], indexes[k]]

    logging.info("Matrix bootstrapped.")
    return this_m


def get_distance_matrix_bootstrap_samples(
    distance_matrix, ensemble_assignment, samples=100, ncores=1
):
    """
    Calculates distance matrices corresponding to bootstrapped ensembles, by
    resampling with replacement.

    Parameters
    ----------

    distance_matrix : encore.utils.TriangularMatrix
        Conformational distance matrix

    ensemble_assignment : str
        Mapping from frames to which ensemble they are from (necessary because
        ensembles are bootstrapped independently)

    samples : int, optional
        How many bootstrap samples to create.

    ncores : int, optional
        Maximum number of cores to be used (default is 1)

    Returns
    -------

    confdistmatrix : list of encore.utils.TriangularMatrix
    """

    bs_args = [
        ([distance_matrix, ensemble_assignment]) for i in range(samples)
    ]

    pc = ParallelCalculation(ncores, bootstrapped_matrix, bs_args)

    pc_results = pc.run()

    bootstrap_matrices = list(zip(*pc_results))[1]

    return bootstrap_matrices


def get_ensemble_bootstrap_samples(ensemble, samples=100):
    """
    Generates a bootstrapped ensemble by resampling with replacement.

    Parameters
    ----------

    ensemble : MDAnalysis.Universe
        Conformational distance matrix

    samples : int, optional
        How many bootstrap samples to create.

    Returns
    -------

    list of MDAnalysis.Universe objects
    """

    ensemble.transfer_to_memory()

    ensembles = []
    for i in range(samples):
        indices = np.random.randint(
            low=0,
            high=ensemble.trajectory.timeseries().shape[1],
            size=ensemble.trajectory.timeseries().shape[1],
        )
        ensembles.append(
            mda.Universe(
                ensemble.filename,
                ensemble.trajectory.timeseries(order="fac")[indices, :, :],
                format=mda.coordinates.memory.MemoryReader,
            )
        )
    return ensembles
