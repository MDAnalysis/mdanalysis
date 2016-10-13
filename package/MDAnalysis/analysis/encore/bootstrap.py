# bootstrap.py --- Bootstrap analysis for ensembles and distance matrices
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
bootstrap procedures --- :mod:`MDAnalysis.analysis.ensemble.bootstrap`
=====================================================================


The module contains functions for bootstrapping either ensembles (Universe
objects) or distance matrices, by resampling with replacement.

:Author: Matteo Tiberti, Wouter Boomsma, Tone Bengtsen
:Year: 2015--2016
:Copyright: GNU Public License v3
:Mantainer: Matteo Tiberti <matteo.tiberti@gmail.com>, mtiberti on github

.. versionadded:: 0.16.0

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
        indexes.append(np.random.randint(low=np.min(old_indexes),
                                         high=np.max(old_indexes) + 1,
                                         size=old_indexes.shape[0]))

    indexes = np.hstack(indexes)
    for j in range(this_m.size):
        for k in range(j):
            this_m[j, k] = matrix[indexes[j], indexes[k]]

    logging.info("Matrix bootstrapped.")
    return this_m


def get_distance_matrix_bootstrap_samples(distance_matrix,
                                          ensemble_assignment,
                                          samples=100,
                                          ncores=1):
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

    bs_args = \
            [([distance_matrix, ensemble_assignment]) for i in range(samples)]

    pc = ParallelCalculation(ncores, bootstrapped_matrix, bs_args)

    pc_results = pc.run()

    bootstrap_matrices = zip(*pc_results)[1]

    return bootstrap_matrices


def get_ensemble_bootstrap_samples(ensemble,
                                   samples=100):
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
            size=ensemble.trajectory.timeseries().shape[1])
        ensembles.append(
            mda.Universe(ensemble.filename,
                        ensemble.trajectory.timeseries(format='afc')[:,indices,:],
                         format=mda.coordinates.memory.MemoryReader))
    return ensembles