# DimensionalityReductionMethod.py --- Interface classes to various
# dimensionality reduction algorithms
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
dimensionality reduction frontend --- :mod:`MDAnalysis.analysis.encore.clustering.DimensionalityReductionMethod`
=====================================================================

The module defines classes for interfacing to various dimensionality reduction
algorithms. One has been implemented natively, and will always be available,
while others are available only if scikit-learn is installed

:Author: Matteo Tiberti, Wouter Boomsma, Tone Bengtsen
:Year: 2015--2016
:Copyright: GNU Public License v3
:Mantainer: Matteo Tiberti <matteo.tiberti@gmail.com>, mtiberti on github

.. versionadded:: 0.16.0

"""

import numpy as np
import logging

# Import native affinity propagation implementation
from . import stochasticproxembed

# Attempt to import scikit-learn clustering algorithms
try:
    import sklearn.decomposition
except ImportError:
   sklearn = None
   msg = "sklearn.decomposition could not be imported: some functionality will " \
         "not be available in encore.dimensionality_reduction()"
   warnings.warn(msg, category=ImportWarning)
   logger.warn(msg)
   del msg


class DimensionalityReductionMethod (object):
    """
    Base class for any Dimensionality Reduction Method
    """

    # Whether the method accepts a distance matrix
    accepts_distance_matrix=True

    def __call__(self, x):
        """
        Parameters
        ----------

        x
            either trajectory coordinate data (np.array) or an
            encore.utils.TriangularMatrix, encoding the conformational
            distance matrix


        Returns
        -------
        numpy.array
            coordinates in reduced space

        """
        raise NotImplementedError("Class {0} doesn't implement __call__()"
                                  .format(self.__class__.__name__))


class StochasticProximityEmbeddingNative(DimensionalityReductionMethod):
    """
    Interface to the natively implemented Affinity propagation procedure.
    """
    def __init__(self,
                 dimension = 2,
                 distance_cutoff = 1.5,
                 min_lam = 0.1,
                 max_lam = 2.0,
                 ncycle = 100,
                 nstep = 10000,
                 stressfreq = -1):
        """
        Parameters
        ----------

        dimension : int
            Number of dimensions to which the conformational space will be reduced
            to (default is 3).

        min_lam : float, optional
            Final lambda learning rate (default is 0.1). Parameter
            for Stochastic Proximity Embedding calculations.

        max_lam : float, optional
            Starting lambda learning rate parameter (default is 2.0). Parameter
            for Stochastic Proximity Embedding calculations.

        ncycle : int, optional
            Number of cycles per run (default is 100). At the end of every
            cycle, lambda is changed.

        nstep : int, optional
            Number of steps per cycle (default is 10000)

        `stressfreq` : int
            calculate and report stress value every stressfreq cycle

        """
        self.dimension = dimension
        self.distance_cutoff = distance_cutoff
        self.min_lam = min_lam
        self.max_lam = max_lam
        self.ncycle = ncycle
        self.nstep = nstep
        self.stressfreq = stressfreq

    def __call__(self, distance_matrix):
        """
        Parameters
        ----------

        distance_matrix : encore.utils.TriangularMatrix
            conformational distance matrix


        Returns
        -------
        numpy.array
            coordinates in reduced space

        """
        final_stress, coordinates = \
            stochasticproxembed.StochasticProximityEmbedding().run(
            s=distance_matrix,
            rco=self.distance_cutoff,
            dim=self.dimension,
            minlam = self.min_lam,
            maxlam = self.max_lam,
            ncycle = self.ncycle,
            nstep = self.nstep,
            stressfreq=-1
        )
        return coordinates, {"final_stress": final_stress}



if sklearn:

    class PrincipleComponentAnalysis(DimensionalityReductionMethod):
        """
        Interface to the PCA dimensionality reduction method implemented in
        sklearn.
        """

        # Whether the method accepts a distance matrix
        accepts_distance_matrix = False

        def __init__(self,
                     dimension = 2,
                     **kwargs):
            """
            Parameters
            ----------

            dimension : int
                Number of dimensions to which the conformational space will be reduced
                to (default is 3).
            """
            self.pca = sklearn.decomposition.PCA(n_components=dimension,
                                                 **kwargs)

        def __call__(self, coordinates):
            """
            Parameters
            ----------

            coordinates : np.array
                trajectory atom coordinates


            Returns
            -------
            numpy.array
                coordinates in reduced space
            """
            coordinates = self.pca.fit_transform(coordinates)
            return coordinates.T, {}
