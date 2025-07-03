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
dimensionality reduction frontend --- :mod:`MDAnalysis.analysis.encore.clustering.DimensionalityReductionMethod`
================================================================================================================

The module defines classes for interfacing to various dimensionality reduction
algorithms. One has been implemented natively, and will always be available,
while others are available only if scikit-learn is installed

:Author: Matteo Tiberti, Wouter Boomsma, Tone Bengtsen

.. versionadded:: 0.16.0

.. deprecated:: 2.8.0
   This module is deprecated in favour of the 
   MDAKit `mdaencore <https://mdanalysis.org/mdaencore/>`_ and will be removed
   in MDAnalysis 3.0.0.

"""
import logging
import warnings

# Import native affinity propagation implementation
from . import stochasticproxembed

# Attempt to import scikit-learn clustering algorithms
try:
    import sklearn.decomposition
except ImportError:
    sklearn = None
    import warnings

    warnings.warn(
        "sklearn.decomposition could not be imported: some "
        "functionality will not be available in "
        "encore.dimensionality_reduction()",
        category=ImportWarning,
    )


class DimensionalityReductionMethod(object):
    """
    Base class for any Dimensionality Reduction Method
    """

    # Whether the method accepts a distance matrix
    accepts_distance_matrix = True

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
        raise NotImplementedError(
            "Class {0} doesn't implement __call__()".format(
                self.__class__.__name__
            )
        )


class StochasticProximityEmbeddingNative(DimensionalityReductionMethod):
    """
    Interface to the natively implemented Affinity propagation procedure.
    """

    def __init__(
        self,
        dimension=2,
        distance_cutoff=1.5,
        min_lam=0.1,
        max_lam=2.0,
        ncycle=100,
        nstep=10000,
    ):
        """
        Parameters
        ----------

        dimension : int
            Number of dimensions to which the conformational space will be
            reduced to (default is 3).

        min_lam : float, optional
            Final lambda learning rate (default is 0.1).

        max_lam : float, optional
            Starting lambda learning rate parameter (default is 2.0).

        ncycle : int, optional
            Number of cycles per run (default is 100). At the end of every
            cycle, lambda is updated.

        nstep : int, optional
            Number of steps per cycle (default is 10000)

        """
        self.dimension = dimension
        self.distance_cutoff = distance_cutoff
        self.min_lam = min_lam
        self.max_lam = max_lam
        self.ncycle = ncycle
        self.nstep = nstep
        self.stressfreq = -1

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
        final_stress, coordinates = (
            stochasticproxembed.StochasticProximityEmbedding(
                s=distance_matrix,
                rco=self.distance_cutoff,
                dim=self.dimension,
                minlam=self.min_lam,
                maxlam=self.max_lam,
                ncycle=self.ncycle,
                nstep=self.nstep,
                stressfreq=self.stressfreq,
            )
        )
        return coordinates, {"final_stress": final_stress}


if sklearn:

    class PrincipalComponentAnalysis(DimensionalityReductionMethod):
        """
        Interface to the PCA dimensionality reduction method implemented in
        sklearn.
        """

        # Whether the method accepts a distance matrix
        accepts_distance_matrix = False

        def __init__(self, dimension=2, **kwargs):
            """
            Parameters
            ----------

            dimension : int
                Number of dimensions to which the conformational space will be
                reduced to (default is 3).
            """
            self.pca = sklearn.decomposition.PCA(
                n_components=dimension, **kwargs
            )

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
