# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""
Diffusion map --- :mod:`MDAnalysis.analysis.pca`
=====================================================================

:Authors: John Detlefs
:Year: 2016
:Copyright: GNU Public License v3

TODO, write documentation

"""
from six.moves import range
import logging

import numpy as np
from MDAnalysis import Universe

from .base import AnalysisBase

logger = logging.getLogger("MDAnalysis.analysis.pca")

class PCA(object):
    """Principal component analysis on an MD trajectory

    Attributes
    ----------
    components: array, (n_components, n_atoms)
        The principal components of the feature space,
        representing the directions of maximum variance in the data.

    explained_variance_ratio : array, (n_components, )
        Percentage of variance explained by each of the selected components.
        If a subset of components is not chosen then all components are stored
        and the sum of explained variances is equal to 1.0




    Methods
    -------
    fit(traj)

    transform(traj, n_components)

    inverse_tranform(pc_space)
    """

    def __init__(self, u, n_components=None
                 **kwargs):
        """
        Parameters
        ----------
        u : MDAnalysis Universe
            The universe containing the trajectory frames for Principal
            Component Analysis.
        """
        self.u = u
        self._calculated = False

    def fit(self, traj=None, n_components = None, start=None, step=None,
            stop=None):
        """ Use a subset of the trajectory to generate principal components """
        if traj is None:
            traj = self.u.trajectory
        elif isinstance(traj, u.trajectory):
            traj = traj
        else:
            raise AttributeError("Trajectory can not be found.")

        self._setup_frames(traj, start, stop, step)



    def transform(self, traj=None, n_components=None):
        """Apply the dimensionality reduction on a trajectory

        Parameters
        ----------
        traj : MDAnalysis Trajectory
            Trajectory for PCA transformation

        Returns
        -------
        pca_space : array, shape (number of atoms, number of components)
        """
        # apply Transform
        #
        pass

    def inverse_tranform(self, pca_space):
        """ Transform PCA-transformed data back to original configuration space.

        Parameters
        ----------
        pca_space : array
            The space corresponding to the trajectory coordinates after the PCA
            transformation.

        Returns
        -------
        original_traj : array, shape (number of frames, number of atoms)
        """
        pass
