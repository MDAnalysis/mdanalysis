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

"""\
Set box dimensions --- :mod: `MDAnalysis.transformations.boxdimensions`
=======================================================================

Set dimensions of the simulation box to a constant vector across all timesteps.


..autocalss:: setdimensions
"""
import numpy as np


class SetDimensions:
    """
    Set simulation box dimensions.

    Timestep dimensions are modified in place.

    Example
    -------

    e.g. set simulation box dimensions to a vector containing unit cell
    dimensions [*a*, *b*, *c*, *alpha*, *beta*, *gamma*], lengths *a*,
    *b*, *c* are in the MDAnalysis length unit (Å), and angles are in degrees.

    dim = [2, 2, 2, 90, 90, 90]
    transform = mda.transformations.boxdimensions.setdimensions(dim)
    u.trajectory.add_transformations(transform)

    Parameters
    ----------
    dimensions: array-like
        vector that contains unit cell lengths and angles.
        Expected shapes are (6, 0) or (1, 6)

    Returns
    -------
    MDAnalysis.coordinates.base.Timestep
    """

    def __init__(self, dimensions):
        self.dimensions = dimensions

        try:
            self.dimensions = np.asarray(self.dimensions, np.float32)
            if self.dimensions.shape != (6, ) and \
               self.dimensions.shape != (1, 6):
                raise ValueError('{} are not valid box dimensions'.format(
                self.dimensions))
            self.dimensions = self.dimensions.reshape(6, )
        except ValueError:
            raise ValueError(f'{self.dimensions} are not valid box dimensions')

    def __call__(self, ts):
        ts.dimensions = self.dimensions
        return ts
