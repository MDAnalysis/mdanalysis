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
Trajectory rotation --- :mod:`MDAnalysis.transformations.rotate`
================================================================

Rotates the coordinates by a given angle arround an axis formed by a direction
and a point.

.. autoclass:: rotateby

"""
import numpy as np
from functools import partial

from ..lib.transformations import rotation_matrix
from ..lib.util import get_weights

from .base import TransformationBase


class rotateby(TransformationBase):
    '''
    Rotates the trajectory by a given angle on a given axis. The axis is defined by
    the user, combining the direction vector and a point. This point can be the center
    of geometry or the center of mass of a user defined AtomGroup, or an array defining
    custom coordinates.

    Note
    ----
    ``max_threads`` is set to 1 for this transformation
    with which it performs better.

    Examples
    --------

    e.g. rotate the coordinates by 90 degrees on a axis formed by the [0,0,1] vector and
    the center of geometry of a given AtomGroup:

    .. code-block:: python

        from MDAnalysis import transformations

        ts = u.trajectory.ts
        angle = 90
        ag = u.atoms
        d = [0,0,1]
        rotated = transformations.rotate.rotateby(angle, direction=d, ag=ag)(ts)

    e.g. rotate the coordinates by a custom axis:

    .. code-block:: python

        from MDAnalysis import transformations

        ts = u.trajectory.ts
        angle = 90
        p = [1,2,3]
        d = [0,0,1]
        rotated = transformations.rotate.rotateby(angle, direction=d, point=p)(ts)

    Parameters
    ----------
    angle: float
        rotation angle in degrees
    direction: array-like
        vector that will define the direction of a custom axis of rotation from the
        provided point. Expected shapes are (3, ) or (1, 3). 
    ag: AtomGroup, optional
        use the weighted center of an AtomGroup as the point from where the rotation axis
        will be defined. If no AtomGroup is given, the `point` argument becomes mandatory
    point: array-like, optional
        list of the coordinates of the point from where a custom axis of rotation will
        be defined. Expected shapes are (3, ) or (1, 3). If no point is given, the
        `ag` argument becomes mandatory.
    weights: {"mass", ``None``} or array_like, optional
        define the weights of the atoms when calculating the center of the AtomGroup.
        With ``"mass"`` uses masses as weights; with ``None`` weigh each atom equally.
        If a float array of the same length as `ag` is provided, use each element of
        the `array_like` as a weight for the corresponding atom in `ag`. Default is 
        None.
    wrap: bool, optional
        If `True`, all the atoms from the given AtomGroup will be moved to the unit cell
        before calculating the center of mass or geometry. Default is `False`, no changes
        to the atom coordinates are done before calculating the center of the AtomGroup. 

    Returns
    -------
    MDAnalysis.coordinates.base.Timestep

    Warning
    -------
    Wrapping/unwrapping the trajectory or performing PBC corrections may not be possible 
    after rotating the trajectory.


    .. versionchanged:: 2.0.0
        The transformation was changed from a function/closure to a class
        with ``__call__``.
    .. versionchanged:: 2.0.0
       The transformation was changed to inherit from the base class for
       limiting threads and checking if it can be used in parallel analysis.
    '''
    def __init__(self,
                 angle,
                 direction,
                 point=None,
                 ag=None,
                 weights=None,
                 wrap=False,
                 max_threads=1,
                 parallelizable=True):
        super().__init__(max_threads=max_threads,
                         parallelizable=parallelizable)

        self.angle = angle
        self.direction = direction
        self.point = point
        self.ag = ag
        self.weights = weights
        self.wrap = wrap

        self.angle = np.deg2rad(self.angle)
        try:
            self.direction = np.asarray(self.direction, np.float32)
            if self.direction.shape != (3, ) and \
               self.direction.shape != (1, 3):
                raise ValueError('{} is not a valid direction'
                                 .format(self.direction))
            self.direction = self.direction.reshape(3, )
        except ValueError:
            raise ValueError(f'{self.direction} is not a valid direction') \
                             from None
        if self.point is not None:
            self.point = np.asarray(self.point, np.float32)
            if self.point.shape != (3, ) and self.point.shape != (1, 3):
                raise ValueError('{} is not a valid point'.format(self.point))
            self.point = self.point.reshape(3, )
        elif self.ag:
            try:
                self.atoms = self.ag.atoms
            except AttributeError:
                raise ValueError(f'{self.ag} is not an AtomGroup object') \
                                from None
            else:
                try:
                    self.weights = get_weights(self.atoms,
                                               weights=self.weights)
                except (ValueError, TypeError):
                    errmsg = ("weights must be {'mass', None} or an iterable "
                              "of the same size as the atomgroup.")
                    raise TypeError(errmsg) from None
            self.center_method = partial(self.atoms.center,
                                         self.weights,
                                         pbc=self.wrap)
        else:
            raise ValueError('A point or an AtomGroup must be specified')

    def _transform(self, ts):
        if self.point is None:
            position = self.center_method()
        else:
            position = self.point
        matrix = rotation_matrix(self.angle, self.direction, position)
        rotation = matrix[:3, :3].T
        translation = matrix[:3, 3]
        ts.positions = np.dot(ts.positions, rotation)
        ts.positions += translation
        return ts
