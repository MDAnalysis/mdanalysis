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
No Jump Trajectory Unwrapping --- :mod:`MDAnalysis.transformations.nojump` 
=======================================================================================

Unwraps the trajectory such that atoms never move more than half a periodic box length.
The algorithm used is based on Kulke & Vermaas 2022, 10.1021/acs.jctc.2c00327

.. autoclass:: NoJump

"""
import numpy as np
import warnings

from .base import TransformationBase


class NoJump(TransformationBase):
    """
    Returns transformed coordinates for the given timestep so that an atom
    does not move more than half the periodic box size, and will not jump
    across the periodic boundary. The algorithm used is based on Kulke &
    Vermaas 2022, 10.1021/acs.jctc.2c00327, equation B8 for non-orthogonal
    systems.

    Example
    -------

    To calculate diffusion, one of the first steps is to compute a trajectory
    where the atoms do not cross the periodic boundary.

    .. code-block:: python

        transformation = NoJump()
        u.trajectory.add_transformations(transformation)
        for ts in u.trajectory:
            print(ts.positions)

    In this case, ``ts.positions`` will return the NoJump unwrapped trajectory.
    To reverse the process, wrap the coordinates again. 

    Returns
    -------
    MDAnalysis.coordinates.timestep.Timestep

    """

    def __init__(
        self,
        check_continuity=True,
        max_threads=None,
        # NoJump transforms are inherently unparallelizable, since
        # it depends on the previous frame's unwrapped coordinates
        parallelizable=False,
    ):
        super().__init__(max_threads=max_threads, parallelizable=parallelizable)
        self.Lprevinverse = np.eye(3)
        self.prevwrap = None
        self.prev = None
        self.old_frame = 0
        self.check_c = check_continuity

    def _transform(self, ts):
        if self.prev is None:
            self.prev = ts.positions.copy()
            self.prevwrap = ts.positions.copy()
            self.Lprevinverse = np.linalg.inv(ts.triclinic_dimensions)
            self.old_frame = ts.frame
            return ts
        if self.check_c and np.abs(self.old_frame - ts.frame) != 1:
            warnings.warn(
                "NoJump transform is only accurate when positions"
                "do not move by more than half a box length."
                "Currently jumping between frames with a step of more than 1."
                "This might be fine, but depending on the trajectory stride,"
                "this might be inaccurate.",
                Warning,
            )
        # Remember that @ is a shorthand for matrix multiplication.
        # np.matmul(a, b) is equivalent to a @ b
        alpha = ts.triclinic_dimensions @ self.Lprevinverse
        Linverse = np.linalg.inv(ts.triclinic_dimensions)
        current = ts.positions.copy()
        ts.positions = (alpha @ self.prev.T).T + current - (alpha @ self.prevwrap.T).T - (ts.triclinic_dimensions @ np.round(Linverse @ current.T - self.Lprevinverse @ self.prevwrap.T)).T
        self.Lprevinverse = Linverse
        self.prevwrap = current
        self.prev = ts.positions.copy()
        self.old_frame = ts.frame

        return ts
