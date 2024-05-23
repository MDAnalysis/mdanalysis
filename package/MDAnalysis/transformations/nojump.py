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
The consequence of this is that particles will diffuse across periodic boundaries when
needed. This unwrapping method is suitable as a preprocessing step to calculate
molecular diffusion, or more commonly to keep multi-domain proteins whole during trajectory
analysis.
The algorithm used is based on :footcite:p:`Kulke2022`.

.. autoclass:: NoJump

"""
import numpy as np
import warnings

from .base import TransformationBase
from ..due import due, Doi
from ..exceptions import NoDataError


class NoJump(TransformationBase):
    """
    Returns transformed coordinates for the given timestep so that an atom
    does not move more than half the periodic box size between two timesteps, and will move
    across periodic boundary edges. The algorithm used is based on :footcite:p:`Kulke2022`,
    equation B6 for non-orthogonal systems, so it is general to most applications where
    molecule trajectories should not "jump" from one side of a periodic box to another.
    
    Note that this transformation depends on a periodic box dimension being set for every
    frame in the trajectory, and that this box dimension can be transformed to an orthonormal
    unit cell. If not, an error is emitted. Since it is typical to transform all frames
    sequentially when unwrapping trajectories, warnings are emitted if non-sequential timesteps
    are transformed using NoJump.

    Example
    -------

    Suppose you had a molecular trajectory with a protein dimer that moved across a periodic
    boundary. This will normally appear to make the dimer split apart. This transformation
    uses the molecular motions between frames in a wrapped trajectory to create an unwrapped
    trajectory where molecules never move more than half a periodic box dimension between subsequent
    frames. Thus, the normal use case is to apply the transformation to every frame of a trajectory,
    and to do so sequentially.

    .. code-block:: python

        transformation = NoJump()
        u.trajectory.add_transformations(transformation)
        for ts in u.trajectory:
            print(ts.positions)

    In this case, ``ts.positions`` will return the NoJump unwrapped trajectory, which would keep
    protein dimers together if they are together in the inital timestep. The unwrapped trajectory may
    also be useful to compute diffusion coefficients via the Einstein relation.
    To reverse the process, wrap the coordinates again.

    Returns
    -------
    MDAnalysis.coordinates.timestep.Timestep

    References
    ----------
    .. footbibliography::

    """

    @due.dcite(
        Doi("10.1021/acs.jctc.2c00327"),
        description="Works through the orthogonal case for unwrapping, "
        "and proposes the non-orthogonal approach.",
        path=__name__,
    )
    def __init__(
        self,
        check_continuity=True,
        max_threads=None,
    ):
        # NoJump transforms are inherently unparallelizable, since
        # it depends on the previous frame's unwrapped coordinates
        super().__init__(max_threads=max_threads, parallelizable=False)
        self.prev = None
        self.old_frame = 0
        self.older_frame = "A"
        self.check_c = check_continuity

    def _transform(self, ts):
        L = ts.triclinic_dimensions
        if L is None:
            msg = f"Periodic box dimensions not provided at step {ts.frame}"
            raise NoDataError(msg)
        try:
            Linverse = np.linalg.inv(L)
        except np.linalg.LinAlgError:
            msg = f"Periodic box dimensions are not invertible at step {ts.frame}"
            raise NoDataError(msg)
        if ts.frame == 0:
            # We don't need to apply the transformation here. However, we need to
            # ensure we have the 0th frame coordinates in reduced form. We also need to
            # set an an appropriate value for self.older frame. This is so that on the
            # following frame we don't return early when we check
            # `self.older_frame != "A"`. If we return early, then the transformation is
            # not applied, and any jumps across boundaries that occur at that frame will
            # not be accounted for.
            self.prev = ts.positions @ Linverse
            self.old_frame = 0
            self.older_frame = -1
            return ts
        if (
            self.check_c
            and self.older_frame != "A"
            and (self.old_frame - self.older_frame) != (ts.frame - self.old_frame)
        ):
            warnings.warn(
                "NoJump detected that the interval between frames is unequal."
                "We are resetting the reference frame to the current timestep.",
                UserWarning,
            )
            self.prev = ts.positions @ Linverse
            self.old_frame = ts.frame
            self.older_frame = "A"
            return ts
        if self.check_c and np.abs(self.old_frame - ts.frame) != 1:
            warnings.warn(
                "NoJump transform is only accurate when positions"
                "do not move by more than half a box length."
                "Currently jumping between frames with a step of more than 1."
                "This might be fine, but depending on the trajectory stride,"
                "this might be inaccurate.",
                UserWarning,
            )
        # Convert into reduced coordinate space
        fcurrent = ts.positions @ Linverse
        fprev = self.prev  # Previous unwrapped coordinates in reduced box coordinates.
        # Calculate the new positions in reduced coordinate space (Equation B6 from
        # 10.1021/acs.jctc.2c00327). As it turns out, the displacement term can
        # be moved inside the round function in this coordinate space, as the
        # difference between wrapped and unwrapped coordinates is an integer.
        newpositions = fcurrent - np.round(fcurrent - fprev)
        # Convert back into real space
        ts.positions = newpositions @ L
        # Set things we need to save for the next frame.
        self.prev = newpositions  # Note that this is in reduced coordinate space.
        self.older_frame = self.old_frame
        self.old_frame = ts.frame

        return ts
