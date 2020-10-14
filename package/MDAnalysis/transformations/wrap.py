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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""\
Wrap/unwrap transformations --- :mod:`MDAnalysis.transformations.wrap`
======================================================================

Wrap/unwrap the atoms of a given AtomGroup in the unit cell. :class:`wrap`
translates the atoms back in the unit cell. :class:`unwrap` translates the
atoms of each molecule so that bons don't split over images. :class:`nojump`
removes jumps across the periodic boundaries when iterating along a trajectory.

.. autoclass:: wrap

.. autoclass:: unwrap

.. autoclass:: nojump
"""

from ..lib.distances import apply_PBC
from ..lib.mdamath import triclinic_vectors
from ..exceptions import NoDataError
import numpy as np


class wrap(object):
    """
    Shift the contents of a given AtomGroup back into the unit cell. ::

    Wrapping with `compound`='atoms':
       +-----------+              +-----------+
       |           |              |           |
       |           | 3 6          | 3 6       |
       |           | ! !          | ! !       |
       |         1-|-2-5-8  ->    |-2-5-8   1-|
       |           | ! !          | ! !       |
       |           | 4 7          | 4 7       |
       |           |              |           |
       +-----------+              +-----------+

    Wrapping with `compound`='fragments':
       +-----------+              +-----------+
       |           |              |           |
       |           | 3 6          | 3 6       |
       |           | ! !          | ! !       |
       |         1-|-2-5-8  ->  1-|-2-5-8     |
       |           | ! !          | ! !       |
       |           | 4 7          | 4 7       |
       |           |              |           |
       +-----------+              +-----------+

    Example
    -------

    .. code-block:: python

        ag = u.atoms
        # Wrapping to the compact unit cell representation
        transform = mda.transformations.wrap(ag, compact=True)
        u.trajectory.add_transformations(transform)

    Parameters
    ----------

    ag: Atomgroup
        Atomgroup to be wrapped in the unit cell
    kwargs : optional
        Arguments to be passed to
        :meth:`MDAnalysis.core.groups.AtomGroup.wrap`.

    Notes
    -----
    This is a wrapper around :meth:`MDAnalysis.core.groups.AtomGroup.wrap`.
    Refer to its documentation on additional `kwargs` parameters. Most
    usefully, `compound` allows the selection of how different sub-groups are
    wrapped (per atom, residue, fragment, molecule, etc.);
    `compact` sets wrapping to a compact unit cell representation;
    `center` allows the choice of the center-of-mass vs the center-of-geometry
    as the reference center for each compound. The same translation is applied
    to all atoms within a compound, meaning it will not be broken by the shift.
    This might however mean that not all atoms from the compound are
    inside the unit cell, but rather the center of the compound is.

    Returns
    -------
    MDAnalysis.coordinates.base.Timestep

    See Also
    --------
    :meth:`MDAnalysis.core.groups.AtomGroup.wrap`


    .. versionchanged:: 2.0.0
        The transformation was changed from a function/closure to a class
        with ``__call__``. It now also fully exposes the interface of
        :meth:`MDAnalysis.core.groups.AtomGroup.wrap`
    """
    def __init__(self, ag, **kwargs):
        self.ag = ag
        self.kwargs = kwargs

    def __call__(self, ts):
        self.ag.wrap(**self.kwargs)
        return ts


class unwrap(object):
    """
    Move all atoms in an AtomGroup so that bonds don't split over images

    Atom positions are modified in place.

    This function is most useful when atoms have been packed into the primary
    unit cell, causing breaks mid molecule, with the molecule then appearing
    on either side of the unit cell. This is problematic for operations
    such as calculating the center of mass of the molecule. ::

       +-----------+     +-----------+
       |           |     |           |
       | 6       3 |     |         3 | 6
       | !       ! |     |         ! | !
       |-5-8   1-2-| ->  |       1-2-|-5-8
       | !       ! |     |         ! | !
       | 7       4 |     |         4 | 7
       |           |     |           |
       +-----------+     +-----------+

    Example
    -------

    .. code-block:: python

        ag = u.atoms
        transform = mda.transformations.unwrap(ag)
        u.trajectory.add_transformations(transform)

    Parameters
    ----------
    atomgroup : AtomGroup
        The :class:`MDAnalysis.core.groups.AtomGroup` to work with.
        The positions of this are modified in place.
    kwargs : optional
        Arguments to be passed to
        :meth:`MDAnalysis.core.groups.AtomGroup.unwrap`.

    Notes
    -----
    This is a wrapper around :meth:`MDAnalysis.core.groups.AtomGroup.unwrap`.
    Refer to its documentation on additional `kwargs` parameters. Most
    usefully, `compound` allows the selection of how different sub-groups are
    wrapped (as a whole, per residue, fragment, molecule, etc.);
    `reference` allows the choice of the center-of-mass vs the
    center-of-geometry as the reference center for each compound.

    Returns
    -------
    :class:`MDAnalysis.coordinates.base.Timestep`

    See Also
    --------
    :func:`MDAnalysis.core.groups.AtomGroup.unwrap`


    .. versionchanged:: 2.0.0
        The transformation was changed from a function/closure to a class
        with ``__call__``. It now also fully exposes the interface of
        :meth:`MDAnalysis.core.groups.AtomGroup.unwrap`
    """
    def __init__(self, ag, **kwargs):
        self.ag = ag
        self.kwargs = kwargs

    def __call__(self, ts):
        self.ag.unwrap(**self.kwargs)
        return ts


class nojump(object):
    """
    Remove jumps across the periodic boundaries in consecutive frames.

    Jumps are identified from translations greater than half the length of each
    box vector. Selective jump removal for x, y, and z is possible.

    This function is based on consecutive frame displacement vectors. It will
    keep molecules whole provided they were whole in the first iteration frame.

    Example
    -------

    .. code-block:: python

        ag = u.atoms
        # starts from a configuration with whole fragments.
        transform = mda.transformations.nojump(ag, initial_workflow=unwrap(ats))
        u.trajectory.add_transformations(transform)

    Parameters
    ----------
    atomgroup : AtomGroup
        The :class:`MDAnalysis.core.groups.AtomGroup` to work with.
        The positions of this are modified in place.
    x, y, z : bool, optional
        Over which dimension(s) to remove jumps.
    initial_workflow : callable or list, optional
        Transform(s) to apply when starting or restarting iteration. Most
        usefully, passing :class:`unwrap` makes compounds whole for the first
        frame, and the nojump operation takes care of keeping them whole for
        the remaining frames.
    refocus : bool or int
        Whether to correct each frame for position drift due to successive
        vector addition. If set to a nonzero int, defines every how many frames
        to correct, representing a tradeoff between performance and accuracy.
    backend : {'serial', 'OpenMP'}, optional
        Keyword selecting the type of acceleration in subcalls to
        :func:`MDAnalysis.lib.distances.apply_PBC`.

    Notes
    -----
    Because a no-jump treatment preserves atom neighborhoods, an unwrap
    transformation at the first frame (specified with `initial_workflow`)
    ensures that fragments are kept whole throughout the transformed
    trajectory, and avoids the need for additional unwrapping transformations.
    However, if :class:`nojump` only operates over some of the dimensions,
    fragments can become broken over the untreated ones.
    Because it needs frame-sequential position information, :class:`nojump`
    will reset the transformation if it detects a backward jump in frames.
    Multi-frame forward jumps do not trigger such a transformation reset, but
    at very large frame skips displacements between treated frames may become
    larger than half the box vectors, yielding wrong results.

    Returns
    -------
    :class:`MDAnalysis.coordinates.base.Timestep`

    See Also
    --------
    :func:`MDAnalysis.core.groups.AtomGroup.unwrap`


    .. versionadded:: 2.0.0
    """
    def __init__(self, ag, x=True, y=True, z=True,
                 initial_workflow=None, refocus=True, backend='serial'):
        self.ag = ag
        self.prev_frame = None
        self.prev_cdx = None
        self.transformed_cdx = None
        self.dims = [i for i, d in enumerate((x, y, z)) if d]
        self.initial_workflow = []
        if initial_workflow:
            if callable(initial_workflow):
                self.initial_workflow = [initial_workflow]
            else:
                self.initial_workflow = initial_workflow
        self.refocus = int(refocus)
        self.backend = backend

    def __call__(self, ts):
        # If all dims are False, this is a null transform
        if not self.dims:
            return ts

        # Shortcut when the same frame is loaded repeatedly
        if ts.frame == self.prev_frame:
            self.ag.positions = self.transformed_cdx
            return ts

        # The case of frame 0 or an iteration reset
        if self.prev_cdx is None or ts.frame - self.prev_frame < 0:
            for transform in self.initial_workflow:
                ts = transform(ts)
            self.prev_cdx = self.ag.positions
            self.transformed_cdx = self.ag.positions
            self.prev_frame = ts.frame
            self.accumulated_frames = 0
            return ts

        self.accumulated_frames += 1
        vecs = self.ag.positions - self.prev_cdx
        self.prev_cdx = self.ag.positions

        vecs_fixed = apply_PBC(vecs, box=ts.dimensions, center=[0., 0., 0.],
                               backend=self.backend)

        vecs[:, self.dims] = vecs_fixed[:, self.dims]
        self.transformed_cdx += vecs
        if self.refocus and self.accumulated_frames >= self.refocus:
            # refocuses self.transformed_cdx to actual shifts of primary
            # cell coordinates.
            self._refocus()
        self.ag.positions = self.transformed_cdx
        self.prev_frame = ts.frame
        return ts

    def _refocus(self):
        corrective_shift = apply_PBC(self.ag.positions - self.transformed_cdx,
                                     box=self.ag.dimensions,
                                     center=[0., 0., 0.],
                                     backend=self.backend)
        self.transformed_cdx += corrective_shift
        self.accumulated_frames = 0
