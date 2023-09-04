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
Fitting transformations --- :mod:`MDAnalysis.transformations.fit`
=================================================================

Translate and/or rotates the coordinates of a given trajectory to align
a given AtomGroup to a reference structure.

.. autoclass:: fit_translation

.. autoclass:: fit_rot_trans

"""
import numpy as np

from ..analysis import align
from ..lib.transformations import euler_from_matrix, euler_matrix

from .base import TransformationBase


class fit_translation(TransformationBase):
    """Translates a given AtomGroup so that its center of geometry/mass matches
    the respective center of the given reference. A plane can be given by the
    user using the option `plane`, and will result in the removal of
    the translation motions of the AtomGroup over that particular plane.

    Example
    -------
    Removing the translations of a given AtomGroup `ag` on the XY plane by
    fitting its center of mass to the center of mass of a reference `ref`:

    .. code-block:: python

        ag = u.select_atoms("protein")
        ref = mda.Universe("reference.pdb")
        transform = mda.transformations.fit_translation(ag, ref, plane="xy",
                                                        weights="mass")
        u.trajectory.add_transformations(transform)

    Parameters
    ----------
    ag : Universe or AtomGroup
       structure to translate, a
       :class:`~MDAnalysis.core.groups.AtomGroup` or a whole
       :class:`~MDAnalysis.core.universe.Universe`
    reference : Universe or AtomGroup
       reference structure, a :class:`~MDAnalysis.core.groups.AtomGroup` or a
       whole :class:`~MDAnalysis.core.universe.Universe`
    plane: str, optional
        used to define the plane on which the translations will be removed.
        Defined as a string of the plane.
        Suported planes are yz, xz and xy planes.
    weights : {"mass", ``None``} or array_like, optional
       choose weights. With ``"mass"`` uses masses as weights; with ``None``
       weigh each atom equally. If a float array of the same length as
       `ag` is provided, use each element of the `array_like` as a
       weight for the corresponding atom in `ag`.

    Returns
    -------
    MDAnalysis.coordinates.timestep.Timestep


    .. versionchanged:: 2.0.0
       The transformation was changed from a function/closure to a class
       with ``__call__``.
    .. versionchanged:: 2.0.0
       The transformation was changed to inherit from the base class for
       limiting threads and checking if it can be used in parallel analysis.
    """
    def __init__(self, ag, reference, plane=None, weights=None,
                 max_threads=None, parallelizable=True):
        super().__init__(max_threads=max_threads,
                         parallelizable=parallelizable)

        self.ag = ag
        self.reference = reference
        self.plane = plane
        self.weights = weights

        if self.plane is not None:
            axes = {'yz': 0, 'xz': 1, 'xy': 2}
            try:
                self.plane = axes[self.plane]
            except (TypeError, KeyError):
                raise ValueError(f'{self.plane} is not a valid plane') \
                                 from None
        try:
            if self.ag.atoms.n_residues != self.reference.atoms.n_residues:
                errmsg = (
                          f"{self.ag} and {self.reference} have mismatched"
                          f"number of residues"
                )

                raise ValueError(errmsg)
        except AttributeError:
            errmsg = (
                      f"{self.ag} or {self.reference} is not valid"
                      f"Universe/AtomGroup"
            )
            raise AttributeError(errmsg) from None
        self.ref, self.mobile = align.get_matching_atoms(self.reference.atoms,
                                                         self.ag.atoms)
        self.weights = align.get_weights(self.ref.atoms, weights=self.weights)
        self.ref_com = self.ref.center(self.weights)

    def _transform(self, ts):
        mobile_com = np.asarray(self.mobile.atoms.center(self.weights),
                                np.float32)
        vector = self.ref_com - mobile_com
        if self.plane is not None:
            vector[self.plane] = 0
        ts.positions += vector

        return ts


class fit_rot_trans(TransformationBase):
    """Perform a spatial superposition by minimizing the RMSD.

    Spatially align the group of atoms `ag` to `reference` by doing a RMSD
    fit.

    This fit works as a way to remove translations and rotations of a given
    AtomGroup in a trajectory. A plane can be given using the flag `plane`
    so that only translations and rotations in that particular plane are
    removed. This is useful for protein-membrane systems to where the membrane
    must remain in the same orientation.

    Note
    ----
    ``max_threads`` is set to 1 for this transformation
    with which it performs better.

    Example
    -------
    Removing the translations and rotations of a given AtomGroup `ag` on the XY plane
    by fitting it to a reference `ref`, using the masses as weights for the RMSD fit:

    .. code-block:: python

        ag = u.select_atoms("protein")
        ref = mda.Universe("reference.pdb")
        transform = mda.transformations.fit_rot_trans(ag, ref, plane="xy",
                                                      weights="mass")
        u.trajectory.add_transformations(transform)

    Parameters
    ----------
    ag : Universe or AtomGroup
       structure to translate and rotate, a
       :class:`~MDAnalysis.core.groups.AtomGroup` or a whole
       :class:`~MDAnalysis.core.universe.Universe`
    reference : Universe or AtomGroup
       reference structure, a :class:`~MDAnalysis.core.groups.AtomGroup` or a whole
       :class:`~MDAnalysis.core.universe.Universe`
    plane: str, optional
        used to define the plane on which the rotations and translations will be removed.
        Defined as a string of the plane. Supported planes are "yz", "xz" and "xy" planes.
    weights : {"mass", ``None``} or array_like, optional
       choose weights. With ``"mass"`` uses masses as weights; with ``None``
       weigh each atom equally. If a float array of the same length as
       `ag` is provided, use each element of the `array_like` as a
       weight for the corresponding atom in `ag`.

    Returns
    -------
    MDAnalysis.coordinates.timestep.Timestep


    .. versionchanged:: 2.0.0
       The transformation was changed from a function/closure to a class
       with ``__call__``.
    .. versionchanged:: 2.0.0
       The transformation was changed to inherit from the base class for
       limiting threads and checking if it can be used in parallel analysis.
    """
    def __init__(self, ag, reference, plane=None, weights=None,
                 max_threads=1, parallelizable=True):
        super().__init__(max_threads=max_threads,
                         parallelizable=parallelizable)

        self.ag = ag
        self.reference = reference
        self.plane = plane
        self.weights = weights

        if self.plane is not None:
            axes = {'yz': 0, 'xz': 1, 'xy': 2}
            try:
                self.plane = axes[self.plane]
            except (TypeError, KeyError):
                raise ValueError(f'{self.plane} is not a valid plane') \
                                 from None
        try:
            if self.ag.atoms.n_residues != self.reference.atoms.n_residues:
                errmsg = (
                    f"{self.ag} and {self.reference} have mismatched "
                    f"number of residues"
                )
                raise ValueError(errmsg)
        except AttributeError:
            errmsg = (
                      f"{self.ag} or {self.reference} is not valid "
                      f"Universe/AtomGroup"
            )
            raise AttributeError(errmsg) from None
        self.ref, self.mobile = align.get_matching_atoms(self.reference.atoms,
                                                         self.ag.atoms)
        self.weights = align.get_weights(self.ref.atoms, weights=self.weights)
        self.ref_com = self.ref.center(self.weights)
        self.ref_coordinates = self.ref.atoms.positions - self.ref_com

    def _transform(self, ts):
        mobile_com = self.mobile.atoms.center(self.weights)
        mobile_coordinates = self.mobile.atoms.positions - mobile_com
        rotation, dump = align.rotation_matrix(mobile_coordinates,
                                               self.ref_coordinates,
                                               weights=self.weights)
        vector = self.ref_com
        if self.plane is not None:
            matrix = np.r_[rotation, np.zeros(3).reshape(1, 3)]
            matrix = np.c_[matrix, np.zeros(4)]
            euler_angs = np.asarray(euler_from_matrix(matrix, axes='sxyz'),
                                    np.float32)
            for i in range(0, euler_angs.size):
                euler_angs[i] = (euler_angs[self.plane] if i == self.plane
                                 else 0)
            rotation = euler_matrix(euler_angs[0],
                                    euler_angs[1],
                                    euler_angs[2],
                                    axes='sxyz')[:3, :3]
            vector[self.plane] = mobile_com[self.plane]
        ts.positions = ts.positions - mobile_com
        ts.positions = np.dot(ts.positions, rotation.T)
        ts.positions = ts.positions + vector
        return ts
