# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2020 The MDAnalysis Development Team and contributors
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

"""
HELANAL --- analysis of protein helices
=======================================

:Author: Lily Wang
:Year: 2020
:Copyright: GNU Public License v3

.. versionadded:: 1.0.0

This module contains code to analyse protein helices using the 
`HELANAL <http://nucleix.mbu.iisc.ernet.in/helanal/helanal/helanal.html>`_ algorithm
([Bansal2000]_ , [Sugeta1967]_ ).

`HELANAL <http://nucleix.mbu.iisc.ernet.in/helanal/helanal/helanal.html>`_ 
quantifies the geometry of helices in proteins on the basis of their Cα
atoms. It can determine local structural features such as the local helical twist and
rise, virtual torsion angle, local helix origins and bending angles between
successive local helix axes. 

Example use
-----------

You can pass in a single selection::

    import MDAnalysis as mda
    from MDAnalysis.tests.datafiles import PSF, DCD
    from MDAnalysis.analysis import helix_analysis as hel
    u = mda.Universe(PSF, DCD)
    helanal = hel.HELANAL(u, select='name CA and resnum 161-187')
    helanal.run()

    print(helanal.summary)

Alternatively, you can analyse several helices at once by passing 
in multiple selection strings::

    helanal2 = hel.HELANAL(u, select=('name CA and resnum 100-160',
                                      'name CA and resnum 200-230'))

"""
from __future__ import division, absolute_import

import warnings
import numpy as np

import MDAnalysis as mda
from ..lib import util, mdamath
from .base import AnalysisBase


def pdot(a, b):
    """Pairwise dot product.

    ``a`` must be the same shape as ``b``.

    Parameters
    ----------
    a: :class:`numpy.ndarray` of shape (N, M)
    b: :class:`numpy.ndarray` of shape (N, M)

    Returns
    -------
    :class:`numpy.ndarray` of shape (N,)
    """
    return np.einsum('ij,ij->i', a, b)


def pnorm(a):
    """Euclidean norm of each vector in a matrix

    Parameters
    ----------
    a: :class:`numpy.ndarray` of shape (N, M)

    Returns
    -------
    :class:`numpy.ndarray` of shape (N,)
    """
    return pdot(a, a)**0.5


def vector_of_best_fit(coordinates):
    """Fit vector through the centered coordinates.

    Parameters
    ----------
    coordinates: :class:`numpy.ndarray` of shape (N, 3)

    Returns
    -------
    :class:`numpy.ndarray` of shape (3,)
        Vector of best fit.
    """
    centered = coordinates - coordinates.mean(axis=0)
    Mt_M = np.matmul(centered.T, centered)
    u, s, vh = np.linalg.linalg.svd(Mt_M)
    vector = vh[0]

    # does vector face first residue?
    angle = mdamath.angle(coordinates[0], vector)
    if angle > np.pi/2 or angle < -np.pi/2:
        vector *= -1
    return vector


def local_screw(best_fit, ref_axis, rotation_vectors):
    """
    Get local screw angles.

    Parameters
    ----------
    best_fit: :class:`numpy.ndarray` of shape (3,)
        Vector of best fit. Tilts are calculated perpendicular to 
        this axis.
    ref_axis: :class:`numpy.ndarray` of shape (3,)
        Reference length-wise axis. One of the reference vectors is 
        orthogonal to this axis.
    rotation_vectors: :class:`numpy.ndarray` of shape (N, 3)
        array of vectors representating the local rotation of each helix.

    Returns
    -------
    :class:`numpy.ndarray` of shape (N,)
        Array of screw angles.
    """
    ref_2 = np.cross(ref_axis, best_fit)
    ref_1 = np.cross(-ref_2, best_fit)  # orthogonal to best_fit and ref_2
    refs = np.array([ref_1, ref_2])  # (2, 3)

    ref_norms = pnorm(refs)
    rot_norms = pnorm(rotation_vectors)

    # angles
    cos = np.matmul(refs, rotation_vectors.T)/np.outer(ref_norms, rot_norms)
    screw_angles, alt_angles = np.arccos(np.clip(cos, -1, 1))  # (2, n_vec)

    # the effect of the below is to move angles into the (pi/4, 3pi/4 domain)
    quarter_pi, half_pi = np.pi/4, np.pi/2
    # replace low angles
    low = screw_angles < quarter_pi
    screw_angles[low] = np.where(alt_angles[low] < half_pi,
                                 half_pi - alt_angles[low],
                                 alt_angles[low] - half_pi)
    # replace high angles
    high = screw_angles > 3*quarter_pi
    screw_angles[high] = np.where(alt_angles[high] < half_pi,
                                  half_pi + alt_angles[high],
                                  3*half_pi - alt_angles[high])

    ortho = np.cross(ref_1, rotation_vectors)
    cos_theta = np.matmul(best_fit, ortho.T)/(pnorm(ortho)*ref_norms[0])
    angle_to_best_fit = np.arccos(np.clip(cos_theta, -1, 1))
    screw_angles[angle_to_best_fit < half_pi] *= -1
    return np.rad2deg(screw_angles)  # (n_vec,)


def helix_analysis(positions, ref_axis=[0, 0, 1]):
    r"""
    Calculate helix properties from atomic coordinates.

    Each property is calculated from a sliding window of 4 atoms, 
    from i to i+3.

    Parameters
    ----------
    positions: :class:`numpy.ndarray` of shape (N, 3)
        Atomic coordinates.
    ref_axis: array-like of length 3, optional
        The reference axis used to calculate the tilt of the vector of best fit, 
        and the local screw angles.

    Returns
    -------
    dict with the following keys:
        local_twists: array, shape (N-3,)
            local twist angle from atom i+1 to i+2
        residues_per_turn: array, shape (N-3,)
            number of residues per turn, based on local_twist
        local_axes:  array, shape (N-3, 3)
            the length-wise helix axis of the local window
        local_bends: array, shape (N-6,)
            the angles between local helix angles, 3 windows apart
        local_heights: array, shape (N-3,)
            the rise of each local helix
        local_rotation_vectors: array, shape (N-2, 3)
            the unit vector from each local origin to atom i+1
        local_origins: array, shape (N-2, 3)
            the projected origin for each helix
        all_bends: array, shape (N-3, N-3)
            angles between each local axis
        global_axes: array, shape (3,)
            vector of best fit through origins
        local_screw: array, shape (N-2,)
            angles of each rotation vector to normal plane of global axes
    """

    #          ^               ^
    #           \             / bi
    #            \           /
    #         CA_i+2 <----- CA_i+1
    #         /    \       /   ^
    #        /    r \     /     \
    #     V /        \ θ /       \
    #      /          \ /       CA_i
    #     v           origin
    #   CA_i+3
    #
    # V: vectors
    # bi: bisectors
    # θ: local_twists
    # origin: origins
    # local_axes: in the plane of the screen. Orthogonal to the bisectors

    vectors = positions[1:] - positions[:-1]  # (n_res-1, 3)
    bisectors = vectors[:-1] - vectors[1:]  # (n_res-2, 3)
    binorms = pnorm(bisectors)  # (n_res-2,)
    adjacent_mag = binorms[:-1] * binorms[1:]  # (n_res-3,)

    # find angle between bisectors for twist and n_residue/turn
    cos_theta = pdot(bisectors[:-1], bisectors[1:])/adjacent_mag
    cos_theta = np.clip(cos_theta, -1, 1)
    twists = np.arccos(cos_theta)  # (n_res-3,)
    local_twists = np.rad2deg(twists)
    residues_per_turn = 2*np.pi / twists

    # find normal to bisectors for local axes
    cross_bi = np.cross(bisectors[:-1], bisectors[1:])  # (n_res-3, 3)
    local_axes = (cross_bi.T / pnorm(cross_bi)).T  # (n_res-3, 3)

    # find angles between axes for bends
    bend_theta = np.matmul(local_axes, local_axes.T)  # (n_res-3, n_res-3)
    bend_matrix = np.rad2deg(np.arccos(np.clip(bend_theta, -1, 1)))
    # local bends are between axes 3 windows apart
    local_bends = np.diagonal(bend_matrix, offset=3)  # (n_res-6,)

    # radius of local cylinder
    radii = (adjacent_mag**0.5) / (2*(1.0-cos_theta))  # (n_res-3,)
    # height of local cylinder
    heights = pdot(vectors[1:-1], local_axes)  # (n_res-3,)

    local_rotation_vectors = (bisectors.T/binorms).T  # (n_res-2, 3)

    # get origins by subtracting radius from atom i+1
    origins = positions[1:-1][:]  # (n_res-2, 3)
    origins[:-1] -= (radii*local_rotation_vectors[:-1].T).T
    # subtract radius from atom i+2 in last one
    origins[-1] -= radii[-1]*local_rotation_vectors[-1]

    helix_axes = vector_of_best_fit(origins)
    screw = local_screw(helix_axes, np.asarray(ref_axis),
                        local_rotation_vectors)

    results = {'local_twists': local_twists,
               'residues_per_turn': residues_per_turn,
               'local_axes': local_axes,
               'local_bends': local_bends,
               'local_heights': heights,
               'local_rotation_vectors': local_rotation_vectors,
               'local_origins': origins,
               'all_bends': bend_matrix,
               'global_axes': helix_axes,
               'local_screw': screw}
    return results


class HELANAL(AnalysisBase):
    r"""
    Perform HELANAL helix analysis on your trajectory.

    Parameters
    ----------
    universe: Universe or AtomGroup
        The Universe or AtomGroup to apply the analysis to.
    select: str or iterable of str, optional
        The selection string to create an atom selection that the HELANAL
        analysis is applied to. Note that HELANAL is designed to work on the 
        alpha-carbon atoms of protein residues. If you pass in multiple 
        selections, the selections will be analysed separately.
    ref_axis: array-like of length 3, optional
        The reference axis used to calculate the tilt of the vector of best fit, 
        and the local screw angles.
    flatten_single_helix: bool, optional
        Whether to flatten results if only one selection is passed.
    verbose : bool, optional
        Turn on more logging and debugging.

    Attributes
    ----------
    local_twists: array or list of arrays
        The local twist angle from atom i+1 to i+2. 
        Each array has shape (n_frames, n_residues-3)
    residues_per_turn: array or list of arrays
        Number of residues per turn, based on local_twist. 
        Each array has shape (n_frames, n_residues-3)
    local_axes: array or list of arrays
        The length-wise helix axis of the local window. 
        Each array has shape (n_frames, n_residues-3, 3)
    local_heights: array or list of arrays
        The rise of each local helix. 
        Each array has shape (n_frames, n_residues-3)
    local_rotation_vectors: array or list of arrays
        The unit vector from each local origin to atom i+1.
        Each array has shape (n_frames, n_residues-2, 3)
    local_origins: array or list of arrays
        The projected origin for each helix.
        Each array has shape (n_frames, n_residues-2, 3)
    local_screw: array or list of arrays
        The local screw angle for each helix.
        Each array has shape (n_frames, n_residues-2)
    local_bends: array or list of arrays
        The angles between local helix axes, 3 windows apart. 
        Each array has shape (n_frames, n_residues-6)
    all_bends: array or list of arrays
        The angles between local helix axes.
        Each array has shape (n_frames, n_residues-3, n_residues-3)
    global_axes: array or list of arrays
        The length-wise axis for the overall helix. 
        Each array has shape (n_frames, 3)
    global_tilts: array or list of arrays
        The angle between the global axis and the reference axis.
        Each array has shape (n_frames,)
    summary: dict or list of dicts
        Summary of stats for each property: the mean, the sample 
        standard deviation, and the mean absolute deviation.
    """

    # shapes of properties from each frame, relative to n_residues
    attr_shapes = {
        'local_twists': (-3,),
        'local_bends': (-6,),
        'local_heights': (-3,),
        'residues_per_turn': (-3,),
        'local_origins': (-2, 3),
        'local_axes': (-3, 3),
        'local_rotation_vectors': (-2, 3),
        'local_screw': (-2,),
    }

    @classmethod
    def from_SecondaryStructure(cls, secondary_structure, select='name CA',
                                ref_axis=[0, 0, 1], verbose=False,
                                flatten_single_helix=True):
        """Construct a HELANAL instance from a SecondaryStructure calculation (e.g. DSSP).

        Parameters
        ----------
        secondary_structure: instance of SecondaryStructureBase class
            The secondary structure calculation to use. Must have been run()
            already.
        select: str, optional
            The selection string to select atoms from each residue found. Note 
            that unlike normal HELANAL creation, this cannot be an iterable of 
            multiple selections.
        ref_axis: array-like of length 3, optional
            The reference axis used to calculate the tilt of the vector of best fit, 
            and the local screw angles.
        flatten_single_helix: bool, optional
            Whether to flatten results if only one selection is passed.
        verbose : bool, optional
            Turn on more logging and debugging.

        Returns
        -------
        :class:`~MDAnalysis.analysis.helix_analysis.HELANAL`
        """
        try:
            ss = secondary_structure.simple_mode
        except AttributeError:
            raise ValueError('The SecondaryStructure instance must be run '
                             'before passing to HELANAL')

        helices = secondary_structure.residues[ss == 'Helix']
        continuous = util.group_consecutive_integers(helices.resindices)

        # must have at least 9 continuous residues for helanal
        selections = []
        universe = secondary_structure._universe
        all_res = universe.residues
        for resindices in continuous:
            if len(hel) >= 9:
                res = all_res[resindices]
                sel = 'resid {}:{} and {}'.format(res[0].resid,
                                                  res[-1].resid,
                                                  select)
                selections.append(sel)
        if not selections:
            raise ValueError('Could not find any helices with at least '
                             '9 continuous residues. Consider running '
                             'HELANAL manually')

        return cls(universe, select=selections, ref_axis=ref_axis,
                   verbose=verbose, flatten_single_helix=flatten_single_helix)

    def __init__(self, universe, select='name CA', ref_axis=[0, 0, 1],
                 verbose=False, flatten_single_helix=True):
        super(HELANAL, self).__init__(universe.universe.trajectory,
                                      verbose=verbose)
        selections = util.asiterable(select)
        self.atomgroups = [universe.select_atoms(s) for s in selections]
        for s, ag in zip(selections, self.atomgroups):
            ids, counts = np.unique(ag.resindices, return_counts=True)
            if np.any(counts > 1):
                dup = ', '.join(map(str, ids[counts > 1]))
                warnings.warn('Your selection {} includes multiple atoms '
                              'for residues with these resindices: {}.'
                              'HELANAL is designed to work on one carbon-alpha '
                              'per residue.'.format(s, dup))
        self.ref_axis = np.asarray(ref_axis)
        self._flatten = flatten_single_helix

    def _zeros_per_frame(self, *dims, n_positions=0):
        """Create zero arrays where first 2 dims are n_frames, n_values"""
        first = dims[0]
        return np.zeros((self.n_frames, first+n_positions, *dims[1:]),
                         dtype=np.float64)

    def _prepare(self):
        n_res = [len(ag) for ag in self.atomgroups]

        for key, dims in self.attr_shapes.items():
            empty = [self._zeros_per_frame(
                *dims, n_positions=n) for n in n_res]
            setattr(self, key, empty)

        self.global_axes = [self._zeros_per_frame(3) for n in n_res]
        self.all_bends = [self._zeros_per_frame(n-3, n-3) for n in n_res]

    def _single_frame(self):
        _f = self._frame_index
        for i, ag in enumerate(self.atomgroups):
            results = helix_analysis(ag.positions)
            for key, value in results.items():
                attr = getattr(self, key)
                attr[i][_f] = value

    def _conclude(self):
        # compute tilt of global axes
        self.global_tilts = []
        norm_ref = (self.ref_axis**2).sum() ** 0.5
        for axes in self.global_axes:
            cos = np.matmul(self.ref_axis, axes.T) / (pnorm(axes)*norm_ref)
            cos = np.clip(cos, -1.0, 1.0)
            self.global_tilts.append(np.rad2deg(np.arccos(cos)))

        global_attrs = ['global_axes', 'global_tilts', 'all_bends']
        attrnames = list(self.attr_shapes.keys()) + global_attrs
        # summarise
        self.summary = []
        for i in range(len(self.atomgroups)):
            stats = {}
            for name in attrnames:
                attr = getattr(self, name)
                mean = attr[i].mean(axis=0)
                dev = np.abs(attr[i]-mean)
                stats[name] = {'mean': mean,
                               'sample_sd': attr[i].std(axis=0, ddof=1),
                               'abs_dev': dev.mean(axis=0)}
            self.summary.append(stats)

        # flatten?
        if len(self.atomgroups) == 1 and self._flatten:
            for name in attrnames + ['summary']:
                attr = getattr(self, name)
                setattr(self, name, attr[0])

    def universe_from_origins(self):
        """
        Create MDAnalysis Universe from the local origins.

        Returns
        -------
        Universe or list of Universes
        """
        try:
            origins = self.local_origins
        except AttributeError:
            raise ValueError('Call run() before universe_from_origins')

        if not isinstance(origins, list):
            origins = [origins]

        universe = []
        for xyz in origins:
            n_res = xyz.shape[1]
            u = mda.Universe.empty(n_res, n_residues=n_res,
                                   atom_resindex=np.arange(n_res),
                                   trajectory=True).load_new(xyz)
            universe.append(u)
        if not isinstance(self.local_origins, list):
            universe = universe[0]
        return universe
