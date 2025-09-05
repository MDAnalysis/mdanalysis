# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2020 The MDAnalysis Development Team and contributors
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
HELANAL --- analysis of protein helices
=======================================

:Author: Lily Wang
:Year: 2020
:Copyright: Lesser GNU Public License v2.1+

.. versionadded:: 2.0.0

This module contains code to analyse protein helices using the
HELANAL_ algorithm
([Bansal2000]_ , [Sugeta1967]_ ).

HELANAL_ quantifies the geometry of helices in proteins on the basis of their
Cα atoms. It can determine local structural features such as the local
helical twist and rise, virtual torsion angle, local helix origins and
bending angles between successive local helix axes.

.. _HELANAL: https://pubmed.ncbi.nlm.nih.gov/10798526/

.. [Sugeta1967] Sugeta, H. and Miyazawa, T. 1967. General method for
   calculating helical parameters of polymer chains from bond lengths, bond
   angles and internal rotation angles. *Biopolymers* 5 673 - 679

.. [Bansal2000] Bansal M, Kumar S, Velavan R. 2000.
   HELANAL - A program to characterise helix geometry in proteins.
   *J Biomol Struct Dyn.*  17(5):811-819.


Example use
-----------

You can pass in a single selection::

    import MDAnalysis as mda
    from MDAnalysis.tests.datafiles import PSF, DCD
    from MDAnalysis.analysis import helix_analysis as hel
    u = mda.Universe(PSF, DCD)
    helanal = hel.HELANAL(u, select='name CA and resnum 161-187')
    helanal.run()

All computed properties are available in ``.results``::

    print(helanal.results.summary)

Alternatively, you can analyse several helices at once by passing
in multiple selection strings::

    helanal2 = hel.HELANAL(u, select=('name CA and resnum 100-160',
                                      'name CA and resnum 200-230'))

The :func:`helix_analysis` function will carry out helix analysis on
atom positions, treating each row of coordinates as an alpha-carbon
equivalent::

    hel_xyz = hel.helix_analysis(u.atoms.positions, ref_axis=[0, 0, 1])

    
Classes
-------

.. autoclass:: HELANAL


Functions
---------

.. autofunction:: helix_analysis

.. autofunction:: vector_of_best_fit

.. autofunction:: local_screw_angles


"""

import warnings
import numpy as np

import MDAnalysis as mda
from ..lib import util, mdamath
from .base import AnalysisBase


def vector_of_best_fit(coordinates):
    """Fit vector through the centered coordinates,
    pointing to the first coordinate (i.e. upside-down).

    Parameters
    ----------
    coordinates : :class:`numpy.ndarray` of shape (N, 3)

    Returns
    -------
    :class:`numpy.ndarray` of shape (3,)
        Vector of best fit.
    """
    centered = coordinates - coordinates.mean(axis=0)
    Mt_M = np.matmul(centered.T, centered)
    u, s, vh = np.linalg.svd(Mt_M)
    vector = vh[0]

    # does vector face first local helix origin?
    angle = mdamath.angle(centered[0], vector)
    if angle > np.pi / 2:
        vector *= -1
    return vector


def local_screw_angles(global_axis, ref_axis, helix_directions):
    """
    Cylindrical azimuth angles between the local direction vectors,
    as projected onto the cross-section of the helix, from (-pi, pi].
    The origin (angle=0) is set to the plane of global_axis and ref_axis.

    Parameters
    ----------
    global_axis : :class:`numpy.ndarray` of shape (3,)
        Vector of best fit. Screw angles are calculated perpendicular to
        this axis.
    ref_axis : :class:`numpy.ndarray` of shape (3,)
        Reference length-wise axis. One of the reference vectors is
        orthogonal to this axis.
    helix_directions : :class:`numpy.ndarray` of shape (N, 3)
        array of vectors representing the local direction of each
        helix window.

    Returns
    -------
    :class:`numpy.ndarray` of shape (N,)
        Array of screw angles.
    """
    global_axis = np.asarray(global_axis)
    # normal to the plane of `ref_axis` & `global_axis`
    perp = np.cross(ref_axis, global_axis)
    if not np.any(perp):  # zero when ref_axis, global_axis parallel
        # use random orthogonal vector
        new_ref = [[1, 0, 0], [0, 0, 1]]
        while not np.any(perp) and new_ref:
            perp = np.cross(new_ref.pop(), global_axis)

    # normal for angle to plane of perp and global_axis
    ortho = np.cross(-perp, global_axis)

    # project helix_directions onto global to remove contribution
    norm_global_sq = np.dot(global_axis, global_axis)
    mag_g = np.matmul(global_axis, helix_directions.T) / norm_global_sq
    # projection onto global_axis
    proj_g = mag_g.reshape(-1, 1) @ global_axis.reshape(1, -1)
    # projection onto plane w/o global_axis contribution
    proj_plane = helix_directions - proj_g

    # angles from projection to perp
    refs = np.array([perp, ortho])  # (2, 3)
    norms = _, ortho_norm = np.outer(
        mdamath.pnorm(refs), mdamath.pnorm(proj_plane)
    )
    cos = cos_perp, cos_ortho = np.matmul(refs, proj_plane.T) / norms
    to_perp, to_ortho = np.arccos(np.clip(cos, -1, 1))  # (2, n_vec)
    to_ortho[ortho_norm == 0] = 0  # ?
    to_ortho[cos_perp < 0] *= -1
    to_ortho[to_ortho == -np.pi] = np.pi  # leave 180 alone
    return np.rad2deg(to_ortho)


from typing import Dict, Tuple, Union, Optional, List
import numpy as np
from numpy.typing import NDArray

def helix_analysis(
    positions: NDArray[np.float64], 
    ref_axis: Union[Tuple[float, float, float], NDArray[np.float64]] = (0, 0, 1)
) -> Dict[str, NDArray]:
    r"""
    Calculate detailed helix properties from atomic coordinates using a sliding window approach.

    This function implements the HELANAL algorithm to analyze protein helices based on
    their Cα atom coordinates. It calculates various geometric properties including
    local twists, bends, and screw angles that characterize the helix structure.

    Parameters
    ----------
    positions : NDArray[float64], shape (N, 3)
        Atomic coordinates of Cα atoms in the helix. Must contain at least 9 atoms
        for meaningful analysis.
    ref_axis : array-like, shape (3,), optional
        Reference axis used for calculating tilt and screw angles. Default is z-axis (0, 0, 1).

    Returns
    -------
    dict
        Dictionary containing the following keys:
        - 'local_twists' : NDArray[float64], shape (N-3,)
            Local twist angle (in degrees) between consecutive residues (i+1 to i+2).
        - 'local_nres_per_turn' : NDArray[float64], shape (N-3,)
            Number of residues per turn, calculated as 360° / local_twist.
        - 'local_axes' : NDArray[float64], shape (N-3, 3)
            Local helix axis vectors for each window.
        - 'local_bends' : NDArray[float64], shape (N-6,)
            Angles (in degrees) between local helix axes 3 windows apart.
        - 'local_heights' : NDArray[float64], shape (N-3,)
            Rise per residue (in Å) along the local helix axis.
        - 'local_helix_directions' : NDArray[float64], shape (N-2, 3)
            Unit vectors from each local origin to the next Cα atom.
        - 'local_origins' : NDArray[float64], shape (N-2, 3)
            Projected origins for each local helix window.
        - 'all_bends' : NDArray[float64], shape (N-3, N-3)
            Matrix of angles between all pairs of local helix axes.
        - 'global_axis' : NDArray[float64], shape (3,)
            Best-fit axis through all helix origins.
        - 'local_screw_angles' : NDArray[float64], shape (N-2,)
            Azimuthal angles (in degrees) in the plane perpendicular to global_axis.

    Notes
    -----
    The algorithm uses a sliding window of 4 consecutive Cα atoms to calculate
    local helix properties. For a helix with N residues, the function returns
    properties for N-3 windows.

    The reference axis is used to define the zero angle for screw angle calculations.
    By default, this is along the z-axis (0, 0, 1).

    References
    ----------
    .. [1] Bansal M, Kumar S, Velavan R. (2000). HELANAL - A program to
       characterise helix geometry in proteins. J Biomol Struct Dyn, 17(5):811-819.
    .. [2] Sugeta H, Miyazawa T. (1967). General method for calculating helical
       parameters of polymer chains from bond lengths, bond angles and internal
       rotation angles. Biopolymers, 5:673-679.
    
    Examples
    --------
    >>> import numpy as np
    >>> from MDAnalysis.analysis import helix_analysis
    >>> # Generate a simple alpha-helix
    >>> t = np.linspace(0, 4*np.pi, 12)  # 12 residues
    >>> x = np.cos(t)
    >>> y = np.sin(t)
    >>> z = t / (4*np.pi) * 16  # ~1.5Å rise per residue
    >>> positions = np.column_stack((x, y, z))
    >>> # Analyze helix properties
    >>> results = helix_analysis(positions)
    >>> # Get average twist angle
    >>> avg_twist = results['local_twists'].mean()
    >>> print(f"Average twist angle: {avg_twist:.1f}°")
    """
    # Input validation
    positions = np.asarray(positions, dtype=np.float64)
    if positions.ndim != 2 or positions.shape[1] != 3:
        raise ValueError("positions must be a 2D array of shape (N, 3)")
    if len(positions) < 9:
        warnings.warn(
            f"Only {len(positions)} atoms provided. For reliable helix analysis, "
            "at least 9 consecutive Cα atoms are recommended.",
            UserWarning
        )
    
    ref_axis = np.asarray(ref_axis, dtype=np.float64)
    if ref_axis.shape != (3,):
        raise ValueError("ref_axis must be a vector of length 3")

    # Optimized vector calculations using numpy array operations
    # Pre-allocate arrays to avoid repeated memory allocation
    n_res = len(positions)
    
    # Calculate vectors between consecutive Cα atoms
    vectors = np.subtract(positions[1:], positions[:-1], dtype=np.float32)  # (n_res-1, 3)
    
    # Calculate bisectors (approximate, not normalized)
    bisectors = np.subtract(vectors[:-1], vectors[1:], dtype=np.float32)  # (n_res-2, 3)
    
    # Pre-compute magnitudes for normalization
    vector_norms = np.linalg.norm(vectors, axis=1)
    bisector_norms = np.linalg.norm(bisectors, axis=1)
    
    # Avoid division by zero for zero-length vectors
    vector_norms[vector_norms == 0] = np.finfo(np.float32).eps
    bisector_norms[bisector_norms == 0] = np.finfo(np.float32).eps
    
    # Normalize vectors for subsequent calculations
    vectors_normed = vectors / vector_norms[:, np.newaxis]
    bisectors_normed = bisectors / bisector_norms[:, np.newaxis]
    
    # Calculate local twists and residues per turn using optimized vector operations
    cos_theta = np.einsum('ij,ij->i', bisectors_normed[:-1], bisectors_normed[1:], optimize=True)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)  # Ensure valid range for arccos
    local_twists = np.rad2deg(np.arccos(cos_theta))
    
    # Avoid division by zero in nres_per_turn calculation
    twist_mask = np.abs(cos_theta) < 1.0 - 1e-10  # Avoid near-1 values that cause numerical instability
    local_nres_per_turn = np.full_like(cos_theta, np.nan)  # Initialize with NaN
    local_nres_per_turn[twist_mask] = 2 * np.pi / np.arccos(cos_theta[twist_mask])
    
    # Calculate local helix axes using optimized cross product
    local_axes = np.cross(vectors_normed[:-1], bisectors_normed, axis=1)
    local_axes_norm = np.linalg.norm(local_axes, axis=1, keepdims=True)
    local_axes_norm[local_axes_norm < np.finfo(np.float32).eps] = 1.0  # Avoid division by zero
    local_axes = np.divide(local_axes, local_axes_norm, out=local_axes, where=local_axes_norm!=0)
    
    # Optimized calculation of local origins using vectorized operations
    v1 = positions[:-2]  # First atom in window (n-2, 3)
    v2 = positions[1:-1]  # Middle atom (n-2, 3)
    v3 = positions[2:]    # Last atom (n-2, 3)
    
    v13 = v3 - v1  # Vector from first to last atom (n-2, 3)
    v13_sq = np.einsum('ij,ij->i', v13, v13, optimize=True)  # Squared lengths (n-2,)
    
    # Calculate projection parameter t, handling division by zero
    t = np.zeros_like(v13_sq)
    mask = v13_sq > np.finfo(np.float32).eps
    t[mask] = np.einsum('ij,ij->i', v2 - v1, v13, optimize=True)[mask] / v13_sq[mask]
    t = np.clip(t, 0.0, 1.0)  # Clamp t to [0,1]
    
    # Calculate local origins using broadcasting
    local_origins = v1 + t[:, np.newaxis] * v13
    
    # Calculate and normalize local helix directions
    local_helix_directions = v2 - local_origins
    dir_norms = np.linalg.norm(local_helix_directions, axis=1, keepdims=True)
    dir_norms[dir_norms < np.finfo(np.float32).eps] = 1.0  # Avoid division by zero
    local_helix_directions = np.divide(local_helix_directions, dir_norms, 
                                      out=local_helix_directions, 
                                      where=dir_norms!=0)
    
    # Calculate local heights using optimized dot product
    local_heights = np.abs(np.einsum('ij,ij->i', v2 - local_origins, local_axes, optimize=True))
    
    # Calculate local bends between axes 3 windows apart
    if len(local_axes) >= 6:
        cos_bends = np.einsum('ij,ij->i', local_axes[3:], local_axes[:-3], optimize=True)
        cos_bends = np.clip(cos_bends, -1.0, 1.0)  # Ensure valid range
        local_bends = np.rad2deg(np.arccos(cos_bends))
    else:
        local_bends = np.array([], dtype=np.float32)
    
    # Calculate radius of local cylinder
    radii = (np.linalg.norm(bisectors, axis=1)**0.5) / (2 * (1.0 - cos_theta))  # (n_res-3,)
    # special case: angle b/w bisectors is 0 (should virtually never happen)
    # guesstimate radius = half bisector magnitude
    radii = np.where(cos_theta != 1, radii, (np.linalg.norm(bisectors, axis=1)**0.5) / 2)
    
    # get origins by subtracting radius from atom i+1
    origins = positions[1:-1].copy()  # (n_res-2, 3)
    origins[:-1] -= (radii * local_helix_directions[:-1].T).T
    # subtract radius from atom i+2 in last one
    origins[-1] -= radii[-1] * local_helix_directions[-1]

    helix_axes = vector_of_best_fit(origins)
    screw = local_screw_angles(
        helix_axes, np.asarray(ref_axis), local_helix_directions
    )

    results = {
        "local_twists": local_twists,
        "local_nres_per_turn": local_nres_per_turn,
        "local_axes": local_axes,
        "local_bends": local_bends,
        "local_heights": heights,
        "local_helix_directions": local_helix_directions,
        "local_origins": origins,
        "all_bends": bend_matrix,
        "global_axis": helix_axes,
        "local_screw_angles": screw,
    }
    return results


class HELANAL(AnalysisBase):
    r"""
    Perform HELANAL helix analysis on your trajectory.

    Parameters
    ----------
    universe : Universe or AtomGroup
        The Universe or AtomGroup to apply the analysis to.
    select : str or iterable of str, optional
        The selection string to create an atom selection that the HELANAL
        analysis is applied to. Note that HELANAL is designed to work on the
        alpha-carbon atoms of protein residues. If you pass in multiple
        selections, the selections will be analysed separately.
    ref_axis : array-like of length 3, optional
        The reference axis used to calculate the tilt of the vector
        of best fit, and the local screw angles.
    flatten_single_helix : bool, optional
        Whether to flatten results if only one selection is passed.
    split_residue_sequences : bool, optional
        Whether to split the residue sequence into individual helices.
        This keyword only applies if a residue gap is present in the
        AtomGroup generated by a ``select`` string.
        If ``False``, the residues will be analysed as a single helix.
        If ``True``, each group of consecutive residues will be treated
        as a separate helix.
    verbose : bool, optional
        Turn on more logging and debugging.

    Attributes
    ----------
    results.local_twists : array or list of arrays
        The local twist angle from atom i+1 to i+2.
        Each array has shape (n_frames, n_residues-3)
    results.local_nres_per_turn : array or list of arrays
        Number of residues per turn, based on local_twist.
        Each array has shape (n_frames, n_residues-3)
    results.local_axes : array or list of arrays
        The length-wise helix axis of the local window.
        Each array has shape (n_frames, n_residues-3, 3)
    results.local_heights : array or list of arrays
        The rise of each local helix.
        Each array has shape (n_frames, n_residues-3)
    results.local_helix_directions : array or list of arrays
        The unit vector from each local origin to atom i+1.
        Each array has shape (n_frames, n_residues-2, 3)
    results.local_origins :array or list of arrays
        The projected origin for each helix.
        Each array has shape (n_frames, n_residues-2, 3)
    results.local_screw_angles : array or list of arrays
        The local screw angle for each helix.
        Each array has shape (n_frames, n_residues-2)
    results.local_bends : array or list of arrays
        The angles between local helix axes, 3 windows apart.
        Each array has shape (n_frames, n_residues-6)
    results.all_bends : array or list of arrays
        The angles between local helix axes.
        Each array has shape (n_frames, n_residues-3, n_residues-3)
    results.global_axis : array or list of arrays
        The length-wise axis for the overall helix. This points at
        the first helix window in the helix, so it runs opposite to
        the direction of the residue numbers.
        Each array has shape (n_frames, 3)
    results.global_tilts : array or list of arrays
        The angle between the global axis and the reference axis.
        Each array has shape (n_frames,)
    results.summary : dict or list of dicts
        Summary of stats for each property: the mean, the sample
        standard deviation, and the mean absolute deviation.
    """

    # shapes of properties from each frame, relative to n_residues
    attr_shapes = {
        "local_twists": (-3,),
        "local_bends": (-6,),
        "local_heights": (-3,),
        "local_nres_per_turn": (-3,),
        "local_origins": (-2, 3),
        "local_axes": (-3, 3),
        "local_helix_directions": (-2, 3),
        "local_screw_angles": (-2,),
    }

    def __init__(
        self,
        universe,
        select="name CA",
        ref_axis=(0, 0, 1),
        verbose=False,
        flatten_single_helix=True,
        split_residue_sequences=True,
    ):
        super(HELANAL, self).__init__(
            universe.universe.trajectory, verbose=verbose
        )
        selections = util.asiterable(select)
        atomgroups = [universe.select_atoms(s) for s in selections]
        consecutive = []
        # check that residues are consecutive and long enough sequence
        for s, ag in zip(selections, atomgroups):
            groups = util.group_same_or_consecutive_integers(ag.resindices)
            counter = 0
            if len(groups) > 1:
                msg = "Your selection {} has gaps in the residues.".format(s)
                if split_residue_sequences:
                    msg += " Splitting into {} helices.".format(len(groups))
                else:
                    groups = [ag.resindices]
                warnings.warn(msg)

            for g in groups:
                ng = len(g)
                counter += ng
                if ng < 9:
                    warnings.warn(
                        "Fewer than 9 atoms found for helix in "
                        "selection {} with these resindices: {}. "
                        "This sequence will be skipped. HELANAL "
                        "is designed to work on at sequences of "
                        "≥9 residues.".format(s, g)
                    )
                    continue

                ids, counts = np.unique(g, return_counts=True)
                if np.any(counts > 1):
                    dup = ", ".join(map(str, ids[counts > 1]))
                    warnings.warn(
                        "Your selection {} includes multiple atoms "
                        "for residues with these resindices: {}."
                        "HELANAL is designed to work on one alpha-"
                        "carbon per residue.".format(s, dup)
                    )

                consecutive.append(ag[counter - ng : counter])

        self.atomgroups = consecutive
        self.ref_axis = np.asarray(ref_axis)
        self._flatten = flatten_single_helix

    def _zeros_per_frame(self, dims, n_positions=0):
        """Create zero arrays where first 2 dims are n_frames, n_values"""
        first = dims[0] + n_positions
        npdims = (
            self.n_frames,
            first,
        ) + dims[
            1:
        ]  # py27 workaround
        return np.zeros(npdims, dtype=np.float64)

    def _prepare(self):
        n_res = [len(ag) for ag in self.atomgroups]

        for key, dims in self.attr_shapes.items():
            empty = [self._zeros_per_frame(dims, n_positions=n) for n in n_res]
            self.results[key] = empty

        self.results.global_axis = [self._zeros_per_frame((3,)) for n in n_res]
        self.results.all_bends = [
            self._zeros_per_frame((n - 3, n - 3)) for n in n_res
        ]

    def _single_frame(self):
        _f = self._frame_index
        for i, ag in enumerate(self.atomgroups):
            results = helix_analysis(ag.positions, ref_axis=self.ref_axis)
            for key, value in results.items():
                attr = self.results[key]
                attr[i][_f] = value

    def _conclude(self):
        # compute tilt of global axes
        self.results.global_tilts = tilts = []
        norm_ref = (self.ref_axis**2).sum() ** 0.5
        for axes in self.results.global_axis:
            cos = np.matmul(self.ref_axis, axes.T) / (
                mdamath.pnorm(axes) * norm_ref
            )
            cos = np.clip(cos, -1.0, 1.0)
            tilts.append(np.rad2deg(np.arccos(cos)))

        global_attrs = ["global_axis", "global_tilts", "all_bends"]
        attrnames = list(self.attr_shapes.keys()) + global_attrs
        # summarise
        self.results.summary = []
        for i in range(len(self.atomgroups)):
            stats = {}
            for name in attrnames:
                attr = self.results[name]
                mean = attr[i].mean(axis=0)
                dev = np.abs(attr[i] - mean)
                stats[name] = {
                    "mean": mean,
                    "sample_sd": attr[i].std(axis=0, ddof=1),
                    "abs_dev": dev.mean(axis=0),
                }
            self.results.summary.append(stats)

        # flatten?
        if len(self.atomgroups) == 1 and self._flatten:
            for name in attrnames + ["summary"]:
                attr = self.results[name]
                self.results[name] = attr[0]

    def universe_from_origins(self):
        """
        Create MDAnalysis Universe from the local origins.

        Returns
        -------
        Universe or list of Universes
        """
        try:
            origins = self.results.local_origins
        except AttributeError:
            raise ValueError("Call run() before universe_from_origins")

        if not isinstance(origins, list):
            origins = [origins]

        universe = []
        for xyz in origins:
            n_res = xyz.shape[1]
            u = mda.Universe.empty(
                n_res,
                n_residues=n_res,
                atom_resindex=np.arange(n_res),
                trajectory=True,
            ).load_new(xyz)
            universe.append(u)
        if not isinstance(self.results.local_origins, list):
            universe = universe[0]
        return universe
