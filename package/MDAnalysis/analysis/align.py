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

"""Coordinate fitting and alignment --- :mod:`MDAnalysis.analysis.align`
=====================================================================

:Author: Oliver Beckstein, Joshua Adelman
:Year: 2010--2013
:Copyright: GNU Public License v3

The module contains functions to fit a target structure to a reference
structure. They use the fast QCP algorithm to calculate the root mean
square distance (RMSD) between two coordinate sets [Theobald2005]_ and
the rotation matrix *R* that minimizes the RMSD [Liu2010]_. (Please
cite these references when using this module.).

Typically, one selects a group of atoms (such as the C-alphas),
calculates the RMSD and transformation matrix, and applys the
transformation to the current frame of a trajectory to obtain the
rotated structure. The :func:`alignto` and :class:`AlignTraj`
functions can be used to do this for individual frames and
trajectories respectively.

The :ref:`RMS-fitting-tutorial` shows how to do the individual steps
manually and explains the intermediate steps.

See Also
--------
:mod:`MDAnalysis.analysis.rms`
     contains functions to compute RMSD (when structural alignment is not
     required)
:mod:`MDAnalysis.lib.qcprot`
     implements the fast RMSD algorithm.


.. _RMS-fitting-tutorial:

RMS-fitting tutorial
--------------------

The example uses files provided as part of the MDAnalysis test suite
(in the variables :data:`~MDAnalysis.tests.datafiles.PSF`,
:data:`~MDAnalysis.tests.datafiles.DCD`, and
:data:`~MDAnalysis.tests.datafiles.PDB_small`). For all further
examples execute first ::

   >>> import MDAnalysis as mda
   >>> from MDAnalysis.analysis import align
   >>> from MDAnalysis.analysis.rms import rmsd
   >>> from MDAnalysis.tests.datafiles import PSF, DCD, PDB_small


In the simplest case, we can simply calculate the C-alpha RMSD between
two structures, using :func:`rmsd`::

   >>> ref = mda.Universe(PDB_small)
   >>> mobile = mda.Universe(PSF, DCD)
   >>> rmsd(mobile.select_atoms('name CA').positions, ref.select_atoms('name CA').positions)
   28.20178579474479

Note that in this example translations have not been removed. In order
to look at the pure rotation one needs to superimpose the centres of
mass (or geometry) first::

   >>> rmsd(mobile.select_atoms('name CA').positions, ref.select_atoms('name CA').positions, center=True)
   21.892591663632704

This has only done a translational superposition. If you want to also do a
rotational superposition use the superposition keyword. This will calculate a
minimized RMSD between the reference and mobile structure::

   >>> rmsd(mobile.select_atoms('name CA').positions, ref.select_atoms('name CA').positions, 
   ...      superposition=True)
   6.809396586471815

The rotation matrix that superimposes *mobile* on *ref* while
minimizing the CA-RMSD is obtained with the :func:`rotation_matrix`
function ::

   >>> mobile0 = mobile.select_atoms('name CA').positions - mobile.select_atoms('name CA').center_of_mass()
   >>> ref0 = ref.select_atoms('name CA').positions - ref.select_atoms('name CA').center_of_mass()
   >>> R, rmsd = align.rotation_matrix(mobile0, ref0)
   >>> rmsd
   6.809396586471805
   >>> R
   array([[ 0.14514539, -0.27259113,  0.95111876],
   ...    [ 0.88652593,  0.46267112, -0.00268642],
   ...    [-0.43932289,  0.84358136,  0.30881368]])

Putting all this together one can superimpose all of *mobile* onto *ref*::

   >>> mobile.atoms.translate(-mobile.select_atoms('name CA').center_of_mass())
   <AtomGroup with 3341 atoms>
   >>> mobile.atoms.rotate(R)
   <AtomGroup with 3341 atoms>
   >>> mobile.atoms.translate(ref.select_atoms('name CA').center_of_mass())
   <AtomGroup with 3341 atoms>
   >>> mobile.atoms.write("mobile_on_ref.pdb")


Common usage
------------

To **fit a single structure** with :func:`alignto`::

   >>> ref = mda.Universe(PSF, PDB_small)
   >>> mobile = mda.Universe(PSF, DCD)     # we use the first frame
   >>> align.alignto(mobile, ref, select="protein and name CA", weights="mass")
   (21.892591663632704, 6.809396586471809)

This will change *all* coordinates in *mobile* so that the protein
C-alpha atoms are optimally superimposed (translation and rotation).

To **fit a whole trajectory** to a reference structure with the
:class:`AlignTraj` class::

   >>> ref = mda.Universe(PSF, PDB_small)   # reference structure 1AKE
   >>> trj = mda.Universe(PSF, DCD)         # trajectory of change 1AKE->4AKE
   >>> alignment = align.AlignTraj(trj, ref, filename='rmsfit.dcd')
   >>> alignment.run()
   <MDAnalysis.analysis.align.AlignTraj object at ...> 

It is also possible to align two arbitrary structures by providing a
mapping between atoms based on a sequence alignment. This allows
fitting of structural homologs or wild type and mutant.

If a alignment was provided as "sequences.aln" one would first produce
the appropriate MDAnalysis selections with the :func:`fasta2select`
function and then feed the resulting dictionary to :class:`AlignTraj`::

   >>> seldict = align.fasta2select('sequences.aln') # doctest: +SKIP
   >>> alignment = align.AlignTraj(trj, ref, filename='rmsfit.dcd', select=seldict) # doctest: +SKIP
   >>> alignment.run() # doctest: +SKIP

(See the documentation of the functions for this advanced usage.)


Functions and Classes
---------------------

.. versionchanged:: 0.10.0
   Function :func:`~MDAnalysis.analysis.rms.rmsd` was removed from
   this module and is now exclusively accessible as
   :func:`~MDAnalysis.analysis.rms.rmsd`.

.. versionchanged:: 0.16.0
   Function :func:`~MDAnalysis.analysis.align.rms_fit_trj` deprecated
   in favor of :class:`AlignTraj` class.

.. versionchanged:: 0.17.0
   removed deprecated :func:`~MDAnalysis.analysis.align.rms_fit_trj`

.. autofunction:: alignto
.. autoclass:: AlignTraj
.. autoclass:: AverageStructure
.. autofunction:: rotation_matrix
.. autofunction:: iterative_average


Helper functions
----------------

The following functions are used by the other functions in this
module. They are probably of more interest to developers than to
normal users.

.. autofunction:: _fit_to
.. autofunction:: fasta2select
.. autofunction:: sequence_alignment
.. autofunction:: get_matching_atoms

"""
import os.path
import warnings
import logging
import collections

import numpy as np

try:
    import Bio.AlignIO
    import Bio.Align
    import Bio.Align.Applications
except ImportError:
    HAS_BIOPYTHON = False
else:
    HAS_BIOPYTHON = True

import MDAnalysis as mda
import MDAnalysis.lib.qcprot as qcp
from MDAnalysis.exceptions import SelectionError, SelectionWarning
import MDAnalysis.analysis.rms as rms
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.lib.util import get_weights
from MDAnalysis.lib.util import deprecate   # remove 3.0
from MDAnalysis.lib.log import ProgressBar
from ..due import due, Doi

from .base import AnalysisBase

logger = logging.getLogger('MDAnalysis.analysis.align')


def rotation_matrix(a, b, weights=None):
    r"""Returns the 3x3 rotation matrix `R` for RMSD fitting coordinate
    sets `a` and `b`.

    The rotation matrix `R` transforms vector `a` to overlap with
    vector `b` (i.e., `b` is the reference structure):

    .. math::
       \mathbf{b} = \mathsf{R} \cdot \mathbf{a}

    Parameters
    ----------
    a : array_like
        coordinates that are to be rotated ("mobile set"); array of N atoms
        of shape N*3 as generated by, e.g.,
        :attr:`MDAnalysis.core.groups.AtomGroup.positions`.
    b : array_like
        reference coordinates; array of N atoms of shape N*3 as generated by,
        e.g., :attr:`MDAnalysis.core.groups.AtomGroup.positions`.
    weights : array_like (optional)
        array of floats of size N for doing weighted RMSD fitting (e.g. the
        masses of the atoms)

    Returns
    -------
    R : ndarray
        rotation matrix
    rmsd : float
        RMSD between `a` and `b` before rotation
    ``(R, rmsd)`` rmsd and rotation matrix *R*

    Example
    -------
    `R` can be used as an argument for
    :meth:`MDAnalysis.core.groups.AtomGroup.rotate` to generate a rotated
    selection, e.g. ::

        >>> from MDAnalysisTests.datafiles import TPR, TRR
        >>> from MDAnalysis.analysis import align
        >>> A = mda.Universe(TPR,TRR)
        >>> B = A.copy()
        >>> R = rotation_matrix(A.select_atoms('backbone').positions,
        ...                     B.select_atoms('backbone').positions)[0]
        >>> A.atoms.rotate(R)
        <AtomGroup with 47681 atoms>
        >>> A.atoms.write("rotated.pdb")

    Notes
    -----
    The function does *not* shift the centers of mass or geometry;
    this needs to be done by the user.

    See Also
    --------
    MDAnalysis.analysis.rms.rmsd: Calculates the RMSD between *a* and *b*.
    alignto: A complete fit of two structures.
    AlignTraj: Fit a whole trajectory.
    """


    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)

    if a.shape != b.shape:
        raise ValueError("'a' and 'b' must have same shape")

    if np.allclose(a, b) and weights is None:
        return np.eye(3, dtype=np.float64), 0.0

    N = b.shape[0]

    if weights is not None:
        # qcp does NOT divide weights relative to the mean
        weights = np.asarray(weights, dtype=np.float64) / np.mean(weights)

    rot = np.zeros(9, dtype=np.float64)

    # Need to transpose coordinates such that the coordinate array is
    # 3xN instead of Nx3. Also qcp requires that the dtype be float64
    # (I think we swapped the position of ref and traj in CalcRMSDRotationalMatrix
    # so that R acts **to the left** and can be broadcasted; we're saving
    # one transpose. [orbeckst])
    rmsd = qcp.CalcRMSDRotationalMatrix(a, b, N, rot, weights)
    return rot.reshape(3, 3), rmsd


def _fit_to(mobile_coordinates, ref_coordinates, mobile_atoms,
            mobile_com, ref_com, weights=None):
    r"""Perform an rmsd-fitting to determine rotation matrix and align atoms

    Parameters
    ----------
    mobile_coordinates : ndarray
        Coordinates of atoms to be aligned
    ref_coordinates : ndarray
        Coordinates of atoms to be fit against
    mobile_atoms : AtomGroup
        Atoms to be translated
    mobile_com: ndarray
        array of xyz coordinate of mobile center of mass
    ref_com : ndarray
        array of xyz coordinate of reference center of mass
    weights : array_like (optional)
       choose weights. With ``None`` weigh each atom equally. If a float array
       of the same length as `mobile_coordinates` is provided, use each element
       of the `array_like` as a weight for the corresponding atom in
       `mobile_coordinates`.

    Returns
    -------
    mobile_atoms : AtomGroup
        AtomGroup of translated and rotated atoms
    min_rmsd : float
        Minimum rmsd of coordinates

    Notes
    -----
    This function assumes that `mobile_coordinates` and `ref_coordinates` have
    already been shifted so that their centers of geometry (or centers of mass,
    depending on `weights`) coincide at the origin. `mobile_com` and `ref_com`
    are the centers *before* this shift.

    1. The rotation matrix :math:`\mathsf{R}` is determined with
       :func:`rotation_matrix` directly from `mobile_coordinates` and
       `ref_coordinates`.
    2. `mobile_atoms` :math:`X` is rotated according to the
       rotation matrix and the centers according to

       .. math::

           X' = \mathsf{R}(X - \bar{X}) + \bar{X}_{\text{ref}}

       where :math:`\bar{X}` is the center.

    """
    R, min_rmsd = rotation_matrix(mobile_coordinates, ref_coordinates,
                                  weights=weights)

    mobile_atoms.translate(-mobile_com)
    mobile_atoms.rotate(R)
    mobile_atoms.translate(ref_com)

    return mobile_atoms, min_rmsd


def alignto(mobile, reference, select=None, weights=None,
            subselection=None, tol_mass=0.1, strict=False,
            match_atoms=True):
    """Perform a spatial superposition by minimizing the RMSD.

    Spatially align the group of atoms `mobile` to `reference` by
    doing a RMSD fit on `select` atoms.

    The superposition is done in the following way:

    1. A rotation matrix is computed that minimizes the RMSD between
       the coordinates of `mobile.select_atoms(sel1)` and
       `reference.select_atoms(sel2)`; before the rotation, `mobile` is
       translated so that its center of geometry (or center of mass)
       coincides with the one of `reference`. (See below for explanation of
       how *sel1* and *sel2* are derived from `select`.)

    2. All atoms in :class:`~MDAnalysis.core.universe.Universe` that
       contain `mobile` are shifted and rotated. (See below for how
       to change this behavior through the `subselection` keyword.)

    The `mobile` and `reference` atom groups can be constructed so that they
    already match atom by atom. In this case, `select` should be set to "all"
    (or ``None``) so that no further selections are applied to `mobile` and
    `reference`, therefore preserving the exact atom ordering (see
    :ref:`ordered-selections-label`).

    .. Warning:: The atom order for `mobile` and `reference` is *only*
       preserved when `select` is either "all" or ``None``. In any other case,
       a new selection will be made that will sort the resulting AtomGroup by
       index and therefore destroy the correspondence between the two groups.
       **It is safest not to mix ordered AtomGroups with selection strings.**

    Parameters
    ----------
    mobile : Universe or AtomGroup
       structure to be aligned, a
       :class:`~MDAnalysis.core.groups.AtomGroup` or a whole
       :class:`~MDAnalysis.core.universe.Universe`
    reference : Universe or AtomGroup
       reference structure, a :class:`~MDAnalysis.core.groups.AtomGroup`
       or a whole :class:`~MDAnalysis.core.universe.Universe`
    select : str or dict or tuple (optional)
       The selection to operate on; can be one of:

       1. any valid selection string for
          :meth:`~MDAnalysis.core.groups.AtomGroup.select_atoms` that
          produces identical selections in `mobile` and `reference`; or

       2. a dictionary ``{'mobile': sel1, 'reference': sel2}`` where *sel1*
          and *sel2* are valid selection strings that are applied to
          `mobile` and `reference` respectively (the
          :func:`MDAnalysis.analysis.align.fasta2select` function returns such
          a dictionary based on a ClustalW_ or STAMP_ sequence alignment); or

       3. a tuple ``(sel1, sel2)``

       When using 2. or 3. with *sel1* and *sel2* then these selection strings
       are applied to `atomgroup` and `reference` respectively and should
       generate *groups of equivalent atoms*.  *sel1* and *sel2* can each also
       be a *list of selection strings* to generate a
       :class:`~MDAnalysis.core.groups.AtomGroup` with defined atom order as
       described under :ref:`ordered-selections-label`).
    match_atoms : bool (optional)
        Whether to match the mobile and reference atom-by-atom. Default ``True``.
    weights : {"mass", ``None``} or array_like (optional)
       choose weights. With ``"mass"`` uses masses as weights; with ``None``
       weigh each atom equally. If a float array of the same length as
       `mobile` is provided, use each element of the `array_like` as a
       weight for the corresponding atom in `mobile`.
    tol_mass: float (optional)
       Reject match if the atomic masses for matched atoms differ by more than
       `tol_mass`, default [0.1]
    strict: bool (optional)
       ``True``
           Will raise :exc:`SelectionError` if a single atom does not
           match between the two selections.
       ``False`` [default]
           Will try to prepare a matching selection by dropping
           residues with non-matching atoms. See :func:`get_matching_atoms`
           for details.
    subselection : str or AtomGroup or None (optional)
       Apply the transformation only to this selection.

       ``None`` [default]
           Apply to ``mobile.universe.atoms`` (i.e., all atoms in the
           context of the selection from `mobile` such as the rest of a
           protein, ligands and the surrounding water)
       *selection-string*
           Apply to ``mobile.select_atoms(selection-string)``
       :class:`~MDAnalysis.core.groups.AtomGroup`
           Apply to the arbitrary group of atoms

    Returns
    -------
    old_rmsd : float
        RMSD before spatial alignment
    new_rmsd : float
        RMSD after spatial alignment

    See Also
    --------
    AlignTraj: More efficient method for RMSD-fitting trajectories.


    .. _ClustalW: http://www.clustal.org/
    .. _STAMP: http://www.compbio.dundee.ac.uk/manuals/stamp.4.2/

    .. versionchanged:: 1.0.0
       Added *match_atoms* keyword to toggle atom matching.

    .. versionchanged:: 0.8
       Added check that the two groups describe the same atoms including
       the new *tol_mass* keyword.

    .. versionchanged:: 0.10.0
       Uses :func:`get_matching_atoms` to work with incomplete selections
       and new `strict` keyword. The new default is to be lenient whereas
       the old behavior was the equivalent of ``strict = True``.

    .. versionchanged:: 0.16.0
       new general 'weights' kwarg replace `mass_weighted`, deprecated `mass_weighted`
    .. deprecated:: 0.16.0
       Instead of ``mass_weighted=True`` use new ``weights='mass'``

    .. versionchanged:: 0.17.0
       Deprecated keyword `mass_weighted` was removed.
    """
    if select in ('all', None):
        # keep the EXACT order in the input AtomGroups; select_atoms('all')
        # orders them by index, which can lead to wrong results if the user
        # has crafted mobile and reference to match atom by atom
        mobile_atoms = mobile.atoms
        ref_atoms = reference.atoms
    else:
        select = rms.process_selection(select)
        mobile_atoms = mobile.select_atoms(*select['mobile'])
        ref_atoms = reference.select_atoms(*select['reference'])


    ref_atoms, mobile_atoms = get_matching_atoms(ref_atoms, mobile_atoms,
                                                tol_mass=tol_mass,
                                                strict=strict,
                                                match_atoms=match_atoms)

    weights = get_weights(ref_atoms, weights)

    mobile_com = mobile_atoms.center(weights)
    ref_com = ref_atoms.center(weights)

    ref_coordinates = ref_atoms.positions - ref_com
    mobile_coordinates = mobile_atoms.positions - mobile_com

    old_rmsd = rms.rmsd(mobile_coordinates, ref_coordinates, weights)

    if subselection is None:
        # mobile_atoms is Universe
        mobile_atoms = mobile.universe.atoms
    elif isinstance(subselection, str):
        # select mobile_atoms from string
        mobile_atoms = mobile.select_atoms(subselection)
    else:
        try:
            # treat subselection as AtomGroup
            mobile_atoms = subselection.atoms
        except AttributeError:
            err = ("subselection must be a selection string, an"
                   " AtomGroup or Universe or None")
            raise TypeError(err) from None


    # _fit_to DOES subtract center of mass, will provide proper min_rmsd
    mobile_atoms, new_rmsd = _fit_to(mobile_coordinates, ref_coordinates,
                                     mobile_atoms, mobile_com, ref_com,
                                     weights=weights)
    return old_rmsd, new_rmsd


@due.dcite(
        Doi("10.1021/acs.jpcb.7b11988"),
        description="Iterative Calculation of Opimal Reference",
        path="MDAnalysis.analysis.align.iterative_average"
)
def iterative_average(
    mobile, reference=None, select='all', weights=None, niter=100,
    eps=1e-6, verbose=False, **kwargs
):
    """Iteratively calculate an optimal reference that is also the average
    structure after an RMSD alignment.

    The optimal reference is defined as average
    structure of a trajectory, with the optimal reference used as input.
    This function computes the optimal reference by using a starting
    reference for the average structure, which is used as the reference
    to calculate the average structure again. This is repeated until the
    reference structure has converged. :footcite:p:`Linke2018`

    Parameters
    ----------
    mobile : mda.Universe
        Universe containing trajectory to be fitted to reference.
    reference : mda.Universe (optional)
        Universe containing the initial reference structure.
    select : str or tuple or dict (optional)
        Atom selection for fitting a substructue. Default is set to all.
        Can be tuple or dict to define different selection strings for
        mobile and target.
    weights : str, array_like (optional)
        Weights that can be used. If `None` use equal weights, if `'mass'`
        use masses of ref as weights or give an array of arbitrary weights.
    niter : int (optional)
        Maximum number of iterations.
    eps : float (optional)
        RMSD distance at which reference and average are assumed to be
        equal.
    verbose : bool (optional)
        Verbosity.
    **kwargs : dict (optional)
        AverageStructure kwargs.

    Returns
    -------
    avg_struc : AverageStructure
        AverageStructure result from the last iteration.

    Example
    -------
    `iterative_average` can be used to obtain a :class:`MDAnalysis.Universe`
    with the optimal reference structure.

    ::

        import MDAnalysis as mda
        from MDAnalysis.analysis import align
        from MDAnalysisTests.datafiles import PSF, DCD

        u = mda.Universe(PSF, DCD)
        av = align.iterative_average(u, u, verbose=True)

        averaged_universe = av.results.universe

    References
    ----------

    .. footbibliography::

    .. versionadded:: 2.8.0
    """
    if not reference:
        reference = mobile

    select = rms.process_selection(select)
    ref = mda.Merge(reference.select_atoms(*select['reference']))
    sel_mobile = select['mobile'][0]

    weights = get_weights(ref.atoms, weights)

    drmsd = np.inf
    for i in ProgressBar(range(niter)):
        # found a converged structure
        if drmsd < eps:
            break

        avg_struc = AverageStructure(
            mobile, reference=ref, select={
                'mobile': sel_mobile, 'reference': 'all'
                },
            weights=weights, **kwargs
        ).run()
        drmsd = rms.rmsd(ref.atoms.positions, avg_struc.results.positions,
                         weights=weights)
        ref = avg_struc.results.universe

        if verbose:
            logger.debug(
                f"iterative_average(): i = {i}, "
                f"rmsd-change = {drmsd:.5f}, "
                f"ave-rmsd = {avg_struc.results.rmsd:.5f}"
            )

    else:
        errmsg = (
            "iterative_average(): Did not converge in "
            f"{niter} iterations to DRMSD < {eps}. "
            f"Final average RMSD = {avg_struc.results.rmsd:.5f}"
        )
        logger.error(errmsg)
        raise RuntimeError(errmsg)

    logger.info(
        f"iterative_average(): Converged to DRMSD < {eps}. "
        f"Final average RMSD = {avg_struc.results.rmsd:.5f}"
    )

    return avg_struc


class AlignTraj(AnalysisBase):
    """RMS-align trajectory to a reference structure using a selection.

    Both the reference `reference` and the trajectory `mobile` must be
    :class:`MDAnalysis.Universe` instances. If they contain a trajectory then
    it is used. The output file format is determined by the file extension of
    `filename`. One can also use the same universe if one wants to fit to the
    current frame.

    .. versionchanged:: 1.0.0
       ``save()`` has now been removed, as an alternative use ``np.savetxt()``
       on :attr:`results.rmsd`.

    """

    def __init__(self, mobile, reference, select='all', filename=None,
                 prefix='rmsfit_', weights=None,
                 tol_mass=0.1, match_atoms=True, strict=False, force=True, in_memory=False,
                 **kwargs):
        """Parameters
        ----------
        mobile : Universe
            Universe containing trajectory to be fitted to reference
        reference : Universe
            Universe containing trajectory frame to be used as reference
        select : str (optional)
            Set as default to all, is used for Universe.select_atoms to choose
            subdomain to be fitted against
        filename : str (optional)
            Provide a filename for results to be written to
        prefix : str (optional)
            Provide a string to prepend to filename for results to be written
            to
        weights : {"mass", ``None``} or array_like (optional)
            choose weights. With ``"mass"`` uses masses of `reference` as
            weights; with ``None`` weigh each atom equally. If a float array of
            the same length as the selection is provided, use each element of
            the `array_like` as a weight for the corresponding atom in the
            selection.
        tol_mass : float (optional)
            Tolerance given to `get_matching_atoms` to find appropriate atoms
        match_atoms : bool (optional)
            Whether to match the mobile and reference atom-by-atom. Default ``True``.
        strict : bool (optional)
            Force `get_matching_atoms` to fail if atoms can't be found using
            exact methods
        force : bool (optional)
            Force overwrite of filename for rmsd-fitting
        in_memory : bool (optional)
            *Permanently* switch `mobile` to an in-memory trajectory
            so that alignment can be done in-place, which can improve
            performance substantially in some cases. In this case, no file
            is written out (`filename` and `prefix` are ignored) and only
            the coordinates of `mobile` are *changed in memory*.
        verbose : bool (optional)
             Set logger to show more information and show detailed progress of
             the calculation if set to ``True``; the default is ``False``.


        Attributes
        ----------
        reference_atoms : AtomGroup
            Atoms of the reference structure to be aligned against
        mobile_atoms : AtomGroup
            Atoms inside each trajectory frame to be rmsd_aligned
        results.rmsd : :class:`numpy.ndarray`
            Array of the rmsd values of the least rmsd between the mobile_atoms
            and reference_atoms after superposition and minimimization of rmsd

            .. versionadded:: 2.0.0

        rmsd : :class:`numpy.ndarray`
            Alias to the :attr:`results.rmsd` attribute.

            .. deprecated:: 2.0.0
               Will be removed in MDAnalysis 3.0.0. Please use
               :attr:`results.rmsd` instead.

        filename : str
            String reflecting the filename of the file where the aligned
            positions will be written to upon running RMSD alignment


        Notes
        -----
        - If set to ``verbose=False``, it is recommended to wrap the statement
          in a ``try ...  finally`` to guarantee restoring of the log level in
          the case of an exception.
        - The ``in_memory`` option changes the `mobile` universe to an
          in-memory representation (see :mod:`MDAnalysis.coordinates.memory`)
          for the remainder of the Python session. If ``mobile.trajectory`` is
          already a :class:`MemoryReader` then it is *always* treated as if
          ``in_memory`` had been set to ``True``.

        .. versionchanged:: 1.0.0
           Default ``filename`` has now been changed to the current directory.

        .. deprecated:: 0.19.1
           Default ``filename`` directory will change in 1.0 to the current directory.

        .. versionchanged:: 0.16.0
           new general ``weights`` kwarg replace ``mass_weights``

        .. deprecated:: 0.16.0
           Instead of ``mass_weighted=True`` use new ``weights='mass'``

        .. versionchanged:: 0.17.0
           removed deprecated `mass_weighted` keyword

        .. versionchanged:: 1.0.0
           Support for the ``start``, ``stop``, and ``step`` keywords has been
           removed. These should instead be passed to :meth:`AlignTraj.run`.

        .. versionchanged:: 2.0.0
           :attr:`rmsd` results are now stored in a
           :class:`MDAnalysis.analysis.base.Results` instance.

        """
        select = rms.process_selection(select)
        self.ref_atoms = reference.select_atoms(*select['reference'])
        self.mobile_atoms = mobile.select_atoms(*select['mobile'])
        if in_memory or isinstance(mobile.trajectory, MemoryReader):
            mobile.transfer_to_memory()
            filename = None
            logger.info("Moved mobile trajectory to in-memory representation")
        else:
            if filename is None:
                fn = os.path.split(mobile.trajectory.filename)[1]
                filename = prefix + fn
                logger.info('filename of rms_align with no filename given'
                            ': {0}'.format(filename))

            if os.path.exists(filename) and not force:
                raise IOError(
                    'Filename already exists in path and force is not set'
                    ' to True')

        # do this after setting the memory reader to have a reference to the
        # right reader.
        super(AlignTraj, self).__init__(mobile.trajectory, **kwargs)
        if not self._verbose:
            logging.disable(logging.WARN)

        # store reference to mobile atoms
        self.mobile = mobile.atoms

        self.filename = filename

        natoms = self.mobile.n_atoms
        self.ref_atoms, self.mobile_atoms = get_matching_atoms(
            self.ref_atoms, self.mobile_atoms, tol_mass=tol_mass,
            strict=strict, match_atoms=match_atoms)

        # with self.filename == None (in_memory), the NullWriter is chosen
        # (which just ignores input) and so only the in_memory trajectory is
        # retained
        self._writer = mda.Writer(self.filename, natoms, **kwargs)

        self._weights = get_weights(self.ref_atoms, weights)

        logger.info("RMS-fitting on {0:d} atoms.".format(len(self.ref_atoms)))

    def _prepare(self):
        # reference centre of mass system
        self._ref_com = self.ref_atoms.center(self._weights)
        self._ref_coordinates = self.ref_atoms.positions - self._ref_com
        # allocate the array for selection atom coords
        self.results.rmsd = np.zeros((self.n_frames,))

    def _single_frame(self):
        index = self._frame_index
        mobile_com = self.mobile_atoms.center(self._weights)
        mobile_coordinates = self.mobile_atoms.positions - mobile_com
        mobile_atoms, self.results.rmsd[index] = _fit_to(
                mobile_coordinates,
                self._ref_coordinates,
                self.mobile,
                mobile_com,
                self._ref_com, self._weights)
        # write whole aligned input trajectory system
        self._writer.write(mobile_atoms)

    def _conclude(self):
        self._writer.close()
        if not self._verbose:
            logging.disable(logging.NOTSET)

    @property
    def rmsd(self):
        wmsg = ("The `rmsd` attribute was deprecated in MDAnalysis 2.0.0 and "
                "will be removed in MDAnalysis 3.0.0. Please use "
                "`results.rmsd` instead.")
        warnings.warn(wmsg, DeprecationWarning)
        return self.results.rmsd


class AverageStructure(AnalysisBase):
    """RMS-align trajectory to a reference structure using a selection,
    and calculate the average coordinates of the trajectory.

    Both the reference `reference` and the trajectory `mobile` must be
    :class:`MDAnalysis.Universe` instances. If they contain a trajectory, then
    it is used. You can also use the same universe if you want to fit to the
    current frame.

    The output file format is determined by the file extension of
    `filename`.

    Example
    -------

    ::

        import MDAnalysis as mda
        from MDAnalysis.tests.datafiles import PSF, DCD
        from MDAnalysis.analysis import align

        u = mda.Universe(PSF, DCD)

        # align to the third frame and average structure
        av = align.AverageStructure(u, ref_frame=3).run()
        averaged_universe = av.results.universe

    """

    def __init__(self, mobile, reference=None, select='all', filename=None,
                weights=None,
                 tol_mass=0.1, match_atoms=True, strict=False, force=True, in_memory=False,
                 ref_frame=0, **kwargs):
        """Parameters
        ----------
        mobile : Universe
            Universe containing trajectory to be fitted to reference
        reference : Universe (optional)
            Universe containing trajectory frame to be used as reference
        select : str (optional)
            Set as default to all, is used for Universe.select_atoms to choose
            subdomain to be fitted against
        filename : str (optional)
            Provide a filename for results to be written to
        weights : {"mass", ``None``} or array_like (optional)
            choose weights. With ``"mass"`` uses masses of `reference` as
            weights; with ``None`` weigh each atom equally. If a float array of
            the same length as the selection is provided, use each element of
            the `array_like` as a weight for the corresponding atom in the
            selection.
        tol_mass : float (optional)
            Tolerance given to `get_matching_atoms` to find appropriate atoms
        match_atoms : bool (optional)
            Whether to match the mobile and reference atom-by-atom. Default ``True``.
        strict : bool (optional)
            Force `get_matching_atoms` to fail if atoms can't be found using
            exact methods
        force : bool (optional)
            Force overwrite of filename for rmsd-fitting
        in_memory : bool (optional)
            *Permanently* switch `mobile` to an in-memory trajectory
            so that alignment can be done in-place, which can improve
            performance substantially in some cases. In this case, no file
            is written out (`filename` and `prefix` are ignored) and only
            the coordinates of `mobile` are *changed in memory*.
        ref_frame : int (optional)
            frame index to select frame from `reference`
        verbose : bool (optional)
            Set logger to show more information and show detailed progress of
            the calculation if set to ``True``; the default is ``False``.


        Attributes
        ----------
        reference_atoms : AtomGroup
            Atoms of the reference structure to be aligned against
        mobile_atoms : AtomGroup
            Atoms inside each trajectory frame to be rmsd_aligned
        results.universe : :class:`MDAnalysis.Universe`
            New Universe with average positions

            .. versionadded:: 2.0.0

        universe : :class:`MDAnalysis.Universe`
            Alias to the :attr:`results.universe` attribute.

            .. deprecated:: 2.0.0
               Will be removed in MDAnalysis 3.0.0. Please use
               :attr:`results.universe` instead.

        results.positions : np.ndarray(dtype=float)
            Average positions

            .. versionadded:: 2.0.0

        positions : np.ndarray(dtype=float)
            Alias to the :attr:`results.positions` attribute.

            .. deprecated:: 2.0.0
               Will be removed in MDAnalysis 3.0.0. Please use
               :attr:`results.positions` instead.

        results.rmsd : float
            Average RMSD per frame

            .. versionadded:: 2.0.0

        rmsd : float
            Alias to the :attr:`results.rmsd` attribute.

            .. deprecated:: 2.0.0
               Will be removed in MDAnalysis 3.0.0. Please use
               :attr:`results.rmsd` instead.

        filename : str
            String reflecting the filename of the file where the average
            structure is written


        Notes
        -----
        - If set to ``verbose=False``, it is recommended to wrap the statement
          in a ``try ...  finally`` to guarantee restoring of the log level in
          the case of an exception.
        - The ``in_memory`` option changes the `mobile` universe to an
          in-memory representation (see :mod:`MDAnalysis.coordinates.memory`)
          for the remainder of the Python session. If ``mobile.trajectory`` is
          already a :class:`MemoryReader` then it is *always* treated as if
          ``in_memory`` had been set to ``True``.


        .. versionadded:: 1.0.0
        .. versionchanged:: 2.0.0
           :attr:`universe`, :attr:`positions`, and :attr:`rmsd` are now
           stored in a :class:`MDAnalysis.analysis.base.Results` instance.
        """
        if in_memory or isinstance(mobile.trajectory, MemoryReader):
            mobile.transfer_to_memory()
            filename = None
            logger.info("Moved mobile trajectory to in-memory representation")

        # do this after setting the memory reader to have a reference to the
        # right reader.
        super(AverageStructure, self).__init__(mobile.trajectory, **kwargs)
        if not self._verbose:
            logging.disable(logging.WARN)

        self.reference = reference if reference is not None else mobile

        select = rms.process_selection(select)
        self.ref_atoms = self.reference.select_atoms(*select['reference'])
        self.mobile_atoms = mobile.select_atoms(*select['mobile'])

        if len(self.ref_atoms) != len(self.mobile_atoms):
            err = ("Reference and trajectory atom selections do "
                   "not contain the same number of atoms: "
                   "N_ref={0:d}, N_traj={1:d}".format(self.ref_atoms.n_atoms,
                                                      self.mobile_atoms.n_atoms))
            logger.exception(err)
            raise SelectionError(err)
        logger.info("RMS calculation "
                    "for {0:d} atoms.".format(len(self.ref_atoms)))

        # store reference to mobile atoms
        self.mobile = mobile.atoms
        self.ref_frame = ref_frame

        self.filename = filename
        self.results.universe = mda.Merge(self.mobile_atoms)

        natoms = len(self.results.universe.atoms)
        self.ref_atoms, self.mobile_atoms = get_matching_atoms(
            self.ref_atoms, self.mobile_atoms, tol_mass=tol_mass,
            strict=strict, match_atoms=match_atoms)

        # with self.filename == None (in_memory), the NullWriter is chosen
        # (which just ignores input) and so only the in_memory trajectory is
        # retained
        self._writer = mda.Writer(self.filename, natoms)
        self._weights = get_weights(self.ref_atoms, weights)

        logger.info("RMS-fitting on {0:d} atoms.".format(len(self.ref_atoms)))

    def _prepare(self):
        current_frame = self.reference.universe.trajectory.ts.frame
        try:
            # Move to the ref_frame
            # (coordinates MUST be stored in case the ref traj is advanced
            # elsewhere or if ref == mobile universe)
            self.reference.universe.trajectory[self.ref_frame]
            self._ref_com = self.ref_atoms.center(self._weights)
            # makes a copy
            self._ref_coordinates = self.ref_atoms.positions - self._ref_com
            self._ref_positions = self.ref_atoms.positions.copy()
        finally:
            # Move back to the original frame
            self.reference.universe.trajectory[current_frame]

        # allocate the array for selection atom coords
        self.results.positions = np.zeros((len(self.mobile_atoms), 3))
        self.results.rmsd = 0

    def _single_frame(self):
        mobile_com = self.mobile_atoms.center(self._weights)
        mobile_coordinates = self.mobile_atoms.positions - mobile_com
        self.results.rmsd += _fit_to(mobile_coordinates,
                                     self._ref_coordinates,
                                     self.mobile,
                                     mobile_com,
                                     self._ref_com, self._weights)[1]
        self.results.positions += self.mobile_atoms.positions

    def _conclude(self):
        self.results.positions /= self.n_frames
        self.results.rmsd /= self.n_frames
        self.results.universe.load_new(
                self.results.positions.reshape((1, -1, 3)))
        self._writer.write(self.results.universe.atoms)
        self._writer.close()
        if not self._verbose:
            logging.disable(logging.NOTSET)

    @property
    def universe(self):
        wmsg = ("The `universe` attribute was deprecated in MDAnalysis 2.0.0 "
                "and will be removed in MDAnalysis 3.0.0. Please use "
                "`results.universe` instead.")
        warnings.warn(wmsg, DeprecationWarning)
        return self.results.universe

    @property
    def positions(self):
        wmsg = ("The `positions` attribute was deprecated in MDAnalysis 2.0.0 "
                "and will be removed in MDAnalysis 3.0.0. Please use "
                "`results.positions` instead.")
        warnings.warn(wmsg, DeprecationWarning)
        return self.results.positions

    @property
    def rmsd(self):
        wmsg = ("The `rmsd` attribute was deprecated in MDAnalysis 2.0.0 "
                "and will be removed in MDAnalysis 3.0.0. Please use "
                "`results.rmsd` instead.")
        warnings.warn(wmsg, DeprecationWarning)
        return self.results.rmsd


@deprecate(release="2.4.0", remove="3.0",
           message="See the documentation under Notes on how to directly use"
                   "Bio.Align.PairwiseAligner with ResidueGroups.")
def sequence_alignment(mobile, reference, match_score=2, mismatch_penalty=-1,
                       gap_penalty=-2, gapextension_penalty=-0.1):
    """Generate a global sequence alignment between two residue groups.

    The residues in `reference` and `mobile` will be globally aligned.
    The global alignment uses the Needleman-Wunsch algorithm as
    implemented in :mod:`Bio.Align.PairwiseAligner`. The parameters of the dynamic
    programming algorithm can be tuned with the keywords. The defaults
    should be suitable for two similar sequences. For sequences with
    low sequence identity, more specialized tools such as clustalw,
    muscle, tcoffee, or similar should be used.

    Parameters
    ----------
    mobile : AtomGroup
        Atom group to be aligned
    reference : AtomGroup
        Atom group to be aligned against
    match_score : float (optional), default 2
         score for matching residues, default 2
    mismatch_penalty : float (optional), default -1
         penalty for residues that do not match , default : -1
    gap_penalty : float (optional), default -2
         penalty for opening a gap; the high default value creates compact
         alignments for highly identical sequences but might not be suitable
         for sequences with low identity, default : -2
    gapextension_penalty : float (optional), default -0.1
         penalty for extending a gap, default: -0.1

    Returns
    -------
    alignment : tuple
        Tuple of top sequence matching output `('Sequence A', 'Sequence B', score,
        begin, end)`

    Raises
    ------
    ImportError
      If optional dependency Biopython is not available.

    Notes
    -----
    If you prefer to work directly with :mod:`Bio.Align` objects then you can
    run your alignment with :class:`Bio.Alig.PairwiseAligner` as ::

      import Bio.Align.PairwiseAligner

      aligner = Bio.Align.PairwiseAligner(
         mode="global",
         match_score=match_score,
         mismatch_score=mismatch_penalty,
         open_gap_score=gap_penalty,
         extend_gap_score=gapextension_penalty)
      aln = aligner.align(reference.residues.sequence(format="Seq"),
                          mobile.residues.sequence(format="Seq"))

      # choose top alignment with highest score
      topalignment = aln[0]

    The ``topalignment`` is a :class:`Bio.Align.PairwiseAlignment` instance
    that can be used in your bioinformatics workflows.

    See Also
    --------
    BioPython documentation for `PairwiseAligner`_. Alternatively, use
    :func:`fasta2select` with :program:`clustalw2` and the option
    ``is_aligned=False``.


    .. _`PairwiseAligner`:
       https://biopython.org/docs/latest/api/Bio.Align.html#Bio.Align.PairwiseAligner


    .. versionadded:: 0.10.0

    .. versionchanged:: 2.4.0
       Replace use of deprecated :func:`Bio.pairwise2.align.globalms` with
       :class:`Bio.Align.PairwiseAligner`.

    .. versionchanged:: 2.7.0
       Biopython is now an optional dependency which this method requires.

    """
    if not HAS_BIOPYTHON:
        errmsg = ("The `sequence_alignment` method requires an installation "
                  "of `Biopython`. Please install `Biopython` to use this "
                  "method: https://biopython.org/wiki/Download")
        raise ImportError(errmsg)

    aligner = Bio.Align.PairwiseAligner(
        mode="global",
        match_score=match_score,
        mismatch_score=mismatch_penalty,
        open_gap_score=gap_penalty,
        extend_gap_score=gapextension_penalty)
    aln = aligner.align(reference.residues.sequence(format="Seq"),
                        mobile.residues.sequence(format="Seq"))
    # choose top alignment with highest score
    topalignment = aln[0]

    # reconstruct the results tuple that used to be of type Bio.pairwise2.Alignment
    AlignmentTuple = collections.namedtuple(
        "Alignment",
        ["seqA", "seqB", "score", "start", "end"])
    # start/stop are not particularly meaningful and there's no obvious way to
    # get the old pairwise2 start/stop from the new PairwiseAligner output.
    return AlignmentTuple(topalignment[0], topalignment[1],
                          topalignment.score,
                          0, max(reference.n_residues, mobile.n_residues))



def fasta2select(fastafilename, is_aligned=False,
                 ref_resids=None, target_resids=None,
                 ref_offset=0, target_offset=0, verbosity=3,
                 alnfilename=None, treefilename=None, clustalw="clustalw2"):
    """Return selection strings that will select equivalent residues.

    The function aligns two sequences provided in a FASTA file and
    constructs MDAnalysis selection strings of the common atoms. When
    these two strings are applied to the two different proteins they
    will generate AtomGroups of the aligned residues.

    `fastafilename` contains the two un-aligned sequences in FASTA
    format. The reference is assumed to be the first sequence, the
    target the second. ClustalW_ produces a pairwise
    alignment (which is written to a file with suffix ``.aln``).  The
    output contains atom selection strings that select the same atoms
    in the two structures.

    Unless `ref_offset` and/or `target_offset` are specified, the resids
    in the structure are assumed to correspond to the positions in the
    un-aligned sequence, namely the first residue has resid == 1.

    In more complicated cases (e.g., when the resid numbering in the
    input structure has gaps due to missing parts), simply provide the
    sequence of resids as they appear in the topology in `ref_resids` or
    `target_resids`, e.g. ::

       target_resids = [a.resid for a in trj.select_atoms('name CA')]

    (This translation table *is* combined with any value for
    `ref_offset` or `target_offset`!)

    Parameters
    ----------
    fastafilename : str, path to filename
        FASTA file with first sequence as reference and
        second the one to be aligned (ORDER IS IMPORTANT!)
    is_aligned : bool (optional)
        ``False`` (default)
            run clustalw for sequence alignment;
        ``True``
            use the alignment in the file (e.g. from STAMP) [``False``]
    ref_offset : int (optional)
        add this number to the column number in the FASTA file
        to get the original residue number, default: 0
    target_offset : int (optional)
        add this number to the column number in the FASTA file
        to get the original residue number, default: 0
    ref_resids : str (optional)
        sequence of resids as they appear in the reference structure
    target_resids : str (optional)
        sequence of resids as they appear in the target
    alnfilename : str (optional)
        filename of ClustalW alignment (clustal format) that is
        produced by *clustalw* when *is_aligned* = ``False``.
        default ``None`` uses the name and path of *fastafilename* and
        substitutes the suffix with '.aln'.
    treefilename: str (optional)
        filename of ClustalW guide tree (Newick format);
        if default ``None``  the the filename is generated from *alnfilename*
        with the suffix '.dnd' instead of '.aln'
    clustalw : str (optional)
        path to the ClustalW (or ClustalW2) binary; only
        needed for `is_aligned` = ``False``, default: "ClustalW2"

    Returns
    -------
    select_dict : dict
        dictionary with 'reference' and 'mobile' selection string
        that can be used immediately in :class:`AlignTraj` as
        ``select=select_dict``.


    See Also
    --------
    :func:`sequence_alignment`, which does not require external
    programs.

    
    Raises
    ------
    ImportError
      If optional dependency Biopython is not available.


    .. _ClustalW: http://www.clustal.org/
    .. _STAMP: http://www.compbio.dundee.ac.uk/manuals/stamp.4.2/

    .. versionchanged:: 1.0.0
       Passing `alnfilename` or `treefilename` as `None` will create a file in
       the current working directory.
    .. versionchanged:: 2.7.0
       Biopython is now an optional dependency which this method requires.

    """
    if not HAS_BIOPYTHON:
        errmsg = ("The `fasta2select` method requires an installation "
                  "of `Biopython`. Please install `Biopython` to use this "
                  "method: https://biopython.org/wiki/Download")
        raise ImportError(errmsg)

    if is_aligned:
        logger.info("Using provided alignment {}".format(fastafilename))
        with open(fastafilename) as fasta:
            alignment = Bio.AlignIO.read(
                fasta, "fasta")
    else:
        if alnfilename is None:
            filepath, ext = os.path.splitext(fastafilename)
            alnfilename = os.path.basename(filepath) + '.aln'
        if treefilename is None:
            filepath, ext = os.path.splitext(alnfilename)
            treefilename = os.path.basename(filepath) + '.dnd'
        run_clustalw = Bio.Align.Applications.ClustalwCommandline(
            clustalw,
            infile=fastafilename,
            type="protein",
            align=True,
            outfile=alnfilename,
            newtree=treefilename)
        logger.debug(
            "Aligning sequences in %(fastafilename)r with %(clustalw)r.",
            vars())
        logger.debug("ClustalW commandline: %r", str(run_clustalw))
        try:
            stdout, stderr = run_clustalw()
        except:
            logger.exception("ClustalW %(clustalw)r failed", vars())
            logger.info(
                "(You can get clustalw2 from http://www.clustal.org/clustal2/)")
            raise
        with open(alnfilename) as aln:
            alignment = Bio.AlignIO.read(
                aln, "clustal")
        logger.info(
            "Using clustalw sequence alignment {0!r}".format(alnfilename))
        logger.info(
            "ClustalW Newick guide tree was also produced: {0!r}".format(treefilename))

    nseq = len(alignment)
    if nseq != 2:
        raise ValueError(
            "Only two sequences in the alignment can be processed.")

    # implict assertion that we only have two sequences in the alignment
    orig_resids = [ref_resids, target_resids]
    offsets = [ref_offset, target_offset]
    GAP = "-"
    for iseq, a in enumerate(alignment):
        # need iseq index to change orig_resids
        if orig_resids[iseq] is None:
            # build default: assume consecutive numbering of all
            # residues in the alignment
            length = len(a.seq) - a.seq.count(GAP)
            orig_resids[iseq] = np.arange(1, length + 1)
        else:
            orig_resids[iseq] = np.asarray(orig_resids[iseq])
    # add offsets to the sequence <--> resid translation table
    seq2resids = [resids + offset for resids, offset in zip(
        orig_resids, offsets)]
    del orig_resids
    del offsets

    def resid_factory(alignment, seq2resids):
        """Return a function that gives the resid for a position ipos in
        the nseq'th alignment.

        resid = resid_factory(alignment,seq2resids)
        r = resid(nseq,ipos)

        It is based on a look up table that translates position in the
        alignment to the residue number in the original
        sequence/structure.

        The first index of resid() is the alignmment number, the
        second the position in the alignment.

        seq2resids translates the residues in the sequence to resid
        numbers in the psf. In the simplest case this is a linear map
        but if whole parts such as loops are ommitted from the protein
        the seq2resids may have big gaps.

        Format: a tuple of two numpy arrays; the first array is for
        the reference, the second for the target, The index in each
        array gives the consecutive number of the amino acid in the
        sequence, the value the resid in the structure/psf.

        Note: assumes that alignments have same length and are padded if
        necessary.
        """
        # could maybe use Bio.PDB.StructureAlignment instead?
        nseq = len(alignment)
        t = np.zeros((nseq, alignment.get_alignment_length()), dtype=int)
        for iseq, a in enumerate(alignment):
            GAP = "-"
            t[iseq, :] = seq2resids[iseq][np.cumsum(np.where(
                np.array(list(a.seq)) == GAP, 0, 1)) - 1]
            # -1 because seq2resid is index-1 based (resids start at 1)

        def resid(nseq, ipos, t=t):
            return t[nseq, ipos]

        return resid

    resid = resid_factory(alignment, seq2resids)

    res_list = []  # collect individual selection string
    # could collect just resid and type (with/without CB) and
    # then post-process and use ranges for continuous stretches, eg
    # ( resid 1:35 and ( backbone or name CB ) ) or ( resid 36 and backbone )
    for ipos in range(alignment.get_alignment_length()):
        aligned = list(alignment[:, ipos])
        if GAP in aligned:
            continue  # skip residue
        template = "resid %i"
        if 'G' not in aligned:
            # can use CB
            template += " and ( backbone or name CB )"
        else:
            template += " and backbone"
        template = "( " + template + " )"

        res_list.append([template % resid(iseq, ipos) for iseq in range(nseq)])

    sel = np.array(res_list).transpose()

    ref_selection = " or ".join(sel[0])
    target_selection = " or ".join(sel[1])
    return {'reference': ref_selection, 'mobile': target_selection}


def get_matching_atoms(ag1, ag2, tol_mass=0.1, strict=False, match_atoms=True):
    """Return two atom groups with one-to-one matched atoms.

    The function takes two :class:`~MDAnalysis.core.groups.AtomGroup`
    instances `ag1` and `ag2` and returns two atom groups `g1` and `g2` that
    consist of atoms so that the mass of atom ``g1[0]`` is the same as the mass
    of atom ``g2[0]``, ``g1[1]`` and ``g2[1]`` etc.

    The current implementation is very simplistic and works on a per-residue basis:

    1. The two groups must contain the same number of residues.
    2. Any residues in each group that have differing number of atoms are discarded.
    3. The masses of corresponding atoms are compared. and if any masses differ
       by more than `tol_mass` the test is considered failed and a
       :exc:`SelectionError` is raised.

    The log file (see :func:`MDAnalysis.start_logging`) will contain detailed
    information about mismatches.

    Parameters
    ----------
    ag1 : AtomGroup
        First :class:`~MDAnalysis.core.groups.AtomGroup` instance that is
        compared
    ag2 : AtomGroup
        Second :class:`~MDAnalysis.core.groups.AtomGroup` instance that is
        compared
    tol_mass : float (optional)
         Reject if the atomic masses for matched atoms differ by more than
         `tol_mass` [0.1]
    strict : bool (optional)
        ``True``
            Will raise :exc:`SelectionError` if a single atom does not
            match between the two selections.
        ``False`` [default]
            Will try to prepare a matching selection by dropping
            residues with non-matching atoms. See :func:`get_matching_atoms`
            for details.
    match_atoms : bool (optional)
        ``True``
            Will attempt to match atoms based on mass
        ``False``
            Will not attempt to match atoms based on mass

    Returns
    -------
    (g1, g2) : tuple
        Tuple with :class:`~MDAnalysis.core.groups.AtomGroup`
        instances that match, atom by atom. The groups are either the
        original groups if all matched or slices of the original
        groups.

    Raises
    ------
    :exc:`SelectionError`
        Error raised if the number of residues does not match or if in the final
        matching masses differ by more than *tol*.

    Notes
    -----
    The algorithm could be improved by using e.g. the Needleman-Wunsch
    algorithm in :mod:`Bio.profile2` to align atoms in each residue (doing a
    global alignment is too expensive).

    .. versionadded:: 0.8

    .. versionchanged:: 0.10.0
       Renamed from :func:`check_same_atoms` to
       :func:`get_matching_atoms` and now returns matching atomgroups
       (possibly with residues removed)

    """

    if ag1.n_atoms != ag2.n_atoms:
        if not match_atoms:
            errmsg = ("Mobile and reference atom selections do not "
                      "contain the same number of atoms and atom "
                      "matching is turned off. To match atoms based "
                      "on residue and mass, try match_atoms=True")
            logger.error(errmsg)
            raise SelectionError(errmsg)
        if ag1.n_residues != ag2.n_residues:
            errmsg = ("Reference and trajectory atom selections do not contain "
                    "the same number of atoms: \n"
                    "atoms:    N_ref={0}, N_traj={1}\n"
                    "and also not the same number of residues:\n"
                    "residues: N_ref={2}, N_traj={3}").format(
                        ag1.n_atoms, ag2.n_atoms,
                        ag1.n_residues, ag2.n_residues)
            logger.error(errmsg)
            raise SelectionError(errmsg)
        else:
            msg = ("Reference and trajectory atom selections do not contain "
                   "the same number of atoms: \n"
                   "atoms:    N_ref={0}, N_traj={1}").format(
                       ag1.n_atoms, ag2.n_atoms)
            if strict:
                logger.error(msg)
                raise SelectionError(msg)

            # continue with trying to create a valid selection
            msg += ("\nbut we attempt to create a valid selection " +
                    "(use strict=True to disable this heuristic).")
            logger.info(msg)
            warnings.warn(msg, category=SelectionWarning)

        # continue with trying to salvage the selection:
        # - number of atoms is different
        # - number of residues is the same
        # We will remove residues with mismatching number of atoms (e.g. not resolved
        # in an X-ray structure)
        assert ag1.n_residues == ag2.n_residues

        # Alternatively, we could align all atoms but Needleman-Wunsch
        # pairwise2 consumes too much memory for thousands of characters in
        # each sequence. Perhaps a solution would be pairwise alignment per residue.
        #
        # aln_elem = Bio.pairwise2.align.globalms("".join([MDAnalysis.topology.
        # core.guess_atom_element(n) for n in gref.atoms.names]),
        # "".join([MDAnalysis.topology.core.guess_atom_element(n)
        # for n in models[0].atoms.names]),
        # 2, -1, -1, -0.1,
        # one_alignment_only=True)

        # For now, just remove the residues that don't have matching numbers
        # NOTE: This can create empty selections, e.g., when comparing a structure
        #       with hydrogens to a PDB structure without hydrogens.
        rsize1 = np.array([r.atoms.n_atoms for r in ag1.residues])
        rsize2 = np.array([r.atoms.n_atoms for r in ag2.residues])
        rsize_mismatches = np.absolute(rsize1 - rsize2)
        mismatch_mask = (rsize_mismatches > 0)
        if np.any(mismatch_mask):
            def get_atoms_byres(g, match_mask=np.logical_not(mismatch_mask)):
                # not pretty... but need to do things on a per-atom basis in
                # order to preserve original selection
                ag = g.atoms
                good = ag.residues.resids[match_mask]  # resid for each residue
                resids = ag.resids                     # resid for each atom
                # boolean array for all matching atoms
                ix_good = np.isin(resids, good)
                return ag[ix_good]

            _ag1 = get_atoms_byres(ag1)
            _ag2 = get_atoms_byres(ag2)

            assert _ag1.atoms.n_atoms == _ag2.atoms.n_atoms

            # diagnostics
            mismatch_resindex = np.arange(ag1.n_residues)[mismatch_mask]
            logger.warning("Removed {0} residues with non-matching numbers of atoms"
                           .format(mismatch_mask.sum()))
            logger.debug("Removed residue ids: group 1: {0}"
                         .format(ag1.residues.resids[mismatch_resindex]))
            logger.debug("Removed residue ids: group 2: {0}"
                         .format(ag2.residues.resids[mismatch_resindex]))
            # replace after logging (still need old ag1 and ag2 for
            # diagnostics)
            ag1 = _ag1
            ag2 = _ag2
            del _ag1, _ag2

            # stop if we created empty selections (by removing ALL residues...)
            if ag1.n_atoms == 0 or ag2.n_atoms == 0:
                errmsg = ("Failed to automatically find matching atoms: created empty selections. "
                          "Try to improve your selections for mobile and reference.")
                logger.error(errmsg)
                raise SelectionError(errmsg)

    if match_atoms:
        # check again because the residue matching heuristic is not very
        # good and can easily be misled (e.g., when one of the selections
        # had fewer atoms but the residues in mobile and reference have
        # each the same number)
        if (not hasattr(ag1, 'masses') or not hasattr(ag2, 'masses')):
            msg = "Atoms could not be matched since they don't contain masses."
            logger.info(msg)
            warnings.warn(msg, category=SelectionWarning)
        else:
            try:
                mass_mismatches = (np.absolute(ag1.masses - ag2.masses) > tol_mass)
            except ValueError:
                errmsg = ("Failed to find matching atoms: len(reference) = {}, len(mobile) = {} "
                          "Try to improve your selections for mobile and reference.").format(
                            ag1.n_atoms, ag2.n_atoms)
                logger.error(errmsg)
                raise SelectionError(errmsg) from None

            if np.any(mass_mismatches):
                # Test 2 failed.
                # diagnostic output:
                logger.error("Atoms: reference | trajectory")
                for ar, at in zip(ag1[mass_mismatches], ag2[mass_mismatches]):
                    logger.error(
                        "{0!s:>4} {1:3d} {2!s:>3} {3!s:>3} {4:6.3f}  |  {5!s:>4} {6:3d} {7!s:>3} {8!s:>3} {9:6.3f}".format(
                            ar.segid,
                            ar.resid,
                            ar.resname,
                            ar.name,
                            ar.mass,
                            at.segid,
                            at.resid,
                            at.resname,
                            at.name,
                            at.mass))
                errmsg = ("Inconsistent selections, masses differ by more than {0}; "
                           "mis-matching atoms are shown above.").format(tol_mass)
                logger.error(errmsg)
                raise SelectionError(errmsg)

    return ag1, ag2
