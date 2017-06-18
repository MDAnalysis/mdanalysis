# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.mdanalysis.org
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
"""
Calculating root mean square quantities --- :mod:`MDAnalysis.analysis.rms`
==========================================================================

:Author: Oliver Beckstein, David L. Dotson, John Detlefs
:Year: 2016
:Copyright: GNU Public License v2

.. versionadded:: 0.7.7
.. versionchanged:: 0.11.0
   Added :class:`RMSF` analysis.
.. versionchanged:: 0.16.0
   Refactored RMSD to fit AnalysisBase API

The module contains code to analyze root mean square quantities such
as the coordinat root mean square distance (:class:`RMSD`) or the
per-residue root mean square fluctuations (:class:`RMSF`).

This module uses the fast QCP algorithm [Theobald2005]_ to calculate
the root mean square distance (RMSD) between two coordinate sets (as
implemented in
:func:`MDAnalysis.lib.qcprot.CalcRMSDRotationalMatrix`).

When using this module in published work please cite [Theobald2005]_.


See Also
--------
:mod:`MDAnalysis.analysis.align`
   aligning structures based on RMSD
:mod:`MDAnalysis.lib.qcprot`
   implements the fast RMSD algorithm.


Example applications
--------------------

Calculating RMSD for multiple domains
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this example we will globally fit a protein to a reference
structure and investigate the relative movements of domains by
computing the RMSD of the domains to the reference. The example is a
DIMS trajectory of adenylate kinase, which samples a large
closed-to-open transition. The protein consists of the CORE, LID, and
NMP domain.

* superimpose on the closed structure (frame 0 of the trajectory),
  using backbone atoms

* calculate the backbone RMSD and RMSD for CORE, LID, NMP (backbone atoms)

The trajectory is included with the test data files. The data in
:attr:`RMSD.rmsd` is plotted with :func:`matplotlib.pyplot.plot`::

   import MDAnalysis
   from MDAnalysis.tests.datafiles import PSF,DCD,CRD
   u = MDAnalysis.Universe(PSF,DCD)
   ref = MDAnalysis.Universe(PSF,DCD)     # reference closed AdK (1AKE) (with the default ref_frame=0)
   #ref = MDAnalysis.Universe(PSF,CRD)    # reference open AdK (4AKE)

   import MDAnalysis.analysis.rms

   R = MDAnalysis.analysis.rms.RMSD(u, ref,
              select="backbone",             # superimpose on whole backbone of the whole protein
              groupselections=["backbone and (resid 1-29 or resid 60-121 or resid 160-214)",   # CORE
                               "backbone and resid 122-159",                                   # LID
                               "backbone and resid 30-59"],                                    # NMP
              filename="rmsd_all_CORE_LID_NMP.dat")
   R.run()
   R.save()

   import matplotlib.pyplot as plt
   rmsd = R.rmsd.T   # transpose makes it easier for plotting
   time = rmsd[1]
   fig = plt.figure(figsize=(4,4))
   ax = fig.add_subplot(111)
   ax.plot(time, rmsd[2], 'k-',  label="all")
   ax.plot(time, rmsd[3], 'k--', label="CORE")
   ax.plot(time, rmsd[4], 'r--', label="LID")
   ax.plot(time, rmsd[5], 'b--', label="NMP")
   ax.legend(loc="best")
   ax.set_xlabel("time (ps)")
   ax.set_ylabel(r"RMSD ($\AA$)")
   fig.savefig("rmsd_all_CORE_LID_NMP_ref1AKE.pdf")


Functions
---------

.. autofunction:: rmsd

Analysis classes
----------------

.. autoclass:: RMSD
   :members:
   :inherited-members:

   .. attribute:: rmsd

       Contains the time series of the RMSD as an N×3 :class:`numpy.ndarray`
       array with content ``[[frame, time (ps), RMSD (A)], [...], ...]``.

.. autoclass:: RMSF
   :members:
   :inherited-members:

   .. attribute:: rmsf

      Results are stored in this N-length :class:`numpy.ndarray` array,
      giving RMSFs for each of the given atoms.

"""
from __future__ import division, absolute_import

from six.moves import zip
import numpy as np
import logging
import warnings


import MDAnalysis.lib.qcprot as qcp
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.exceptions import SelectionError, NoDataError
from MDAnalysis.lib.log import ProgressMeter, _set_verbose
from MDAnalysis.lib.util import asiterable


logger = logging.getLogger('MDAnalysis.analysis.rmsd')


def rmsd(a, b, weights=None, center=False, superposition=False):
    """Returns RMSD between two coordinate sets `a` and `b`.

    `a` and `b` are arrays of the coordinates of N atoms of shape N*3
    as generated by, e.g.,
    :meth:`MDAnalysis.core.groups.AtomGroup.coordinates`.

    .. note::
       If you use trajectory data from simulations performed under **periodic
       boundary conditions** then you *must make your molecules whole* before
       performing RMSD calculations so that the centers of mass of the mobile
       and reference structure are properly superimposed.


    Parameters
    ----------
    a, b : array_like
        coordinates to align
    weights : array_like (optional)
        1D array with weights, use to compute weighted average
    center : bool (optional)
        subtract center of geometry before calculation. With weights given
        compute weighted average as center.
    superposition : bool (optional)
        perform a rotational and translational superposition with the fast QCP
        algorithm [Theobald2005]_ before calculating the RMSD

    Returns
    -------
    rmsd : float
        RMSD between a and b

    Example
    -------
    >>> u = Universe(PSF,DCD)
    >>> bb = u.select_atoms('backbone')
    >>> A = bb.positions.copy()  # coordinates of first frame
    >>> u.trajectory[-1]         # forward to last frame
    >>> B = bb.positions.copy()  # coordinates of last frame
    >>> rmsd(A, B, center=True)
    3.9482355416565049

    .. versionchanged: 0.8.1
       *center* keyword added
    .. versionchanged: 0.14.0
       *superposition* keyword added

    """

    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    N = b.shape[0]

    if a.shape != b.shape:
        raise ValueError('a and b must have same shape')

    # superposition only works if structures are centered
    if center or superposition:
        # make copies (do not change the user data!)
        # weights=None is equivalent to all weights 1
        a = a - np.average(a, axis=0, weights=weights)
        b = b - np.average(b, axis=0, weights=weights)

    if weights is not None:
        if len(weights) != len(a):
            raise ValueError('weights must have same length as a/b')
        # weights are constructed as relative to the mean
        weights = np.asarray(weights, dtype=np.float64) / np.mean(weights)

    if superposition:
        return qcp.CalcRMSDRotationalMatrix(a, b, N, None, weights)
    else:
        if weights is not None:
            return np.sqrt(np.sum(weights[:, np.newaxis]
                                  * ((a - b) ** 2)) / N)
        else:
            return np.sqrt(np.sum((a - b) ** 2) / N)


def process_selection(select):
    """Return a canonical selection dictionary.

    Parameters
    ----------
    select : str / tuple / dict

        - `str` -> Any valid string selection
        - `dict` -> ``{'mobile':sel1, 'reference':sel2}``
        - `tuple` -> ``(sel1, sel2)``

    Returns
    -------
    dict
        selections for 'reference' and 'mobile'. Values are guarenteed to be
        iterable (so that one can provide selections to retain order)

    Notes
    -----
    The dictionary input for `select` can be generated by
    :func:`fasta2select` based on a ClustalW_ or STAMP_ sequence alignment.
    """
    if type(select) is str:
        select = {'reference': select, 'mobile': select}
    elif type(select) is tuple:
        try:
            select = {'mobile': select[0], 'reference': select[1]}
        except IndexError:
            raise IndexError("select must contain two selection strings "
                             "(reference, mobile)")
    elif type(select) is dict:
        # compatability hack to use new nomenclature
        try:
            select['mobile'] = select['target']
            warnings.warn("use key 'mobile' instead of deprecated 'target'; "
                          "'target' will be removed in 0.8",
                          DeprecationWarning)
        except KeyError:
            pass
        try:
            select['mobile']
            select['reference']
        except KeyError:
            raise KeyError("select dictionary must contain entries for keys "
                           "'mobile' and 'reference'.")
    else:
        raise TypeError("'select' must be either a string, 2-tuple, or dict")
    select['mobile'] = asiterable(select['mobile'])
    select['reference'] = asiterable(select['reference'])
    return select


class RMSD(AnalysisBase):
    r"""Class to perform RMSD analysis on a trajectory.

    The RMSD will be computed between `select` and `reference` for
    all frames in the trajectory belonging to `atomgroup`.

    .. Note::
       If you use trajectory data from simulations performed under **periodic
       boundary conditions** then you *must make your molecules whole* before
       performing RMSD calculations so that the centers of mass of the mobile
       and reference structure are properly superimposed.

    Run the analysis with :meth:`RMSD.run`, which stores the results
    in the array :attr:`RMSD.rmsd`.


    Parameters
    ----------
    atomgroup : :class:`~MDAnalysis.core.groups.AtomGroup` or
                :class:`~MDAnalysis.core.universe.Universe`
        MDAnalysis `AtomGroup` or `Universe` with an associated trajectory
        to be RMSD fit.
    reference : :class:`~MDAnalysis.core.groups.AtomGroup` or
                :class:`~MDAnalysis.core.universe.Universe`, optional
        Reference coordinates, if ``None`` current frame of `atomgroup` is
        used
    select : str / dict / tuple (optional)
        The selection to operate on; can be one of:

        1. any valid selection string for
           :meth:`~MDAnalysis.core.groups.AtomGroup.select_atoms` that
           produces identical selections in *mobile* and *reference*; or

        2. a dictionary ``{'mobile':sel1, 'reference':sel2}`` (the
           :func:`MDAnalysis.analysis.align.fasta2select` function returns
           such a dictionary based on a ClustalW_ or STAMP_ sequence
           alignment); or
        3. a tuple ``(sel1, sel2)``

        When using 2. or 3. with *sel1* and *sel2* then these selections
        can also each be a list of selection strings (to generate a
        AtomGroup with defined atom order as described under
        :ref:`ordered-selections-label`).
    groupselections : list (optional)
        A list of selections as described for `select` Each selection
        describes additional RMSDs to be computed *after the structures
        have be superpositioned* according to `select`. The output contains
        one additional column for each selection. [``None``]

        .. Note:: Experimental feature. Only limited error checking
                  implemented.

    filename : str, optional
        write RSMD into file file :meth:`RMSD.save`
    mass_weighted : bool (deprecated)
         do a mass-weighted RMSD fit
    weights : str/array_like (optional)
         choose weights. If 'str' uses masses as weights
    tol_mass : float (optional)
         Reject match if the atomic masses for matched atoms differ by more
         than `tol_mass`
    ref_frame : int (optional)
         frame index to select frame from `reference`

    Notes
    -----
    The root mean square deviation of a group of :math:`N` atoms relative to a
    reference structure as a function of time is calculated as

    .. math::

       \rho(t) = \sqrt{\frac{1}{N} \sum_{i=1}^N \left(\mathbf{x}_i(t)
                                - \mathbf{x}_i^{\text{ref}}\right)^2}

    This class uses Douglas Theobald's fast QCP algorithm [Theobald2005]_ to
    calculate the RMSD (see :mod:`MDAnalysis.lib.qcprot` for implementation
    details).

    The class runs various checks on the input to ensure that the two atom
    groups can be compared. This includes a comparison of atom masses (i.e.,
    only the positions of atoms of the same mass will be considered to be
    correct for comparison). If masses should not be checked, just set
    `tol_mass` to a large value such as 1000.

    .. _ClustalW: http://www.clustal.org/
    .. _STAMP: http://www.compbio.dundee.ac.uk/manuals/stamp.4.2/

    .. versionadded:: 0.7.7
    .. versionchanged:: 0.8
       `groupselections` added
    .. versionchanged:: 0.16.0
       Flexible weighting scheme with new `weights` keyword.
    .. deprecated:: 0.16.0
       Instead of ``mass_weighted=True`` use new ``weights='mass'``;
       refactored to fit with AnalysisBase API

    """

    def __init__(self, atomgroup, reference=None, select='all',
                 groupselections=None, filename="rmsd.dat",
                 mass_weighted=None,
                 weights=None, tol_mass=0.1, ref_frame=0, **kwargs):
        super(RMSD, self).__init__(atomgroup.universe.trajectory,
                                   **kwargs)
        self.universe = atomgroup.universe
        self.reference = reference if reference is not None else self.universe

        select = process_selection(select)
        self.groupselections = ([process_selection(s) for s in groupselections]
                                if groupselections is not None else [])
        if mass_weighted is not None:
            warnings.warn("mass weighted is deprecated argument. Please use "
                          " 'weights=\"mass\" instead. Will be removed in 0.17.0",
                          category=DeprecationWarning)
            if mass_weighted:
                weights = 'mass'
        self.weights = weights
        self.tol_mass = tol_mass
        self.ref_frame = ref_frame
        self.filename = filename

        self.ref_atoms = self.reference.select_atoms(*select['reference'])
        self.mobile_atoms = self.universe.select_atoms(*select['mobile'])

        if len(self.ref_atoms) != len(self.mobile_atoms):
            err = ("Reference and trajectory atom selections do "
                   "not contain the same number of atoms: "
                   "N_ref={0:d}, N_traj={1:d}".format(self.ref_atoms.n_atoms,
                                                      self.mobile_atoms.n_atoms))
            logger.exception(err)
            raise SelectionError(err)
        logger.info("RMS calculation "
                    "for {0:d} atoms.".format(len(self.ref_atoms)))
        mass_mismatches = (np.absolute((self.ref_atoms.masses -
                                        self.mobile_atoms.masses)) >
                           self.tol_mass)

        if np.any(mass_mismatches):
            # diagnostic output:
            logger.error("Atoms: reference | mobile")
            for ar, at in zip(self.ref_atoms, self.mobile_atoms):
                if ar.name != at.name:
                    logger.error("{0!s:>4} {1:3d} {2!s:>3} {3!s:>3} {4:6.3f}"
                                 "|  {5!s:>4} {6:3d} {7!s:>3} {8!s:>3}"
                                 "{9:6.3f}".format(ar.segid, ar.resid,
                                                   ar.resname, ar.name,
                                                   ar.mass, at.segid, at.resid,
                                                   at.resname, at.name,
                                                   at.mass))
            errmsg = ("Inconsistent selections, masses differ by more than"
                      "{0:f}; mis-matching atoms"
                      "are shown above.".format(self.tol_mass))
            logger.error(errmsg)
            raise SelectionError(errmsg)
        del mass_mismatches

        # TODO:
        # - make a group comparison a class that contains the checks above
        # - use this class for the *select* group and the additional
        #   *groupselections* groups each a dict with reference/mobile
        self._groupselections_atoms = [
            {
                'reference': self.reference.select_atoms(*s['reference']),
                'mobile': self.universe.select_atoms(*s['mobile']),
            }
            for s in self.groupselections]
        # sanity check
        for igroup, (sel, atoms) in enumerate(zip(self.groupselections,
                                                  self._groupselections_atoms)):
            if len(atoms['mobile']) != len(atoms['reference']):
                logger.exception('SelectionError: Group Selection')
                raise SelectionError(
                    "Group selection {0}: {1} | {2}: Reference and trajectory "
                    "atom selections do not contain the same number of atoms: "
                    "N_ref={3}, N_traj={4}".format(
                        igroup, sel['reference'], sel['mobile'],
                        len(atoms['reference']), len(atoms['mobile'])))
        # initialized to note for testing the save function
        self.rmsd = None


    def _prepare(self):
        self._n_atoms = self.mobile_atoms.n_atoms

        if not isinstance(self.weights, (list, tuple, np.ndarray)) and self.weights == 'mass':
            self.weights = self.ref_atoms.masses
        if self.weights is not None:
            self.weights = (self.weights / self.weights.mean()).astype(np.float64)

        current_frame = self.reference.trajectory.ts.frame

        try:
            # Move to the ref_frame
            # (coordinates MUST be stored in case the ref traj is advanced
            # elsewhere or if ref == mobile universe)
            self.reference.trajectory[self.ref_frame]
            self._ref_com = self.ref_atoms.center(self.weights)
            # makes a copy
            self._ref_coordinates = self.ref_atoms.positions - self._ref_com
            if self._groupselections_atoms:
                self._groupselections_ref_coords_64 = [(self.reference.
                    select_atoms(*s['reference']).
                    positions.astype(np.float64)) for s in
                    self.groupselections]
        finally:
            # Move back to the original frame
            self.reference.trajectory[current_frame]

        self._ref_coordinates_64 = self._ref_coordinates.astype(np.float64)

        if self._groupselections_atoms:
            # Only carry out a rotation if we want to calculate secondary
            # RMSDs.
            # R: rotation matrix that aligns r-r_com, x~-x~com
            #    (x~: selected coordinates, x: all coordinates)
            # Final transformed traj coordinates: x' = (x-x~_com)*R + ref_com
            self._rot = np.zeros(9, dtype=np.float64)  # allocate space
            self._R = np.matrix(self._rot.reshape(3, 3))
        else:
            self._rot = None

        self.rmsd = np.zeros((self.n_frames,
                              3 + len(self._groupselections_atoms)))

        self._pm.format = ("RMSD {rmsd:5.2f} A at frame "
                           "{step:5d}/{numsteps}  [{percentage:5.1f}%]\r")
        self._mobile_coordinates64 = self.mobile_atoms.positions.copy().astype(np.float64)

    def _single_frame(self):
        mobile_com = self.mobile_atoms.center(self.weights).astype(np.float64)
        self._mobile_coordinates64[:] = self.mobile_atoms.positions
        self._mobile_coordinates64 -= mobile_com

        self.rmsd[self._frame_index, :2] = self._ts.frame, self._trajectory.time

        if self._groupselections_atoms:
            # superimpose structures: MDAnalysis qcprot needs Nx3 coordinate
            # array with float64 datatype (float32 leads to errors up to 1e-3 in
            # RMSD). Note that R is defined in such a way that it acts **to the
            # left** so that we can easily use broadcasting and save one
            # expensive numpy transposition.

            self.rmsd[self._frame_index, 2] = qcp.CalcRMSDRotationalMatrix(
                self._ref_coordinates_64, self._mobile_coordinates64,
                self._n_atoms, self._rot, self.weights)

            self._R[:, :] = self._rot.reshape(3, 3)
            # Transform each atom in the trajectory (use inplace ops to
            # avoid copying arrays) (Marginally (~3%) faster than
            # "ts.positions[:] = (ts.positions - x_com) * R + ref_com".)
            self._ts.positions[:] -= mobile_com

            # R acts to the left & is broadcasted N times.
            self._ts.positions[:,:] = (self._mobile_coordinates64[:] *
                                       self._R)
            self._ts.positions[:] += self._ref_com

            # 2) calculate secondary RMSDs
            for igroup, (refpos, atoms) in enumerate(
                    zip(self._groupselections_ref_coords_64,
                        self._groupselections_atoms), 3):
                self.rmsd[self._frame_index, igroup] = qcp.CalcRMSDRotationalMatrix(
                    refpos, atoms['mobile'].positions.astype(np.float64),
                    atoms['mobile'].n_atoms, None, self.weights)
        else:
            # only calculate RMSD by setting the Rmatrix to None (no need
            # to carry out the rotation as we already get the optimum RMSD)
            self.rmsd[self._frame_index, 2] = qcp.CalcRMSDRotationalMatrix(
                self._ref_coordinates_64, self._mobile_coordinates64,
                self._n_atoms, None, self.weights)

        self._pm.rmsd = self.rmsd[self._frame_index, 2]

    def save(self, filename=None):
        """Save RMSD from :attr:`RMSD.rmsd` to text file *filename*.

        Parameters
        ----------
        filename : str (optional)
            if no filename is given the default provided to the constructor is
            used.
        """
        filename = filename or self.filename
        if filename is not None:
            if self.rmsd is None:
                raise NoDataError("rmsd has not been calculated yet")
            np.savetxt(filename, self.rmsd)
            logger.info("Wrote RMSD timeseries  to file %r", filename)
        return filename


class RMSF(AnalysisBase):
    r"""Calculate RMSF of given atoms across a trajectory.

    Parameters
    ----------
    atomgroup : AtomGroup
        Atoms for which RMSF is calculated
    start : int (optional)
        starting frame, default None becomes 0.
    stop : int (optional)
        Frame index to stop analysis. Default: None becomes
        n_frames. Iteration stops *before* this frame number,
        which means that the trajectory would be read until the end.
    step : int (optional)
        step between frames, default None becomes 1.
    weights : str / array_like (optional)
        used weights. If ``'mass'`` use masses of atomgroup, it ``None``
        use uniform weights.
    verbose : bool (optional)
        if ``False``, suppress all output

    Note
    ----
    No RMSD-superposition is performed; it is assumed that the user is
    providing a trajectory where the protein of interest has been structurally
    aligned to a reference structure.


    Notes
    -----
    The root mean square fluctuation of an atom :math:`i` is computed as the
    time average

    .. math::

      \rho_i = \sqrt{\left\langle (\mathbf{x}_i - \langle\mathbf{x}_i\rangle)^2 \right\rangle}

    This method implements an algorithm for computing sums of squares while
    avoiding overflows and underflows [Welford1962]_.

    References
    ----------
    .. [Welford1962] B. P. Welford (1962). "Note on a Method for
       Calculating Corrected Sums of Squares and Products." Technometrics
       4(3):419-420.

    .. versionadded:: 0.11.0
    .. versionchanged:: 0.16.0
       Flexible weighting scheme with new `weights` keyword.
    .. deprecated:: 0.16.0
       Instead of ``mass_weighted=True`` use new ``weights='mass'``;
       refactored to fit with AnalysisBase API;
       the keyword argument `quiet` is deprecated in favor of `verbose`.

    """
    def __init__(self, atomgroup, weights=None, **kwargs):
        super(RMSF, self).__init__(atomgroup.universe.trajectory, **kwargs)
        self.atomgroup = atomgroup
        if not isinstance(weights, (list, tuple, np.ndarray)) and weights == 'mass':
            weights = self.atomgroup.masses
        self.weights = weights

    def run(self, start=None, stop=None, step=None, progout=None,
            verbose=None, quiet=None):
        if any([el is not None for el in (start, stop, step, progout, quiet)]):
            warnings.warn("run arguments are deprecated. Please pass them at "
                          "class construction. These options will be removed in 0.17.0",
                          category=DeprecationWarning)
            verbose = _set_verbose(verbose, quiet, default=False)
            # regenerate class with correct args
            super(RMSF, self).__init__(self.atomgroup.universe.trajectory,
                                       start=start, stop=stop, step=step,
                                       verbose=verbose)
        return super(RMSF, self).run()

    def _prepare(self):
        self.sumsquares = np.zeros((self.atomgroup.n_atoms, 3))
        self.mean = self.sumsquares.copy()

    def _single_frame(self):
        k = self._frame_index
        self.sumsquares += (k / (k+1.0)) * (self.atomgroup.positions - self.mean) ** 2
        self.mean = (k * self.mean + self.atomgroup.positions) / (k + 1)

    def _conclude(self):
        k = self._frame_index
        self.rmsf = np.sqrt(self.sumsquares.sum(axis=1) / (k + 1))

        if not (self.rmsf >= 0).all():
            raise ValueError("Some RMSF values negative; overflow " +
                             "or underflow occurred")
