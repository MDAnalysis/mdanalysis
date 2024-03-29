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
:attr:`RMSD.results.rmsd` is plotted with :func:`matplotlib.pyplot.plot` (see Figure :ref:`RMSD plot figure <figure-RMSD>`)::

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
                               "backbone and resid 30-59"])                                    # NMP
   R.run()

   import matplotlib.pyplot as plt
   rmsd = R.results.rmsd.T   # transpose makes it easier for plotting
   time = rmsd[1]
   fig = plt.figure(figsize=(4,4))
   ax = fig.add_subplot(111)
   ax.plot(time, rmsd[2], 'k-',  label="all")
   ax.plot(time, rmsd[3], 'k--', label="CORE")
   ax.plot(time, rmsd[4], 'r--', label="LID")
   ax.plot(time, rmsd[5], 'b--', label="NMP")
   ax.legend(loc="best")
   ax.set_xlabel("time (ps)")
   ax.set_ylabel(r"RMSD ($\\AA$)")
   fig.savefig("rmsd_all_CORE_LID_NMP_ref1AKE.pdf")
   
.. _figure-RMSD:

.. figure:: /images/RSMD_plot.png
      :scale: 50 %
      :alt: RMSD plot
      
      RMSD plot for backbone and CORE, LID, NMP domain of the protein.

      
Functions
---------

.. autofunction:: rmsd

Analysis classes
----------------

.. autoclass:: RMSD
   :members:
   :inherited-members:

   .. attribute:: results.rmsd

       Contains the time series of the RMSD as an N×3 :class:`numpy.ndarray`
       array with content ``[[frame, time (ps), RMSD (A)], [...], ...]``.

       .. versionadded:: 2.0.0

   .. attribute:: rmsd

       Alias to the :attr:`results.rmsd` attribute.

       .. deprecated:: 2.0.0
          Will be removed in MDAnalysis 3.0.0. Please use :attr:`results.rmsd`
          instead.


.. autoclass:: RMSF
   :members:
   :inherited-members:

   .. attribute:: results.rmsf

      Results are stored in this N-length :class:`numpy.ndarray` array,
      giving RMSFs for each of the given atoms.

      .. versionadded:: 2.0.0

   .. attribute:: rmsf

      Alias to the :attr:`results.rmsf` attribute.

      .. deprecated:: 2.0.0
         Will be removed in MDAnalysis 3.0.0. Please use :attr:`results.rmsf`
         instead.

"""
import numpy as np

import logging
import warnings

import MDAnalysis as mda
import MDAnalysis.lib.qcprot as qcp
from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.exceptions import SelectionError, NoDataError
from MDAnalysis.lib.util import asiterable, iterable, get_weights


logger = logging.getLogger('MDAnalysis.analysis.rmsd')


def rmsd(a, b, weights=None, center=False, superposition=False):
    r"""Returns RMSD between two coordinate sets `a` and `b`.

    `a` and `b` are arrays of the coordinates of N atoms of shape
    :math:`N times 3` as generated by, e.g.,
    :meth:`MDAnalysis.core.groups.AtomGroup.positions`.

    Note
    ----
    If you use trajectory data from simulations performed under **periodic
    boundary conditions** then you *must make your molecules whole* before
    performing RMSD calculations so that the centers of mass of the mobile and
    reference structure are properly superimposed.


    Parameters
    ----------
    a : array_like
        coordinates to align to `b`
    b : array_like
        coordinates to align to (same shape as `a`)
    weights : array_like (optional)
        1D array with weights, use to compute weighted average
    center : bool (optional)
        subtract center of geometry before calculation. With weights given
        compute weighted average as center.
    superposition : bool (optional)
        perform a rotational and translational superposition with the fast QCP
        algorithm [Theobald2005]_ before calculating the RMSD; implies
        ``center=True``.

    Returns
    -------
    rmsd : float
        RMSD between `a` and `b`

    Notes
    -----
    The RMSD :math:`\rho(t)` as a function of time is calculated as

    .. math::

       \rho(t) = \sqrt{\frac{1}{N} \sum_{i=1}^N w_i \left(\mathbf{x}_i(t)
                                - \mathbf{x}_i^{\text{ref}}\right)^2}

    It is the Euclidean distance in configuration space of the current
    configuration (possibly after optimal translation and rotation) from a
    reference configuration divided by :math:`1/\sqrt{N}` where :math:`N` is
    the number of coordinates.

    The weights :math:`w_i` are calculated from the input weights
    `weights` :math:`w'_i` as relative to the mean:

    .. math::

       w_i = \frac{w'_i}{\langle w' \rangle}


    Example
    -------
    >>> import MDAnalysis as mda
    >>> from MDAnalysis.analysis.rms import rmsd
    >>> from MDAnalysis.tests.datafiles import PSF, DCD
    >>> u = mda.Universe(PSF, DCD)
    >>> bb = u.select_atoms('backbone')
    >>> A = bb.positions.copy()  # coordinates of first frame
    >>> _ = u.trajectory[-1]  # forward to last frame
    >>> B = bb.positions.copy()  # coordinates of last frame
    >>> rmsd(A, B, center=True)
    6.838544558398293


    .. versionchanged:: 0.8.1
       *center* keyword added
    .. versionchanged:: 0.14.0
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
            raise ValueError('weights must have same length as a and b')
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
    select : str or tuple or dict

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

    if isinstance(select, str):
        select = {'reference': str(select), 'mobile': str(select)}
    elif type(select) is tuple:
        try:
            select = {'mobile': select[0], 'reference': select[1]}
        except IndexError:
            raise IndexError(
                "select must contain two selection strings "
                "(reference, mobile)") from None
    elif type(select) is dict:
        # compatability hack to use new nomenclature
        try:
            select['mobile']
            select['reference']
        except KeyError:
            raise KeyError(
                    "select dictionary must contain entries for keys "
                    "'mobile' and 'reference'."
                    ) from None
    else:
        raise TypeError("'select' must be either a string, 2-tuple, or dict")
    select['mobile'] = asiterable(select['mobile'])
    select['reference'] = asiterable(select['reference'])
    return select


class RMSD(AnalysisBase):
    r"""Class to perform RMSD analysis on a trajectory.

    The RMSD will be computed for two groups of atoms and all frames in the
    trajectory belonging to `atomgroup`. The groups of atoms are obtained by
    applying the selection selection `select` to the changing `atomgroup` and
    the fixed `reference`.

    Note
    ----
    If you use trajectory data from simulations performed under **periodic
    boundary conditions** then you *must make your molecules whole* before
    performing RMSD calculations so that the centers of mass of the selected
    and reference structure are properly superimposed.


    Run the analysis with :meth:`RMSD.run`, which stores the results
    in the array :attr:`RMSD.results.rmsd`.


    .. versionchanged:: 1.0.0
       ``save()`` method was removed, use ``np.savetxt()`` on
       :attr:`RMSD.results.rmsd` instead.
    .. versionchanged:: 2.0.0
       :attr:`rmsd` results are now stored in a
       :class:`MDAnalysis.analysis.base.Results` instance.

    """
    def __init__(self, atomgroup, reference=None, select='all',
                 groupselections=None, weights=None, weights_groupselections=False,
                 tol_mass=0.1, ref_frame=0, **kwargs):
        r"""Parameters
        ----------
        atomgroup : AtomGroup or Universe
            Group of atoms for which the RMSD is calculated. If a trajectory is
            associated with the atoms then the computation iterates over the
            trajectory.
        reference : AtomGroup or Universe (optional)
            Group of reference atoms; if ``None`` then the current frame of
            `atomgroup` is used.
        select : str or dict or tuple (optional)
            The selection to operate on; can be one of:

            1. any valid selection string for
               :meth:`~MDAnalysis.core.groups.AtomGroup.select_atoms` that
               produces identical selections in `atomgroup` and `reference`; or

            2. a dictionary ``{'mobile': sel1, 'reference': sel2}`` where *sel1*
               and *sel2* are valid selection strings that are applied to
               `atomgroup` and `reference` respectively (the
               :func:`MDAnalysis.analysis.align.fasta2select` function returns such
               a dictionary based on a ClustalW_ or STAMP_ sequence alignment); or

            3. a tuple ``(sel1, sel2)``

            When using 2. or 3. with *sel1* and *sel2* then these selection strings
            are applied to `atomgroup` and `reference` respectively and should
            generate *groups of equivalent atoms*.  *sel1* and *sel2* can each also
            be a *list of selection strings* to generate a
            :class:`~MDAnalysis.core.groups.AtomGroup` with defined atom order as
            described under :ref:`ordered-selections-label`).

        groupselections : list (optional)
            A list of selections as described for `select`, with the difference
            that these selections are *always applied to the full universes*,
            i.e., ``atomgroup.universe.select_atoms(sel1)`` and
            ``reference.universe.select_atoms(sel2)``. Each selection describes
            additional RMSDs to be computed *after the structures have been
            superimposed* according to `select`. No additional fitting is
            performed.The output contains one additional column for each
            selection.

            .. Note:: Experimental feature. Only limited error checking
                      implemented.

        weights : {"mass", ``None``} or array_like (optional)
             1. "mass" will use masses as weights for both `select` and `groupselections`.

             2. ``None`` will weigh each atom equally for both `select` and `groupselections`.

             3. If 1D float array of the same length as `atomgroup` is provided,
             use each element of the `array_like` as a weight for the
             corresponding atom in `select`, and assumes ``None`` for `groupselections`.

        weights_groupselections : False or list of {"mass", ``None`` or array_like} (optional)
             1. ``False`` will apply imposed weights to `groupselections` from
             ``weights`` option if ``weights`` is either ``"mass"`` or ``None``. 
             Otherwise will assume a list of length equal to length of 
             `groupselections` filled with ``None`` values.

             2. A list of {"mass", ``None`` or array_like} with the length of `groupselections`
             will apply the weights to `groupselections` correspondingly.

        tol_mass : float (optional)
             Reject match if the atomic masses for matched atoms differ by more
             than `tol_mass`.
        ref_frame : int (optional)
             frame index to select frame from `reference`
        verbose : bool (optional)
             Show detailed progress of the calculation if set to ``True``; the
             default is ``False``.

        Raises
        ------
        SelectionError
             If the selections from `atomgroup` and `reference` do not match.
        TypeError
             If `weights` or `weights_groupselections` is not of the appropriate type;
             see also :func:`MDAnalysis.lib.util.get_weights`
        ValueError
             If `weights` are not compatible with `atomgroup` (not the same
             length) or if it is not a 1D array (see
             :func:`MDAnalysis.lib.util.get_weights`).

             A :exc:`ValueError` is also raised if the length of `weights_groupselections`
             are not compatible with `groupselections`.

        Notes
        -----
        The root mean square deviation :math:`\rho(t)` of a group of :math:`N`
        atoms relative to a reference structure as a function of time is
        calculated as

        .. math::

           \rho(t) = \sqrt{\frac{1}{N} \sum_{i=1}^N w_i \left(\mathbf{x}_i(t)
                                    - \mathbf{x}_i^{\text{ref}}\right)^2}

        The weights :math:`w_i` are calculated from the input weights `weights`
        :math:`w'_i` as relative to the mean of the input weights:

        .. math::

           w_i = \frac{w'_i}{\langle w' \rangle}

        The selected coordinates from `atomgroup` are optimally superimposed
        (translation and rotation) on the `reference` coordinates at each time step
        as to minimize the RMSD. Douglas Theobald's fast QCP algorithm
        [Theobald2005]_ is used for the rotational superposition and to calculate
        the RMSD (see :mod:`MDAnalysis.lib.qcprot` for implementation details).

        The class runs various checks on the input to ensure that the two atom
        groups can be compared. This includes a comparison of atom masses (i.e.,
        only the positions of atoms of the same mass will be considered to be
        correct for comparison). If masses should not be checked, just set
        `tol_mass` to a large value such as 1000.

        .. _ClustalW: http://www.clustal.org/
        .. _STAMP: http://www.compbio.dundee.ac.uk/manuals/stamp.4.2/


        See Also
        --------
        rmsd


        .. versionadded:: 0.7.7
        .. versionchanged:: 0.8
           `groupselections` added
        .. versionchanged:: 0.16.0
           Flexible weighting scheme with new `weights` keyword.
        .. deprecated:: 0.16.0
           Instead of ``mass_weighted=True`` (removal in 0.17.0) use new
           ``weights='mass'``; refactored to fit with AnalysisBase API
        .. versionchanged:: 0.17.0
           removed deprecated `mass_weighted` keyword; `groupselections`
           are *not* rotationally superimposed any more.
        .. versionchanged:: 1.0.0
           `filename` keyword was removed.

        """
        super(RMSD, self).__init__(atomgroup.universe.trajectory,
                                   **kwargs)
        self.atomgroup = atomgroup
        self.reference = reference if reference is not None else self.atomgroup

        select = process_selection(select)
        self.groupselections = ([process_selection(s) for s in groupselections]
                                if groupselections is not None else [])
        self.weights = weights
        self.tol_mass = tol_mass
        self.ref_frame = ref_frame
        self.weights_groupselections = weights_groupselections
        self.ref_atoms = self.reference.select_atoms(*select['reference'])
        self.mobile_atoms = self.atomgroup.select_atoms(*select['mobile'])

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
                'reference': self.reference.universe.select_atoms(*s['reference']),
                'mobile': self.atomgroup.universe.select_atoms(*s['mobile']),
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

        # check weights type
        acceptable_dtypes = (np.dtype('float64'), np.dtype('int64'))
        msg = ("weights should only be 'mass', None or 1D float array."
              "For weights on groupselections, "
              "use **weight_groupselections**")

        if iterable(self.weights):
            element_lens = []
            for element in self.weights:
                if iterable(element):
                    element_lens.append(len(element))
                else:
                    element_lens.append(1)
            if np.unique(element_lens).size > 1:
                # jagged data structure
                raise TypeError(msg)
            if np.array(element).dtype not in acceptable_dtypes:
                raise TypeError(msg)

        if iterable(self.weights) or self.weights != "mass":
            get_weights(self.mobile_atoms, self.weights)

        if self.weights_groupselections:
            if len(self.weights_groupselections) != len(self.groupselections):
                raise ValueError("Length of weights_groupselections is not equal to "
                                 "length of groupselections ")
            for weights, atoms, selection in zip(self.weights_groupselections,
                                                 self._groupselections_atoms,
                                                 self.groupselections):
                try:
                    if iterable(weights) or weights != "mass":
                        get_weights(atoms['mobile'], weights)
                except Exception as e:
                    raise type(e)(str(e) + ' happens in selection %s' % selection['mobile'])


    def _prepare(self):
        self._n_atoms = self.mobile_atoms.n_atoms
        if not self.weights_groupselections:
            if not iterable(self.weights):         # apply 'mass' or 'None' to groupselections
                self.weights_groupselections = [self.weights] * len(self.groupselections)
            else:
                self.weights_groupselections = [None] * len(self.groupselections)

        for igroup, (weights, atoms) in enumerate(zip(self.weights_groupselections,
                                                      self._groupselections_atoms)):
            if str(weights) == 'mass':
                self.weights_groupselections[igroup] = atoms['mobile'].masses
            if weights is not None:
                self.weights_groupselections[igroup] = np.asarray(self.weights_groupselections[igroup],
                                                                  dtype=np.float64) /  \
                                             np.mean(self.weights_groupselections[igroup])
        # add the array of weights to weights_select
        self.weights_select = get_weights(self.mobile_atoms, self.weights)
        self.weights_ref = get_weights(self.ref_atoms, self.weights)
        if self.weights_select is not None:
            self.weights_select = np.asarray(self.weights_select, dtype=np.float64) /  \
                                  np.mean(self.weights_select)
            self.weights_ref = np.asarray(self.weights_ref, dtype=np.float64) / \
                               np.mean(self.weights_ref)

        current_frame = self.reference.universe.trajectory.ts.frame

        try:
            # Move to the ref_frame
            # (coordinates MUST be stored in case the ref traj is advanced
            # elsewhere or if ref == mobile universe)
            self.reference.universe.trajectory[self.ref_frame]
            self._ref_com = self.ref_atoms.center(self.weights_ref)
            # makes a copy
            self._ref_coordinates = self.ref_atoms.positions - self._ref_com
            if self._groupselections_atoms:
                self._groupselections_ref_coords64 = [(self.reference.
                                                       select_atoms(*s['reference']).
                                                       positions.astype(np.float64)) for s in
                                                     self.groupselections]
        finally:
            # Move back to the original frame
            self.reference.universe.trajectory[current_frame]

        self._ref_coordinates64 = self._ref_coordinates.astype(np.float64)

        if self._groupselections_atoms:
            # Only carry out a rotation if we want to calculate secondary
            # RMSDs.
            # R: rotation matrix that aligns r-r_com, x~-x~com
            #    (x~: selected coordinates, x: all coordinates)
            # Final transformed traj coordinates: x' = (x-x~_com)*R + ref_com
            self._rot = np.zeros(9, dtype=np.float64)  # allocate space
            self._R = self._rot.reshape(3, 3)
        else:
            self._rot = None

        self.results.rmsd = np.zeros((self.n_frames,
                                      3 + len(self._groupselections_atoms)))

        self._mobile_coordinates64 = self.mobile_atoms.positions.copy().astype(np.float64)

    def _single_frame(self):
        mobile_com = self.mobile_atoms.center(self.weights_select).astype(np.float64)
        self._mobile_coordinates64[:] = self.mobile_atoms.positions
        self._mobile_coordinates64 -= mobile_com

        self.results.rmsd[self._frame_index, :2] = self._ts.frame, self._trajectory.time

        if self._groupselections_atoms:
            # superimpose structures: MDAnalysis qcprot needs Nx3 coordinate
            # array with float64 datatype (float32 leads to errors up to 1e-3 in
            # RMSD). Note that R is defined in such a way that it acts **to the
            # left** so that we can easily use broadcasting and save one
            # expensive numpy transposition.

            self.results.rmsd[self._frame_index, 2] = qcp.CalcRMSDRotationalMatrix(
                self._ref_coordinates64, self._mobile_coordinates64,
                self._n_atoms, self._rot, self.weights_select)

            self._R[:, :] = self._rot.reshape(3, 3)
            # Transform each atom in the trajectory (use inplace ops to
            # avoid copying arrays) (Marginally (~3%) faster than
            # "ts.positions[:] = (ts.positions - x_com) * R + ref_com".)
            self._ts.positions[:] -= mobile_com

            # R acts to the left & is broadcasted N times.
            self._ts.positions[:] = np.dot(self._ts.positions, self._R)
            self._ts.positions += self._ref_com

            # 2) calculate secondary RMSDs (without any further
            #    superposition)
            for igroup, (refpos, atoms) in enumerate(
                    zip(self._groupselections_ref_coords64,
                        self._groupselections_atoms), 3):
                self.results.rmsd[self._frame_index, igroup] = rmsd(
                    refpos, atoms['mobile'].positions,
                    weights=self.weights_groupselections[igroup-3],
                    center=False, superposition=False)
        else:
            # only calculate RMSD by setting the Rmatrix to None (no need
            # to carry out the rotation as we already get the optimum RMSD)
            self.results.rmsd[self._frame_index, 2] = qcp.CalcRMSDRotationalMatrix(
                self._ref_coordinates64, self._mobile_coordinates64,
                self._n_atoms, None, self.weights_select)

    @property
    def rmsd(self):
        wmsg = ("The `rmsd` attribute was deprecated in MDAnalysis 2.0.0 and "
                "will be removed in MDAnalysis 3.0.0. Please use "
                "`results.rmsd` instead.")
        warnings.warn(wmsg, DeprecationWarning)
        return self.results.rmsd


class RMSF(AnalysisBase):
    r"""Calculate RMSF of given atoms across a trajectory.

    Note
    ----
    No RMSD-superposition is performed; it is assumed that the user is
    providing a trajectory where the protein of interest has been structurally
    aligned to a reference structure (see the Examples section below). The
    protein also has be whole because periodic boundaries are not taken into
    account.


    Run the analysis with :meth:`RMSF.run`, which stores the results
    in the array :attr:`RMSF.results.rmsf`.

    """
    def __init__(self, atomgroup, **kwargs):
        r"""Parameters
        ----------
        atomgroup : AtomGroup
            Atoms for which RMSF is calculated
        verbose : bool (optional)
             Show detailed progress of the calculation if set to ``True``; the
             default is ``False``.

        Raises
        ------
        ValueError
             raised if negative values are calculated, which indicates that a
             numerical overflow or underflow occured


        Notes
        -----
        The root mean square fluctuation of an atom :math:`i` is computed as the
        time average

        .. math::

          \rho_i = \sqrt{\left\langle (\mathbf{x}_i - \langle\mathbf{x}_i\rangle)^2 \right\rangle}

        No mass weighting is performed.

        This method implements an algorithm for computing sums of squares while
        avoiding overflows and underflows :footcite:p:`Welford1962`.


        Examples
        --------

        In this example we calculate the residue RMSF fluctuations by analyzing
        the :math:`\text{C}_\alpha` atoms. First we need to fit the trajectory
        to the average structure as a reference. That requires calculating the
        average structure first. Because we need to analyze and manipulate the
        same trajectory multiple times, we are going to load it into memory
        using the :mod:`~MDAnalysis.coordinates.MemoryReader`. (If your
        trajectory does not fit into memory, you will need to :ref:`write out
        intermediate trajectories <writing-trajectories>` to disk or
        :ref:`generate an in-memory universe
        <creating-in-memory-trajectory-label>` that only contains, say, the
        protein)::

           import MDAnalysis as mda
           from MDAnalysis.analysis import align

           from MDAnalysis.tests.datafiles import TPR, XTC

           u = mda.Universe(TPR, XTC, in_memory=True)
           protein = u.select_atoms("protein")

           # 1) the current trajectory contains a protein split across
           #    periodic boundaries, so we first make the protein whole and
           #    center it in the box using on-the-fly transformations
           import MDAnalysis.transformations as trans

           not_protein = u.select_atoms('not protein')
           transforms = [trans.unwrap(protein),
                         trans.center_in_box(protein, wrap=True),
                         trans.wrap(not_protein)]
           u.trajectory.add_transformations(*transforms)

           # 2) fit to the initial frame to get a better average structure
           #    (the trajectory is changed in memory)
           prealigner = align.AlignTraj(u, u, select="protein and name CA",
                                        in_memory=True).run()

           # 3) reference = average structure
           ref_coordinates = u.trajectory.timeseries(asel=protein).mean(axis=1)
           # make a reference structure (need to reshape into a 1-frame
           # "trajectory")
           reference = mda.Merge(protein).load_new(ref_coordinates[:, None, :],
                                                   order="afc")

        We created a new universe ``reference`` that contains a single frame
        with the averaged coordinates of the protein.  Now we need to fit the
        whole trajectory to the reference by minimizing the RMSD. We use
        :class:`MDAnalysis.analysis.align.AlignTraj`::

           aligner = align.AlignTraj(u, reference,
                                     select="protein and name CA",
                                     in_memory=True).run()

        The trajectory is now fitted to the reference (the RMSD is stored as
        `aligner.results.rmsd` for further inspection). Now we can calculate
        the RMSF::

           from MDAnalysis.analysis.rms import RMSF

           calphas = protein.select_atoms("name CA")
           rmsfer = RMSF(calphas, verbose=True).run()

        and plot::

           import matplotlib.pyplot as plt

           plt.plot(calphas.resnums, rmsfer.results.rmsf)



        References
        ----------
        .. footbibliography::


        .. versionadded:: 0.11.0
        .. versionchanged:: 0.16.0
           refactored to fit with AnalysisBase API
        .. deprecated:: 0.16.0
           the keyword argument `quiet` is deprecated in favor of `verbose`.
        .. versionchanged:: 0.17.0
           removed unused keyword `weights`
        .. versionchanged:: 1.0.0
           Support for the ``start``, ``stop``, and ``step`` keywords has been
           removed. These should instead be passed to :meth:`RMSF.run`.

        """
        super(RMSF, self).__init__(atomgroup.universe.trajectory, **kwargs)
        self.atomgroup = atomgroup

    def _prepare(self):
        self.sumsquares = np.zeros((self.atomgroup.n_atoms, 3))
        self.mean = self.sumsquares.copy()

    def _single_frame(self):
        k = self._frame_index
        self.sumsquares += (k / (k+1.0)) * (self.atomgroup.positions - self.mean) ** 2
        self.mean = (k * self.mean + self.atomgroup.positions) / (k + 1)

    def _conclude(self):
        k = self._frame_index
        self.results.rmsf = np.sqrt(self.sumsquares.sum(axis=1) / (k + 1))

        if not (self.results.rmsf >= 0).all():
            raise ValueError("Some RMSF values negative; overflow " +
                             "or underflow occurred")

    @property
    def rmsf(self):
        wmsg = ("The `rmsf` attribute was deprecated in MDAnalysis 2.0.0 and "
                "will be removed in MDAnalysis 3.0.0. Please use "
                "`results.rmsd` instead.")
        warnings.warn(wmsg, DeprecationWarning)
        return self.results.rmsf
