# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
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
Native contacts analysis --- :mod:`MDAnalysis.analysis.contacts`
================================================================

This module contains classes to analyze native contacts *Q* over a
trajectory. Native contacts of a conformation are contacts that exist
in a reference structure and in the conformation. Contacts in the
reference structure are always defined as being closer than a distance
`radius`. The fraction of native contacts for a conformation can be
calculated in different ways. This module supports 3 different metrics
listed below, as well as custom metrics.

1. *Hard Cut*: To count as a contact the atoms *i* and *j* have to be at least
   as close as in the reference structure.

2. *Soft Cut*: The atom pair *i* and *j* is assigned based on a soft potential
   that is 1 if the distance is 0, 1/2 if the distance is the same as in
   the reference and 0 for large distances. For the exact definition of the
   potential and parameters have a look at function :func:`soft_cut_q`.

3. *Radius Cut*: To count as a contact the atoms *i* and *j* cannot be further
   apart than some distance `radius`.

The "fraction of native contacts" *Q(t)* is a number between 0 and 1 and
calculated as the total number of native contacts for a given time frame
divided by the total number of contacts in the reference structure.


Examples for contact analysis
-----------------------------

One-dimensional contact analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As an example we analyze the opening ("unzipping") of salt bridges
when the AdK enzyme opens up; this is one of the example trajectories
in MDAnalysis. ::

    import numpy as np
    import matplotlib.pyplot as plt
    import MDAnalysis as mda
    from MDAnalysis.analysis import contacts
    from MDAnalysis.tests.datafiles import PSF,DCD
    # example trajectory (transition of AdK from closed to open)
    u = mda.Universe(PSF,DCD)
    # crude definition of salt bridges as contacts between NH/NZ in ARG/LYS and
    # OE*/OD* in ASP/GLU. You might want to think a little bit harder about the
    # problem before using this for real work.
    sel_basic = "(resname ARG LYS) and (name NH* NZ)"
    sel_acidic = "(resname ASP GLU) and (name OE* OD*)"
    # reference groups (first frame of the trajectory, but you could also use a
    # separate PDB, eg crystal structure)
    acidic = u.select_atoms(sel_acidic)
    basic = u.select_atoms(sel_basic)
    # set up analysis of native contacts ("salt bridges"); salt bridges have a
    # distance <6 A
    ca1 = contacts.Contacts(u, select=(sel_acidic, sel_basic),
                            refgroup=(acidic, basic), radius=6.0)
    # iterate through trajectory and perform analysis of "native contacts" Q
    ca1.run()
    # print number of averave contacts
    average_contacts = np.mean(ca1.results.timeseries[:, 1])
    print('average contacts = {}'.format(average_contacts))
    # plot time series q(t)
    fig, ax = plt.subplots()
    ax.plot(ca1.results.timeseries[:, 0], ca1.results.timeseries[:, 1])
    ax.set(xlabel='frame', ylabel='fraction of native contacts',
           title='Native Contacts, average = {:.2f}'.format(average_contacts))
    fig.show()


The first graph shows that when AdK opens, about 20% of the salt
bridges that existed in the closed state disappear when the enzyme
opens. They open in a step-wise fashion (made more clear by the movie
`AdK_zipper_cartoon.avi`_).

.. _`AdK_zipper_cartoon.avi`:
   http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2803350/bin/NIHMS150766-supplement-03.avi

.. rubric:: Notes

Suggested cutoff distances for different simulations

* For all-atom simulations, cutoff = 4.5 Å
* For coarse-grained simulations, cutoff = 6.0 Å


Two-dimensional contact analysis (q1-q2)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Analyze a single DIMS transition of AdK between its closed and open
conformation and plot the trajectory projected on q1-q2
:footcite:p:`Franklin2007` ::


    import MDAnalysis as mda
    from MDAnalysis.analysis import contacts
    from MDAnalysisTests.datafiles import PSF, DCD
    u = mda.Universe(PSF, DCD)
    q1q2 = contacts.q1q2(u, 'name CA', radius=8)
    q1q2.run()

    f, ax = plt.subplots(1, 2, figsize=plt.figaspect(0.5))
    ax[0].plot(q1q2.results.timeseries[:, 0], q1q2.results.timeseries[:, 1],
               label='q1')
    ax[0].plot(q1q2.results.timeseries[:, 0], q1q2.results.timeseries[:, 2],
               label='q2')
    ax[0].legend(loc='best')
    ax[1].plot(q1q2.results.timeseries[:, 1],
               q1q2.results.timeseries[:, 2], '.-')
    f.show()

Compare the resulting pathway to the `MinActionPath result for AdK`_
:footcite:p:`Franklin2007`.

.. _MinActionPath result for AdK:
   http://lorentz.dynstr.pasteur.fr/joel/adenylate.php


Writing your own contact analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :class:`Contacts` class has been designed to be extensible for your own
analysis. As an example we will analyze when the acidic and basic groups of AdK
are in contact which each other; this means that at least one of the contacts
formed in the reference is closer than 2.5 Å.

For this we define a new function to determine if any contact is closer than
2.5 Å; this function must implement the API prescribed by :class:`Contacts`::

    def is_any_closer(r, r0, dist=2.5):
        return np.any(r < dist)

The first two parameters `r` and `r0` are provided by :class:`Contacts` when it
calls :func:`is_any_closer` while the others can be passed as keyword args
using the `kwargs` parameter in :class:`Contacts`.

Next we are creating an instance of the :class:`Contacts` class and use the
:func:`is_any_closer` function as an argument to `method` and run the analysis::

    # crude definition of salt bridges as contacts between NH/NZ in ARG/LYS and
    # OE*/OD* in ASP/GLU. You might want to think a little bit harder about the
    # problem before using this for real work.
    sel_basic = "(resname ARG LYS) and (name NH* NZ)"
    sel_acidic = "(resname ASP GLU) and (name OE* OD*)"

    # reference groups (first frame of the trajectory, but you could also use a
    # separate PDB, eg crystal structure)
    acidic = u.select_atoms(sel_acidic)
    basic = u.select_atoms(sel_basic)

    nc = contacts.Contacts(u, select=(sel_acidic, sel_basic),
                           method=is_any_closer,
                           refgroup=(acidic, basic), kwargs={'dist': 2.5})
    nc.run()

    bound = nc.results.timeseries[:, 1]
    frames = nc.results.timeseries[:, 0]

    f, ax = plt.subplots()

    ax.plot(frames, bound, '.')
    ax.set(xlabel='frame', ylabel='is Bound',
           ylim=(-0.1, 1.1))

    f.show()


Functions
---------

.. autofunction:: hard_cut_q
.. autofunction:: soft_cut_q
.. autofunction:: radius_cut_q
.. autofunction:: contact_matrix
.. autofunction:: q1q2

Classes
-------

.. autoclass:: Contacts
   :members:

.. rubric:: References
.. footbibliography::

"""
import os
import errno
import warnings
import bz2
import functools

import numpy as np

import logging

import MDAnalysis
import MDAnalysis.lib.distances
from MDAnalysis.lib.util import openany
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.core.groups import AtomGroup, UpdatingAtomGroup
from .base import AnalysisBase, ResultsGroup

logger = logging.getLogger("MDAnalysis.analysis.contacts")


def soft_cut_q(r, r0, beta=5.0, lambda_constant=1.8):
    r"""Calculate fraction of native contacts *Q* for a soft cut off

    The native contact function is defined as :footcite:p:`Best2013`

    .. math::

        Q(r, r_0) = \frac{1}{1 + e^{\beta (r - \lambda r_0)}}

    Reasonable values for different simulation types are

    - *All Atom*: `lambda_constant = 1.8` (unitless)
    - *Coarse Grained*: `lambda_constant = 1.5` (unitless)

    Parameters
    ----------
    r: array
      Contact distances at time t
    r0: array
      Contact distances at time t=0, reference distances
    beta: float (default 5.0 Angstrom)
      Softness of the switching function
    lambda_constant: float (default 1.8, unitless)
      Reference distance tolerance

    Returns
    -------
    Q : float
      fraction of native contacts
    """
    r = np.asarray(r)
    r0 = np.asarray(r0)
    result = 1 / (1 + np.exp(beta * (r - lambda_constant * r0)))

    return result.sum() / len(r0)


def hard_cut_q(r, cutoff):
    """Calculate fraction of native contacts *Q* for a hard cut off.

    The cutoff can either be a float or a :class:`~numpy.ndarray` of the same
    shape as `r`.

    Parameters
    ----------
    r : ndarray
        distance matrix
    cutoff : ndarray | float
        cut off value to count distances. Can either be a float of a ndarray of
        the same size as distances

    Returns
    -------
    Q : float
        fraction of contacts

    """
    r = np.asarray(r)
    cutoff = np.asarray(cutoff)
    y = r <= cutoff
    return y.sum() / r.size


def radius_cut_q(r, r0, radius):
    """calculate native contacts *Q* based on the single distance radius.

    Parameters
    ----------
    r : ndarray
        distance array between atoms
    r0 : ndarray
        unused to fullfill :class:`Contacts` API
    radius : float
        Distance between atoms at which a contact is formed

    Returns
    -------
    Q : float
        fraction of contacts

    """
    return hard_cut_q(r, radius)


def contact_matrix(d, radius, out=None):
    """calculate contacts from distance matrix

    Parameters
    ----------
    d : array-like
        distance matrix
    radius : float
        distance below which a contact is formed.
    out : array (optional)
        If `out` is supplied as a pre-allocated array of the correct
        shape then it is filled instead of allocating a new one in
        order to increase performance.

    Returns
    -------
    contacts : ndarray
        boolean array of formed contacts
    """
    if out is not None:
        out[:] = d <= radius
    else:
        out = d <= radius
    return out


class Contacts(AnalysisBase):
    """Calculate contacts based observables.

    The standard methods used in this class calculate the fraction of native
    contacts *Q* from a trajectory.


    .. rubric:: Contact API

    By defining your own method it is possible to calculate other observables
    that only depend on the distances and a possible reference distance. The
    **Contact API** prescribes that this method must be a function with call
    signature ``func(r, r0, **kwargs)`` and must be provided in the keyword
    argument `method`.

    Attributes
    ----------
    results.timeseries : numpy.ndarray
        2D array containing *Q* for all refgroup pairs and analyzed frames

    timeseries : numpy.ndarray
        Alias to the :attr:`results.timeseries` attribute.

        .. deprecated:: 2.0.0
           Will be removed in MDAnalysis 3.0.0. Please use
           :attr:`results.timeseries` instead.


    .. versionchanged:: 1.0.0
       ``save()`` method has been removed. Use ``np.savetxt()`` on
       :attr:`Contacts.results.timeseries` instead.
    .. versionchanged:: 1.0.0
        added ``pbc`` attribute to calculate distances using PBC.
    .. versionchanged:: 2.0.0
       :attr:`timeseries` results are now stored in a
       :class:`MDAnalysis.analysis.base.Results` instance.
    .. versionchanged:: 2.2.0
       :class:`Contacts` accepts both AtomGroup and string for `select`
    .. versionchanged:: 2.9.0
       Introduced :meth:`get_supported_backends` allowing
       for parallel execution on :mod:`multiprocessing`
       and :mod:`dask` backends.
    """

    _analysis_algorithm_is_parallelizable = True

    @classmethod
    def get_supported_backends(cls):
        return (
            "serial",
            "multiprocessing",
            "dask",
        )

    def __init__(
        self,
        u,
        select,
        refgroup,
        method="hard_cut",
        radius=4.5,
        pbc=True,
        kwargs=None,
        **basekwargs,
    ):
        """
        Parameters
        ----------
        u : Universe
            trajectory
        select : tuple(AtomGroup, AtomGroup) | tuple(string, string)
            two contacting groups that change over time
        refgroup : tuple(AtomGroup, AtomGroup)
            two contacting atomgroups in their reference conformation. This
            can also be a list of tuples containing different atom groups
        radius : float, optional (4.5 Angstroms)
            radius within which contacts exist in refgroup
        method : string | callable (optional)
            Can either be one of ``['hard_cut' , 'soft_cut', 'radius_cut']`` or a callable
            with call signature ``func(r, r0, **kwargs)`` (the "Contacts API").
        pbc : bool (optional)
            Uses periodic boundary conditions to calculate distances if set to ``True``; the
            default is ``True``.
        kwargs : dict, optional
            dictionary of additional kwargs passed to `method`. Check
            respective functions for reasonable values.
        verbose : bool (optional)
             Show detailed progress of the calculation if set to ``True``; the
             default is ``False``.

        Attributes
        ----------
        n_initial_contacts : int
             Total number of initial contacts.
        r0 : list[numpy.ndarray]
             List of distance arrays between reference groups.

        Notes
        -----

        .. versionchanged:: 1.0.0
           Changed `selection` keyword to `select`
        """
        self.u = u
        super(Contacts, self).__init__(self.u.trajectory, **basekwargs)

        self.fraction_kwargs = kwargs if kwargs is not None else {}

        if method == "hard_cut":
            self.fraction_contacts = hard_cut_q
        elif method == "soft_cut":
            self.fraction_contacts = soft_cut_q
        elif method == "radius_cut":
            self.fraction_contacts = functools.partial(
                radius_cut_q, radius=radius
            )
        else:
            if not callable(method):
                raise ValueError("method has to be callable")
            self.fraction_contacts = method

        self.select = select

        self.grA, self.grB = (self._get_atomgroup(u, sel) for sel in select)

        self.pbc = pbc

        # contacts formed in reference
        self.r0 = []
        self.initial_contacts = []

        # get dimensions via partial for parallelization compatibility
        self._get_box = functools.partial(self._get_box_func, pbc=self.pbc)

        if isinstance(refgroup[0], AtomGroup):
            refA, refB = refgroup
            self.r0.append(
                distance_array(
                    refA.positions,
                    refB.positions,
                    box=self._get_box(refA.universe),
                )
            )
            self.initial_contacts.append(contact_matrix(self.r0[-1], radius))

        else:
            for refA, refB in refgroup:
                self.r0.append(
                    distance_array(
                        refA.positions,
                        refB.positions,
                        box=self._get_box(refA.universe),
                    )
                )
                self.initial_contacts.append(
                    contact_matrix(self.r0[-1], radius)
                )

        self.n_initial_contacts = self.initial_contacts[0].sum()

    @staticmethod
    def _get_atomgroup(u, sel):
        select_error_message = (
            "selection must be either string or a "
            "static AtomGroup. Updating AtomGroups "
            "are not supported."
        )
        if isinstance(sel, str):
            return u.select_atoms(sel)
        elif isinstance(sel, AtomGroup):
            if isinstance(sel, UpdatingAtomGroup):
                raise TypeError(select_error_message)
            else:
                return sel
        else:
            raise TypeError(select_error_message)

    @staticmethod
    def _get_box_func(ts, pbc):
        """Retrieve the dimensions of the simulation box based on PBC.

        Parameters
        ----------
        ts : Timestep
            The current timestep of the simulation, which contains the
            box dimensions.
        pbc : bool
            A flag indicating whether periodic boundary conditions (PBC)
            are enabled. If `True`, the box dimensions are returned,
            else returns `None`.

        Returns
        -------
        box_dimensions : ndarray or None
            The dimensions of the simulation box as a NumPy array if PBC
            is True, else returns `None`.
        """
        return ts.dimensions if pbc else None

    def _prepare(self):
        self.results.timeseries = np.empty((self.n_frames, len(self.r0) + 1))

    def _single_frame(self):
        self.results.timeseries[self._frame_index][0] = self._ts.frame

        # compute distance array for a frame
        d = distance_array(
            self.grA.positions, self.grB.positions, box=self._get_box(self._ts)
        )

        for i, (initial_contacts, r0) in enumerate(
            zip(self.initial_contacts, self.r0), 1
        ):
            # select only the contacts that were formed in the reference state
            r = d[initial_contacts]
            r0 = r0[initial_contacts]
            q = self.fraction_contacts(r, r0, **self.fraction_kwargs)
            self.results.timeseries[self._frame_index][i] = q

    @property
    def timeseries(self):
        wmsg = (
            "The `timeseries` attribute was deprecated in MDAnalysis "
            "2.0.0 and will be removed in MDAnalysis 3.0.0. Please use "
            "`results.timeseries` instead"
        )
        warnings.warn(wmsg, DeprecationWarning)
        return self.results.timeseries

    def _get_aggregator(self):
        return ResultsGroup(lookup={"timeseries": ResultsGroup.ndarray_vstack})


def _new_selections(u_orig, selections, frame):
    """create stand alone AGs from selections at frame"""
    u = MDAnalysis.Universe(u_orig.filename, u_orig.trajectory.filename)
    u.trajectory[frame]
    return [u.select_atoms(s) for s in selections]


def q1q2(u, select="all", radius=4.5):
    """Perform a q1-q2 analysis.

    Compares native contacts between the starting structure and final structure
    of a trajectory :footcite:p:`Franklin2007`.

    Parameters
    ----------
    u : Universe
        Universe with a trajectory
    select : string, optional
        atoms to do analysis on
    radius : float, optional
        distance at which contact is formed

    Returns
    -------
    contacts : :class:`Contacts`
        Contact Analysis that is set up for a q1-q2 analysis


    .. versionchanged:: 1.0.0
       Changed `selection` keyword to `select`
       Support for setting ``start``, ``stop``, and ``step`` has been removed.
       These should now be directly passed to :meth:`Contacts.run`.
    """
    selection = (select, select)
    first_frame_refs = _new_selections(u, selection, 0)
    last_frame_refs = _new_selections(u, selection, -1)
    return Contacts(
        u,
        selection,
        (first_frame_refs, last_frame_refs),
        radius=radius,
        method="radius_cut",
    )
