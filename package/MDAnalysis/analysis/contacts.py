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

"""
Native contacts analysis --- :mod:`MDAnalysis.analysis.contacts`
================================================================

This module contains classes to analyze native contacts *Q* over a
trajectory. Native contacts of a conformation are contacts that exist
in a reference structure and in the conformation. Contacts in the
reference structure are always defined as being closer then a distance
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

    import MDAnalysis as mda
    from MDAnalysis.analysis import contacts
    from MDAnalysis.tests.datafiles import PSF,DCD
    import matplotlib.pyplot as plt
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
    ca1 = contacts.Contacts(u, selection=(sel_acidic, sel_basic),
                            refgroup=(acidic, basic), radius=6.0)
    # iterate through trajectory and perform analysis of "native contacts" Q
    ca1.run()
    # print number of averave contacts
    average_contacts = np.mean(ca1.timeseries[:, 1])
    print('average contacts = {}'.format(average_contacts))
    # plot time series q(t)
    f, ax = plt.subplots()
    ax.plot(ca1.timeseries[:, 0], ca1.timeseries[:, 1])
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
conformation and plot the trajectory projected on q1-q2 [Franklin2007]_ ::


    import MDAnalysis as mda
    from MDAnalysis.analysis import contacts
    from MDAnalysisTests.datafiles import PSF, DCD
    u = mda.Universe(PSF, DCD)
    q1q2 = contacts.q1q2(u, 'name CA', radius=8)
    q1q2.run()

    f, ax = plt.subplots(1, 2, figsize=plt.figaspect(0.5))
    ax[0].plot(q1q2.timeseries[:, 0], q1q2.timeseries[:, 1], label='q1')
    ax[0].plot(q1q2.timeseries[:, 0], q1q2.timeseries[:, 2], label='q2')
    ax[0].legend(loc='best')
    ax[1].plot(q1q2.timeseries[:, 1], q1q2.timeseries[:, 2], '.-')
    f.show()

Compare the resulting pathway to the `MinActionPath result for AdK`_
[Franklin2007]_.

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

    nc = contacts.Contacts(u, selection=(sel_acidic, sel_basic),
                           method=is_any_closer,
                           refgroup=(acidic, basic), kwargs={'dist': 2.5})
    nc.run()

    bound = nc.timeseries[:, 1]
    frames = nc.timeseries[:, 0]

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

"""
from __future__ import division, absolute_import
from six.moves import zip

import os
import errno
import warnings
import bz2

import numpy as np

import logging

import MDAnalysis
import MDAnalysis.lib.distances
from MDAnalysis.lib.util import openany, deprecate
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.core.groups import AtomGroup
from .base import AnalysisBase

logger = logging.getLogger("MDAnalysis.analysis.contacts")


def soft_cut_q(r, r0, beta=5.0, lambda_constant=1.8):
    r"""Calculate fraction of native contacts *Q* for a soft cut off

    The native contact function is defined as [Best2013]_

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

    References
    ----------
    .. [Best2013] RB Best, G Hummer, and WA Eaton, "Native contacts determine protein
       folding mechanisms in atomistic simulations" _PNAS_ **110** (2013),
       17874–17879. doi: `10.1073/pnas.1311599110
       <http://doi.org/10.1073/pnas.1311599110>`_.

    """
    r = np.asarray(r)
    r0 = np.asarray(r0)
    result = 1/(1 + np.exp(beta*(r - lambda_constant * r0)))

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

    References
    ----------
    .. [Franklin2007] Franklin, J., Koehl, P., Doniach, S., & Delarue,
       M. (2007).  MinActionPath: Maximum likelihood trajectory for large-scale
       structural transitions in a coarse-grained locally harmonic energy
       landscape.  Nucleic Acids Research, 35(SUPPL.2), 477–482.
       doi: `10.1093/nar/gkm342 <http://doi.org/10.1093/nar/gkm342>`_

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
    out: array (optional)
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
    timeseries : list
        list containing *Q* for all refgroup pairs and analyzed frames

    """
    def __init__(self, u, selection, refgroup, method="hard_cut", radius=4.5,
                 kwargs=None, **basekwargs):
        """
        Parameters
        ----------
        u : Universe
            trajectory
        selection : tuple(string, string)
            two contacting groups that change over time
        refgroup : tuple(AtomGroup, AtomGroup)
            two contacting atomgroups in their reference conformation. This
            can also be a list of tuples containing different atom groups
        radius : float, optional (4.5 Angstroms)
            radius within which contacts exist in refgroup
        method : string | callable (optional)
            Can either be one of ``['hard_cut' , 'soft_cut']`` or a callable
            with call signature ``func(r, r0, **kwargs)`` (the "Contacts API").
        kwargs : dict, optional
            dictionary of additional kwargs passed to `method`. Check
            respective functions for reasonable values.
        start : int, optional
            First frame of trajectory to analyse, Default: None becomes 0.
        stop : int, optional
            Frame index to stop analysis. Default: None becomes
            n_frames. Iteration stops *before* this frame number,
            which means that the trajectory would be read until the end.
        step : int, optional
            Step between frames to analyse, Default: None becomes 1.
        verbose : bool (optional)
             Show detailed progress of the calculation if set to ``True``; the
             default is ``False``.
        """
        self.u = u
        super(Contacts, self).__init__(self.u.trajectory, **basekwargs)

        if method == 'hard_cut':
            self.fraction_contacts = hard_cut_q
        elif method == 'soft_cut':
            self.fraction_contacts = soft_cut_q
        else:
            if not callable(method):
                raise ValueError("method has to be callable")
            self.fraction_contacts = method

        self.selection = selection
        self.grA = u.select_atoms(selection[0])
        self.grB = u.select_atoms(selection[1])

        # contacts formed in reference
        self.r0 = []
        self.initial_contacts = []

        if isinstance(refgroup[0], AtomGroup):
            refA, refB = refgroup
            self.r0.append(distance_array(refA.positions, refB.positions))
            self.initial_contacts.append(contact_matrix(self.r0[-1], radius))
        else:
            for refA, refB in refgroup:
                self.r0.append(distance_array(refA.positions, refB.positions))
                self.initial_contacts.append(contact_matrix(self.r0[-1],
                                                            radius))

        self.fraction_kwargs = kwargs if kwargs is not None else {}
        self.timeseries = []

    def _single_frame(self):
        # compute distance array for a frame
        d = distance_array(self.grA.positions, self.grB.positions)

        y = np.empty(len(self.r0) + 1)
        y[0] = self._ts.frame
        for i, (initial_contacts, r0) in enumerate(zip(self.initial_contacts,
                                                       self.r0)):
            # select only the contacts that were formed in the reference state
            r = d[initial_contacts]
            r0 = r0[initial_contacts]
            y[i + 1] = self.fraction_contacts(r, r0, **self.fraction_kwargs)

        if len(y) == 1:
            y = y[0]
        self.timeseries.append(y)

    def _conclude(self):
        self.timeseries = np.array(self.timeseries, dtype=float)

    @deprecate(release="0.19.0", remove="1.0.0")
    def save(self, outfile):
        """save contacts timeseries

        Parameters
        ----------
        outfile : str
            file to save contacts

        """
        np.savetxt(outfile, self.timeseries,
                   header="# q1 analysis\n", comments='')


def _new_selections(u_orig, selections, frame):
    """create stand alone AGs from selections at frame"""
    u = MDAnalysis.Universe(u_orig.filename, u_orig.trajectory.filename)
    u.trajectory[frame]
    return [u.select_atoms(s) for s in selections]


def q1q2(u, selection='all', radius=4.5,
         start=None, stop=None, step=None):
    """Perform a q1-q2 analysis.

    Compares native contacts between the starting structure and final structure
    of a trajectory [Franklin2007]_.

    Parameters
    ----------
    u : Universe
        Universe with a trajectory
    selection : string, optional
        atoms to do analysis on
    radius : float, optional
        distance at which contact is formed
    start : int, optional
        First frame of trajectory to analyse, Default: 0
    stop : int, optional
        Last frame of trajectory to analyse, Default: -1
    step : int, optional
        Step between frames to analyse, Default: 1

    Returns
    -------
    contacts : :class:`Contacts`
        Contact Analysis that is set up for a q1-q2 analysis

    """
    selection = (selection, selection)
    first_frame_refs = _new_selections(u, selection, 0)
    last_frame_refs = _new_selections(u, selection, -1)
    return Contacts(u, selection,
                    (first_frame_refs, last_frame_refs),
                    radius=radius, method=radius_cut_q,
                    start=start, stop=stop, step=step,
                    kwargs={'radius': radius})
