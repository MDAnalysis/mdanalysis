# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
# Beckstein and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""Native contacts analysis --- :mod:`MDAnalysis.analysis.contacts`
================================================================


Analysis of native contacts *Q* over a trajectory. Native contacts of a
conformation are contacts that exist in a reference structure and in the
conformation. Contacts in the reference structure are always defined as being
closer then a distance `radius`. The fraction of native contacts for a
conformation can be calculated in different ways. This module supports 3
different metrics liseted below, as wel as custom metrics.

1. *Hard Cut*: To count as a contact the atoms *i* and *j* have to be at least
   as close as in the reference structure.

2. *Soft Cut*: The atom pair *i* and *j* is assigned based on a soft potential
   that is 1 for if the distance is 0, 1./2 if the distance is the same as in
   the reference and 0 for large distances. For the exact definition of the
   potential and parameters have a look at `soft_cut_q`.

3. *Radius Cut*: To count as a contact the atoms *i* and *j* cannot be further
   apart then some distance `radius`.

The "fraction of native contacts" *Q(t)* is a number between 0 and 1 and
calculated as the total number of native contacts for a given time frame
divided by the total number of contacts in the reference structure.

Examples
--------

One-dimensional contact analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As an example we analyze the opening ("unzipping") of salt bridges
when the AdK enzyme opens up; this is one of the example trajectories
in MDAnalysis. ::

>>> import MDAnalysis as mda
>>> from MDAnalysis.analysis import contacts
>>> from MDAnalysis.tests.datafiles import PSF,DCD
>>> import matplotlib.pyplot as plt
>>> # example trajectory (transition of AdK from closed to open)
>>> u = mda.Universe(PSF,DCD)
>>> # crude definition of salt bridges as contacts between NH/NZ in ARG/LYS and
>>> # OE*/OD* in ASP/GLU. You might want to think a little bit harder about the
>>> # problem before using this for real work.
>>> sel_basic = "(resname ARG LYS) and (name NH* NZ)"
>>> sel_acidic = "(resname ASP GLU) and (name OE* OD*)"
>>> # reference groups (first frame of the trajectory, but you could also use a
>>> # separate PDB, eg crystal structure)
>>> acidic = u.select_atoms(sel_acidic)
>>> basic = u.select_atoms(sel_basic)
>>> # set up analysis of native contacts ("salt bridges"); salt bridges have a
>>> # distance <6 A
>>> ca1 = contacts.Contacts(u, selection=(sel_acidic, sel_basic),
>>>                         refgroup=(acidic, basic), radius=6.0)
>>> # iterate through trajectory and perform analysis of "native contacts" Q
>>> ca1.run()
>>> # print number of averave contacts
>>> average_contacts = np.mean(ca1.timeseries[:, 1])
>>> print('average contacts = {}'.format(average_contacts))
>>> # plot time series q(t)
>>> f, ax = plt.subplots()
>>> ax.plot(ca1.timeseries[:, 0], ca1.timeseries[:, 1])
>>> ax.set(xlabel='frame', ylabel='fraction of native contacts',
           title='Native Contacts, average = {:.2f}'.format(average_contacts))
>>> fig.show()


The first graph shows that when AdK opens, about 20% of the salt
bridges that existed in the closed state disappear when the enzyme
opens. They open in a step-wise fashion (made more clear by the movie
`AdK_zipper_cartoon.avi`_).

.. AdK_zipper_cartoon.avi:
   http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2803350/bin/NIHMS150766-supplement-03.avi

Notes
-----
Suggested cutoff distances for different simulations
* For all-atom simulations, cutoff = 4.5 A
* For coarse-grained simulations, cutoff = 6.0 A

Two-dimensional contact analysis (q1-q2)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Analyze a single DIMS transition of AdK between its closed and open
conformation and plot the trajectory projected on q1-q2::


>>> import MDAnalysis as mda
>>> from MDAnalysis.analysis import contacts
>>> from MDAnalysisTests.datafiles import PSF, DCD
>>> u = mda.Universe(PSF, DCD)
>>> q1q2 = contacts.q1q2(u, 'name CA', radius=8)
>>> q1q2.run()
>>>
>>> f, ax = plt.subplots(1, 2, figsize=plt.figaspect(0.5))
>>> ax[0].plot(q1q2.timeseries[:, 0], q1q2.timeseries[:, 1], label='q1')
>>> ax[0].plot(q1q2.timeseries[:, 0], q1q2.timeseries[:, 2], label='q2')
>>> ax[0].legend(loc='best')
>>> ax[1].plot(q1q2.timeseries[:, 1], q1q2.timeseries[:, 2], '.-')
>>> f.show()

Compare the resulting pathway to the `MinActionPath result for AdK`_.

.. _MinActionPath result for AdK:
   http://lorentz.dynstr.pasteur.fr/joel/adenylate.php

Writing your own contact analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ..class:`Contacts` has been designed to be extensible for your own
analysis. As an example we will analysis when the acidic and basic groups of
are in contact which each other, this means that at least one of the contacts
formed in the reference is closer then 2.5 Angstrom. For this we define a new
method to determine if any contact is closer then 2.5 Angström that implements
the API ..class:`Contacts` except.

The first to parameters `r` and `r0` are provided by ..class:`Contacts` the
others can be passed as keyword args using the `kwargs` parameter in
..class:`Contacts`.

>>> def is_any_closer(r, r0, dist=2.5):
>>>     return np.any(r < dist)

Next we are creating an instance of the Constants class and use the
`is_any_closer` function as an argument to `method` and run the analysus

>>> # crude definition of salt bridges as contacts between NH/NZ in ARG/LYS and
>>> # OE*/OD* in ASP/GLU. You might want to think a little bit harder about the
>>> # problem before using this for real work.
>>> sel_basic = "(resname ARG LYS) and (name NH* NZ)"
>>> sel_acidic = "(resname ASP GLU) and (name OE* OD*)"
>>> # reference groups (first frame of the trajectory, but you could also use a
>>> # separate PDB, eg crystal structure)
>>> acidic = u.select_atoms(sel_acidic)
>>> basic = u.select_atoms(sel_basic)
>>> nc = contacts.Contacts(u, selection=(sel_acidic, sel_basic),
>>>                        method=is_any_closer,
>>>                        refgroup=(acidic, basic), kwargs={'dist': 2.5})
>>>
>>> nc.run()
>>>
>>> bound = nc.timeseries[:, 1]
>>> frames = nc.timeseries[:, 0]
>>>
>>> f, ax = plt.subplots()
>>>
>>> ax.plot(frames, bound, '.')
>>> ax.set(xlabel='frame', ylabel='is Bound',
>>>        ylim=(-0.1, 1.1))
>>>
>>> f.show()

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


Deprecated
----------

.. autoclass:: ContactAnalysis1
   :members:
.. autoclass:: ContactAnalysis
   :members:

"""
from __future__ import division

import os
import errno
import warnings
import bz2
from six.moves import zip

import numpy as np
from numpy.lib.utils import deprecate

import logging

import MDAnalysis
import MDAnalysis.lib.distances
from MDAnalysis.lib.util import openany
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.core.AtomGroup import AtomGroup
from .base import AnalysisBase

logger = logging.getLogger("MDAnalysis.analysis.contacts")


def soft_cut_q(r, r0, beta=5.0, lambda_constant=1.8):
    r"""Calculate fraction of native contacts *Q* for a soft cut off

    ..math::
        Q(r, r_0) = \frac{1}{1 + e^{\beta (r - \lambda r_0)}}

    Reasonable values for different simulation types are

    *All Atom*: lambda_constant = 1.8 (unitless)
    *Coarse Grained*: lambda_constant = 1.5 (unitless)

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
    .. [1] RB Best, G Hummer, and WA Eaton, "Native contacts determine protein
       folding mechanisms in atomistic simulations" _PNAS_ **110** (2013),
       17874–17879. `10.1073/pnas.1311599110
       <http://doi.org/10.1073/pnas.1311599110>`_.

    """
    r = np.asarray(r)
    r0 = np.asarray(r0)
    result = 1/(1 + np.exp(beta*(r - lambda_constant * r0)))

    return result.sum() / len(r0)


def hard_cut_q(r, cutoff):
    """Calculate fraction of native contacts *Q* for a hard cut off. The cutoff
    can either be a float a ndarray the same shape as `r`

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
    """calculate native contacts *Q* based on the single sistance radius

    Parameters
    ----------
    r : ndarray
        distance array between atoms
    r0 : ndarray
        unused to fullfill Contacts API
    radius : float
        Distance between atoms at which a contact is formed

    Returns
    -------
    Q : float
        fraction of contacts

    References
    ----------
    .. [1] Franklin, J., Koehl, P., Doniach, S., & Delarue, M. (2007).
       MinActionPath: Maximum likelihood trajectory for large-scale structural
       transitions in a coarse-grained locally harmonic energy landscape.
       Nucleic Acids Research, 35(SUPPL.2), 477–482.
       http://doi.org/10.1093/nar/gkm342

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
    contacts *Q* from a trajectory. By defining your own method it is possible
    to calculate other observables that only depend on the distances and a
    possible reference distance.

    Attributes
    ----------
    timeseries : list
        list containing *Q* for all refgroup pairs and analyzed frames

    """
    def __init__(self, u, selection, refgroup, method="hard_cut", radius=4.5,
                 kwargs=None, start=None, stop=None, step=None,):
        """Initialization

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
            Can either be one of ['hard_cut' , 'soft_cut'] or a callable that
            implements a API (r, r0, **kwargs).
        kwargs : dict, optional
            dictionary of additional kwargs passed to `method`. Check
            respective functions for reasonable values.
        start : int, optional
            First frame of trajectory to analyse, Default: 0
        stop : int, optional
            Last frame of trajectory to analyse, Default: -1
        step : int, optional
            Step between frames to analyse, Default: 1

        """
        if method == 'hard_cut':
            self.fraction_contacts = hard_cut_q
        elif method == 'soft_cut':
            self.fraction_contacts = soft_cut_q
        else:
            if not callable(method):
                raise ValueError("method has to be callable")
            self.fraction_contacts = method

        # setup boilerplate
        self.u = u
        self._setup_frames(self.u.trajectory, start, stop, step)

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

    def save(self, outfile):
        """save contacts timeseries

        Parameter
        ---------
        outfile : str
            file to save contacts

        """
        with open(outfile, "w") as f:
            f.write("# q1 analysis\n")
            np.savetxt(f, self.timeseries)


def _new_selections(u_orig, selections, frame):
    """create stand alone AGs from selections at frame"""
    u = MDAnalysis.Universe(u_orig.filename, u_orig.trajectory.filename)
    u.trajectory[frame]
    return [u.select_atoms(s) for s in selections]


def q1q2(u, selection='all', radius=4.5,
         start=None, stop=None, step=None):
    """Do a q1-q2 analysis to compare native contacts between the starting
    structure and final structure of a trajectory.

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
    contacts
        Contact Analysis setup for a q1-q2 analysis

    """
    selection = (selection, selection)
    first_frame_refs = _new_selections(u, selection, 0)
    last_frame_refs = _new_selections(u, selection, -1)
    return Contacts(u, selection,
                    (first_frame_refs, last_frame_refs),
                    radius=radius, method=radius_cut_q,
                    start=start, stop=stop, step=step,
                    kwargs={'radius': radius})

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################


# What comes now are old deprecated contact Analysis classes


# ContactAnalysis needs to be cleaned up and possibly renamed but
# until then it remains because we don't have the functionality
# elsewhere.

@deprecate(new_name="Contacts", message="This class will be removed in 0.17")
class ContactAnalysis(object):
    """Perform a native contact analysis ("q1-q2").

    The analysis of the trajectory is performed with the
    :meth:`ContactAnalysis.run` method. The result is stored in
    :attr:`ContactAnalysis.timeseries`. It is a numpy array which
    contains the frame number at index 0, q1 and q2 at index 1 and 2,
    and the total number of contacts in 3 and 4. ::

        frame  q1 q2  n1 n2

    The total number of contacts in the reference states 1 and 2 are
    stored in :attr:`ContactAnalysis.nref` (index 0 and 1).

    The :meth:`ContactAnalysis.run` method calculates the percentage of native
    contacts *q1* and *q2* along a trajectory. "Contacts" are defined as the
    number of Ca atoms (or per-residue *centroids* of a user defined
    *selection*) within *radius* of a primary Ca. *q1* is the fraction of
    contacts relative to the reference state 1 (typically the starting
    conformation of the trajectory) and *q2* is the fraction of contacts
    relative to the conformation 2.

    The timeseries is written to a bzip2-compressed file in `targetdir`
    named "basename(trajectory)infix_q1q2.dat.bz2" and is also
    accessible as the attribute
    :attr:`ContactAnalysis.timeseries`.

    """

    def __init__(self, topology, trajectory, ref1=None, ref2=None, radius=8.0,
                 targetdir=os.path.curdir, infix="", force=False,
                 selection="name CA", centroids=False):
        """Calculate native contacts from two reference structures.

        Parameters
        ----------
        topology : filename
            topology file
        trajectory : filename
            trajectory
        ref1 : filename or ``None``, optional
            structure of the reference conformation 1 (pdb); if ``None`` the
            *first* frame of the trajectory is chosen
        ref2 : filename or ``None``, optional
            structure of the reference conformation 2 (pdb); if ``None`` the
            *last* frame of the trajectory is chosen
        radius : float, optional, default 8 A
            contacts are deemed any Ca within radius
        targetdir : path, optional, default ``.``
            output files are saved in this directory
        infix : string, optional
            additional tag string that is inserted into the output filename of
            the data file
         selection : string, optional, default ``"name CA"``
            MDAnalysis selection string that selects the particles of
            interest; the default is to only select the C-alpha atoms
            in `ref1` and `ref2`

            .. Note:: If `selection` produces more than one atom per
                      residue then you will get multiple contacts per
                      residue unless you also set `centroids` = ``True``
         centroids : bool
            If set to ``True``, use the centroids for the selected atoms on a
            per-residue basis to compute contacts. This allows, for instance
            defining the sidechains as `selection` and then computing distances
            between sidechain centroids.

        """

        self.topology = topology
        self.trajectory = trajectory
        self.radius = radius
        self.targetdir = targetdir
        self.force = force
        self.selection = selection
        self.centroids = centroids

        trajectorybase = os.path.splitext(os.path.basename(trajectory))[0]
        output = trajectorybase + infix + '_q1q2.dat'
        self.output = os.path.join(self.targetdir, output)
        self.output_bz2 = self.output + '.bz2'

        self.timeseries = None  # final result

        # short circuit if output file already exists: skip everything
        if self.output_exists():
            self._skip = True
            # do not bother reading any data or initializing arrays... !!
            return

        # don't bother if trajectory is empty (can lead to segfaults so better
        # catch it)
        stats = os.stat(trajectory)
        if stats.st_size == 0:
            warnings.warn('trajectory = {trajectory!s} is empty, '
                          'skipping...'.format(**vars()))
            self._skip = True
            return
        # under normal circumstances we do not skip
        self._skip = False

        # expensive initialization starts with building Universes :-)
        self.u = MDAnalysis.Universe(topology, trajectory)

        if ref1 is None:
            ref1 = os.path.join(self.targetdir, trajectorybase + '_first.pdb')
            self.u.trajectory[0]  # extract first frame
            self.u.atoms.write(ref1)
        self.ref1 = ref1
        if ref2 is None:
            ref2 = os.path.join(self.targetdir, trajectorybase + '_last.pdb')
            self.u.trajectory[-1]  # extract last frame
            self.u.atoms.write(ref2)
            self.u.trajectory[0]  # rewind, just in case...
        self.ref2 = ref2

        r1 = MDAnalysis.Universe(topology, self.ref1)
        r2 = MDAnalysis.Universe(topology, self.ref2)

        self.ca = self.u.select_atoms(self.selection)
        ca1 = r1.select_atoms(self.selection)
        ca2 = r2.select_atoms(self.selection)

        # NOTE: self_distance_array() produces a 1D array; this works here
        #       but is not the same as the 2D output from distance_array()!
        #       See the docs for self_distance_array().
        dref = [self.get_distance_array(ca1), self.get_distance_array(ca2)]
        self.qref = [self.qarray(dref[0]), self.qarray(dref[1])]
        self.nref = [self.qref[0].sum(), self.qref[1].sum()]

        self.d = np.zeros_like(dref[0])
        self.q = self.qarray(self.d)
        self._qtmp = np.zeros_like(self.q)  # pre-allocated array

    def get_distance_array(self, g, **kwargs):
        """Calculate the self_distance_array for atoms in group *g*.

        Parameters
        ----------
        g : AtomGroup
              group of atoms to calculate distance array for
        results : array, optional
              passed on to :func:`MDAnalysis.lib.distances.self_distance_array`
              as a preallocated array
        centroids : bool, optional, default ``None``
              ``True``: calculate per-residue centroids from the selected
              atoms; ``False``: consider each atom separately; ``None``: use
              the class default for *centroids* [``None``]

        """
        centroids = kwargs.pop("centroids", None)
        centroids = self.centroids if centroids is None else centroids
        if not centroids:
            coordinates = g.positions
        else:
            # centroids per residue (but only including the selected atoms)
            coordinates = np.array([residue.centroid()
                                    for residue in g.split("residue")])
        return MDAnalysis.lib.distances.self_distance_array(coordinates,
                                                            **kwargs)

    def output_exists(self, force=False):
        """Return True if default output file already exists.

        Disable with force=True (will always return False)
        """
        return (os.path.isfile(self.output) or
                os.path.isfile(self.output_bz2)) and not (self.force or force)

    def run(self, store=True, force=False):
        """Analyze trajectory and produce timeseries.

        Stores results in :attr:`ContactAnalysis.timeseries` (if
        store=True) and writes them to a bzip2-compressed data file.
        """
        if self._skip or self.output_exists(force=force):
            warnings.warn("File {output!r} or {output_bz2!r} already exists, "
                          "loading {trajectory!r}.".format(**vars(self)))
            try:
                self.load(self.output)
            except IOError:
                self.load(self.output_bz2)
            return None

        outbz2 = bz2.BZ2File(self.output_bz2, mode='w', buffering=8192)
        try:
            outbz2.write("# q1-q2 analysis\n"
                         "# nref1 = {0:d}\n"
                         "# nref2 = {1:d}\n".format(self.nref[0],
                                                    self.nref[1]))
            outbz2.write("# frame  q1  q2   n1  n2\n")
            records = []
            for ts in self.u.trajectory:
                frame = ts.frame
                # use pre-allocated distance array to save a little bit of time
                self.get_distance_array(self.ca, result=self.d)
                self.qarray(self.d, out=self.q)
                n1, q1 = self.qN(self.q, 0, out=self._qtmp)
                n2, q2 = self.qN(self.q, 1, out=self._qtmp)

                if store:
                    records.append((frame, q1, q2, n1, n2))
                outbz2.write("{frame:4d}  {q1:8.6f} {q2:8.6f}  {n1:5d} {n2:5d}\n".format(**vars()))
        finally:
            outbz2.close()
        if store:
            self.timeseries = np.array(records).T
        return self.output_bz2

    def qarray(self, d, out=None):
        """Return array with ``True`` for contacts.

        Note
        ----
        This method is typically only used internally.

        Arguments
        ---------
        d : array
          2D array of distances. The method uses the value of
          :attr:`radius` to determine if a ``distance < radius``
          is considered a contact.
        out : array, optional
          If `out` is supplied as a pre-allocated array of the correct
          shape then it is filled instead of allocating a new one in
          order to increase performance.

        Returns
        -------
        array
           contact matrix

        """
        if out is None:
            out = (d <= self.radius)
        else:
            out[:] = (d <= self.radius)
        return out

    def qN(self, q, n, out=None):
        """Calculate native contacts relative to reference state.

        Note
        ----
        This method is typically only used internally.

        Arguments
        ---------
        q : array
          contact matrix (see :meth:`Contacts.qarray`)
        out : array, optional
          If `out` is supplied as a pre-allocated array of the correct
          shape then it will contain the contact matrix relative
          to the reference state, i.e. only those contacts that
          are also seen in the reference state.

        Returns
        -------
        contacts : integer
           total number of contacts
        fraction : float
           fraction of contacts relative to the reference state

        """

        if out is None:
            out = np.logical_and(q, self.qref[n])
        else:
            np.logical_and(q, self.qref[n], out)
        contacts = out.sum()
        return contacts, float(contacts) / self.nref[n]

    def load(self, filename):
        """Load the data file."""
        records = []
        with openany(filename) as data:
            for line in data:
                if line.startswith('#'):
                    continue
                records.append(map(float, line.split()))
        self.timeseries = np.array(records).T

    def plot(self, **kwargs):
        """Plot q1-q2."""
        from pylab import plot, xlabel, ylabel

        kwargs.setdefault('color', 'black')
        if self.timeseries is None:
            raise ValueError("No timeseries data; do "
                             "'ContactAnalysis.run(store=True)' first.")
        t = self.timeseries
        plot(t[1], t[2], **kwargs)
        xlabel(r"$q_1$")
        ylabel(r"$q_2$")


@deprecate(new_name="Contacts", message="This class will be removed in 0.17")
class ContactAnalysis1(object):
    """Perform a very flexible native contact analysis with respect to a single
    reference.

    This analysis class allows one to calculate the fraction of native contacts
    *q* between two arbitrary groups of atoms with respect to an arbitrary
    reference structure. For instance, as a reference one could take a crystal
    structure of a complex, and as the two groups atoms one selects two
    molecules A and B in the complex. Then the question to be answered by *q*
    is, is which percentage of the contacts between A and B persist during the
    simulation.

    First prepare :class:`~MDAnalysis.core.AtomGroup.AtomGroup` selections for
    the reference atoms; this example uses some arbitrary selections::

      ref = Universe('crystal.pdb')
      refA = ref.select_atoms('name CA and segid A and resid 6:100')
      refB = ref.select_atoms('name CA and segid B and resid 1:40')

    Load the trajectory::

      u = Universe(topology, trajectory)

    We then need two selection strings *selA* and *selB* that, when applied as
    ``u.select_atoms(selA)`` produce a list of atoms that is equivalent to the
    reference (i.e. ``u.select_atoms(selA)`` must select the same atoms as
    ``refA`` in this example)::

      selA = 'name CA and resid 1:95'     # corresponds to refA
      selB = 'name CA and resid 150:189'  # corresponds to refB

    .. Note::

       It is the user's responsibility to provide a reference group
       (or groups) that describe equivalent atoms to the ones selected
       by *selection*.

    Now we are ready to set up the analysis::

      CA1 = ContactAnalysis1(u, selection=(selA,selB), refgroup=(refA,refB),
                             radius=8.0, outfile="q.dat")

    If the groups do not match in length then a :exc:`ValueError` is raised.

    The analysis across the whole trajectory is performed with ::

      CA1.run()

    Results are saved to *outfile* (``framenumber q N`` per line) and
    can also be plotted with ::

      CA1.plot()        # plots the time series q(t)
      CA1.plot_qavg()   # plots the matrix of average contacts <q>

    **Description of computed values** in the output file:

    *N*
         number of native contacts

    *q*
         fraction of native contacts relative to the reference

    """

    def __init__(self, *args, **kwargs):
        """Calculate native contacts within a group or between two groups.

        :Arguments:
          *topology*
            psf or pdb file
          *trajectory*
            dcd or xtc/trr file
          *universe*
            instead of a topology/trajectory combination, one can also supply
            a :class:`MDAnalysis.Universe`

        :Keywords:
          *selection*
            selection string that determines which distances are calculated; if
            this is a tuple or list with two entries then distances are
            calculated between these two different groups ["name CA or name
            B*"]
          *refgroup*
            reference group, either a single
            :class:`~MDAnalysis.core.AtomGroup.AtomGroup` (if there is only a
            single *selection*) or a list of two such groups. The reference
            contacts are directly computed from *refgroup* and hence the atoms
            in the reference group(s) must be equivalent to the ones produced
            by the *selection* on the input trajectory.
          *radius*
            contacts are deemed any atoms within radius [8.0 A]
          *outfile*
            name of the output file; with the gz or bz2 suffix, a compressed
            file is written. The average <q> is written to a second, gzipped
            file that has the same name with 'array' included. E.g. for the
            default name "q1.dat.gz" the <q> file will be "q1.array.gz". The
            format is the matrix in column-row format, i.e. selection 1
            residues are the columns and selection 2 residues are rows. The
            file can be read with :func:`np.loadtxt`.  ["q1.dat.gz"]

        The function calculates the percentage of native contacts q1
        along a trajectory. "Contacts" are defined as the number of atoms
        within *radius* of a given other atom. *q1* is the fraction of contacts
        relative to the reference state 1 (typically the starting conformation
        of the trajectory).

        The timeseries is written to a file *outfile* and is also accessible as
        the attribute :attr:`ContactAnalysis1.timeseries`.

        .. deprecated: 0.14.0

        """

        # XX or should I use as input
        #   sel = (group1, group2), ref = (refgroup1, refgroup2)
        # and get the universe from sel?
        # Currently it's a odd hybrid.
        #
        # Enhancements:
        # - select contact pairs to write out as a timecourse
        # - make this selection based on qavg
        from os.path import splitext

        warnings.warn("ContactAnalysis1 is deprecated and will be removed "
                      "in 1.0. Use Contacts instead.",
                      category=DeprecationWarning)

        self.selection_strings = self._return_tuple2(kwargs.pop(
            'selection', "name CA or name B*"), "selection")
        self.references = self._return_tuple2(kwargs.pop('refgroup', None),
                                              "refgroup")
        self.radius = kwargs.pop('radius', 8.0)
        self.targetdir = kwargs.pop('targetdir', os.path.curdir)
        self.output = kwargs.pop('outfile', "q1.dat.gz")
        self.outarray = splitext(splitext(self.output)[0])[0] + ".array.gz"
        self.force = kwargs.pop('force', False)

        self.timeseries = None  # final result

        self.filenames = args
        self.universe = MDAnalysis.as_Universe(*args, **kwargs)

        self.selections = [self.universe.select_atoms(s)
                           for s in self.selection_strings]

        # sanity checkes
        for x in self.references:
            if x is None:
                raise ValueError("a reference AtomGroup must be supplied")
        for ref, sel, s in zip(self.references,
                               self.selections,
                               self.selection_strings):
            if ref.atoms.n_atoms != sel.atoms.n_atoms:
                raise ValueError("selection=%r: Number of atoms differ "
                                 "between reference (%d) and trajectory (%d)" %
                                 (s, ref.atoms.n_atoms, sel.atoms.n_atoms))

        # compute reference contacts
        dref = MDAnalysis.lib.distances.distance_array(
            self.references[0].positions, self.references[1].positions)
        self.qref = self.qarray(dref)
        self.nref = self.qref.sum()

        # setup arrays for the trajectory
        self.d = np.zeros_like(dref)
        self.q = self.qarray(self.d)
        self._qtmp = np.zeros_like(self.q)  # pre-allocated array

        self.qavg = np.zeros(shape=self.q.shape, dtype=np.float64)

    def _return_tuple2(self, x, name):
        if not isinstance(x, (tuple, list, np.ndarray)):
            t = (x,)
        else:
            t = x
        if len(t) == 2:
            return t
        elif len(t) == 1:
            return (x, x)
        else:
            raise ValueError("%(name)s must be a single object or a "
                             "tuple/list with two objects and not %(x)r" % vars())

    def output_exists(self, force=False):
        """Return True if default output file already exists.

        Disable with force=True (will always return False)
        """
        return os.path.isfile(self.output) and not (self.force or force)

    def run(self, store=True, force=False, start=0, stop=None, step=1,
            **kwargs):
        """Analyze trajectory and produce timeseries.

        Stores results in :attr:`ContactAnalysis1.timeseries` (if store=True)
        and writes them to a data file. The average q is written to a second
        data file.
        *start*
            The value of the first frame index in the trajectory to be used
            (default: index 0)
        *stop*
            The value of the last frame index in the trajectory to be used
            (default: None -- use all frames)
        *step*
            The number of frames to skip during trajectory iteration (default:
            use every frame)

        """

        if 'start_frame' in kwargs:
            warnings.warn("start_frame argument has been deprecated, use "
                          "start instead -- removal targeted for version "
                          "0.15.0", DeprecationWarning)
            start = kwargs.pop('start_frame')

        if 'end_frame' in kwargs:
            warnings.warn("end_frame argument has been deprecated, use "
                          "stop instead -- removal targeted for version "
                          "0.15.0", DeprecationWarning)
            stop = kwargs.pop('end_frame')

        if 'step_value' in kwargs:
            warnings.warn("step_value argument has been deprecated, use "
                          "step instead -- removal targeted for version "
                          "0.15.0", DeprecationWarning)
            step = kwargs.pop('step_value')

        if self.output_exists(force=force):
            warnings.warn("File %r already exists, loading it INSTEAD of "
                          "trajectory %r. Use force=True to overwrite "
                          "the output file. " %
                          (self.output, self.universe.trajectory.filename))
            self.load(self.output)
            return None

        with openany(self.output, 'w') as out:
            out.write("# q1 analysis\n# nref = {0:d}\n".format((self.nref)))
            out.write("# frame  q1  n1\n")
            records = []
            self.qavg *= 0  # average contact existence
            A, B = self.selections
            for ts in self.universe.trajectory[start:stop:step]:
                frame = ts.frame
                # use pre-allocated distance array to save a little bit of time
                MDAnalysis.lib.distances.distance_array(A.coordinates(),
                                                        B.coordinates(),
                                                        result=self.d)
                self.qarray(self.d, out=self.q)
                n1, q1 = self.qN(self.q, out=self._qtmp)
                self.qavg += self.q
                if store:
                    records.append((frame, q1, n1))
                out.write("{frame:4d}  {q1:8.6f} {n1:5d}\n".format(**vars()))
        if store:
            self.timeseries = np.array(records).T
        n_frames = len(np.arange(
            self.universe.trajectory.n_frames)[start:stop:step])
        if n_frames > 0:
            self.qavg /= n_frames
        else:
            logger.warn("No frames were analyzed. "
                        "Check values of start, stop, step.")
            logger.debug("start={start} stop={stop} step={step}".format(**vars()))
        np.savetxt(self.outarray, self.qavg, fmt="%8.6f")
        return self.output

    def qarray(self, d, out=None):
        """Return distance array with True for contacts.

        *d* is the matrix of distances. The method uses the value of
        :attr:`ContactAnalysis1.radius` to determine if a ``distance < radius``
        is considered a contact.

        If *out* is supplied as a pre-allocated array of the correct
        shape then it is filled instead of allocating a new one in
        order to increase performance.

        This method is typically only used internally.
        """
        if out is None:
            out = (d <= self.radius)
        else:
            out[:] = (d <= self.radius)
        return out

    def qN(self, q, out=None):
        """Calculate native contacts relative to reference state.

        *q* is the matrix of contacts (e.g. :attr:`~ContactAnalysis1.q`).

        If *out* is supplied as a pre-allocated array of the correct
        shape then it is filled instead of allocating a new one in
        order to increase performance.

        This method is typically only used internally.
        """
        if out is None:
            out = np.logical_and(q, self.qref)
        else:
            np.logical_and(q, self.qref, out)
        contacts = out.sum()
        return contacts, float(contacts) / self.nref

    def load(self, filename):
        """Load the data file."""
        records = []
        with openany(filename) as data:
            for line in data:
                if line.startswith('#'):
                    continue
                records.append(map(float, line.split()))
        self.timeseries = np.array(records).T
        try:
            self.qavg = np.loadtxt(self.outarray)
        except IOError as err:
            if err.errno != errno.ENOENT:
                raise

    def plot(self, filename=None, **kwargs):
        """Plot q(t).

        .. function:: ContactAnalysis1.plot([filename, ...])

        If *filename* is supplied then the figure is also written to file (the
        suffix determines the file type, e.g. pdf, png, eps, ...). All other
        keyword arguments are passed on to :func:`pylab.plot`.
        """
        from pylab import plot, xlabel, ylabel, savefig

        kwargs.setdefault('color', 'black')
        kwargs.setdefault('linewidth', 2)
        if self.timeseries is None:
            raise ValueError("No timeseries data; "
                             "do 'ContactAnalysis.run(store=True)' first.")
        t = self.timeseries
        plot(t[0], t[1], **kwargs)
        xlabel(r"frame number $t$")
        ylabel(r"native contacts $q_1$")

        if filename is not None:
            savefig(filename)

    def _plot_qavg_pcolor(self, filename=None, **kwargs):
        """Plot :attr:`ContactAnalysis1.qavg`, the matrix of average native
        contacts."""
        from pylab import (pcolor, gca, meshgrid, xlabel, ylabel, xlim, ylim,
                           colorbar, savefig)

        x, y = self.selections[0].resids, self.selections[1].resids
        X, Y = meshgrid(x, y)

        pcolor(X, Y, self.qavg.T, **kwargs)
        gca().set_aspect('equal')

        xlim(min(x), max(x))
        ylim(min(y), max(y))

        xlabel("residues")
        ylabel("residues")

        colorbar()

        if filename is not None:
            savefig(filename)

    def plot_qavg(self, filename=None, **kwargs):
        """Plot :attr:`ContactAnalysis1.qavg`, the matrix of average native contacts.

        .. function:: ContactAnalysis1.plot_qavg([filename, ...])

        If *filename* is supplied then the figure is also written to file (the
        suffix determines the file type, e.g. pdf, png, eps, ...). All other
        keyword arguments are passed on to :func:`pylab.imshow`.
        """
        from pylab import (imshow, xlabel, ylabel, xlim, ylim, colorbar, cm,
                           clf, savefig)

        x, y = self.selections[0].resids, self.selections[1].resids

        kwargs['origin'] = 'lower'
        kwargs.setdefault('aspect', 'equal')
        kwargs.setdefault('interpolation', 'nearest')
        kwargs.setdefault('vmin', 0)
        kwargs.setdefault('vmax', 1)
        kwargs.setdefault('cmap', cm.hot)
        kwargs.setdefault('extent', (min(x), max(x), min(y), max(y)))

        clf()
        imshow(self.qavg.T, **kwargs)

        xlim(min(x), max(x))
        ylim(min(y), max(y))

        xlabel("residue from {0!r}".format(self.selection_strings[0]))
        ylabel("residue from {0!r}".format(self.selection_strings[1]))

        colorbar()

        if filename is not None:
            savefig(filename)
