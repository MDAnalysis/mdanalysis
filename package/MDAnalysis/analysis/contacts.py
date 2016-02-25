# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""
Native contacts analysis --- :mod:`MDAnalysis.analysis.contacts`
================================================================

:Author: Oliver Beckstein
:Year: 2010
:Copyright: GNU Public License v3

Analysis of native contacts *q* over a trajectory.

* a "contact" exists between two atoms *i* and *j* if the distance between them is
  smaller than a given *radius*

* a "native contact" exists between *i* and *j* if a contact exists and if the
  contact also exists between the equivalent atoms in a reference structure or
  conformation

The "fraction of native contacts" *q(t)* is a number between 0 and 1 and
calculated as the total number of native contacts for a given time frame
divided by the total number of contacts in the reference structure.

Classes are available for two somewhat different ways to perform a contact
analysis:

1. Contacts between two groups of atoms are defined with
   :class:`ContactAnalysis1`), which allows one to calculate *q(t)* over
   time. This is especially useful in order to look at native contacts during
   an equilibrium simulation where one can also look at the average matrix of
   native contacts (see :meth:`ContactAnalysis1.plot_qavg`).

2. Contacts are defined within one group in a protein (e.g. all C-alpha atoms)
   but relative to *two different conformations* 1 and 2, using
   :class:`ContactAnalysis`. This allows one to do a *q1-q2* analysis that
   shows how native contacts of state 1 change in comparison to native contacts
   of state 2.  Transition pathways have been analyzed in terms of these two
   variables q1 and q2 that relate to the native contacts in the end states of
   the transition.

.. SeeAlso:: See http://lorentz.dynstr.pasteur.fr/joel/adenylate.php for an
   example of contact analysis applied to MinActionPath trajectories of AdK
   (although this was *not* performed with MDAnalysis --- it's provided as a
   very good illustrative example).


Examples
--------

One-dimensional contact analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As an example we analyze the opening ("unzipping") of salt bridges
when the AdK enzyme opens up; this is one of the example trajectories
in MDAnalysis. ::

    import MDAnalysis
    import MDAnalysis.analysis.contacts
    from MDAnalysis.tests.datafiles import PSF,DCD

    # example trajectory (transition of AdK from closed to open)
    u = MDAnalysis.Universe(PSF,DCD)

    # crude definition of salt bridges as contacts between NH/NZ in ARG/LYS and OE*/OD* in ASP/GLU.
    # You might want to think a little bit harder about the problem before using this for real work.
    sel_basic = "(resname ARG or resname LYS) and (name NH* or name NZ)"
    sel_acidic = "(resname ASP or resname GLU) and (name OE* or name OD*)"

    # reference groups (first frame of the trajectory, but you could also use a separate PDB, eg crystal structure)
    acidic = u.select_atoms(sel_acidic)
    basic = u.select_atoms(sel_basic)

    # set up analysis of native contacts ("salt bridges"); salt bridges have a distance <6 A
    CA1 = MDAnalysis.analysis.contacts.ContactAnalysis1(u, selection=(sel_acidic, sel_basic), refgroup=(acidic,
    basic), radius=6.0, outfile="qsalt.dat")

    # iterate through trajectory and perform analysis of "native contacts" q
    # (force=True ignores any previous results, force=True is useful when testing)
    CA1.run(force=True)

    # plot time series q(t) [possibly do "import pylab; pylab.clf()" do clear the figure first...]
    CA1.plot(filename="adk_saltbridge_contact_analysis1.pdf", linewidth=3, color="blue")

    # or plot the data in qsalt.dat yourself.
    CA1.plot_qavg(filename="adk_saltbridge_contact_analysis1_matrix.pdf")

The first graph shows that when AdK opens, about 20% of the salt
bridges that existed in the closed state disappear when the enzyme
opens. They open in a step-wise fashion (made more clear by the movie
`AdK_zipper_cartoon.avi`_).

The output graphs can be made prettier but if you look at the code
itself then you'll quickly figure out what to do. The qavg plot is the
matrix of all contacts, averaged over the trajectory. This plot makes
more sense for an equilibrium trajectory than for the example above
but is is included for illustration.

See the docs for :class:`ContactAnalysis1` for another example.


.. AdK_zipper_cartoon.avi:
   http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2803350/bin/NIHMS150766-supplement-03.avi

Two-dimensional contact analysis (q1-q2)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Analyze a single DIMS transition of AdK between its closed and open
conformation and plot the trajectory projected on q1-q2::

  import MDAnalysis.analysis.contacts
  from MDAnalysis.tests.datafiles import PSF, DCD
  C = MDAnalysis.analysis.contacts.ContactAnalysis(PSF, DCD)
  C.run()
  C.plot()

Compare the resulting pathway to the `MinActionPath result for AdK`_.

.. _MinActionPath result for AdK:
   http://lorentz.dynstr.pasteur.fr/joel/adenylate.php

Classes
-------

.. autoclass:: ContactAnalysis
   :members:
.. autoclass:: ContactAnalysis1
   :members:

"""

import os
import errno
import warnings
import bz2
from six.moves import zip
import numpy as np
import logging

import MDAnalysis
import MDAnalysis.lib.distances
from MDAnalysis.lib.util import openany
from MDAnalysis.analysis.distances import distance_array
from .base import AnalysisBase

logger = logging.getLogger("MDAnalysis.analysis.contacts")

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
    """

    def __init__(self, topology, trajectory, ref1=None, ref2=None, radius=8.0,
                 targetdir=os.path.curdir, infix="", force=False,
                 selection="name CA", centroids=False):
        """Calculate native contacts from two reference structures.

        :Arguments:
          *topology*
            psf or pdb file
          *trajectory*
            dcd or xtc/trr file
          *ref1*
            structure of the reference conformation 1 (pdb); if ``None`` the *first*
            frame of the trajectory is chosen
          *ref2*
            structure of the reference conformation 2 (pdb); if ``None`` the *last*
            frame of the trajectory is chosen
          *radius*
            contacts are deemed any Ca within radius [8 A]
          *targetdir*
            output files are saved there [.]
          *infix*
            additional tag string that is inserted into the output filename of the
            data file [""]
          *selection*
            MDAnalysis selection string that selects the particles of
            interest; the default is to only select the C-alpha atoms
            in *ref1* and *ref*2 ["name CA"]

            .. Note:: If *selection* produces more than one atom per
                      residue then you will get multiple contacts per
                      residue unless you also set *centroids* = ``True``
          *centroids*
            If set to ``True``, use the centroids for the selected atoms on a
            per-residue basis to compute contacts. This allows, for instance
            defining the sidechains as *selection* and then computing distances
            between sidechain centroids.

        The function calculates the percentage of native contacts *q1* and *q2*
        along a trajectory. "Contacts" are defined as the number of Ca atoms (or
        per-residue *centroids* of a user defined *selection*) within *radius* of
        a primary Ca. *q1* is the fraction of contacts relative to the reference
        state 1 (typically the starting conformation of the trajectory) and *q2*
        is the fraction of contacts relative to the conformation 2.

        The timeseries is written to a bzip2-compressed file in *targetdir*
        named "basename(*trajectory*)*infix*_q1q2.dat.bz2" and is also
        accessible as the attribute
        :attr:`ContactAnalysis.timeseries`.
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
            return  # do not bother reading any data or initializing arrays... !!
        # don't bother if trajectory is empty (can lead to segfaults so better catch it)
        stats = os.stat(trajectory)
        if stats.st_size == 0:
            warnings.warn('trajectory = {trajectory!s} is empty, skipping...'.format(**vars()))
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

        :Keywords:
           *results*
              passed on to :func:`MDAnalysis.lib.distances.self_distance_array`
              as a preallocated array
           *centroids*
              ``True``: calculate per-residue centroids from the selected atoms;
              ``False``: consider each atom separately; ``None``: use the class
              default for *centroids* [``None``]

        """
        centroids = kwargs.pop("centroids", None)
        centroids = self.centroids if centroids is None else centroids
        if not centroids:
            coordinates = g.positions
        else:
            # centroids per residue (but only including the selected atoms)
            coordinates = np.array([residue.centroid() for residue in g.split("residue")])
        return MDAnalysis.lib.distances.self_distance_array(coordinates, **kwargs)

    def output_exists(self, force=False):
        """Return True if default output file already exists.

        Disable with force=True (will always return False)
        """
        return (os.path.isfile(self.output) or os.path.isfile(self.output_bz2)) and not (self.force or force)

    def run(self, store=True, force=False):
        """Analyze trajectory and produce timeseries.

        Stores results in :attr:`ContactAnalysis.timeseries` (if
        store=True) and writes them to a bzip2-compressed data file.
        """
        if self._skip or self.output_exists(force=force):
            warnings.warn("File {output!r} or {output_bz2!r} already exists, loading {trajectory!r}.".format(**vars(self)))
            try:
                self.load(self.output)
            except IOError:
                self.load(self.output_bz2)
            return None

        outbz2 = bz2.BZ2File(self.output_bz2, mode='w', buffering=8192)
        try:
            outbz2.write("# q1-q2 analysis\n# nref1 = {0:d}\n# nref2 = {1:d}\n".format(self.nref[0], self.nref[1]))
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
        """Return distance array with True for contacts.

        If *out* is supplied as a pre-allocated array of the correct
        shape then it is filled instead of allocating a new one in
        order to increase performance.
        """
        if out is None:
            out = (d <= self.radius)
        else:
            out[:] = (d <= self.radius)
        return out

    def qN(self, q, n, out=None):
        """Calculate native contacts relative to state n.

        If *out* is supplied as a pre-allocated array of the correct
        shape then it is filled instead of allocating a new one in
        order to increase performance.
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
            raise ValueError("No timeseries data; do 'ContactAnalysis.run(store=True)' first.")
        t = self.timeseries
        plot(t[1], t[2], **kwargs)
        xlabel(r"$q_1$")
        ylabel(r"$q_2$")



# ContactAnalysis1 is a (hopefully) temporary hack. It should be unified with ContactAnalysis
# or either should be derived from a base class because many methods are copy&paste with
# minor changes (mostly for going from q1q2 -> q1 only).
# If ContactAnalysis is enhanced to accept two references then this should be even easier.
# It might also be worthwhile making a simpler class that just does the q calculation
# and use it for both reference and trajectory data.
class ContactAnalysis1(object):
    """Perform a very flexible native contact analysis with respect to a single reference.

    This analysis class allows one to calculate the fraction of native contacts
    *q* between two arbitrary groups of atoms with respect to an arbitrary
    reference structure. For instance, as a reference one could take a crystal
    structure of a complex, and as the two groups atoms one selects two
    molecules A and B in the complex. Then the question to be answered by *q*
    is, is which percentage of the contacts between A and B persist during the simulation.

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

      CA1 = ContactAnalysis1(u, selection=(selA,selB), refgroup=(refA,refB), radius=8.0, outfile="q.dat")

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
            selection string that determines which distances are calculated; if this
            is a tuple or list with two entries then distances are calculated between
            these two different groups ["name CA or name B*"]
          *refgroup*
            reference group, either a single :class:`~MDAnalysis.core.AtomGroup.AtomGroup`
            (if there is only a single *selection*) or a list of two such groups.
            The reference contacts are directly computed from *refgroup* and hence
            the atoms in the reference group(s) must be equivalent to the ones produced
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
        
        warnings.warn("ContactAnalysis1 is deprecated and will be removed in 1.0/. Use Contacts instead.", DeprecationWarning)

        self.selection_strings = self._return_tuple2(kwargs.pop('selection', "name CA or name B*"), "selection")
        self.references = self._return_tuple2(kwargs.pop('refgroup', None), "refgroup")
        self.radius = kwargs.pop('radius', 8.0)
        self.targetdir = kwargs.pop('targetdir', os.path.curdir)
        self.output = kwargs.pop('outfile', "q1.dat.gz")
        self.outarray = splitext(splitext(self.output)[0])[0] + ".array.gz"
        self.force = kwargs.pop('force', False)

        self.timeseries = None  # final result

        self.filenames = args
        self.universe = MDAnalysis.as_Universe(*args, **kwargs)

        self.selections = [self.universe.select_atoms(s) for s in self.selection_strings]

        # sanity checkes
        for x in self.references:
            if x is None:
                raise ValueError("a reference AtomGroup must be supplied")
        for ref, sel, s in zip(self.references, self.selections, self.selection_strings):
            if ref.atoms.n_atoms != sel.atoms.n_atoms:
                raise ValueError("selection=%r: Number of atoms differ between "
                                 "reference (%d) and trajectory (%d)" %
                                 (s, ref.atoms.n_atoms, sel.atoms.n_atoms))

        # compute reference contacts
        dref = MDAnalysis.lib.distances.distance_array(
            self.references[0].coordinates(), self.references[1].coordinates())
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
            raise ValueError("%(name)s must be a single object or a tuple/list with two objects "
                             "and not %(x)r" % vars())

    def output_exists(self, force=False):
        """Return True if default output file already exists.

        Disable with force=True (will always return False)
        """
        return os.path.isfile(self.output) and not (self.force or force)

    def run(self, store=True, force=False, start=0, stop=None, step=1, **kwargs):
        """Analyze trajectory and produce timeseries.

        Stores results in :attr:`ContactAnalysis1.timeseries` (if store=True)
        and writes them to a data file. The average q is written to a second
        data file.
        *start*
            The value of the first frame index in the trajectory to be used (default: index 0)
        *stop*
            The value of the last frame index in the trajectory to be used (default: None -- use all frames)
        *step*
            The number of frames to skip during trajectory iteration (default: use every frame)
        """

        if 'start_frame' in kwargs:
            warnings.warn("start_frame argument has been deprecated, use start instead --"
                           "removal targeted for version 0.15.0", DeprecationWarning)
            start = kwargs.pop('start_frame')

        if 'end_frame' in kwargs:
            warnings.warn("end_frame argument has been deprecated, use stop instead --"
                           "removal targeted for version 0.15.0", DeprecationWarning)
            stop = kwargs.pop('end_frame')

        if 'step_value' in kwargs:
            warnings.warn("step_value argument has been deprecated, use step instead --"
                           "removal targeted for version 0.15.0", DeprecationWarning)
            step = kwargs.pop('step_value')

        if self.output_exists(force=force):
            warnings.warn("File %r already exists, loading it INSTEAD of trajectory %r. "
                          "Use force=True to overwrite the output file. " %
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
                MDAnalysis.lib.distances.distance_array(A.coordinates(), B.coordinates(), result=self.d)
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
            logger.warn("No frames were analyzed. Check values of start, stop, step.") 
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
            raise ValueError("No timeseries data; do 'ContactAnalysis.run(store=True)' first.")
        t = self.timeseries
        plot(t[0], t[1], **kwargs)
        xlabel(r"frame number $t$")
        ylabel(r"native contacts $q_1$")

        if filename is not None:
            savefig(filename)

    def _plot_qavg_pcolor(self, filename=None, **kwargs):
        """Plot :attr:`ContactAnalysis1.qavg`, the matrix of average native contacts."""
        from pylab import pcolor, gca, meshgrid, xlabel, ylabel, xlim, ylim, colorbar, savefig

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
        from pylab import imshow, xlabel, ylabel, xlim, ylim, colorbar, cm, clf, savefig

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


def best_hummer_q(r, r0, beta=5.0, lambda_constant=1.8):
    """Calculate the Best-Hummer fraction of native contacts (Q)

    A soft-cutoff contacts analysis

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
      Fraction of native contacts
    result : array
      Intermediate, r-r0 array transformed by the switching function
    """
    result = 1/(1 + np.exp(beta*(r - lambda_constant * r0)))

    return result.sum() / len(r0), result


class Contacts(AnalysisBase):
    """Calculate fraction of native contacts (Q) from a trajectory

    Inputs
    ------
    Two string selections for the contacting AtomGroups, 
    the groups could be protein-lipid or protein-protein. 

    Use two reference AtomGroups to obtain reference distances (r0)
    for the cotacts. 

    Methods available
    -----------------
    Supports either hard-cutoff or soft-cutoff (Best-Hummer like) contacts. 

    Returns
    -------
    list 
        Returns a list of following structure::
            {
                [[t1, q1], [t2, q2], ... [tn, qn]]
            }
        where t is time in ps and q is the fraction of native contacts

    Examples
    --------

    1. Protein folding

        ref = Universe("villin.gro")
        u = Universe("conf_protein.gro", "traj_protein.xtc")
        Q = calculate_contacts(u, ref, "protein and not name H*", "protein and not name H*")

    2. A pair of helices

        ref = Universe("glycophorin_dimer.pdb")
        u = Universe("conf_protein.gro", "traj_protein.xtc")
        Q = calculate_contacts(u, ref, "protein and resid 75-92 and not name H* and segid A", "protein and resid 75-92 and not name H* and segid B")

    Parater choices
    ---------------
    There are recommendations and not strict orders. 
    These parameters should be insensitive to small changes. 
    * For all-atom simulations, radius = 4.5 A and lambda_constant = 1.8 (unitless)
    * For coarse-grained simulations, radius = 6.0 A and lambda_constant = 1.5 (unitless)

    Additional
    ----------
    Supports writing and reading the analysis results to and from a text file. 
    Supports simple plotting operations, for exploratory data analysis. 

    Reference
    ---------
    Using the definition from Best, Hummer, and Eaton, "Native contacts determine protein folding mechanisms in atomistic simulations" PNAS (2013) 10.1073/pnas.1311599110        

    Eq. (1) of the SI defines the expression for the fraction of native contacts, $Q(X)$:
    .. math:: 

        Q(X) = \frac{1}{|S|} \sum_{(i,j) \in S} \frac{1}{1 + \exp[\beta(r_{ij}(X) - \lambda r_{ij}^0)]}

    where:

    * :math:`X` is a conformation,
    * :math:`r_{ij}(X)` is the distance between atoms $i$ and $j$ in conformation $X$,
    * :math:`r^0_{ij}` is the distance from heavy atom i to j in the native state conformation,
    * :math:`S` is the set of all pairs of heavy atoms $(i,j)$ belonging to residues $\theta_i$ and $\theta_j$ such that $|\theta_i - \theta_j| > 3$ and $r^0_{i,} < 4.5 \unicode{x212B}$,
    * :math:`\beta=5 \unicode{x212B}^{-1},
    * :math:`\lambda=1.8` for all-atom simulations


    """
    def __init__(self, u, selection, refgroup, method="cutoff", radius=4.5, outfile=None,
                 start=None, stop=None, step=None, **kwargs):
        """Calculate the persistence length for polymer chains

        Parameters
        ----------
        u: Universe
            trajectory
        selection: tuple(string, string)
            two contacting groups that change over time
        refgroup: tuple(AtomGroup, AtomGroup)
            two contacting groups in their reference conformation
        radius: float, optional (4.5 Angstroms)
            radius within which contacts exist
        method: string
            either 'cutoff' or 'best-hummer'

        start : int, optional
            First frame of trajectory to analyse, Default: 0
        stop : int, optional
            Last frame of trajectory to analyse, Default: -1
        step : int, optional
            Step between frames to analyse, Default: 1

        Parameters for 'best-hummer' method
        -----------------------------------
        lambda_constant: float, optional (1.8 unitless)
            contact is considered formed between (lambda*r0,r0)
        beta: float, optional (5 Angstroms^-1)
            softness of the switching function, the lower the softer

        Attributes
        ----------
        results: list 
            Fraction of native contacts for each frame
        """

        # check method
        if not method in ("cutoff", "best-hummer"): 
            raise ValueError("method has to be 'cutoff' or 'best-hummer'")
        self._method = method

        # method-specific parameters
        if method == "best-hummer":
            self.beta = kwargs.get('beta', 5.0)
            self.lambda_constant  = kwargs.get('lambda_constant', 1.8)

        # steup boilerplate
        self.u = u
        self._setup_frames(self.u.trajectory,
                           start=start,
                           stop=stop,
                           step=step)

        self.selection = selection
        grA, grB = u.select_atoms(selection[0]), u.select_atoms(selection[1])
        self.grA, self.grB = grA, grB        
        refA, refB = refgroup

        # contacts formed in reference
        r0 = distance_array(refA.positions, refB.positions)
        self.r0 = r0
        self.mask = r0 < radius

        self.contact_matrix = []
        self.timeseries = []
        self.outfile = outfile

    @staticmethod
    def load(filename):
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

    def _single_frame(self):
        grA, grB, r0, mask = self.grA, self.grB, self.r0, self.mask

        # compute distance array for a frame
        d = distance_array(grA.positions, grB.positions)

        # select only the contacts that were formed in the reference state
        # r, r0 are 1D array
        r, r0 = d[mask], r0[mask]

        if self._method == "cutoff":
            y = r <= r0
            y = float(y.sum())/mask.sum()
        elif self._method == "best-hummer":
            y, _ = best_hummer_q(r, r0, self.beta, self.lambda_constant)
        else:
            raise ValueError("Unknown method type, has to be 'cutoff' or 'best-hummer'")

        cm = np.zeros((grA.positions.shape[0], grB.positions.shape[0]))
        cm[mask] = y
        self.contact_matrix.append(cm)
        self.timeseries.append((self._ts.frame , y, mask.sum()))

    def _conclude(self):
        """Finalise the timeseries you've gathered.

        Called at the end of the run() method to finish everything up.
        """
        # write output
        if not self.outfile: return
        with open(self.outfile, "w") as f:
            f.write("# q1 analysis\n# nref = {0:d}\n".format(self.mask.sum()))
            f.write("# frame  q1  n1\n")
            for frame, q1, n1 in self.timeseries:
                f.write("{frame:4d}  {q1:8.6f} {n1:5d}\n".format(**vars()))            

    def contact_matrix(self, d, out=None):
        """Return distance array with True for contacts.

        Parameters
        ----------
        d : array
            is the matrix of distances. The method uses the value of
            `ContactAnalysis1.radius` to determine if a ``distance < radius``
            is considered a contact.
        out: array (optional)
            If `out` is supplied as a pre-allocated array of the correct
            shape then it is filled instead of allocating a new one in
            order to increase performance.

        Returns
        -------
        array 
            boolean array of which contacts are formed

        Note
        ----
        This method is typically only used internally.
        """
        if out:
            out[:] = (d <= self.radius)
        else:
            out = (d <= self.radius)            
        return out

    def fraction_native(q, out=None):
        """Calculate native contacts relative to reference state.

        Parameters
        ----------
        q: array
            is the matrix of contacts (e.g. `ContactAnalysis1.q`).
        out: array
            If *out* is supplied as a pre-allocated array of the correct
            shape then it is filled instead of allocating a new one in
            order to increase performance.

        Returns
        -------
        array
            Fraction of native contacts (Q) calculated from a contact matrix

        Note
        ----
        This method is typically only used internally.
        """
        if out:
            np.logical_and(q, self.mask, out)
        else:
            out = np.logical_and(q, self.mask)
        contacts = out.sum()
        return contacts, float(contacts) / self.mask.sum()

    def plot(self, filename=None, **kwargs):
        """Plot q(t).

        Parameters
        ----------
        filename : str
            If `filename` is supplied then the figure is also written to file (the
            suffix determines the file type, e.g. pdf, png, eps, ...). All other
            keyword arguments are passed on to `pylab.plot`.
        **kwargs
            Arbitrary keyword arguments for the plotting function            
        """
        if not self.timeseries :
            raise ValueError("No timeseries data; do 'Contacts.run()' first.")
        x, y, _ = zip(*self.timeseries)

        import matplotlib.pyplot as plt
        kwargs.setdefault('color', 'black')
        kwargs.setdefault('linewidth', 2)       

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(x, y, **kwargs)
        ax.set_xlabel(r"frame number $t$")
        ax.set_ylabel(r"contacts $q_1$")

        if filename:
            fig.savefig(filename)            
        else:
            fig.show()


    def plot_qavg(self, filename=None, **kwargs):
        """Plot `Contacts.qavg`, the matrix of average contacts.

        Parameters
        ----------
        filename : str
            If `filename` is supplied then the figure is also written to file (the
            suffix determines the file type, e.g. pdf, png, eps, ...). All other
            keyword arguments are passed on to `pylab.imshow`.
        **kwargs
            Arbitrary keyword arguments for the plotting function
        """
        if not self.contact_matrix :
            raise ValueError("No timeseries data; do 'Contacts.run()' first.")
        # collapse on the time-axis
        data = np.array(self.contact_matrix)
        data = data.mean(axis=0)

        import matplotlib.pyplot as plt
        import matplotlib.cm as cm

        kwargs['origin'] = 'lower'
        kwargs.setdefault('aspect', 'equal')
        kwargs.setdefault('interpolation', 'nearest')
        kwargs.setdefault('vmin', 0)
        kwargs.setdefault('vmax', 1)
        kwargs.setdefault('cmap', cm.hot)

        fig = plt.figure()
        ax = fig.add_subplot(111)
        cax = ax.imshow(data, **kwargs)

        cbar = fig.colorbar(cax)

        if filename:
            fig.savefig(filename)            
        else: 
            fig.show()