# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2011 Naveen Michaud-Agrawal,
#               Elizabeth J. Denning, Oliver Beckstein,
#               and contributors (see website for details)
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
#     N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
#     O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
#     Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
#     doi:10.1002/jcc.21787
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
    acidic = u.selectAtoms(sel_acidic)
    basic = u.selectAtoms(sel_basic)

    # set up analysis of native contacts ("salt bridges"); salt bridges have a distance <6 A
    CA1 = MDAnalysis.analysis.contacts.ContactAnalysis1(u, selection=(sel_acidic, sel_basic), refgroup=(acidic, basic), radius=6.0, outfile="qsalt.dat")

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
http://sbcb.bioch.ox.ac.uk/oliver/Movies/AdK/AdK_zipper_cartoon.avi
(divx, on Mac use http://perian.org)).

The output graphs can be made prettier but if you look at the code
itself then you'll quickly figure out what to do. The qavg plot is the
matrix of all contacts, averaged over the trajectory. This plot makes
more sense for an equilibrium trajectory than for the example above
but is is included for illustration.

See the docs for :class:`ContactAnalysis1` for another example.


Two-dimensional contact analysis (q1-q2)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Analyze a single DIMS transition of AdK between its closed and open
conformation and plot the trajectory projected on q1-q2::

  import MDAnalysis.analysis.contacts
  from MDAnalysis.tests.datafiles import *
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
from __future__ import with_statement

import os, errno
import warnings
import bz2
from itertools import izip
import numpy
import MDAnalysis
from MDAnalysis.core.distances import distance_array, self_distance_array
from MDAnalysis.core.util import openany

import logging
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
                 targetdir=os.path.curdir, infix="", force=False):
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

        The function calculates the percentage of native contacts q1 and q2
        along a trajectory. "Contacts" are defined as the number of Ca atoms
        within *radius* of a primary Ca. *q1* is the fraction of contacts
        relative to the reference state 1 (typically the starting conformation
        of the trajectory) and *q2* is the fraction of contacts relative to the
        conformation 2.

        The timeseries is written to a bzip2-compressed file in *targetdir*
        named "basename(*trajectory*)*infix*_q1q2.dat.bz2" and is also
        accessible as the attribute :attr:`ContactAnalysis.timeseries`.
        """

        self.topology = topology
        self.trajectory = trajectory
        self.radius = radius
        self.targetdir = targetdir
        self.force = force

        trajectorybase = os.path.splitext(os.path.basename(trajectory))[0]
        output = trajectorybase + infix + '_q1q2.dat'
        self.output = os.path.join(self.targetdir, output)
        self.output_bz2 = self.output + '.bz2'

        self.timeseries = None  # final result

        # short circuit if output file already exists: skip everything
        if self.output_exists():
            self._skip = True
            return   # do not bother reading any data or initializing arrays... !!
        # don't bother if trajectory is empty (can lead to segfaults so better catch it)
        stats = os.stat(trajectory)
        if stats.st_size == 0:
            warnings.warn('trajectory = %(trajectory)s is empty, skipping...' % vars())
            self._skip = True
            return
        # under normal circumstances we do not skip
        self._skip = False

        # expensive initialization starts with building Universes :-)
        self.u = MDAnalysis.Universe(topology, trajectory)

        if ref1 is None:
            ref1 = os.path.join(self.targetdir, trajectorybase + '_first.pdb')
            self.u.trajectory[0]      # extract first frame
            self.u.atoms.write(ref1)
        self.ref1 = ref1
        if ref2 is None:
            ref2 = os.path.join(self.targetdir, trajectorybase + '_last.pdb')
            self.u.trajectory[-1]     # extract last frame
            self.u.atoms.write(ref2)
            self.u.trajectory[0]      # rewind, just in case...
        self.ref2 = ref2

        r1 = MDAnalysis.Universe(topology, self.ref1)
        r2 = MDAnalysis.Universe(topology, self.ref2)

        self.ca = self.u.selectAtoms('name CA')
        ca1 = r1.selectAtoms('name CA')
        ca2 = r2.selectAtoms('name CA')

        # NOTE: self_distance_array() produces a 1D array; this works here
        #       but is not the same as the 2D output from distance_array()!
        #       See the docs for self_distance_array().
        dref =  [self_distance_array(ca1.coordinates()), self_distance_array(ca2.coordinates())]
        self.qref = [self.qarray(dref[0]), self.qarray(dref[1])]
        self.nref = [self.qref[0].sum(), self.qref[1].sum()]

        self.d = numpy.zeros_like(dref[0])
        self.q = self.qarray(self.d)
        self._qtmp = numpy.zeros_like(self.q)  # pre-allocated array

    def output_exists(self, force=False):
        """Return True if default output file already exists.

        Disable with force=True (will always return False)
        """
        return (os.path.isfile(self.output) or os.path.isfile(self.output_bz2)) \
               and not (self.force or force)

    def run(self, store=True, force=False):
        """Analyze trajectory and produce timeseries.

        Stores results in :attr:`ContactAnalysis.timeseries` (if
        store=True) and writes them to a bzip2-compressed data file.
        """
        if self._skip or self.output_exists(force=force):
            import warnings
            warnings.warn("File %(output)r or %(output_bz2)r already exists, loading %(trajectory)r." % vars(self))
            try:
                self.load(self.output)
            except IOError:
                self.load(self.output_bz2)
            return None

        outbz2 = bz2.BZ2File(self.output_bz2, mode='w', buffering=8192)
        try:
            outbz2.write("# q1-q2 analysis\n# nref1 = %d\n# nref2 = %d\n"
                         % (self.nref[0], self.nref[1]))
            outbz2.write("# frame  q1  q2   n1  n2\n")
            records = []
            for ts in self.u.trajectory:
                frame = ts.frame
                # use pre-allocated distance array to save a little bit of time
                self_distance_array(self.ca.coordinates(), result=self.d)
                self.qarray(self.d, out=self.q)
                n1, q1 = self.qN(self.q, 0, out=self._qtmp)
                n2, q2 = self.qN(self.q, 1, out=self._qtmp)

                if store:
                    records.append((frame, q1, q2, n1, n2))
                outbz2.write("%(frame)4d  %(q1)8.6f %(q2)8.6f  %(n1)5d %(n2)5d\n" % vars())
        finally:
            outbz2.close()
        if store:
            self.timeseries = numpy.array(records).T
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
            out = numpy.logical_and(q, self.qref[n])
        else:
            numpy.logical_and(q, self.qref[n], out)
        contacts = out.sum()
        return contacts, float(contacts)/self.nref[n]

    def load(self, filename):
        """Load the data file."""
        from MDAnalysis.core.util import openany
        records = []
        with openany(filename) as data:
            for line in data:
                if line.startswith('#'):
                    continue
                records.append(map(float, line.split()))
        self.timeseries = numpy.array(records).T

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

    .. class:: ContactAnalysis1(topology, trajectory[,selection[,refgroup[,radius[,outfile]]]])

    .. class:: ContactAnalysis1(universe[,selection[,refgroup[,radius[,outfile]]]])


    This analysis class allows one to calculate the fraction of native contacts
    *q* between two arbitrary groups of atoms with respect to an arbitrary
    reference structure. For instance, as a reference one could take a crystal
    structure of a complex, and as the two groups atoms one selects two
    molecules A and B in the complex. Then the question to be answered by *q*
    is, is which percentage of the contacts between A and B persist during the simulation.

    First prepare :class:`~MDAnalysis.core.AtomGroup.AtomGroup` selections for
    the reference atoms; this example uses some arbitrary selections::

      ref = Universe('crystal.pdb')
      refA = re.selectAtoms('name CA and segid A and resid 6:100')
      refB = re.selectAtoms('name CA and segid B and resid 1:40')

    Load the trajectory::

      u = Universe(topology, trajectory)

    We then need two selection strings *selA* and *selB* that, when applied as
    ``u.selectAtoms(selA)`` produce a list of atoms that is equivalent to the
    reference (i.e. ``u.selectAtoms(selA)`` must select the same atoms as
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
            file can be read with :func:`numpy.loadtxt`.  ["q1.dat.gz"]

        The function calculates the percentage of native contacts q1
        along a trajectory. "Contacts" are defined as the number of atoms
        within *radius* of a given other atom. *q1* is the fraction of contacts
        relative to the reference state 1 (typically the starting conformation
        of the trajectory).

        The timeseries is written to a file *outfile* and is also accessible as
        the attribute :attr:`ContactAnalysis1.timeseries`.
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

        self.selection_strings = self._return_tuple2(kwargs.pop('selection', "name CA or name B*"), "selection")
        self.references = self._return_tuple2(kwargs.pop('refgroup', None), "refgroup")
        self.radius = kwargs.pop('radius', 8.0)
        self.targetdir = kwargs.pop('targetdir', os.path.curdir)
        self.output = kwargs.pop('outfile', "q1.dat.gz")
        self.outarray = splitext(splitext(self.output)[0])[0]+".array.gz"
        self.force = kwargs.pop('force', False)

        self.timeseries = None  # final result

        self.filenames = args
        self.universe = MDAnalysis.asUniverse(*args, **kwargs)

        self.selections = [self.universe.selectAtoms(s) for s in self.selection_strings]

        # sanity checkes
        for x in self.references:
            if x is None:
                raise ValueError("a reference AtomGroup must be supplied")
        for ref,sel,s in izip(self.references, self.selections, self.selection_strings):
            if ref.atoms.numberOfAtoms() != sel.atoms.numberOfAtoms():
                raise ValueError("selection=%r: Number of atoms differ between "
                                 "reference (%d) and trajectory (%d)" %
                                 (s, ref.atoms.numberOfAtoms(), sel.atoms.numberOfAtoms()))

        # compute reference contacts
        dref = distance_array(self.references[0].coordinates(), self.references[1].coordinates())
        self.qref = self.qarray(dref)
        self.nref = self.qref.sum()

        # setup arrays for the trajectory
        self.d = numpy.zeros_like(dref)
        self.q = self.qarray(self.d)
        self._qtmp = numpy.zeros_like(self.q)  # pre-allocated array

        self.qavg = numpy.zeros(shape=self.q.shape, dtype=numpy.float64)

    def _return_tuple2(self, x, name):
        if not isinstance(x, (tuple, list, numpy.ndarray)):
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

    def run(self, store=True, force=False):
        """Analyze trajectory and produce timeseries.

        Stores results in :attr:`ContactAnalysis1.timeseries` (if store=True)
        and writes them to a data file. The average q is written to a second
        data file.
        """
        if self.output_exists(force=force):
            import warnings
            warnings.warn("File %r already exists, loading it INSTEAD of trajectory %r. "
                          "Use force=True to overwrite the output file. " %
                          (self.output, self.universe.trajectory.filename))
            self.load(self.output)
            return None

        with openany(self.output, 'w') as out:
            out.write("# q1 analysis\n# nref = %d\n" % (self.nref))
            out.write("# frame  q1  n1\n")
            records = []
            self.qavg *= 0  # average contact existence
            A,B = self.selections
            for ts in self.universe.trajectory:
                frame = ts.frame
                # use pre-allocated distance array to save a little bit of time
                distance_array(A.coordinates(), B.coordinates(), result=self.d)
                self.qarray(self.d, out=self.q)
                n1, q1 = self.qN(self.q, out=self._qtmp)
                self.qavg += self.q
                if store:
                    records.append((frame, q1, n1))
                out.write("%(frame)4d  %(q1)8.6f %(n1)5d\n" % vars())
        if store:
            self.timeseries = numpy.array(records).T
        self.qavg /= self.universe.trajectory.numframes
        numpy.savetxt(self.outarray, self.qavg, fmt="%8.6f")
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
            out = numpy.logical_and(q, self.qref)
        else:
            numpy.logical_and(q, self.qref, out)
        contacts = out.sum()
        return contacts, float(contacts)/self.nref

    def load(self, filename):
        """Load the data file."""
        from MDAnalysis.core.util import openany
        records = []
        with openany(filename) as data:
            for line in data:
                if line.startswith('#'):
                    continue
                records.append(map(float, line.split()))
        self.timeseries = numpy.array(records).T
        try:
            self.qavg = numpy.loadtxt(self.outarray)
        except IOError, err:
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

        if not filename is None:
            savefig(filename)

    def _plot_qavg_pcolor(self, filename=None, **kwargs):
        """Plot :attr:`ContactAnalysis1.qavg`, the matrix of average native contacts."""
        from pylab import pcolor, gca, meshgrid, xlabel, ylabel, xlim, ylim, colorbar, savefig

        x,y = self.selections[0].resids(), self.selections[1].resids()
        X,Y = meshgrid(x,y)

        pcolor(X,Y,self.qavg.T, **kwargs)
        gca().set_aspect('equal')

        xlim(min(x),max(x))
        ylim(min(y),max(y))

        xlabel("residues")
        ylabel("residues")

        colorbar()

        if not filename is None:
            savefig(filename)


    def plot_qavg(self, filename=None, **kwargs):
        """Plot :attr:`ContactAnalysis1.qavg`, the matrix of average native contacts.

        .. function:: ContactAnalysis1.plot_qavg([filename, ...])

        If *filename* is supplied then the figure is also written to file (the
        suffix determines the file type, e.g. pdf, png, eps, ...). All other
        keyword arguments are passed on to :func:`pylab.imshow`.
        """
        from pylab import imshow, xlabel, ylabel, xlim, ylim, colorbar, cm, clf, savefig

        x,y = self.selections[0].resids(), self.selections[1].resids()

        kwargs['origin'] = 'lower'
        kwargs.setdefault('aspect', 'equal')
        kwargs.setdefault('interpolation', 'nearest')
        kwargs.setdefault('vmin', 0)
        kwargs.setdefault('vmax', 1)
        kwargs.setdefault('cmap', cm.hot)
        kwargs.setdefault('extent', (min(x), max(x), min(y), max(y)))

        clf()
        imshow(self.qavg.T, **kwargs)

        xlim(min(x),max(x))
        ylim(min(y),max(y))

        xlabel("residue from %r" % self.selection_strings[0])
        ylabel("residue from %r" % self.selection_strings[1])

        colorbar()

        if not filename is None:
            savefig(filename)
