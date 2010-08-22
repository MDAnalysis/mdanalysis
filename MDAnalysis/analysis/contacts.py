# native contact analysis
# Part of MDAnalysis http://mdanalysis.googlecode.com
# Copyright (c) 2006-2010 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# Released under the GNU Public Licence, v2
"""
:Author: Oliver Beckstein
:Year: 2010
:Copyright: GNU Public License v3

Contact analysis ("q1-q2")
==========================

See http://lorentz.dynstr.pasteur.fr/joel/adenylate.php for an example of
contact analysis applied to MinActionPath trajectories of AdK (although this
was *not* performed with MDAnalysis --- it is just to give an idea what it is
about).

Example
-------

Analyze a single DIMS transition of AdK between its closed and open
conformation and plot the trajectory projected on q1-q2::

  import MDAnalysis.analysis.contacts
  from MDAnalysis.tests.datafiles import *
  C = MDAnalysis.analysis.contacts.ContactAnalysis(PSF, DCD)
  C.run()
  C.plot()

"""

import os
import warnings
import bz2
import numpy
import MDAnalysis
from MDAnalysis.core.distances import self_distance_array


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

        r1 = MDAnalysis.Universe(topology, pdbfilename=self.ref1)
        r2 = MDAnalysis.Universe(topology, pdbfilename=self.ref2)
    
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
        self.qtmp = numpy.zeros_like(self.q)  # pre-allocated array

    def output_exists(self, force=False):
        """Return True if default output file already exists.

        Disable with force=True (will always return False)
        """
        return (os.path.isfile(self.output) or os.path.isfile(self.output_bz2)) \
               and not (self.force or force)

    def run(self, store=True, force=False):
        """Analyze trajectory and produce timeseries.

        Stores results in :attr:`ContactAnalysis.timeseries` (if
        store=True) and write it to a bzip2-compressed data file.
        """
        if self._skip or self.output_exists(force=force):
            print "File %(output)r or %(output_bz2)r already exists, skipping %(trajectory)r." % vars(self)
            return None
        
        try:
            outbz2 = bz2.BZ2File(self.output_bz2, mode='w', buffering=8192)
            outbz2.write("# q1-q2 analysis\n# nref1 = %d\n# nref2 = %d\n" 
                         % (self.nref[0], self.nref[1]))
            outbz2.write("# frame  q1  q2   n1  n2\n")
            records = []
            for ts in self.u.trajectory:
                frame = ts.frame
                # use pre-allocated distance array to save a little bit of time
                self_distance_array(self.ca.coordinates(), result=self.d)
                self.qarray(self.d, out=self.q)
                n1, q1 = self.qN(self.q, 0, out=self.qtmp)
                n2, q2 = self.qN(self.q, 1, out=self.qtmp)

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


