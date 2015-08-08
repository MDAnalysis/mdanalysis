# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
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
Generation and Analysis of X3DNA helicoidal parameter profiles --- :mod:`MDAnalysis.analysis.x3dna`
===================================================================================================

:Author: Elizabeth Denning
:Year: 2013-2014
:Copyright: GNU Public License v2

.. versionadded: 0.8

With the help of this module, X3DNA_ can be run on frames in a trajectory. Data
can be combined and analyzed. X3DNA_ [Lu2003]_ [Lu2008]_ must be installed
separately.


.. rubric:: References

.. [Lu2003] Xiang-Jun Lu & Wilma K. Olson (2003).
            3DNA: a software package for the analysis, rebuilding and visualization
            for three-dimensional nucleic acid structure
            Nucleic Acids Res. 31(17), 5108-21.

.. [Lu2008] Xiang-Jun Lu & Wilma K. Olson (2008).
            3DNA: a versatile, integrated software system for the analysis, rebuilding
            and visualization of three-dimensional nucleic-acid structures.
            Nat Protoc. 3(7), 1213-27.

.. _X3DNA: http://x3dna.org/


Examples
--------

Single structure
~~~~~~~~~~~~~~~~

B-DNA structure::

   from MDAnalysis.analysis.x3dna import X3DNA, X3DNAtraj
   from MDAnalysis.tests.datafiles import PDB_X3DNA

   # set path to your x3dna binary in bashrc file
   H = X3DNA(PDB_X3DNA, executable="x3dna_ensemble analyze -b 355d.bps -p pdbfile")
   H.run()
   H.collect()
   H.plot()


Trajectory
~~~~~~~~~~

Analyzing a trajectory::

  u = MDAnalysis.Universe(psf, trajectory)
  H = X3DNAtraj(u, ...)
  H.run()
  H.plot()
  H.save()

The profiles are available as the attribute :attr:`X3DNAtraj.profiles`
(``H.profiles`` in the example) and are indexed by frame number but
can also be indexed by an arbitrary order parameter as shown in the
next example.



Analysis
--------

.. autoclass:: X3DNA
   :members:
   :inherited-members:

   .. attribute:: profiles

      ``x3dna_ensemble analyze -b 355d.bps -p pdbfile attribute``:
      After running :meth:`X3DNA.collect`, this dict contains all the
      X3DNA profiles, indexed by the frame number. If only a single
      frame was analyzed then this will be ``X3DNA.profiles[0]``. Note
      that the order is random; one needs to sort the keys first.

.. autoclass:: X3DNAtraj
   :members:
   :inherited-members:

   .. attribute:: profiles

      After running :meth:`X3DNA.collect`, this dict contains all the
      X3DNA profiles, indexed by the frame number.



Utilities
---------

.. autoclass:: ApplicationError


"""

import glob
import os
import errno
import shutil
import warnings
import numpy as np
import os.path
import subprocess
import tempfile
import textwrap

from MDAnalysis import ApplicationError
from MDAnalysis.lib.util import which, realpath, asiterable

import logging

logger = logging.getLogger("MDAnalysis.analysis.x3dna")


def mean_std_from_x3dnaPickle(profile):
    """Get mean and standard deviation of helicoidal parameters from a saved *profile*.

    The *profile* should have been saved with :meth:`BaseX3DNA.save`. Then
    load it with ::

      profile = cPickle.load(open("x3dna.pickle"))
      h_mean, h_std = mean_std_from_x3dnaPickle(profile)
    """
    if profile.x3dna_param is False:
        bp_shear, bp_stretch, bp_stagger, bp_rise, bp_shift, bp_slide, bp_buckle, bp_prop, bp_open, bp_tilt, bp_roll,\
            bp_twist = [], [], [], [], [], [], [], [], [], [], [], []
        for i in range(len(profile)):
            bp_shear.append(profile.values()[i].Shear)
            bp_stretch.append(profile.values()[i].Stretch)
            bp_stagger.append(profile.values()[i].Stagger)
            bp_buckle.append(profile.values()[i].Buckle)
            bp_prop.append(profile.values()[i].Propeller)
            bp_open.append(profile.values()[i].Opening)
            bp_rise.append(profile.values()[i].Rise)
            bp_shift.append(profile.values()[i].Shift)
            bp_slide.append(profile.values()[i].Slide)
            bp_tilt.append(profile.values()[i].Tilt)
            bp_roll.append(profile.values()[i].Roll)
            bp_twist.append(profile.values()[i].Twist)
        bp_shear, bp_stretch, bp_stagger, bp_rise, bp_shift, bp_slide, bp_buckle, bp_prop, bp_open, bp_tilt, bp_roll,\
            bp_twist = np.array(bp_shear), np.array(bp_stretch), np.array(bp_stagger), np.array(bp_rise),\
            np.array(bp_shift), np.array(bp_slide), np.array(bp_buckle), np.array(bp_prop), \
            np.array(bp_open), np.array(bp_tilt), np.array(bp_roll), np.array(bp_twist)
        na_avg, na_std = [], []
        for j in range(len(bp_shear[0])):
            na_avg.append([
                np.mean(bp_shear[:, j]), np.mean(bp_stretch[:, j]), np.mean(bp_stagger[:, j]),
                np.mean(bp_buckle[:, j]), np.mean(bp_prop[:, j]), np.mean(bp_open[:, j]),
                np.mean(bp_shift[:, j]), np.mean(bp_slide[:, j]), np.mean(bp_rise[:, j]),
                np.mean(bp_tilt[:, j]), np.mean(bp_roll[:, j]), np.mean(bp_twist[:, j])])
            na_std.append([
                np.std(bp_shear[:, j]), np.std(bp_stretch[:, j]), np.std(bp_stagger[:, j]),
                np.std(bp_buckle[:, j]), np.std(bp_prop[:, j]), np.std(bp_open[:, j]),
                np.std(bp_shift[:, j]), np.std(bp_slide[:, j]), np.std(bp_rise[:, j]),
                np.std(bp_tilt[:, j]), np.std(bp_roll[:, j]), np.std(bp_twist[:, j])])
    else:
        bp_rise, bp_shift, bp_slide, bp_tilt, bp_roll, bp_twist = [], [], [], [], [], [], [], [], [], [], [], []
        for i in range(len(profile)):
            #print i
            bp_rise.append(profile.values()[i].Rise)
            bp_shift.append(profile.values()[i].Shift)
            bp_slide.append(profile.values()[i].Slide)
            bp_tilt.append(profile.values()[i].Tilt)
            bp_roll.append(profile.values()[i].Roll)
            bp_twist.append(profile.values()[i].Twist)
        bp_rise, bp_shift, bp_slide, bp_tilt, bp_roll, bp_twist = np.array(bp_shear),np.array(bp_stretch),\
            np.array(bp_stagger), np.array(bp_rise), np.array(bp_shift), np.array(bp_slide),\
            np.array(bp_buckle), np.array(bp_prop), np.array(bp_open), np.array(bp_tilt),\
            np.array(bp_roll), np.array(bp_twist)
        na_avg, na_std = [], []
        for j in range(len(bp_shift[0])):
            na_avg.append([
                np.mean(bp_shift[:, j]), np.mean(bp_slide[:, j]), np.mean(bp_rise[:, j]),
                np.mean(bp_tilt[:, j]), np.mean(bp_roll[:, j]), np.mean(bp_twist[:, j])])
            na_std.append([
                np.std(bp_shift[:, j]), np.std(bp_slide[:, j]), np.std(bp_rise[:, j]),
                np.std(bp_tilt[:, j]), np.std(bp_roll[:, j]), np.std(bp_twist[:, j])])

    na_avg, na_std = np.array(na_avg), np.array(na_std)
    return na_avg, na_std


class BaseX3DNA(object):
    """Baseclass for X3DNA_ analysis, providing plotting and utility functions.


    .. _X3DNA: http://x3dna.org
    """

    def save(self, filename="x3dna.pickle"):
        """Save :attr:`profiles` as a Python pickle file *filename*.

        Load profiles dictionary with ::

           import cPickle
           profiles = cPickle.load(open(filename))

        """
        import cPickle

        cPickle.dump(self.profiles, open(filename, "wb"), cPickle.HIGHEST_PROTOCOL)

    def mean_std(self):
        """H.mean_std() returns the mean and standard deviation of base parameters (order of base parameters:Shear,
        Stretch,Stagger,Buckle,Propeller,Opening,Shift,Slide,Rise,Tilt,Roll,Twist) for each nucleic acid pair.

        """

        bp_shear, bp_stretch, bp_stagger, bp_rise, bp_shift, bp_slide, bp_buckle, bp_prop, bp_open, bp_tilt, bp_roll,\
            bp_twist = [], [], [], [], [], [], [], [], [], [], [], []
        for i in range(len(self.profiles)):
            bp_shear.append(self.profiles.values()[i].Shear)
            bp_stretch.append(self.profiles.values()[i].Stretch)
            bp_stagger.append(self.profiles.values()[i].Stagger)
            bp_buckle.append(self.profiles.values()[i].Buckle)
            bp_prop.append(self.profiles.values()[i].Propeller)
            bp_open.append(self.profiles.values()[i].Opening)
            bp_rise.append(self.profiles.values()[i].Rise)
            bp_shift.append(self.profiles.values()[i].Shift)
            bp_slide.append(self.profiles.values()[i].Slide)
            bp_tilt.append(self.profiles.values()[i].Tilt)
            bp_roll.append(self.profiles.values()[i].Roll)
            bp_twist.append(self.profiles.values()[i].Twist)
        bp_shear, bp_stretch, bp_stagger, bp_rise, bp_shift, bp_slide, bp_buckle, bp_prop, bp_open, bp_tilt, bp_roll,\
            bp_twist = np.array(bp_shear), np.array(bp_stretch), np.array(bp_stagger), np.array(bp_rise),\
            np.array(bp_shift), np.array(bp_slide), np.array(bp_buckle), np.array(bp_prop),\
            np.array(bp_open), np.array(bp_tilt), np.array(bp_roll), np.array(bp_twist)
        na_avg, na_std = [], []
        for j in range(len(bp_shear[0])):
            na_avg.append([
                np.mean(bp_shear[:, j]), np.mean(bp_stretch[:, j]), np.mean(bp_stagger[:, j]),
                np.mean(bp_buckle[:, j]), np.mean(bp_prop[:, j]), np.mean(bp_open[:, j]),
                np.mean(bp_shift[:, j]), np.mean(bp_slide[:, j]), np.mean(bp_rise[:, j]),
                np.mean(bp_tilt[:, j]), np.mean(bp_roll[:, j]), np.mean(bp_twist[:, j])])
            na_std.append([
                np.std(bp_shear[:, j]), np.std(bp_stretch[:, j]), np.std(bp_stagger[:, j]),
                np.std(bp_buckle[:, j]), np.std(bp_prop[:, j]), np.std(bp_open[:, j]),
                np.std(bp_shift[:, j]), np.std(bp_slide[:, j]), np.std(bp_rise[:, j]),
                np.std(bp_tilt[:, j]), np.std(bp_roll[:, j]), np.std(bp_twist[:, j])])
        na_avg, na_std = np.array(na_avg), np.array(na_std)
        return na_avg, na_std

    def mean(self):
        """H.mean() returns the mean value for the base parameters (order of base parameters:Shear,Stretch,Stagger,
        Buckle,Propeller,Opening,Shift,Slide,Rise,Tilt,Roll,Twist) for each nucleic acid pair.

        """
        bp_shear, bp_stretch, bp_stagger, bp_rise, bp_shift, bp_slide, bp_buckle, bp_prop, bp_open, bp_tilt, bp_roll,\
            bp_twist = [], [], [], [], [], [], [], [], [], [], [], []
        for i in range(len(self.profiles)):
            bp_shear.append(self.profiles.values()[i].Shear)
            bp_stretch.append(self.profiles.values()[i].Stretch)
            bp_stagger.append(self.profiles.values()[i].Stagger)
            bp_buckle.append(self.profiles.values()[i].Buckle)
            bp_prop.append(self.profiles.values()[i].Propeller)
            bp_open.append(self.profiles.values()[i].Opening)
            bp_rise.append(self.profiles.values()[i].Rise)
            bp_shift.append(self.profiles.values()[i].Shift)
            bp_slide.append(self.profiles.values()[i].Slide)
            bp_tilt.append(self.profiles.values()[i].Tilt)
            bp_roll.append(self.profiles.values()[i].Roll)
            bp_twist.append(self.profiles.values()[i].Twist)
        bp_shear, bp_stretch, bp_stagger, bp_rise, bp_shift, bp_slide, bp_buckle, bp_prop, bp_open, bp_tilt, bp_roll,\
            bp_twist = np.array(bp_shear), np.array(bp_stretch), np.array(bp_stagger), np.array(bp_rise),\
            np.array(bp_shift), np.array(bp_slide), np.array(bp_buckle), np.array(bp_prop),\
            np.array(bp_open), np.array(bp_tilt), np.array(bp_roll), np.array(bp_twist)

        na_avg = []
        for j in range(len(bp_shear[0])):
            na_avg.append([
                np.mean(bp_shear[:, j]), np.mean(bp_stretch[:, j]), np.mean(bp_stagger[:, j]),
                np.mean(bp_buckle[:j]), np.mean(bp_prop[:, j]), np.mean(bp_open[:, j]),
                np.mean(bp_shift[:, j]), np.mean(bp_slide[:, j]), np.mean(bp_rise[:, j]),
                np.mean(bp_tilt[:, j]), np.mean(bp_roll[:, j]), np.mean(bp_twist[:, j])])
        na_avg = np.array(na_avg)
        return na_avg

    def std(self):
        """H.std() returns the standard deviation of base parameters (order of base parameters:Shear,Stretch,Stagger,
        Buckle,Propeller,Opening,Shift,Slide,Rise,Tilt,Roll,Twist) for each nucleic acid pair.

        """

        bp_shear, bp_stretch, bp_stagger, bp_rise, bp_shift, bp_slide, bp_buckle, bp_prop, bp_open, bp_tilt, bp_roll,\
            bp_twist = [], [], [], [], [], [], [], [], [], [], [], []
        for i in range(len(self.profiles)):
            bp_shear.append(self.profiles.values()[i].Shear)
            bp_stretch.append(self.profiles.values()[i].Stretch)
            bp_stagger.append(self.profiles.values()[i].Stagger)
            bp_buckle.append(self.profiles.values()[i].Buckle)
            bp_prop.append(self.profiles.values()[i].Propeller)
            bp_open.append(self.profiles.values()[i].Opening)
            bp_rise.append(self.profiles.values()[i].Rise)
            bp_shift.append(self.profiles.values()[i].Shift)
            bp_slide.append(self.profiles.values()[i].Slide)
            bp_tilt.append(self.profiles.values()[i].Tilt)
            bp_roll.append(self.profiles.values()[i].Roll)
            bp_twist.append(self.profiles.values()[i].Twist)
        bp_shear, bp_stretch, bp_stagger, bp_rise, bp_shift, bp_slide, bp_buckle, bp_prop, bp_open, bp_tilt, bp_roll,\
            bp_twist = np.array(bp_shear), np.array(bp_stretch), np.array(bp_stagger), np.array(bp_rise),\
            np.array(bp_shift), np.array(bp_slide), np.array(bp_buckle), np.array(bp_prop),\
            np.array(bp_open), np.array(bp_tilt), np.array(bp_roll), np.array(bp_twist)

        na_std = []
        for j in range(len(bp_shear[0])):
            na_std.append([
                np.std(bp_shear[:, j]), np.std(bp_stretch[:, j]), np.std(bp_stagger[:, j]),
                np.std(bp_buckle[:j]), np.std(bp_prop[:, j]), np.std(bp_open[:, j]), np.std(bp_shift[:, j]),
                np.std(bp_slide[:, j]), np.std(bp_rise[:, j]), np.std(bp_tilt[:, j]), np.std(bp_roll[:, j]),
                np.std(bp_twist[:, j])])
        na_std = np.array(na_std)
        return na_std

    def plot(self, **kwargs):
        """Plot base parameter profiles  in a 1D graph
           (plotting of the mean and standard deviation parameter value for each individual base pair).

        :Keywords:
           *Not available at this time*

        """
        import matplotlib.pyplot as plt
        from itertools import izip

        na_avg, na_std = self.mean_std()
        for k in range(len(na_avg[0])):
            ax = kwargs.pop('ax', plt.subplot(111))
            x = range(1, len(na_avg[:, k]) + 1)
            ax.errorbar(x, na_avg[:, k], yerr=na_std[:, k], fmt='-o')
            ax.set_xlim(0, len(na_avg[:, k]) + 1)
            ax.set_xlabel(r"Nucleic Acid Number")
            param = self.profiles.values()[0].dtype.names[k]
            if param in ["Shear", "Stretch", "Stagger", "Rise", "Shift", "Slide"]:
                ax.set_ylabel("%s ($\AA$)" % (param))
            else:
                ax.set_ylabel("%s (deg)" % (param))
            ax.figure.savefig("%s.png" % (param))
            ax.figure.clf()

    def sorted_profiles_iter(self):
        """Return an iterator over profiles sorted by frame/order parameter *q*.

        The iterator produces tuples ``(q, profile)``.
        """
        if self.profiles is None:
            raise StopIteration
        for q in sorted(self.profiles):
            yield (q, self.profiles[q])

    __iter__ = sorted_profiles_iter


class X3DNA(BaseX3DNA):
    """Run X3DNA_ on a single frame or a DCD trajectory.

    Only a subset of all X3DNA control parameters is supported and can be set
    with keyword arguments. For further details on X3DNA_ see the `X3DNA docs`_.

    Running X3DNA with the :class:`X3DNA` class is a 3-step process:

     1. set up the class with all desired parameters
     2. run X3DNA with :meth:`X3DNA.run`
     3. collect the data from the output file with :meth:`X3DNA.collect`

    The class also provides some simple plotting functions of the collected
    data such as :meth:`X3DNA.plot` or :meth:`X3DNA.plot3D`.

    .. versionadded:: 0.8

    .. _`X3DNA docs`: http://forum.x3dna.org/
    """

    def __init__(self, filename, **kwargs):
        """Set up parameters to run X3DNA_ on PDB *filename*.

        :Arguments:

          *filename*

               The *filename* is used as input for X3DNA in the "xdna_ensemble"
               command.  It specifies the name of a PDB coordinate file
               to be used. This must be in Brookhaven protein databank format
               or something closely approximating this. Both ATOM and HETATM
               records are read. Note that if water molecules or ions are
               present in the channel these can be ignored on read by the use
               of the *ignore_residues* keyword.

        :Keywords:

          *executable*

               Path to the :program:`xdna_ensemble` executable directories
               (e.g. ``/opt/x3dna/2.1 and /opt/x3dna/2.1/bin``) must be set
               and then added to export in bashrc file. See X3DNA
               documentation for set-up instructions.

          *x3dna_param*

               Determines whether base step or base pair parameters will be
               calculated. If True then stacked base step parameters will be
               analyzed [Default is True].  If False then stacked base pair
               parameters will be analyzed.
        """
        # list of temporary files, to be cleaned up on __del__
        self.tempfiles = [
            "auxiliary.par", "bestpairs.pdb", "bp_order.dat", "bp_helical.par", "cf_7methods.par",
            "col_chains.scr", "col_helices.scr", "hel_regions.pdb", "ref_frames.dat", "hstacking.pdb", "stacking.pdb"
        ]
        self.tempdirs = []
        self.filename = filename

        logger.info("Setting up X3DNA analysis for %(filename)r", vars(self))

        # guess executables
        self.exe = {}
        x3dna_exe_name = kwargs.pop('executable', 'xdna_ensemble')
        self.x3dna_param = kwargs.pop('x3dna_param', True)
        self.exe['xdna_ensemble'] = which(x3dna_exe_name)
        if self.exe['xdna_ensemble'] is None:
            errmsg = "X3DNA binary %(x3dna_exe_name)r not found." % vars()
            logger.fatal(errmsg)
            logger.fatal("%(x3dna_exe_name)r must be on the PATH or provided as keyword argument 'executable'.",
                         vars())
            raise OSError(errno.ENOENT, errmsg)
        x3dnapath = os.path.dirname(self.exe['xdna_ensemble'])
        self.logfile = kwargs.pop("logfile", "bp_step.par")

        if self.x3dna_param is False:
            self.template = textwrap.dedent("""x3dna_ensemble analyze -b 355d.bps --one %(filename)r """)
        else:
            self.template = textwrap.dedent("""find_pair -s %(filename)r stdout |analyze stdin """)

        # sanity checks
        for program, path in self.exe.items():
            if path is None or which(path) is None:
                logger.error("Executable %(program)r not found, should have been %(path)r.",
                             vars())
        # results
        self.profiles = {}

    def run(self, **kwargs):
        """Run X3DNA on the input file."""
        inpname = kwargs.pop("inpfile", None)
        outname = kwargs.pop("outfile", self.logfile)

        x3dnaargs = vars(self).copy()
        x3dnaargs.update(kwargs)
        x3dna_param = kwargs.pop('x3dna_param', self.x3dna_param)

        inp = self.template % x3dnaargs
        if inpname:
            with open(inpname, "w") as f:
                f.write(inp)
            logger.debug("Wrote X3DNA input file %r for inspection", inpname)

        logger.info("Starting X3DNA on %(filename)r (trajectory: %(dcd)r)", x3dnaargs)
        logger.debug("%s", self.exe['xdna_ensemble'])
        with open(outname, "w") as output:
            x3dna = subprocess.call([inp], shell=True)
        with open(outname, "r") as output:
            # X3DNA is not very good at setting returncodes so check ourselves
            for line in output:
                if line.strip().startswith(('*** ERROR ***', 'ERROR')):
                    x3dna.returncode = 255
                    break
        if x3dna.bit_length != 0:
            logger.fatal("X3DNA Failure (%d). Check output %r", x3dna.bit_length, outname)
        logger.info("X3DNA finished: output file %(outname)r", vars())

    def collect(self, **kwargs):
        """Parse the output from a X3DNA run into numpy recarrays.

        Can deal with outputs containing multiple frames. Output format:

        The method saves the result as :attr:`X3DNA.profiles`, a dictionary
        indexed by the frame number. Each entry is a
        :class:`np.recarray`.

        If the keyword *outdir* is supplied (e.g. ".") then each profile is
        saved to a gzipped data file.

        :Keywords:
           *run*
              identifier, free form [1]
           *outdir*
              save output data under *outdir*/*run* if set to any other
              value but ``None`` [``None``]

        """
        #        Shear    Stretch   Stagger   Buckle   Prop-Tw   Opening     Shift     Slide     Rise      Tilt
        # Roll      Twist
        #0123456789.0123456789.0123456789.0123456789.0123456789.0123456789.
        #            11          22          33          44
        #123456789.123456789.123456789.123456789.123456789.123456789.123456789.
        #1           13
        #T-A     -0.033    -0.176     0.158   -12.177    -8.979     1.440     0.000     0.000     0.000     0.000
        #  0.000     0.000
        #C-G     -0.529     0.122    -0.002    -7.983   -10.083    -0.091    -0.911     1.375     3.213    -0.766
        # -4.065    41.492
        # only parse bp_step.par
        x3dna_output = kwargs.pop("x3dnaout", self.logfile)

        run = kwargs.pop("run", 1)  # id number
        outdir = kwargs.pop("outdir", os.path.curdir)

        logger.info("Collecting X3DNA profiles for run with id %s", run)
        length = 1  # length of trajectory --- is this really needed?? No... just for info
        if '*' in self.filename:
            import glob

            filenames = glob.glob(self.filename)
            length = len(filenames)
            if length == 0:
                logger.error("Glob pattern %r did not find any files.", self.filename)
                raise ValueError("Glob pattern %r did not find any files." % (self.filename,))
            logger.info("Found %d input files based on glob pattern %s", length, self.filename)

        # one recarray for each frame, indexed by frame number
        # TODO: use ordered dict
        self.profiles = {}

        logger.info("Run %s: Reading %d X3DNA profiles from %r", run, length, x3dna_output)
        x3dna_profile_no = 0
        records = []
        with open(x3dna_output, "r") as x3dna:
            read_data = False
            for line in x3dna:
                line = line.rstrip()  # preserve columns (FORTRAN output...)
                if self.x3dna_param is False:
                    if line.startswith("#        Shear"):
                        read_data = True
                        logger.debug("Started reading data")
                        fields = line.split()
                        x3dna_profile_no = int(1)  # useless int value code based off hole plugin
                        records = []
                        continue
                    if read_data:
                        if len(line.strip()) != 0:
                            try:
                                Sequence, Shear, Stretch, Stagger, Buckle, Propeller, Opening, Shift, Slide, Rise, \
                                    Tilt, Roll, Twist = line.split()
                            except:
                                logger.critical("Run %d: Problem parsing line %r", run, line.strip())
                                logger.exception("Check input file %r.", x3dna_output)
                                raise
                            records.append(
                                [float(Shear), float(Stretch), float(Stagger), float(Buckle), float(Propeller),
                                    float(Opening), float(Shift), float(Slide), float(Rise), float(Tilt), float(Roll),
                                    float(Twist)])
                            continue
                        else:
                            # end of records (empty line)
                            read_data = False
                else:
                    if line.startswith("#      Shift"):
                        read_data = True
                        logger.debug("Started reading data")
                        fields = line.split()
                        x3dna_profile_no = int(1)  # useless int value code based off hole plugin
                        records = []
                        continue
                    if read_data:
                        if len(line.strip()) != 0:
                            try:
                                Sequence, Shift, Slide, Rise, Tilt, Roll, Twist = line.split()
                            except:
                                logger.critical("Run %d: Problem parsing line %r", run, line.strip())
                                logger.exception("Check input file %r.", x3dna_output)
                                raise
                            records.append(
                                [float(Shift), float(Slide), float(Rise), float(Tilt), float(Roll), float(Twist)])
                            continue
                        else:
                            # end of records (empty line)
                            read_data = False
            if self.x3dna_param is False:
                frame_x3dna_output = np.rec.fromrecords(records, formats="f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8",
                                                           names="Shear,Stretch,Stagger,Buckle,Propeller,Opening,"
                                                                 "Shift,Slide,Rise,Tilt,Roll,Twist")
            else:
                frame_x3dna_output = np.rec.fromrecords(records, formats="f8,f8,f8,f8,f8,f8",
                                                           names="Shift,Slide,Rise,Tilt,Roll,Twist")
            # store the profile
            self.profiles[x3dna_profile_no] = frame_x3dna_output
            logger.debug("Collected X3DNA profile for frame %d (%d datapoints)",
                         x3dna_profile_no, len(frame_x3dna_output))
            # save a profile for each frame (for debugging and scripted processing)
            # a tmp folder for each trajectory
            if outdir is not None:
                rundir = os.path.join(outdir, "run_" + str(run))
                os.system("rm -f tmp*.out")
                if not os.path.exists(rundir):
                    os.makedirs(rundir)
                frame_x3dna_txt = os.path.join(rundir, "bp_step_%s_%04d.dat.gz" % (run, x3dna_profile_no))
                np.savetxt(frame_x3dna_txt, frame_x3dna_output)
                logger.debug("Finished with frame %d, saved as %r", x3dna_profile_no, frame_x3dna_txt)
                # if we get here then we haven't found anything interesting
        if len(self.profiles) == length:
            logger.info("Collected X3DNA profiles for %d frames", len(self.profiles))
        else:
            logger.warn("Missing data: Found %d X3DNA profiles from %d frames.", len(self.profiles), length)

    def __del__(self):
        for f in self.tempfiles:
            try:
                os.unlink(f)
            except OSError:
                pass
        for d in self.tempdirs:
            shutil.rmtree(d, ignore_errors=True)


class X3DNAtraj(BaseX3DNA):
    """Analyze all frames in a trajectory.

    The :class:`X3DNA` class provides a direct interface to X3DNA. X3DNA itself
    has limited support for analysing trajectories but cannot deal with all the
    trajectory formats understood by MDAnalysis. This class can take any
    universe and feed it to X3DNA. By default it sequentially creates a PDB for
    each frame and runs X3DNA on the frame.
    """

    def __init__(self, universe, **kwargs):
        """Set up the class.

        :Arguments:

          *universe*
               The input trajectory as part of a
               :class:`~MDAnalysis.core.AtomGroup.Universe`. trajectory is
               converted to a sequence of PDB files and X3DNA is run on each
               individual file. (Use the *start*, *stop*, and *step* keywords
               to slice the trajectory.)

          *start*, *stop*, *step*
               frame indices to slice the trajectory as
               ``universe.trajectory[start, stop, step]``

          *x3dna_param*
               indicates whether stacked bases or stacked base-pairs will be
               analyzed. True is bases or False is stacked base-pairs [Default
               is True].

          *kwargs*
               All other keywords are passed on to :class:`X3DNA` (see there for description).

        """
        self.universe = universe
        self.selection = kwargs.pop("selection", "nucleic")

        self.x3dna_param = kwargs.pop('x3dna_param', True)
        self.start = kwargs.pop('start', None)
        self.stop = kwargs.pop('stop', None)
        self.step = kwargs.pop('step', None)

        self.x3dna_kwargs = kwargs

        # processing

    def run(self, **kwargs):
        """Run X3DNA on the whole trajectory and collect profiles.

        Keyword arguments *start*, *stop*, and *step* can be used to only
        analyse part of the trajectory.
        """
        from itertools import izip

        start = kwargs.pop('start', self.start)
        stop = kwargs.pop('stop', self.stop)
        step = kwargs.pop('step', self.step)
        x3dna_param = kwargs.pop('x3dna_param', self.x3dna_param)
        x3dna_kw = self.x3dna_kwargs.copy()
        x3dna_kw.update(kwargs)

        profiles = {}

        nucleic = self.universe.select_atoms(self.selection)
        for ts in self.universe.trajectory[start:stop:step]:
            print ts.frame
            print nucleic.atoms
            logger.info("X3DNA analysis frame %4d ", ts.frame)
            fd, pdbfile = tempfile.mkstemp(suffix=".pdb")
            os.close(fd)
            nucleic.write(pdbfile)
            os.system("find_pair %s 355d.bps" % pdbfile)
            try:
                nucleic.write(pdbfile)
                '''
                if int(ts.frame) in [int(start+1), 1]:
                    #print start
                    os.system("find_pair %s 355d.bps" %(pdbfile))
                    #print "complete"
                '''
                x3dna_profiles = self.run_x3dna(pdbfile, **x3dna_kw)
            finally:
                try:
                    os.unlink(pdbfile)
                except OSError:
                    pass
            if len(x3dna_profiles) != 1:
                err_msg = "Got {0} profiles ({1}) --- should be 1 (time step {2})".format(
                    len(x3dna_profiles), x3dna_profiles.keys(), ts)
                logger.error(err_msg)
                warnings.warn(err_msg)
            profiles[ts.frame] = x3dna_profiles.values()[0]
        self.profiles = profiles

    def run_x3dna(self, pdbfile, **kwargs):
        """Run X3DNA on a single PDB file *pdbfile*."""
        kwargs['x3dna_param'] = self.x3dna_param
        H = X3DNA(pdbfile, **kwargs)
        H.run()
        H.collect()
        return H.profiles
