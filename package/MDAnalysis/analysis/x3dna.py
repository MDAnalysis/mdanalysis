# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; encoding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2012 Naveen Michaud-Agrawal,
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
Generation and Analysis of X3DNA helicoidal parameter profiles --- :mod:`MDAnalysis.analysis.x3dna`
=================================================================================

:Author: Elizabeth Denning 
:Year: 2013-2014
:Copyright: GNU Public License v2

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

.. _X3DNA: http://http://x3dna.org/


Examples
--------

Single structure
~~~~~~~~~~~~~~~~

B-DNA structure::

   from MDAnalysis.analysis.x3dna import X3DNA, X3DNAtraj
   from MDAnalysis.tests.datafiles import PDB_X3DNA

   H = X3DNA(PDB_X3DNA, executable="find_pair pdbfile stdout | analyze stdin ")  # set path to your x3dna binary in bashrc file
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

The profiles are available as the attribute :attr:`X3DNAtraj.profiles`
(``H.profiles`` in the example) and are indexed by frame number but
can also be indexed by an arbitrary order parameter as shown in the
next example.



Analysis
--------

.. autoclass:: X3DNA
   :members:
   :inherited-members:


   find_pair pdbfile stdout | analyze stdin ") attribute:: profiles
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
from __future__ import with_statement

import numpy

import matplotlib
import pylab
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
from MDAnalysis.core.util import which, realpath, asiterable

import logging
logger = logging.getLogger("MDAnalysis.analysis.x3dna")


class BaseX3DNA(object):
    """Baseclass for X3DNA analysis, providing plotting and utility functions"""

    def save(self, filename="x3dna.pickle"):
        """Save :attr:`profiles` as a Python pickle file *filename*.

        Load profiles dictionary with ::

           import cPickle
           profiles = cPickle.load(open(filename))

        """
        import cPickle
        cPickle.dump(self.profiles, open(filename, "wb"), cPickle.HIGHEST_PROTOCOL)

    def mean_std(self):
        """H.mean_std() returns the mean and standard deviation of base parameters (order of base parameters:Shear,Stretch,Stagger,Buckle,Propeller,Opening,Shift,Slide,Rise,Tilt,Roll,Twist) for each nucleic acid pair. 

        """

        bp_shear,bp_stretch,bp_stagger,bp_rise,bp_shift,bp_slide,bp_buckle,bp_prop,bp_open,bp_tilt,bp_roll,bp_twist = [],[],[],[],[],[],[],[],[],[],[],[]
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
        bp_shear,bp_stretch,bp_stagger,bp_rise,bp_shift,bp_slide,bp_buckle,bp_prop,bp_open,bp_tilt,bp_roll,bp_twist = numpy.array(bp_shear),numpy.array(bp_stretch),numpy.array(bp_stagger),numpy.array(bp_rise),numpy.array(bp_shift),numpy.array(bp_slide),numpy.array(bp_buckle),numpy.array(bp_prop),numpy.array(bp_open),numpy.array(bp_tilt),numpy.array(bp_roll),numpy.array(bp_twist)
        na_avg,na_std = [], []
        for j in range(len(bp_shear[0])):
            na_avg.append([numpy.mean(bp_shear[:,j]),numpy.mean(bp_stretch[:,j]),numpy.mean(bp_stagger[:,j]),numpy.mean(bp_buckle[:j]),numpy.mean(bp_prop[:,j]),numpy.mean(bp_open[:,j]),numpy.mean(bp_shift[:,j]),numpy.mean(bp_slide[:,j]),numpy.mean(bp_rise[:,j]),numpy.mean(bp_tilt[:,j]),numpy.mean(bp_roll[:,j]),numpy.mean(bp_twist[:,j])])
            na_std.append([numpy.std(bp_shear[:,j]),numpy.std(bp_stretch[:,j]),numpy.std(bp_stagger[:,j]),numpy.std(bp_buckle[:j]),numpy.std(bp_prop[:,j]),numpy.std(bp_open[:,j]),numpy.std(bp_shift[:,j]),numpy.std(bp_slide[:,j]),numpy.std(bp_rise[:,j]),numpy.std(bp_tilt[:,j]),numpy.std(bp_roll[:,j]),numpy.std(bp_twist[:,j])])
        na_avg,na_std = numpy.array(na_avg),numpy.array(na_std)
        return na_avg,na_std

    def mean(self):
        """H.mean() returns the mean value for the base parameters (order of base parameters:Shear,Stretch,Stagger,Buckle,Propeller,Opening,Shift,Slide,Rise,Tilt,Roll,Twist) for each nucleic acid pair.

        """
        bp_shear,bp_stretch,bp_stagger,bp_rise,bp_shift,bp_slide,bp_buckle,bp_prop,bp_open,bp_tilt,bp_roll,bp_twist = [],[],[],[],[],[],[],[],[],[],[],[]
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
        bp_shear,bp_stretch,bp_stagger,bp_rise,bp_shift,bp_slide,bp_buckle,bp_prop,bp_open,bp_tilt,bp_roll,bp_twist = numpy.array(bp_shear),numpy.array(bp_stretch),numpy.array(bp_stagger),numpy.array(bp_rise),numpy.array(bp_shift),numpy.array(bp_slide),numpy.array(bp_buckle),numpy.array(bp_prop),numpy.array(bp_open),numpy.array(bp_tilt),numpy.array(bp_roll),numpy.array(bp_twist)

        na_avg = []
        for j in range(len(bp_shear[0])):
            na_avg.append([numpy.mean(bp_shear[:,j]),numpy.mean(bp_stretch[:,j]),numpy.mean(bp_stagger[:,j]),numpy.mean(bp_buckle[:j]),numpy.mean(bp_prop[:,j]),numpy.mean(bp_open[:,j]),numpy.mean(bp_shift[:,j]),numpy.mean(bp_slide[:,j]),numpy.mean(bp_rise[:,j]),numpy.mean(bp_tilt[:,j]),numpy.mean(bp_roll[:,j]),numpy.mean(bp_twist[:,j])])
        na_avg = numpy.array(na_avg)
        return na_avg

    def std(self):
        """H.std() returns the standard deviation of base parameters (order of base parameters:Shear,Stretch,Stagger,Buckle,Propeller,Opening,Shift,Slide,Rise,Tilt,Roll,Twist) for each nucleic acid pair.

        """

        bp_shear,bp_stretch,bp_stagger,bp_rise,bp_shift,bp_slide,bp_buckle,bp_prop,bp_open,bp_tilt,bp_roll,bp_twist = [],[],[],[],[],[],[],[],[],[],[],[]
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
        bp_shear,bp_stretch,bp_stagger,bp_rise,bp_shift,bp_slide,bp_buckle,bp_prop,bp_open,bp_tilt,bp_roll,bp_twist = numpy.array(bp_shear),numpy.array(bp_stretch),numpy.array(bp_stagger),numpy.array(bp_rise),numpy.array(bp_shift),numpy.array(bp_slide),numpy.array(bp_buckle),numpy.array(bp_prop),numpy.array(bp_open),numpy.array(bp_tilt),numpy.array(bp_roll),numpy.array(bp_twist)

        na_std = []
        for j in range(len(bp_shear[0])):
            na_std.append([numpy.std(bp_shear[:,j]),numpy.std(bp_stretch[:,j]),numpy.std(bp_stagger[:,j]),numpy.std(bp_buckle[:j]),numpy.std(bp_prop[:,j]),numpy.std(bp_open[:,j]),numpy.std(bp_shift[:,j]),numpy.std(bp_slide[:,j]),numpy.std(bp_rise[:,j]),numpy.std(bp_tilt[:,j]),numpy.std(bp_roll[:,j]),numpy.std(bp_twist[:,j])])
        na_std = numpy.array(na_std)
        return na_std

    def plot(self, **kwargs):
        """Plot base parameter profiles  in a 1D graph 
           (plotting of the mean and standard deviation parameter value for each individual base pair).

        :Keywords:
           *Not available at this time*

        """
        import matplotlib.pyplot as plt
        from itertools import izip

        #kw, kwargs = self._process_plot_kwargs(kwargs)

        #ax = kwargs.pop('ax', plt.subplot(111))

        na_avg,na_std = self.mean_std()

        for k in range(len(na_avg[0])):
            ax = kwargs.pop('ax', plt.subplot(111))
            x = range(1,len(na_avg[:,k])+1)
            ax.errorbar(x,na_avg[:,k],yerr=na_std[:,k],fmt='-o')
            ax.set_xlim(0,len(na_avg[:,k])+1)
            ax.set_xlabel(r"Nucleic Acid Number")
            param = self.profiles.values()[0].dtype.names[k]
            print param
            if param in ["Shear","Stretch","Stagger","Rise","Shift","Slide"]:
                ax.set_ylabel("%s ($\AA$)"%(param))
            else:
                ax.set_ylabel("%s (deg)"%(param))
            ax.figure.savefig("%s.png"%(param))
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
    """Run X3DNA on a single frame or a DCD trajectory.

    Only a subset of all `X3DNA control parameters`_ is supported and can be set
    with keyword arguments.

    Running X3DNA with the :class:`X3DNA` class is a 3-step process:

     1. set up the class with all desired parameters
     2. run X3DNA with :meth:`X3DNA.run`
     3. collect the data from the output file with :meth:`X3DNA.collect`

    The class also provides some simple plotting functions of the collected
    data such as :meth:`X3DNA.plot` or :meth:`X3DNA.plot3D`.

    .. versionadded:: 0.7.7

    .. _`X3DNA docs`:
       http://http://forum.x3dna.org/ 
    """

    #: List of residues that are ignore by default. Can be changed with
    #: the *ignore_residues* keyword.

    def __init__(self, filename, **kwargs):
        """Set up parameters to run X3DNA_ on PDB *filename*.

        :Arguments:

          *filename*

               The *filename* is used as input for X3DNA in the "find_pair"    
               command.  It specifies the name of a PDB coordinate file
               to be used. This must be in Brookhaven protein databank format
               or something closely approximating this. Both ATOM and HETATM
               records are read. Note that if water molecules or ions are
               present in the channel these can be ignored on read by the use
               of the *ignore_residues* keyword.

               **Wildcard pattern**. A new feature (in release 2.1 of X3DNA) was
               the option to include a wild card (``*``) in the filename. e.g.,
               *filename* = `"ab*.pdb"` will apply x3dna to all files in the
               directory whose name starts with ``ab`` and ends with
               ``.pdb``. This is intended to aid the analysis of multiple
               copies of the same molecule - produced during molecular dynamics
               or other method. The x3dna procedure will be applied to each file
               in turn with the same setup conditions (initial point, sampling
               distance etc.). Graphics files will contain a combination of the
               individual runs, one after another. Note that the pdb files are
               read independently so that they need not have an identical
               number of atoms or atom order etc. (though they should be
               sufficiently similar for a X3DNA run from identical starting
               conditions to be useful).

               .. SeeAlso::

                  An alternative way to load in multiple files is a direct read
                  from a CHARMM binary dynamics DCD coordinate file - using the
                  *dcd* keyword or use :class:`X3DNAtraj`.


        :Keywords:

          *dcd*
               DCD trajectory (must be supplied together with a matching
               PDB file *filename*) and then X3DNA runs its analysis on each frame.

               It does multiple X3DNA runs on positions taken from a CHARMM binary
               dynamics format .DCD trajectory file. The *dcd* file must have
               exactly the same number of atoms in exactly the same order as
               the pdb file specified by *filename*. Note that if this option
               is used the pdb file is used as a template only - the
               coordinates are ignored. Note that structural parameters
               determined for each individual structure are written in a tagged
               format so that it is possible to extract the information from
               the text output file using a :program:`grep` command. The
               reading of the file can be controlled by the *step* keyword
               and/or setting :attr:`X3DNA.dcd_iniskip` to the number of frames
               to be skipped initially.

               .. Note::

                  X3DNA is very picky and does not read all DCD-like formats. If
                  in doubt, look into the *logfile* for error diagnostics.

                  At the moment, DCDs generated with MDAnalysis are not
                  accepted by X3DNA — use :class:`X3DNAtraj`, which works with
                  anything that MDAnalysis can read.

          *logfile*

               name of the file collecting X3DNA's output (which can be parsed
               using :meth:`X3DNA.collect` ["tmp*.out"]

          *step*

               step size for going through the trajectory (skips *step* - 1
               frames) [1]

          *executable*

               Path to the :program:`find_pair` executable directories
               (e.g. ``/opt/x3dna/2.1 and /opt/x3dna/2.1/bin``) must be set 
               and then added to export in bashrc file. See X3DNA 
               documentation for set-up instructions.
        """
        # list of temporary files, to be cleaned up on __del__
        self.tempfiles = ["auxiliary.par","bestpairs.pdb","bp_order.dat","bp_helical.par","cf_7methods.par","col_chains.scr","col_helices.scr","hel_regions.pdb","ref_frames.dat","hstacking.pdb","stacking.pdb"]
        self.tempdirs = []
        self.filename = filename
        #self.coordinates = self.check_and_fix_long_filename(self.filename)
        self.dcd = kwargs.pop('dcd', None)
        #if self.dcd:
        #    self.dcd = self.check_and_fix_long_filename(self.dcd)
        self.dcd_step = kwargs.pop("step", 1) - 1   # X3DNA docs description is confusing: step or skip??
        self.dcd_iniskip = 0
        self.shorto = int(kwargs.pop("shorto", 0))     # look at using SHORTO 2 for minimum output

        logger.info("Setting up X3DNA analysis for %(filename)r", vars(self))

        # guess executables
        self.exe = {}
        x3dna_exe_name = kwargs.pop('executable', 'find_pair')
        self.exe['find_pair'] = which(x3dna_exe_name)
        if self.exe['find_pair'] is None:
            errmsg = "X3DNA binary %(x3dna_exe_name)r not found." % vars()
            logger.fatal(errmsg)
            logger.fatal("%(x3dna_exe_name)r must be on the PATH or provided as keyword argument 'executable'.",
                         vars())
            raise OSError(errno.ENOENT, errmsg)
        x3dnapath = os.path.dirname(self.exe['find_pair'])
        self.logfile = kwargs.pop("logfile", "bp_step.par")
       
        self.template = textwrap.dedent("""find_pair %(filename)r stdout |analyze stdin """)

        if self.dcd:
            # CHARMD -- DCD (matches COORD)
            # CHARMS int int -- ignore_first_N_frames   skip_every_X_frames
            # http://d2o.bioch.ox.ac.uk:38080/doc/hole_d03.html#CHARMD
            self.template += "\nCHARMD %(dcd)s\nCHARMS %(dcd_iniskip)d %(dcd_step)d\n"

        # sanity checks
        if self.shorto > 2:
            logger.warn("SHORTO (%d) needs to be < 3 in order to extract a X3DNA profile!",
                        self.shorto)
        for program, path in self.exe.items():
            if path is None or which(path) is None:
                logger.error("Executable %(program)r not found, should have been %(path)r.",
                             vars())
        # results
        self.profiles = {}

    def check_and_fix_long_filename(self, filename, tmpdir=os.path.curdir):
        """Return *filename* suitable for X3DNA.

        X3DNA is limited to filenames <= :attr:`X3DNA.X3DNA_MAX_LENGTH`. This method

        1. returns *filename* if X3DNA can process it
        2. returns a relative path (see :func:`os.path.relpath`) if that shortens the
           path sufficiently
        3. creates a symlink to *filename* (:func:`os.symlink`) in a safe temporary
           directory and returns the path of the symlink. The temporary directory and
           the symlink are stored in :attr:`X3DNA.tempfiles` and :attr:`X3DNA.tempdirs`
           and deleted when the :class:`X3DNA` instance is deleted or garbage collected.

           By default the temporary directory is created inside the current
           directory in order to keep that path name short. This can be changed
           with the *tmpdir* keyword (e.g. one can use "/tmp").
        """
        #if len(filename) <= self.X3DNA_MAX_LENGTH:
        #    return filename

        logger.debug("path check: X3DNA will not read %r because it has more than %d characters.", filename)
        # try a relative path
        newname = filename #os.path.relpath(filename)

        # shorten path by creating a symlink inside a safe temp dir
        root,ext = filename #os.path.splitext(filename)
        dirname = tempfile.mkdtemp(dir=tmpdir)
        newname = os.path.join(dirname, os.path.basename(filename))
        self.tempfiles.append(newname)
        self.tempdirs.append(dirname)
        os.symlink(filename, newname)
        logger.debug("path check: Using symlink: %r --> %r", filename, newname)
        return newname

    def run(self, **kwargs):
        """Run X3DNA on the input file."""
        inpname = kwargs.pop("inpfile", None)
        outname = kwargs.pop("outfile", self.logfile)
        
        x3dnaargs = vars(self).copy()
        x3dnaargs.update(kwargs)

        inp = self.template % x3dnaargs
        if inpname:
            with open(inpname, "w") as f:
                f.write(inp)
            logger.debug("Wrote X3DNA input file %r for inspection", inpname)

        logger.info("Starting X3DNA on %(filename)r (trajectory: %(dcd)r)", x3dnaargs)
        logger.debug("%s", self.exe['find_pair'])
        with open(outname, "w") as output:
            x3dna = subprocess.call([inp],shell=True)
        with open(outname, "r") as output:
            # X3DNA is not very good at setting returncodes so check ourselves
            for line in output:
                if line.strip().startswith(('*** ERROR ***', 'ERROR')):
                    x3dna.returncode = 255
                    break
        if x3dna.bit_length != 0:
            logger.fatal("X3DNA Failure (%d). Check output %r", x3dna.bit_length, outname)
            #if stderr is not None:
            #    logger.fatal(stderr)
            #raise ApplicationError(x3dna.returncode, "X3DNA %r failed. Check output %r." % (self.exe['find_pair'], outname))
        logger.info("X3DNA finished: output file %(outname)r", vars())

    def collect(self, **kwargs):
        """Parse the output from a X3DNA run into numpy recarrays.

        Can deal with outputs containing multiple frames. Output format::


        The method saves the result as :attr:`X3DNA.profiles`, a dictionary
        indexed by the frame number. Each entry is a
        :class:`numpy.recarray`.

        If the keyword *outdir* is supplied (e.g. ".") then each profile is
        saved to a gzipped data file.

        :Keywords:
           *run*
              identifier, free form [1]
           *outdir*
              save output data under *outdir*/*run* if set to any other
              value but ``None`` [``None``]

        """
        #        Shear    Stretch   Stagger   Buckle   Prop-Tw   Opening     Shift     Slide     Rise      Tilt      Roll      Twist
        #0123456789.0123456789.0123456789.0123456789.0123456789.0123456789.
        #            11          22          33          44
        #123456789.123456789.123456789.123456789.123456789.123456789.123456789.
        #1           13
        #T-A     -0.033    -0.176     0.158   -12.177    -8.979     1.440     0.000     0.000     0.000     0.000     0.000     0.000
        #C-G     -0.529     0.122    -0.002    -7.983   -10.083    -0.091    -0.911     1.375     3.213    -0.766    -4.065    41.492
        # only parse bp_step.par      
        x3dna_output = kwargs.pop("x3dnaout", self.logfile)
        run = kwargs.pop("run", 1)     # id number
        outdir = kwargs.pop("outdir", os.path.curdir)

        logger.info("Collecting X3DNA profiles for run with id %s", run)
        length = 1   # length of trajectory --- is this really needed?? No... just for info
        if '*' in self.filename:
            import glob
            filenames = glob.glob(self.filename)
            length = len(filenames)
            if length == 0:
                logger.error("Glob pattern %r did not find any files.", self.filename)
                raise ValueError("Glob pattern %r did not find any files." % (self.filename,))
            logger.info("Found %d input files based on glob pattern %s", length, self.filename)
        if self.dcd:
            from MDAnalysis import Universe
            u = Universe(self.filename, self.dcd)
            length = int((u.trajectory.numframes - self.dcd_iniskip)/(self.dcd_step+1))
            logger.info("Found %d input frames in DCD trajectory %r", length, self.dcd)

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
                #print line
                # not sure if correct startswith
                if line.startswith("#        Shear"):
                    read_data = True
                    logger.debug("Started reading data")
                    fields = line.split()
                    # not sure if correct fields
                    x3dna_profile_no = int(13) #fields[1:])
                    records = []
                    #print fields
                    continue
                if read_data:
                    if len(line.strip()) != 0:
                        try:
                            Sequence, Shear,Stretch,Stagger,Buckle,Propeller,Opening,Shift,Slide,Rise,Tilt,Roll,Twist = line.split() #x3dnaformat.read(line)
                        except:
                            logger.critical("Run %d: Problem parsing line %r", run, line.strip())
                            logger.exception("Check input file %r.", x3dna_output)
                            raise
                        records.append([float(Shear),float(Stretch),float(Stagger),float(Buckle),float(Propeller),float(Opening),float(Shift),float(Slide),float(Rise),float(Tilt),float(Roll),float(Twist)])
                        #print records
                        #print len(records)
                        continue
                    else:
                        # end of records (empty line)

                        read_data = False
            frame_x3dna_output = numpy.rec.fromrecords(records, formats="f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8", names="Shear,Stretch,Stagger,Buckle,Propeller,Opening,Shift,Slide,Rise,Tilt,Roll,Twist")
            # store the profile
            self.profiles[x3dna_profile_no] = frame_x3dna_output
            #print self.profiles
            logger.debug("Collected X3DNA profile for frame %d (%d datapoints)",
                x3dna_profile_no, len(frame_x3dna_output))
            # save a profile for each frame (for debugging and scripted processing)
            # a tmp folder for each trajectory
            if outdir is not None:
                rundir = os.path.join(outdir, "run_"+str(run))
                os.system("mv tmp*.out %s"%rundir)
                if not os.path.exists(rundir):
                    os.makedirs(rundir)
                frame_x3dna_txt = os.path.join(rundir, "bp_step_%s_%04d.dat.gz" % (run, x3dna_profile_no))
                numpy.savetxt(frame_x3dna_txt, frame_x3dna_output)
                logger.debug("Finished with frame %d, saved as %r", x3dna_profile_no, frame_x3dna_txt)
                #continue
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

          *orderparameters*
               Sequence or text file with list of numbers corresponding to the
               frames in the trajectory.

          *start*, *stop*, *step*
               frame indices to slice the trajectory as
               ``universe.trajectory[start, stop, step]``

          *kwargs*
               All other keywords are passed on to :class:`X3DNA` (see there for description).

        """
        self.universe = universe
        self.selection = kwargs.pop("selection", "nucleic")
        self.orderparametersfile = kwargs.pop("orderparameters", None)

        self.start = kwargs.pop('start', None)
        self.stop = kwargs.pop('stop', None)
        self.step = kwargs.pop('step', None)


        self.x3dna_kwargs = kwargs

        # processing
        self.orderparameters = self._process_orderparameters(self.orderparametersfile)

    def _process_orderparameters(self, data):
        """Read orderparameters from *data*

        * If *data* is a string: Read orderparameters from *filename*.
        * If data is a array/list: use as is
        * If ``None``: assign frame numbers from trajectory
        """
        if isinstance(data, basestring):
            q = numpy.loadtxt(data)
        elif data is None:
            # frame numbers
            q = numpy.arange(1, self.universe.trajectory.numframes+1)
        else:
            q = numpy.asarray(data)

        if len(q.shape) != 1:
            raise TypeError("Order parameter array must be 1D.")
        if len(q) != self.universe.trajectory.numframes:
            errmsg = "Not same number of orderparameters ({0}) as trajectory frames ({1})".format(
                len(q), self.universe.trajectory.numframes)
            logger.error(errmsg)
            raise ValueError(errmsg)
        return q
    def run(self, **kwargs):
        """Run X3DNA on the whole trajectory and collect profiles.

        Keyword arguments *start*, *stop*, and *step* can be used to only
        analyse part of the trajectory.
        """
        from itertools import izip

        start = kwargs.pop('start', self.start)
        stop = kwargs.pop('stop', self.stop)
        step = kwargs.pop('step', self.step)
        x3dna_kw = self.x3dna_kwargs.copy()
        x3dna_kw.update(kwargs)

        profiles = {}  # index by orderparameters: NOTE: can overwrite!

        nucleic = self.universe.selectAtoms(self.selection)
        for q,ts in izip(self.orderparameters[start:stop:step], self.universe.trajectory[start:stop:step]):
            logger.info("X3DNA analysis frame %4d (orderparameter %g)", ts.frame, q)
            fd, pdbfile = tempfile.mkstemp(suffix=".pdb")
            try:
                nucleic.write(pdbfile)
                x3dna_profiles = self.run_x3dna(pdbfile, **x3dna_kw)
                print "yes"
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
            profiles[q] = x3dna_profiles.values()[0]
        self.profiles = profiles

    def run_x3dna(self, pdbfile, **kwargs):
        """Run X3DNA on a single PDB file *pdbfile*."""
        H = X3DNA(pdbfile, **kwargs)
        H.run()
        H.collect()
        return H.profiles

