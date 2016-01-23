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
Generation and Analysis of HOLE pore profiles --- :mod:`MDAnalysis.analysis.hole`
=================================================================================

:Author: Lukas Stelzl, Oliver Beckstein
:Year: 2011-2012
:Copyright: GNU Public License v2

With the help of this module, HOLE_ can be run on frames in a trajectory. Data
can be combined and analyzed. HOLE_ [Smart1993]_ [Smart1996]_ must be installed
separately and can be obtained from http://www.smartsci.uk/hole/.


.. rubric:: References

.. [Smart1993] O.S. Smart, J.M. Goodfellow and B.A. Wallace.
               The Pore Dimensions of Gramicidin A. Biophysical Journal 65:2455-2460, 1993.
.. [Smart1996] O.S. Smart, J.G. Neduvelil, X. Wang, B.A. Wallace, and M.S.P. Sansom.
               HOLE: A program for the analysis of the pore dimensions of ion channel
               structural models. J.Mol.Graph., 14:354–360, 1996.
               URL http://www.smartsci.uk/hole/

.. _HOLE: http://www.smartsci.uk/hole/


Examples
--------

Single structure
~~~~~~~~~~~~~~~~

Gramicidin A (gA) channel::

   from MDAnalysis.analysis.hole import HOLE, HOLEtraj
   from MDAnalysis.tests.datafiles import PDB_HOLE

   H = HOLE(PDB_HOLE, executable="~/hole2/bin/hole")  # set path to your hole binary
   H.run()
   H.collect()
   H.plot(linewidth=3, color="black", label=False)


Trajectory
~~~~~~~~~~

Analyzing a trajectory::

  u = MDAnalysis.Universe(psf, trajectory)
  H = HOLEtraj(u, ...)
  H.run()
  H.plot3D()

The profiles are available as the attribute :attr:`HOLEtraj.profiles`
(``H.profiles`` in the example) and are indexed by frame number but
can also be indexed by an arbitrary order parameter as shown in the
next example.


Trajectory with RMSD as order parameter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to classify the HOLE profiles the RMSD to a reference
structure is calculated for each trajectory frame (e.g. using the
:class:`MDAnalysis.analysis.rms.RMSD` analysis class). Then the HOLE
profiles can be ordered by the RMSD, which acts as an order
parameter. ::

  import MDAnalysis.analysis.hole
  import MDAnalysis.analysis.rms
  MDAnalysis.start_logging()
  ref = MDAnalysis.Universe("ref.psf", "ref.pdb")         # reference structure
  u = MDAnalysis.Universe("ref.psf", "traj.xtc")          # trajectory

  # calculate RMSD
  R = MDAnalysis.analysis.rms.RMSD(u, reference=ref, select="protein", mass_weighted=True)
  R.run()

  # HOLE analysis with order parameters
  H = MDAnalysis.analysis.hole.HOLEtraj(u, cvect=[0,0,1], orderparameters=R.rmsd[:,2])
  H.run()

The :attr:`HOLEtraj.profiles` dictionary will have the order parameter as key
for each frame. The plot functions will automatically sort the profiles by
ascending order parameter. To access the individual profiles one can simply
iterate over the sorted profiles (see :meth:`HOLEtraj.sorted_profiles_iter`) ::

  for q, profile in H:
     print("orderparameter = %g" % q)
     print("min(R) = %g" % profile.radius.min())


Data structures
---------------

A profile is stored as a :class:`numpy.recarray`::

   frame   rxncoord   radius

frame
   integer frame number (only important when HOLE itself reads a
   trajectory)
rxncoord
   the distance along the pore axis, in Å
radius
   the pore radius, in Å

The :attr:`HOLE.profiles` or :attr:`HOLEtraj.profiles` dictionary
holds one profile for each key. By default the keys are the frame
numbers but :class:`HOLEtraj` can take the optional *orderparameters*
keyword argument and load an arbitrary order parameter for each
frame. In that case, the key becomes the orderparameter.

.. Note::
   The profiles dict is not ordered and hence one typically needs to manually
   order the keys first. Furthermore, duplicate keys are not possible:
   In the case of *duplicate orderparameters*, the last one read will
   be stored in the dict.



Analysis
--------

.. autoclass:: HOLE
   :members:
   :inherited-members:

   .. attribute:: profiles

      After running :meth:`HOLE.collect`, this dict contains all the
      HOLE profiles, indexed by the frame number. If only a single
      frame was analyzed then this will be ``HOLE.profiles[0]``. Note
      that the order is random; one needs to sort the keys first.

.. autoclass:: HOLEtraj
   :members:
   :inherited-members:

   .. attribute:: profiles

      After running :meth:`HOLE.collect`, this dict contains all the
      HOLE profiles, indexed by the frame number or the order
      parameter (if *orderparameters* was supplied). Note that the
      order is random; one needs to sort the keys first.

.. autodata:: SIMPLE2_RAD


Utilities
---------

.. autofunction:: write_simplerad2
.. autofunction:: seq2str
.. autoclass:: ApplicationError


"""

from six.moves import zip, cPickle
import six

import numpy as np
import glob
import os
import errno
import shutil
import warnings
import os.path
import subprocess
import tempfile
import textwrap
import logging
from itertools import cycle

from MDAnalysis import Universe
from MDAnalysis import ApplicationError
from MDAnalysis.lib.util import which, realpath, asiterable
from MDAnalysis.lib.util import FORTRANReader


logger = logging.getLogger("MDAnalysis.analysis.hole")

#: Built-in HOLE radii (based on ``simple.rad`` from the HOLE_ distribution):
#: van der Waals radii are AMBER united atom from Weiner et al. (1984), JACS, vol 106 pp765-768.
#: *Simple* - Only use one value for each element C O H etc.
#: Added radii for K+, NA+, CL- (Pauling hydration radius from Hille 2002).
#: The data file can be written with the convenience function :func:`write_simplerad2`.
SIMPLE2_RAD = """
remark: Time-stamp: <2005-11-21 13:57:55 oliver> [OB]
remark: van der Waals radii: AMBER united atom
remark: from Weiner et al. (1984), JACS, vol 106 pp765-768
remark: Simple - Only use one value for each element C O H etc.
remark: van der Waals radii
remark: general last
VDWR C??? ??? 1.85
VDWR O??? ??? 1.65
VDWR S??? ??? 2.00
VDWR N??? ??? 1.75
VDWR H??? ??? 1.00
VDWR H?   ??? 1.00
VDWR P??? ??? 2.10
remark: ASN, GLN polar H (odd names for these atoms in xplor)
VDWR E2?  GLN 1.00
VDWR D2?  ASN 1.00
remark: amber lone pairs on sulphurs
VDWR LP?? ??? 0.00
remark: for some funny reason it wants radius for K even though
remark: it is on the exclude list
remark: Use Pauling hydration radius (Hille 2001) [OB]
VDWR K?   ??? 1.33
VDWR NA?  ??? 0.95
VDWR CL?  ??? 1.81
remark: funny hydrogens in gA structure [OB]
VDWR 1H?  ??? 1.00
VDWR 2H?  ??? 1.00
VDWR 3H?  ??? 1.00
remark: need bond rad for molqpt option
BOND C??? 0.85
BOND N??? 0.75
BOND O??? 0.7
BOND S??? 1.1
BOND H??? 0.5
BOND P??? 1.0
BOND ???? 0.85
"""


def write_simplerad2(filename="simple2.rad"):
    """Write the built-in radii in :data:`SIMPLE2_RAD` to *filename*.

    Does nothing if *filename* already exists.
    """
    if not os.path.exists(filename):
        with open(filename, "w") as rad:
            rad.write(SIMPLE2_RAD + "\n")
        logger.debug("Created simple radii file %(filename)r", vars())
    return filename


def seq2str(v):
    """Return sequence *v* as a string of numbers with spaces as separators.

    In the special case of ``None``, the empty string "" is returned.
    """
    if v is None:
        return ""
    return " ".join([str(x) for x in v])


class BaseHOLE(object):
    """Baseclass for HOLE analysis, providing plotting and utility functions"""

    def save(self, filename="hole.pickle"):
        """Save :attr:`profiles` as a Python pickle file *filename*.

        Load profiles dictionary with ::

           import cPickle
           profiles = cPickle.load(open(filename))

        """
        cPickle.dump(self.profiles, open(filename, "wb"), cPickle.HIGHEST_PROTOCOL)

    def _process_plot_kwargs(self, kwargs):
        import matplotlib.colors

        kw = {}
        frames = kwargs.pop('frames', None)
        if frames is None:
            frames = np.sort(self.profiles.keys()[::kwargs.pop('step', 1)])
        else:
            frames = asiterable(frames)
        kw['frames'] = frames
        kw['yshift'] = kwargs.pop('yshift', 0.0)
        kw['rmax'] = kwargs.pop('rmax', None)
        kw['label'] = kwargs.pop('label', True)
        if kw['label'] == "_nolegend_":
            kw['label'] = False
        elif kw['label'] is None:
            kw['label'] = True
        color = kwargs.pop('color', None)
        if color is None:
            cmap = kwargs.pop('cmap', matplotlib.cm.jet)
            normalize = matplotlib.colors.normalize(vmin=np.min(frames), vmax=np.max(frames))
            colors = cmap(normalize(frames))
        else:
            colors = cycle(asiterable(color))
        kw['colors'] = colors
        kw['linestyles'] = cycle(asiterable(kwargs.pop('linestyle', '-')))
        return kw, kwargs

    def plot(self, **kwargs):
        """Plot profiles in a 1D graph.

        Lines are coloured according to the colour map *cmap*.

        :Keywords:
           *step*
              only plot every *step* profile [1]
           *cmap*
              matplotlib color map to continuously color graphs [jet]
           *color*
              color or iterable of colors to specify graph colors;
              ``None`` will use *cmap* instead [``None``]
           *linestyle*
              linestyle supported by :func:`matplotlib.pyplot.plot`;
              an iterable will be applied in turn ['-']
           *yshift*
              displace each R(zeta) profile by *yshift* in the
              y-direction [0]
           *frames*
              integer or list: only plot these specific frame(s);
              the default is to plot everything (see also *step*).
              [``None``]
           *label*
              ignored and replaced with the frame number; if it
              is set to "_nolegend_" or ``False`` then no
              legend is displayed
           *ax*
              :class:`matplotlib.axes.Axes` instance to plot into
              [``None``]

        All other *kwargs* are passed to :func:`matplotlib.pyplot.plot`.
        """
        import matplotlib.pyplot as plt

        kw, kwargs = self._process_plot_kwargs(kwargs)

        ax = kwargs.pop('ax', plt.subplot(111))
        for iplot, (frame, color, linestyle) in enumerate(zip(kw['frames'], kw['colors'], kw['linestyles'])):
            kwargs['color'] = color
            kwargs['linestyle'] = linestyle
            kwargs['zorder'] = -frame
            kwargs['label'] = str(frame)
            dy = iplot * kw['yshift']
            profile = self.profiles[frame]
            ax.plot(profile.rxncoord, profile.radius + dy, **kwargs)
            #logger.debug("Plotted frame %d with color %r and linestyle %r", frame, color, linestyle)
        ax.set_xlabel(r"pore coordinate $\zeta$ ($\AA$)")
        ax.set_ylabel(r"HOLE radius $r$ ($\AA$)")
        if kw['label']:
            ax.legend(loc="best")

    def plot3D(self, **kwargs):
        """Stacked graph of profiles.

        Lines are coloured according to the colour map *cmap*.

        :Keywords:
           *step*
              only plot every *step* profile [1]
           *cmap*
              matplotlib color map [jet]
           *rmax*
              only display radii up to *rmax* [``None``]
           *ylabel*
              label of the reaction coordinate axis ["frames"]

        Based on Stack Overflow post `3d plots using matplotlib`_.

        .. _`3d plots using matplotlib`:
           http://stackoverflow.com/questions/9053255/3d-plots-using-matplotlib

        """
        import matplotlib.pyplot as plt
        import mpl_toolkits.mplot3d.axes3d as axes3d

        kw, kwargs = self._process_plot_kwargs(kwargs)
        rmax = kw.pop('rmax', None)
        ylabel = kwargs.pop('ylabel', "frames")

        fig = plt.figure(figsize=kwargs.pop('figsize', (6, 6)))
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        for frame, color, linestyle in zip(kw['frames'], kw['colors'], kw['linestyles']):
            kwargs['color'] = color
            kwargs['linestyle'] = linestyle
            kwargs['zorder'] = -frame
            profile = self.profiles[frame]
            if rmax is None:
                radius = profile.radius
                rxncoord = profile.rxncoord
            else:
                # does not seem to work with masked arrays but with nan hack!
                # http://stackoverflow.com/questions/4913306/python-matplotlib-mplot3d-how-do-i-set-a-maximum-value
                # -for-the-z-axis
                #radius = np.ma.masked_greater(profile.radius, rmax)
                #rxncoord = np.ma.array(profile.rxncoord, mask=radius.mask)
                rxncoord = profile.rxncoord
                radius = profile.radius.copy()
                radius[radius > rmax] = np.nan
            ax.plot(rxncoord, frame * np.ones_like(rxncoord), radius, **kwargs)

        ax.set_xlabel(r"pore coordinate $z$")
        ax.set_ylabel(ylabel)
        ax.set_zlabel(r"HOLE radius $r$")
        plt.draw()

    def min_radius(self):
        """Return the minimum radius over all profiles as a function of q"""
        f = []
        for q, profile in self:
            f.append([q, profile.radius.min()])
        return np.array(f)

    def sorted_profiles_iter(self):
        """Return an iterator over profiles sorted by frame/order parameter *q*.

        The iterator produces tuples ``(q, profile)``.
        """
        if self.profiles is None:
            raise StopIteration
        for q in sorted(self.profiles):
            yield (q, self.profiles[q])

    __iter__ = sorted_profiles_iter


class HOLE(BaseHOLE):
    """Run HOLE on a single frame or a DCD trajectory.

    Only a subset of all `HOLE control parameters`_ is supported and can be set
    with keyword arguments.

    HOLE (as a FORTRAN77 program) has a number of limitations when it comes to
    filename lengths (must be shorter than the empirically found
    :attr:`HOLE.HOLE_MAX_LENGTH`). This class tries to work around them by
    creating temporary symlinks to files when needed but this can still fail
    when permissions are not correctly set on the current directory.

    Running HOLE with the :class:`HOLE` class is a 3-step process:

     1. set up the class with all desired parameters
     2. run HOLE with :meth:`HOLE.run`
     3. collect the data from the output file with :meth:`HOLE.collect`

    The class also provides some simple plotting functions of the collected
    data such as :meth:`HOLE.plot` or :meth:`HOLE.plot3D`.

    .. versionadded:: 0.7.7

    .. _`HOLE control parameters`:
       http://s3.smartsci.uk/hole2/doc/old/hole_d03.html
    """
    #: Maximum number of characters in a filename (limitation of HOLE)
    HOLE_MAX_LENGTH = 70

    #: List of residues that are ignore by default. Can be changed with
    #: the *ignore_residues* keyword.
    default_ignore_residues = ["SOL", "WAT", "TIP", "HOH", "K  ", "NA ", "CL "]

    def __init__(self, filename, **kwargs):
        """Set up parameters to run HOLE_ on PDB *filename*.

        :Arguments:

          *filename*

               The *filename* is used as input for HOLE in the "COORD" card of
               the input file.  It specifies the name of a PDB co-ordinate file
               to be used. This must be in Brookhaven protein databank format
               or something closely approximating this. Both ATOM and HETATM
               records are read. Note that if water molecules or ions are
               present in the channel these can be ignored on read by the use
               of the *ignore_residues* keyword.

               **Wildcard pattern**. A new feature (in release 2.1 of HOLE) was
               the option to include a wild card (``*``) in the filename. e.g.,
               *filename* = `"ab*.pdb"` will apply hole to all files in the
               directory whose name starts with ``ab`` and ends with
               ``.pdb``. This is intended to aid the analysis of multiple
               copies of the same molecule - produced during molecular dynamics
               or other method. The hole procedure will be applied to each file
               in turn with the same setup conditions (initial point, sampling
               distance etc.). Graphics files will contain a combination of the
               individual runs, one after another. Note that the pdb files are
               read independently so that they need not have an identical
               number of atoms or atom order etc. (though they should be
               sufficiently similar for a HOLE run from identical starting
               conditions to be useful).

               .. SeeAlso::

                  An alternative way to load in multiple files is a direct read
                  from a CHARMM binary dynamics DCD coordinate file - using the
                  *dcd* keyword or use :class:`HOLEtraj`.


        :Keywords:

          *dcd*
               DCD trajectory (must be supplied together with a matching
               PDB file *filename*) and then HOLE runs its analysis on each frame.

               It does multiple HOLE runs on positions taken from a CHARMM binary
               dynamics format .DCD trajectory file. The *dcd* file must have
               exactly the same number of atoms in exactly the same order as
               the pdb file specified by *filename*. Note that if this option
               is used the pdb file is used as a template only - the
               coordinates are ignored. Note that structural parameters
               determined for each individual structure are written in a tagged
               format so that it is possible to extract the information from
               the text output file using a :program:`grep` command. The
               reading of the file can be controlled by the *step* keyword
               and/or setting :attr:`HOLE.dcd_iniskip` to the number of frames
               to be skipped initially.

               .. Note::

                  HOLE is very picky and does not read all DCD-like formats. If
                  in doubt, look into the *logfile* for error diagnostics.

                  At the moment, DCDs generated with MDAnalysis are not
                  accepted by HOLE — use :class:`HOLEtraj`, which works with
                  anything that MDAnalysis can read.

          *logfile*

               name of the file collecting HOLE's output (which can be parsed
               using :meth:`HOLE.collect` ["hole.out"]

          *sphpdb*

               name of the HOLE sph file, a PDB-like file containig the coordinates
               of the pore centers ["hole.sph"]

          *step*

               step size for going through the trajectory (skips *step* - 1
               frames) [1]

          *cpoint*

               coordinates of a point inside the pore, e.g. ``[12.3, 0.7,
               18.55]``. If ``None`` then HOLE's own simple search algorithm is
               used.

               This specifies a point which lies within the channel, for simple
               channels such as gramicidin results do not show great
               sensitivity to the exact point taken. An easy way to produce an
               initial point is to use molecular graphics to find two atoms
               which lie either side of the pore and to average their
               co-ordinates. Or if the channel structure contains water
               molecules or counter ions then take the coordinates of one of
               these (and use the *ignore_residues* keyword to ignore them in
               the pore radius calculation).

               If this card is not specified then HOLE now (from version 2.2)
               attempts to make a guess where the channel will be. The
               procedure assumes the channel is reasonably symmetric. The
               initial guess on cpoint will be the centroid of all alpha carbon
               atoms (name 'CA' in pdb file). This is then refined by a crude
               grid search up to 5 Å from the original position. This procedure
               works most of the time but is clearly far from infallible —
               results should be careful checked (with molecular graphics) if
               it is used. [``None``]

          *cvect*

               Search direction, should be parallel to the pore axis,
               e.g. ``[0,0,1]`` for the z-axis.

               If this keyword is ``None`` then HOLE now attempts to make a
               guess where the channel will be. The procedure assumes the
               channel is reasonably symmetric. The guess will be either along
               the X axis (1,0,0), Y axis (0,1,0) or Z axis (0,0,1). If the
               structure is not aligned on one of these axis the results will
               clearly be approximate. If a guess is used then results should
               be carefully checked. [``None``]

          *sample*
               distance of sample points in Å

               Specifies the distance between the planes used in the HOLE
               procedure. The default value is 0.2 Å, this should be reasonable
               for most purposes. However, if you wish to visualize a very
               tight constriction then specify a smaller value.  [0.2]

          *dotden*
               density of facettes for generating a 3D pore representation

               This number controls the density of dots which will be used by
               the program. A sphere of dots is placed on each centre
               determined in the Monte Carlo procedure. Only dots which do not
               lie within any other sphere are considered. The actual number of
               dots written is therefore controlled by *dotden* and
               *sample*. *dotden* should be set to between 5 (few dots per
               sphere) and 35 (large number of dots per sphere). [15]

          *endrad*
               Radius which is considered to be the end of the pore. This
               keyword can be used to specify the radius above which the
               program regards a result as indicating that the end of the pore
               has been reached. The default value is 22.0 Å. This may need to
               be increased for large channels or reduced for small. [22.0]

          *shorto*

               Determines the output of output in the *logfile*; for automated processing
               this must be < 3.

               - 0: Full text output
               - 1: All text output given except "run in progress" (i.e.,
                 detailed contemporary description of what HOLE is doing).
               - 2: Ditto plus no graph type output - only leaving minimum
                 radius and conductance calculations.
               - 3: All text output other than input card mirroring and error messages
                 turned off.

          *ignore_residues*

               sequence of three-letter residues that are not taken into account during the
               calculation; wildcards are *not* supported
               [ ``["SOL","WAT", "TIP", "HOH", "K  ", "NA ", "CL "]`` ]

          *radius*
               Path to the radii; if set to None then a set of default radii, :data:`SIMPLE2_RAD`,
               is used (an extension of ``simple.rad`` from the HOLE distribution)

               This specifies the name for the file specifying van der Waals
               radii for each atom. A number of files with different values are
               supplied with HOLE.

          *sphpdb*

               This keyword specifies the filename for output of the sphere
               centre information in pdb form. Its typical suffix is
               ".sph". The co-ordinates are set to the sphere centres and the
               occupancies are the sphere radii. All centres are assigned the
               atom name QSS and residue name SPH and the residue number is set
               to the storage number of the centre. The file can be imported
               into molecular graphics programs but are likely to be bonded
               together in a awful manner - as the points are very close to one
               another. In VMD sph objects are best displayed as
               "Points". Displaying .sph objects rather than rendered or dot
               surfaces can be useful to analyze the distance of particular
               atoms from the sphere-centre line.

               Most usefully .sph files can be used to produce molecular
               graphical output from a hole run. This is achieved by using the
               :program:`sph_process` program to read the .sph
               file. ["hole.sph"]

          *executable*

               Path to the :program:`hole` executable
               (e.g. ``/opt/hole/exe/hole``); the other programs
               :program:`sph_process` and :program:`sos_triangle` are assumed
               to live in the same directory as :program:`hole`. If
               :program:`hole` is found on the :envvar:`PATH` then the bare
               executable name is sufficient.  ["hole"]
        """
        # list of temporary files, to be cleaned up on __del__
        self.tempfiles = []
        self.tempdirs = []

        self.filename = filename
        self.coordinates = self.check_and_fix_long_filename(self.filename)
        self.dcd = kwargs.pop('dcd', None)
        if self.dcd:
            self.dcd = self.check_and_fix_long_filename(self.dcd)
        self.dcd_step = kwargs.pop("step", 1) - 1  # HOLE docs description is confusing: step or skip??
        self.dcd_iniskip = 0
        self.cpoint = kwargs.pop("cpoint", None)
        self.cvect = kwargs.pop("cvect", None)
        self.sample = float(kwargs.pop("sample", 0.20))
        self.dotden = int(kwargs.pop("dotden", 15))
        self.endrad = float(kwargs.pop("endrad", 22.))
        self.shorto = int(kwargs.pop("shorto", 0))  # look at using SHORTO 2 for minimum output
        self.ignore_residues = kwargs.pop("ignore_residues", self.default_ignore_residues)
        self.radius = self.check_and_fix_long_filename(
            realpath(kwargs.pop('radius', None) or write_simplerad2()))

        logger.info("Setting up HOLE analysis for %(filename)r", vars(self))
        logger.info("Using radius file %(radius)r", vars(self))

        # guess executables
        self.exe = {}
        hole_exe_name = kwargs.pop('executable', 'hole')
        self.exe['hole'] = which(hole_exe_name)
        if self.exe['hole'] is None:
            errmsg = "HOLE binary {hole_exe_name!r} not found.".format(**vars())
            logger.fatal(errmsg)
            logger.fatal("%(hole_exe_name)r must be on the PATH or provided as keyword argument 'executable'.",
                         vars())
            raise OSError(errno.ENOENT, errmsg)
        holepath = os.path.dirname(self.exe['hole'])
        self.exe['sos_triangle'] = os.path.join(holepath, "sos_triangle")
        self.exe['sph_process'] = os.path.join(holepath, "sph_process")

        self.sphpdb = kwargs.pop("sphpdb", "hole.sph")
        self.logfile = kwargs.pop("logfile", "hole.out")

        self.template = textwrap.dedent("""
            ! Input file for Oliver Smart's HOLE program
            ! written by MDAnalysis.analysis.hole.HOLE
            ! filename = %(filename)s
            COORD  %(coordinates)s
            RADIUS %(radius)s
            SPHPDB %(sphpdb)s
            SAMPLE %(sample)f
            ENDRAD %(endrad)f
            IGNORE %(ignore)s
            SHORTO %(shorto)d
            """)
        if self.cpoint is not None:
            # note: if it is None then we can't change this with a kw for run() !!
            self.template += "CPOINT %(cpoint_xyz)s\n"
        else:
            logger.info("HOLE will guess CPOINT")
        if self.cvect is not None:
            # note: if it is None then we can't change this with a kw for run() !!
            self.template += "CVECT  %(cvect_xyz)s\n"
        else:
            logger.info("HOLE will guess CVECT")

        if self.dcd:
            # CHARMD -- DCD (matches COORD)
            # CHARMS int int -- ignore_first_N_frames   skip_every_X_frames
            #        http://s3.smartsci.uk/hole2/doc/old/hole_d03.html#CHARMD
            self.template += "\nCHARMD %(dcd)s\nCHARMS %(dcd_iniskip)d %(dcd_step)d\n"

        # sanity checks
        if self.shorto > 2:
            logger.warn("SHORTO (%d) needs to be < 3 in order to extract a HOLE profile!",
                        self.shorto)
        for program, path in self.exe.items():
            if path is None or which(path) is None:
                logger.error("Executable %(program)r not found, should have been %(path)r.",
                             vars())
        # results
        self.profiles = {}

    def check_and_fix_long_filename(self, filename, tmpdir=os.path.curdir):
        """Return *filename* suitable for HOLE.

        HOLE is limited to filenames <= :attr:`HOLE.HOLE_MAX_LENGTH`. This method

        1. returns *filename* if HOLE can process it
        2. returns a relative path (see :func:`os.path.relpath`) if that shortens the
           path sufficiently
        3. creates a symlink to *filename* (:func:`os.symlink`) in a safe temporary
           directory and returns the path of the symlink. The temporary directory and
           the symlink are stored in :attr:`HOLE.tempfiles` and :attr:`HOLE.tempdirs`
           and deleted when the :class:`HOLE` instance is deleted or garbage collected.

           By default the temporary directory is created inside the current
           directory in order to keep that path name short. This can be changed
           with the *tmpdir* keyword (e.g. one can use "/tmp").
        """
        if len(filename) <= self.HOLE_MAX_LENGTH:
            return filename

        logger.debug("path check: HOLE will not read %r because it has more than %d characters.", filename,
                     self.HOLE_MAX_LENGTH)
        # try a relative path
        newname = os.path.relpath(filename)
        if len(newname) <= self.HOLE_MAX_LENGTH:
            logger.debug("path check: Using relative path: %r --> %r", filename, newname)
            return newname

        # shorten path by creating a symlink inside a safe temp dir
        root, ext = os.path.splitext(filename)
        dirname = tempfile.mkdtemp(dir=tmpdir)
        newname = os.path.join(dirname, os.path.basename(filename))
        self.tempfiles.append(newname)
        self.tempdirs.append(dirname)
        if len(newname) > self.HOLE_MAX_LENGTH:
            logger.fatal("path check: Failed to shorten filename %r --> %r", filename, newname)
            raise RuntimeError("Failed to shorten filename %r --> %r", filename, newname)
        os.symlink(filename, newname)
        logger.debug("path check: Using symlink: %r --> %r", filename, newname)
        return newname

    def run(self, **kwargs):
        """Run HOLE on the input file."""
        inpname = kwargs.pop("inpfile", None)
        outname = kwargs.pop("outfile", self.logfile)
        # NOTE: If cvect and cpoint had been None in the constructor then they are
        #       ignored here. Arguably a bug... but then again, the keywords for run() are
        #       not even officially documented :-).
        kwargs.setdefault("cvect_xyz", seq2str(kwargs.pop('cvect', self.cvect)))
        kwargs.setdefault("cpoint_xyz", seq2str(kwargs.pop('cpoint', self.cpoint)))
        kwargs.setdefault("ignore", seq2str(kwargs.pop('ignore_residues', self.ignore_residues)))
        holeargs = vars(self).copy()
        holeargs.update(kwargs)

        inp = self.template % holeargs
        if inpname:
            with open(inpname, "w") as f:
                f.write(inp)
            logger.debug("Wrote HOLE input file %r for inspection", inpname)

        logger.info("Starting HOLE on %(filename)r (trajectory: %(dcd)r)", holeargs)
        logger.debug("%s <%s >%s", self.exe['hole'], (inpname if inpname else "(input)"), outname)
        with open(outname, "w") as output:
            hole = subprocess.Popen([self.exe['hole']], stdin=subprocess.PIPE, stdout=output)
            stdout, stderr = hole.communicate(inp)
        with open(outname, "r") as output:
            # HOLE is not very good at setting returncodes so check ourselves
            for line in output:
                if line.strip().startswith(('*** ERROR ***', 'ERROR')):
                    hole.returncode = 255
                    break
        if hole.returncode != 0:
            logger.fatal("HOLE Failure (%d). Check output %r", hole.returncode, outname)
            if stderr is not None:
                logger.fatal(stderr)
            raise ApplicationError(hole.returncode, "HOLE {0!r} failed. Check output {1!r}.".format(self.exe['hole'], outname))
        logger.info("HOLE finished: output file %(outname)r", vars())

    def create_vmd_surface(self, filename="hole.vmd", **kwargs):
        """Process HOLE output to create a smooth pore surface suitable for VMD.

        Takes the :attr:`sphpdb` file and feeds it to :program:`sph_process`
        and :program:`sos_triangle` as described under `Visualization of HOLE
        results`_.

        Load the output file *filename* into VMD by issuing in the tcl console ::

           source hole.vmd

        The level of detail is determined by :attr:`HOLE.dotden` (which can be
        overriden by keyword *dotden*).

        The surface will be colored so that parts that are inaccessible to
        water (pore radius < 1.15 Å) are red, water accessible parts (1.15 Å >
        pore radius < 2.30 Å) are green and wide areas (pore radius > 2.30 Å
        are blue).

        .. _`Visualization of HOLE results`: http://s3.smartsci.uk/hole2/doc/index.html#_producing_a_triangulated_surface_and_visualizing_in_vmd
        .. _`sph_process`: http://s3.smartsci.uk/hole2/doc/old/hole_d04.html#sph_process
        """
        # not sure how this works when run on multiple frames...
        # see http://s3.smartsci.uk/hole2/doc/old/hole_d04.html#sph_process
        kwargs.setdefault("dotden", self.dotden)

        fd, tmp_sos = tempfile.mkstemp(suffix=".sos", text=True)
        os.close(fd)
        try:
            rc = subprocess.call([
                self.exe["sph_process"], "-sos", "-dotden", str(kwargs['dotden']),
                "-color", self.sphpdb, tmp_sos])
            if rc != 0:
                logger.fatal("sph_process failed (%d)", rc)
                raise OSError(rc, "sph_process failed")
            with open(tmp_sos) as sos:
                with open(filename, "w") as triangles:
                    rc = subprocess.call([self.exe["sos_triangle"], "-s"], stdin=sos, stdout=triangles)
            if rc != 0:
                logger.fatal("sos_triangle failed (%d)", rc)
                raise OSError(rc, "sos_triangle failed")
        except:
            os.unlink(tmp_sos)
            raise
        return filename

    def collect(self, **kwargs):
        """Parse the output from a HOLE run into numpy recarrays.

        Can deal with outputs containing multiple frames. Output format::

            frame  rxncoord   radius

        The method saves the result as :attr:`HOLE.profiles`, a dictionary
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
        # cenxyz.cvec      radius  cen_line_D sum{s/(area point sourc
        #0123456789.0123456789.0123456789.0123456789.0123456789.0123456789.
        #            11          22          33          44
        #123456789.123456789.123456789.123456789.123456789.123456789.123456789.
        #1           13
        # 3F12
        #   -27.17082    15.77469   -73.19195     0.00013   (sampled)
        #   -27.07082    12.91103   -69.39840     0.00032 (mid-point)
        #
        # only parse cenxyz.cvec      radius  cen_line_D
        # because the area can be ********
        holeformat = FORTRANReader('3F12')

        hole_output = kwargs.pop("holeout", self.logfile)
        run = kwargs.pop("run", 1)  # id number
        outdir = kwargs.pop("outdir", os.path.curdir)

        logger.info("Collecting HOLE profiles for run with id %s", run)
        length = 1  # length of trajectory --- is this really needed?? No... just for info
        if '*' in self.filename:
            filenames = glob.glob(self.filename)
            length = len(filenames)
            if length == 0:
                logger.error("Glob pattern %r did not find any files.", self.filename)
                raise ValueError("Glob pattern {0!r} did not find any files.".format(self.filename))
            logger.info("Found %d input files based on glob pattern %s", length, self.filename)
        if self.dcd:
            u = Universe(self.filename, self.dcd)
            length = int((u.trajectory.n_frames - self.dcd_iniskip) / (self.dcd_step + 1))
            logger.info("Found %d input frames in DCD trajectory %r", length, self.dcd)

        # one recarray for each frame, indexed by frame number
        # TODO: use ordered dict
        self.profiles = {}

        logger.info("Run %s: Reading %d HOLE profiles from %r", run, length, hole_output)
        hole_profile_no = 0
        records = []
        with open(hole_output, "r") as hole:
            read_data = False
            for line in hole:
                line = line.rstrip()  # preserve columns (FORTRAN output...)
                if line.startswith(" Starting calculation for position number"):
                    fields = line.split()
                    hole_profile_no = int(fields[5])
                    records = []
                    logger.debug("Started reading profile %d", hole_profile_no)
                    continue
                if line.startswith(" cenxyz.cvec"):
                    read_data = True
                    logger.debug("Started reading data")
                    continue
                if read_data:
                    if len(line.strip()) != 0:
                        try:
                            rxncoord, radius, cenlineD = holeformat.read(line)
                        except:
                            logger.critical("Run %d: Problem parsing line %r", run, line.strip())
                            logger.exception("Check input file %r.", hole_output)
                            raise
                        records.append((hole_profile_no, rxncoord, radius))
                        continue
                    else:
                        # end of records (empty line)
                        read_data = False
                        frame_hole_output = np.rec.fromrecords(records, formats="i4,f8,f8",
                                                                  names="frame,rxncoord,radius")
                        # store the profile
                        self.profiles[hole_profile_no] = frame_hole_output
                        logger.debug("Collected HOLE profile for frame %d (%d datapoints)",
                                     hole_profile_no, len(frame_hole_output))
                        # save a profile for each frame (for debugging and scripted processing)
                        # a tmp folder for each trajectory
                        if outdir is not None:
                            rundir = os.path.join(outdir, "run_" + str(run))
                            if not os.path.exists(rundir):
                                os.makedirs(rundir)
                            frame_hole_txt = os.path.join(rundir, "radii_{0!s}_{1:04d}.dat.gz".format(run, hole_profile_no))
                            np.savetxt(frame_hole_txt, frame_hole_output)
                            logger.debug("Finished with frame %d, saved as %r", hole_profile_no, frame_hole_txt)
                        continue
                        # if we get here then we haven't found anything interesting
        if len(self.profiles) == length:
            logger.info("Collected HOLE radius profiles for %d frames", len(self.profiles))
        else:
            logger.warn("Missing data: Found %d HOLE profiles from %d frames.", len(self.profiles), length)

    def __del__(self):
        for f in self.tempfiles:
            try:
                os.unlink(f)
            except OSError:
                pass
        for d in self.tempdirs:
            shutil.rmtree(d, ignore_errors=True)


class HOLEtraj(BaseHOLE):
    """Analyze all frames in a trajectory.

    The :class:`HOLE` class provides a direct interface to HOLE. HOLE itself
    has limited support for analysing trajectories but cannot deal with all the
    trajectory formats understood by MDAnalysis. This class can take any
    universe and feed it to HOLE. By default it sequentially creates a PDB for
    each frame and runs HOLE on the frame.
    """

    def __init__(self, universe, **kwargs):
        """Set up the class.

        :Arguments:

          *universe*
               The input trajectory as part of a
               :class:`~MDAnalysis.core.AtomGroup.Universe`. trajectory is
               converted to a sequence of PDB files and HOLE is run on each
               individual file. (Use the *start*, *stop*, and *step* keywords
               to slice the trajectory.)

          *orderparameters*
               Sequence or text file with list of numbers corresponding to the
               frames in the trajectory.

          *start*, *stop*, *step*
               frame indices to slice the trajectory as
               ``universe.trajectory[start, stop, step]``

          *selection*
               selection string for :meth:`~MDAnalysis.core.AtomGroup.Universe.select_atoms` to select
               the group of atoms that is to be analysed by HOLE ["protein"]

          *cpoint*
               Point inside the pore.

               If set to ``True`` then *cpoint* is guessed as the
               :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.center_of_geometry` of
               the *selection* from the first frame of the trajectory.

               If *cpoint* is not set or set to ``None`` then HOLE guesses it
               with its own algorithm (for each individual frame). [``None``]

          *kwargs*
               All other keywords are passed on to :class:`HOLE` (see there for description).

        """
        self.universe = universe
        self.selection = kwargs.pop("selection", "protein")
        self.orderparametersfile = kwargs.pop("orderparameters", None)

        self.start = kwargs.pop('start', None)
        self.stop = kwargs.pop('stop', None)
        self.step = kwargs.pop('step', None)

        self.cpoint = kwargs.pop('cpoint', None)
        if self.cpoint is True:
            self.cpoint = self.guess_cpoint(selection=self.selection)
            logger.info("Guessed CPOINT = %r from selection %r", self.cpoint, self.selection)
        kwargs['cpoint'] = self.cpoint

        self.hole_kwargs = kwargs

        # processing
        self.orderparameters = self._process_orderparameters(self.orderparametersfile)

    def guess_cpoint(self, **kwargs):
        """Guess a point inside the pore.

        This method simply uses the center of geometry of the protein selection
        as a guess. *selection* is "protein" by default.
        """
        return self.universe.select_atoms(kwargs.get("selection", "protein")).center_of_geometry()

    def _process_orderparameters(self, data):
        """Read orderparameters from *data*

        * If *data* is a string: Read orderparameters from *filename*.
        * If data is a array/list: use as is
        * If ``None``: assign frame numbers from trajectory
        """
        if isinstance(data, six.string_types):
            q = np.loadtxt(data)
        elif data is None:
            # frame numbers
            q = np.arange(1, self.universe.trajectory.n_frames + 1)
        else:
            q = np.asarray(data)

        if len(q.shape) != 1:
            raise TypeError("Order parameter array must be 1D.")
        if len(q) != self.universe.trajectory.n_frames:
            errmsg = "Not same number of orderparameters ({0}) as trajectory frames ({1})".format(
                len(q), self.universe.trajectory.n_frames)
            logger.error(errmsg)
            raise ValueError(errmsg)
        return q

    def run(self, **kwargs):
        """Run HOLE on the whole trajectory and collect profiles.

        Keyword arguments *start*, *stop*, and *step* can be used to only
        analyse part of the trajectory.
        """
        start = kwargs.pop('start', self.start)
        stop = kwargs.pop('stop', self.stop)
        step = kwargs.pop('step', self.step)
        hole_kw = self.hole_kwargs.copy()
        hole_kw.update(kwargs)

        profiles = {}  # index by orderparameters: NOTE: can overwrite!

        # better: if not a DCD, convert to DCD because HOLE can read it (does it need a PSF?)
        # ... maybe faster than splittinginto PDBs...
        # but currently HOLE chokes on MDAnalysis-generated DCDs (header comment wrong/too long?)

        # probably slow (especially with additional linking for filename length) --
        # can't we do this in memory?? Or at least check that we can get a local tmp dir

        # TODO: alternatively, dump all frames with leading framenumber and use a wildcard
        #       (although the file renaming might create problems...)
        protein = self.universe.select_atoms(self.selection)
        for q, ts in zip(self.orderparameters[start:stop:step], self.universe.trajectory[start:stop:step]):
            logger.info("HOLE analysis frame %4d (orderparameter %g)", ts.frame, q)
            fd, pdbfile = tempfile.mkstemp(suffix=".pdb")
            os.close(fd)  # only need an empty file that can be overwritten, close right away (Issue 129)
            try:
                protein.write(pdbfile)
                hole_profiles = self.run_hole(pdbfile, **hole_kw)
            finally:
                try:
                    os.unlink(pdbfile)
                except OSError:
                    pass
            if len(hole_profiles) != 1:
                err_msg = "Got {0} profiles ({1}) --- should be 1 (time step {2})".format(
                    len(hole_profiles), hole_profiles.keys(), ts)
                logger.error(err_msg)
                warnings.warn(err_msg)
            profiles[q] = hole_profiles.values()[0]
        self.profiles = profiles

    def run_hole(self, pdbfile, **kwargs):
        """Run HOLE on a single PDB file *pdbfile*."""
        H = HOLE(pdbfile, **kwargs)
        H.run()
        H.collect()
        return H.profiles
