# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
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
AMBER trajectories --- :mod:`MDAnalysis.coordinates.TRJ`
========================================================

AMBER_ can write :ref:`ASCII trajectories<ascii-trajectories>` ("traj") and
:ref:`binary trajectories<netcdf-trajectories>` ("netcdf"). MDAnalysis supports
reading of both formats and writing for the binary trajectories.

.. Note::

   Support for AMBER is *experimental* and feedback and contributions
   are highly appreciated. Use the `Issue Tracker`_ or get in touch on
   the `MDAnalysis mailinglist`_.

.. rubric:: Units

* lengths in Angstrom (Å)
* time in ps (but see below)

AMBER trajectory coordinate frames are based on a :class:`Timestep`
object.

.. autoclass:: Timestep
   :members:

   .. attribute:: _pos

      coordinates of the atoms as a :class:`numpy.ndarray` of shape `(numatoms, 3)`

   .. attribute:: _velocities

      velocities of the atoms as a :class:`numpy.ndarray` of shape `(numatoms, 3)`;
      only available if the trajectory contains velocities or if the
      *velocities* = ``True`` keyword has been supplied.


.. _ascii-trajectories:

ASCII TRAJ trajectories
-----------------------

ASCII AMBER_ TRJ coordinate files (as defined in `AMBER TRJ format`_)
are handled by the :class:`TRJReader`. It is also possible to directly
read *bzip2* or *gzip* compressed files.

AMBER ASCII trajectories are recognised by the suffix '.trj' or
'.mdcrd' (possibly with an additional '.gz' or '.bz2').

.. rubric:: Limitations

* Periodic boxes are only stored as box lengths A, B, C in an AMBER
  trajectory; the reader always assumes that these are orthorhombic
  boxes.

* The trajectory does not contain time information so we simply set
  the time step to 1 ps (or the user could provide it as kwarg *delta*)

* No direct access of frames is implemented, only iteration through
  the trajectory.

* Trajectories with fewer than 4 atoms probably fail to be read (BUG).

* If the trajectory contains exactly *one* atom then it is always
  assumed to be non-periodic (for technical reasons).


.. autoclass:: TRJReader
   :members:


.. _netcdf-trajectories:

Binary NetCDF trajectories
--------------------------

The `AMBER netcdf`_ format make use of NetCDF_ (Network Common Data
Form) format. Such binary trajectories are recognized in MDAnalysis by
the '.ncdf' suffix and read by the :class:`NCDFReader`.

Binary trajectories can also contain velocities and can record the exact time
step. In principle, the trajectories can be in different units than the AMBER
defaults of ångström and picoseconds but at the moment MDAnalysis only supports
those and will raise a :exc:`NotImplementedError` if anything else is detected.

.. autoclass:: NCDFReader
   :members:

   .. automethod:: __getitem__
   .. automethod:: __iter__

.. autoclass:: NCDFWriter
   :members:


.. Links

.. _AMBER: http://ambermd.org
.. _AMBER TRJ format: http://ambermd.org/formats.html#trajectory
.. _AMBER netcdf format: http://ambermd.org/netcdf/nctraj.html
.. _AMBER netcdf: http://ambermd.org/netcdf/nctraj.html
.. _NetCDF: http://www.unidata.ucar.edu/software/netcdf
.. _Issue Tracker: http://code.google.com/p/mdanalysis/issues/list
.. _MDAnalysis mailinglist: http://groups.google.com/group/mdnalysis-discussion



"""

import numpy
import warnings

import MDAnalysis
import base
import MDAnalysis.core.util as util

import errno
import logging
logger = logging.getLogger("MDAnalysis.coordinates.AMBER")

try:
        import netCDF4 as netcdf
except ImportError:
        # Just to notify the user; the module will still load. However, NCDFReader and NCDFWriter
        # will raise a proper ImportError if they are called without the netCDF4 library present.
        # See Issue 122 for a discussion.
        logger.debug("Failed to import netCDF4; AMBER NETCDFReader/Writer will not work. "
                     "Install netCDF4 from http://code.google.com/p/netcdf4-python/.")
        logger.debug("See also https://code.google.com/p/mdanalysis/wiki/netcdf")


class Timestep(base.Timestep):
        """AMBER trajectory Timestep.

        The Timestep can be initialized with *arg* being

        1. an integer (the number of atoms) and an optional keyword argument *velocities* to allocate
           space for both coordinates and velocities;
        2. another :class:`Timestep` instance, in which case a copy is made (If the copied Timestep
           does not contain velocities but *velocities* = ``True`` is provided, then space for
           velocities is allocated);
        3. a :class:`numpy.ndarray` of shape ``(numatoms, 3)`` (for positions only) or
           ``(numatoms, 6)`` (for positions and velocities): ``positions = arg[:,:3]``,
           ``velocities = arg[:3:6]``.

        """
        # based on TRR Timestep (MDAnalysis.coordinates.xdrfile.TRR.Timestep)
        #
        # NOTE: kwargs is a new thing for Timesteps and not yet in the
        # trajectory API; if this breaks something then we should simply make
        # two different classes, one with the other without velocities. [orbeckst, 2012-05-29]
        def __init__(self, arg, **kwargs):
                velocities = kwargs.pop('velocities', False)
                DIM = 3
                if numpy.dtype(type(arg)) == numpy.dtype(int):
                        self.frame = 0
                        self.step = 0
                        self.time = 0
                        self.numatoms = arg
                        # C floats and C-order for arrays
                        self._pos = numpy.zeros((self.numatoms, DIM), dtype=numpy.float32, order='C')
                        if velocities:
                                self._velocities = numpy.zeros((self.numatoms, DIM), dtype=numpy.float32, order='C')
                        self._unitcell = numpy.zeros(2*DIM, dtype=numpy.float32)  # A,B,C,alpha,beta,gamma
                elif isinstance(arg, Timestep): # Copy constructor
                        # This makes a deepcopy of the timestep
                        self.frame = arg.frame
                        self.numatoms = arg.numatoms
                        self._unitcell = numpy.array(arg._unitcell)
                        self._pos = numpy.array(arg._pos)
                        try:
                                self._velocities = numpy.array(arg._velocities)
                        except AttributeError:
                                pass
                        for attr in ('step', 'time'):
                                if hasattr(arg, attr):
                                        self.__setattr__(attr, arg.__getattribute__(attr))
                elif isinstance(arg, numpy.ndarray):
                        # provide packed array shape == (natoms, 2*DIM)
                        # which contains pos = arg[:,0:3], v = arg[:,3:6]
                        # or just positions: pos = arg[:,0:3] == arg
                        if len(arg.shape) != 2:
                                raise ValueError("packed numpy array (x,v) can only have 2 dimensions")
                        self._unitcell = numpy.zeros(2*DIM, dtype=numpy.float32)
                        self.frame = 0
                        self.step = 0
                        self.time = 0
                        if (arg.shape[0] == 2*DIM and arg.shape[1] != 2*DIM) or \
                                    (arg.shape[0] == DIM and arg.shape[1] != DIM):
                                # wrong order (but need to exclude case where natoms == DIM or natoms == 2*DIM!)
                                raise ValueError("AMBER timestep is to be initialized from (natoms, 2*3) or (natoms, 3) array")
                        self.numatoms = arg.shape[0]
                        self._pos = arg[:,0:DIM].copy('C')                       # C-order
                        if arg.shape[1] == 2*DIM:
                                self._velocities = arg[:,DIM:2*DIM].copy('C')    # C-order
                        elif arg.shape[1] == DIM and velocities:
                                self._velocities = numpy.zeros_like(self._pos)
                        elif arg.shape[1] == DIM:
                                # only positions
                                pass
                        else:
                                raise ValueError("AMBER timestep has no second dimension 3 or 6: shape=%r" % (arg.shape,))
                else:
                        raise ValueError("Cannot create an empty Timestep")
                self._x = self._pos[:,0]
                self._y = self._pos[:,1]
                self._z = self._pos[:,2]

        @property
        def dimensions(self):
                """unitcell dimensions (`A, B, C, alpha, beta, gamma`)

                - `A, B, C` are the lengths of the primitive cell vectors `e1, e2, e3`
                - `alpha` = angle(`e1, e2`)
                - `beta` = angle(`e1, e3`)
                - `gamma` = angle(`e2, e3`)

                .. Note::

                   A :ref:`ASCII AMBER trajectory<ascii-trajectories>` only contains box lengths
                   `A,B,C`; we assume an orthorhombic box and set all
                   angles to 90º.
                """
                # Layout of unitcell is [A,B,C,90,90,90] with the primitive cell vectors
                return self._unitcell

class TRJReader(base.Reader):
        """AMBER trajectory reader.

        Reads the ASCII formatted `AMBER TRJ format`_. Periodic box information
        is auto-detected.

        The number of atoms in a timestep *must* be provided in the `numatoms`
        keyword because it is not stored in the trajectory header and cannot be
        reliably autodetected. The constructor raises a :exc:`ValueError` if
        `numatoms` is left at its default value of ``None``.

        The length of a timestep is not stored in the trajectory itself but can
        be set by passing the `delta` keyword argument to the constructor; it
        is assumed to be in ps. The default value is 1 ps.

        Functionality is currently limited to simple iteration over the
        trajectory.

        .. _AMBER TRJ format: http://ambermd.org/formats.html#trajectory
        """
        format = 'TRJ'
        units = {'time': 'ps', 'length': 'Angstrom'}
        _Timestep = Timestep

        # TODO: implement random access via seek
        #       - compute size of frame
        #       - compute seek offset & go
        #       - check that this works for files >2GB

        def __init__(self, filename, numatoms=None, **kwargs):
                # amber trj REQUIRES the number of atoms from the topology
                if numatoms is None:
                        raise ValueError("AMBER TRJ reader REQUIRES the numatoms keyword")
                self.filename = filename
                self.__numatoms = numatoms
                self.__numframes = None

                self.trjfile = None  # have _read_next_timestep() open it properly!
                self.fixed = 0
                self.skip = 1
                self.skip_timestep = 1   # always 1 for trj at the moment
                self.delta = kwargs.pop("delta", 1.0)   # can set delta manually, default is 1ps
                self.ts = Timestep(self.numatoms)

                # FORMAT(10F8.3)  (X(i), Y(i), Z(i), i=1,NATOM)
                self.default_line_parser = util.FORTRANReader("10F8.3")
                self.lines_per_frame = int(numpy.ceil(3.0 * self.numatoms / len(self.default_line_parser)))
                # The last line per frame might have fewer than 10
                # We determine right away what parser we need for the last
                # line because it will be the same for all frames.
                last_per_line = 3 * self.numatoms % len(self.default_line_parser)
                self.last_line_parser = util.FORTRANReader("%dF8.3" % last_per_line)

                # FORMAT(10F8.3)  BOX(1), BOX(2), BOX(3)
                # is this always on a separate line??
                self.box_line_parser = util.FORTRANReader("3F8.3")

                # Now check for box
                self.periodic = False
                self._detect_amber_box()

                # open file, read first frame
                self._read_next_timestep()

        def _read_next_timestep(self):
                # FORMAT(10F8.3)  (X(i), Y(i), Z(i), i=1,NATOM)
                ts = self.ts
                if self.trjfile is None:
                        self.open_trajectory()

                # Read coordinat frame:
                #coordinates = numpy.zeros(3*self.numatoms, dtype=numpy.float32)
                _coords = []
                for number,line in enumerate(self.trjfile):
                        try:
                                _coords.extend(self.default_line_parser.read(line))
                        except ValueError:
                                # less than 10 entries on the line:
                                _coords.extend(self.last_line_parser.read(line))
                        if number == self.lines_per_frame - 1:
                                # read all atoms that are there in this frame
                                break
                if _coords == []:
                        # at the end of the stream (the loop has not been entered)
                        raise EOFError

                # Read box information
                if self.periodic:
                        line = self.trjfile.next()
                        box = self.box_line_parser.read(line)
                        ts._unitcell[:3] = numpy.array(box, dtype=numpy.float32)
                        ts._unitcell[3:] = [90.,90.,90.]  # assumed

                # probably slow ... could be optimized by storing the coordinates in X,Y,Z
                # lists or directly filling the array; the array/reshape is not good
                # because it creates an intermediate array
                ts._pos[:] = numpy.array(_coords).reshape(self.numatoms, 3)
                ts.frame += 1
                return ts

        def _detect_amber_box(self):
                """Detecting a box in a AMBER trajectory

                Rewind trajectory and check for potential box data
                after the first frame.

                Set :attr:`TRJReader.periodic` to ``True`` if box was
                found, ``False`` otherwise.

                Only run at the beginning as it *rewinds* the trajctory.

                 - see if there's data after the atoms have been read that looks
                   like::

                     FORMAT(10F8.3)  BOX(1), BOX(2), BOX(3)
                     BOX    : size of periodic box

                 - this WILL fail if we have exactly 1 atom in the trajectory because
                   there's no way to distinguish the coordinates from the box
                   so for 1 atom we always assume no box
                 XXX: needs a Timestep that knows about AMBER unitcells!
                """
                if self.numatoms == 1:
                        # for 1 atom we cannot detect the box with the current approach
                        self.periodic = False   # see _read_next_timestep()!
                        wmsg = "Trajectory contains a single atom: assuming periodic=False"
                        warnings.warn(wmsg)
                        return False

                self._reopen()
                self.periodic = False      # make sure that only coordinates are read
                self._read_next_timestep()
                ts = self.ts
                # TODO: what do we do with 1-frame trajectories? Try..except EOFError?
                line = self.trjfile.next()
                nentries = self.default_line_parser.number_of_matches(line)
                if nentries == 3:
                        self.periodic = True
                        ts._unitcell[:3] = self.box_line_parser.read(line)
                        ts._unitcell[3:] = [90.,90.,90.]  # assumed
                else:
                        self.periodic = False
                        ts._unitcell = numpy.zeros(6, numpy.float32)
                self.close()
                return self.periodic

        @property
        def numframes(self):
                """Number of frames (obtained from reading the whole trajectory)."""
                if not self.__numframes is None:   # return cached value
                        return self.__numframes
                try:
                        self.__numframes = self._read_trj_numframes(self.filename)
                except IOError:
                        return 0
                else:
                        return self.__numframes

        def _read_trj_numatoms(self, filename):
                raise NotImplementedError("It is not possible to reliably deduce NATOMS from AMBER trj files")

        def _read_trj_numframes(self, filename):
                self._reopen()
                # the number of lines in the XYZ file will be 2 greater than the number of atoms
                linesPerFrame = self.numatoms * 3. / 10.

                counter = 0
                # step through the file (assuming xyzfile has an iterator)
                for i in self.trjfile:
                        counter = counter + 1
                self.close()
                # need to check this is an integer!
                numframes = int(counter/linesPerFrame)
                return numframes

        @property
        def numatoms(self):
                if not self.__numatoms is None:   # return cached value
                        return self.__numatoms
                try:
                        self.__numatoms = self._read_trj_numatoms(self.filename)
                except IOError:
                        return 0
                else:
                        return self.__numatoms

        def __del__(self):
                if not self.trjfile is None:
                        self.close()

        def __len__(self):
                return self.numframes

        def _reopen(self):
                self.close()
                self.open_trajectory()

        def open_trajectory(self):
                """Open the trajectory for reading and load first frame."""
                self.trjfile, filename = util.anyopen(self.filename, 'r')
                self.header = self.trjfile.readline()  # ignore first line
                if len(self.header.rstrip()) > 80:
                        # Chimera uses this check
                        raise OSError("Header of AMBER formatted trajectory has more than 80 chars. "
                                      "This is probably not a AMBER trajectory.")
                # reset ts
                ts = self.ts
                ts.status = 1
                ts.frame = 0
                ts.step = 0
                ts.time = 0
                return self.trjfile

        def close(self):
                """Close trj trajectory file if it was open."""
                if self.trjfile is None:
                        return
                self.trjfile.close()
                self.trjfile = None

        def rewind(self):
                """Reposition at the beginning of the trajectory"""
                self._reopen()
                self.next()

        def __iter__(self):
                self._reopen()
                self.ts.frame = 0  # start at 0 so that the first frame becomes 1
                while True:
                        try:
                                yield self._read_next_timestep()
                        except EOFError:
                                self.close()
                                raise StopIteration


class NCDFReader(base.Reader):
        """Reader for `AMBER NETCDF format`_ (version 1.0).

        AMBER binary trajectories are automatically recognised by the
        file extension ".ncdf".

        The number of atoms (*numatoms*) does not have to be provided as it can
        be read from the trajectory. The trajectory reader can randomly access
        frames and therefore supports direct indexing (with 0-based frame
        indices) and full-feature trajectory iteration, including slicing.

        Velocities are autodetected and read into the
        :attr:`Timestep._velocities` attribute.

        Periodic unit cell information is detected and used to populate the
        :attr:`Timestep.dimensions` attribute. (If no unit cell is available in
        the trajectory, then :attr:`Timestep.dimensions` will return
        ``[0,0,0,0,0,0]``.)

        Current limitations:

        * only trajectories with time in ps and lengths in Angstroem are processed
        * scale_factors are not supported (and not checked)

        .. _AMBER NETCDF format: http://ambermd.org/netcdf/nctraj.html

        .. SeeAlso:: :class:`NCDFWriter`

        .. versionadded: 0.7.6

        """

        format = 'NCDF'
        version = "1.0"
        units = {'time': 'ps', 'length': 'Angstrom', 'velocity': 'Angstrom/ps'}
        _Timestep = Timestep
        def __init__(self, filename, numatoms=None, **kwargs):
                try:
                        import netCDF4 as netcdf
                except ImportError:
                        logger.fatal("netcdf4-python with the netCDF and HDF5 libraries must be installed for the AMBER ncdf Reader.")
                        logger.fatal("See installation instructions at https://code.google.com/p/mdanalysis/wiki/netcdf")
                        raise ImportError("netCDF4 package missing.\n"
                                          "netcdf4-python with the netCDF and HDF5 libraries must be installed for the AMBER ncdf Reader.\n"
                                          "See installation instructions at https://code.google.com/p/mdanalysis/wiki/netcdf")

                self.filename = filename
                convert_units = kwargs.pop('convert_units', None)
                if convert_units is None:
                        convert_units = MDAnalysis.core.flags['convert_gromacs_lengths']
                        self.convert_units = convert_units  # convert length and time to base units

                self.trjfile = netcdf.Dataset(self.filename)

                if not ('AMBER' in self.trjfile.Conventions.split(',') or 'AMBER' in self.trjfile.Conventions.split()):
                        errmsg = ("NCDF trajectory {0} does not conform to AMBER specifications, "+
                                  "http://ambermd.org/netcdf/nctraj.html ('AMBER' must be one of the tokens "+
                                  "in attribute Conventions)").format(self.filename)
                        logger.fatal(errmsg)
                        raise TypeError(errmsg)
                if not self.trjfile.ConventionVersion == self.version:
                        wmsg = "NCDF trajectory format is %s but the reader implements format %s" % (self.trjfile.ConventionVersion, self.version)
                        warnings.warn(wmsg)
                        logger.warn(wmsg)

                self.numatoms = len(self.trjfile.dimensions['atom'])
                self.numframes = len(self.trjfile.dimensions['frame'])
                # also records time steps in data.variables['time'] and unit
                # but my example only has 0

                try:
                        self.remarks = self.trjfile.title
                except AttributeError:
                        self.remarks = u""
                # other metadata (*= requd):
                # - program*              sander
                # - programVersion*       9.0
                # - application           AMBER
                #

                # checks for not-implemented features (other units would need to be hacked into core.units)
                if self.trjfile.variables['time'].units != "picosecond":
                        raise NotImplementedError("NETCDFReader currently assumes that the trajectory was written with a time unit of picoseconds and not {0}.".format(self.trjfile.variables['time'].units))
                if self.trjfile.variables['coordinates'].units != "angstrom":
                        raise NotImplementedError("NETCDFReader currently assumes that the trajectory was written with a length unit of Angstroem and not {0}.".format(self.trjfile.variables['coordinates'].units))
                if hasattr(self.trjfile.variables['coordinates'], 'scale_factor'):
                        raise NotImplementedError("scale_factors are not implemented")
                if numatoms is not None:
                        if numatoms != self.numatoms:
                                raise ValueError("Supplied numatoms (%d) != natom from ncdf (%d). "
                                                 "Note: numatoms can be None and then the ncdf value is used!" % (numatoms, self.numatoms))

                self.has_velocities =  ('velocities' in self.trjfile.variables)
                self.fixed = 0
                self.skip = 1
                self.skip_timestep = 1   # always 1 for trj at the moment ? CHECK DOCS??
                self.delta = kwargs.pop("delta", 1.0)   # SHOULD GET FROM NCDF (as mean(diff(time)))??
                self.periodic = 'cell_lengths' in self.trjfile.variables
                self._current_frame = 0

                self.ts = self._Timestep(self.numatoms, velocities=self.has_velocities)

                # load first data frame
                self._read_frame(0, self.ts)

        def _read_frame(self, frame, ts):
                if self.trjfile is None:
                        raise IOError("Trajectory is closed")
                if numpy.dtype(type(frame)) != numpy.dtype(int):
                        # convention... for netcdf could also be a slice
                        raise TypeError("frame must be a positive integer")
                if frame >= self.numframes or frame < 0:
                        raise IndexError("frame index must be 0 <= frame < {0}".format(self.numframes))
                # note: self.trjfile.variables['coordinates'].shape == (frames, numatoms, 3)
                ts._pos[:] = self.trjfile.variables['coordinates'][frame]
                ts.time = self.trjfile.variables['time'][frame]
                if self.has_velocities:
                        ts._velocities[:] = self.trjfile.variables['velocities'][frame]
                if self.periodic:
                        ts._unitcell[:3] = self.trjfile.variables['cell_lengths'][frame]
                        ts._unitcell[3:] = self.trjfile.variables['cell_angles'][frame]
                if self.convert_units:
                        self.convert_pos_from_native(ts._pos)                       # in-place !
                        self.convert_time_from_native(ts.time)                      # in-place ! (hope this works...)
                        if self.has_velocities:
                                self.convert_velocities_from_native(ts._velocities) # in-place !
                        if self.periodic:
                                self.convert_pos_from_native(ts._unitcell[:3])      # in-place ! (only lengths)
                ts.frame = frame+1   # frame labels are 1-based
                self._current_frame = frame
                return ts

        def _read_next_timestep(self, ts=None):
                if ts is None:
                        ts = self.ts
                try:
                        return self._read_frame(self._current_frame+1, ts)
                except IndexError:
                        raise IOError


        def __getitem__(self, frame):
                """Return the Timestep corresponding to *frame*.

                If *frame* is a integer then the corresponding frame is
                returned. Negative numbers are counted from the end.

                If frame is a :class:`slice` then an iterator is returned that
                allows iteration over that part of the trajectory.

                .. Note:: *frame* is a 0-based frame index.
                """
                if numpy.dtype(type(frame)) != numpy.dtype(int) and type(frame) != slice:
                        raise TypeError("Can only index NETCDF trajectory with int or a slice.")
                if numpy.dtype(type(frame)) == numpy.dtype(int):
                        if frame < 0:
                                # Interpret similar to a sequence
                                frame = len(self) + frame
                                if frame < 0 or frame >= len(self):
                                        raise IndexError
                        return self._read_frame(frame, self.ts)
                elif type(frame) == slice: # if frame is a slice object
                        if not (((type(frame.start) == int) or (frame.start == None)) and
                                ((type(frame.stop) == int) or (frame.stop == None)) and
                                ((type(frame.step) == int) or (frame.step == None))):
                                raise TypeError("Slice indices are not integers")
                        def iterNETCDF(start=frame.start, stop=frame.stop, step=frame.step):
                                start, stop, step = self._check_slice_indices(start, stop, step)
                                for i in xrange(start, stop, step):
                                        yield self._read_frame(i, self.ts)
                        return iterNETCDF()
                raise ValueError("Type {0} of argument {1} not supported".format(type(i), i))

        def __iter__(self):
                """Iterate over the whole trajectory"""
                for i in xrange(0, self.numframes):
                        try:
                                yield self._read_frame(i, self.ts)
                        except IndexError:
                                raise StopIteration

        def close(self):
                """Close trajectory; any further access will raise an :exc:`IOError`"""
                if not self.trjfile is None:
                        self.trjfile.close()
                        self.trjfile = None

        __del__ = close

        def Writer(self, filename, **kwargs):
                """Returns a NCDFWriter for *filename* with the same parameters as this NCDF.

                All values can be changed through keyword arguments.

                :Arguments:
                  *filename*
                      filename of the output NCDF trajectory
                :Keywords:
                  *numatoms*
                      number of atoms
                  *delta*
                      length of one timestep in picoseconds
                  *remarks*
                      string that is stored in the title field

                :Returns: :class:`NCDFWriter`
                """
                numatoms = kwargs.pop('numatoms', self.numatoms)
                kwargs.setdefault('remarks', self.remarks)
                kwargs.setdefault('delta', self.delta)
                return NCDFWriter(filename, numatoms, **kwargs)


class NCDFWriter(base.Writer):
        """Writer for `AMBER NETCDF format`_ (version 1.0).

        AMBER binary trajectories are automatically recognised by the
        file extension ".ncdf".

        Velocities are written out if they are detected in the input
        :class:`Timestep`. The trajectories are always written with ångström
        for the lengths and picoseconds for the time (and hence Å/ps for
        velocities).

        Unit cell information is written if available.

        .. _AMBER NETCDF format: http://ambermd.org/netcdf/nctraj.html

        .. SeeAlso:: :class:`NCDFReader`

        .. versionadded: 0.7.6

        """

        format = 'NCDF'
        version = "1.0"
        units = {'time': 'ps', 'length': 'Angstrom', 'velocity': 'Angstrom/ps'}

        def __init__(self, filename, numatoms, start=0, step=1, delta=1.0, remarks=None,
                     convert_units=None, zlib=False, cmplevel=1):
                '''Create a new NCDFWriter

                :Arguments:
                 *filename*
                    name of output file
                 *numatoms*
                    number of atoms in trajectory file

                :Keywords:
                  *start*
                    starting timestep
                  *step*
                    skip between subsequent timesteps
                  *delta*
                    timestep
                  *convert_units*
                    ``True``: units are converted to the AMBER base format; ``None`` selects
                    the value of :data:`MDAnalysis.core.flags` ['convert_gromacs_lengths'].
                    (see :ref:`flags-label`)
                  *zlib*
                    compress data [``False``]
                  *cmplevel*
                    compression level (1-9) [1]
                '''
                self.filename = filename
                if numatoms == 0:
                        raise ValueError("NCDFWriter: no atoms in output trajectory")
                self.numatoms = numatoms
                if convert_units is None:
                        convert_units = MDAnalysis.core.flags['convert_gromacs_lengths']
                self.convert_units = convert_units    # convert length and time to base units on the fly?

                self.start = start    # do we use those?
                self.step = step      # do we use those?
                self.delta = delta
                self.remarks = remarks or u"AMBER NetCDF format (MDAnalysis.coordinates.trj.NCDFWriter)"

                self.zlib = zlib
                self.cmplevel = cmplevel

                self.ts = None                # when/why would this be assigned??
                self.__first_frame = True     # signals to open trajectory
                self.trjfile = None           # open on first write with _init_netcdf()
                self.periodic = None          # detect on first write
                self.has_velocities = False   # velocities disabled for the moment
                self.curr_frame = 0

        def _init_netcdf(self, periodic=True, velocities=False):
                """Initialize netcdf AMBER 1.0 trajectory.

                The trajectory is opened when the first frame is written
                because that is the earlies time that we can detect if the
                output should contain periodicity information (i.e. the unit
                cell dimensions).

                Based on Joshua Adelman's `netcdf4storage.py`_ in `Issue 109`_.

                .. _`Issue 109`:
                   http://code.google.com/p/mdanalysis/issues/detail?id=109#c2
                .. _`netcdf4storage.py`:
                   http://code.google.com/p/mdanalysis/issues/attachmentText?id=109&aid=1090002000&name=netcdf4storage.py
                """
                try:
                        import netCDF4 as netcdf
                except ImportError:
                        logger.fatal("netcdf4-python with the netCDF and HDF5 libraries must be installed for the AMBER ncdf Writer.")
                        logger.fatal("See installation instructions at https://code.google.com/p/mdanalysis/wiki/netcdf")
                        raise ImportError("netCDF4 package missing.\n"
                                          "netcdf4-python with the netCDF and HDF5 libraries must be installed for the AMBER ncdf Writer.\n"
                                          "See installation instructions at https://code.google.com/p/mdanalysis/wiki/netcdf")

                if not self.__first_frame:
                        raise IOError(errno.EIO, "Attempt to write to closed file {0}".format(self.filename))

                ncfile = netcdf.Dataset(self.filename, clobber=True, mode='w', format='NETCDF3_64BIT')

                # Set global attributes.
                setattr(ncfile, 'program', 'MDAnalysis.coordinates.TRJ.NCDFWriter')
                setattr(ncfile, 'programVersion', MDAnalysis.__version__)
                setattr(ncfile, 'Conventions', 'AMBER')
                setattr(ncfile, 'ConventionVersion', '1.0')
                setattr(ncfile, 'application', 'MDAnalysis')

                # Create dimensions
                ncfile.createDimension('frame', None)         # unlimited number of steps (can append)
                ncfile.createDimension('atom', self.numatoms) # number of atoms in system
                ncfile.createDimension('spatial', 3)          # number of spatial dimensions
                ncfile.createDimension('cell_spatial', 3)     # unitcell lengths
                ncfile.createDimension('cell_angular', 3)     # unitcell angles

                # Create variables.
                coords = ncfile.createVariable('coordinates','f8', ('frame','atom','spatial'),
                                                     zlib=self.zlib, complevel=self.cmplevel)
                setattr(coords, 'units', 'angstrom')
                time = ncfile.createVariable('time','f8', ('frame',),
                                                     zlib=self.zlib, complevel=self.cmplevel)
                setattr(time, 'units', 'picosecond')

                self.periodic = periodic
                if self.periodic:
                        cell_lengths = ncfile.createVariable('cell_lengths','f8', ('frame', 'cell_spatial'),
                                                             zlib=self.zlib, complevel=self.cmplevel)
                        setattr(cell_lengths, 'units', 'angstrom')
                        cell_angles = ncfile.createVariable('cell_angles','f8', ('frame', 'cell_angular'),
                                                            zlib=self.zlib, complevel=self.cmplevel)
                        setattr(cell_angles, 'units', 'degrees')

                self.has_velocities = velocities
                if self.has_velocities:
                        velocs = ncfile.createVariable('velocities','f8', ('frame','atom','spatial'),
                                                             zlib=self.zlib, complevel=self.cmplevel)
                        setattr(velocs, 'units', 'angstrom/picosecond')

                ncfile.sync()
                self.__first_frame = False
                self.trjfile = ncfile

        def is_periodic(self, ts=None):
                """Return ``True`` if :class:`Timestep` *ts* contains a valid simulation box"""
                ts = ts if ts is not None else self.ts
                return numpy.all(ts.dimensions > 0)

        def write_next_timestep(self, ts=None):
                '''write a new timestep to the trj file

                *ts* is a :class:`Timestep` instance containing coordinates to
                be written to trajectory file
                '''
                if ts is None:
                        if not hasattr(self, "ts") or self.ts is None:
                                raise IOError("NCDFWriter: no coordinate data to write to trajectory file")
                        else:
                                ts = self.ts   # self.ts would have to be assigned manually!
                elif ts.numatoms != self.numatoms:
                        raise IOError("NCDFWriter: Timestep does not have the correct number of atoms")

                if self.trjfile is None:
                        # first time step: analyze data and open trajectory accordingly
                        # (also sets self.periodic and self.has_velocities)
                        self._init_netcdf(periodic=self.is_periodic(ts), velocities=hasattr(ts, "_velocities"))

                return self._write_next_timestep(ts)


        def _write_next_timestep(self, ts):
                """Write coordinates and unitcell information to NCDF file.

                Do not call this method directly; instead use
                :meth:`write_next_timestep` because some essential setup is done
                there before writing the first frame.

                Based on Joshua Adelman's `netcdf4storage.py`_ in `Issue 109`_.

                .. _`Issue 109`:
                   http://code.google.com/p/mdanalysis/issues/detail?id=109#c2
                .. _`netcdf4storage.py`:
                   http://code.google.com/p/mdanalysis/issues/attachmentText?id=109&aid=1090002000&name=netcdf4storage.py
                """
                assert self.trjfile is not None, "trjfile must be open in order to write to it"

                if self.convert_units:
                        # make a copy of the scaled positions so that the in-memory
                        # timestep is not changed (would have lead to wrong results if
                        # analysed *after* writing a time step to disk). The new
                        # implementation could lead to memory problems and/or slow-down for
                        # very big systems because we temporarily create a new array pos
                        # for each frame written
                        pos = self.convert_pos_to_native(ts._pos, inplace=False)
                        try:
                                time = self.convert_time_to_native(ts.time, inplace=False)
                        except AttributeError:
                                time = ts.frame * self.convert_time_to_native(self.delta, inplace=False)
                        unitcell = self.convert_dimensions_to_unitcell(ts)
                else:
                        pos = ts._pos
                        try:
                                time = ts.time
                        except AttributeError:
                                time = ts.frame * self.delta
                        unitcell = ts.dimensions

                if not hasattr(ts, 'step'):
                        # bogus, should be actual MD step number, i.e. frame * delta/dt
                        ts.step = ts.frame

                # write step
                self.trjfile.variables['coordinates'][self.curr_frame,:,:] = pos
                self.trjfile.variables['time'][self.curr_frame] = time
                if self.periodic:
                        self.trjfile.variables['cell_lengths'][self.curr_frame,:] = unitcell[:3]
                        self.trjfile.variables['cell_angles'][self.curr_frame,:] = unitcell[3:]
                if self.has_velocities:
                        if self.convert_units:
                                velocities = self.convert_velocities_to_native(ts._velocities, inplace=False)
                        else:
                                velocities = ts._velocities
                        self.trjfile.variables['velocities'][self.curr_frame,:,:] = velocities
                self.trjfile.sync()
                self.curr_frame += 1

        def close(self):
                if not self.trjfile is None:
                        self.trjfile.close()
                        self.trjfile = None

        __del__ = close
