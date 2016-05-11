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

      coordinates of the atoms as a :class:`numpy.ndarray` of shape `(n_atoms, 3)`

   .. attribute:: _velocities

      velocities of the atoms as a :class:`numpy.ndarray` of shape `(n_atoms, 3)`;
      only available if the trajectory contains velocities or if the
      *velocities* = ``True`` keyword has been supplied.

   .. attribute:: _forces

      forces of the atoms as a :class:`numpy.ndarray` of shape `(n_atoms, 3)`;
      only available if the trajectory contains forces or if the
      *forces* = ``True`` keyword has been supplied.


.. _ascii-trajectories:

ASCII TRAJ trajectories
-----------------------

ASCII AMBER_ TRJ coordinate files (as defined in `AMBER TRJ format`_)
are handled by the :class:`TRJReader`. It is also possible to directly
read *bzip2* or *gzip* compressed files.

AMBER ASCII trajectories are recognised by the suffix '.trj',
'.mdcrd' or '.crdbox (possibly with an additional '.gz' or '.bz2').

.. rubric:: Limitations

* Periodic boxes are only stored as box lengths A, B, C in an AMBER
  trajectory; the reader always assumes that these are orthorhombic
  boxes.

* The trajectory does not contain time information so we simply set
  the time step to 1 ps (or the user could provide it as kwarg *dt*)

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

Binary trajectories can also contain velocities and forces, and can record the
exact time
step. In principle, the trajectories can be in different units than the AMBER
defaults of ångström and picoseconds but at the moment MDAnalysis only supports
those and will raise a :exc:`NotImplementedError` if anything else is detected.

.. autoclass:: NCDFReader
   :members:

.. autoclass:: NCDFWriter
   :members:


.. Links

.. _AMBER: http://ambermd.org
.. _AMBER TRJ format: http://ambermd.org/formats.html#trajectory
.. _AMBER netcdf format: http://ambermd.org/netcdf/nctraj.html
.. _AMBER netcdf: http://ambermd.org/netcdf/nctraj.html
.. _NetCDF: http://www.unidata.ucar.edu/software/netcdf
.. _Issue Tracker: https://github.com/MDAnalysis/mdanalysis/issues
.. _MDAnalysis mailinglist: http://groups.google.com/group/mdnalysis-discussion



"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
import warnings
import errno
import logging

import MDAnalysis
from ..core import flags
from . import base
from ..lib import util

logger = logging.getLogger("MDAnalysis.coordinates.AMBER")

try:
    import netCDF4 as netcdf
except ImportError:
    # Just to notify the user; the module will still load. However, NCDFReader
    # and NCDFWriter will raise a proper ImportError if they are called without
    # the netCDF4 library present. See Issue 122 for a discussion.
    logger.debug(
        "Failed to import netCDF4; AMBER NETCDFReader/Writer will not work. "
        "Install netCDF4 from https://github.com/Unidata/netcdf4-python.")
    logger.debug(
        "See also https://github.com/MDAnalysis/mdanalysis/wiki/netcdf")


class Timestep(base.Timestep):
    """AMBER trajectory Timestep.

    The Timestep can be initialized with *arg* being an integer
    (the number of atoms) and an optional keyword argument *velocities* to
    allocate space for both coordinates and velocities;

    .. versionchanged:: 0.10.0
       Added ability to contain Forces
    """
    order = 'C'


class TRJReader(base.Reader):
    """AMBER trajectory reader.

    Reads the ASCII formatted `AMBER TRJ format`_. Periodic box information
    is auto-detected.

    The number of atoms in a timestep *must* be provided in the `n_atoms`
    keyword because it is not stored in the trajectory header and cannot be
    reliably autodetected. The constructor raises a :exc:`ValueError` if
    `n_atoms` is left at its default value of ``None``.

    The length of a timestep is not stored in the trajectory itself but can
    be set by passing the `dt` keyword argument to the constructor; it
    is assumed to be in ps. The default value is 1 ps.

    Functionality is currently limited to simple iteration over the
    trajectory.

    .. _AMBER TRJ format: http://ambermd.org/formats.html#trajectory

    .. versionchanged:: 0.11.0
       Frames now 0-based instead of 1-based
       kwarg 'delta' renamed to 'dt', for uniformity with other Readers
    """
    format = ['TRJ', 'MDCRD', 'CRDBOX']
    units = {'time': 'ps', 'length': 'Angstrom'}
    _Timestep = Timestep

    def __init__(self, filename, n_atoms=None, **kwargs):
        super(TRJReader, self).__init__(filename, **kwargs)
        if n_atoms is None:
            raise ValueError("AMBER TRJ reader REQUIRES the n_atoms keyword")
        self._n_atoms = n_atoms
        self._n_frames = None

        self.trjfile = None  # have _read_next_timestep() open it properly!
        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)

        # FORMAT(10F8.3)  (X(i), Y(i), Z(i), i=1,NATOM)
        self.default_line_parser = util.FORTRANReader("10F8.3")
        self.lines_per_frame = int(np.ceil(3.0 * self.n_atoms / len(
            self.default_line_parser)))
        # The last line per frame might have fewer than 10
        # We determine right away what parser we need for the last
        # line because it will be the same for all frames.
        last_per_line = 3 * self.n_atoms % len(self.default_line_parser)
        self.last_line_parser = util.FORTRANReader("{0:d}F8.3".format(
            last_per_line))

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
        # coordinates = numpy.zeros(3*self.n_atoms, dtype=np.float32)
        _coords = []
        for number, line in enumerate(self.trjfile):
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
            line = next(self.trjfile)
            box = self.box_line_parser.read(line)
            ts._unitcell[:3] = np.array(box, dtype=np.float32)
            ts._unitcell[3:] = [90., 90., 90.]  # assumed

        # probably slow ... could be optimized by storing the coordinates in
        # X,Y,Z lists or directly filling the array; the array/reshape is not
        # good because it creates an intermediate array
        ts._pos[:] = np.array(_coords).reshape(self.n_atoms, 3)
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
        if self.n_atoms == 1:
            # for 1 atom we cannot detect the box with the current approach
            self.periodic = False  # see _read_next_timestep()!
            wmsg = "Trajectory contains a single atom: assuming periodic=False"
            warnings.warn(wmsg)
            return False

        self._reopen()
        self.periodic = False  # make sure that only coordinates are read
        self._read_next_timestep()
        ts = self.ts
        # TODO: what do we do with 1-frame trajectories? Try..except EOFError?
        line = next(self.trjfile)
        nentries = self.default_line_parser.number_of_matches(line)
        if nentries == 3:
            self.periodic = True
            ts._unitcell[:3] = self.box_line_parser.read(line)
            ts._unitcell[3:] = [90., 90., 90.]  # assumed
        else:
            self.periodic = False
            ts._unitcell = np.zeros(6, np.float32)
        self.close()
        return self.periodic

    @property
    def n_frames(self):
        """Number of frames (obtained from reading the whole trajectory)."""
        if self._n_frames is not None:  # return cached value
            return self._n_frames
        try:
            self._n_frames = self._read_trj_n_frames(self.filename)
        except IOError:
            return 0
        else:
            return self._n_frames

    def _read_trj_n_frames(self, filename):
        self._reopen()

        counter = 0
        try:
            while True:
                next(self)
                counter += 1
        except EOFError:
            self.rewind()

        return counter

    @property
    def n_atoms(self):
        return self._n_atoms

    def _reopen(self):
        self.close()
        self.open_trajectory()

    def open_trajectory(self):
        """Open the trajectory for reading and load first frame."""
        self.trjfile = util.anyopen(self.filename, 'r')
        self.header = self.trjfile.readline()  # ignore first line
        if len(self.header.rstrip()) > 80:
            # Chimera uses this check
            raise OSError(
                "Header of AMBER formatted trajectory has more than 80 chars. "
                "This is probably not a AMBER trajectory.")
        # reset ts
        ts = self.ts
        ts.frame = -1

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
        next(self)


class NCDFReader(base.Reader):
    """Reader for `AMBER NETCDF format`_ (version 1.0).

    AMBER binary trajectories are automatically recognised by the
    file extension ".ncdf".

    The number of atoms (*n_atoms*) does not have to be provided as it can
    be read from the trajectory. The trajectory reader can randomly access
    frames and therefore supports direct indexing (with 0-based frame
    indices) and full-feature trajectory iteration, including slicing.

    Velocities are autodetected and read into the
    :attr:`Timestep._velocities` attribute.

    Forces are autodetected and read into the
    :attr:`Timestep._forces` attribute.

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
    .. versionchanged:: 0.10.0
       Added ability to read Forces
    .. versionchanged:: 0.11.0
       Frame labels now 0-based instead of 1-based
       kwarg 'delta' renamed to 'dt', for uniformity with other Readers
    """

    format = ['NCDF', 'NC']
    multiframe = True
    version = "1.0"
    units = {'time': 'ps',
             'length': 'Angstrom',
             'velocity': 'Angstrom/ps',
             'force': 'kcal/(mol*Angstrom)'}
    _Timestep = Timestep

    def __init__(self, filename, n_atoms=None, **kwargs):
        try:
            import netCDF4 as netcdf
        except ImportError:
            errmsg = ("netCDF4 package missing.\n"
                      "netcdf4-python with the netCDF and HDF5 libraries must"
                      " be installed for the AMBER ncdf Reader.\n"
                      "See installation instructions at "
                      "https://github.com/MDAnalysis/mdanalysis/wiki/netcdf")
            logger.fatal(errmsg)
            raise ImportError(errmsg)

        super(NCDFReader, self).__init__(filename, **kwargs)

        self.trjfile = netcdf.Dataset(self.filename)

        if not ('AMBER' in self.trjfile.Conventions.split(',') or
                'AMBER' in self.trjfile.Conventions.split()):
            errmsg = ("NCDF trajectory {0} does not conform to AMBER "
                      "specifications, http://ambermd.org/netcdf/nctraj.html "
                      "('AMBER' must be one of the tokens in attribute "
                      "Conventions)".format(self.filename))
            logger.fatal(errmsg)
            raise TypeError(errmsg)
        if not self.trjfile.ConventionVersion == self.version:
            wmsg = ("NCDF trajectory format is {0!s} but the reader "
                    "implements format {1!s}".format(
                        self.trjfile.ConventionVersion, self.version))
            warnings.warn(wmsg)
            logger.warn(wmsg)

        self.n_atoms = len(self.trjfile.dimensions['atom'])
        self.n_frames = len(self.trjfile.dimensions['frame'])
        # also records time steps in data.variables['time'] and unit
        # but my example only has 0

        try:
            self.remarks = self.trjfile.title
        except AttributeError:
            self.remarks = ""
        # other metadata (*= requd):
        # - program*              sander
        # - programVersion*       9.0
        # - application           AMBER
        #

        # checks for not-implemented features (other units would need to be
        # hacked into MDAnalysis.units)
        if self.trjfile.variables['time'].units != "picosecond":
            raise NotImplementedError(
                "NETCDFReader currently assumes that the trajectory was "
                "written with a time unit of picoseconds and "
                "not {0}.".format(self.trjfile.variables['time'].units))
        if self.trjfile.variables['coordinates'].units != "angstrom":
            raise NotImplementedError(
                "NETCDFReader currently assumes that the trajectory was "
                "written with a length unit of Angstroem and "
                "not {0}.".format(self.trjfile.variables['coordinates'].units))
        if hasattr(self.trjfile.variables['coordinates'], 'scale_factor'):
            raise NotImplementedError("scale_factors are not implemented")
        if n_atoms is not None:
            if n_atoms != self.n_atoms:
                raise ValueError(
                    "Supplied n_atoms ({0}) != natom from ncdf ({1}). Note: "
                    "n_atoms can be None and then the ncdf value is used!"
                    "".format(n_atoms, self.n_atoms))

        self.has_velocities = 'velocities' in self.trjfile.variables
        self.has_forces = 'forces' in self.trjfile.variables

        self.periodic = 'cell_lengths' in self.trjfile.variables
        self._current_frame = 0

        self.ts = self._Timestep(self.n_atoms,
                                 velocities=self.has_velocities,
                                 forces=self.has_forces,
                                 reader=self,  # for dt
                                 **self._ts_kwargs)

        # load first data frame
        self._read_frame(0)

    def _read_frame(self, frame):
        ts = self.ts

        if self.trjfile is None:
            raise IOError("Trajectory is closed")
        if frame >= self.n_frames or frame < 0:
            raise IndexError("frame index must be 0 <= frame < {0}"
                             "".format(self.n_frames))
        # note: self.trjfile.variables['coordinates'].shape == (frames, n_atoms, 3)
        ts._pos[:] = self.trjfile.variables['coordinates'][frame]
        ts.time = self.trjfile.variables['time'][frame]
        if self.has_velocities:
            ts._velocities[:] = self.trjfile.variables['velocities'][frame]
        if self.has_forces:
            ts._forces[:] = self.trjfile.variables['forces'][frame]
        if self.periodic:
            ts._unitcell[:3] = self.trjfile.variables['cell_lengths'][frame]
            ts._unitcell[3:] = self.trjfile.variables['cell_angles'][frame]
        if self.convert_units:
            self.convert_pos_from_native(ts._pos)  # in-place !
            self.convert_time_from_native(
                ts.time)  # in-place ! (hope this works...)
            if self.has_velocities:
                self.convert_velocities_from_native(ts._velocities,
                                                    inplace=True)
            if self.has_forces:
                self.convert_forces_from_native(ts._forces, inplace=True)
            if self.periodic:
                self.convert_pos_from_native(
                    ts._unitcell[:3])  # in-place ! (only lengths)
        ts.frame = frame  # frame labels are 0-based
        self._current_frame = frame
        return ts

    def _reopen(self):
        self._current_frame = -1

    def _read_next_timestep(self, ts=None):
        if ts is None:
            ts = self.ts
        try:
            return self._read_frame(self._current_frame + 1)
        except IndexError:
            raise IOError

    def _get_dt(self):
        t1 = self.trjfile.variables['time'][1]
        t0 = self.trjfile.variables['time'][0]
        return t1 - t0

    def close(self):
        """Close trajectory; any further access will raise an :exc:`IOError`"""
        if self.trjfile is not None:
            self.trjfile.close()
            self.trjfile = None

    def Writer(self, filename, **kwargs):
        """Returns a NCDFWriter for *filename* with the same parameters as this NCDF.

        All values can be changed through keyword arguments.

        :Arguments:
          *filename*
              filename of the output NCDF trajectory
        :Keywords:
          *n_atoms*
              number of atoms
          *dt*
              length of one timestep in picoseconds
          *remarks*
              string that is stored in the title field

        :Returns: :class:`NCDFWriter`
        """
        n_atoms = kwargs.pop('n_atoms', self.n_atoms)
        kwargs.setdefault('remarks', self.remarks)
        kwargs.setdefault('dt', self.dt)
        return NCDFWriter(filename, n_atoms, **kwargs)


class NCDFWriter(base.Writer):
    """Writer for `AMBER NETCDF format`_ (version 1.0).

    AMBER binary trajectories are automatically recognised by the
    file extension ".ncdf" or ".nc".

    Velocities are written out if they are detected in the input
    :class:`Timestep`. The trajectories are always written with ångström
    for the lengths and picoseconds for the time (and hence Å/ps for
    velocities).

    Unit cell information is written if available.

    .. _AMBER NETCDF format: http://ambermd.org/netcdf/nctraj.html

    .. SeeAlso:: :class:`NCDFReader`

    .. versionadded: 0.7.6

    .. versionchanged:: 0.10.0
       Added ability to write velocities and forces
    .. versionchanged:: 0.11.0
       kwarg 'delta' renamed to 'dt', for uniformity with other Readers
    """

    format = 'NCDF'
    version = "1.0"
    units = {'time': 'ps',
             'length': 'Angstrom',
             'velocity': 'Angstrom/ps',
             'force': 'kcal/(mol*Angstrom)'}

    def __init__(self,
                 filename,
                 n_atoms,
                 start=0,
                 step=1,
                 dt=1.0,
                 remarks=None,
                 convert_units=None,
                 zlib=False,
                 cmplevel=1,
                 **kwargs):
        """Create a new NCDFWriter

        :Arguments:
         *filename*
            name of output file
         *n_atoms*
            number of atoms in trajectory file

        :Keywords:
          *start*
            starting timestep
          *step*
            skip between subsequent timesteps
          *dt*
            timestep
          *convert_units*
            ``True``: units are converted to the AMBER base format; ``None``
            selects the value of :data:`MDAnalysis.core.flags`
            ['convert_lengths'] (see :ref:`flags-label`).
          *zlib*
            compress data [``False``]
          *cmplevel*
            compression level (1-9) [1]
          *velocities*
            Write velocities into the trajectory [``False``]
          *forces*
            Write forces into the trajectory [``False``]
        """
        self.filename = filename
        if n_atoms == 0:
            raise ValueError("NCDFWriter: no atoms in output trajectory")
        self.n_atoms = n_atoms
        if convert_units is None:
            convert_units = flags['convert_lengths']
        # convert length and time to base units on the fly?
        self.convert_units = convert_units

        self.start = start  # do we use those?
        self.step = step  # do we use those?
        self.dt = dt
        self.remarks = remarks or "AMBER NetCDF format (MDAnalysis.coordinates.trj.NCDFWriter)"

        self.zlib = zlib
        self.cmplevel = cmplevel

        self.ts = None  # when/why would this be assigned??
        self._first_frame = True  # signals to open trajectory
        self.trjfile = None  # open on first write with _init_netcdf()
        self.periodic = None  # detect on first write
        self.has_velocities = kwargs.get('velocities', False)
        self.has_forces = kwargs.get('forces', False)
        self.curr_frame = 0

    def _init_netcdf(self, periodic=True):
        """Initialize netcdf AMBER 1.0 trajectory.

        The trajectory is opened when the first frame is written
        because that this is the earliest time that we can detect if the
        output should contain periodicity information (i.e. the unit
        cell dimensions).

        Based on Joshua Adelman's `netcdf4storage.py`_ in `Issue 109`_.

        .. _`Issue 109`:
           https://github.com/MDAnalysis/mdanalysis/issues/109
        .. _`netcdf4storage.py`:
           https://storage.googleapis.com/google-code-attachments/mdanalysis/issue-109/comment-2/netcdf4storage.py
        """
        try:
            import netCDF4 as netcdf
        except ImportError:
            logger.fatal("netcdf4-python with the netCDF and HDF5 libraries "
                         "must be installed for the AMBER ncdf Writer."
                         "See installation instructions at "
                         "https://github.com/MDAnalysis/mdanalysis/wiki/netcdf")
            raise ImportError(
                "netCDF4 package missing.\n"
                "netcdf4-python with the netCDF and HDF5 libraries must be "
                "installed for the AMBER ncdf Writer.\n"
                "See installation instructions at "
                "https://github.com/MDAnalysis/mdanalysis/wiki/netcdf")

        if not self._first_frame:
            raise IOError(
                errno.EIO,
                "Attempt to write to closed file {0}".format(self.filename))

        ncfile = netcdf.Dataset(self.filename,
                                clobber=True,
                                mode='w',
                                format='NETCDF3_64BIT')

        # Set global attributes.
        setattr(ncfile, 'program', 'MDAnalysis.coordinates.TRJ.NCDFWriter')
        setattr(ncfile, 'programVersion', MDAnalysis.__version__)
        setattr(ncfile, 'Conventions', 'AMBER')
        setattr(ncfile, 'ConventionVersion', '1.0')
        setattr(ncfile, 'application', 'MDAnalysis')

        # Create dimensions
        ncfile.createDimension('frame',
                               None)  # unlimited number of steps (can append)
        ncfile.createDimension('atom',
                               self.n_atoms)  # number of atoms in system
        ncfile.createDimension('spatial', 3)  # number of spatial dimensions
        ncfile.createDimension('cell_spatial', 3)  # unitcell lengths
        ncfile.createDimension('cell_angular', 3)  # unitcell angles
        ncfile.createDimension('label', 5)  # needed for cell_angular

        # Create variables.
        coords = ncfile.createVariable('coordinates',
                                       'f4', ('frame', 'atom', 'spatial'),
                                       zlib=self.zlib,
                                       complevel=self.cmplevel)
        setattr(coords, 'units', 'angstrom')

        spatial = ncfile.createVariable('spatial', 'c', ('spatial', ))
        spatial[:] = np.asarray(list('xyz'))

        time = ncfile.createVariable('time',
                                     'f4', ('frame', ),
                                     zlib=self.zlib,
                                     complevel=self.cmplevel)
        setattr(time, 'units', 'picosecond')

        self.periodic = periodic
        if self.periodic:
            cell_lengths = ncfile.createVariable('cell_lengths',
                                                 'f8',
                                                 ('frame', 'cell_spatial'),
                                                 zlib=self.zlib,
                                                 complevel=self.cmplevel)
            setattr(cell_lengths, 'units', 'angstrom')

            cell_spatial = ncfile.createVariable('cell_spatial', 'c',
                                                 ('cell_spatial', ))
            cell_spatial[:] = np.asarray(list('abc'))

            cell_angles = ncfile.createVariable('cell_angles',
                                                'f8',
                                                ('frame', 'cell_angular'),
                                                zlib=self.zlib,
                                                complevel=self.cmplevel)
            setattr(cell_angles, 'units', 'degrees')

            cell_angular = ncfile.createVariable('cell_angular', 'c',
                                                 ('cell_angular', 'label'))
            cell_angular[:] = np.asarray([list('alpha'), list('beta '), list(
                'gamma')])

        # These properties are optional, and are specified on Writer creation
        if self.has_velocities:
            velocs = ncfile.createVariable('velocities',
                                           'f8', ('frame', 'atom', 'spatial'),
                                           zlib=self.zlib,
                                           complevel=self.cmplevel)
            setattr(velocs, 'units', 'angstrom/picosecond')
        if self.has_forces:
            forces = ncfile.createVariable('forces',
                                           'f8', ('frame', 'atom', 'spatial'),
                                           zlib=self.zlib,
                                           complevel=self.cmplevel)
            setattr(forces, 'units', 'kilocalorie/mole/angstrom')

        ncfile.sync()
        self._first_frame = False
        self.trjfile = ncfile

    def is_periodic(self, ts=None):
        """Return ``True`` if :class:`Timestep` *ts* contains a valid
        simulation box
        """
        ts = ts if ts is not None else self.ts
        return np.all(ts.dimensions > 0)

    def write_next_timestep(self, ts=None):
        '''write a new timestep to the trj file

        *ts* is a :class:`Timestep` instance containing coordinates to
        be written to trajectory file
        '''
        if ts is None:
            if not hasattr(self, "ts") or self.ts is None:
                raise IOError(
                    "NCDFWriter: no coordinate data to write to trajectory file")
            else:
                ts = self.ts  # self.ts would have to be assigned manually!
        elif ts.n_atoms != self.n_atoms:
            raise IOError(
                "NCDFWriter: Timestep does not have the correct number of atoms")

        if self.trjfile is None:
            # first time step: analyze data and open trajectory accordingly
            self._init_netcdf(periodic=self.is_periodic(ts))

        return self._write_next_timestep(ts)

    def _write_next_timestep(self, ts):
        """Write coordinates and unitcell information to NCDF file.

        Do not call this method directly; instead use
        :meth:`write_next_timestep` because some essential setup is done
        there before writing the first frame.

        Based on Joshua Adelman's `netcdf4storage.py`_ in `Issue 109`_.

        .. _`Issue 109`:
           https://github.com/MDAnalysis/mdanalysis/issues/109
        .. _`netcdf4storage.py`:
           https://storage.googleapis.com/google-code-attachments/mdanalysis/issue-109/comment-2/netcdf4storage.py
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
                time = ts.frame * self.convert_time_to_native(self.dt,
                                                              inplace=False)
            unitcell = self.convert_dimensions_to_unitcell(ts)
        else:
            pos = ts._pos
            time = ts.time

            unitcell = ts.dimensions

        # write step
        self.trjfile.variables['coordinates'][self.curr_frame, :, :] = pos
        self.trjfile.variables['time'][self.curr_frame] = time
        if self.periodic:
            self.trjfile.variables['cell_lengths'][
                self.curr_frame, :] = unitcell[:3]
            self.trjfile.variables['cell_angles'][
                self.curr_frame, :] = unitcell[3:]
        if self.has_velocities:
            if self.convert_units:
                velocities = self.convert_velocities_to_native(ts._velocities,
                                                               inplace=False)
            else:
                velocities = ts._velocities
            self.trjfile.variables['velocities'][
                self.curr_frame, :, :] = velocities
        if self.has_forces:
            if self.convert_units:
                forces = self.convert_forces_to_native(ts._forces,
                                                       inplace=False)
            else:
                forces = ts._velocities
            self.trjfile.variables['forces'][self.curr_frame, :, :] = forces
        self.trjfile.sync()
        self.curr_frame += 1

    def close(self):
        if self.trjfile is not None:
            self.trjfile.close()
            self.trjfile = None
