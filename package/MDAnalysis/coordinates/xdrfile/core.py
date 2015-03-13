# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://mdanalysis.googlecode.com
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
Common high-level Gromacs XDR functionality --- :mod:`MDAnalysis.coordinates.xdrfile.core`
==========================================================================================

The :mod:`MDAnalysis.coordinates.xdrfile.core` module contains generic
classes to access Gromacs_ XDR-encoded trajectory formats such as TRR
and XTC.

A generic Gromacs_ trajectory is simply called "trj" within this
module.

.. SeeAlso:: :mod:`MDAnalysis.coordinates.base` for the generic MDAnalysis base
             classes and :mod:`MDAnalysis.coordinates.xdrfile.libxdrfile2` for
             the low-level bindings to the XDR trajectories.

.. _Gromacs: http://www.gromacs.org

Generic xdr trj classes
-----------------------

The generic classes are subclassed to generate the specific classes
for the XTC and TRR format.

.. versionchanged:: 0.8.0
   The XTC/TRR I/O interface now uses
   :mod:`~MDAnalysis.coordinates.xdrfile.libxdrfile2`, which has seeking and
   indexing capabilities. Note that unlike
   :mod:`~MDAnalysis.coordinates.xdrfile.libxdrfile` before it,
   :mod:`~MDAnalysis.coordinates.xdrfile.libxdrfile2` is distributed under the
   GNU GENERAL PUBLIC LICENSE, version 2 (or higher).

.. versionchanged:: 0.9.0
   TrjReader now stores the offsets used for frame seeking automatically as a
   hidden file in the same directory as the source trajectory. These offsets
   are automatically retrieved upon TrjReader instantiation, resulting in
   substantially quicker initialization times for long trajectories. The
   ctime and filesize of the trajectory are stored with the offsets, and these
   are checked against the trajectory on load to ensure the offsets aren't
   stale. The offsets are automatically regenerated if they are stale or
   missing.

.. autoclass:: Timestep
   :members:
.. autoclass:: TrjReader
   :members:
.. autoclass:: TrjWriter
   :members:

"""

import os
import errno
import numpy
import sys
import cPickle
import warnings

import libxdrfile2
import statno
from MDAnalysis.coordinates import base
from MDAnalysis.coordinates.core import triclinic_box, triclinic_vectors
import MDAnalysis.core


# This is the XTC class. The TRR overrides with it's own.
class Timestep(base.Timestep):
    """Timestep for a Gromacs trajectory."""

    def __init__(self, arg, **kwargs):
        DIM = libxdrfile2.DIM  # compiled-in dimension (most likely 3)
        if numpy.dtype(type(arg)) == numpy.dtype(int):
            self.frame = 0
            self.numatoms = arg
            # C floats and C-order for arrays (see libxdrfile2.i)
            self._pos = numpy.zeros((self.numatoms, DIM), dtype=numpy.float32, order='C')
            self._unitcell = numpy.zeros((DIM, DIM), dtype=numpy.float32)
            # additional data for xtc
            self.status = libxdrfile2.exdrOK
            self.step = 0
            self.time = 0
            self.prec = 0
        elif isinstance(arg, Timestep):  # Copy constructor
            # This makes a deepcopy of the timestep
            self.frame = arg.frame
            self.numatoms = arg.numatoms
            self._unitcell = numpy.array(arg._unitcell)
            try:
                self._pos = numpy.array(arg._pos, dtype=numpy.float32)
            except ValueError as err:
                raise ValueError("Attempted to create a Timestep with invalid coordinate data: " + err.message)
            for attr in ('status', 'step', 'time', 'prec', 'lmbda'):
                if hasattr(arg, attr):
                    self.__setattr__(attr, arg.__getattribute__(attr))
        elif isinstance(arg, numpy.ndarray):
            if len(arg.shape) != 2:
                raise ValueError("numpy array can only have 2 dimensions")
            self._unitcell = numpy.zeros((DIM, DIM), dtype=numpy.float32)
            self.frame = 0
            if arg.shape[0] == DIM and arg.shape[1] != DIM:  # wrong order (need to exclude natoms == DIM!)
                raise ValueError("Coordinate array must have shape (natoms, %(DIM)d)" % vars())
            if arg.shape[1] == DIM:
                self.numatoms = arg.shape[0]
            else:
                raise ValueError("XDR timestep has not second dimension 3: shape=%r" % (arg.shape,))
            self._pos = arg.astype(numpy.float32).copy('C')  # C-order
            # additional data for xtc
            self.status = libxdrfile2.exdrOK
            self.step = 0
            self.time = 0
            self.prec = 0
        else:
            raise ValueError("Cannot create an empty Timestep")
        self._x = self._pos[:, 0]
        self._y = self._pos[:, 1]
        self._z = self._pos[:, 2]

    @property
    def dimensions(self):
        """unitcell dimensions (A, B, C, alpha, beta, gamma)

        - A, B, C are the lengths of the primitive cell vectors e1, e2, e3
        - alpha = angle(e1, e2)
        - beta = angle(e1, e3)
        - gamma = angle(e2, e3)
        """
        # Layout of unitcell is [X, Y, Z] with the primitive cell vectors
        x = self._unitcell[0]
        y = self._unitcell[1]
        z = self._unitcell[2]
        return triclinic_box(x, y, z)

    @dimensions.setter
    def dimensions(self, box):
        self._unitcell = triclinic_vectors(box)


class TrjWriter(base.Writer):
    """Writes to a Gromacs trajectory file

    (Base class)
    """
    #: units of time (ps) and length (nm) in Gromacs
    units = {'time': 'ps', 'length': 'nm'}
    #: override to define trajectory format of the reader (XTC or TRR)
    format = None

    def __init__(self, filename, numatoms, start=0, step=1, delta=None, precision=1000.0, remarks=None,
                 convert_units=None):
        """ Create a new TrjWriter

        :Arguments:
          *filename*
             name of output file
          *numatoms*
             number of atoms in trajectory file

        :Keywords:
          *start*
             starting timestep; only used when *delta* is set.
          *step*
             skip between subsequent timesteps; only used when *delta* is set.
          *delta*
             timestep to use. If set will override any time information contained in the
             passed :class:`Timestep` objects; otherwise that will be used. If in the
             latter case :attr:`~Timestep.time` is unavailable the TrjWriter will default
             to setting the trajectory time at 1 MDAnalysis unit (typically 1ps) per step.
          *precision*
              accuracy for lossy XTC format as a power of 10 (ignored
              for TRR) [1000.0]
          *convert_units*
             ``True``: units are converted to the MDAnalysis base format; ``None`` selects
             the value of :data:`MDAnalysis.core.flags` ['convert_lengths'].
             (see :ref:`flags-label`)

        .. versionchanged:: 0.8.0
           The TRR writer is now able to write TRRs without coordinates/velocities/forces,
           depending on the properties available in the :class:`Timestep` objects passed to
           :meth:`~TRRWriter.write`.
        """
        assert self.format in ('XTC', 'TRR')

        if numatoms == 0:
            raise ValueError("TrjWriter: no atoms in output trajectory")
        self.filename = filename
        # Convert filename to ascii because of SWIG bug.
        # See: http://sourceforge.net/p/swig/feature-requests/75
        # Only needed for Python < 3
        if sys.version_info[0] < 3:
            if isinstance(filename, unicode):
                self.filename = filename.encode("UTF-8")

        if convert_units is None:
            convert_units = MDAnalysis.core.flags['convert_lengths']
        self.convert_units = convert_units  # convert length and time to base units on the fly?
        self.numatoms = numatoms

        self.frames_written = 0
        self.start = start
        self.step = step
        self.delta = delta
        self.remarks = remarks
        self.precision = precision  # only for XTC
        self.xdrfile = libxdrfile2.xdrfile_open(self.filename, 'w')

        self.ts = None
        # To flag empty properties to be skipped when writing a TRR it suffices to pass an empty 2D array with shape(
        # natoms,0)
        if self.format == 'TRR':
            self._emptyarr = numpy.array([], dtype=numpy.float32).reshape(self.numatoms, 0)

    def write_next_timestep(self, ts=None):
        """ write a new timestep to the trj file

        *ts* is a :class:`Timestep` instance containing coordinates to
        be written to trajectory file
        """
        if self.xdrfile is None:
            raise IOError("Attempted to write to closed file %r", self.filename)
        if ts is None:
            if not hasattr(self, "ts"):
                raise IOError("TrjWriter: no coordinate data to write to trajectory file")
            else:
                ts = self.ts
        elif not ts.numatoms == self.numatoms:
            # Check to make sure Timestep has the correct number of atoms
            raise IOError("TrjWriter: Timestep does not have the correct number of atoms")

        status = self._write_next_timestep(ts)

        if status != libxdrfile2.exdrOK:
            raise IOError(errno.EIO, "Error writing %s file (status %d)" % (self.format, status), self.filename)
        self.frames_written += 1

    def _write_next_timestep(self, ts):
        """Generic writer for XTC and TRR with minimum intelligence; override if necessary."""

        # (1) data common to XTC and TRR

        # Time-writing logic: if the writer was created with a delta parameter,
        #  use delta*(start+step*frames_written)
        #  otherwise use the provided Timestep obj time attribute. Fallback to 1 if not present.
        if self.delta is None:
            try:
                time = ts.time
            except AttributeError:
                # Default to 1 time unit.
                time = self.start + self.step * self.frames_written
        else:
            time = (self.start + self.step * self.frames_written) * self.delta

        if self.convert_units:
            time = self.convert_time_to_native(time, inplace=False)

        if not hasattr(ts, 'step'):
            # bogus, should be actual MD step number, i.e. frame * delta/dt
            ts.step = ts.frame
        unitcell = self.convert_dimensions_to_unitcell(ts).astype(numpy.float32)  # must be float32 (!)

        # make a copy of the scaled positions so that the in-memory
        # timestep is not changed (would have lead to wrong results if
        # analysed *after* writing a time step to disk). The new
        # implementation could lead to memory problems and/or slow-down for
        # very big systems because we temporarily create a new array pos
        # for each frame written
        #
        # For TRR only go through the trouble if the frame actually has valid
        #  coords/vels/forces; otherwise they won't be written anyway (pointers
        #  set to an empty array that libxdrfile2.py knows it should set to NULL).
        #
        # (2) have to treat XTC and TRR somewhat differently
        if self.format == 'XTC':
            if self.convert_units:
                pos = self.convert_pos_to_native(ts._pos, inplace=False)
            else:
                pos = ts._pos
            status = libxdrfile2.write_xtc(self.xdrfile, int(ts.step), float(time), unitcell, pos, self.precision)
        #
        elif self.format == 'TRR':
            if not hasattr(ts, 'lmbda'):
                ts.lmbda = 1.0
            # Same assignment logic as for TRR Timestep creation (because we might be getting here frames from
            #  other formats' Timesteps that don't have 'has_' flags).
            has_x = ts.__dict__.get("has_x", True)
            has_v = ts.__dict__.get("has_v", hasattr(ts, "_velocities"))
            has_f = ts.__dict__.get("has_f", hasattr(ts, "_forces"))
            # COORDINATES
            if has_x:
                if self.convert_units:
                    pos = self.convert_pos_to_native(ts._pos, inplace=False)
                else:
                    pos = ts._pos
            else:
                pos = self._emptyarr
            #VELOCITIES
            if has_v:
                if self.convert_units:
                    velocities = self.convert_velocities_to_native(ts._velocities, inplace=False)
                else:
                    velocities = ts._velocities
            else:
                velocities = self._emptyarr
            # FORCES
            if has_f:
                if self.convert_units:
                    forces = self.convert_forces_to_native(ts._forces, inplace=False)
                else:
                    forces = ts._forces
            else:
                forces = self._emptyarr
            #
            status = libxdrfile2.write_trr(self.xdrfile, int(ts.step), float(time), float(ts.lmbda), unitcell,
                                           pos, velocities, forces)
        else:
            raise NotImplementedError("Gromacs trajectory format %s not known." % self.format)
        return status

    def close(self):
        status = libxdrfile2.exdrCLOSE
        if not self.xdrfile is None:
            status = libxdrfile2.xdrfile_close(self.xdrfile)
            self.xdrfile = None
        return status

    def convert_dimensions_to_unitcell(self, ts):
        """Read dimensions from timestep *ts* and return Gromacs box vectors"""
        return self.convert_pos_to_native(triclinic_vectors(ts.dimensions))


class TrjReader(base.Reader):
    """Generic base class for reading Gromacs trajectories inside MDAnalysis.

    Derive classes and set :attr:`TrjReader.format`, :attr:`TrjReader._read_trj` and :attr:`TrjReader._read_trj_atoms`.

    Example::
       reader = TrjReader("file.trj")
       for ts in reader:
          print ts

    """
    #: units of time (ps) and length (nm) in Gromacs
    units = {'time': 'ps', 'length': 'nm'}
    #: override to define trajectory format of the reader (XTC or TRR)
    format = None
    #: supply the appropriate Timestep class, e.g.
    #: :class:`MDAnalysis.coordinates.xdrfile.XTC.Timestep` for XTC
    _Timestep = Timestep
    #: writer class that matches this reader (override appropriately)
    _Writer = TrjWriter

    def __init__(self, filename, convert_units=None, sub=None, **kwargs):
        """
        :Arguments:
            *filename*
                the name of the trr file.

        :Keywords:
            *sub*
                an numpy integer array of what subset of trajectory atoms to load into
                the timestep. Intended to work similarly to the 'sub' argument to Gromacs_' trjconv.

                This is usefull when one has a Universe loaded with only an unsolvated protein, and
                wants to read a solvated trajectory.

                The length of this array must be <= to the actual number of atoms in the trajectory, and
                equal to number of atoms in the Universe.
            *refresh_offsets*
                if ``True``, do not retrieve stored offsets, but instead generate
                new ones; if ``False``, use retrieved offsets if available [``False``]

        .. versionchanged:: 0.9.0
           New keyword *refresh_offsets*

        """
        self.filename = filename
        # Convert filename to ascii because of SWIG bug.
        # See: http://sourceforge.net/p/swig/feature-requests/75
        # Only needed for Python < 3
        if sys.version_info[0] < 3:
            if isinstance(filename, unicode):
                self.filename = filename.encode("UTF-8")

        if convert_units is None:
            convert_units = MDAnalysis.core.flags['convert_lengths']
        self.convert_units = convert_units  # convert length and time to base units on the fly?
        self.xdrfile = None

        self.__numframes = None  # takes a long time, avoid accessing self.numframes
        self.skip_timestep = 1  # always 1 for xdr files
        self.__delta = None  # compute from time in first two frames!
        self.fixed = 0  # not relevant for Gromacs xtc/trr
        self.__offsets = None  # storage of offsets in the file
        self.skip = 1
        self.periodic = False

        # actual number of atoms in the trr file
        # first time file is opened, exception should be thrown if bad file
        self.__trr_numatoms = self._read_trj_natoms(self.filename)

        # logic for handling sub sections of trr:
        # this class has some tmp buffers into which the libxdrfile2 functions read the
        # entire frame, then this class copies the relevant sub section into the timestep.
        # the sub section logic is contained entierly within this class.

        # check to make sure sub is valid (if given)
        if sub is not None:
            # only valid type
            if not isinstance(sub, numpy.ndarray) or len(sub.shape) != 1 or sub.dtype.kind != 'i':
                raise TypeError("sub MUST be a single dimensional numpy array of integers")
            if len(sub) > self.__trr_numatoms:
                raise ValueError("sub MUST be less than or equal to the number of actual trr atoms, {} in this case".
                                 format(self.__trr_numatoms))
            if numpy.max(sub) >= self.__trr_numatoms or numpy.min(sub) < 0:
                raise IndexError("sub contains out-of-range elements for the given trajectory")
            # sub appears to be valid
            self.__sub = sub

            # make tmp buffers
            # C floats and C-order for arrays (see libxdrfile2.i)
            DIM = libxdrfile2.DIM  # compiled-in dimension (most likely 3)
            # both xdr and trr have positions
            self.__pos_buf = numpy.zeros((self.__trr_numatoms, DIM), dtype=numpy.float32, order='C')
            if self.format == "TRR":
                # only trr has velocities and forces
                self.__velocities_buf = numpy.zeros((self.__trr_numatoms, DIM), dtype=numpy.float32, order='C')
                self.__forces_buf = numpy.zeros((self.__trr_numatoms, DIM), dtype=numpy.float32, order='C')
            else:
                # xdr does not have vel or forces.
                self.__velocities_buf = None
                self.__forces_buf = None
        else:
            self.__sub = None
            self.__pos_buf = None
            self.__velocities_buf = None
            self.__forces_buf = None

        # make the timestep, this is ALWAYS the used the public number of atoms
        # (same as the calling Universe)
        # at this time, __trr_numatoms and __sub are set, so self.numatoms has all it needs
        # to determine number of atoms.
        self.ts = self._Timestep(self.numatoms)

        # Read in the first timestep
        self._read_next_timestep()

        # try retrieving stored offsets
        if not kwargs.pop('refresh_offsets', False):
            self._retrieve_offsets()

    @property
    def numatoms(self):
        """The number of publically available atoms that this reader will store in the timestep.

        If 'sub' was not given in the ctor, then this value will just be the actual
        number of atoms in the underlying trajectory file. If however 'sub' was given, then
        this value is the number specified by the 'sub' sub-selection.

        If for any reason the trajectory cannot be read then a negative value is returned.
        """
        return len(self.__sub) if self.__sub is not None else self.__trr_numatoms

    @property
    def numframes(self):
        """Read the number of frames from the trajectory.

        The result is cached. If for any reason the trajectory cannot
        be read then 0 is returned.

        This takes a long time because the frames are counted by
        iterating through the whole trajectory. If the trajectory
        was previously loaded and saved offsets exist, then
        loading will be significantly faster.

        .. SeeAlso:: :meth:`TrjReader.load_offsets` and :meth:`TrjReader.save_offsets`
        """
        if not self.__numframes is None:  # return cached value
            return self.__numframes
        try:
            self._read_trj_numframes(self.filename)
        except IOError:
            self.__numframes = 0
            return 0
        else:
            return self.__numframes

    @property
    def offsets(self):
        if self.__offsets is not None:
            return self.__offsets
        try:
            self._read_trj_numframes(self.filename)
        except IOError:
            self.__offsets = []
            return 0
        else:
            return self.__offsets

    @property
    def delta(self):
        """Time step length in ps.

        The result is computed from the trajectory and cached. If for
        any reason the trajectory cannot be read then 0 is returned.
        """
        # no need for conversion: it's alread in our base unit ps
        if not self.__delta is None:  # return cached value
            return self.__delta
        try:
            t0 = self.ts.time
            self.next()
            t1 = self.ts.time
            self.__delta = t1 - t0
        except IOError:
            return 0
        finally:
            self.rewind()
        return self.__delta

    def _read_trj_natoms(self, filename):
        """Generic number-of-atoms extractor with minimum intelligence. Override if necessary."""
        if self.format == 'XTC':
            numatoms = libxdrfile2.read_xtc_natoms(filename)
        elif self.format == 'TRR':
            numatoms = libxdrfile2.read_trr_natoms(filename)
        else:
            raise NotImplementedError("Gromacs trajectory format %s not known." % self.format)
        return numatoms

    def _read_trj_numframes(self, filename):
        """Number-of-frames extractor/indexer for XTC and TRR."""
        if self.format == 'XTC':
            self.__numframes, self.__offsets = libxdrfile2.read_xtc_numframes(filename)
            self._store_offsets()
        elif self.format == 'TRR':
            self.__numframes, self.__offsets = libxdrfile2.read_trr_numframes(filename)
            self._store_offsets()
        else:
            raise NotImplementedError("Gromacs trajectory format %s not known." % self.format)
        return

    def _offset_filename(self):
        head, tail = os.path.split(self.filename)
        return os.path.join(head, '.{}_offsets.pkl'.format(tail))

    def _store_offsets(self):
        """Stores offsets for trajectory as a hidden file in the same directory
            as the trajectory itself.

        .. versionadded: 0.9.0
        """
        # try to store offsets; if fails (due perhaps to permissions), then
        # don't bother
        try:
            self.save_offsets(self._offset_filename())
        except IOError:
            warnings.warn("Offsets could not be stored; they will rebuilt when needed next.")

    def _retrieve_offsets(self):
        """Attempts to retrieve previously autosaved offsets for trajectory.

        .. versionadded: 0.9.0
        """
        try:
            self.load_offsets(self._offset_filename(), check=True)
        except IOError:
            warnings.warn("Offsets could not be retrieved; they will be rebuilt instead.")

    def save_offsets(self, filename):
        """Saves current trajectory offsets into *filename*, as a pickled object.

        Along with the offsets themselves, the ctime and file size of the
        trajectory file are also saved.  These are used upon load as a check to
        ensure the offsets still match the trajectory they are being applied
        to.

        The offset file is a pickled dictionary with keys/values::
          *ctime*
             the ctime of the trajectory file
          *size*
             the size of the trajectory file
          *offsets*
             a numpy array of the offsets themselves

        :Arguments:
          *filename*
              filename in which to save the frame offsets

        .. versionadded: 0.8.0
        .. versionchanged: 0.9.0
           Format of the offsets file has changed. It is no longer a pickled numpy
           array, but now a pickled dictionary. See details above. Old offset files
           can no longer be loaded.
            
        """
        if self.__offsets is None:
            self._read_trj_numframes(self.filename)

        output = {'ctime': os.path.getctime(self.filename),
                  'size': os.path.getsize(self.filename),
                  'offsets': self.__offsets}

        with open(filename, 'wb') as f:
            cPickle.dump(output, f)

    def load_offsets(self, filename, check=False):
        """Loads current trajectory offsets from pickled *filename*. 

        Checks if ctime and size of trajectory file matches that stored in
        pickled *filename*.  If either one does not match (and *check* == ``True``)
        then the offsets are not loaded.  This is intended to conservatively
        avoid loading out-of-date offsets.

        The offset file is expected to be a pickled dictionary with keys/values::
          *ctime*
             the ctime of the trajectory file
          *size*
             the size of the trajectory file
          *offsets*
             a numpy array of the offsets themselves

        :Arguments:
          *filename*
              filename of pickle file saved with :meth:`~TrjReader.save_offsets` 
              with the frame offsets for the loaded trajectory

        :Keywords:
          *check*
              if False, ignore ctime and size check of trajectory file
        
        :Raises: :exc:`IOError` if the file cannot be read (see :func:`open`).

        .. versionadded: 0.8.0
        .. versionchanged: 0.9.0
           Format of the offsets file has changed. It is no longer a pickled numpy
           array, but now a pickled dictionary. See details above. Old offset files
           can no longer be loaded.

        """
        if not os.path.isfile(filename):  # Return silently if the offset file is not present
            return

        with open(filename, 'rb') as f:
            offsets = cPickle.load(f)

        if check:
            conditions = False
            try:
                ## ensure all conditions are met
                # ctime of file must match that stored
                key = 'ctime'
                conditions = (os.path.getctime(self.filename) == offsets[key])

                # file size must also match
                key = 'size'
                conditions = (os.path.getsize(self.filename) == offsets[key]) and conditions
            except KeyError:
                warnings.warn("Offsets in file '{}' not suitable; missing {}.".format(filename, key))
                return

            # if conditions not met, abort immediately
            if not conditions:
                warnings.warn("Aborted loading offsets from file; ctime or size did not match.")
                return

        # try to load offsets
        try:
            self.__offsets = offsets['offsets']
        except KeyError:
            warnings.warn("Missing key 'offsets' in file '{}'; aborting load of offsets.".format(filename))
            return
        self.__numframes = len(self.__offsets)

        # finally, check that loaded offsets appear to work by trying
        # to load last frame; otherwise, dump them so they get regenerated
        # on next call to ``self.numframes``

        #store current frame
        frame = self.frame
        try:
            self.__getitem__(-1)
            # ensure we return to the frame we started with
            self.__getitem__(frame - 1)
        except (IndexError, IOError):
            warnings.warn("Could not access last frame with loaded offsets; will rebuild offsets instead.")
            self.__offsets = None
            self.__numframes = None

    def open_trajectory(self):
        """Open xdr trajectory file.

        :Returns: pointer to XDRFILE (and sets self.xdrfile)
        :Raises:  :exc:`IOError` with code EALREADY if file was already opened or
                  ENOENT if the file cannot be found
        """
        if not self.xdrfile is None:
            raise IOError(errno.EALREADY, 'XDR file already opened', self.filename)
        if not os.path.exists(self.filename):
            # must check; otherwise might segmentation fault
            raise IOError(errno.ENOENT, 'XDR file not found', self.filename)
        self.xdrfile = libxdrfile2.xdrfile_open(self.filename, 'r')
        # reset ts
        ts = self.ts
        ts.status = libxdrfile2.exdrOK
        ts.frame = 0
        ts.step = 0
        ts.time = 0
        # additional data for XTC
        ts.prec = 0
        # additional data for TRR
        ts.lmbda = 0
        return self.xdrfile

    def close(self):
        """Close xdr trajectory file if it was open."""
        if self.xdrfile is None:
            return
        libxdrfile2.xdrfile_close(self.xdrfile)
        self.xdrfile = None  # guard against  crashing with a double-free pointer

    def Writer(self, filename, **kwargs):
        """Returns a Gromacs TrjWriter for *filename* with the same parameters as this trajectory.

        All values can be changed through keyword arguments.

        :Arguments:
          *filename*
              filename of the output trajectory

        :Keywords:
          *numatoms*
              number of atoms
          *delta*
              Time interval between frames.
          *precision*
              accuracy for lossy XTC format as a power of 10 (ignored
              for TRR) [1000.0]

        :Returns: appropriate :class:`TrjWriter`
        """
        numatoms = kwargs.pop('numatoms', self.numatoms)
        kwargs['step'] = self.skip_timestep
        kwargs.setdefault('delta', self.delta)
        try:
            kwargs['start'] = self[0].time / kwargs['delta']
        except (AttributeError, ZeroDivisionError):
            kwargs['start'] = 0
        try:
            kwargs.setdefault('precision', self.precision)
        except AttributeError:
            pass  # not needed for TRR
        return self._Writer(filename, numatoms, **kwargs)

    def __iter__(self):
        self.ts.frame = 0  # start at 0 so that the first frame becomes 1
        self._reopen()
        while True:
            try:
                ts = self._read_next_timestep()
            except IOError as err:
                if err.errno == errno.EIO:
                    break
                else:
                    self.close()
                    raise
            except:
                self.close()
                raise
            else:
                yield ts

    def _read_next_timestep(self, ts=None):
        """Generic ts reader with minimum intelligence. Override if necessary."""
        if ts is None:
            ts = self.ts
        elif self.format == 'TRR' and not hasattr(ts,
                                                  '_tpos'):  # If a foreign Timestep is passed as the receptacle of
                                                  # the data
            ts = TRR.Timestep(ts)  # we must make sure the access-checking stuff gets set up.

        if self.xdrfile is None:
            self.open_trajectory()

        if self.format == 'XTC':
            if self.__sub is None:
                ts.status, ts.step, ts.time, ts.prec = libxdrfile2.read_xtc(self.xdrfile, ts._unitcell, ts._pos)
            else:
                ts.status, ts.step, ts.time, ts.prec = libxdrfile2.read_xtc(self.xdrfile, ts._unitcell, self.__pos_buf)
                ts._pos[:] = self.__pos_buf[self.__sub]
        elif self.format == 'TRR':
            if self.__sub is None:
                ts.status, ts.step, ts.time, ts.lmbda,\
                    ts.has_x, ts.has_v, ts.has_f = libxdrfile2.read_trr(self.xdrfile,
                                                                        ts._unitcell,
                                                                        ts._tpos,
                                                                        ts._tvelocities,
                                                                        ts._tforces)
            else:
                ts.status, ts.step, ts.time, ts.lmbda,\
                    ts.has_x, ts.has_v, ts.has_f = libxdrfile2.read_trr(self.xdrfile,
                                                                        ts._unitcell,
                                                                        self.__pos_buf,
                                                                        self.__velocities_buf,
                                                                        self.__forces_buf)
                ts._tpos[:] = self.__pos_buf[self.__sub]
                ts._tvelocities[:] = self.__velocities_buf[self.__sub]
                ts._tforces[:] = self.__forces_buf[self.__sub]
        else:
            raise NotImplementedError("Gromacs trajectory format %s not known." % self.format)

        if (ts.status == libxdrfile2.exdrENDOFFILE) or \
                (ts.status == libxdrfile2.exdrINT and self.format == 'TRR'):
            # seems that trr files can get a exdrINT when reaching EOF (??)
            raise IOError(errno.EIO, "End of file reached for %s file" % self.format,
                          self.filename)
        elif not ts.status == libxdrfile2.exdrOK:
            raise IOError(errno.EBADF, "Problem with %s file, status %s" %
                                       (self.format, statno.errorcode[ts.status]), self.filename)
        if self.convert_units:
            # TRRs have the annoying possibility of frames without coordinates/velocities/forces...
            if self.format == 'XTC' or ts.has_x:
                self.convert_pos_from_native(ts._pos)  # in-place !
            self.convert_pos_from_native(ts._unitcell)  # in-place ! (note: xtc/trr contain unit vecs!)
            ts.time = self.convert_time_from_native(ts.time)  # in-place does not work with scalars
            if self.format == 'TRR':
                if ts.has_v:
                    self.convert_velocities_from_native(ts._velocities)  # in-place
                if ts.has_f:
                    self.convert_forces_from_native(ts._forces)  # in-place
        ts.frame += 1
        return ts

    def rewind(self):
        """Position at beginning of trajectory"""
        self._reopen()
        self.next()  # read first frame

    def _reopen(self):
        self.close()
        self.open_trajectory()

    def _goto_frame(self, frame):
        """ Fast, index-based, random frame access

        :TODO: Should we treat frame as 0-based or 1-based??  Right
               now: 0-based (which is inconsistent) but analogous to
               DCDReader
        """
        if self.__offsets is None:
            self._read_trj_numframes(self.filename)
        self._seek(self.__offsets[frame])
        self._read_next_timestep()
        self.ts.frame = frame + 1
        return self.ts

    def __getitem__(self, frame):
        if (numpy.dtype(type(frame)) != numpy.dtype(int)) and (type(frame) != slice):
            raise TypeError
        if numpy.dtype(type(frame)) == numpy.dtype(int):
            if frame < 0:
                # Interpret similar to a sequence
                frame += len(self)
            if (frame < 0) or (frame >= len(self)):
                raise IndexError("0 <= frame < len(traj) is outside of trajectory boundaries")
            return self._goto_frame(frame)
        elif type(frame) == slice:  # if frame is a slice object
            if not (((type(frame.start) == int) or (frame.start is None)) and
               ((type(frame.stop) == int) or (frame.stop is None)) and
               ((type(frame.step) == int) or (frame.step is None))):
                raise TypeError("Slice indices are not integers")

            def _iter(start=frame.start, stop=frame.stop, step=frame.step):
                # check_slices implicitly calls len(self), thus setting numframes
                start, stop, step = self._check_slice_indices(start, stop, step)
                if step == 1:
                    rng = xrange(start + 1, stop, step)
                    yield self._goto_frame(start)
                    for framenr in rng:
                        yield self._read_next_timestep()
                else:
                    rng = xrange(start, stop, step)
                    for framenr in rng:
                        yield self._goto_frame(framenr)

            return _iter()

    def _check_slice_indices(self, start, stop, step):
        if step is None:
            step = 1
        if (start < 0) and start is not None:
            start += len(self)
        if (stop < 0) and stop is not None:
            stop += len(self)
        if step > 0:
            if start is None:
                start = 0
            if stop is None:
                stop = len(self)
        else:  # Fixes reverse iteration with default slice. (eg. trajectory[::-1])
            if start is None:
                start = len(self) - 1
            if stop is None:
                stop = -1
        if stop > len(self):
            stop = len(self)

        if step > 0 and stop <= start:
            raise IndexError("Stop frame is lower than start frame")
        if step == 0:
            raise IndexError("Iteration step 0.")
        if start >= len(self):
            raise IndexError("Frame start outside of the range of the trajectory.")
        return start, stop, step

    def timeseries(self, asel, start=0, stop=-1, skip=1, format='afc'):
        raise NotImplementedError("timeseries not available for Gromacs trajectories")

    def correl(self, timeseries, start=0, stop=-1, skip=1):
        raise NotImplementedError("correl not available for Gromacs trajectories")

    def __del__(self):
        self.close()

    def _seek(self, pos, rel=False):
        """Traj seeker"""
        if self.format == 'XTC' or self.format == 'TRR':
            if rel:
                status = libxdrfile2.xdr_seek(self.xdrfile, long(pos), libxdrfile2.SEEK_CUR)
            else:
                status = libxdrfile2.xdr_seek(self.xdrfile, long(pos), libxdrfile2.SEEK_SET)
            if status != libxdrfile2.exdrOK:
                raise IOError(errno.EIO,
                              "Problem seeking to offset %d (relative to current = %s) on file %s, status %s.\n"
                              "Perhaps you are trying to read a file >2GB and your system does not have large file"
                              "support?" % (pos, rel, self.filename, status))
        else:
            raise NotImplementedError("Gromacs trajectory format %s not known." % self.format)

    def _tell(self):
        """Traj pos getter"""
        if self.format == 'XTC' or self.format == 'TRR':
            offset = libxdrfile2.xdr_tell(self.xdrfile)
        else:
            raise NotImplementedError("Gromacs trajectory format %s not known." % self.format)
        return offset
