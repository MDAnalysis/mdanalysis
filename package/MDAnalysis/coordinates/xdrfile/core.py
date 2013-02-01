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
Common high-level functionality for accessing Gromacs trajectories
==================================================================

The :mod:`MDAnalysis.coordinates.xdrfile.core` module contains generic
classes to access Gromacs_ XDR-encoded trajectory formats such as TRR
and XTC.

A generic Gromacs_ trajectory is simply called "trj" within this
module.

.. SeeAlso:: :mod:`MDAnalysis.coordinates.base` for the generic MDAnalysis base
             classes and :mod:`MDAnalysis.coordinates.xdrfile.libxdrfile` for
             the low-level bindings to the XDR trajectories.

.. _Gromacs: http://www.gromacs.org

Generic xdr trj classes
-----------------------

The generic classes are subclassed to generate the specific classes
for the XTC and TRR format.

.. autoclass:: Timestep
   :members:
.. autoclass:: TrjReader
   :members:
.. autoclass:: TrjWriter
   :members:

"""

import os.path
import errno
import numpy

import libxdrfile, statno
from MDAnalysis.coordinates import base
from MDAnalysis.coordinates.core import triclinic_box, triclinic_vectors

import MDAnalysis.core

class Timestep(base.Timestep):
    """Timestep for a Gromacs trajectory."""
    def __init__(self, arg):
        DIM = libxdrfile.DIM    # compiled-in dimension (most likely 3)
        if numpy.dtype(type(arg)) == numpy.dtype(int):
            self.frame = 0
            self.numatoms = arg
            # C floats and C-order for arrays (see libxdrfile.i)
            self._pos = numpy.zeros((self.numatoms, DIM), dtype=numpy.float32, order='C')
            self._unitcell = numpy.zeros((DIM,DIM), dtype=numpy.float32)
            # additional data for xtc
            self.status = libxdrfile.exdrOK
            self.step = 0
            self.time = 0
            self.prec = 0
        elif isinstance(arg, Timestep): # Copy constructor
            # This makes a deepcopy of the timestep
            self.frame = arg.frame
            self.numatoms = arg.numatoms
            self._unitcell = numpy.array(arg._unitcell)
            self._pos = numpy.array(arg._pos)
            for attr in ('status', 'step', 'time', 'prec', 'lmbda'):
                if hasattr(arg, attr):
                    self.__setattr__(attr, arg.__getattribute__(attr))
        elif isinstance(arg, numpy.ndarray):
            if len(arg.shape) != 2:
                raise ValueError("numpy array can only have 2 dimensions")
            self._unitcell = numpy.zeros((DIM,DIM), dtype=numpy.float32)
            self.frame = 0
            if arg.shape[0] == DIM and arg.shape[1] != DIM:    ## wrong order (need to exclude natoms == DIM!)
                raise ValueError("Coordinate array must have shape (natoms, %(DIM)d)" % vars())
            if arg.shape[1] == DIM:
                self.numatoms = arg.shape[0]
            else:
                raise ValueError("XDR timestep has not second dimension 3: shape=%r" % (arg.shape,))
            self._pos = arg.astype(numpy.float32).copy('C')  ## C-order
            # additional data for xtc
            self.status = libxdrfile.exdrOK
            self.step = 0
            self.time = 0
            self.prec = 0
        else:
            raise ValueError("Cannot create an empty Timestep")
        self._x = self._pos[:,0]
        self._y = self._pos[:,1]
        self._z = self._pos[:,2]

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
        return triclinic_box(x,y,z)

class TrjWriter(base.Writer):
    """Writes to a Gromacs trajectory file

    (Base class)
    """
    #: units of time (ps) and length (nm) in Gromacs
    units = {'time': 'ps', 'length': 'nm'}
    #: override to define trajectory format of the reader (XTC or TRR)
    format = None

    def __init__(self, filename, numatoms, start=0, step=1, delta=1.0, precision=1000.0, remarks=None,
                  convert_units=None):
        ''' Create a new TrjWriter

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
          *precision*
             accuracy for lossy XTC format [1000]
          *convert_units*
             ``True``: units are converted to the MDAnalysis base format; ``None`` selects
             the value of :data:`MDAnalysis.core.flags` ['convert_gromacs_lengths'].
             (see :ref:`flags-label`)
        '''
        assert self.format in ('XTC', 'TRR')

        if numatoms == 0:
            raise ValueError("TrjWriter: no atoms in output trajectory")
        self.filename = filename
        if convert_units is None:
            convert_units = MDAnalysis.core.flags['convert_gromacs_lengths']
        self.convert_units = convert_units    # convert length and time to base units on the fly?
        self.numatoms = numatoms

        self.frames_written = 0
        self.start = start
        self.step = step
        self.delta = delta
        self.remarks = remarks
        self.precision = precision  # only for XTC
        self.xdrfile = libxdrfile.xdrfile_open(filename, 'w')

        self.ts = None

    def write_next_timestep(self, ts=None):
        ''' write a new timestep to the trj file

        *ts* is a :class:`Timestep` instance containing coordinates to
        be written to trajectory file
        '''
        if self.xdrfile is None:
            raise IOError("Attempted to write to closed file %r", self.filename)
        if ts is None:
            if not hasattr(self, "ts"):
                raise IOError("TrjWriter: no coordinate data to write to trajectory file")
            else:
                ts=self.ts
        elif not ts.numatoms == self.numatoms:
            # Check to make sure Timestep has the correct number of atoms
            raise IOError("TrjWriter: Timestep does not have the correct number of atoms")

        status = self._write_next_timestep(ts)

        if status != libxdrfile.exdrOK:
            raise IOError(errno.EIO, "Error writing %s file (status %d)" % (self.format, status), self.filename)
        self.frames_written += 1

    def _write_next_timestep(self, ts):
        """Generic writer for XTC and TRR with minimum intelligence; override if necessary."""

        # (1) data common to XTC and TRR
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
        else:
            pos = ts._pos
            try:
                time = ts.time
            except AttributeError:
                time = ts.frame * self.delta
        if not hasattr(ts, 'step'):
            # bogus, should be actual MD step number, i.e. frame * delta/dt
            ts.step = ts.frame
        unitcell = self.convert_dimensions_to_unitcell(ts).astype(numpy.float32)  # must be float32 (!)

        # (2) have to treat XTC and TRR somewhat differently
        if self.format == 'XTC':
            status = libxdrfile.write_xtc(self.xdrfile, int(ts.step), float(time), unitcell, pos, self.precision)
        elif self.format == 'TRR':
            if not hasattr(ts, 'lmbda'):
                ts.lmbda = 1.0
            if hasattr(ts, '_velocities'):
                if self.convert_units:
                    velocities = self.convert_velocities_to_native(ts._velocities, inplace=False)
                else:
                    velocities = ts._velocities
            else:
                # 0-velocities (no need to convert those); add it to ts as a sideeffect
                # so that we don't have to reallocate next time
                velocities = ts._velocities = numpy.zeros((3,ts.numatoms), dtype=numpy.float32)
            if hasattr(ts, '_forces'):
                if self.convert_units:
                    forces = self.convert_forces_to_native(ts._forces, inplace=False)
                else:
                    forces = ts._velocities                
            else:
                # 0-forces (no need to convert those); add it to ts as a sideeffect
                # so that we don't have to reallocate next time
                forces = ts._forces = numpy.zeros((3,ts.numatoms), dtype=numpy.float32)
            status = libxdrfile.write_trr(self.xdrfile, int(ts.step), float(time), float(ts.lmbda), unitcell,
                                          pos, velocities, forces)
        else:
            raise NotImplementedError("Gromacs trajectory format %s not known." % self.format)
        return status

    def close(self):
        status = libxdrfile.exdrCLOSE
        if not self.xdrfile is None:
            status = libxdrfile.xdrfile_close(self.xdrfile)
            self.xdrfile = None
        return status

    def convert_dimensions_to_unitcell(self, ts):
        """Read dimensions from timestep *ts* and return Gromacs box vectors"""
        return self.convert_pos_to_native(triclinic_vectors(ts.dimensions))

class TrjReader(base.Reader):
    """Generic base class for reading Gromacs trajectories inside MDAnalysis.

    Derive classes and set :attr:`TrjReader.format`,
    :attr:`TrjReader._read_trj` and :attr:`TrjReader._read_trj_atoms`.

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
            *sub*
            
                an numpy integer array of what subset of trajectory atoms to load into
                the timestep. Intended to work similarly to the 'sub' argument to gromac's trjconv.

                This is usefull when one has a Universe loaded with only an unsolvated protein, and
                wants to read a solvated trajectory. 
                
                The length of this array must be <= to the actual number of atoms in the trajectory, and
                equal to number of atoms in the Universe.
            
        """
        self.filename = filename
        if convert_units is None:
            convert_units = MDAnalysis.core.flags['convert_gromacs_lengths']
        self.convert_units = convert_units  # convert length and time to base units on the fly?
        self.xdrfile = None

        self.__numframes = None  # takes a long time, avoid accessing self.numframes
        self.skip_timestep = 1   # always 1 for xdr files
        self.__delta = None      # compute from time in first two frames!
        self.fixed = 0           # not relevant for Gromacs xtc/trr
        self.skip = 1
        self.periodic = False
        
        # actual number of atoms in the trr file
        # first time file is opened, exception should be thrown if bad file
        self.__trr_numatoms = self._read_trj_natoms(self.filename)
        
        # logic for handling sub sections of trr: 
        # this class has some tmp buffers into which the libxdrfile functions read the 
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
                raise IndexError("sub contains out of range elements for the given trajectory")
            # sub appears to be valid
            self.__sub = sub
            
            # make tmp buffers
            # C floats and C-order for arrays (see libxdrfile.i)
            DIM = libxdrfile.DIM    # compiled-in dimension (most likely 3)
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

        This  takes a  long time  because  the frames  are counted  by
        iterating through the whole trajectory.
        """
        if not self.__numframes is None:   # return cached value
            return self.__numframes
        try:
            self.__numframes = self._read_trj_numframes(self.filename)
        except IOError:
            return 0
        else:
            return self.__numframes

    @property
    def delta(self):
        """Time step length in ps.

        The result is computed from the trajectory and cached. If for
        any reason the trajectory cannot be read then 0 is returned.
        """
        # no need for conversion: it's alread in our base unit ps
        if not self.__delta is None:   # return cached value
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
            numatoms = libxdrfile.read_xtc_natoms(self.filename)
        elif self.format == 'TRR':
            numatoms = libxdrfile.read_trr_natoms(self.filename)
        else:
            raise NotImplementedError("Gromacs trajectory format %s not known." % self.format)
        return numatoms

    def _read_trj_numframes(self, filename):
        """Generic numer-of-frames extractor with minimum intelligence. Override if necessary."""
        if self.format == 'XTC':
            numframes = libxdrfile.read_xtc_numframes(self.filename)
        elif self.format == 'TRR':
            numframes = libxdrfile.read_trr_numframes(self.filename)
        else:
            raise NotImplementedError("Gromacs trajectory format %s not known." % self.format)
        return numframes

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
        self.xdrfile = libxdrfile.xdrfile_open(self.filename, 'r')
        # reset ts
        ts = self.ts
        ts.status = libxdrfile.exdrOK
        ts.frame = 0
        ts.step = 0
        ts.time = 0
        # additional data for xtc
        ts.prec = 0
        # additional data for TRR
        ts.lmbda = 0
        return self.xdrfile

    def close(self):
        """Close xdr trajectory file if it was open."""
        if self.xdrfile is None:
            return
        libxdrfile.xdrfile_close(self.xdrfile)
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
              Time interval between frames in ps.
          *precision*
              accuracy for lossy XTC format as a power of 10 (ignored
              for TRR) [1000.0]

        :Returns: appropriate :class:`TrjWriter`
        """
        numatoms = kwargs.pop('numatoms', self.numatoms)
        kwargs.pop('start', None)            # ignored by TrjWriter
        kwargs['step'] = self.skip_timestep  # ignored/fixed to 1
        kwargs.setdefault('delta', self.delta)
        try:
            kwargs.setdefault('precision', self.precision)
        except AttributeError:
            pass                             # not needed for TRR
        return self._Writer(filename, numatoms, **kwargs)

    def __iter__(self):
        self.ts.frame = 0  # start at 0 so that the first frame becomes 1
        self._reopen()
        while True:
            try:
                ts = self._read_next_timestep()
            except IOError, err:
                if err.errno == errno.ENODATA:
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
        if self.xdrfile is None:
            self.open_trajectory()

        if self.format == 'XTC':
            if self.__sub is None:            
                ts.status, ts.step, ts.time, ts.prec = libxdrfile.read_xtc(self.xdrfile, ts._unitcell, ts._pos)
            else:
                ts.status, ts.step, ts.time, ts.prec = libxdrfile.read_xtc(self.xdrfile, ts._unitcell, self.__pos_buf)
                ts._pos[:] = self.__pos_buf[self.__sub]
        elif self.format == 'TRR':
            if self.__sub is None:
                ts.status, ts.step, ts.time, ts.lmbda = libxdrfile.read_trr(self.xdrfile, ts._unitcell, ts._pos,
                                                                            ts._velocities, ts._forces)
            else:
                ts.status, ts.step, ts.time, ts.lmbda = libxdrfile.read_trr(self.xdrfile, ts._unitcell, self.__pos_buf,
                                                                            self.__velocities_buf, self.__forces_buf)
                ts._pos[:] = self.__pos_buf[self.__sub]
                ts._velocities[:] = self.__velocities_buf[self.__sub]
                ts._forces[:] = self.__forces_buf[self.__sub]
        else:
            raise NotImplementedError("Gromacs trajectory format %s not known." % self.format)
        
        if (ts.status == libxdrfile.exdrENDOFFILE) or \
                (ts.status == libxdrfile.exdrINT and self.format == 'TRR'):
            # seems that trr files can get a exdrINT when reaching EOF (??)
            raise IOError(errno.ENODATA, "End of file reached for %s file" % self.format,
                          self.filename)
        elif not ts.status == libxdrfile.exdrOK:
            raise IOError(errno.EFAULT, "Problem with %s file, status %s" %
                          (self.format, statno.errorcode[ts.status]), self.filename)
        if self.convert_units:
            self.convert_pos_from_native(ts._pos)             # in-place !
            self.convert_pos_from_native(ts._unitcell)        # in-place ! (note: xtc/trr contain unit vecs!)
            ts.time = self.convert_time_from_native(ts.time)  # in-place does not work with scalars
            if self.format == 'TRR':
                self.convert_velocities_from_native(ts._velocities) # in-place
                self.convert_forces_from_native(ts._forces)     # in-place
        ts.frame += 1
        return ts

    def rewind(self):
        """Position at beginning of trajectory"""
        self._reopen()
        self.next()   # read first frame

    def _reopen(self):
        self.close()
        self.open_trajectory()

    def _forward_to_frame(self, frameindex):
        """Slow implementation: must read sequentially.

        .. Note::

           *frameindex* starts from 0; i.e. *frameindex* = frame - 1.

        :TODO: Should we treat frame as 0-based or 1-based??  Right
               now: 0-based (which is inconsistent) but analogous to
               DCDReader
        """
        self.rewind()
        if frameindex == 0:
            return
        for ts in self:
            if ts.frame-1 >= frameindex:
                break

    def __getitem__(self, frame):
        if (numpy.dtype(type(frame)) != numpy.dtype(int)) and (type(frame) != slice):
            raise TypeError
        if (numpy.dtype(type(frame)) == numpy.dtype(int)):
            if (frame < 0):
                # Interpret similar to a sequence
                frame = len(self) + frame
            if (frame < 0) or (frame >= len(self)):
                raise IndexError("0 <= frame < len(traj) is outside of trajectory boundaries")
            self._forward_to_frame(frame)
            return self.ts
        elif type(frame) == slice: # if frame is a slice object
            if not (((type(frame.start) == int) or (frame.start == None)) and
                    ((type(frame.stop) == int) or (frame.stop == None)) and
                    ((type(frame.step) == int) or (frame.step == None))):
                raise TypeError("Slice indices are not integers")
            def iterDCD(start=frame.start, stop=frame.stop, step=frame.step):
                start, stop, step = self._check_slice_indices(start, stop, step)
                if step < 0:
                    raise NotImplementedError("XTC/TRR do not support reverse iterating for performance reasons")
                for ts in self:
                    # simple sequential plodding through trajectory --- SLOW!
                    frameindex = ts.frame - 1
                    if frameindex < start:  continue
                    if frameindex >= stop:  break
                    if (frameindex - start) % step != 0: continue
                    yield self.ts
            return iterDCD()

    def timeseries(self, asel, start=0, stop=-1, skip=1, format='afc'):
        raise NotImplementedError("timeseries not available for Gromacs trajectories")
    def correl(self, timeseries, start=0, stop=-1, skip=1):
        raise NotImplementedError("correl not available for Gromacs trajectories")

    def __del__(self):
        self.close()
