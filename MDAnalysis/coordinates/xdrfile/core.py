"""Reading of Gromacs trajectories."""

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
            if len(arg.shape) != 2: raise Exception("numpy array can only have 2 dimensions")
            self._unitcell = numpy.zeros((DIM,DIM), dtype=numpy.float32)
            self.frame = 0
            if arg.shape[0] == DIM:    ## wrong order
                self.numatoms = arg.shape[-1]
            else: 
                self.numatoms = arg.shape[0]
            self._pos = arg.copy('C')  ## C-order ! (?) -- does this work or do I have to transpose?
            # additional data for xtc
            self.status = libxdrfile.exdrOK
            self.step = 0
            self.time = 0
            self.prec = 0
        else: 
            raise Exception("Cannot create an empty Timestep")
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

    def __init__(self, filename, convert_units=None, **kwargs):
        self.filename = filename
        if convert_units is None:
            convert_units = MDAnalysis.core.flags['convert_gromacs_lengths']
        self.convert_units = convert_units  # convert length and time to base units on the fly?
        self.xdrfile = None
        self.__numatoms = None
        self.__numframes = None  # takes a long time, avoid accessing self.numframes
        self.skip_timestep = 1   # always 1 for xdr files
        self.__delta = None      # compute from time in first two frames!
        self.fixed = 0           # not relevant for Gromacs xtc/trr
        self.skip = 1
        self.periodic = False
        self.ts = self._Timestep(self.numatoms)
        # Read in the first timestep
        self._read_next_timestep()

    @property
    def numatoms(self):
        """Read the number of atoms from the trajectory.

        The result is cached. If for any reason the trajectory cannot
        be read then 0 is returned.
        """
        if not self.__numatoms is None:   # return cached value
            return self.__numatoms
        try:
            self.__numatoms = self._read_trj_natoms(self.filename)
        except IOError:
            return 0
        else:
            return self.__numatoms

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

    def close_trajectory(self):
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
          *start*
              number of the first recorded MD step
          *step*
              indicate that *step* MD steps (!) make up one trajectory frame
          *delta*
              MD integrator time step (!), in AKMA units
          *precision*
              accuracy for lossy XTC format

        :Returns: :class:`DCDWriter`

        .. Note:: The keyword arguments set the low-level attributes of the DCD
                  according to the CHARMM format. The time between two frames
                  would be *delta* * *step* !
        """
        numatoms = kwargs.pop('numatoms', self.numatoms)
        kwargs.setdefault('start', self.start_timestep)
        kwargs.setdefault('step', self.skip_timestep)
        kwargs.setdefault('delta', self.delta)
        kwargs.setdefault('precision', self.precision)
        return TrjWriter(filename, numatoms, **kwargs)

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
                    self.close_trajectory()
                    raise
            except:
                self.close_trajectory()
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
            ts.status, ts.step, ts.time, ts.prec = libxdrfile.read_xtc(self.xdrfile, ts._unitcell, ts._pos)
        elif self.format == 'TRR':
            ts.status, ts.step, ts.time, ts.lmbda = libxdrfile.read_trr(self.xdrfile, ts._unitcell, ts._pos,
                                                                        ts._velocities, ts._forces)
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
        ts.frame += 1
        return ts

    def rewind(self):
        """Position at beginning of trajectory"""
        self._reopen()
        self.next()   # read first frame

    def _reopen(self):
        self.close_trajectory()
        self.open_trajectory()        

    def _forward_to_frame(self, frameindex):
        """Slow implementation: must read sequentially.

        .. Note:: *frameindex* starts from 0; i.e. *frameindex* =
                  frame - 1.

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
        self.close_trajectory()
        

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
        filename - name of output file
        numatoms - number of atoms in trajectory file
        start - starting timestep
        step  - skip between subsequent timesteps
        delta - timestep
        precision - accuracy for lossy XTC format [1000]
        convert_units - units are converted to the MDAnalysis base format; ``None`` selects
                        the value of :data:`MDAnalysis.core.flags`['convert_gromacs_lengths']
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
            ts - timestep object containing coordinates to be written to dcd file
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
        """Generic writer with minimum intelligence; override if necessary."""
        if self.convert_units:
            self.convert_pos_to_native(ts._pos)             # in-place !
            try:
                ts.time = self.convert_time_to_native(ts.time)
            except AttributeError:
                ts.time = ts.frame * self.convert_time_to_native(self.delta)
        if not hasattr(ts, 'step'):
            # bogus, should be actual MD step number, i.e. frame * delta/dt
            ts.step = ts.frame
        unitcell = self.convert_dimensions_to_unitcell(ts).astype(numpy.float32)  # must be float32 (!)
        if self.format == 'XTC':
            status = libxdrfile.write_xtc(self.xdrfile, ts.step, ts.time, unitcell, ts._pos, self.precision)
        elif self.format == 'TRR':
            if not hasattr(ts, 'lmbda'):
                ts.lmbda = 1.0
            if not hasattr(ts, '_velocities'):
                ts._velocities = numpy.zeros((3,ts.numatoms), dtype=numpy.float32)
            if not hasattr(ts, '_forces'):
                ts._forces = numpy.zeros((3,ts.numatoms), dtype=numpy.float32)
            status = libxdrfile.write_trr(self.xdrfile, ts.step, ts.time, ts.lmbda, unitcell, 
                                           ts._pos, ts._velocities, ts._forces)
        else:
            raise NotImplementedError("Gromacs trajectory format %s not known." % self.format)
        return status

    def close_trajectory(self):
        status = libxdrfile.exdrCLOSE
        if not self.xdrfile is None:
            status = libxdrfile.xdrfile_close(self.xdrfile)
            self.xdrfile = None
        return status

    def convert_dimensions_to_unitcell(self, ts):
        """Read dimensions from timestep *ts* and return Gromacs box vectors"""
        return self.convert_pos_to_native(triclinic_vectors(ts.dimensions))
