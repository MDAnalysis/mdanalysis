"""DCD Hierarchy

"""
import numpy

class Timestep:
    """Timestep data for one frame

Data:     numatoms                   - number of atoms
          frame                      - frame number
          dimensions                 - system box dimensions (x, y, z, alpha, beta, gamma)

Methods:  t = Timestep(numatoms) - create a timestep object with space for numatoms (done automatically)
          t[i]                   - return coordinates for the i'th atom (0-based)
          t[start:stop:skip]     - return an array of coordinates, where start, stop and skip correspond to atom indices (0-based)
"""
    def __init__(self, arg):
        if numpy.dtype(type(arg)) == numpy.dtype(int):
            self.frame = 0
            self.numatoms = arg
            self._pos = numpy.zeros((self.numatoms, 3), dtype=numpy.float32, order='F')
            #self._pos = numpy.zeros((3, self.numatoms), numpy.float32)
            self._unitcell = numpy.zeros((6), numpy.float32)
        elif isinstance(arg, Timestep): # Copy constructor
            # This makes a deepcopy of the timestep
            self.frame = arg.frame
            self.numatoms = arg.numatoms
            self._unitcell = numpy.array(arg._unitcell)
            self._pos = numpy.array(arg._pos)
        elif isinstance(arg, numpy.ndarray):
            if len(arg.shape) != 2: raise Exception("numpy array can only have 2 dimensions")
            self._unitcell = numpy.zeros((6), numpy.float32)
            self.frame = 0
            if arg.shape[0] == 3: self.numatoms = arg.shape[0]
            else: self.numatoms = arg.shape[-1]
            self._pos = arg.copy('Fortran')
        else: raise Exception("Cannot create an empty Timestep")
        self._x = self._pos[:,0]
        self._y = self._pos[:,1]
        self._z = self._pos[:,2]
    def __getitem__(self, atoms):
        if numpy.dtype(type(atoms)) == numpy.dtype(int):
            if (atoms < 0):
                atoms = self.numatoms + atoms
            if (atoms < 0) or (atoms >= self.numatoms):
                raise IndexError
            return self._pos[atoms]
        elif type(atoms) == slice or type(atoms) == numpy.ndarray:
            return self._pos[atoms]
        else: raise TypeError
    def __len__(self):
        return self.numatoms
    def __iter__(self):
        def iterTS():
            for i in xrange(self.numatoms):
                yield self[i]
        return iterTS()
    def __repr__(self):
        return "< Timestep "+ repr(self.frame) + " with unit cell dimensions " + repr(self.dimensions) + " >"
    def copy(self):
        return self.__deepcopy__()
    def __deepcopy__(self):
        # Is this the best way?
        return Timestep(self)

    # Properties
    def dimensions():
        def fget(self):
            # Layout of unitcell is [A, alpha, B, beta, gamma, C]
            return numpy.take(self._unitcell, [0,2,5,1,3,4])
        return locals()
    dimensions = property(**dimensions())

class DCDWriter:
    """Writes to a DCD file

Data:

Methods:
    d = DCDWriter(dcdfilename, numatoms, start, step, delta, remarks)
"""
    def __init__(self, dcdfilename, numatoms, start=0, step=1, delta=1.0, remarks="Created by DCDWriter"):
        ''' Create a new DCDWriter
        dcdfilename - name of output file
        numatoms - number of atoms in dcd file
        start - starting timestep
        step  - skip between subsequent timesteps
        delta - timestep
        remarks - comments to annotate dcd file
        '''
        if numatoms == 0:
            raise Exception("DCDWriter: no atoms in output trajectory")
        self.dcdfilename = dcdfilename
        self.numatoms = numatoms

        self.frames_written = 0
        self.start = start
        self.step = step
        self.delta = delta
        self.dcdfile = file(dcdfilename, 'wb')
        self.remarks = remarks
        self._write_dcd_header(numatoms, start, step, delta, remarks)
    def dcd_header(self):
        import struct
        desc = ['file_desc', 'header_size', 'natoms', 'nsets', 'setsread', 'istart', 'nsavc', 'delta', 'nfixed', 'freeind_ptr', 'fixedcoords_ptr', 'reverse', 'charmm', 'first', 'with_unitcell']
        return dict(zip(desc, struct.unpack("iiiiiiidiPPiiii",self._dcd_C_str)))
    def write_next_timestep(self, ts=None):
        ''' write a new timestep to the dcd file
            ts - timestep object containing coordinates to be written to dcd file
        '''
        if ts is None:
            if not hasattr(self, "ts"):
                raise Exception("DCDWriter: no coordinate data to write to trajectory file")
            else:
                ts=self.ts
        # Check to make sure Timestep has the correct number of atoms
        elif not ts.numatoms == self.numatoms:
            raise Exception("DCDWriter: Timestep does not have the correct number of atoms")
        self._write_next_frame(ts._x, ts._y, ts._z, ts._unitcell)
        self.frames_written += 1
    def close_trajectory(self):
        # Do i really need this?
        self._finish_dcd_write()
        self.dcdfile.close()
        self.dcdfile = None
    def __del__(self):
        if self.dcdfile is not None:
            self.close_trajectory()
    def __repr__(self):
        return "< DCDWriter '"+ self.dcdfilename + "' for " + repr(self.numatoms) + " atoms >"

class DCDReader:
    """Reads from a DCD file
Data:
    ts                     - Timestep object containing coordinates of current frame

Methods:
    dcd = DCD(dcdfilename)             - open dcd file and read header
    len(dcd)                           - return number of frames in dcd
    for ts in dcd:                     - iterate through trajectory
    for ts in dcd[start:stop:skip]:    - iterate through a trajectory
    dcd[i]                             - random access into the trajectory (i corresponds to frame number)
    data = dcd.timeseries(...)         - retrieve a subset of coordinate information for a group of atoms
    data = dcd.correl(...)             - populate a Timeseries.Collection object with computed timeseries
"""
    def __init__(self, dcdfilename):
        self.dcdfilename = dcdfilename
        self.dcdfile = file(dcdfilename, 'rb')
        self.numatoms = 0
        self.numframes = 0
        self.fixed = 0
        self.skip = 1
        self.periodic = False
        self._read_dcd_header()
        self.ts = Timestep(self.numatoms)
        # Read in the first timestep
        self._read_next_timestep()
    def dcd_header(self):
        import struct
        desc = ['file_desc', 'header_size', 'natoms', 'nsets', 'setsread', 'istart', 'nsavc', 'delta', 'nfixed', 'freeind_ptr', 'fixedcoords_ptr', 'reverse', 'charmm', 'first', 'with_unitcell']
        return dict(zip(desc, struct.unpack("iiiiiiidiPPiiii",self._dcd_C_str)))
    def __len__(self):
        return self.numframes
    def __iter__(self):
        # Reset the trajectory file, read from the start
        # usage is "from ts in dcd:" where dcd does not have indexes
        self._reset_dcd_read()
        def iterDCD():
            for i in xrange(0, self.numframes, self.skip):
                try: yield self._read_next_timestep()
                except IOError: raise StopIteration
        return iterDCD()
    def __getitem__(self, frame):
        if (numpy.dtype(type(frame)) != numpy.dtype(int)) and (type(frame) != slice):
            raise TypeError
        if (numpy.dtype(type(frame)) == numpy.dtype(int)):
            if (frame < 0):
                # Interpret similar to a sequence
                frame = len(self) + frame
            if (frame < 0) or (frame >= len(self)):
                raise IndexError
            self._jump_to_frame(frame)
            ts = self.ts
            ts.frame = self._read_next_frame(ts._x, ts._y, ts._z, ts._unitcell, 1)
            return ts
        elif type(frame) == slice: # if frame is a slice object
            if not (((type(frame.start) == int) or (frame.start == None)) and
                    ((type(frame.stop) == int) or (frame.stop == None)) and
                    ((type(frame.step) == int) or (frame.step == None))):
                raise TypeError("Slice indices are not integers")
            def iterDCD(start=frame.start, stop=frame.stop, step=frame.step):
                start, stop, step = self._check_slice_indices(start, stop, step)
                for i in xrange(start, stop, step):
                    yield self[i]
            return iterDCD()
    def _check_slice_indices(self, start, stop, skip):
        if start == None: start = 0
        if stop == None: stop = len(self)
        if skip == None: skip = 1
        if (start < 0): start += len(self)
        if (stop < 0): stop += len(self)
        elif (stop > len(self)): stop = len(self)
        if (stop <= start): raise Exception("Stop frame is lower than start frame")
        if ((start < 0) or (start >= len(self)) or
           (stop < 0) or (stop > len(self))):
               raise IndexError
        return (start, stop, skip)
    def _read_next_timestep(self, ts=None):
        if ts is None: ts = self.ts
        ts.frame = self._read_next_frame(ts._x, ts._y, ts._z, ts._unitcell, self.skip)
        return ts
    def timeseries(self, asel, start=0, stop=-1, skip=1, format='afc'):
        ''' Return a subset of coordinate data for an AtomGroup
            asel - AtomGroup object
            start, stop, skip - range of trajectory to access, start and stop are inclusive
            format - the order/shape of the return data array, corresponding to (a)tom, (f)rame, (c)oordinates
                     all six combinations of 'a', 'f', 'c' are allowed
                     ie "fac" - return array where the shape is (frame, number of atoms, coordinates)
        '''
        start, stop, skip = self._check_slice_indices(start, stop, skip)
        if len(asel) == 0:
            raise Exception("Timeseries requires at least one atom to analyze")
        if len(format) != 3 and format not in ['afc', 'acf', 'caf', 'cfa', 'fac', 'fca']:
            raise Exception("Invalid timeseries format")
        atom_numbers = list(asel.indices())
        # Check if the atom numbers can be grouped for efficiency, then we can read partial buffers
        # from trajectory file instead of an entire timestep
        # XXX needs to be implemented
        return self._read_timeseries(atom_numbers, start, stop, skip, format)
    def correl(self, timeseries, start=0, stop=-1, skip=1):
        ''' Populate a TimeseriesCollection object with timeseries computed from the trajectory
            timeseries - TimeseriesCollection
            start, stop, skip - subset of trajectory to use, with start and stop being inclusive
        '''
        start, stop, skip = self._check_slice_indices(start, stop, skip)
        atomlist = timeseries._getAtomList()
        format = timeseries._getFormat()
        lowerb, upperb = timeseries._getBounds()
        sizedata = timeseries._getDataSize()
        atomcounts = timeseries._getAtomCounts()
        auxdata = timeseries._getAuxData()
        return self._read_timecorrel(atomlist, atomcounts, format, auxdata, sizedata, lowerb, upperb, start, stop, skip)
    def close_trajectory(self):
        self._finish_dcd_read()
        self.dcdfile.close()
        self.dcdfile = None
    def __del__(self):
        if self.dcdfile is not None:
            self.close_trajectory()
    def __repr__(self):
            return "< DCDReader '"+ self.dcdfilename + "' with " + repr(self.numframes) + " frames of " + repr(self.numatoms) + " atoms (" + repr(self.fixed) + " fixed) >"

# Add the c functions to their respective classes so they act as class methods
import _dcdmodule
import new
DCDReader._read_dcd_header = new.instancemethod(_dcdmodule.__read_dcd_header, None, DCDReader)
DCDReader._read_next_frame = new.instancemethod(_dcdmodule.__read_next_frame, None, DCDReader)
DCDReader._jump_to_frame = new.instancemethod(_dcdmodule.__jump_to_frame, None, DCDReader)
DCDReader._reset_dcd_read = new.instancemethod(_dcdmodule.__reset_dcd_read, None, DCDReader)
DCDReader._finish_dcd_read = new.instancemethod(_dcdmodule.__finish_dcd_read, None, DCDReader)
DCDReader._read_timeseries = new.instancemethod(_dcdmodule.__read_timeseries, None, DCDReader)

DCDWriter._write_dcd_header = new.instancemethod(_dcdmodule.__write_dcd_header, None, DCDWriter)
DCDWriter._write_next_frame = new.instancemethod(_dcdmodule.__write_next_frame, None, DCDWriter)
DCDWriter._finish_dcd_write = new.instancemethod(_dcdmodule.__finish_dcd_write, None, DCDWriter)
del(_dcdmodule)

# _dcdtest is implemented with Pyrex - hopefully all dcd reading functionality can move to pyrex
import _dcdtest
#DCDReader._read_timeseries = new.instancemethod(_dcdtest.__read_timeseries, None, DCDReader)
DCDReader._read_timecorrel = new.instancemethod(_dcdtest.__read_timecorrel, None, DCDReader)
del(_dcdtest)
del(new)
