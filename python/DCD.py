"""DCD Hierarchy

"""
import Numeric

class Timestep:
    """Timestep data for one frame

Data: numatoms, frame, x, y, z, unitcell

Methods:

"""
    def __init__(self, arg):
        if type(arg) == int:
            self.frame = 0
            self.numatoms = arg
            self._pos = Numeric.zeros((3, self.numatoms), Numeric.Float32)
            self._unitcell = Numeric.zeros((6), Numeric.Float32)
        elif isinstance(arg, Timestep): # Copy constructor
            # This makes a deepcopy of the timestep passed in, transposing it if necessary
            self.frame = arg.frame
            self.numatoms = arg.numatoms
            self._unitcell = Numeric.array(arg._unitcell)
            if arg._pos.shape[0] != 3:
                self._pos = Numeric.array(Numeric.transpose(arg._pos))
            else:
                self._pos = Numeric.array(arg._pos)
        elif isinstance(arg, Numeric.arraytype):
            self._unitcell = Numeric.zeros((6), Numeric.Float32)
            self.frame = 0
            if arg.shape[0] != 3:
                self.numatoms = arg.shape[0]
                self._pos = Numeric.asarray(Numeric.transpose(arg))
            else:
                self.numatoms = arg.shape[-1]
                self._pos = Numeric.asarray(arg)
        else: raise Exception("Can't have an empty Timestep")
        self._x = self._pos[0]
        self._y = self._pos[1]
        self._z = self._pos[2]
    def __getattr__(self, name):
        if (name == "dimensions"):
            # Layout of unitcell is [A, alpha, B, beta, gamma, C]
            uc = self._unitcell
            return Numeric.take(uc, [0,2,5,1,3,4])
        else: raise AttributeError("class "+repr(self.__class__.__name__)+" has no attribute "+ name)
    def __getitem__(self, atomno):
        if (type(atomno) != int):
            raise TypeError
        if (atomno < 0):
            atomno = self.numatoms + atomno
        if (atomno < 0) or (atomno >= self.numatoms):
            raise IndexError
        return self._pos[:,atomno]
    def __len__(self):
        return self.numatoms
    def __iter__(self):
        def iterTS():
            for i in xrange(self.numatoms):
                yield self[i]
        return iterTS()
    def __str__(self):
        return repr(self)
    def __repr__(self):
        return "< Timestep "+ repr(self.frame) + " with unit cell dimensions " + repr(self.dimensions) + " >"
    def copy(self):
        return self.__deepcopy__()
    def __deepcopy__(self):
        # Is this the best way?
        return Timestep(self)

class DCDWriter:
    """Writes to a DCD file

Data:

Methods:
    d = DCDWriter(dcdfile)
"""
    def __init__(self, dcdfilename, numatoms, start=0, step=1, delta=1.0, remarks="Created by DCDWriter"):
        if numatoms == 0:
            raise Exception("DCDWriter: no atoms in trajectory file")
        self.dcdfilename = dcdfilename
        self.numatoms = numatoms

        self.numframes = 0
        self.start = start
        self.step = step
        self.delta = delta
        self.dcdfile = file(dcdfilename, 'wb')
        self.remarks = remarks
        self._write_dcd_header(numatoms, start, step, delta, remarks)
    def _dcdhandle(self):
        import struct
        desc = ['file_desc', 'header_size', 'natoms', 'nsets', 'setsread', 'istart', 'nsavc', 'delta', 'nfixed', 'freeind_ptr', 'fixedcoords_ptr', 'reverse', 'charmm', 'first', 'with_unitcell']
        return dict(zip(desc, struct.unpack("iiiiiiidiPPiiii",self._dcd_C_str)))
    def write_next_timestep(self, ts=None):
        if ts is None:
            if not hasattr(self, "ts"):
                raise Exception("DCDWriter: no coordinate data to write to trajectory file")
            else:
                ts=self.ts
        # Check to make sure Timestep has the correct number of atoms
        elif not ts.numatoms == self.numatoms:
            raise Exception("Timestep does not have the correct number of atoms")
        self._write_next_frame(ts._x, ts._y, ts._z, ts._unitcell)
        self.numframes += 1
    def close_trajectory(self):
        # Do i really need this?
        self._finish_dcd_write()
        self.dcdfile.close()
        self.dcdfile = None
    def __del__(self):
        if self.dcdfile is not None:
            self.close_trajectory()
    def __repr__(self):
        return "< DCDWriter '"+ self.dcdfilename + "' with " + repr(self.numframes) + " frames of " + repr(self.numatoms) + " atoms >"

class DCDReader:
    """Reads from a DCD file

Data:

Methods:
    d = DCD(dcdfile) - open dcdfile and read header
    len(d) - return number of frames

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
        self.read_next_timestep()
    def __dcdhandle(self):
        import struct
        desc = ['file_desc', 'header_size', 'natoms', 'nsets', 'setsread', 'istart', 'nsavc', 'delta', 'nfixed', 'freeind_ptr', 'fixedcoords_ptr', 'reverse', 'charmm', 'first', 'with_unitcell']
        return dict(zip(desc, struct.unpack("iiiiiiidiPPiiii",self._dcd_C_str)))
    def __getitem__(self, frame):
        if (type(frame) != int) and (type(frame) != slice):
            raise TypeError
        if (type(frame) == int):
            if (frame < 0):
                # Interpret similar to a sequence
                frame = len(self) + frame
            if (frame < 0) or (frame >= len(self)):
                raise IndexError
            self._jump_to_frame(frame)
            ts = self.ts
            ts.frame = self._read_next_frame(ts._x, ts._y, ts._z, ts._unitcell, 1)
            return ts
        else: # if frame is a slice object
            if not (((type(frame.start) == int) or (frame.start == None)) and
                    ((type(frame.stop) == int) or (frame.stop == None)) and
                    ((type(frame.step) == int) or (frame.step == None))):
                raise TypeError("Slicing indices are not integers or None")
            # Does it make sense to support slice objects?
            # XXX Maybe one day as a generator
            def iterDCD(start=frame.start, stop=frame.stop, step = frame.step):
                if (start < 0): start += len(self)
                if (stop < 0): stop += len(self)
                if (stop <= start): raise Exception("Stop frame is lower than start frame")
                if ((start < 0) or (start >= len(self)) or
                   (stop < 0) or (stop > len(self))):
                       raise IndexError
                if (step == None): step = 1
                for i in xrange(start, stop, step):
                    yield self[i]
            return iterDCD()
    #def reset_dcd_read(self):
    #    self._reset_dcd_read()
    def __iter__(self):
        # Reset the trajectory file
        self._reset_dcd_read()
        def iterDCD():
            for i in xrange(0, self.numframes, self.skip):
                try:
                    #self._jump_to_frame(i)
                    yield self.read_next_timestep()
                except IOError:
                    raise StopIteration
        return iterDCD()
    def read_next_timestep(self, ts=None):
        if (ts==None):
            ts = self.ts
        ts.frame = self._read_next_frame(ts._x, ts._y, ts._z, ts._unitcell, self.skip)
        return ts
    def timeseries(self, asel, start=0, stop=-1, skip=1, format='afc'):
        if (start < 0): start += len(self)
        if (stop < 0): stop += len(self)
        if (stop <= start): raise Exception("Stop frame is lower than start frame")
        if ((start < 0) or (start >= len(self)) or
           (stop < 0) or (stop >= len(self))):
               raise IndexError
        if len(asel) == 0:
            raise Exception("timeseries requires at least one atom to analyze")
        if len(asel) == 0:
            raise Exception("timeseries requires at least one atom to analyze")
        atom_numbers = list(asel.indices())
        # Check if the atom numbers can be grouped for efficiency
        # XXX needs to be implemented
        #results = []
        #for a in atom_numbers:
        #    results.append(self._read_correl([a], skip))
        #import Numeric
        #return Numeric.concatenate(tuple(results), 0)
        return self._read_timeseries(atom_numbers, start, stop, skip, format)
    def correl(self, timeseries, start=0, stop=-1, skip=1):
        if (start < 0): start += len(self)
        if (stop < 0): stop += len(self)
        if (stop <= start): raise Exception("Stop frame is lower than start frame")
        if ((start < 0) or (start >= len(self)) or
           (stop < 0) or (stop >= len(self))):
               raise IndexError
        atomlist = timeseries.getAtomList()
        format = timeseries.getFormat()
        lowerb, upperb = timeseries.getBounds()
        sizedata = timeseries.getDataSize()
        atomcounts = timeseries.getAtomCounts()
        auxdata = timeseries.getAuxData()
        return self._read_timecorrel(atomlist, atomcounts, format, auxdata, sizedata, lowerb, upperb, start, stop, skip)
    def __len__(self):
        return self.numframes
    def close_trajectory(self):
        self._finish_dcd_read()
        self.dcdfile.close()
        self.dcdfile = None
    def __del__(self):
        if self.dcdfile is not None:
            self.close_trajectory()
    def __repr__(self):
            return "< DCDReader '"+ self.dcdfilename + "' with " + repr(self.numframes) + " frames of " + repr(self.numatoms) + " atoms (" + repr(self.fixed) + " fixed) >"

import _dcd
import new
DCDReader._read_dcd_header = new.instancemethod(_dcd.__read_dcd_header, None, DCDReader)
DCDReader._read_next_frame = new.instancemethod(_dcd.__read_next_frame, None, DCDReader)
DCDReader._jump_to_frame = new.instancemethod(_dcd.__jump_to_frame, None, DCDReader)
DCDReader._reset_dcd_read = new.instancemethod(_dcd.__reset_dcd_read, None, DCDReader)
DCDReader._finish_dcd_read = new.instancemethod(_dcd.__finish_dcd_read, None, DCDReader)
DCDReader._read_timeseries = new.instancemethod(_dcd.__read_timeseries, None, DCDReader)

DCDWriter._write_dcd_header = new.instancemethod(_dcd.__write_dcd_header, None, DCDWriter)
DCDWriter._write_next_frame = new.instancemethod(_dcd.__write_next_frame, None, DCDWriter)
DCDWriter._finish_dcd_write = new.instancemethod(_dcd.__finish_dcd_write, None, DCDWriter)
del(_dcd)

import _dcdtest
#DCDReader._read_timeseries = new.instancemethod(_dcdtest.__read_timeseries, None, DCDReader)
DCDReader._read_timecorrel = new.instancemethod(_dcdtest.__read_timecorrel, None, DCDReader)
del(_dcdtest)
del(new)
