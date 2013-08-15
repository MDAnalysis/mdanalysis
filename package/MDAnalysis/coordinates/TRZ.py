import os, errno
import base
from base import Timestep
import MDAnalysis.core
import numpy
import struct

class Timestep(base.Timestep):
    """ TRZ custom Timestep
    """
    @property
    def dimensions(self):
        """Unit cell dimensions, (A, B, C, alpha, beta, gamma)
        """
        return self._unitcell

class TRZReader(base.Reader):
    """ Reads an IBIsCO or YASP trajectory file
    TRZ format detailed below, each line is a single fortran write statement, so is surrounded by 4 bytes when written by fortran
    In brackets after each entry is the size of the content of each line
    Header :title(80c)
            nrec (int4)
    Frame  :nframe, ntrj*nframe, natoms, treal (3*int4, real8)
            boxx, 0.0, 0.0, 0.0, boxy, 0.0, 0.0, 0.0, boxz (real8 * 9)
            pressure, pt11, pt12, pt22, pt13, pt23, pt33 (real8 *7)
            6, etot, ptot, ek, t, 0.0, 0.0 (int4, real8 * 6)
            rx (real4 * natoms)
            ry 
            rz
            vx
            vy
            vz
    """

    format = "TRZ"

    units = {'time':'ps', 'length':'nm'}

    def __init__(self, trzfilename, numatoms=None, **kwargs):
        if numatoms is None:
            raise ValueError('TRZReader requires the numatoms keyword')

        self.filename = trzfilename
        self.trzfile = file(self.filename, 'rb') #file() is like open()

        self.__numatoms = numatoms
        self.__numframes = None

        self.fixed = 0 #Are any atoms fixed in place? Not used in trz files
        self.skip = 1 #Step size for iterating through trajectory
        self.periodic = False # Box info for PBC
        self.delta = kwargs.pop("delta", 1.0) # Length of timestep
        self.skip_timestep = 1 # Number of steps between frames, can be found in trz files

        self._read_trz_header()
        self.ts = Timestep(self.numatoms)
        self._read_next_timestep()

    def _read_trz_header(self):
        #Read the header of the file
        #file.read(4)
        #title = struct.unpack('80c',file.read(80))
        #file.read(4)
        #file.read(4)
        #nrec = struct.unpack('i',file.read(4))
        #file.read(4)
        self.trzfile.seek(100,1) # Header is 100 bytes in total, but contains nothing "useful"

    def _read_next_timestep(self, ts=None): # self.next() is from base Reader class and calls this
        #Read a timestep from binary file
        if ts is None:
            ts = self.ts

        try:
            #Skip into frame
            self.trzfile.seek(12,1) # Seek forward 12 bytes from current position (1)
            natoms = struct.unpack('i',self.trzfile.read(4))[0]
            ts.time = struct.unpack('d',self.trzfile.read(8)) # Real time of the system
            #Read box data
            self.trzfile.seek(8,1) #Includes 4 from previous write statement
            box = list(struct.unpack('9d',self.trzfile.read(72)))
            ts.dimensions[0] = box[0] # lx, 0.0, 0.0
            ts.dimensions[1] = box[4] # 0.0, ly, 0.0
            ts.dimensions[2] = box[-1]# 0.0, 0.0, lz
            ts.dimensions[3:] = [90., 90., 90.] #Assume box is orthagonal
            self.trzfile.seek(4,1)
            #Skip some more stuff
            midskip = (4 + 8*7 + 4) + (4 + 4 + 8*6 + 4) #Skip energy and pressure entries
            self.trzfile.seek(midskip,1)
            #Read coordinate data
            readarg = str(natoms) + 'f'
            self.trzfile.seek(4,1)
            ts._x[:] = struct.unpack(readarg,self.trzfile.read(4*natoms)) # ts._pos[:,0] = x coord
            self.trzfile.seek(8,1)
            ts._y[:] = struct.unpack(readarg,self.trzfile.read(4*natoms))
            self.trzfile.seek(8,1)
            ts._z[:] = struct.unpack(readarg,self.trzfile.read(4*natoms))
            self.trzfile.seek(4,1)
            #Skip velocities at end of frame
            # ts._velocities[:,0] is vx
            endskip = 3*(4 + natoms*4 + 4) # 3 real4 arrays natoms in length bookended by 4 bytes, velocity info
            self.trzfile.seek(endskip,1)

            ts.frame += 1
            return ts
        except struct.error: #End of file reached if struct fails to unpack
            raise IOError

    @property
    def numatoms(self):
        """Number of atoms in a frame"""
        if not self.__numatoms is None:
            return  self.__numatoms
        try:
            self._reopen()
            self.__numatoms = self._read_trz_natoms(self.trzfile)
        except IOError:
            return 0
        else:
            return self.__numatoms

    def _read_trz_natoms(self, f):
        #Read start of next frame and reopen file
        try:
            f.seek(12,1) #Reads 4 bytes at start, then nframe, ntrj*nframe
            natoms = struct.unpack('i',f.read(4))
        except struct.error:
            raise IOError
        else:
            self._reopen()
            return natoms

    @property
    def numframes(self):
        """Total number of frames in a trajectory"""
        if not self.__numframes is None:
            return self.__numframes
        try:
            self.__numframes = self._read_trz_numframes(self.trzfile)
        except IOError:
            return 0
        else:
            return self.__numframes

    def _read_trz_numframes(self, trzfile):
        framecounter = 0
        self._reopen()
        while True:
            try:
                self._read_next_timestep()
                framecounter += 1
            except IOError:
                self.rewind()
                return framecounter

    def __iter__(self):
        self.ts.frame = 0
        self._reopen()
        while True:
            try:
                yield self._read_next_timestep()
            except IOError:
                self.rewind()
                raise StopIteration

    def rewind(self):
        """Reposition onto first frame"""
        self._reopen()
        self.next()

    def _reopen(self):
        self.close()
        self.open_trajectory()
        self._read_trz_header() # Moves to start of first frame

    def open_trajectory(self):
        """Open the trajectory file"""
        if not self.trzfile is None:
            raise IOError(errno.EALREADY, 'TRZ file already opened', self.filename)
        if not os.path.exists(self.filename):
            raise IOError(errno.ENOENT, 'TRZ file not found', self.filename)

        self.trzfile = file(self.filename, 'rb')

        #Reset ts
        ts = self.ts
        ts.status = 1
        ts.frame =0
        ts.step = 0
        ts.time = 0
        return self.trzfile

    def close(self):
        """Close trz file if it was open"""
        if self.trzfile is None:
            return
        self.trzfile.close()
        self.trzfile = None

    def __del__(self):
        if not self.trzfile is None:
            self.close()

