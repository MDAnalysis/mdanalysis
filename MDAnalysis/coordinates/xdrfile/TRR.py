"""Reading of Gromacs trr trajectories."""

import numpy

import core
import libxdrfile

class Timestep(core.Timestep):
    """Timestep for a Gromacs TRR trajectory."""
    def __init__(self, arg):
        DIM = libxdrfile.DIM    # compiled-in dimension (most likely 3)
        if numpy.dtype(type(arg)) == numpy.dtype(int):
            self.frame = 0
            self.numatoms = arg
            # C floats and C-order for arrays (see libxdrfile.i)
            self._pos = numpy.zeros((self.numatoms, DIM), dtype=numpy.float32, order='C')            
            self._velocities = numpy.zeros((self.numatoms, DIM), dtype=numpy.float32, order='C')            
            self._forces = numpy.zeros((self.numatoms, DIM), dtype=numpy.float32, order='C')            
            self._unitcell = numpy.zeros((DIM,DIM), dtype=numpy.float32)
            # additional data for xtc
            self.status = libxdrfile.exdrOK
            self.step = 0
            self.time = 0
            self.lmbda = 0
        elif isinstance(arg, Timestep): # Copy constructor
            # This makes a deepcopy of the timestep
            self.frame = arg.frame
            self.numatoms = arg.numatoms
            self._unitcell = numpy.array(arg._unitcell)
            self._pos = numpy.array(arg._pos)
            self._velocities = numpy.array(arg._velocities)
            self._forces = numpy.array(arg._forces)
            for attr in ('status', 'step', 'time', 'lmbda'):
                if hasattr(arg, attr):
                    self.__setattr__(attr, arg.__getattribute__(attr))
        elif isinstance(arg, numpy.ndarray):
            # provide packed array shape == (natoms, 3*DIM)
            # which contains pos = arg[:,0:3], v = arg[:,3:6], f = arg[:, 6:9]
            if len(arg.shape) != 2: raise Exception("packed numpy array (x,v,f) can only have 2 dimensions")
            self._unitcell = numpy.zeros((DIM,DIM), dtype=numpy.float32)
            self.frame = 0
            if arg.shape[0] == 3*DIM: 
                self.numatoms = arg.shape[-1]   # C-order
            elif arg.shape[1] == 3*DIM: 
                self.numatoms = arg.shape[0]    # C-order
            else:
                raise ValueError("TRR timestep is initialized from (natoms, 3*3) array")
            # OB -- not sure if this is really doing what i's supposed to, especially
            # when the order is really swapped in arg
            self._pos = arg[:,0:DIM].copy('C')               # C-order
            self._velocities = arg[:,DIM:2*DIM].copy('C')    # C-order
            self._pos = arg[:,2*DIM:3*DIM].copy('C')         # C-order
            # additional data for trr
            self.status = libxdrfile.exdrOK
            self.step = 0
            self.time = 0
            self.lmbda = 0
        else: 
            raise Exception("Cannot create an empty Timestep")
        self._x = self._pos[:,0]
        self._y = self._pos[:,1]
        self._z = self._pos[:,2]


class TRRReader(core.TrjReader):
    """Read a `Gromacs <www.gromacs.org>` TRR trajectory."""
    format = "TRR"
    _Timestep = Timestep

class TRRWriter(core.TrjWriter):
    """Write a `Gromacs <www.gromacs.org>` TRR trajectory."""
    format = "TRR"
