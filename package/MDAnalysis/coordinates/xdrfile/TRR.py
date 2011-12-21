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
Reading of `Gromacs TRR trajectories`_.

.. _Gromacs TRR trajectories: http://www.gromacs.org/Documentation/File_Formats/.trr_File
.. _Gromacs: http://www.gromacs.org


.. SeeAlso:: :mod:`MDAnalysis.coordinates.xdrfile.libxdrfile` for low-level
   bindings to the Gromacs trajectory file formats

Classes
-------

.. autoclass:: Timestep
   :members:
   :inherited-members:
.. autoclass:: TRRReader
   :members:
   :inherited-members:
.. autoclass:: TRRWriter
   :members:
   :inherited-members:
"""

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
            # or just positions: pos = arg[:,0:3] == arg
            if len(arg.shape) != 2:
                raise ValueError("packed numpy array (x,v,f) can only have 2 dimensions")
            self._unitcell = numpy.zeros((DIM,DIM), dtype=numpy.float32)
            self.frame = 0
            if (arg.shape[0] == 3*DIM and arg.shape[1] != 3*DIM) or \
                    (arg.shape[0] == DIM and arg.shape[1] != DIM):
                # wrong order (but need to exclude case where natoms == DIM or natoms == 3*DIM!)
                raise ValueError("TRR timestep is to be initialized from (natoms, 3*3) or (natoms, 3) array")
            self.numatoms = arg.shape[0]
            self._pos = arg[:,0:DIM].copy('C')                   # C-order
            if arg.shape[1] == 3*DIM:
                self._velocities = arg[:,DIM:2*DIM].copy('C')    # C-order
                self._forces = arg[:,2*DIM:3*DIM].copy('C')      # C-order
            elif arg.shape[1] == DIM:
                self._velocities = numpy.zeros_like(self._pos)
                self._forces = numpy.zeros_like(self._pos)
            else:
                raise ValueError("TRR timestep has not second dimension 3 or 9: shape=%r" % (arg.shape,))
            # additional data for trr
            self.status = libxdrfile.exdrOK
            self.step = 0
            self.time = 0
            self.lmbda = 0
        else:
            raise ValueError("Cannot create an empty Timestep")
        self._x = self._pos[:,0]
        self._y = self._pos[:,1]
        self._z = self._pos[:,2]

class TRRWriter(core.TrjWriter):
    """Write a Gromacs_ TRR trajectory."""
    format = "TRR"
    units = {'time': 'ps', 'length':'nm', 'velocity':'nm/ps', 'force':'kJ/(mol*nm)'}

class TRRReader(core.TrjReader):
    """Read a Gromacs_ TRR trajectory."""
    format = "TRR"
    _Timestep = Timestep
    _Writer = TRRWriter
    units = {'time': 'ps', 'length':'nm', 'velocity':'nm/ps', 'force':'kJ/(mol*nm)'}

