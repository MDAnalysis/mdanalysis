import sys
import os
import pyh5md
import numpy as np
import MDAnalysis as mda
from MDAnalysis.coordinates import base, core

class Timestep(base.Timestep):
    """H5MD Timestep

    """
    order = 'C'

    def _init_unitcell(self):
        return np.zeros((3, 3), dtype=np.float32)

    @property
    def dimensions(self):
        """unitcell dimensions (*A*, *B*, *C*, *alpha*, *beta*, *gamma*)

        lengths *a*, *b*, *c* are in the MDAnalysis length unit (Å), and
        angles are in degrees.

        Setting dimensions will populate the underlying native format description (triclinic box vectors)
        If `edges<https://nongnu.org/h5md/h5md.html#simulation-box>`_ is a matrix,
        the box is of triclinic shape with the edge vectors given by the rows of the matrix.)
        """
        return core.triclinic_box(*self._unitcell) # vectors are the rows

    @dimensions.setter
    def dimensions(self, box):
        self._unitcell[:] = core.triclinic_vectors(box)

class H5MDReader(base.ReaderBase):
    """Reader for the H5MD format.
    
    At the moment, all data is read directly from file without unit conversion.
    
    Data that is currently read from an H5MD file includes: n_frames, dimensions,
    positions, velocities, forces, data['step']
    
    Data that is not currently read from an H5MD file includes: masses and others
    """
    format = 'H5MD'
    units = {'time': None, 'length': None}
    # need to add units
    _Timestep = Timestep

    def __init__(self, filename, **kwargs):
        """
        Parameters
        ----------
        filename : str
            trajectory filename
        **kwargs : dict
            General reader arguments.

        """
        super(H5MDReader, self).__init__(filename, **kwargs)
        self.filename = filename
        self.open_trajectory() # also initializes particles_group('trajectory') as self._trajectory
        self.n_atoms = pyh5md.element(self._trajectory, 'n_atoms').value[()]
        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
        self._read_next_timestep()


    def open_trajectory(self) :
        """opens the trajectory file using pyh5md module"""
        self._frame = -1
        self._file = pyh5md.File(self.filename, 'r')
        self._trajectory = self._file.particles_group('trajectory')

    def close(self):
        """close reader"""
        self._file.close()

    @property
    def n_frames(self):
        """number of frames in trajectory"""
        return pyh5md.element(self._trajectory, 'positions').value.shape[0]

    def _reopen(self):
        """reopen trajectory"""
        self.close()
        self.open_trajectory()

    def _read_frame(self, frame):
        try:
            myframe = pyh5md.element(self._trajectory, 'positions').step[frame]
        except ValueError:
            raise IOError from None

        # set frame number
        self._frame = frame
        ts = self.ts

        # sets the Timestep object
        ts.frame = frame
        ts.data['step'] = pyh5md.element(self._trajectory, 'positions').step[frame]
        # will add more data

        # set frame box dimensions
        # set triclinic box vectors
        ts._unitcell[:] = pyh5md.element(self._trajectory['box'], 'edges').value[frame, :]

        # set particle positions
        frame_positions = pyh5md.element(self._trajectory, 'positions').value[frame, :]
        n_atoms_now = frame_positions.shape[0]
        if n_atoms_now != self.n_atoms :
            raise ValueError("Frame %d has %d atoms but the initial frame has %d"
                " atoms. MDAnalysis in unable to deal with variable"
                " topology!"%(frame, n_atoms_now, self.n_atoms))
        else :
            ts.positions = frame_positions

        # set particle velocities
        frame_velocities = pyh5md.element(self._trajectory, 'velocities').value[frame, :]
        n_atoms_now = frame_velocities.shape[0]
        if n_atoms_now != self.n_atoms :
            raise ValueError("Frame %d has %d atoms but the initial frame has %d"
                " atoms. MDAnalysis in unable to deal with variable"
                " topology!"%(frame, n_atoms_now, self.n_atoms))
        else :
            ts.velocities = frame_velocities

        # set particle forces
        frame_forces = pyh5md.element(self._trajectory, 'forces').value[frame, :]
        n_atoms_now = frame_forces.shape[0]
        if n_atoms_now != self.n_atoms :
            raise ValueError("Frame %d has %d atoms but the initial frame has %d"
                " atoms. MDAnalysis in unable to deal with variable"
                " topology!"%(frame, n_atoms_now, self.n_atoms))
        else :
            ts.forces = frame_forces

        return ts

    def _read_next_timestep(self) :
        """read next frame in trajectory"""
        return self._read_frame(self._frame + 1)
