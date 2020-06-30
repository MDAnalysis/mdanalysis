import numpy as np
import sys
import os
import pyh5md
import MDAnalysis as mda
from MDAnalysis.coordinates import base

class H5MDReader(base.ReaderBase):
    """Reader for the H5MD format.

    """
    format = 'H5MD'
    units = {'time': None, 'length': None}
    
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
        self.open_trajectory()
        self.n_atoms = self._file['particles']['atoms']['n_atoms'][()]
        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
        self._read_next_timestep()
        
    def open_trajectory(self) :
        """opens the trajectory file using pyh5md module"""
        self._frame = -1 # why last frame?
        self._file = pyh5md.File(self.filename, 'r')

    def close(self):
        """close reader"""
        self._file.close()
        
    @property
    def n_frames(self):
        """number of frames in trajectory"""
        return pyh5md.element(atoms,'positions').value.shape[0]

    def _reopen(self):
        """reopen trajectory"""
        self.close()
        self.open_trajectory()

    def _read_frame(self, frame):

        # set frame number
        self._frame = frame

        # sets the Timestep object
        self.ts.frame = frame
        #self.ts.data['step'] = myframe.configuration.step

        # set frame box dimensions
        self.ts.dimensions = self._file['particles']['atoms']['box']['edges']['value'][frame,:,:] #triclinic dimension
        #for i in range(3,6) :
            #self.ts.dimensions[i] = np.arccos(self.ts.dimensions[i]) * 180.0 / np.pi

        # set particle positions
        frame_positions = self._file['particles']['atoms']['positions']['value'][frame,:,:]
        n_atoms_now = frame_positions.shape[0]
        if n_atoms_now != self.n_atoms :
            raise ValueError("Frame %d has %d atoms but the initial frame has %d"
                " atoms. MDAnalysis in unable to deal with variable"
                " topology!"%(frame, n_atoms_now, self.n_atoms))
        else :
            self.ts.positions = frame_positions
            
        # set particle velocities
        frame_velocities = self._file['particles']['atoms']['velocities']['value'][frame,:,:]
        n_atoms_now = frame_velocities.shape[0]
        if n_atoms_now != self.n_atoms :
            raise ValueError("Frame %d has %d atoms but the initial frame has %d"
                " atoms. MDAnalysis in unable to deal with variable"
                " topology!"%(frame, n_atoms_now, self.n_atoms))
        else :
            self.ts.velocities = frame_velocities
            
        # set particle forces
        frame_forces = self._file['particles']['atoms']['forces']['value'][frame,:,:]
        n_atoms_now = frame_forces.shape[0]
        if n_atoms_now != self.n_atoms :
            raise ValueError("Frame %d has %d atoms but the initial frame has %d"
                " atoms. MDAnalysis in unable to deal with variable"
                " topology!"%(frame, n_atoms_now, self.n_atoms))
        else :
            self.ts.forces = frame_forces
            
        # set particle masses
        frame_masses = self._file['particles']['atoms']['masses'][()]
        n_atoms_now = frame_masses.shape
        if n_atoms_now != self.n_atoms :
            raise ValueError("Frame %d has %d atoms but the initial frame has %d"
                " atoms. MDAnalysis in unable to deal with variable"
                " topology!"%(frame, n_atoms_now, self.n_atoms))
        else :
            self.ts.masses = frame_masses
            
        return self.ts

    def _read_next_timestep(self) :
        """read next frame in trajectory"""
        return self._read_frame(self._frame + 1)