import h5py
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

    Currently reads .h5md files with the following HDF5 hierarchy:
    Notation - (name) is an HDF5 group and [variable] is an HDF5 dataset,
               <dtype> is dataset datatype

    H5MD root
     \-- (h5md)
     \-- (particles)
        \-- (trajectory)
            \-- (box)
                \-- (edges)
                    \-- [step] <int32>, gives frame#
                    \-- [value] <float>, gives box dimensions
            \-- (positions)
                \-- [step] <int32>, gives frame #
                \-- [time] <float>, gives time
                \-- [value] <float>, gives trajectory positions
            \-- (velocities)
                \-- [step] <int32>, gives frame #
                \-- [time] <float>, gives time
                \-- [value] <float>, gives trajectory velocities
            \-- (forces)
                \-- [step] <int32>, gives frame #
                \-- [time] <float>, gives time
                \-- [value] <float>, gives trajectory forces
            \-- (data)
                \-- (dt)
                    \-- [step] <int32>, gives frame #
                    \-- [value] <float>, gives dt
                \-- (lambda)
                    \-- [step] <int32>, gives frame #
                    \-- [value] <float>, gives lambda
                \-- (step)
                    \-- [step] <int32>, gives frame #
                    \-- [value] <int32>, gives step
            \-- [n_atoms] <int>, gives # of atoms in trajectory

    Data that is not currently read from an H5MD file includes: masses
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
        self.open_trajectory()
        self.n_atoms = self._file['particles']['trajectory']['n_atoms'][()]
        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
        self._read_next_timestep()


    def open_trajectory(self) :
        """opens the trajectory file using h5py package"""
        self._frame = -1
        self._file = h5py.File(self.filename, 'r')

    def close(self):
        """close reader"""
        self._file.close()

    @property
    def n_frames(self):
        """number of frames in trajectory"""
        return self._file['particles']['trajectory']['positions']['value'].shape[0]

    def _reopen(self):
        """reopen trajectory"""
        self.close()
        self.open_trajectory()

    def _read_frame(self, frame):
        try:
            myframe = self._file['particles']['trajectory']['positions']['step'][frame]
        except ValueError:
            raise IOError from None

        # set frame number
        self._frame = frame
        ts = self.ts

        # sets the Timestep object
        ts.frame = frame

        # set data dictionary values
        ts.data['time'] = self._file['particles']['trajectory']['positions']['time'][frame]
        ts.data['step'] = self._file['particles']['trajectory']['data']['step']['value'][frame]
        ts.data['lambda'] = self._file['particles']['trajectory']['data']['lambda']['value'][frame]
        ts.data['dt'] = self._file['particles']['trajectory']['data']['dt']['value'][frame]

        # set frame box dimensions
        # set triclinic box vectors
        ts._unitcell[:] = self._file['particles']['trajectory']['box']['edges']['value'][frame, :]

        # set particle positions
        frame_positions = self._file['particles']['trajectory']['positions']['value'][frame, :]
        n_atoms_now = frame_positions.shape[0]
        if n_atoms_now != self.n_atoms :
            raise ValueError("Frame %d has %d atoms but the initial frame has %d"
                " atoms. MDAnalysis in unable to deal with variable"
                " topology!"%(frame, n_atoms_now, self.n_atoms))
        else :
            ts.positions = frame_positions

        # set particle velocities
        frame_velocities = self._file['particles']['trajectory']['velocities']['value'][frame, :]
        n_atoms_now = frame_velocities.shape[0]
        if n_atoms_now != self.n_atoms :
            raise ValueError("Frame %d has %d atoms but the initial frame has %d"
                " atoms. MDAnalysis in unable to deal with variable"
                " topology!"%(frame, n_atoms_now, self.n_atoms))
        else :
            ts.velocities = frame_velocities

        # set particle forces
        frame_forces = self._file['particles']['trajectory']['forces']['value'][frame, :]
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
