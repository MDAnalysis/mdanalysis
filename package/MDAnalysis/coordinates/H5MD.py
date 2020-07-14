import numpy as np
import MDAnalysis as mda
from MDAnalysis.coordinates import base, core
try:
    import h5py
except ImportError:
    raise ImportError("Cannot import h5py")

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

    Notation
    --------
    (name) is an HDF5 group that the reader recognizes
    {name} is an HDF5 group with arbitrary name
    [variable] is an HDF5 dataset
    <dtype> is dataset datatype
    +-- is an attribute of a group or dataset
    See `h5md documentation<https://nongnu.org/h5md/h5md.html#>`_
    for detailed overview

    H5MD root
     \-- (h5md)
     \-- (particles)
        \-- {group1}
            \-- (box)
                +-- dimension <int>, gives box spatial dimension
                +-- boundary <str>, gives boundary conditions
                \-- (edges)
                    \-- [step] <int32>, gives frame#
                    \-- [value] <float>, gives box dimensions
            \-- (position)
                \-- [step] <int>, gives frame #
                \-- [time] <float>, gives time
                    +-- units <str>
                \-- [value] <float>, gives trajectory positions
                    +-- units <str>
            \-- (velocity)
                \-- [step] <int>, gives frame #
                \-- [time] <float>, gives time
                    +-- units <str>
                \-- [value] <float>, gives trajectory velocities
                    +-- units <str>
            \-- (force)
                \-- [step] <int>, gives frame #
                \-- [time] <float>, gives time
                    +-- units <str>
                \-- [value] <float>, gives trajectory forces
                    +-- units <str>
            \-- (data)
                \-- (dt)
                    \-- [step] <int>, gives frame #
                    \-- [value] <float>, gives dt
                \-- (lambda)
                    \-- [step] <int>, gives frame #
                    \-- [value] <float>, gives lambda
                \-- (step)
                    \-- [step] <int>, gives frame #
                    \-- [value] <int>, gives step

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
        self.n_atoms = self._particle_group['position/value'].shape[1]
        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
        self._read_next_timestep()


    def open_trajectory(self) :
        """opens the trajectory file using h5py library"""
        self._frame = -1
        self._file = h5py.File(self.filename, 'r')
        # pulls first key out of 'particles'
        # allows for arbitrary name of group1 in 'particles'
        self._particle_group = self._file['particles'][list(self._file['particles'])[0]]

    def close(self):
        """close reader"""
        self._file.close()

    @property
    def n_frames(self):
        """number of frames in trajectory"""
        return self._particle_group['position/value'].shape[0]

    def _reopen(self):
        """reopen trajectory"""
        self.close()
        self.open_trajectory()

    def _read_frame(self, frame):
        try:
            myframe = self._particle_group['position/step'][frame]
        except ValueError:
            raise IOError from None

        self._frame = frame
        ts = self.ts
        particle_group = self._particle_group
        ts.frame = frame

        # set data dictionary values
        ts.data['time'] = particle_group['position/time'][frame]
        ts.data['step'] = particle_group['data/step/value'][frame]
        ts.data['lambda'] = particle_group['data/lambda/value'][frame]
        ts.data['dt'] = particle_group['data/dt/value'][frame]

        # set frame box dimensions
        # set triclinic box vectors
        ts._unitcell[:] = particle_group['box/edges/value'][frame, :]

        # set particle positions
        frame_positions = particle_group['position/value'][frame, :]
        n_atoms_now = frame_positions.shape[0]
        if n_atoms_now != self.n_atoms :
            raise ValueError("Frame %d has %d atoms but the initial frame has %d"
                             " atoms. MDAnalysis in unable to deal with variable"
                             " topology!"%(frame, n_atoms_now, self.n_atoms))
        else :
            ts.positions = frame_positions

        # set particle velocities
        frame_velocities = particle_group['velocity/value'][frame, :]
        n_atoms_now = frame_velocities.shape[0]
        if n_atoms_now != self.n_atoms :
            raise ValueError("Frame %d has %d atoms but the initial frame has %d"
                             " atoms. MDAnalysis in unable to deal with variable"
                             " topology!"%(frame, n_atoms_now, self.n_atoms))
        else :
            ts.velocities = frame_velocities

        # set particle forces
        frame_forces = particle_group['force/value'][frame, :]
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

    def _read_h5md_units(self):
        """converts units from H5MD to MDAnalysis notation"""
        velocity_dict = {
        'nm ps-1': 'nm/ps'
        }
        force_dict = {
        'kJ mol-1 nm-1':'kJ/(mol*nm)'
        }
