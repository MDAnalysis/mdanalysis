# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
"""
H5MD trajectories --- :mod:`MDAnalysis.coordinates.H5MD`
========================================================

Classes
-------

.. autoclass:: Timestep
   :members:

   .. attribute:: _pos

      coordinates of the atoms as a :class:`numpy.ndarray` of shape `(n_atoms, 3)`

   .. attribute:: _velocities

      velocities of the atoms as a :class:`numpy.ndarray` of shape `(n_atoms, 3)`;
      only available if the trajectory contains velocities or if the
      *velocities* = ``True`` keyword has been supplied.

   .. attribute:: _forces

      forces of the atoms as a :class:`numpy.ndarray` of shape `(n_atoms, 3)`;
      only available if the trajectory contains forces or if the
      *forces* = ``True`` keyword has been supplied.


.. autoclass:: H5MDReader
   :members:


"""

import numpy as np
import MDAnalysis as mda
from MDAnalysis.coordinates import base, core
try:
    import h5py
except ImportError:
    HAS_H5PY = False
else:
    HAS_H5PY = True

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
        If `edges <https://nongnu.org/h5md/h5md.html#simulation-box>`_ is a matrix,
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
    See `h5md documentation <https://nongnu.org/h5md/h5md.html#>`_
    for detailed overview

    H5MD root
     \-- (h5md)
        +-- version
        \-- author
            +-- name: <str>, author's name
            +-- email: <str>, option email
        \-- creator
            +-- name: <str>, gives file that created h5md file
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
    units = {'time': None,
             'length': None,
             'velocity': None,
             'force': None}
    # units are added from h5md file
    _Timestep = Timestep

    def __init__(self, filename, convert_units=True, **kwargs):
        """
        Parameters
        ----------
        filename : str or h5py.File
            trajectory filename or open h5py file
        **kwargs : dict
            General reader arguments.
        """
        if not HAS_H5PY:
            raise RuntimeError("Please install h5py")
        super(H5MDReader, self).__init__(filename, **kwargs)
        self.filename = filename
        self.open_trajectory()
        self.n_atoms = self._particle_group['position/value'].shape[1]

        self.has_positions = 'position' in self._particle_group
        #if not self.has_positions:
            #raise ValueError("'position' group must be in 'particles' group")
        self.has_velocities = 'velocity' in self._particle_group
        self.has_forces = 'force' in self._particle_group

        self.ts = self._Timestep(self.n_atoms,
                                 velocities=self.has_velocities,
                                 forces=self.has_forces,
                                 **self._ts_kwargs)
        self._read_next_timestep()


    @staticmethod
    def _format_hint(thing):
        """Can this Reader read *thing*"""
        # nb, filename strings can still get passed through if
        # format='H5MD' is used
        return HAS_H5PY and isinstance(thing, h5py.File)


    def open_trajectory(self) :
        """opens the trajectory file using h5py library"""
        self._frame = -1
        if isinstance(self.filename, h5py.File):
            self._file = self.filename
        else:
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
        if 'data' in self._particle_group:
            if 'time' in self._particle_group['data']:
                ts.data['time'] = particle_group['position/time'][frame]
            if 'step' in self._particle_group['data']:
                ts.data['step'] = particle_group['data/step/value'][frame]
            if 'lambda' in self._particle_group['data']:
                ts.data['lambda'] = particle_group['data/lambda/value'][frame]
            if 'dt' in self._particle_group['data']:
                ts.data['dt'] = particle_group['data/dt/value'][frame]

        # set frame box dimensions
        # set triclinic box vectors
        ts._unitcell[:] = particle_group['box/edges/value'][frame, :]

        # set particle positions
        frame_positions = particle_group['position/value'][frame, :]
        n_atoms_now = frame_positions.shape[0]
        if n_atoms_now != self.n_atoms :
            raise ValueError("Frame %d has %d atoms but the initial frame has %d"
                             " atoms. MDAnalysis is unable to deal with variable"
                             " topology!"%(frame, n_atoms_now, self.n_atoms))
        else :
            ts.positions = frame_positions

        # set particle velocities
        if self.has_velocities:
            frame_velocities = particle_group['velocity/value'][frame, :]
            n_atoms_now = frame_velocities.shape[0]
            if n_atoms_now != self.n_atoms :
                raise ValueError("Frame %d has %d atoms but the initial frame has %d"
                                 " atoms. MDAnalysis is unable to deal with variable"
                                 " topology!"%(frame, n_atoms_now, self.n_atoms))
            else :
                ts.velocities = frame_velocities

        # set particle forces
        if self.has_forces:
            frame_forces = particle_group['force/value'][frame, :]
            n_atoms_now = frame_forces.shape[0]
            if n_atoms_now != self.n_atoms :
                raise ValueError("Frame %d has %d atoms but the initial frame has %d"
                                 " atoms. MDAnalysis is unable to deal with variable"
                                 " topology!"%(frame, n_atoms_now, self.n_atoms))
            else :
                ts.forces = frame_forces

        # unit conversion
        self._translate_h5md_units()  # fills units dictionary
        if self.convert_units:
            #self.ts.time = self.convert_time_to_native(self.ts.time)
            self.convert_pos_from_native(self.ts.dimensions[:3])
            self.convert_pos_from_native(self.ts.positions)
            if self.has_velocities:
                self.convert_velocities_from_native(self.ts.velocities)
            if self.has_forces:
                self.convert_forces_from_native(self.ts.forces)

        return ts


    def _read_next_timestep(self) :
        """read next frame in trajectory"""
        return self._read_frame(self._frame + 1)


    def _translate_h5md_units(self):
        """converts units from H5MD to MDAnalysis notation
        and fills units dictionary"""

        velocity_dict = {
        'nm ps-1': 'nm/ps',
        'Angstrom ps-1': 'Angstrom/ps'
        }
        force_dict = {
        'kJ mol-1 nm-1': 'kJ/(mol*nm)',
        'kJ mol-1 Angstrom-1': 'kJ/(mol*Angstrom)'
        }

        if 'units' in self._particle_group['position/time'].attrs:
            try:
                self.units['time'] = self._particle_group['position/time'].attrs['units']
            except KeyError:
                print("Unit %d not recognized by H5MDReader. Please "
                      "raise an issue on `MDAnalysis github <https://github.com/MDAnalysis/mdanalysis>`_"
                      ""%(self._particle_group['position/time'].attrs['units']))

        if 'units' in self._particle_group['position'].attrs:
            try:
                self.units['length'] = self._particle_group['position'].attrs['units']
            except KeyError:
                print("Unit %d not recognized by H5MDReader. Please "
                      "raise an issue on `MDAnalysis github <https://github.com/MDAnalysis/mdanalysis>`_"
                      ""%(self._particle_group['position'].attrs['units']))

        # looks for unit in h5md file and passes it to velocity_dict for translation
        if self.has_velocities:
            if 'units' in self._particle_group['velocity'].attrs:
                try:
                    self.units['velocity'] = velocity_dict[self._particle_group['velocity'].attrs['units']]
                except KeyError:
                    print("Unit %d not recognized by H5MDReader. Please "
                          "raise an issue on `MDAnalysis github <https://github.com/MDAnalysis/mdanalysis>`_"
                          ""%(self._particle_group['velocity'].attrs['units']))

        # looks for unit in h5md file and passes it to force_dict for translation
        if self.has_forces:
            if 'units' in self._particle_group['force'].attrs:
                try:
                    self.units['force'] = force_dict[self._particle_group['force'].attrs['units']]
                except KeyError:
                    print("Unit %d not recognized by H5MDReader. Please "
                          "raise an issue on `MDAnalysis github <https://github.com/MDAnalysis/mdanalysis>`_"
                          ""%(self._particle_group['force'].attrs['units']))
