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

The `H5MD`_ file format is based upon `HDF5`_, which makes use of parallel file
system features through the MPI-IO interface of the HDF5 library.

The reader currently uses the `H5PY`_ library to access data from an H5MD file.

.. _`H5MD`: https://nongnu.org/h5md/index.html
.. _`HDF5`: https://www.hdfgroup.org/solutions/hdf5/
.. _`H5PY`: http://docs.h5py.org/en/stable/


Classes
-------

.. autoclass:: Timestep
   :members:

   .. attribute:: positions

      coordinates of the atoms as a :class:`numpy.ndarray` of shape
      `(n_atoms, 3)`

   .. attribute:: velocities

      velocities of the atoms as a :class:`numpy.ndarray` of shape
      `(n_atoms, 3)`; only available if the trajectory contains velocities
       or if the *velocities* = ``True`` keyword has been supplied.

   .. attribute:: forces

      forces of the atoms as a :class:`numpy.ndarray` of shape
      `(n_atoms, 3)`; only available if the trajectory contains forces
      or if the *forces* = ``True`` keyword has been supplied.


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

        lengths *A*, *B*, *C* are in the MDAnalysis length unit (Å), and
        angles are in degrees.

        Setting dimensions will populate the underlying native format
        description (triclinic box vectors). If `edges
        <https://nongnu.org/h5md/h5md.html#simulation-box>`_ is a matrix,
        the box is of triclinic shape with the edge vectors given by
        the rows of the matrix.
        """
        return core.triclinic_box(*self._unitcell)

    @dimensions.setter
    def dimensions(self, box):
        self._unitcell[:] = core.triclinic_vectors(box)


class H5MDReader(base.ReaderBase):
    """Reader for the H5MD format.
    (See `h5md documentation <https://nongnu.org/h5md/h5md.html#>`_
    for a detailed overview)

    Currently reads .h5md files with the following HDF5 hierarchy:

    .. code-block:: text

        Notation:
        (name) is an HDF5 group that the reader recognizes
        {name} is an HDF5 group with arbitrary name
        [variable] is an HDF5 dataset
        <dtype> is dataset datatype
        +-- is an attribute of a group or dataset

        H5MD root
         \-- (h5md)
            +-- version <int>
            \-- author
                +-- name <str>, author's name
                +-- email <str>, optional email address
            \-- creator
                +-- name <str>, file that created .h5md file
                +-- version
         \-- (particles)
            \-- {group1}
                \-- (box)
                    +-- dimension : <int>, number of spatial dimensions
                    +-- boundary : <str>, boundary conditions of unit cell
                    \-- (edges)
                        \-- [step] <int>, gives frame
                        \-- [value] <float>, gives box dimensions
                \-- (position)
                    \-- [step] <int>, gives frame
                    \-- [time] <float>, gives time
                        +-- units <str>
                    \-- [value] <float>, gives numpy arrary of positions
                                         with shape (n_atoms, 3)
                        +-- units <str>
                \-- (velocity)
                    \-- [step] <int>, gives frame
                    \-- [time] <float>, gives time
                        +-- units <str>
                    \-- [value] <float>, gives numpy arrary of velocities
                                         with shape (n_atoms, 3)
                        +-- units <str>
                \-- (force)
                    \-- [step] <int>, gives frame
                    \-- [time] <float>, gives time
                        +-- units <str>
                    \-- [value] <float>, gives numpy arrary of forces
                                         with shape (n_atoms, 3)
                        +-- units <str>
                \-- (data)
                    \-- (dt)
                        \-- [step] <int>, gives frame number
                        \-- [value] <float>, gives dt
                    \-- (lambda)
                        \-- [step] <int>, gives frame
                        \-- [value] <float>, gives lambda
                    \-- (step)
                        \-- [step] <int>, gives frame
                        \-- [value] <int>, gives step


    .. note::
        The reader does not currently read mass data.

    """

    format = 'H5MD'
    # units are added from h5md file
    units = {'time': None,
             'length': None,
             'velocity': None,
             'force': None}
    # Translate H5MD units (https://nongnu.org/h5md/modules/units.html)
    # to MDAnalysis units.
    _unit_translation = {
        'time': {
            'ps': 'ps',
            'fs': 'fs',
            'ns': 'ns',
            'second': 'second',
            'sec': 'sec',
            's': 's',
            'AKMA': 'AKMA',
            },
        'length': {
            'Angstrom': 'Angstrom',
            'angstrom': 'Angstrom',
            'A': 'Angstrom',
            'nm': 'nm',
            'pm': 'pm',
            'fm': 'fm',
            },
        'velocity': {
            'Angstrom ps-1': 'Angstrom/ps',
            'A ps-1': 'Angstrom/ps',
            'Angstrom fs-1': 'Angstrom/fs',
            'A fs-1': 'Angstrom/fs',
            'Angstrom AKMA-1': 'Angstrom/AKMA',
            'A AKMA-1': 'Angstrom/AKMA',
            'nm ps-1': 'nm/ps',
            'nm ns-1': 'nm/ns',
            'pm ps-1': 'pm/ps',
            'm s-1': 'm/s'
            },
        'force':  {
            'kJ mol-1 Angstrom-1': 'kJ/(mol*Angstrom)',
            'kJ mol-1 nm-1': 'kJ/(mol*nm)',
            'Newton': 'Newton',
            'N': 'N',
            'J m-1': 'J/m',
            'kcal mol-1 Angstrom-1': 'kcal/(mol*Angstrom)',
            'kcal mol-1 A-1': 'kcal/(mol*Angstrom)'
            }
    }
    _Timestep = Timestep

    def __init__(self, filename, convert_units=True, **kwargs):
        """
        Parameters
        ----------
        filename : str or h5py.File
            trajectory filename or open h5py file
        convert_units : bool (optional)
            convert units to MDAnalysis units
        **kwargs : dict
            General reader arguments.
        """
        if not HAS_H5PY:
            raise RuntimeError("Please install h5py")
        super(H5MDReader, self).__init__(filename, **kwargs)
        self.filename = filename
        self.open_trajectory()
        self.has_positions = 'position' in self._particle_group
        if not self.has_positions:
            raise ValueError("'position' group must be in 'particles' group")
        self.n_atoms = self._particle_group['position/value'].shape[1]
        self.has_velocities = 'velocity' in self._particle_group
        self.has_forces = 'force' in self._particle_group
        self.ts = self._Timestep(self.n_atoms,
                                 velocities=self.has_velocities,
                                 forces=self.has_forces,
                                 **self._ts_kwargs)
        self._translate_h5md_units()  # fills units dictionary
        self._read_next_timestep()

    @staticmethod
    def _format_hint(thing):
        """Can this Reader read *thing*"""
        # nb, filename strings can still get passed through if
        # format='H5MD' is used
        return HAS_H5PY and isinstance(thing, h5py.File)

    def open_trajectory(self):
        """opens the trajectory file using h5py library"""
        self._frame = -1
        if isinstance(self.filename, h5py.File):
            self._file = self.filename
        else:
            self._file = h5py.File(self.filename, 'r')
        # pulls first key out of 'particles'
        # allows for arbitrary name of group1 in 'particles'
        self._particle_group = self._file['particles'][
                               list(self._file['particles'])[0]]

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
        """reads data from h5md file and copies to current timestep"""
        try:
            myframe = self._particle_group['position/step'][frame]
        except ValueError:
            raise IOError from None

        self._frame = frame
        ts = self.ts
        particle_group = self._particle_group
        ts.frame = frame

        # set data dictionary values
        if 'data' in particle_group:
            if 'time' in particle_group['data']:
                ts.data['time'] = particle_group['data/time/value'][frame]
            if 'step' in particle_group['data']:
                ts.data['step'] = particle_group['data/step/value'][frame]
            if 'lambda' in particle_group['data']:
                ts.data['lambda'] = particle_group['data/lambda/value'][frame]
            if 'dt' in particle_group['data']:
                ts.data['dt'] = particle_group['data/dt/value'][frame]

        # set frame box dimensions
        # set triclinic box vectors
        ts._unitcell[:] = particle_group['box/edges/value'][frame, :]

        # set particle positions
        frame_positions = particle_group['position/value'][frame, :]
        n_atoms_now = frame_positions.shape[0]
        if n_atoms_now != self.n_atoms:
            raise ValueError("Frame {} has {} atoms but the initial frame"
                             " has {} atoms. MDAnalysis is unable to deal"
                             " with variable topology!"
                             "".format(frame, n_atoms_now, self.n_atoms))
        else:
            ts.positions = frame_positions

        # set particle velocities
        if self.has_velocities:
            frame_velocities = particle_group['velocity/value'][frame, :]
            n_atoms_now = frame_velocities.shape[0]
            if n_atoms_now != self.n_atoms:
                raise ValueError("Frame {} has {} atoms but the initial frame"
                                 " has {} atoms. MDAnalysis is unable to deal"
                                 " with variable topology!"
                                 "".format(frame, n_atoms_now, self.n_atoms))
            else:
                ts.velocities = frame_velocities

        # set particle forces
        if self.has_forces:
            frame_forces = particle_group['force/value'][frame, :]
            n_atoms_now = frame_forces.shape[0]
            if n_atoms_now != self.n_atoms:
                raise ValueError("Frame {} has {} atoms but the initial frame"
                                 " has {} atoms. MDAnalysis is unable to deal"
                                 " with variable topology!"
                                 "".format(frame, n_atoms_now, self.n_atoms))
            else:
                ts.forces = frame_forces

        # unit conversion
        if self.convert_units:
            ts.time = self.convert_time_from_native(ts.time)
            self.convert_pos_from_native(ts.dimensions[:3])
            self.convert_pos_from_native(ts.positions)
            if self.has_velocities:
                self.convert_velocities_from_native(ts.velocities)
            if self.has_forces:
                self.convert_forces_from_native(ts.forces)

        return ts

    def _read_next_timestep(self):
        """read next frame in trajectory"""
        return self._read_frame(self._frame + 1)

    def _translate_h5md_units(self):
        """converts units from H5MD to MDAnalysis notation
        and fills units dictionary"""

        if 'units' in self._particle_group['position/time'].attrs:
            try:
                self.units['time'] = self._unit_translation['time'][
                                     self._particle_group
                                     ['position/time'].attrs['units']]
            except KeyError:
                raise RuntimeError("Time unit '{}' is not recognized by"
                                   " H5MDReader. Please raise an issue in"
                                   " https://github.com/MDAnalysis/"
                                   "mdanalysis/issues".format(
                                    self._particle_group['position/time'].
                                    attrs['units'])) from None

        if 'units' in self._particle_group['position'].attrs:
            try:
                self.units['length'] = self._unit_translation['length'][
                                       self._particle_group
                                       ['position'].attrs['units']]
            except KeyError:
                raise RuntimeError("Length unit '{}' is not recognized by"
                                   " H5MDReader. Please raise an issue in"
                                   " https://github.com/MDAnalysis/"
                                   "mdanalysis/issues".format(
                                    self._particle_group['position'].
                                    attrs['units'])) from None

        if self.has_velocities:
            if 'units' in self._particle_group['velocity'].attrs:
                try:
                    self.units['velocity'] = self._unit_translation[
                                             'velocity'][self._particle_group[
                                             'velocity'].attrs['units']]
                except KeyError:
                    raise RuntimeError("Velocity unit '{}' is not recognized"
                                       " by H5MDReader. Please raise an issue"
                                       " in https://github.com/MDAnalysis/"
                                       "mdanalysis/issues".format(
                                        self._particle_group['velocity'].
                                        attrs['units'])) from None

        if self.has_forces:
            if 'units' in self._particle_group['force'].attrs:
                try:
                    self.units['force'] = self._unit_translation['force'][
                                          self._particle_group
                                          ['force'].attrs['units']]
                except KeyError:
                    raise RuntimeError("Force unit '{}' is not recognized by"
                                       " H5MDReader. Please raise an issue in"
                                       " https://github.com/MDAnalysis/"
                                       "mdanalysis/issues".format(
                                        self._particle_group['force'].
                                        attrs['units'])) from None
