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
"""H5MD trajectories --- :mod:`MDAnalysis.coordinates.H5MD`
========================================================

The `H5MD`_ trajectory file format is based upon the general, high performance
`HDF5`_ file format.
HDF5 files are self documenting and can be accessed with the `h5py`_ library.
The HDF5 library (and `H5PY`_) must be installed; otherwise, H5MD files
cannot be read by MDAnalysis. If `H5PY`_ is not installed, a ``RuntimeError``
is raised.

HDF5 can make use of parallel file system features through the MPI-IO
interface of the HDF5 library to improve parallel reads and writes.


The `H5MD`_ file format is based upon `HDF5`_, which makes use of parallel
file system features through the MPI-IO interface of the HDF5 library.
The reader currently uses the `H5PY`_ library to access data from an H5MD file.

.. _`H5MD`: https://nongnu.org/h5md/index.html
.. _`HDF5`: https://www.hdfgroup.org/solutions/hdf5/
.. _`H5PY`: http://docs.h5py.org/
.. _`H5MD notation`: https://nongnu.org/h5md/modules/units.html
.. _`MDAnalysis notation`: https://userguide.mdanalysis.org/1.0.0/units.html


Classes
-------

.. autoclass:: H5MDWriter
   :members:

"""

import numpy as np
import MDAnalysis as mda
from MDAnalysis.coordinates import base, core
from MDAnalysis.exceptions import NoDataError
try:
    import h5py
except ImportError:
    HAS_H5PY = False
else:
    HAS_H5PY = True


class H5MDWriter(base.WriterBase):
    """Writer for the H5MD format.

    Classes
    -------

    .. autoclass:: Timestep
       :members:

       .. attribute:: positions

          coordinates of the atoms as a :class:`numpy.ndarray` of shape
          `(n_atoms, 3)`; only available if the trajectory contains positions
          or if the ``positions = True`` keyword has been supplied.

       .. attribute:: velocities

          velocities of the atoms as a :class:`numpy.ndarray` of shape
          `(n_atoms, 3)`; only available if the trajectory contains velocities
          or if the ``velocities = True`` keyword has been supplied.

       .. attribute:: forces

          forces of the atoms as a :class:`numpy.ndarray` of shape
          `(n_atoms, 3)`; only available if the trajectory contains forces
          or if the ``forces = True`` keyword has been supplied.


    .. autoclass:: H5MDReader
       :members:

    .. versionadded:: 2.0.0

    """

    format = 'H5MD'
    multiframe = True
    _unit_translation = {
        'time': {
            'ps': 'ps',
            'fs': 'fs',
            'ns': 'ns',
            'second': 's',
            'sec': 's',
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
            'Angstrom/ps': 'Angstrom ps-1',
            'A/ps': 'Angstrom ps-1',
            'Angstrom/fs': 'Angstrom fs-1',
            'A/fs': 'Angstrom fs-1',
            'Angstrom/AKMA': 'Angstrom AKMA-1',
            'A/AKMA': 'Angstrom AKMA-1',
            'nm/ps': 'nm ps-1',
            'nm/ns': 'nm ns-1',
            'pm/ps': 'pm ps-1',
            'm/s': 'm s-1'
            },
        'force':  {
            'kJ/(mol*Angstrom)': 'kJ mol-1 Angstrom-1',
            'kJ/(mol*nm)': 'kJ mol-1 nm-1',
            'Newton': 'Newton',
            'N': 'N',
            'J/m': 'J m-1',
            'kcal/(mol*Angstrom)': 'kcal mol-1 Angstrom-1',
            'kcal/(mol*A)': 'kcal mol-1 Angstrom-1'
            }
    }

    def __init__(self,
                 filename,
                 n_atoms,
                 convert_units=True,
                 units=None,
                 **kwargs):
        """
        Parameters
        ----------
        filename : str or :class:`h5py.File`
            trajectory filename or open h5py file
        n_atoms : int
            number of atoms in trajectory
        convert_units : bool (optional)
            convert units to MDAnalysis units
        **kwargs : dict

        """
        self.filename = filename
        if n_atoms == 0:
            raise ValueError("H5MDWriter: no atoms in output trajectory")
        self.n_atoms = n_atoms
        self.convert_units = convert_units

        self.h5mdfile = None
        self.has_positions = kwargs.get('positions', False)
        self.has_velocities = kwargs.get('velocities', False)
        self.has_forces = kwargs.get('forces', False)

        # Pull out various keywords to store in 'h5md' group
        self.author = kwargs.pop('author', 'N/A')
        self.author_email = kwargs.pop('author_email', None)
        self.creator = kwargs.pop('creator', 'MDAnalysis')
        self.creator_version = kwargs.pop('creator_version', mda.__version__)

        self.units = units

    def _init_h5md(self, periodic=True):
        """Initializes H5MD trajectory.

        The `H5MD`_ file is opened using the `H5PY`_ library. """

        h5md_file = h5py.File(self.filename, 'w')
        h5md_group = h5md_file.create_group('h5md')
        h5md_group.attrs['version'] = np.array([1,1])
        author_group = h5md_file.create_group('h5md/author')
        author_group.attrs['name'] = self.author
        if self.author_email is not None:
            author_group.attrs['email'] = self.author_email
        creator_group = h5md_file.create_group('h5md/creator')
        creator_group.attrs['name'] = self.creator
        creator_group.attrs['version'] = self.creator_version

        self.h5mdfile = h5md_file
        self._trajectory_group = h5mdfile.require_group('particles').require_group('trajectory')

    def is_periodic(self, ts):
        """Test if timestep ``ts`` contains a periodic box.

        Parameters
        ----------
        ts : :class:`Timestep`
             :class:`Timestep` instance containing coordinates to
             be written to trajectory file

        Returns
        -------
        bool
            Return ``True`` if `ts` contains a valid simulation box
        """
        return np.all(ts.dimensions > 0)

    def _write_next_frame(self, ag):
        """Write information associated with ``ag`` at current frame
        into trajectory

        Parameters
        ----------
        ag : AtomGroup or Universe

        .. versionadded:: 2.0.0

        """
        try:
            # Atomgroup?
            ts = ag.ts
        except AttributeError:
            try:
                # Universe?
                ts = ag.trajectory.ts
            except AttributeError:
                errmsg = "Input obj is neither an AtomGroup or Universe"
                raise TypeError(errmsg) from None

        if ts.n_atoms != self.n_atoms:
            raise IOError("H5MDWriter: Timestep does not have"
                          " the correct number of atoms")

        # Opens H5MD file to write with H5PY library and initializes
        # 'h5md' group and checks if unitcell is periodic
        if self.h5mdfile is None:
            self._init_h5md(periodic = self.is_periodic(ts))

        trajectory = self._trajectory_group  # conventient assignment

        # Checks to see if writer given boolean arguments for position,
        # velocity, and force, then instantiates H5MD dataset with
        # _create_element() method. Also adds translated units to each group.
        # Box dimensions are created at first available boolean
        if self.has_positions:
            self.position = self._create_element(parent_group=trajectory,
                                                 name='position',
                                                 data=ts.positions,
                                                 dataset_type='time')
            trajectory['position'].attrs['units'] = self._unit_translation['length'][self.units['length']]
            self.box = self._create_box(dimension=3,
                                        boundary=['periodic',
                                                  'periodic',
                                                  'periodic'],
                                        data=ts.triclinic_dimensions,
                                        dataset_type='time',
                                        step_from=self.position)

        if self.has_velocities:
            if self.has_positions:
                self.velocity = self._create_element(parent_group=trajectory,
                                                     name='velocity',
                                                     data=ts.velocities,
                                                     step_from=self.position,
                                                     dataset_type='time')
                trajectory['velocity'].attrs['units'] = self._unit_translation['velocity'][self.units['velocity']]
            else:
                self.velocity = self._create_element(parent_group=trajectory,
                                                     name='velocity',
                                                     data=ts.velocities,
                                                     dataset_type='time')
                trajectory['velocity'].attrs['units'] = self._unit_translation['velocity'][self.units['velocity']]
                self.box = self._create_box(dimension=3,
                                            boundary=['periodic',
                                                      'periodic',
                                                      'periodic'],
                                            data=ts.triclinic_dimensions,
                                            step_from=self.velocity)

        if self.has_forces:
            if self.has_positions:
                self.force = self._create_element(parent_group=trajectory,
                                                  name='force',
                                                  data=ts.forces,
                                                  step_from=self.position,
                                                  dataset_type='time')
                trajectory['force'].attrs['units'] = self._unit_translation['force'][self.units['force']]
            elif self.has_velocities:
                self.force = self._create_element(parent_group=trajectory,
                                                  name='force',
                                                  data=ts.forces,
                                                  step_from=self.velocity,
                                                  dataset_type='time')
                trajectory['force'].attrs['units'] = self._unit_translation['force'][self.units['force']]
            else:
                self.force = self._create_element(parent_group=trajectory,
                                                  name='force',
                                                  data=ts.forces,
                                                  dataset_type='time')
                trajectory['force'].attrs['units'] = self._unit_translation['force'][self.units['force']]
                self.box = self._create_box(dimension=3,
                                            boundary=['periodic',
                                                      'periodic',
                                                      'periodic'],
                                            data=ts.triclinic_dimensions,
                                            step_from=self.force)


        return self._write_next_timestep(ts)

    def _write_next_timestep(self, ts):
        """Write coordinates and unitcell information to H5MD file.

        Do not call this method directly; instead use
        :meth:`write` because some essential setup is done
        there before writing the first frame.

        """

        if self.convert_units:
            self.convert_pos_to_native(ts.positions)
            self.convert_pos_to_native(ts.triclinic_dimensions)
            self.convert_velocities_to_native(ts.velocities)
            self.convert_forces_to_native(ts.forces)
            ts.time = self.convert_time_to_native(ts.time)

        if self.has_positions:
            self.position.append(value=ts.positions,
                                 step=ts.frame,
                                 time=ts.time)
        if self.has_velocities:
            self.velocity.append(value=ts.velocities,
                                 step=ts.frame,
                                 time=ts.time)
        if self.has_forces:
            self.force.append(value=ts.forces,
                              step=ts.frame,
                              time=time)
        self.box.edges.append(value=ts.triclinic_dimensions,
                              step=ts.frame,
                              time=ts.time)
