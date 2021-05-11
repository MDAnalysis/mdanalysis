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

The `H5MD`_ trajectory file format is based upon the general, high performance
`HDF5`_ file format.
HDF5 files are self documenting and can be accessed with the `h5py`_ library.
HDF5 can make use of parallel file system features through the MPI-IO
interface of the HDF5 library to improve parallel reads and writes.

The HDF5 library and `h5py`_ must be installed; otherwise, H5MD files
cannot be read by MDAnalysis. If `h5py`_ is not installed, a
:exc:`RuntimeError` is raised.

Units
-----

H5MD files are very flexible and can store data in a wide range of physical
units. The :class:`H5MDReader` will attempt to match the units in order to
convert all data to the standard MDAnalysis units (see
:mod:`MDAnalysis.units`).

Units are read from the attributes of the position, velocity, force, and time
datasets provided by the H5MD file. The unit string is translated from `H5MD
notation`_ to `MDAnalysis notation`_. If MDAnalysis does not recognize the unit
(likely because that unit string is not defined in :mod:`MDAnalysis.units`)
provided, a :exc:`RuntimeError` is raised.  If no units are provided,
MDAnalysis stores a value of ``None`` for each unit.  If the H5MD file does not
contain units and ``convert_units=True``, MDAnalysis will raise a
:exc`ValueError`. To load a universe from an H5MD file with no units, set
``convert_units=False``.


Example: Loading an H5MD simulation
-----------------------------------

To load an H5MD simulation from an H5MD trajectory data file (using the
:class:`~MDAnalysis.coordinates.H5MD.H5MDReader`), pass the topology
and trajectory files to :class:`~MDAnalysis.core.universe.Universe`::

    import MDAnalysis as mda
    u = mda.Universe("topology.tpr", "trajectory.h5md")

It is also possible to pass an open :class:`h5py.File` file stream
into the reader::

    import MDAnalysis as mda
    with h5py.File("trajectory.h5md", 'r') as f:
         u = mda.Universe("topology.tpr", f)

.. Note:: Directly using a `h5py.File` does not work yet.
   See issue `#2884 <https://github.com/MDAnalysis/mdanalysis/issues/2884>`_.

Example: Opening an H5MD file in parallel
-----------------------------------------

The parallel features of HDF5 can be accessed through h5py
(see `parallel h5py docs`_ for more detail) by using the `mpi4py`_ Python
package with a Parallel build of HDF5. To load a an H5MD simulation with
parallel HDF5, pass `driver` and `comm` arguments to
:class:`~MDAnalysis.core.universe.Universe`::

    import MDAnalysis as mda
    from mpi4py import MPI
    u = mda.Universe("topology.tpr", "trajectory.h5md",
                     driver="mpio", comm=MPI.COMM_WORLD)

.. Note::
   :mod:`h5py` must be built with parallel features enabled on top of a parallel
   HDF5 build, and HDF5 and :mod:`mpi4py` must be built with a working MPI
   implementation. See instructions below.

Building parallel h5py and HDF5 on Linux
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Building a working parallel HDF5/h5py/mpi4py environment can be
challenging and is often specific to your local computing resources,
e.g., the supercomputer that you're running on typically already has
its preferred MPI installation. As a starting point we provide
instructions that worked in a specific, fairly generic environment.

These instructions successfully built parallel HDF5/h5py
with *OpenMPI 4.0.4*, *HDF5 1.10.6*, *h5py 2.9.0*, and *mpi4py 3.0.3*
on *Ubuntu 16.0.6*. You may have to play around with different combinations of
versions of h5py/HDF5 to get a working parallel build.

    1. `Build MPI from sources`_
    2. `Build HDF5 from sources`_ with parallel settings enabled:

       .. code-block:: bash

          ./configure --enable-parallel --enable-shared
          make
          make install

    3. `Install mpi4py`_, making sure to point `mpicc` to where you've
       installed your MPI implemenation:

       .. code-block:: bash

          env MPICC=/path/to/mpicc pip install mpi4py

    4. `Build h5py from sources`_, making sure to enable mpi and to point
       to your parallel build of HDF5:

       .. code-block:: bash

          export HDF5_PATH=path-to-parallel-hdf5
          python setup.py clean --all
          python setup.py configure -r --hdf5-version=X.Y.Z --mpi --hdf5=$HDF5_PATH
          export gcc=gcc
          CC=mpicc HDF5_DIR=$HDF5_PATH python setup.py build
          python setup.py install

If you have questions or want to share how you managed to build
parallel hdf5/h5py/mpi4py please let everyone know on the
`MDAnalysis forums`_.

.. _`H5MD`: https://nongnu.org/h5md/index.html
.. _`HDF5`: https://www.hdfgroup.org/solutions/hdf5/
.. _`H5PY`: http://docs.h5py.org/
.. _`parallel h5py docs`: https://docs.h5py.org/en/stable/mpi.html
.. _`mpi4py`: https://mpi4py.readthedocs.io/en/stable/index.html
.. _`Build MPI from sources`: https://mpi4py.readthedocs.io/en/stable/appendix.html#building-mpi-from-sources
.. _`Build HDF5 from sources`: https://support.hdfgroup.org/ftp/HDF5/current/src/unpacked/release_docs/INSTALL_parallel
.. _`Install mpi4py`: https://mpi4py.readthedocs.io/en/stable/install.html#requirements
.. _`Build h5py from sources`: https://docs.h5py.org/en/stable/mpi.html#building-against-parallel-hdf5
.. _`H5MD notation`: https://nongnu.org/h5md/modules/units.html
.. _`MDAnalysis notation`: https://userguide.mdanalysis.org/units.html
.. _`MDAnalysis forums`: https://www.mdanalysis.org/#participating


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
      or if the `velocities` = ``True`` keyword has been supplied.

   .. attribute:: forces

      forces of the atoms as a :class:`numpy.ndarray` of shape
      `(n_atoms, 3)`; only available if the trajectory contains forces
      or if the `forces` = ``True`` keyword has been supplied.


.. autoclass:: H5MDReader
   :members:

   .. automethod:: H5MDReader._reopen

.. autoclass:: H5PYPicklable
   :members:

"""

import numpy as np
import MDAnalysis as mda
from . import base, core
from ..exceptions import NoDataError
try:
    import h5py
except ImportError:
    HAS_H5PY = False

    # Allow building documentation even if h5py is not installed
    import types

    class MockH5pyFile:
        pass
    h5py = types.ModuleType("h5py")
    h5py.File = MockH5pyFile

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
        if self._unitcell is not None:
            return core.triclinic_box(*self._unitcell)

    @dimensions.setter
    def dimensions(self, box):
        self._unitcell[:] = core.triclinic_vectors(box)


class H5MDReader(base.ReaderBase):
    r"""Reader for the H5MD format.

    See `h5md documentation <https://nongnu.org/h5md/h5md.html>`_
    for a detailed overview of the H5MD file format.

    The reader attempts to convert units in the trajectory file to
    the standard MDAnalysis units (:mod:`MDAnalysis.units`) if
    `convert_units` is set to ``True``.

    Additional data in the *observables* group of the H5MD file are
    loaded into the :attr:`Timestep.data` dictionary.

    Only 3D-periodic boxes or no periodicity are supported; for no
    periodicity, :attr:`Timestep.dimensions` will return ``None``.

    Although H5MD can store varying numbers of particles per time step
    as produced by, e.g., GCMC simulations, MDAnalysis can currently
    only process a fixed number of particles per step. If the number
    of particles changes a :exc:`ValueError` is raised.

    The :class:`H5MDReader` reads .h5md files with the following
    HDF5 hierarchy:

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
                            +-- unit <str>
                \-- (position)
                    \-- [step] <int>, gives frame
                    \-- [time] <float>, gives time
                        +-- unit <str>
                    \-- [value] <float>, gives numpy arrary of positions
                                         with shape (n_atoms, 3)
                        +-- unit <str>
                \-- (velocity)
                    \-- [step] <int>, gives frame
                    \-- [time] <float>, gives time
                        +-- unit <str>
                    \-- [value] <float>, gives numpy arrary of velocities
                                         with shape (n_atoms, 3)
                        +-- unit <str>
                \-- (force)
                    \-- [step] <int>, gives frame
                    \-- [time] <float>, gives time
                        +-- unit <str>
                    \-- [value] <float>, gives numpy arrary of forces
                                         with shape (n_atoms, 3)
                        +-- unit <str>
         \-- (observables)
            \-- (lambda)
                \-- [step] <int>, gives frame
                \-- [time] <float>, gives time
                \-- [value] <float>
            \-- (step)
                \-- [step] <int>, gives frame
                \-- [time] <float>, gives time
                \-- [value] <int>, gives integration step


    .. note::
        The reader does not currently read mass or charge data.

    .. note::
        If the `driver` and `comm` arguments were used to open the
        hdf5 file (specifically, ``driver="mpio"``) then the :meth:`_reopen`
        method does *not* close and open the file like most readers because
        the information about the MPI communicator would be lost; instead
        it rewinds the trajectory back to the first timestep.


    .. versionadded:: 2.0.0

    """

    format = 'H5MD'
    # units is defined as instance-level variable and set from the
    # H5MD file in __init__() below

    # This dictionary is used to translate H5MD units to MDAnalysis units.
    # (https://nongnu.org/h5md/modules/units.html)
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

    def __init__(self, filename,
                 convert_units=True,
                 driver=None,
                 comm=None,
                 **kwargs):
        """
        Parameters
        ----------
        filename : str or :class:`h5py.File`
            trajectory filename or open h5py file
        convert_units : bool (optional)
            convert units to MDAnalysis units
        driver : str (optional)
            H5PY file driver used to open H5MD file
        comm : :class:`MPI.Comm` (optional)
            MPI communicator used to open H5MD file
            Must be passed with `'mpio'` file driver
        **kwargs : dict
            General reader arguments.

        Raises
        ------
        RuntimeError
            when `H5PY`_ is not installed
        RuntimeError
            when a unit is not recognized by MDAnalysis
        ValueError
            when ``n_atoms`` changes values between timesteps
        ValueError
            when ``convert_units=True`` but the H5MD file contains no units
        ValueError
            when dimension of unitcell is not 3
        ValueError
            when an MPI communicator object is passed to the reader
            but ``driver != 'mpio'``
        NoDataError
            when the H5MD file has no 'position', 'velocity', or
            'force' group

        """
        if not HAS_H5PY:
            raise RuntimeError("Please install h5py")
        super(H5MDReader, self).__init__(filename, **kwargs)
        self.filename = filename
        self.convert_units = convert_units

        # if comm is provided, driver must be 'mpio' and file will be
        # opened with parallel h5py/hdf5 enabled
        self._driver = driver
        self._comm = comm
        if (self._comm is not None) and (self._driver != 'mpio'):
            raise ValueError("If MPI communicator object is used to open"
                             " h5md file, ``driver='mpio'`` must be passed.")

        self.open_trajectory()
        if self._particle_group['box'].attrs['dimension'] != 3:
            raise ValueError("MDAnalysis only supports 3-dimensional"
                             " simulation boxes")

        # _has dictionary used for checking whether h5md file has
        # 'position', 'velocity', or 'force' groups in the file
        self._has = {name: name in self._particle_group for
                     name in ('position', 'velocity', 'force')}

        # Gets n_atoms from first available group
        for name, value in self._has.items():
            if value:
                self.n_atoms = self._particle_group[name]['value'].shape[1]
                break
        else:
            raise NoDataError("Provide at least a position, velocity"
                              " or force group in the h5md file.")

        self.ts = self._Timestep(self.n_atoms,
                                 positions=self.has_positions,
                                 velocities=self.has_velocities,
                                 forces=self.has_forces,
                                 **self._ts_kwargs)

        self.units = {'time': None,
                      'length': None,
                      'velocity': None,
                      'force': None}
        self._set_translated_units()  # fills units dictionary
        self._read_next_timestep()

    def _set_translated_units(self):
        """converts units from H5MD to MDAnalysis notation
        and fills units dictionary"""

        # need this dictionary to associate 'position': 'length'
        _group_unit_dict = {'time': 'time',
                            'position': 'length',
                            'velocity': 'velocity',
                            'force': 'force'
                            }

        for group, unit in _group_unit_dict.items():
            self._translate_h5md_units(group, unit)
            self._check_units(group, unit)

    def _translate_h5md_units(self, group, unit):
        """stores the translated unit string into the units dictionary"""

        errmsg = "{} unit '{}' is not recognized by H5MDReader. Please raise"
        " an issue in https://github.com/MDAnalysis/mdanalysis/issues"

        # doing time unit separately because time has to fish for
        # first available parent group - either position, velocity, or force
        if unit == 'time':
            for name, value in self._has.items():
                if value:
                    if 'unit' in self._particle_group[name]['time'].attrs:
                        try:
                            self.units['time'] = self._unit_translation[
                                'time'][self._particle_group[name][
                                    'time'].attrs['unit']]
                            break
                        except KeyError:
                            raise RuntimeError(errmsg.format(
                                               unit, self._particle_group[
                                               name]['time'].attrs['unit'])
                                               ) from None

        else:
            if self._has[group]:
                if 'unit' in self._particle_group[group]['value'].attrs:
                    try:
                        self.units[unit] = self._unit_translation[unit][
                            self._particle_group[group]['value'].attrs['unit']]
                    except KeyError:
                        raise RuntimeError(errmsg.format(
                                           unit, self._particle_group[group][
                                           'value'].attrs['unit'])
                                           ) from None

            # if position group is not provided, can still get 'length' unit
            # from unitcell box
            if (not self._has['position']) and ('edges' in self._particle_group['box']):
                if 'unit' in self._particle_group['box/edges/value'].attrs:
                    try:
                        self.units['length'] = self._unit_translation[
                                               'length'][self._particle_group[
                                               'box/edges/value'
                                               ].attrs['unit']]
                    except KeyError:
                        raise RuntimeError(errmsg.format(unit,
                                           self._particle_group[
                                           'box/edges/value'].attrs[
                                           'unit'])) from None

    def _check_units(self, group, unit):
        """Raises error if no units are provided from H5MD file
        and convert_units=True"""

        if not self.convert_units:
            return

        errmsg = "H5MD file must have readable units if ``convert_units`` is"
        " set to ``True``. MDAnalysis sets ``convert_units=True`` by default."
        " Set ``convert_units=False`` to load Universe without units."

        if unit == 'time':
            if self.units['time'] is None:
                raise ValueError(errmsg)

        else:
            if self._has[group]:
                if self.units[unit] is None:
                    raise ValueError(errmsg)

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
            self._driver = self._file.driver
        else:
            if self._comm is not None:
                # can only pass comm argument to h5py.File if driver='mpio'
                assert self._driver == 'mpio'
                self._file = H5PYPicklable(name=self.filename,  # pragma: no cover
                                           mode='r',
                                           driver=self._driver,
                                           comm=self._comm)
            elif self._driver is not None:
                self._file = H5PYPicklable(name=self.filename, mode='r',
                                           driver=self._driver)
            else:
                self._file = H5PYPicklable(name=self.filename, mode='r')
        # pulls first key out of 'particles'
        # allows for arbitrary name of group1 in 'particles'
        self._particle_group = self._file['particles'][
            list(self._file['particles'])[0]]

    @property
    def n_frames(self):
        """number of frames in trajectory"""
        for name, value in self._has.items():
            if value:
                return self._particle_group[name]['value'].shape[0]

    def _read_frame(self, frame):
        """reads data from h5md file and copies to current timestep"""
        try:
            for name, value in self._has.items():
                if value:
                    _ = self._particle_group[name]['step'][frame]
                    break
            else:
                raise NoDataError("Provide at least a position, velocity"
                                  " or force group in the h5md file.")
        except (ValueError, IndexError):
            raise IOError from None

        self._frame = frame
        ts = self.ts
        particle_group = self._particle_group
        ts.frame = frame

        # fills data dictionary from 'observables' group
        # Note: dt is not read into data as it is not decided whether
        # Timestep should have a dt attribute (see Issue #2825)
        self._copy_to_data()

        # Sets frame box dimensions
        # Note: H5MD files must contain 'box' group in each 'particles' group
        if 'edges' in particle_group['box'] and ts._unitcell is not None:
            ts._unitcell[:] = particle_group['box/edges/value'][frame, :]
        else:
            # sets ts.dimensions = None
            ts._unitcell = None

        # set the timestep positions, velocities, and forces with
        # current frame dataset
        if self._has['position']:
            ts.positions = self._get_frame_dataset('position')
        if self._has['velocity']:
            ts.velocities = self._get_frame_dataset('velocity')
        if self._has['force']:
            ts.forces = self._get_frame_dataset('force')

        if self.convert_units:
            self._convert_units()

        return ts

    def _copy_to_data(self):
        """assigns values to keys in data dictionary"""

        if 'observables' in self._file:
            for key in self._file['observables'].keys():
                self.ts.data[key] = self._file['observables'][key][
                    'value'][self._frame]

        # pulls 'time' out of first available parent group
        for name, value in self._has.items():
            if value:
                if 'time' in self._particle_group[name]:
                    self.ts.data['time'] = self._particle_group[name][
                        'time'][self._frame]
                    break

    def _get_frame_dataset(self, dataset):
        """retrieves dataset array at current frame"""

        frame_dataset = self._particle_group[
            dataset]['value'][self._frame, :]
        n_atoms_now = frame_dataset.shape[0]
        if n_atoms_now != self.n_atoms:
            raise ValueError("Frame {} has {} atoms but the initial frame"
                             " has {} atoms. MDAnalysis is unable to deal"
                             " with variable topology!"
                             "".format(self._frame,
                                       n_atoms_now,
                                       self.n_atoms))
        return frame_dataset

    def _convert_units(self):
        """converts time, position, velocity, and force values if they
        are not given in MDAnalysis standard units

        See https://userguide.mdanalysis.org/1.0.0/units.html
        """

        self.ts.time = self.convert_time_from_native(self.ts.time)

        if 'edges' in self._particle_group['box']:
            self.convert_pos_from_native(self.ts.dimensions[:3])

        if self._has['position']:
            self.convert_pos_from_native(self.ts.positions)

        if self._has['velocity']:
            self.convert_velocities_from_native(self.ts.velocities)

        if self._has['force']:
            self.convert_forces_from_native(self.ts.forces)

    def _read_next_timestep(self):
        """read next frame in trajectory"""
        return self._read_frame(self._frame + 1)

    def close(self):
        """close reader"""
        self._file.close()

    def _reopen(self):
        """reopen trajectory

        Note
        ----

        If the `driver` and `comm` arguments were used to open the
        hdf5 file (specifically, ``driver="mpio"``) then this method
        does *not* close and open the file like most readers because
        the information about the MPI communicator would be lost; instead
        it rewinds the trajectory back to the first timstep.

        """
        if self._driver == "mpio":  # pragma: no cover
            self._read_frame(-1)
            return

        self.close()
        self.open_trajectory()

    @property
    def has_positions(self):
        """``True`` if 'position' group is in trajectory."""
        return self._has['position']

    @has_positions.setter
    def has_positions(self, value: bool):
        self._has['position'] = value

    @property
    def has_velocities(self):
        """``True`` if 'velocity' group is in trajectory."""
        return self._has['velocity']

    @has_velocities.setter
    def has_velocities(self, value: bool):
        self._has['velocity'] = value

    @property
    def has_forces(self):
        """``True`` if 'force' group is in trajectory."""
        return self._has['force']

    @has_forces.setter
    def has_forces(self, value: bool):
        self._has['force'] = value

    def __getstate__(self):
        state = self.__dict__.copy()
        del state['_particle_group']
        return state

    def __setstate__(self, state):
        self.__dict__ = state
        self._particle_group = self._file['particles'][
                               list(self._file['particles'])[0]]
        self[self.ts.frame]


class H5PYPicklable(h5py.File):
    """H5PY file object (read-only) that can be pickled.

    This class provides a file-like object (as returned by
    :class:`h5py.File`) that,
    unlike standard Python file objects,
    can be pickled. Only read mode is supported.

    When the file is pickled, filename, mode, driver, and comm of
    :class:`h5py.File` in the file are saved. On unpickling, the file
    is opened by filename, mode, driver. This means that for a successful
    unpickle, the original file still has to be accessible with its filename.

    Parameters
    ----------
    filename : str or file-like
        a filename given a text or byte string.
    driver : str (optional)
        H5PY file driver used to open H5MD file

    Example
    -------
    ::

        f = H5PYPicklable('filename', 'r')
        print(f['particles/trajectory/position/value'][0])
        f.close()

    can also be used as context manager::

        with H5PYPicklable('filename', 'r'):
            print(f['particles/trajectory/position/value'][0])

    Note
    ----
    Pickling of an `h5py.File` opened with `driver="mpio"` and an MPI
    communicator is currently not supported

    See Also
    ---------
    :class:`MDAnalysis.lib.picklable_file_io.FileIOPicklable`
    :class:`MDAnalysis.lib.picklable_file_io.BufferIOPicklable`
    :class:`MDAnalysis.lib.picklable_file_io.TextIOPicklable`
    :class:`MDAnalysis.lib.picklable_file_io.GzipPicklable`
    :class:`MDAnalysis.lib.picklable_file_io.BZ2Picklable`


    .. versionadded:: 2.0.0
    """

    def __getstate__(self):
        driver = self.driver
        # Current issues: Need a way to retrieve MPI communicator object
        # from self and pickle MPI.Comm object. Parallel driver is excluded
        # from test because h5py calls for an MPI configuration when driver is
        # 'mpio', so this will need to be patched in the test function.
        if driver == 'mpio':  # pragma: no cover
            raise TypeError("Parallel pickling of `h5py.File` with"  # pragma: no cover
                            " 'mpio' driver is currently not supported.")

        return {'name': self.filename,
                'mode': self.mode,
                'driver': driver}

    def __setstate__(self, state):
        self.__init__(name=state['name'],
                      mode=state['mode'],
                      driver=state['driver'])

    def __getnewargs__(self):
        """Override the h5py getnewargs to skip its error message"""
        return ()
