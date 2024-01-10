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
:exc:`ValueError`. To load a universe from an H5MD file with no units defined,
set ``convert_units=False``.

:class:`H5MDWriter` detects the native units of the parent trajectory and
writes the trajectory with those units, unless one of `timeunit`,
`lengthunit`, `velocityunit`, `forceunit` arugments are supplied. In
this case, the writer will write the corresponding dataset with the selected
unit only if it is recognized by `MDAnalysis units`_.


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


Example: Writing an H5MD file
-----------------------------

To write to an H5MD file from a trajectory loaded with MDAnalysis, do:

.. code-block:: python

    import MDAnalysis as mda
    u = mda.Universe("topology.tpr", "trajectory.h5md")
    with mda.Writer("output.h5md", n_atoms=u.trajectory.n_atoms) as W:
        for ts in u.trajectory:
            W.write(u)

To write an H5MD file with contiguous datasets, you must specify the
number of frames to be written and set ``chunks=False``:

.. code-block:: python

    with mda.Writer("output_contigous.h5md", n_atoms=u.trajectory.n_atoms,
                    n_frames=3, chunks=False) as W:
        for ts in u.trajectory[:3]:
            W.write(u)

The writer also supports writing directly from an :class:`~MDAnalysis.core.groups.AtomGroup`::

    ag = u.select_atoms("protein and name CA")
    ag.write("alpha_carbons.h5md", frames='all')


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
.. _`MDAnalysis units`: https://userguide.mdanalysis.org/units.html
.. _`MDAnalysis forums`: https://www.mdanalysis.org/#participating


Classes
-------

.. autoclass:: H5MDReader
   :members:
   :inherited-members:

   .. automethod:: H5MDReader._reopen

.. autoclass:: H5MDWriter
   :members:
   :inherited-members:

.. autoclass:: H5PYPicklable
   :members:

"""

import numpy as np
import MDAnalysis as mda
from . import base, core
from ..exceptions import NoDataError
from ..due import due, Doi
from MDAnalysis.lib.util import store_init_arguments
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


class H5MDReader(base.ReaderBase):
    r"""Reader for the H5MD format.

    See `h5md documentation <https://nongnu.org/h5md/h5md.html>`_
    for a detailed overview of the H5MD file format.

    The reader attempts to convert units in the trajectory file to
    the standard MDAnalysis units (:mod:`MDAnalysis.units`) if
    `convert_units` is set to ``True``.

    Additional data in the *observables* group of the H5MD file are
    loaded into the :attr:`Timestep.data
    <MDAnalysis.coordinates.timestep.Timestep.data>` dictionary.

    Only 3D-periodic boxes or no periodicity are supported; for no
    periodicity, :attr:`Timestep.dimensions
    <MDAnalysis.coordinates.timestep.Timestep.dimensions>` will return ``None``.

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
    .. versionchanged:: 2.1.0
       Adds :meth:`parse_n_atoms` method to obtain the number of atoms directly
       from the trajectory by evaluating the shape of the ``position``,
       ``velocity``, or ``force`` groups.
    .. versionchanged:: 2.5.0
       Add correct handling of simple cuboid boxes

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

    @due.dcite(Doi("10.25080/majora-1b6fd038-005"),
               description="MDAnalysis trajectory reader/writer of the H5MD"
               "format", path=__name__)
    @due.dcite(Doi("10.1016/j.cpc.2014.01.018"),
               description="Specifications of the H5MD standard",
               path=__name__, version='1.1')
    @store_init_arguments
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

        # Gets some info about what settings the datasets were created with
        # from first available group
        for name, value in self._has.items():
            if value:
                dset = self._particle_group[f'{name}/value']
                self.n_atoms = dset.shape[1]
                self.compression = dset.compression
                self.compression_opts = dset.compression_opts
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

    @staticmethod
    def parse_n_atoms(filename):
        with h5py.File(filename, 'r') as f:
            for group in f['particles/trajectory']:
                if group in ('position', 'velocity', 'force'):
                    n_atoms = f[f'particles/trajectory/{group}/value'].shape[1]
                    return n_atoms

            raise NoDataError("Could not construct minimal topology from the "
                              "H5MD trajectory file, as it did not contain a "
                              "'position', 'velocity', or 'force' group. "
                              "You must include a topology file.")

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
            else:
                self._file = H5PYPicklable(name=self.filename,
                                           mode='r',
                                           driver=self._driver)

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
        if "edges" in particle_group["box"]:
            edges = particle_group["box/edges/value"][frame, :]
            # A D-dimensional vector or a D × D matrix, depending on the
            # geometry of the box, of Float or Integer type. If edges is a
            # vector, it specifies the space diagonal of a cuboid-shaped box.
            # If edges is a matrix, the box is of triclinic shape with the edge
            # vectors given by the rows of the matrix.
            if edges.shape == (3,):
                ts.dimensions = [*edges, 90, 90, 90]
            else:
                ts.dimensions = core.triclinic_box(*edges)
        else:
            ts.dimensions = None

        # set the timestep positions, velocities, and forces with
        # current frame dataset
        if self._has['position']:
            self._read_dataset_into_ts('position', ts.positions)
        if self._has['velocity']:
            self._read_dataset_into_ts('velocity', ts.velocities)
        if self._has['force']:
            self._read_dataset_into_ts('force', ts.forces)

        if self.convert_units:
            self._convert_units()

        return ts

    def _copy_to_data(self):
        """assigns values to keys in data dictionary"""

        if 'observables' in self._file:
            for key in self._file['observables'].keys():
                self.ts.data[key] = self._file['observables'][key][
                    'value'][self._frame]

        # pulls 'time' and 'step' out of first available parent group
        for name, value in self._has.items():
            if value:
                if 'time' in self._particle_group[name]:
                    self.ts.time = self._particle_group[name][
                        'time'][self._frame]
                    break
        for name, value in self._has.items():
            if value:
                if 'step' in self._particle_group[name]:
                    self.ts.data['step'] = self._particle_group[name][
                        'step'][self._frame]
                    break

    def _read_dataset_into_ts(self, dataset, attribute):
        """reads position, velocity, or force dataset array at current frame
        into corresponding ts attribute"""

        n_atoms_now = self._particle_group[f'{dataset}/value'][
                                           self._frame].shape[0]
        if n_atoms_now != self.n_atoms:
            raise ValueError(f"Frame {self._frame} of the {dataset} dataset"
                             f" has {n_atoms_now} atoms but the initial frame"
                             " of either the postion, velocity, or force"
                             f" dataset had {self.n_atoms} atoms."
                             " MDAnalysis is unable to deal"
                             " with variable topology!")

        self._particle_group[f'{dataset}/value'].read_direct(
                             attribute, source_sel=np.s_[self._frame, :])

    def _convert_units(self):
        """converts time, position, velocity, and force values if they
        are not given in MDAnalysis standard units

        See https://userguide.mdanalysis.org/1.0.0/units.html
        """

        self.ts.time = self.convert_time_from_native(self.ts.time)

        if 'edges' in self._particle_group['box'] and self.ts.dimensions is not None:
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

    def Writer(self, filename, n_atoms=None, **kwargs):
        """Return writer for trajectory format

        Note
        ----
        The chunk shape of the input file will not be copied to the output
        file, as :class:`H5MDWriter` uses a chunk shape of ``(1, n_atoms, 3)``
        by default. To use a custom chunk shape, you must specify the
        `chunks` argument. If you would like to copy an existing chunk
        format from a dataset (positions, velocities, or forces), do
        the following::

            chunks = u.trajectory._particle_group['position/value'].chunks

        Note that the writer will set the same layout for all particle groups.

        See Also
        --------
        :class:`H5MDWriter`  Output class for the H5MD format


        .. versionadded:: 2.0.0

        """
        if n_atoms is None:
            n_atoms = self.n_atoms
        kwargs.setdefault('driver', self._driver)
        kwargs.setdefault('compression', self.compression)
        kwargs.setdefault('compression_opts', self.compression_opts)
        kwargs.setdefault('positions', self.has_positions)
        kwargs.setdefault('velocities', self.has_velocities)
        kwargs.setdefault('forces', self.has_forces)
        return H5MDWriter(filename, n_atoms, **kwargs)

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


class H5MDWriter(base.WriterBase):
    """Writer for `H5MD`_ format (version 1.1).

    H5MD trajectories are automatically recognised by the
    file extension ".h5md".

    All data from the input :class:`~MDAnalysis.coordinates.timestep.Timestep` is
    written by default. For detailed information on how :class:`H5MDWriter`
    handles units, compression, and chunking, see the Notes section below.

    Note
    ----
    Parellel writing with the use of a MPI communicator and the ``'mpio'``
    HDF5 driver is currently not supported.

    Note
    ----
    :exc:`NoDataError` is raised if no positions, velocities, or forces are
    found in the input trajectory. While the H5MD standard allows for this
    case, :class:`H5MDReader` cannot currently read files without at least
    one of these three groups.

    Note
    ----
    Writing H5MD files with fancy trajectory slicing where the Timestep
    does not increase monotonically such as ``u.trajectory[[2,1,0]]``
    or ``u.trajectory[[0,1,2,0,1,2]]`` raises a :exc:`ValueError` as this
    violates the rules of the step dataset in the H5MD standard.

    Parameters
    ----------
    filename : str or :class:`h5py.File`
        trajectory filename or open h5py file
    n_atoms : int
        number of atoms in trajectory
    n_frames : int (optional)
        number of frames to be written in trajectory
    driver : str (optional)
        H5PY file driver used to open H5MD file. See `H5PY drivers`_ for
        list of available drivers.
    convert_units : bool (optional)
        Convert units from MDAnalysis to desired units
    chunks : tuple (optional)
        Custom chunk layout to be applied to the position,
        velocity, and force datasets. By default, these datasets
        are chunked in ``(1, n_atoms, 3)`` blocks
    compression : str or int (optional)
        HDF5 dataset compression setting to be applied
        to position, velocity, and force datasets. Allowed
        settings are 'gzip', 'szip', 'lzf'. If an integer
        in range(10), this indicates gzip compression level.
        Otherwise, an integer indicates the number of a
        dynamically loaded compression filter.
    compression_opts : int or tup (optional)
        Compression settings.  This is an integer for gzip, 2-tuple for
        szip, etc. If specifying a dynamically loaded compression filter
        number, this must be a tuple of values. For gzip, 1 indicates
        the lowest level of compression and 9 indicates maximum compression.
    positions : bool (optional)
        Write positions into the trajectory [``True``]
    velocities : bool (optional)
        Write velocities into the trajectory [``True``]
    forces : bool (optional)
        Write forces into the trajectory [``True``]
    timeunit : str (optional)
        Option to convert values in the 'time' dataset to a custom unit,
        must be recognizable by MDAnalysis
    lengthunit : str (optional)
        Option to convert values in the 'position/value' dataset to a
        custom unit, must be recognizable by MDAnalysis
    velocityunit : str (optional)
        Option to convert values in the 'velocity/value' dataset to a
        custom unit, must be recognizable by MDAnalysis
    forceunit : str (optional)
        Option to convert values in the 'force/value' dataset to a
        custom unit, must be recognizable by MDAnalysis
    author : str (optional)
        Name of the author of the file
    author_email : str (optional)
        Email of the author of the file
    creator : str (optional)
        Software that wrote the file [``MDAnalysis``]
    creator_version : str (optional)
        Version of software that wrote the file
        [:attr:`MDAnalysis.__version__`]

    Raises
    ------
    RuntimeError
        when `H5PY`_ is not installed
    ValueError
        when `n_atoms` is 0
    ValueError
        when ``chunks=False`` but the user did not specify `n_frames`
    ValueError
        when `positions`, `velocities`, and `forces` are all
        set to ``False``
    TypeError
        when the input object is not a :class:`Universe` or
        :class:`AtomGroup`
    IOError
        when `n_atoms` of the :class:`Universe` or :class:`AtomGroup`
        being written does not match `n_atoms` passed as an argument
        to the writer
    ValueError
        when any of the optional `timeunit`, `lengthunit`,
        `velocityunit`, or `forceunit` keyword arguments are
        not recognized by MDAnalysis

    Notes
    -----

    By default, the writer will write all available data (positions,
    velocities, and forces) if detected in the input
    :class:`~MDAnalysis.coordinates.timestep.Timestep`. In addition, the settings
    for `compression` and `compression_opts` will be read from
    the first available group of positions, velocities, or forces and used as
    the default value. To write a file without any one of these datsets,
    set `positions`, `velocities`, or `forces` to ``False``.

    .. rubric:: Units

    The H5MD format is very flexible with regards to units, as there is no
    standard defined unit for the format. For this reason, :class:`H5MDWriter`
    does not enforce any units. The units of the written trajectory can be set
    explicitly with the keyword arguments `lengthunit`, `velocityunit`,
    and `forceunit`. If units are not explicitly specified, they are set to
    the native units of the trajectory that is the source of the coordinates.
    For example, if one converts a DCD trajectory, then positions are written
    in ångstrom and time in AKMA units. A GROMACS XTC will be written in nm and
    ps. The units are stored in the metadata of the H5MD file so when
    MDAnalysis loads the H5MD trajectory, the units will be automatically
    set correctly.

    .. rubric:: Compression

    HDF5 natively supports various compression modes. To write the trajectory
    with compressed datasets, set ``compression='gzip'``, ``compression='lzf'``
    , etc. See `H5PY compression options`_ for all supported modes of
    compression. An additional argument, `compression_opts`, can be used to
    fine tune the level of compression. For example, for GZIP compression,
    `compression_opts` can be set to 1 for minimum compression and 9 for
    maximum compression.

    .. rubric:: HDF5 Chunking

    HDF5 datasets can be *chunked*, meaning the dataset can be split into equal
    sized pieces and stored in different, noncontiguous places on disk.
    If HDF5 tries to read an element from a chunked dataset, the *entire*
    dataset must be read, therefore an ill-thought-out chunking scheme can
    drastically effect file I/O performance. In the case of all MDAnalysis
    writers, in general, the number of frames being written is not known
    apriori by the writer, therefore the HDF5 must be extendable. However, the
    allocation of diskspace is defined when the dataset is created, therefore
    extendable HDF5 datasets *must* be chunked so as to allow dynamic storage
    on disk of any incoming data to the writer. In such cases where chunking
    isn't explicity defined by the user, H5PY automatically selects a chunk
    shape via an algorithm that attempts to make mostly square chunks between
    1 KiB - 1 MiB, however this can lead to suboptimal I/O performance.
    :class:`H5MDWriter` uses a default chunk shape of ``(1, n_atoms, 3)`` so
    as to mimic the typical access pattern of a trajectory by MDAnalysis. In
    our tests ([Jakupovic2021]_), this chunk shape led to a speedup on the
    order of 10x versus H5PY's auto-chunked shape. Users can set a custom
    chunk shape with the `chunks` argument. Additionaly, the datasets in a
    file can be written with a contiguous layout by setting ``chunks=False``,
    however this must be accompanied by setting `n_frames` equal to the
    number of frames being written, as HDF5 must know how much space to
    allocate on disk when creating the dataset.

    .. _`H5PY compression options`: https://docs.h5py.org/en/stable/high/dataset.html#filter-pipeline
    .. _`H5PY drivers`: https://docs.h5py.org/en/stable/high/file.html#file-drivers


    .. versionadded:: 2.0.0

    """

    format = 'H5MD'
    multiframe = True
    #: These variables are not written from :attr:`Timestep.data`
    #: dictionary to the observables group in the H5MD file
    data_blacklist = ['step', 'time', 'dt']

    #: currently written version of the file format
    H5MD_VERSION = (1, 1)

    # This dictionary is used to translate MDAnalysis units to H5MD units.
    # (https://nongnu.org/h5md/modules/units.html)
    _unit_translation_dict = {
        'time': {
            'ps': 'ps',
            'fs': 'fs',
            'ns': 'ns',
            'second': 's',
            'sec': 's',
            's': 's',
            'AKMA': 'AKMA'},
        'length': {
            'Angstrom': 'Angstrom',
            'angstrom': 'Angstrom',
            'A': 'Angstrom',
            'nm': 'nm',
            'pm': 'pm',
            'fm': 'fm'},
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
            'm/s': 'm s-1'},
        'force':  {
            'kJ/(mol*Angstrom)': 'kJ mol-1 Angstrom-1',
            'kJ/(mol*nm)': 'kJ mol-1 nm-1',
            'Newton': 'Newton',
            'N': 'N',
            'J/m': 'J m-1',
            'kcal/(mol*Angstrom)': 'kcal mol-1 Angstrom-1',
            'kcal/(mol*A)': 'kcal mol-1 Angstrom-1'}}

    @due.dcite(Doi("10.25080/majora-1b6fd038-005"),
               description="MDAnalysis trajectory reader/writer of the H5MD"
               "format", path=__name__)
    @due.dcite(Doi("10.1016/j.cpc.2014.01.018"),
               description="Specifications of the H5MD standard",
               path=__name__, version='1.1')
    def __init__(self, filename, n_atoms, n_frames=None, driver=None,
                 convert_units=True, chunks=None, compression=None,
                 compression_opts=None, positions=True, velocities=True,
                 forces=True, timeunit=None, lengthunit=None,
                 velocityunit=None, forceunit=None, author='N/A',
                 author_email=None, creator='MDAnalysis',
                 creator_version=mda.__version__, **kwargs):

        if not HAS_H5PY:
            raise RuntimeError("H5MDWriter: Please install h5py")
        self.filename = filename
        if n_atoms == 0:
            raise ValueError("H5MDWriter: no atoms in output trajectory")
        self._driver = driver
        if self._driver == 'mpio':
            raise ValueError("H5MDWriter: parallel writing with MPI I/O "
                             "is not currently supported.")
        self.n_atoms = n_atoms
        self.n_frames = n_frames
        self.chunks = (1, n_atoms, 3) if chunks is None else chunks
        if self.chunks is False and self.n_frames is None:
            raise ValueError("H5MDWriter must know how many frames will be "
                             "written if ``chunks=False``.")
        self.contiguous = self.chunks is False and self.n_frames is not None
        self.compression = compression
        self.compression_opts = compression_opts
        self.convert_units = convert_units
        self.h5md_file = None

        # The writer defaults to writing all data from the parent Timestep if
        # it exists. If these are True, the writer will check each
        # Timestep.has_*  value and fill the self._has dictionary accordingly
        # in _initialize_hdf5_datasets()
        self._write_positions = positions
        self._write_velocities = velocities
        self._write_forces = forces
        if not any([self._write_positions,
                    self._write_velocities,
                    self._write_velocities]):
            raise ValueError("At least one of positions, velocities, or "
                             "forces must be set to ``True``.")

        self._new_units = {'time': timeunit,
                           'length': lengthunit,
                           'velocity': velocityunit,
                           'force': forceunit}

        # Pull out various keywords to store metadata in 'h5md' group
        self.author = author
        self.author_email = author_email
        self.creator = creator
        self.creator_version = creator_version

    def _write_next_frame(self, ag):
        """Write information associated with ``ag`` at current frame
        into trajectory

        Parameters
        ----------
        ag : AtomGroup or Universe

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

        # This should only be called once when first timestep is read.
        if self.h5md_file is None:
            self._determine_units(ag)
            self._open_file()
            self._initialize_hdf5_datasets(ts)

        return self._write_next_timestep(ts)

    def _determine_units(self, ag):
        """determine which units the file will be written with

        By default, it fills the :attr:`self.units` dictionary by copying the
        units dictionary of the parent file. Because H5MD files have no
        standard unit restrictions, users may pass a kwarg in ``(timeunit,
        lengthunit, velocityunit, forceunit)`` to the writer so long as
        MDAnalysis has a conversion factor for it (:exc:`ValueError` raised if
        it does not). These custom unit arguments must be in
        `MDAnalysis notation`_. If custom units are supplied from the user,
        :attr`self.units[unit]` is replaced with the corresponding
        `unit` argument.

        """

        self.units = ag.universe.trajectory.units.copy()

        # set user input units
        for key, value in self._new_units.items():
            if value is not None:
                if value not in self._unit_translation_dict[key]:
                    raise ValueError(f"{value} is not a unit recognized by"
                                     " MDAnalysis. Allowed units are:"
                                     f" {self._unit_translation_dict.keys()}"
                                     " For more information on units, see"
                                     " `MDAnalysis units`_.")
                else:
                    self.units[key] = self._new_units[key]

        if self.convert_units:
            # check if all units are None
            if not any(self.units.values()):
                raise ValueError("The trajectory has no units, but "
                                 "`convert_units` is set to ``True`` by "
                                 "default in MDAnalysis. To write the file "
                                 "with no units, set ``convert_units=False``.")

    def _open_file(self):
        """Opens file with `H5PY`_ library and fills in metadata from kwargs.

        :attr:`self.h5md_file` becomes file handle that links to root level.

        """

        self.h5md_file = h5py.File(name=self.filename,
                                   mode='w',
                                   driver=self._driver)

        # fill in H5MD metadata from kwargs
        self.h5md_file.require_group('h5md')
        self.h5md_file['h5md'].attrs['version'] = np.array(self.H5MD_VERSION)
        self.h5md_file['h5md'].require_group('author')
        self.h5md_file['h5md/author'].attrs['name'] = self.author
        if self.author_email is not None:
            self.h5md_file['h5md/author'].attrs['email'] = self.author_email
        self.h5md_file['h5md'].require_group('creator')
        self.h5md_file['h5md/creator'].attrs['name'] = self.creator
        self.h5md_file['h5md/creator'].attrs['version'] = self.creator_version

    def _initialize_hdf5_datasets(self, ts):
        """initializes all datasets that will be written to by
        :meth:`_write_next_timestep`

        Note
        ----
        :exc:`NoDataError` is raised if no positions, velocities, or forces are
        found in the input trajectory. While the H5MD standard allows for this
        case, :class:`H5MDReader` cannot currently read files without at least
        one of these three groups. A future change to both the reader and
        writer will allow this case.


        """

        # for keeping track of where to write in the dataset
        self._counter = 0

        # ask the parent file if it has positions, velocities, and forces
        # if prompted by the writer with the self._write_* attributes
        self._has = {group: getattr(ts, f'has_{attr}')
                     if getattr(self, f'_write_{attr}')
                     else False for group, attr in zip(
                     ('position', 'velocity', 'force'),
                     ('positions', 'velocities', 'forces'))}

        # initialize trajectory group
        self.h5md_file.require_group('particles').require_group('trajectory')
        self._traj = self.h5md_file['particles/trajectory']
        self.data_keys = [
            key for key in ts.data.keys() if key not in self.data_blacklist]
        if self.data_keys:
            self._obsv = self.h5md_file.require_group('observables')

        # box group is required for every group in 'particles'
        self._traj.require_group('box')
        self._traj['box'].attrs['dimension'] = 3
        if ts.dimensions is not None and np.all(ts.dimensions > 0):
            self._traj['box'].attrs['boundary'] = 3*['periodic']
            self._traj['box'].require_group('edges')
            self._edges = self._traj.require_dataset('box/edges/value',
                                                     shape=(0, 3, 3),
                                                     maxshape=(None, 3, 3),
                                                     dtype=np.float32)
            self._step = self._traj.require_dataset('box/edges/step',
                                                    shape=(0,),
                                                    maxshape=(None,),
                                                    dtype=np.int32)
            self._time = self._traj.require_dataset('box/edges/time',
                                                    shape=(0,),
                                                    maxshape=(None,),
                                                    dtype=np.float32)
            self._set_attr_unit(self._edges, 'length')
            self._set_attr_unit(self._time, 'time')
        else:
            # if no box, boundary attr must be "none" according to H5MD
            self._traj['box'].attrs['boundary'] = 3*['none']
            self._create_step_and_time_datasets()

        if self.has_positions:
            self._create_trajectory_dataset('position')
            self._pos = self._traj['position/value']
            self._set_attr_unit(self._pos, 'length')
        if self.has_velocities:
            self._create_trajectory_dataset('velocity')
            self._vel = self._traj['velocity/value']
            self._set_attr_unit(self._vel, 'velocity')
        if self.has_forces:
            self._create_trajectory_dataset('force')
            self._force = self._traj['force/value']
            self._set_attr_unit(self._force, 'force')

        # intialize observable datasets from ts.data dictionary that
        # are NOT in self.data_blacklist
        if self.data_keys:
            for key in self.data_keys:
                self._create_observables_dataset(key, ts.data[key])

    def _create_step_and_time_datasets(self):
        """helper function to initialize a dataset for step and time

        Hunts down first available location to create the step and time
        datasets. This should only be called if the trajectory has no
        dimension, otherwise the 'box/edges' group creates step and time
        datasets since 'box' is the only required group in 'particles'.

        :attr:`self._step` and :attr`self._time` serve as links to the created
        datasets that other datasets can also point to for their step and time.
        This serves two purposes:
            1. Avoid redundant writing of multiple datasets that share the
               same step and time data.
            2. In HDF5, each chunked dataset has a cache (default 1 MiB),
               so only 1 read is required to access step and time data
               for all datasets that share the same step and time.

        """

        for group, value in self._has.items():
            if value:
                self._step = self._traj.require_dataset(f'{group}/step',
                                                        shape=(0,),
                                                        maxshape=(None,),
                                                        dtype=np.int32)
                self._time = self._traj.require_dataset(f'{group}/time',
                                                        shape=(0,),
                                                        maxshape=(None,),
                                                        dtype=np.float32)
                self._set_attr_unit(self._time, 'time')
                break

    def _create_trajectory_dataset(self, group):
        """helper function to initialize a dataset for
        position, velocity, and force"""

        if self.n_frames is None:
            shape = (0, self.n_atoms, 3)
            maxshape = (None, self.n_atoms, 3)
        else:
            shape = (self.n_frames, self.n_atoms, 3)
            maxshape = None

        chunks = None if self.contiguous else self.chunks

        self._traj.require_group(group)
        self._traj.require_dataset(f'{group}/value',
                                   shape=shape,
                                   maxshape=maxshape,
                                   dtype=np.float32,
                                   chunks=chunks,
                                   compression=self.compression,
                                   compression_opts=self.compression_opts)
        if 'step' not in self._traj[group]:
            self._traj[f'{group}/step'] = self._step
        if 'time' not in self._traj[group]:
            self._traj[f'{group}/time'] = self._time

    def _create_observables_dataset(self, group, data):
        """helper function to initialize a dataset for each observable"""

        self._obsv.require_group(group)
        # guarantee ints and floats have a shape ()
        data = np.asarray(data)
        self._obsv.require_dataset(f'{group}/value',
                                   shape=(0,) + data.shape,
                                   maxshape=(None,) + data.shape,
                                   dtype=data.dtype)
        if 'step' not in self._obsv[group]:
            self._obsv[f'{group}/step'] = self._step
        if 'time' not in self._obsv[group]:
            self._obsv[f'{group}/time'] = self._time

    def _set_attr_unit(self, dset, unit):
        """helper function to set a 'unit' attribute for an HDF5 dataset"""

        if self.units[unit] is None:
            return

        dset.attrs['unit'] = self._unit_translation_dict[unit][self.units[unit]]

    def _write_next_timestep(self, ts):
        """Write coordinates and unitcell information to H5MD file.

        Do not call this method directly; instead use
        :meth:`write` because some essential setup is done
        there before writing the first frame.

        The first dimension of each dataset is extended by +1 and
        then the data is written to the new slot.

        Note
        ----
        Writing H5MD files with fancy trajectory slicing where the Timestep
        does not increase monotonically such as ``u.trajectory[[2,1,0]]``
        or ``u.trajectory[[0,1,2,0,1,2]]`` raises a :exc:`ValueError` as this
        violates the rules of the step dataset in the H5MD standard.

        """

        i = self._counter

        # H5MD step refers to the integration step at which the data were
        # sampled, therefore ts.data['step'] is the most appropriate value
        # to use. However, step is also necessary in H5MD to allow
        # temporal matching of the data, so ts.frame is used as an alternative
        self._step.resize(self._step.shape[0]+1, axis=0)
        try:
            self._step[i] = ts.data['step']
        except(KeyError):
            self._step[i] = ts.frame
        if len(self._step) > 1 and self._step[i] < self._step[i-1]:
            raise ValueError("The H5MD standard dictates that the step "
                             "dataset must increase monotonically in value.")

        # the dataset.resize() method should work with any chunk shape
        self._time.resize(self._time.shape[0]+1, axis=0)
        self._time[i] = ts.time

        if 'edges' in self._traj['box']:
            self._edges.resize(self._edges.shape[0]+1, axis=0)
            self._edges.write_direct(ts.triclinic_dimensions,
                                     dest_sel=np.s_[i, :])
        # These datasets are not resized if n_frames was provided as an
        # argument, as they were initialized with their full size.
        if self.has_positions:
            if self.n_frames is None:
                self._pos.resize(self._pos.shape[0]+1, axis=0)
            self._pos.write_direct(ts.positions, dest_sel=np.s_[i, :])
        if self.has_velocities:
            if self.n_frames is None:
                self._vel.resize(self._vel.shape[0]+1, axis=0)
            self._vel.write_direct(ts.velocities, dest_sel=np.s_[i, :])
        if self.has_forces:
            if self.n_frames is None:
                self._force.resize(self._force.shape[0]+1, axis=0)
            self._force.write_direct(ts.forces, dest_sel=np.s_[i, :])

        if self.data_keys:
            for key in self.data_keys:
                obs = self._obsv[f'{key}/value']
                obs.resize(obs.shape[0]+1, axis=0)
                obs[i] = ts.data[key]

        if self.convert_units:
            self._convert_dataset_with_units(i)

        self._counter += 1

    def _convert_dataset_with_units(self, i):
        """convert values in the dataset arrays with self.units dictionary"""

        # Note: simply doing convert_pos_to_native(self._pos[-1]) does not
        # actually change the values in the dataset, so assignment required
        if self.units['time'] is not None:
            self._time[i] = self.convert_time_to_native(self._time[i])
        if self.units['length'] is not None:
            if self._has['position']:
                self._pos[i] = self.convert_pos_to_native(self._pos[i])
            if 'edges' in self._traj['box']:
                self._edges[i] = self.convert_pos_to_native(self._edges[i])
        if self._has['velocity']:
            if self.units['velocity'] is not None:
                self._vel[i] = self.convert_velocities_to_native(self._vel[i])
        if self._has['force']:
            if self.units['force'] is not None:
                self._force[i] = self.convert_forces_to_native(self._force[i])

    @property
    def has_positions(self):
        """``True`` if writer is writing positions from Timestep."""
        return self._has['position']

    @property
    def has_velocities(self):
        """``True`` if writer is writing velocities from Timestep."""
        return self._has['velocity']

    @property
    def has_forces(self):
        """``True`` if writer is writing forces from Timestep."""
        return self._has['force']


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
