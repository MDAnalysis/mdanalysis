# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#


"""
Gromacs TRR trajectory I/O --- :mod:`MDAnalysis.coordinates.xdrfile.TRR`
========================================================================

Classes for reading and writing of `Gromacs TRR trajectories`_
together with supporting code.

.. note:: Users should access classes from :mod:`MDAnalysis.coordinates.TRR`.

.. _Gromacs TRR trajectories: http://www.gromacs.org/Documentation/File_Formats/.trr_File
.. _Gromacs: http://www.gromacs.org


.. SeeAlso:: :mod:`MDAnalysis.coordinates.xdrfile.libxdrfile2` for low-level
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

import numpy as np
import errno

from . import statno
from . import core
from . import libxdrfile2
from MDAnalysis import NoDataError


class Timestep(core.Timestep):
    """Timestep for a Gromacs_ TRR trajectory.

    TRR Timestep always has positions, velocities and forces allocated, meaning
    that the `_pos`, `_velocities` and `_forces` attributes will always exist.
    Whether or not this data is valid for the current frame is determined by the
    flags `has_positions`, `has_velocities` and `has_forces`.  These are controlled
    by the entries in the data dictionary: `position_source`, `velocity_source`
    and `force_source`.
    When accessing the `.positions` attribute, this will only be returned if
    `position_source` matches the current `frame`, otherwise a :class:`~MDAnalysis.NoDataError`
    will be raised.  This scheme applies to both velocities and forces.

    When writing to the Timestep, it is important that the :attr:`~Timestep.frame`
    attribute is updated first!  This is because the source of the data is
    taken as the current frame of the Timestep.

    For example::

        ts = MDAnalysis.coordinates.TRR.Timestep(N)     # N being the number of atoms
        ts.frame = 0  # The first frame
        ts.velocities = vel_array   # Where vel_array is an existing array of shape (N, DIM)
                                    # This will also automatically set the `has_velocities` flag
                                    # And the velocity source will be "0"

    Attempting to populate the array before setting the frame instead will cause the 
    velocity data to be considered "out of date"::

        ts = MDAnalysis.coordinates.TRR.Timestep(N)     # N being the number of atoms
        ts.velocities = vel_array
        ts.frame = 0  # This frame number doesn't match the frame number of when
                      # velocities was populated

    .. versionchanged:: 0.8.0
       TRR :class:`Timestep` objects are now fully aware of the existence or
       not of coordinate/velocity/force information in frames.
    .. versionchanged:: 0.11.0
       Velocities and forces can now be directly accessed via the `velocities` and
       `forces` attributes, with a :class:`~MDAnalysis.NoDataError` being raised
       if no data was given for that frame.
       Management redone to use `has_positions`, `has_velocities` and
       `has_forces`.  Setting these will update the data dictionary entries
       which track the data.
       Timestep lmbda now stored in the data dictionary 

    .. _Gromacs: http://www.gromacs.org
    """
    # The exception error
    _nodataerr = ("You are accessing the {attrs} of a Timestep but there are none in"
                  " this frame. It might be the case that your trajectory doesn't"
                  " have {attrs}, or not every frame. In the latter case you can "
                  "(1) get rid of {attr}-less frames before passing the trajectory"
                  " to MDAnalysis, or "
                  "(2) reflow your code to not access {attrs} when they're not there,"
                  " by making use of the '{flag}' flag of Timestep objects.")

    def __init__(self, n_atoms, **kwargs):
        kwargs.update(positions=True, velocities=True, forces=True)
        super(Timestep, self).__init__(n_atoms, **kwargs)

        # Set initial sources to None, so they are never valid
        # _source is compared against current frame when accessing
        # if they do not match, then the Timestep returns _nodataerr
        self.data['position_source'] = None
        self.data['velocity_source'] = None
        self.data['force_source'] = None

        # TRR always has pos vel & force allocated
        self._pos = np.zeros((self.n_atoms, 3), dtype=np.float32,
                             order=self.order)
        self._velocities = np.zeros((self.n_atoms, 3), dtype=np.float32,
                                    order=self.order)
        self._forces = np.zeros((self.n_atoms, 3), dtype=np.float32,
                                order=self.order)

        self.data['lmbda'] = 0

    @property
    def has_positions(self):
        # When setting position the frame is recorded
        # if frame gets changed (by reading next frame)
        # then this data is "out of date", so has_positions
        # is False
        return self.data['position_source'] == self.frame

    @has_positions.setter
    def has_positions(self, val):
        # If val evaluates to True, say that position data
        # comes from this frame
        # If val evaluates to False, say that position data
        # comes from -1 (ie. never)
        self.data['position_source'] = self.frame if val else None

    @property
    def positions(self):
        # If the position info came from this frame, yield it
        # else play dumb
        if self.has_positions:
            return self._pos
        else:
            raise NoDataError(self._nodataerr.format(
                attr='position', attrs='positions', flag='has_positions'))

    @positions.setter
    def positions(self, new):
        self._pos[:] = new
        self.has_positions = True

    @property
    def has_velocities(self):
        return self.data['velocity_source'] == self.frame

    @has_velocities.setter
    def has_velocities(self, val):
        self.data['velocity_source'] = self.frame if val else None

    @property
    def velocities(self):
        if self.has_velocities:
            return self._velocities
        else:
            raise NoDataError(self._nodataerr.format(
                attr='velocity', attrs='velocities', flag='has_velocities'))

    @velocities.setter
    def velocities(self, new):
        self._velocities[:] = new
        self.has_velocities = True

    @property
    def has_forces(self):
        return self.data['force_source'] == self.frame

    @has_forces.setter
    def has_forces(self, val):
        self.data['force_source'] = self.frame if val else None

    @property
    def forces(self):
        if self.has_forces:
            return self._forces
        else:
            raise NoDataError(self._nodataerr.format(
                attr='force', attrs='forces', flag='has_forces'))

    @forces.setter
    def forces(self, new):
        self._forces[:] = new
        self.has_forces = True


class TRRWriter(core.TrjWriter):
    """Write a Gromacs_ TRR trajectory.

    .. _Gromacs: http://www.gromacs.org
    """
    format = "TRR"
    units = {'time': 'ps', 'length': 'nm', 'velocity': 'nm/ps', 'force': 'kJ/(mol*nm)'}


class TRRReader(core.TrjReader):
    """Read a Gromacs_ TRR trajectory.

    .. versionchanged:: 0.8.0
       :class:`Timestep` objects returned from TRR files now have
       :attr:`~Timestep.has_x`, :attr:`~Timestep.has_velocities`, and :attr:`~Timestep.has_forces`
       flags reflecting whether coordinates/velocities/forces were read.
       Attempting to access such data when the corresponding flag is set to ``False``
       will raise a :exc:`NoDataError`.
    .. versionchanged:: 0.11.0
       Returned Timesteps will always be a :class:`Timestep` instance.
       Attributes lmbda and status now stored in the TS.data dictionary

    .. _Gromacs: http://www.gromacs.org
    """
    format = "TRR"
    _Timestep = Timestep
    _Writer = TRRWriter
    units = {'time': 'ps', 'length': 'nm', 'velocity': 'nm/ps', 'force': 'kJ/(mol*nm)'}

    def _allocate_sub(self, DIM):
        self._pos_buf = np.zeros((self._trr_n_atoms, DIM), dtype=np.float32, order='C')
        self._velocities_buf = np.zeros((self._trr_n_atoms, DIM), dtype=np.float32, order='C')
        self._forces_buf = np.zeros((self._trr_n_atoms, DIM), dtype=np.float32, order='C')
    
    def _read_trj_natoms(self, filename):
        return libxdrfile2.read_trr_natoms(filename)
    
    def _read_trj_n_frames(self, filename):
        self._n_frames, self._offsets = libxdrfile2.read_trr_n_frames(filename)
        self._store_offsets()

    def _read_next_timestep(self, ts=None):
        if ts is None:
            ts = self.ts
        elif not isinstance(ts, Timestep):
            # Requires using a TRR Timestep
            ts = Timestep.from_timestep(ts)

        if self.xdrfile is None:
            self.open_trajectory()

        if self._sub is None:
            ts.data['status'], ts._frame, ts.time, ts.data['lmbda'],\
                has_x, has_v, has_f = libxdrfile2.read_trr(
                    self.xdrfile, ts._unitcell, ts._pos, ts._velocities, ts._forces)
        else:
            ts.data['status'], ts._frame, ts.time, ts.data['lmbda'],\
                has_x, has_v, has_f = libxdrfile2.read_trr(
                    self.xdrfile, ts._unitcell, self._pos_buf, self._velocities_buf, self._forces_buf)
            ts._pos[:] = self._pos_buf[self._sub]
            ts._velocities[:] = self._velocities_buf[self._sub]
            ts._forces[:] = self._forces_buf[self._sub]

        # Updated frame number first!
        ts.frame += 1
        # Update the sources of data
        if has_x:
            ts.has_positions = True
        if has_v:
            ts.has_velocities = True
        if has_f:
            ts.has_forces = True

        if ((ts.data['status'] == libxdrfile2.exdrENDOFFILE) or 
            (ts.data['status'] == libxdrfile2.exdrINT)):
            # seems that trr files can get a exdrINT when reaching EOF (??)
            raise IOError(errno.EIO, "End of file reached for {0} file".format(self.format),
                          self.filename)
        elif not ts.data['status'] == libxdrfile2.exdrOK:
            raise IOError(errno.EBADF,("Problem with {0} file, status {1}"
                                       "".format((self.format, statno.ERRORCODE[ts.data['status']]), self.filename)))

        if self.convert_units:
            # TRRs have the annoying possibility of frames without coordinates/velocities/forces...
            if ts.has_positions:
                self.convert_pos_from_native(ts._pos)  # in-place !
            self.convert_pos_from_native(ts._unitcell)  # in-place ! (note: trr contain unit vecs!)
            ts.time = self.convert_time_from_native(ts.time)  # in-place does not work with scalars
            if ts.has_velocities:
                self.convert_velocities_from_native(ts._velocities)  # in-place
            if ts.has_forces:
                self.convert_forces_from_native(ts._forces)  # in-place

        return ts
