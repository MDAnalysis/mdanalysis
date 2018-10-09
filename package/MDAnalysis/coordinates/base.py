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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#


"""\
Base classes --- :mod:`MDAnalysis.coordinates.base`
===================================================

Derive other Timestep, Reader and Writer classes from the classes in
this module. The derived classes must follow the :ref:`Trajectory API`
in :mod:`MDAnalysis.coordinates.__init__`.

Timestep
--------

A :class:`Timestep` holds information for the current time frame in
the trajectory. It is one of the central data structures in
MDAnalysis.

.. class:: Timestep

   .. automethod:: __init__
   .. automethod:: from_coordinates
   .. automethod:: from_timestep
   .. automethod:: _init_unitcell
   .. autoattribute:: n_atoms
   .. attribute::`frame`

      frame number (0-based)

      .. versionchanged:: 0.11.0
         Frames now 0-based; was 1-based

   .. autoattribute:: time
   .. autoattribute:: dt
   .. autoattribute:: positions
   .. autoattribute:: velocities
   .. autoattribute:: forces
   .. autoattribute:: has_positions
   .. autoattribute:: has_velocities
   .. autoattribute:: has_forces
   .. attribute:: _pos

      :class:`numpy.ndarray` of dtype :class:`~numpy.float32` of shape
      (*n_atoms*, 3) and internal FORTRAN order, holding the raw
      cartesian coordinates (in MDAnalysis units, i.e. Å).

      .. Note::

         Normally one does not directly access :attr:`_pos` but uses
         the :meth:`~MDAnalysis.core.groups.AtomGroup.coordinates`
         method of an :class:`~MDAnalysis.core.groups.AtomGroup` but
         sometimes it can be faster to directly use the raw
         coordinates. Any changes to this array are immediately
         reflected in atom positions. If the frame is written to a new
         trajectory then the coordinates are changed. If a new
         trajectory frame is loaded, then *all* contents of
         :attr:`_pos` are overwritten.

   .. attribute:: _velocities

      :class:`numpy.ndarray` of dtype :class:`~numpy.float32`. of shape
      (*n_atoms*, 3), holding the raw velocities (in MDAnalysis
      units, i.e. typically Å/ps).

      .. Note::

         Normally velocities are accessed through the
         :attr:`velocities` or the
         :meth:`~MDAnalysis.core.groups.AtomGroup.velocities`
         method of an :class:`~MDAnalysis.core.groups.AtomGroup`

         :attr:`~Timestep._velocities` only exists if the :attr:`has_velocities`
         flag is True

      .. versionadded:: 0.7.5

   .. attribute:: _forces

      :class:`numpy.ndarray` of dtype :class:`~numpy.float32`. of shape
      (*n_atoms*, 3), holding the forces

      :attr:`~Timestep._forces` only exists if :attr:`has_forces`
      is True

      .. versionadded:: 0.11.0
         Added as optional to :class:`Timestep`

   .. autoattribute:: dimensions
   .. autoattribute:: triclinic_dimensions
   .. autoattribute:: volume
   .. automethod:: __getitem__
   .. automethod:: __eq__
   .. automethod:: __iter__
   .. automethod:: copy
   .. automethod:: copy_slice


Readers
-------

Readers know how to take trajectory data in a given format and present it in a
common API to the user in MDAnalysis. There are two types of readers:

1. Readers for *multi frame trajectories*, i.e., file formats that typically
   contain many frames. These readers are typically derived from
   :class:`ReaderBase`.

2. Readers for *single frame formats*: These file formats only contain a single
   coordinate set. These readers are derived from
   :class`:SingleFrameReaderBase`.

The underlying low-level readers handle closing of files in different
ways. Typically, the MDAnalysis readers try to ensure that files are always
closed when a reader instance is garbage collected, which relies on
implementing a :meth:`~ReaderBase.__del__` method. However, in some cases, this
is not necessary (for instance, for the single frame formats) and then such a
method can lead to undesirable side effects (such as memory leaks). In this
case, :class:`ProtoReader` should be used.


.. autoclass:: ReaderBase
   :members:
   :inherited-members:

.. autoclass:: SingleFrameReaderBase
   :members:
   :inherited-members:

.. autoclass:: ProtoReader
   :members:



Writers
-------

Writers know how to write information in a :class:`Timestep` to a trajectory
file.

.. autoclass:: WriterBase
   :members:
   :inherited-members:


Helper classes
--------------

The following classes contain basic functionality that all readers and
writers share.

.. autoclass:: IOBase
   :members:

"""
from __future__ import absolute_import
import six
from six.moves import range

import numpy as np
import numbers
import copy
import warnings
import weakref

from . import core
from .. import NoDataError
from .. import (
    _READERS,
    _SINGLEFRAME_WRITERS,
    _MULTIFRAME_WRITERS,
)
from .. import units
from ..auxiliary.base import AuxReader
from ..auxiliary.core import auxreader
from ..core import flags
from ..lib.util import asiterable, Namespace


class Timestep(object):
    """Timestep data for one frame

    :Methods:

      ``ts = Timestep(n_atoms)``

         create a timestep object with space for n_atoms

    .. versionchanged:: 0.11.0
       Added :meth:`from_timestep` and :meth:`from_coordinates` constructor
       methods.
       :class:`Timestep` init now only accepts integer creation.
       :attr:`n_atoms` now a read only property.
       :attr:`frame` now 0-based instead of 1-based.
       Attributes `status` and `step` removed.
    """
    order = 'F'

    def __init__(self, n_atoms, **kwargs):
        """Create a Timestep, representing a frame of a trajectory

        Parameters
        ----------
        n_atoms : int
          The total number of atoms this Timestep describes
        positions : bool, optional
          Whether this Timestep has position information [``True``]
        velocities : bool (optional)
          Whether this Timestep has velocity information [``False``]
        forces : bool (optional)
          Whether this Timestep has force information [``False``]
        reader : Reader (optional)
          A weak reference to the owning Reader.  Used for
          when attributes require trajectory manipulation (e.g. dt)
        dt : float (optional)
          The time difference between frames (ps).  If :attr:`time`
          is set, then `dt` will be ignored.
        time_offset : float (optional)
          The starting time from which to calculate time (in ps)


        .. versionchanged:: 0.11.0
           Added keywords for `positions`, `velocities` and `forces`.
           Can add and remove position/velocity/force information by using
           the ``has_*`` attribute.
        """
        # readers call Reader._read_next_timestep() on init, incrementing
        # self.frame to 0
        self.frame = -1
        self._n_atoms = n_atoms

        self.data = {}

        for att in ('dt', 'time_offset'):
            try:
                self.data[att] = kwargs[att]
            except KeyError:
                pass
        try:
            # do I have a hook back to the Reader?
            self._reader = weakref.ref(kwargs['reader'])
        except KeyError:
            pass

        # Stupid hack to make it allocate first time round
        # ie we have to go from not having, to having positions
        # to make the Timestep allocate
        self._has_positions = False
        self._has_velocities = False
        self._has_forces = False

        # These will allocate the arrays if the has flag
        # gets set to True
        self.has_positions = kwargs.get('positions', True)
        self.has_velocities = kwargs.get('velocities', False)
        self.has_forces = kwargs.get('forces', False)

        self._unitcell = self._init_unitcell()

        # set up aux namespace for adding auxiliary data
        self.aux = Namespace()


    @classmethod
    def from_timestep(cls, other, **kwargs):
        """Create a copy of another Timestep, in the format of this Timestep

        .. versionadded:: 0.11.0
        """
        ts = cls(other.n_atoms,
                 positions=other.has_positions,
                 velocities=other.has_velocities,
                 forces=other.has_forces,
                 **kwargs)
        ts.frame = other.frame
        ts.dimensions = other.dimensions
        try:
            ts.positions = other.positions.copy(order=cls.order)
        except NoDataError:
            pass
        try:
            ts.velocities = other.velocities.copy(order=cls.order)
        except NoDataError:
            pass
        try:
            ts.forces = other.forces.copy(order=cls.order)
        except NoDataError:
            pass

        # Optional attributes that don't live in .data
        # should probably iron out these last kinks
        for att in ('_frame',):
            try:
                setattr(ts, att, getattr(other, att))
            except AttributeError:
                pass

        if hasattr(ts, '_reader'):
            other._reader = weakref.ref(ts._reader())

        ts.data = copy.deepcopy(other.data)

        return ts

    @classmethod
    def from_coordinates(cls,
                         positions=None,
                         velocities=None,
                         forces=None,
                         **kwargs):
        """Create an instance of this Timestep, from coordinate data

        Can pass position, velocity and force data to form a Timestep.

        .. versionadded:: 0.11.0
        """
        has_positions = positions is not None
        has_velocities = velocities is not None
        has_forces = forces is not None

        lens = [len(a) for a in [positions, velocities, forces]
                if a is not None]
        if not lens:
            raise ValueError("Must specify at least one set of data")
        n_atoms = max(lens)
        # Check arrays are matched length?
        if not all(val == n_atoms for val in lens):
            raise ValueError("Lengths of input data mismatched")

        ts = cls(n_atoms,
                 positions=has_positions,
                 velocities=has_velocities,
                 forces=has_forces,
                 **kwargs)
        if has_positions:
            ts.positions = positions
        if has_velocities:
            ts.velocities = velocities
        if has_forces:
            ts.forces = forces

        return ts

    def _init_unitcell(self):
        """Create custom datastructure for :attr:`_unitcell`."""
        # override for other Timesteps
        return np.zeros((6), np.float32)

    def __eq__(self, other):
        """Compare with another Timestep

        .. versionadded:: 0.11.0
        """
        if not isinstance(other, Timestep):
            return False

        if not self.frame == other.frame:
            return False

        if not self.n_atoms == other.n_atoms:
            return False

        if not self.has_positions == other.has_positions:
            return False
        if self.has_positions:
            if not (self.positions == other.positions).all():
                return False

        if not self.has_velocities == other.has_velocities:
            return False
        if self.has_velocities:
            if not (self.velocities == other.velocities).all():
                return False

        if not self.has_forces == other.has_forces:
            return False
        if self.has_forces:
            if not (self.forces == other.forces).all():
                return False

        return True

    def __ne__(self, other):
        return not self == other

    def __getitem__(self, atoms):
        """Get a selection of coordinates

        ``ts[i]``

           return coordinates for the i'th atom (0-based)

        ``ts[start:stop:skip]``

           return an array of coordinates, where start, stop and skip
           correspond to atom indices,
           :attr:`MDAnalysis.core.groups.Atom.index` (0-based)
        """
        if isinstance(atoms, numbers.Integral):
            return self._pos[atoms]
        elif isinstance(atoms, (slice, np.ndarray)):
            return self._pos[atoms]
        else:
            raise TypeError

    def __len__(self):
        return self.n_atoms

    def __iter__(self):
        """Iterate over coordinates

        ``for x in ts``

            iterate of the coordinates, atom by atom
        """
        for i in range(self.n_atoms):
            yield self[i]

    def __repr__(self):
        desc = "< Timestep {0}".format(self.frame)
        try:
            tail = " with unit cell dimensions {0} >".format(self.dimensions)
        except NotImplementedError:
            tail = " >"
        return desc + tail

    def copy(self):
        """Make an independent ("deep") copy of the whole :class:`Timestep`."""
        return self.__deepcopy__()

    def __deepcopy__(self):
        return self.from_timestep(self)

    def copy_slice(self, sel):
        """Make a new `Timestep` containing a subset of the original `Timestep`.

        Parameters
        ----------
        sel : array_like or slice
            The underlying position, velocity, and force arrays are sliced
            using a :class:`list`, :class:`slice`, or any array-like.

        Returns
        -------
        :class:`Timestep`
            A `Timestep` object of the same type containing all header
            information and all atom information relevant to the selection.

        Note
        ----
        The selection must be a 0 based :class:`slice` or array of the atom indices
        in this :class:`Timestep`

        Example
        -------
        Using a Python :class:`slice` object::

           new_ts = ts.copy_slice(slice(start, stop, step))

        Using a list of indices::

           new_ts = ts.copy_slice([0, 2, 10, 20, 23])


        .. versionadded:: 0.8
        .. versionchanged:: 0.11.0
           Reworked to follow new Timestep API.  Now will strictly only
           copy official attributes of the Timestep.

        """
        # Detect the size of the Timestep by doing a dummy slice
        try:
            pos = self.positions[sel]
        except NoDataError:
            # It's cool if there's no Data, we'll live
            pos = None
        except:
            raise TypeError("Selection type must be compatible with slicing"
                            " the coordinates")
        try:
            vel = self.velocities[sel]
        except NoDataError:
            vel = None
        except:
            raise TypeError("Selection type must be compatible with slicing"
                            " the coordinates")
        try:
            force = self.forces[sel]
        except NoDataError:
            force = None
        except:
            raise TypeError("Selection type must be compatible with slicing"
                            " the coordinates")

        new_TS = self.__class__.from_coordinates(
            positions=pos,
            velocities=vel,
            forces=force)

        new_TS._unitcell = self._unitcell.copy()

        new_TS.frame = self.frame

        for att in ('_frame',):
            try:
                setattr(new_TS, att, getattr(self, att))
            except AttributeError:
                pass

        if hasattr(self, '_reader'):
            new_TS._reader = weakref.ref(self._reader())

        new_TS.data = copy.deepcopy(self.data)

        return new_TS

    @property
    def n_atoms(self):
        """A read only view of the number of atoms this Timestep has

        .. versionchanged:: 0.11.0
           Changed to read only property
        """
        # In future could do some magic here to make setting n_atoms
        # resize the coordinate arrays, but
        # - not sure if that is ever useful
        # - not sure how to manage existing data upon extension
        return self._n_atoms

    @property
    def has_positions(self):
        """A boolean of whether this Timestep has position data

        This can be changed to ``True`` or ``False`` to allocate space for
        or remove the data.

        .. versionadded:: 0.11.0
        """
        return self._has_positions

    @has_positions.setter
    def has_positions(self, val):
        if val and not self._has_positions:
            # Setting this will always reallocate position data
            # ie
            # True -> False -> True will wipe data from first True state
            self._pos = np.zeros((self.n_atoms, 3), dtype=np.float32,
                                 order=self.order)
            self._has_positions = True
        elif not val:
            # Unsetting val won't delete the numpy array
            self._has_positions = False

    @property
    def positions(self):
        """A record of the positions of all atoms in this Timestep

        Setting this attribute will add positions to the Timestep if they
        weren't originally present.

        Returns
        -------
        positions : numpy.ndarray with dtype numpy.float32
               position data of shape ``(n_atoms, 3)`` for all atoms

        Raises
        ------
        :exc:`MDAnalysis.exceptions.NoDataError`
               if the Timestep has no position data


        .. versionchanged:: 0.11.0
           Now can raise :exc:`NoDataError` when no position data present
        """
        if self.has_positions:
            return self._pos
        else:
            raise NoDataError("This Timestep has no positions")

    @positions.setter
    def positions(self, new):
        self.has_positions = True
        self._pos[:] = new

    @property
    def _x(self):
        """A view onto the x dimension of position data

        .. versionchanged:: 0.11.0
           Now read only
        """
        return self.positions[:, 0]

    @property
    def _y(self):
        """A view onto the y dimension of position data

        .. versionchanged:: 0.11.0
           Now read only
        """
        return self.positions[:, 1]

    @property
    def _z(self):
        """A view onto the z dimension of position data

        .. versionchanged:: 0.11.0
           Now read only
        """
        return self.positions[:, 2]

    @property
    def has_velocities(self):
        """A boolean of whether this Timestep has velocity data

        This can be changed to ``True`` or ``False`` to allocate space for
        or remove the data.

        .. versionadded:: 0.11.0
        """
        return self._has_velocities

    @has_velocities.setter
    def has_velocities(self, val):
        if val and not self._has_velocities:
            self._velocities = np.zeros((self.n_atoms, 3), dtype=np.float32,
                                        order=self.order)
            self._has_velocities = True
        elif not val:
            self._has_velocities = False

    @property
    def velocities(self):
        """A record of the velocities of all atoms in this Timestep

        Setting this attribute will add velocities to the Timestep if they
        weren't originally present.

        Returns
        -------
        velocities : numpy.ndarray with dtype numpy.float32
               velocity data of shape ``(n_atoms, 3)`` for all atoms

        Raises
        ------
        :exc:`MDAnalysis.exceptions.NoDataError`
               if the Timestep has no velocity data


        .. versionadded:: 0.11.0
        """
        if self.has_velocities:
            return self._velocities
        else:
            raise NoDataError("This Timestep has no velocities")

    @velocities.setter
    def velocities(self, new):
        self.has_velocities = True
        self._velocities[:] = new

    @property
    def has_forces(self):
        """A boolean of whether this Timestep has force data

        This can be changed to ``True`` or ``False`` to allocate space for
        or remove the data.

        .. versionadded:: 0.11.0
        """
        return self._has_forces

    @has_forces.setter
    def has_forces(self, val):
        if val and not self._has_forces:
            self._forces = np.zeros((self.n_atoms, 3), dtype=np.float32,
                                    order=self.order)
            self._has_forces = True
        elif not val:
            self._has_forces = False

    @property
    def forces(self):
        """A record of the forces of all atoms in this Timestep

        Setting this attribute will add forces to the Timestep if they
        weren't originally present.

        Returns
        -------
        forces : numpy.ndarray with dtype numpy.float32
               force data of shape ``(n_atoms, 3)`` for all atoms

        Raises
        ------
        :exc:`MDAnalysis.exceptions.NoDataError`
               if the Timestep has no force data


        .. versionadded:: 0.11.0
        """
        if self.has_forces:
            return self._forces
        else:
            raise NoDataError("This Timestep has no forces")

    @forces.setter
    def forces(self, new):
        self.has_forces = True
        self._forces[:] = new

    @property
    def dimensions(self):
        """unitcell dimensions (*A*, *B*, *C*, *alpha*, *beta*, *gamma*)

        lengths *a*, *b*, *c* are in the MDAnalysis length unit (Å), and
        angles are in degrees.

        Setting dimensions will populate the underlying native format
        description of box dimensions
        """
        # The actual Timestep._unitcell depends on the underlying
        # trajectory format. It can be e.g. six floats representing
        # the box edges and angles or the 6 unique components of the
        # box matrix or the full box matrix.
        return self._unitcell

    @dimensions.setter
    def dimensions(self, box):
        self._unitcell[:] = box

    @property
    def volume(self):
        """volume of the unitcell"""
        return core.box_volume(self.dimensions)

    @property
    def triclinic_dimensions(self):
        """The unitcell dimensions represented as triclinic vectors

        Returns
        -------
        numpy.ndarray
             A (3, 3) numpy.ndarray of unit cell vectors

        Examples
        --------
        The unitcell for a given system can be queried as either three
        vectors lengths followed by their respective angle, or as three
        triclinic vectors.

          >>> ts.dimensions
          array([ 13.,  14.,  15.,  90.,  90.,  90.], dtype=float32)
          >>> ts.triclinic_dimensions
          array([[ 13.,   0.,   0.],
                 [  0.,  14.,   0.],
                 [  0.,   0.,  15.]], dtype=float32)

        Setting the attribute also works::

          >>> ts.triclinic_dimensions = [[15, 0, 0], [5, 15, 0], [5, 5, 15]]
          >>> ts.dimensions
          array([ 15.        ,  15.81138802,  16.58312416,  67.58049774,
                  72.45159912,  71.56504822], dtype=float32)

        See Also
        --------
        :func:`MDAnalysis.lib.mdamath.triclinic_vectors`


        .. versionadded:: 0.11.0
        """
        return core.triclinic_vectors(self.dimensions)

    @triclinic_dimensions.setter
    def triclinic_dimensions(self, new):
        """Set the unitcell for this Timestep as defined by triclinic vectors

        .. versionadded:: 0.11.0
        """
        self.dimensions = core.triclinic_box(*new)

    @property
    def dt(self):
        """The time difference in ps between timesteps

        Note
        ----
        This defaults to 1.0 ps in the absence of time data


        .. versionadded:: 0.11.0
        """
        try:
            return self.data['dt']
        except KeyError:
            pass
        try:
            dt = self.data['dt'] = self._reader()._get_dt()
            return dt
        except AttributeError:
            pass
        warnings.warn("Reader has no dt information, set to 1.0 ps")
        return 1.0

    @dt.setter
    def dt(self, new):
        self.data['dt'] = new

    @dt.deleter
    def dt(self):
        del self.data['dt']

    @property
    def time(self):
        """The time in ps of this timestep

        This is calculated as::

          time = ts.data['time_offset'] + ts.time

        Or, if the trajectory doesn't provide time information::

          time = ts.data['time_offset'] + ts.frame * ts.dt

        .. versionadded:: 0.11.0
        """
        offset = self.data.get('time_offset', 0)
        try:
            return self.data['time'] + offset
        except KeyError:
            return self.dt * self.frame + offset

    @time.setter
    def time(self, new):
        self.data['time'] = new

    @time.deleter
    def time(self):
        del self.data['time']


class FrameIteratorBase(object):
    """
    Base iterable over the frames of a trajectory.

    A frame iterable has a length that can be accessed with the :func:`len`
    function, and can be indexed similarly to a full trajectory. When indexed,
    indices are resolved relative to the iterable and not relative to the
    trajectory.

    Parameters
    ----------
    trajectory: ProtoReader
        The trajectory over which to iterate.

    .. versionadded:: 0.19.0

    """
    def __init__(self, trajectory):
        self._trajectory = trajectory

    def __len__(self):
        raise NotImplementedError()

    @staticmethod
    def _avoid_bool_list(frames):
        if isinstance(frames, list) and frames and isinstance(frames[0], bool):
            return np.array(frames, dtype=bool)
        return frames

    @property
    def trajectory(self):
        return self._trajectory


class FrameIteratorSliced(FrameIteratorBase):
    """
    Iterable over the frames of a trajectory on the basis of a slice.

    Parameters
    ----------
    trajectory: ProtoReader
        The trajectory over which to iterate.
    frames: slice
        A slice to select the frames of interest.

    See Also
    --------
    FrameIteratorBase

    .. versionadded:: 0.19.0

    """
    def __init__(self, trajectory, frames):
        # It would be easier to store directly a range object, as it would
        # store its parameters in a single place, calculate its length, and
        # take care of most the indexing. Though, doing so is not compatible
        # with python 2 where xrange (or range with six) is only an iterator.
        super(FrameIteratorSliced, self).__init__(trajectory)
        self._start, self._stop, self._step = trajectory.check_slice_indices(
            frames.start, frames.stop, frames.step,
        )

    def __len__(self):
        start, stop, step = self.start, self.stop, self.step
        if (step > 0 and start < stop):
            # We go from a lesser number to a larger one.
            return int(1 + (stop - 1 - start) // step)
        elif (step < 0 and start > stop):
            # We count backward from a larger number to a lesser one.
            return int(1 + (start - 1 - stop) // (-step))
        else:
            # The range is empty.
            return 0

    def __iter__(self):
        for i in range(self.start, self.stop, self.step):
            yield self.trajectory[i]
        self.trajectory.rewind()

    def __getitem__(self, frame):
        if isinstance(frame, numbers.Integral):
            length = len(self)
            if not -length < frame < length:
                raise IndexError('Index {} is out of range of the range of length {}.'
                                 .format(frame, length))
            if frame < 0:
                frame = len(self) + frame
            frame = self.start + frame * self.step
            return self.trajectory._read_frame_with_aux(frame)
        elif isinstance(frame, slice):
            start = self.start + (frame.start or 0) * self.step
            if frame.stop is None:
                stop = self.stop
            else:
                stop = self.start + (frame.stop or 0) * self.step
            step = (frame.step or 1) * self.step

            if step > 0:
                start = max(0, start)
            else:
                stop = max(0, stop)

            new_slice = slice(start, stop, step)
            return FrameIteratorSliced(self.trajectory, new_slice)
        else:
            # Indexing with a lists of bools does not behave the same in all
            # version of numpy.
            frame = self._avoid_bool_list(frame)
            frames = np.array(list(range(self.start, self.stop, self.step)))[frame]
            return FrameIteratorIndices(self.trajectory, frames)

    @property
    def start(self):
        return self._start

    @property
    def stop(self):
        return self._stop

    @property
    def step(self):
        return self._step


class FrameIteratorAll(FrameIteratorBase):
    """
    Iterable over all the frames of a trajectory.

    Parameters
    ----------
    trajectory: ProtoReader
        The trajectory over which to iterate.

    See Also
    --------
    FrameIteratorBase

    .. versionadded:: 0.19.0

    """
    def __init__(self, trajectory):
        super(FrameIteratorAll, self).__init__(trajectory)

    def __len__(self):
        return self.trajectory.n_frames

    def __iter__(self):
        return iter(self.trajectory)

    def __getitem__(self, frame):
        return self.trajectory[frame]


class FrameIteratorIndices(FrameIteratorBase):
    """
    Iterable over the frames of a trajectory listed in a sequence of indices.

    Parameters
    ----------
    trajectory: ProtoReader
        The trajectory over which to iterate.
    frames: sequence
        A sequence of indices.

    See Also
    --------
    FrameIteratorBase
    """
    def __init__(self, trajectory, frames):
        super(FrameIteratorIndices, self).__init__(trajectory)
        self._frames = []
        for frame in frames:
            if not isinstance(frame, numbers.Integral):
                raise TypeError("Frames indices must be integers.")
            frame = trajectory._apply_limits(frame)
            self._frames.append(frame)
        self._frames = tuple(self._frames)

    def __len__(self):
        return len(self.frames)

    def __iter__(self):
        for frame in self.frames:
            yield self.trajectory._read_frame_with_aux(frame)

    def __getitem__(self, frame):
        if isinstance(frame, numbers.Integral):
            frame = self.frames[frame]
            return self.trajectory._read_frame_with_aux(frame)
        else:
            frame = self._avoid_bool_list(frame)
            frames = np.array(self.frames)[frame]
            return FrameIteratorIndices(self.trajectory, frames)

    @property
    def frames(self):
        return self._frames


class IOBase(object):
    """Base class bundling common functionality for trajectory I/O.

    .. versionchanged:: 0.8
       Added context manager protocol.
    """

    #: dict with units of of *time* and *length* (and *velocity*, *force*,
    #: ... for formats that support it)
    units = {'time': None, 'length': None, 'velocity': None}

    def convert_pos_from_native(self, x, inplace=True):
        """Conversion of coordinate array x from native units to base units.

        Parameters
        ----------
        x : array_like
          Positions to transform
        inplace : bool (optional)
          Whether to modify the array inplace, overwriting previous data

        Note
        ----
        By default, the input `x` is modified in place and also returned.
        In-place operations improve performance because allocating new arrays
        is avoided.


        .. versionchanged:: 0.7.5
           Keyword `inplace` can be set to ``False`` so that a
           modified copy is returned *unless* no conversion takes
           place, in which case the reference to the unmodified `x` is
           returned.

        """
        f = units.get_conversion_factor(
            'length', self.units['length'], flags['length_unit'])
        if f == 1.:
            return x
        if not inplace:
            return f * x
        x *= f
        return x

    def convert_velocities_from_native(self, v, inplace=True):
        """Conversion of velocities array *v* from native to base units

        Parameters
        ----------
        v : array_like
          Velocities to transform
        inplace : bool (optional)
          Whether to modify the array inplace, overwriting previous data

        Note
        ----
        By default, the input *v* is modified in place and also returned.
        In-place operations improve performance because allocating new arrays
        is avoided.


        .. versionadded:: 0.7.5
        """
        f = units.get_conversion_factor(
            'speed', self.units['velocity'], flags['speed_unit'])
        if f == 1.:
            return v
        if not inplace:
            return f * v
        v *= f
        return v

    def convert_forces_from_native(self, force, inplace=True):
        """Conversion of forces array *force* from native to base units

        Parameters
        ----------
        force : array_like
          Forces to transform
        inplace : bool (optional)
          Whether to modify the array inplace, overwriting previous data

        Note
        ----
        By default, the input *force* is modified in place and also returned.
        In-place operations improve performance because allocating new arrays
        is avoided.

        .. versionadded:: 0.7.7
        """
        f = units.get_conversion_factor(
            'force', self.units['force'], flags['force_unit'])
        if f == 1.:
            return force
        if not inplace:
            return f * force
        force *= f
        return force

    def convert_time_from_native(self, t, inplace=True):
        """Convert time *t* from native units to base units.

        Parameters
        ----------
        t : array_like
          Time values to transform
        inplace : bool (optional)
          Whether to modify the array inplace, overwriting previous data

        Note
        ----
        By default, the input `t` is modified in place and also returned
        (although note that scalar values `t` are passed by value in Python and
        hence an in-place modification has no effect on the caller.)  In-place
        operations improve performance because allocating new arrays is
        avoided.


        .. versionchanged:: 0.7.5
           Keyword `inplace` can be set to ``False`` so that a
           modified copy is returned *unless* no conversion takes
           place, in which case the reference to the unmodified `x` is
           returned.

        """
        f = units.get_conversion_factor(
            'time', self.units['time'], flags['time_unit'])
        if f == 1.:
            return t
        if not inplace:
            return f * t
        t *= f
        return t

    def convert_pos_to_native(self, x, inplace=True):
        """Conversion of coordinate array `x` from base units to native units.

        Parameters
        ----------
        x : array_like
          Positions to transform
        inplace : bool (optional)
          Whether to modify the array inplace, overwriting previous data

        Note
        ----
        By default, the input `x` is modified in place and also returned.
        In-place operations improve performance because allocating new arrays
        is avoided.


        .. versionchanged:: 0.7.5
           Keyword `inplace` can be set to ``False`` so that a
           modified copy is returned *unless* no conversion takes
           place, in which case the reference to the unmodified `x` is
           returned.

        """
        f = units.get_conversion_factor(
            'length', flags['length_unit'], self.units['length'])
        if f == 1.:
            return x
        if not inplace:
            return f * x
        x *= f
        return x

    def convert_velocities_to_native(self, v, inplace=True):
        """Conversion of coordinate array `v` from base to native units

        Parameters
        ----------
        v : array_like
          Velocities to transform
        inplace : bool (optional)
          Whether to modify the array inplace, overwriting previous data

        Note
        ----
        By default, the input `v` is modified in place and also returned.
        In-place operations improve performance because allocating new arrays
        is avoided.


        .. versionadded:: 0.7.5
        """
        f = units.get_conversion_factor(
            'speed', flags['speed_unit'], self.units['velocity'])
        if f == 1.:
            return v
        if not inplace:
            return f * v
        v *= f
        return v

    def convert_forces_to_native(self, force, inplace=True):
        """Conversion of force array *force* from base to native units.

        Parameters
        ----------
        force : array_like
          Forces to transform
        inplace : bool (optional)
          Whether to modify the array inplace, overwriting previous data

        Note
        ----
        By default, the input `force` is modified in place and also returned.
        In-place operations improve performance because allocating new arrays
        is avoided.


        .. versionadded:: 0.7.7
        """
        f = units.get_conversion_factor(
            'force', flags['force_unit'], self.units['force'])
        if f == 1.:
            return force
        if not inplace:
            return f * force
        force *= f
        return force

    def convert_time_to_native(self, t, inplace=True):
        """Convert time *t* from base units to native units.

        Parameters
        ----------
        t : array_like
          Time values to transform
        inplace : bool, optional
          Whether to modify the array inplace, overwriting previous data

        Note
        ----
        By default, the input *t* is modified in place and also
        returned. (Also note that scalar values *t* are passed by
        value in Python and hence an in-place modification has no
        effect on the caller.)

        .. versionchanged:: 0.7.5
           Keyword *inplace* can be set to ``False`` so that a
           modified copy is returned *unless* no conversion takes
           place, in which case the reference to the unmodified *x* is
           returned.

        """
        f = units.get_conversion_factor(
            'time', flags['time_unit'], self.units['time'])
        if f == 1.:
            return t
        if not inplace:
            return f * t
        t *= f
        return t

    def close(self):
        """Close the trajectory file."""
        pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # see http://docs.python.org/2/library/stdtypes.html#typecontextmanager
        self.close()
        return False  # do not suppress exceptions


class _Readermeta(type):
    # Auto register upon class creation
    def __init__(cls, name, bases, classdict):
        type.__init__(type, name, bases, classdict)
        try:
            fmt = asiterable(classdict['format'])
        except KeyError:
            pass
        else:
            for f in fmt:
                f = f.upper()
                _READERS[f] = cls


class ProtoReader(six.with_metaclass(_Readermeta, IOBase)):
    """Base class for Readers, without a :meth:`__del__` method.

    Extends :class:`IOBase` with most attributes and methods of a generic
    Reader, with the exception of a :meth:`__del__` method. It should be used
    as base for Readers that do not need :meth:`__del__`, especially since
    having even an empty :meth:`__del__` might lead to memory leaks.

    See the :ref:`Trajectory API` definition in
    :mod:`MDAnalysis.coordinates.__init__` for the required attributes and
    methods.

    See Also
    --------
    :class:`ReaderBase`


    .. versionchanged:: 0.11.0
       Frames now 0-based instead of 1-based
    """

    #: The appropriate Timestep class, e.g.
    #: :class:`MDAnalysis.coordinates.xdrfile.XTC.Timestep` for XTC.
    _Timestep = Timestep

    def __init__(self):
        # initialise list to store added auxiliary readers in
        # subclasses should now call super
        self._auxs = {}
        self._transformations=[]

    def __len__(self):
        return self.n_frames

    @classmethod
    def parse_n_atoms(cls, filename, **kwargs):
        """Read the coordinate file and deduce the number of atoms

        Returns
        -------
        n_atoms : int
          the number of atoms in the coordinate file

        Raises
        ------
        NotImplementedError
          when the number of atoms can't be deduced
        """
        raise NotImplementedError("{} cannot deduce the number of atoms"
                                  "".format(cls.__name__))

    def next(self):
        """Forward one step to next frame."""
        try:
            ts = self._read_next_timestep()
        except (EOFError, IOError):
            self.rewind()
            raise StopIteration
        else:
            for auxname in self.aux_list:
                ts = self._auxs[auxname].update_ts(ts)

            ts = self._apply_transformations(ts)

        return ts

    def __next__(self):
        """Forward one step to next frame when using the `next` builtin."""
        return self.next()

    def rewind(self):
        """Position at beginning of trajectory"""
        self._reopen()
        self.next()

    @property
    def dt(self):
        """Time between two trajectory frames in picoseconds."""
        return self.ts.dt

    @property
    def totaltime(self):
        """Total length of the trajectory

        The time is calculated as ``(n_frames - 1) * dt``, i.e., we assume that
        the first frame no time as elapsed. Thus, a trajectory with two frames will
        be considered to have a length of a single time step `dt` and a
        "trajectory" with a single frame will be reported as length 0.

        """
        return (self.n_frames - 1) * self.dt

    @property
    def frame(self):
        """Frame number of the current time step.

        This is a simple short cut to :attr:`Timestep.frame`.
        """
        return self.ts.frame

    @property
    def time(self):
        """Time of the current frame in MDAnalysis time units (typically ps).

        This is either read straight from the Timestep, or calculated as
        time = :attr:`Timestep.frame` * :attr:`Timestep.dt`
        """
        return self.ts.time

    @property
    def trajectory(self):
        # Makes a reader effectively commpatible with a FrameIteratorBase
        return self

    def Writer(self, filename, **kwargs):
        """A trajectory writer with the same properties as this trajectory."""
        raise NotImplementedError(
            "Sorry, there is no Writer for this format in MDAnalysis. "
            "Please file an enhancement request at "
            "https://github.com/MDAnalysis/mdanalysis/issues")

    def OtherWriter(self, filename, **kwargs):
        """Returns a writer appropriate for *filename*.

        Sets the default keywords *start*, *step* and *dt* (if
        available). *n_atoms* is always set from :attr:`Reader.n_atoms`.


        See Also
        --------
        :meth:`Reader.Writer` and :func:`MDAnalysis.Writer`

        """
        kwargs['n_atoms'] = self.n_atoms  # essential
        kwargs.setdefault('start', self.frame)
        try:
            kwargs.setdefault('dt', self.dt)
        except KeyError:
            pass
        return core.writer(filename, **kwargs)

    def _read_next_timestep(self, ts=None):  # pragma: no cover
        # Example from DCDReader:
        #     if ts is None:
        #         ts = self.ts
        #     ts.frame = self._read_next_frame(etc)
        #     return ts
        raise NotImplementedError(
            "BUG: Override _read_next_timestep() in the trajectory reader!")

    def __iter__(self):
        """ Iterate over trajectory frames. """
        self._reopen()
        return self

    def _reopen(self):
        """Should position Reader to just before first frame

        Calling next after this should return the first frame
        """
        pass

    def _apply_limits(self, frame):
        if frame < 0:
            frame += len(self)
        if frame < 0 or frame >= len(self):
            raise IndexError("Index {} exceeds length of trajectory ({})."
                             "".format(frame, len(self)))
        return frame

    def __getitem__(self, frame):
        """Return the Timestep corresponding to *frame*.

        If `frame` is a integer then the corresponding frame is
        returned. Negative numbers are counted from the end.

        If frame is a :class:`slice` then an iterator is returned that
        allows iteration over that part of the trajectory.

        Note
        ----
        *frame* is a 0-based frame index.
        """
        if isinstance(frame, numbers.Integral):
            frame = self._apply_limits(frame)
            return self._read_frame_with_aux(frame)
        elif isinstance(frame, (list, np.ndarray)):
            if len(frame) != 0 and isinstance(frame[0], (bool, np.bool_)):
                # Avoid having list of bools
                frame = np.asarray(frame, dtype=np.bool)
                # Convert bool array to int array
                frame = np.arange(len(self))[frame]
            return FrameIteratorIndices(self, frame)
        elif isinstance(frame, slice):
            start, stop, step = self.check_slice_indices(
                frame.start, frame.stop, frame.step)
            if start == 0 and stop == len(self) and step == 1:
                return FrameIteratorAll(self)
            else:
                return FrameIteratorSliced(self, frame)
        else:
            raise TypeError("Trajectories must be an indexed using an integer,"
                            " slice or list of indices")

    def _read_frame(self, frame):
        """Move to *frame* and fill timestep with data."""
        raise TypeError("{0} does not support direct frame indexing."
                        "".format(self.__class__.__name__))
        # Example implementation in the DCDReader:
        # self._jump_to_frame(frame)
        # ts = self.ts
        # ts.frame = self._read_next_frame(ts._x, ts._y, ts._z,
        #                                  ts._unitcell, 1)
        # return ts

    def _read_frame_with_aux(self, frame):
        """Move to *frame*, updating ts with trajectory, transformations and auxiliary data."""
        ts = self._read_frame(frame)  # pylint: disable=assignment-from-no-return
        for aux in self.aux_list:
            ts = self._auxs[aux].update_ts(ts)

        ts = self._apply_transformations(ts)

        return ts

    def _sliced_iter(self, start, stop, step):
        """Generator for slicing a trajectory.

        *start* *stop* and *step* are 3 integers describing the slice.
        Error checking is not done past this point.

        A :exc:`NotImplementedError` is raised if random access to
        frames is not implemented.
        """
        # override with an appropriate implementation e.g. using self[i] might
        # be much slower than skipping steps in a next() loop
        try:
            for i in range(start, stop, step):
                yield self._read_frame_with_aux(i)
            self.rewind()
        except TypeError:  # if _read_frame not implemented
            raise TypeError("{0} does not support slicing."
                            "".format(self.__class__.__name__))

    def check_slice_indices(self, start, stop, step):
        """Check frame indices are valid and clip to fit trajectory.

        The usage follows standard Python conventions for :func:`range` but see
        the warning below.

        Parameters
        ----------
        start : int or None
          Starting frame index (inclusive). ``None`` corresponds to the default
          of 0, i.e., the initial frame.
        stop : int or None
          Last frame index (exclusive). ``None`` corresponds to the default
          of n_frames, i.e., it includes the last frame of the trajectory.
        step : int or None
          step size of the slice, ``None`` corresponds to the default of 1, i.e,
          include every frame in the range `start`, `stop`.

        Returns
        -------
        start, stop, step : tuple (int, int, int)
          Integers representing the slice

        Warning
        -------
        The returned values `start`, `stop` and `step` give the expected result
        when passed in :func:`range` but gives unexpected behavior when passed
        in a :class:`slice` when ``stop=None`` and ``step=-1``

        This can be a problem for downstream processing of the output from this
        method. For example, slicing of trajectories is implemented by passing
        the values returned by :meth:`check_slice_indices` to :func:`range` ::

          range(start, stop, step)

        and using them as the indices to randomly seek to. On the other hand,
        in :class:`MDAnalysis.analysis.base.AnalysisBase` the values returned
        by :meth:`check_slice_indices` are used to splice the trajectory by
        creating a :class:`slice` instance ::

          slice(start, stop, step)

        This creates a discrepancy because these two lines are not equivalent::

            range(10, -1, -1)             # [10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0]
            range(10)[slice(10, -1, -1)]  # []

        """

        slice_dict = {'start': start, 'stop': stop, 'step': step}
        for varname, var in slice_dict.items():
            if isinstance(var, numbers.Integral):
                slice_dict[varname] = int(var)
            elif (var is None):
                pass
            else:
                raise TypeError("{0} is not an integer".format(varname))

        start = slice_dict['start']
        stop = slice_dict['stop']
        step = slice_dict['step']

        if step == 0:
            raise ValueError("Step size is zero")

        nframes = len(self)
        step = step or 1

        if start is None:
            start = 0 if step > 0 else nframes - 1
        elif start < 0:
            start += nframes
        if start < 0:
            start = 0

        if step < 0 and start >= nframes:
            start = nframes - 1

        if stop is None:
            stop = nframes if step > 0 else -1
        elif stop < 0:
            stop += nframes

        if step > 0 and stop > nframes:
            stop = nframes

        return start, stop, step

    def __repr__(self):
        return ("<{cls} {fname} with {nframes} frames of {natoms} atoms>"
                "".format(
            cls=self.__class__.__name__,
            fname=self.filename,
            nframes=self.n_frames,
            natoms=self.n_atoms
        ))

    def add_auxiliary(self, auxname, auxdata, format=None, **kwargs):
        """Add auxiliary data to be read alongside trajectory.

        Auxiliary data may be any data timeseries from the trajectory additional
        to that read in by the trajectory reader. *auxdata* can be an
        :class:`~MDAnalysis.auxiliary.base.AuxReader` instance, or the data
        itself as e.g. a filename; in the latter case an appropriate
        :class:`~MDAnalysis.auxiliary.base.AuxReader` is guessed from the
        data/file format. An appropriate `format` may also be directly provided
        as a key word argument.

        On adding, the AuxReader is initially matched to the current timestep
        of the trajectory, and will be updated when the trajectory timestep
        changes (through a call to :meth:`next()` or jumping timesteps with
        ``trajectory[i]``).

        The representative value(s) of the auxiliary data for each timestep (as
        calculated by the :class:`~MDAnalysis.auxiliary.base.AuxReader`) are
        stored in the current timestep in the ``ts.aux`` namespace under *auxname*;
        e.g. to add additional pull force data stored in pull-force.xvg::

            u = MDAnalysis.Universe(PDB, XTC)
            u.trajectory.add_auxiliary('pull', 'pull-force.xvg')

        The representative value for the current timestep may then be accessed
        as ``u.trajectory.ts.aux.pull`` or ``u.trajectory.ts.aux['pull']``.

        See Also
        --------
        :meth:`remove_auxiliary`

        Note
        ----
        Auxiliary data is assumed to be time-ordered, with no duplicates. See
        the :ref:`Auxiliary API`.
        """
        if auxname in self.aux_list:
            raise ValueError("Auxiliary data with name {name} already "
                             "exists".format(name=auxname))
        if isinstance(auxdata, AuxReader):
            aux = auxdata
            aux.auxname = auxname
        else:
            aux = auxreader(auxdata, format=format, auxname=auxname, **kwargs)
        self._auxs[auxname] = aux
        self.ts = aux.update_ts(self.ts)

    def remove_auxiliary(self, auxname):
        """Clear data and close the :class:`~MDAnalysis.auxiliary.base.AuxReader`
        for the auxiliary *auxname*.

        See Also
        --------
        :meth:`add_auxiliary`
        """
        aux = self._check_for_aux(auxname)
        aux.close()
        del aux
        delattr(self.ts.aux, auxname)

    @property
    def aux_list(self):
        """ Lists the names of added auxiliary data. """
        return self._auxs.keys()

    def _check_for_aux(self, auxname):
        """ Check for the existance of an auxiliary *auxname*. If present,
        return the AuxReader; if not, raise ValueError
        """
        if auxname in self.aux_list:
            return self._auxs[auxname]
        else:
            raise ValueError("No auxiliary named {name}".format(name=auxname))

    def next_as_aux(self, auxname):
        """ Move to the next timestep for which there is at least one step from
        the auxiliary *auxname* within the cutoff specified in *auxname*.

        This allows progression through the trajectory without encountering
        ``NaN`` representative values (unless these are specifically part of the
        auxiliary data).

        If the auxiliary cutoff is not set, where auxiliary steps are less frequent
        (``auxiliary.dt > trajectory.dt``), this allows progression at the
        auxiliary pace (rounded to nearest timestep); while if the auxiliary
        steps are more frequent, this will work the same as calling
        :meth:`next()`.

        See the :ref:`Auxiliary API`.

        See Also
        --------
        :meth:`iter_as_aux`
        """

        aux = self._check_for_aux(auxname)
        ts = self.ts
        # catch up auxiliary if it starts earlier than trajectory
        while aux.step_to_frame(aux.step + 1, ts) is None:
            next(aux)
        # find the next frame that'll have a representative value
        next_frame = aux.next_nonempty_frame(ts)
        if next_frame is None:
            # no more frames with corresponding auxiliary values; stop iteration
            raise StopIteration
        # some readers set self._frame to -1, rather than self.frame, on
        # _reopen; catch here or doesn't read first frame
        while self.frame != next_frame or getattr(self, '_frame', 0) == -1:
            # iterate trajectory until frame is reached
            ts = self.next()
        return ts

    def iter_as_aux(self, auxname):
        """Iterate through timesteps for which there is at least one assigned
        step from the auxiliary *auxname* within the cutoff specified in *auxname*.

        See Also
        --------
        :meth:`next_as_aux`
        :meth:`iter_auxiliary`
        """
        aux = self._check_for_aux(auxname)
        self._reopen()
        aux._restart()
        while True:
            try:
                yield self.next_as_aux(auxname)
            except StopIteration:
                return

    def iter_auxiliary(self, auxname, start=None, stop=None, step=None,
                       selected=None):
        """ Iterate through the auxiliary *auxname* independently of the trajectory.

        Will iterate over the specified steps of the auxiliary (defaults to all
        steps). Allows to access all values in an auxiliary, including those out
        of the time range of the trajectory, without having to also iterate
        through the trajectory.

        After interation, the auxiliary will be repositioned at the current step.

        Parameters
        ----------
        auxname : str
            Name of the auxiliary to iterate over.
        (start, stop, step) : optional
            Options for iterating over a slice of the auxiliary.
        selected : lst | ndarray, optional
            List of steps to iterate over.

        Yields
        ------
        :class:`~MDAnalysis.auxiliary.base.AuxStep` object

        See Also
        --------
        :meth:`iter_as_aux`
        """
        aux = self._check_for_aux(auxname)
        if selected is not None:
            selection = selected
        else:
            selection = slice(start, stop, step)
        for i in aux[selection]:
            yield i
        aux.read_ts(self.ts)

    def get_aux_attribute(self, auxname, attrname):
        """Get the value of *attrname* from the auxiliary *auxname*

        Parameters
        ----------
        auxname : str
            Name of the auxiliary to get value for
        attrname : str
            Name of gettable attribute in the auxiliary reader

        See Also
        --------
        :meth:`set_aux_attribute`
        """
        aux = self._check_for_aux(auxname)
        return getattr(aux, attrname)

    def set_aux_attribute(self, auxname, attrname, new):
        """ Set the value of *attrname* in the auxiliary *auxname*.

        Parameters
        ----------
        auxname : str
            Name of the auxiliary to alter
        attrname : str
            Name of settable attribute in the auxiliary reader
        new
            New value to try set *attrname* to

        See Also
        --------
        :meth:`get_aux_attribute`
        :meth:`rename_aux` - to change the *auxname* attribute
        """
        aux = self._check_for_aux(auxname)
        if attrname == 'auxname':
            self.rename_aux(auxname, new)
        else:
            setattr(aux, attrname, new)

    def rename_aux(self, auxname, new):
        """ Change the name of the auxiliary *auxname* to *new*.

        Provided there is not already an auxiliary named *new*, the auxiliary
        name will be changed in ts.aux namespace, the trajectory's
        list of added auxiliaries, and in the auxiliary reader itself.

        Parameters
        ----------
        auxname : str
             Name of the auxiliary to rename
        new : str
             New name to try set

        Raises
        ------
        ValueError
             If the name *new* is already in use by an existing auxiliary.
        """
        aux = self._check_for_aux(auxname)
        if new in self.aux_list:
            raise ValueError("Auxiliary data with name {name} already "
                             "exists".format(name=new))
        aux.auxname = new
        self._auxs[new] = self._auxs.pop(auxname)
        setattr(self.ts.aux, new, self.ts.aux[auxname])
        delattr(self.ts.aux, auxname)

    def get_aux_descriptions(self, auxnames=None):
        """Get descriptions to allow reloading the specified auxiliaries.

        If no auxnames are provided, defaults to the full list of added
        auxiliaries.

        Passing the resultant description to ``add_auxiliary()`` will allow
        recreation of the auxiliary. e.g., to duplicate all auxiliaries into a
        second trajectory::

           descriptions = trajectory_1.get_aux_descriptions()
           for aux in descriptions:
               trajectory_2.add_auxiliary(**aux)


        Returns
        -------
        list
            List of dictionaries of the args/kwargs describing each auxiliary.

        See Also
        --------
        :meth:`MDAnalysis.auxiliary.base.AuxReader.get_description`
        """
        if not auxnames:
            auxnames = self.aux_list
        descriptions = [self._auxs[aux].get_description() for aux in auxnames]
        return descriptions

    @property
    def transformations(self):
        """ Returns the list of transformations"""
        return self._transformations

    @transformations.setter
    def transformations(self, transformations):
        if not self._transformations:
            self._transformations = transformations
        else:
            raise ValueError("Transformations are already set")

    def add_transformations(self, *transformations):
        """ Add all transformations to be applied to the trajectory.

        This function take as list of transformations as an argument. These
        transformations are functions that will be called by the Reader and given
        a :class:`Timestep` object as argument, which will be transformed and returned
        to the Reader.
        The transformations can be part of the :mod:`~MDAnalysis.transformations`
        module, or created by the user, and are stored as a list `transformations`.
        This list can only be modified once, and further calls of this function will
        raise an exception.

        .. code-block:: python

          u = MDAnalysis.Universe(topology, coordinates)
          workflow = [some_transform, another_transform, this_transform]
          u.trajectory.add_transformations(*workflow)

        Parameters
        ----------
        transform_list : list
            list of all the transformations that will be applied to the coordinates

        See Also
        --------
        :mod:`MDAnalysis.transformations`
        """

        try:
            self.transformations = transformations
        except ValueError:
            raise ValueError("Can't add transformations again. Please create new Universe object")
        else:
            self.ts = self._apply_transformations(self.ts)


        # call reader here to apply the newly added transformation on the
        # current loaded frame?

    def _apply_transformations(self, ts):
        """Applies all the transformations given by the user """

        for transform in self.transformations:
            ts = transform(ts)

        return ts



class ReaderBase(ProtoReader):
    """Base class for trajectory readers that extends :class:`ProtoReader` with a
    :meth:`__del__` method.

    New Readers should subclass :class:`ReaderBase` and properly implement a
    :meth:`close` method, to ensure proper release of resources (mainly file
    handles). Readers that are inherently safe in this regard should subclass
    :class:`ProtoReader` instead.

    See the :ref:`Trajectory API` definition in
    :mod:`MDAnalysis.coordinates.__init__` for the required attributes and
    methods.

    See Also
    --------
    :class:`ProtoReader`


    .. versionchanged:: 0.11.0
       Most of the base Reader class definitions were offloaded to
       :class:`ProtoReader` so as to allow the subclassing of ReaderBases without a
       :meth:`__del__` method.  Created init method to create common
       functionality, all ReaderBase subclasses must now :func:`super` through this
       class.  Added attribute :attr:`_ts_kwargs`, which is created in init.
       Provides kwargs to be passed to :class:`Timestep`

    """

    def __init__(self, filename, convert_units=None, **kwargs):
        super(ReaderBase, self).__init__()

        self.filename = filename

        if convert_units is None:
            convert_units = flags['convert_lengths']
        self.convert_units = convert_units

        ts_kwargs = {}
        for att in ('dt', 'time_offset'):
            try:
                val = kwargs[att]
            except KeyError:
                pass
            else:
                ts_kwargs[att] = val

        self._ts_kwargs = ts_kwargs

    def copy(self):
        """Return independent copy of this Reader.

        New Reader will have its own file handle and can seek/iterate
        independently of the original.

        Will also copy the current state of the Timestep held in
        the original Reader
        """
        new = self.__class__(self.filename,
                             n_atoms=self.n_atoms)
        if self.transformations:
            new.add_transformations(*self.transformations)
        # seek the new reader to the same frame we started with
        new[self.ts.frame]
        # then copy over the current Timestep in case it has
        # been modified since initial load
        new.ts = self.ts.copy()
        for auxname, auxread in self._auxs.items():
            new.add_auxiliary(auxname, auxread.copy())
        return new

    def __del__(self):
        for aux in self.aux_list:
            self._auxs[aux].close()
        self.close()


class _Writermeta(type):
    # Auto register this format upon class creation
    def __init__(cls, name, bases, classdict):
        type.__init__(type, name, bases, classdict)
        try:
            # grab the string which describes this format
            # could be either 'PDB' or ['PDB', 'ENT'] for multiple formats
            fmt = asiterable(classdict['format'])
        except KeyError:
            # not required however
            pass
        else:
            # does the Writer support single and multiframe writing?
            single = classdict.get('singleframe', True)
            multi = classdict.get('multiframe', False)

            if single:
                for f in fmt:
                    f = f.upper()
                    _SINGLEFRAME_WRITERS[f] = cls
            if multi:
                for f in fmt:
                    f = f.upper()
                    _MULTIFRAME_WRITERS[f] = cls


class WriterBase(six.with_metaclass(_Writermeta, IOBase)):
    """Base class for trajectory writers.

    See Trajectory API definition in :mod:`MDAnalysis.coordinates.__init__` for
    the required attributes and methods.
    """

    def convert_dimensions_to_unitcell(self, ts, inplace=True):
        """Read dimensions from timestep *ts* and return appropriate unitcell.

        The default is to return ``[A,B,C,alpha,beta,gamma]``; if this
        is not appropriate then this method has to be overriden.
        """
        # override if the native trajectory format does NOT use
        # [A,B,C,alpha,beta,gamma]
        lengths, angles = ts.dimensions[:3], ts.dimensions[3:]
        if not inplace:
            lengths = lengths.copy()
        lengths = self.convert_pos_to_native(lengths)
        return np.concatenate([lengths, angles])

    def write(self, obj):
        """Write current timestep, using the supplied `obj`.

        Parameters
        ----------
        obj : :class:`~MDAnalysis.core.groups.AtomGroup` or :class:`~MDAnalysis.core.universe.Universe` or a :class:`Timestep`
            write coordinate information associate with `obj`

        Note
        ----
        The size of the `obj` must be the same as the number of atoms provided
        when setting up the trajectory.
        """
        if isinstance(obj, Timestep):
            ts = obj
        else:
            try:
                ts = obj.ts
            except AttributeError:
                try:
                    # special case: can supply a Universe, too...
                    ts = obj.trajectory.ts
                except AttributeError:
                    raise TypeError("No Timestep found in obj argument")
        return self.write_next_timestep(ts)

    def __del__(self):
        self.close()

    def __repr__(self):
        try:
            return "< {0!s} {1!r} for {2:d} atoms >".format(self.__class__.__name__, self.filename, self.n_atoms)
        except (TypeError, AttributeError):
            # no trajectory loaded yet or a Writer that does not need e.g.
            # self.n_atoms
            return "< {0!s} {1!r} >".format(self.__class__.__name__, self.filename)

    def has_valid_coordinates(self, criteria, x):
        """Returns ``True`` if all values are within limit values of their formats.

        Due to rounding, the test is asymmetric (and *min* is supposed to be negative):

           min < x <= max

        Parameters
        ----------
        criteria : dict
            dictionary containing the *max* and *min* values in native units
        x : numpy.ndarray
            ``(x, y, z)`` coordinates of atoms selected to be written out

        Returns
        -------
        bool
        """
        x = np.ravel(x)
        return np.all(criteria["min"] < x) and np.all(x <= criteria["max"])

        # def write_next_timestep(self, ts=None)


class SingleFrameReaderBase(ProtoReader):
    """Base class for Readers that only have one frame.

    To use this base class, define the method :meth:`_read_first_frame` to
    read from file `self.filename`.  This should populate the attribute
    `self.ts` with a :class:`Timestep` object.

    .. versionadded:: 0.10.0
    .. versionchanged:: 0.11.0
       Added attribute "_ts_kwargs" for subclasses
       Keywords "dt" and "time_offset" read into _ts_kwargs
    """
    _err = "{0} only contains a single frame"

    def __init__(self, filename, convert_units=None, n_atoms=None, **kwargs):
        super(SingleFrameReaderBase, self).__init__()

        self.filename = filename
        if convert_units is None:
            convert_units = flags['convert_lengths']
        self.convert_units = convert_units

        self.n_frames = 1
        self.n_atom = n_atoms

        ts_kwargs = {}
        for att in ('dt', 'time_offset'):
            try:
                val = kwargs[att]
            except KeyError:
                pass
            else:
                ts_kwargs[att] = val

        self._ts_kwargs = ts_kwargs
        self._read_first_frame()

    def copy(self):
        """Return independent copy of this Reader.

        New Reader will have its own file handle and can seek/iterate
        independently of the original.

        Will also copy the current state of the Timestep held in
        the original Reader
        """
        new = self.__class__(self.filename,
                             n_atoms=self.n_atoms)
        new.ts = self.ts.copy()
        for auxname, auxread in self._auxs.items():
            new.add_auxiliary(auxname, auxread.copy())
        # since the transformations have already been applied to the frame
        # simply copy the property
        new.transformations = self.transformations

        return new

    def _read_first_frame(self):  # pragma: no cover
        # Override this in subclasses to create and fill a Timestep
        pass

    def rewind(self):
        self._read_first_frame()
        for auxname, auxread in self._auxs.items():
            self.ts = auxread.update_ts(self.ts)
        super(SingleFrameReaderBase, self)._apply_transformations(self.ts)

    def _reopen(self):
        pass

    def next(self):
        raise StopIteration(self._err.format(self.__class__.__name__))

    def __iter__(self):
        yield self.ts
        return

    def _read_frame(self, frame):
        if frame != 0:
            raise IndexError(self._err.format(self.__class__.__name__))

        return self.ts

    def close(self):
        # all single frame readers should use context managers to access
        # self.filename. Explicitly setting it to the null action in case
        # the IOBase.close method is ever changed from that.
        pass

    def add_transformations(self, *transformations):
        """ Add all transformations to be applied to the trajectory.

        This function take as list of transformations as an argument. These
        transformations are functions that will be called by the Reader and given
        a :class:`Timestep` object as argument, which will be transformed and returned
        to the Reader.
        The transformations can be part of the :mod:`~MDAnalysis.transformations`
        module, or created by the user, and are stored as a list `transformations`.
        This list can only be modified once, and further calls of this function will
        raise an exception.

        .. code-block:: python

          u = MDAnalysis.Universe(topology, coordinates)
          workflow = [some_transform, another_transform, this_transform]
          u.trajectory.add_transformations(*workflow)

        Parameters
        ----------
        transform_list : list
            list of all the transformations that will be applied to the coordinates

        See Also
        --------
        :mod:`MDAnalysis.transformations`
        """
        #Overrides :meth:`~MDAnalysis.coordinates.base.ProtoReader.add_transformations`
        #to avoid unintended behaviour where the coordinates of each frame are transformed
        #multiple times when iterating over the trajectory.
        #In this method, the trajectory is modified all at once and once only.

        super(SingleFrameReaderBase, self).add_transformations(*transformations)
        for transform in self.transformations:
            self.ts = transform(self.ts)

    def _apply_transformations(self, ts):
        """ Applies the transformations to the timestep."""
        # Overrides :meth:`~MDAnalysis.coordinates.base.ProtoReader.add_transformations`
        # to avoid applying the same transformations multiple times on each frame

        return ts
