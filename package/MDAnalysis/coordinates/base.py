# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2011 Naveen Michaud-Agrawal,
#               Elizabeth J. Denning, Oliver Beckstein,
#               and contributors (see website for details)
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
#     N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
#     O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
#     Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
#     doi:10.1002/jcc.21787
#

"""
Base classes --- :mod:`MDAnalysis.coordinates.base`
====================================================

Derive other Reader and Writer classes from the classes in this
module. The derived classes must follow the Trajectory API in
:mod:`MDAnalysis.coordinates.__init__`.

.. autoclass:: Timestep
   :members:

.. autoclass:: IObase
   :members:

.. autoclass:: Reader
   :members:

.. autoclass:: ChainReader
   :members:

.. autoclass:: Writer
   :members:

"""

import itertools
import os.path
import warnings

import bisect
import numpy

import MDAnalysis.core
from MDAnalysis.core import units, flags
from MDAnalysis.core.util import iterable, asiterable
import core

class Timestep(object):
    """Timestep data for one frame

    :Attributes:
      :attr:`numatoms`
        number of atoms
      :attr:`frame`
        frame number
      :attr:`_pos`
        coordinates as a (*numatoms*,3) :class:`numpy.ndarray` of dtype
        :data:`~numpy.float32`.

        .. Note::

           normally one does not directly access :attr:`_pos` but uses the
           :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.coordinates` method of
           an :meth:`~MDAnalysis.core.AtomGroup.AtomGroup` but sometimes it can
           be faster to directly use the raw coordinates. Any changes to this
           array are imediately reflected in atom positions. If the frame is
           written to a new trajectory then the coordinates are changed. If a
           new trajectory frame is loaded, then *all* contents of :attr:`_pos`
           are overwritten.

      :attr:`dimensions`:
         system box dimensions (A, B, C, alpha, beta, gamma); lengths
         are in the MDAnalysis length unit, and angles are in degrees.

    :Methods:
      ``t = Timestep(numatoms)``
         create a timestep object with space for numatoms (done automatically)
      ``t[i]``
         return coordinates for the i'th atom (0-based)
      ``t[start:stop:skip]``
         return an array of coordinates, where start, stop and skip correspond to atom indices (0-based)
    """
    def __init__(self, arg):
        if numpy.dtype(type(arg)) == numpy.dtype(int):
            self.frame = 0
            self.numatoms = arg
            self._pos = numpy.zeros((self.numatoms, 3), dtype=numpy.float32, order='F')
            #self._pos = numpy.zeros((3, self.numatoms), numpy.float32)
            self._unitcell = numpy.zeros((6), numpy.float32)
        elif isinstance(arg, Timestep): # Copy constructor
            # This makes a deepcopy of the timestep
            self.frame = arg.frame
            self.numatoms = arg.numatoms
            self._unitcell = numpy.array(arg._unitcell)
            self._pos = numpy.array(arg._pos, order='F')
        elif isinstance(arg, numpy.ndarray):
            if len(arg.shape) != 2: raise ValueError("numpy array can only have 2 dimensions")
            self._unitcell = numpy.zeros((6), numpy.float32)
            self.frame = 0
            #if arg.shape[0] == 3: self.numatoms = arg.shape[0]  # ??? is this correct ??? [OB]  # Nope, not sure what the aim was so i've left the lines in as comments [DP]
            #else: self.numatoms = arg.shape[-1]                 # ??? reverse ??? [OB]
            if arg.shape[1] == 3: self.numatoms = arg.shape[0]
            else: self.numatoms = arg.shape[0]   # Or should an exception be raised if coordinate structure is not 3-dimensional? Maybe velocities could be read one day... [DP]

            self._pos = arg.astype(numpy.float32).copy('Fortran',)
        else: raise ValueError("Cannot create an empty Timestep")
        self._x = self._pos[:,0]
        self._y = self._pos[:,1]
        self._z = self._pos[:,2]
    def __getitem__(self, atoms):
        if numpy.dtype(type(atoms)) == numpy.dtype(int):
            if (atoms < 0):
                atoms = self.numatoms + atoms
            if (atoms < 0) or (atoms >= self.numatoms):
                raise IndexError
            return self._pos[atoms]
        elif type(atoms) == slice or type(atoms) == numpy.ndarray:
            return self._pos[atoms]
        else: raise TypeError
    def __len__(self):
        return self.numatoms
    def __iter__(self):
        def iterTS():
            for i in xrange(self.numatoms):
                yield self[i]
        return iterTS()
    def __repr__(self):
        return "< Timestep "+ repr(self.frame) + " with unit cell dimensions " + repr(self.dimensions) + " >"
    def copy(self):
        return self.__deepcopy__()
    def __deepcopy__(self):
        # Is this the best way?
        return Timestep(self)

    @property
    def dimensions(self):
        """unitcell dimensions (A, B, C, alpha, beta, gamma)"""
        # Layout of unitcell is [A, alpha, B, beta, gamma, C] --- (originally CHARMM DCD)
        # override for other formats; this strange ordering is kept for historical reasons
        # (the user should not need concern themselves with this)
        return numpy.take(self._unitcell, [0,2,5,1,3,4])

    @property
    def volume(self):
        """volume of the unitcell"""
        return core.box_volume(self.dimensions)


class IObase(object):
    """Base class bundling common functionality for trajectory I/O.
    """
    #: override to define trajectory format of the reader/writer (DCD, XTC, ...)
    format = None

    #: dict with units of of *time* and *length* (and *velocity*, *force*,
    #: ... for formats that support it)
    units = {'time': None, 'length': None}

    def convert_pos_from_native(self, x):
        """In-place conversion of coordinate array x from native units to base units."""
        f = units.get_conversion_factor('length', self.units['length'], MDAnalysis.core.flags['length_unit'])
        x *= f
        return x

    def convert_time_from_native(self, t):
        """Convert time *t* from native units to base units."""
        f = units.get_conversion_factor('time', self.units['time'], MDAnalysis.core.flags['time_unit'])
        t *= f
        return t

    def convert_pos_to_native(self, x):
        """In-place conversion of coordinate array x from base units to native units."""
        f = units.get_conversion_factor('length', MDAnalysis.core.flags['length_unit'], self.units['length'])
        x *= f
        return x

    def convert_time_to_native(self, t):
        """Convert time *t* from base units to native units."""
        f = units.get_conversion_factor('time', MDAnalysis.core.flags['time_unit'], self.units['time'])
        t *= f
        return t

    def close(self):
        """Close the trajectory file."""
        pass

    def close_trajectory(self):
        """Specific implementation of trajectory closing."""
        # close_trajectory() was the pre-0.7.0 way of closing a
        # trajectory; it is being kept around but user code should
        # use close() and not rely on close_trajectory()
        warnings.warn("close_trajectory() will be removed in MDAnalysis 0.8. "
                      "Use close() instead.", DeprecationWarning)
        self.close()

class Reader(IObase):
    """Base class for trajectory readers.

    See the :ref:`Trajectory API` definition in
    :mod:`MDAnalysis.coordinates.__init__` for the required attributes and methods.
    """

    #: The appropriate Timestep class, e.g.
    #: :class:`MDAnalysis.coordinates.xdrfile.XTC.Timestep` for XTC.
    _Timestep = Timestep

    def __len__(self):
        return self.numframes

    def next(self):
        """Forward one step to next frame."""
        return self._read_next_timestep()

    def rewind(self):
        """Position at beginning of trajectory"""
        self[0]

    @property
    def dt(self):
        """Time between two trajectory frames in picoseconds."""
        return self.skip_timestep * self.convert_time_from_native(self.delta)

    @property
    def totaltime(self):
        """Total length of the trajectory numframes * dt."""
        return self.numframes * self.dt

    @property
    def frame(self):
        """Frame number of the current time step.

        This is a simple short cut to :attr:`Timestep.frame`.
        """
        return self.ts.frame

    @property
    def time(self):
        """Time of the current frame in MDAnalysis time units (typically ps).

        time = :attr:`Timestep.frame` * :attr:`Reader.dt`
        """
        try:
            return self.ts.frame * self.dt
        except KeyError:
            # single frame formats fail with KeyError because they do not define
            # a unit for time so we just always return 0 because time is a required
            # attribute of the Reader
            return 0.0

    def Writer(self, filename, **kwargs):
        """Returns a trajectory writer with the same properties as this trajectory."""
        raise NotImplementedError("Sorry, there is no Writer for this format in MDAnalysis. "
                                  "Please file an enhancement request at http://code.google.com/p/mdanalysis/issues/")

    def OtherWriter(self, filename, **kwargs):
        """Returns a writer appropriate for *filename*.

        Sets the default keywords *start*, *step* and *delta* (if
        available). *numatoms* is always set from :attr:`Reader.numatoms`.

        .. SeeAlso:: :meth:`Reader.Writer` and :func:`MDAnalysis.Writer`
        """
        kwargs['numatoms'] = self.numatoms            # essential
        kwargs.setdefault('start', self.frame - 1)    # -1 should be correct... [orbeckst] ?!?
        kwargs.setdefault('step', self.skip_timestep)
        try:
            kwargs.setdefault('delta', self.dt)
        except KeyError:
            pass
        return core.writer(filename, **kwargs)

    def _read_next_timestep(self, ts=None):
        # Example from DCDReader:
        #     if ts is None: ts = self.ts
        #     ts.frame = self._read_next_frame(ts._x, ts._y, ts._z, ts._unitcell, self.skip)
        #     return ts
        raise NotImplementedError("BUG: Override _read_next_timestep() in the trajectory reader!")

    def __del__(self):
        pass

    def __iter__(self):
        raise NotImplementedError("BUG: Override __iter__() in the trajectory reader!")

    def _check_slice_indices(self, start, stop, skip):
        if start == None: start = 0
        if stop == None: stop = len(self)
        if skip == None: skip = 1
        if (start < 0): start += len(self)
        if (stop < 0): stop += len(self)
        elif (stop > len(self)): stop = len(self)
        if skip > 0 and stop <= start:
            raise IndexError("Stop frame is lower than start frame")
        if ((start < 0) or (start >= len(self)) or (stop < 0) or (stop > len(self))):
            raise IndexError("Frame start/stop outside of the range of the trajectory.")
        return (start, stop, skip)

    def __repr__(self):
        return "< %s %r with %d frames of %d atoms (%d fixed) >" % \
            (self.__class__.__name__, self.filename, self.numframes, self.numatoms, self.fixed)


class ChainReader(Reader):
    """Reader that concatenates multiple trajectories on the fly.

    **Known issues**

    - Trajectory API attributes exist but most of them only reflect
      the first trajectory in the list; :attr:`ChainReader.numframes`,
      :attr:`ChainReader.numatoms`, and :attr:`ChainReader.fixed` are
      properly set, though

    - slicing not implemented

    - :attr:`time` not implemented yet
    """
    format = 'CHAIN'

    def __init__(self, filenames, **kwargs):
        """Set up the chain reader.

        :Arguments:
           *filenames*
               file name or list of file names; the reader will open
               all file names and provide frames in the order of
               trajectories from the list. Each trajectory must
               contain the same number of atoms in the same order
               (i.e. they all must belong to the same topology). The trajectory
               format is deduced from the extension of *filename*.
           *skip*
               skip step (also passed on to the individual trajectory
               readers); must be same for all trajectories
           *kwargs*
               all other keyword arguments are passed on to each
               trajectory reader unchanged
        """
        self.filenames = asiterable(filenames)
        self.readers = [core.reader(filename, **kwargs) for filename in self.filenames]
        self.__active_reader_index = 0   # pointer to "active" trajectory index into self.readers

        self.skip = kwargs.get('skip', 1)
        self.numatoms = self._get_same('numatoms')
        self.fixed = self._get_same('fixed')

        # Translation between virtual frames and frames in individual
        # trajectories.
        # Assumes that individual trajectories i contain frames that can
        # be addressed with an index 0 <= f < numframes[i]

        # Build a map of frames: ordered list of starting virtual
        # frames; the index i into this list corresponds to the index
        # into self.readers
        #
        # For virtual frame k (1...sum(numframes)) find corresponding
        # trajectory i and local frame f (i.e. readers[i][f] will
        # correspond to ChainReader[k]).

        # build map 'start_frames', which is used by _get_local_frame()
        numframes = self._get('numframes')
        # [0]: frames are 0-indexed internally
        # (see Timestep._check_slice_indices())
        self.__start_frames = numpy.cumsum([0] + numframes)

        self.numframes = numpy.sum(numframes)

        #: source for trajectories frame (fakes trajectory)
        self.__chained_trajectories_iter = None

        # make sure that iteration always yields frame 1
        # rewind() also sets self.ts
        self.ts = None
        self.rewind()

    def _get_local_frame(self, k):
        """Find trajectory index and trajectory frame for chained frame k.

        Frame *k* in the chained trajectory can be found in the
        trajectory at index *i* and frame index *f*.

        Frames are internally treated as 0-based indices into the
        trajectory. (This might not be fully consistent across
        MDAnalysis at the moment!)

        :Returns: **local frame** tuple `(i, f)`

        :Raises: :exc:`IndexError` for `k<0` or `i<0`.

        .. Note::

           Does not check if *k* is larger than the maximum number of frames in
           the chained trajectory.
        """
        if k < 0:
            raise IndexError("Virtual (chained) frames must be >= 0")
        # trajectory index i
        i = bisect.bisect_right(self.__start_frames, k) - 1
        if i < 0:
            raise IndexError("Cannot find trajectory for virtual frame %d" % k)
        # local frame index f in trajectory i (frame indices are 0-based)
        f = k - self.__start_frames[i]
        return i, f

    # methods that can change with the current reader
    def convert_time_from_native(self, t):
        return self.active_reader.convert_time_from_native(t)
    def convert_time_to_native(self, t):
        return self.active_reader.convert_time_to_native(t)
    def convert_pos_from_native(self, x):
        return self.active_reader.convert_from_native(x)
    def convert_pos_to_native(self, x):
        return self.active_reader.convert_pos_to_native(x)

    # attributes that can change with the current reader
    @property
    def filename(self):
        return self.active_reader.filename
    @property
    def skip_timestep(self):
        return self.active_reader.skip_timestep
    @property
    def delta(self):
        return self.active_reader.delta
    @property
    def periodic(self):
        return self.active_reader.periodic
    @property
    def units(self):
        return self.active_reader.units
    @property
    def compressed(self):
        try:
            return self.active_reader.compressed
        except AttributeError:
            return None

    @property
    def frame(self):
        """Cumulative frame number of the current time step.

        .. Note::

           The frame number is 1-based, i.e. the first frame has frame number
           1. However, frame indices (used for indexing and slicing with the
           `trajectory[frame_index]` notation use a 0-based index, i.e. *frame*
           - 1.
        """
        return self.ts.frame

    @property
    def time(self):
        """Cumulative time of the current frame in MDAnalysis time units (typically ps)."""
        raise NotImplementedError

    def _apply(self, method, **kwargs):
        """Execute *method* with *kwargs* for all readers."""
        return [reader.__getattribute__(method)(**kwargs) for reader in self.readers]

    def _get(self, attr):
        """Get value of *attr* for all readers."""
        return [reader.__getattribute__(attr) for reader in self.readers]

    def _get_same(self, attr):
        """Verify that *attr* has the same value for all readers and return value.

        :Arguments: *attr* attribute name
        :Returns: common value of the attribute
        :Raises: :Exc:`ValueError` if not all readers have the same value
        """
        values = numpy.array(self._get(attr))
        value = values[0]
        if not numpy.all(values == value):
            bad_traj = numpy.array(self.filenames)[values != value]
            raise ValueError("The following trajectories do not have the correct %s "
                             " (%d):\n%r" % (attr, value, bad_traj))
        return value

    def __activate_reader(self, i):
        """Make reader *i* the active reader."""
        # private method, not to be used by user to avoid a total mess
        if i < 0 or i >= len(self.readers):
            raise IndexError("Reader index must be 0 <= i < %d" % len(self.readers))
        self.__active_reader_index = i

    @property
    def active_reader(self):
        """Reader instance from which frames are being read."""
        return self.readers[self.__active_reader_index]

    def _jump_to_frame(self, frame):
        """Position trajectory at frame index *frame*.

        The frame is translated to the corresponding reader and local
        frame index and the Timestep instance in
        :attr:`ChainReader.ts` is updated.

        .. Note::

           *frame* is 0-based, i.e. the first frame in the trajectory is
           accessed with *frame* = 0.

        .. SeeAlso:: :meth:`~ChainReader._get_local_frame`.
        """
        i,f = self._get_local_frame(frame)
        # seek to (1) reader i and (2) frame f in trajectory i
        self.__activate_reader(i)
        self.active_reader[f]      # rely on reader to implement __getitem__()
        # update Timestep
        self.ts = self.active_reader.ts
        self.ts.frame = frame + 1  # continuous frames, 1-based

    def _chained_iterator(self):
        """Iterator that presents itself as a chained trajectory."""
        self._rewind()  # must rewind all readers
        readers = itertools.chain(*self.readers)
        for frame, ts in enumerate(readers):
            ts.frame = frame+1  # fake continuous frames, 1-based
            self.ts = ts
            # make sure that the active reader is in sync
            i,f = self._get_local_frame(frame)  # uses 0-based frames!
            self.__activate_reader(i)
            yield ts

    def _read_next_timestep(self, ts=None):
        self.ts = self.__chained_trajectories_iter.next()
        return self.ts

    def rewind(self):
        """Set current frame to the beginning."""
        self._rewind()
        self.__chained_trajectories_iter = self._chained_iterator()
        self.ts = self.__chained_trajectories_iter.next()  # set time step to frame 1

    def _rewind(self):
        """Internal method: Rewind trajectories themselves and trj pointer."""
        self._apply('rewind')
        self.__activate_reader(0)

    def close(self):
        self._apply('close')

    def __iter__(self):
        """Generator for all frames, starting at frame 1."""
        self._rewind()
        self.__chained_trajectories_iter = self._chained_iterator()  # start from first frame
        for ts in self.__chained_trajectories_iter:
            yield ts

    def __getitem__(self, frame):
        """Index with 0-based frame index to jump to a frame.

        Updates timestep; slices will be implemented in the future.
        """
        if numpy.dtype(type(frame)) != numpy.dtype(int) and type(frame) != slice:
            raise TypeError("frames (0-based) must be int or a slice")
        if numpy.dtype(type(frame)) == numpy.dtype(int):
            if frame < 0:
                # Interpret similar to a sequence
                frame = len(self) + frame
            if frame < 0 or frame >= len(self):
                raise IndexError("frame index out of bounds")
            self._jump_to_frame(frame)  # updates ts
            return self.ts
        elif type(frame) == slice: # if frame is a slice object
            raise NotImplementedError("Slices not implemented for ChainReader")

    def __repr__(self):
        return "< %s %r with %d frames of %d atoms (%d fixed) >" % \
            (self.__class__.__name__,
             [os.path.basename(fn) for fn in self.filenames],
             self.numframes, self.numatoms, self.fixed)



class Writer(IObase):
    """Base class for trajectory writers.

    See Trajectory API definition in :mod:`MDAnalysis.coordinates.__init__` for
    the required attributes and methods.
    """

    def convert_dimensions_to_unitcell(self, ts):
        """Read dimensions from timestep *ts* and return appropriate unitcell"""
        # override if the native trajectory format does NOT use [A,B,C,alpha,beta,gamma]
        lengths, angles = ts.dimensions[:3], ts.dimensions[3:]
        self.convert_pos_to_native(lengths)
        return numpy.concatenate([lengths, angles])

    def write(self, obj):
        """Write current timestep, using the supplied *obj*.

        The argument should be a :class:`~MDAnalysis.core.AtomGroup.AtomGroup` or
        a :class:`~MDAnalysis.Universe` or a :class:`Timestep` instance.

        .. Note::

           The size of the *obj* must be the same as the number of atom
           provided when setting up the trajectory.
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
            return "< %s %r for %d atoms >" % (self.__class__.__name__, self.filename, self.numatoms)
        except TypeError:
            # no trajectory loaded yet
            return "< %s %r >" % (self.__class__.__name__, self.filename)

    def has_valid_coordinates(self, criteria, x):
        """Returns ``True`` if all values are within limit values of their formats.

        Due to rounding, the test is asymmetric (and *min* is supposed to be negative):

           min < x <= max

        :Arguments:
            *criteria*
               dictionary containing the *max* and *min* values in native units
            *x*
               :class:`numpy.ndarray` of ``(x, y, z)`` coordinates of atoms selected to be written out.
        :Returns: boolean
        """
        x = numpy.ravel(x)
        return numpy.all(criteria["min"] < x) and numpy.all(x <= criteria["max"])

    # def write_next_timestep(self, ts=None)
