# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; encoding: utf-8 -*-
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

   .. attribute:: _pos

      :class:`numpy.ndarray` of dtype :class:`~numpy.float32` of shape
      (*numatoms*, 3) and internal FORTRAN order, holding the raw
      cartesian coordinates (in MDAnalysis units, i.e. Å).

      .. Note::

         Normally one does not directly access :attr:`_pos` but uses
         the :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.coordinates`
         method of an :class:`~MDAnalysis.core.AtomGroup.AtomGroup` but
         sometimes it can be faster to directly use the raw
         coordinates. Any changes to this array are immediately
         reflected in atom positions. If the frame is written to a new
         trajectory then the coordinates are changed. If a new
         trajectory frame is loaded, then *all* contents of
         :attr:`_pos` are overwritten.

   .. attribute:: _velocities

      :class:`numpy.ndarray` of dtype :class:`~numpy.float32`. of shape
      (*numatoms*, 3), holding the raw velocities (in MDAnalysis
      units, i.e. typically Å/ps).

      .. Note::

         Normally velocities are accessed through the
         :meth:`~MDAnalysis.core.AtomGroup.AtomGroup.velocities()`
         method of an :class:`~MDAnalysis.core.AtomGroup.AtomGroup`
         but this attribute is documented as there can be occasions
         when it is required (e.g. in order to *change* velocities) or
         much more convenient or faster to access the raw velocities
         directly.

         :attr:`~Timestep._velocities` only exists if the underlying
         trajectory format supports velocities. Your code should check
         for its existence or handle an :exc:`AttributeError`
         gracefully.


      .. versionadded:: 0.7.5

   .. attribute:: numatoms

      number of atoms

   .. attribute::`frame`

      frame number

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

    :Methods:

      ``ts = Timestep(numatoms)``

         create a timestep object with space for numatoms (done
         automatically)

      ``ts[i]``

         return coordinates for the i'th atom (0-based)

      ``ts[start:stop:skip]``

         return an array of coordinates, where start, stop and skip
         correspond to atom indices,
         :attr:`MDAnalysis.core.AtomGroup.Atom.number` (0-based)

      ``for x in ts``

         iterate of the coordinates, atom by atom
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
            if len(arg.shape) != 2:
                raise ValueError("numpy array can only have 2 dimensions")
            self._unitcell = numpy.zeros((6), numpy.float32)
            self.frame = 0
            if arg.shape[1] == 3:
                self.numatoms = arg.shape[0]
            else:
                self.numatoms = arg.shape[0]
                # Or should an exception be raised if coordinate
                # structure is not 3-dimensional? Maybe velocities
                # could be read one day... [DP]
            self._pos = arg.astype(numpy.float32).copy('Fortran',)
        else:
            raise ValueError("Cannot create an empty Timestep")
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
        """Make an independent ("deep") copy of the whole :class:`Timestep`."""
        return self.__deepcopy__()
    def __deepcopy__(self):
        # Is this the best way?
        return Timestep(self)

    @property
    def dimensions(self):
        """unitcell dimensions (*A*, *B*, *C*, *alpha*, *beta*, *gamma*)

        lengths *a*, *b*, *c* are in the MDAnalysis length unit (Å), and
        angles are in degrees.

        :attr:`dimensions` is read-only because it transforms the
        actual format of the unitcell (which differs between different
        trajectory formats) to the representation described here,
        which is used everywhere in MDAnalysis.
        """

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
    units = {'time': None, 'length': None, 'velocity': None}

    def convert_pos_from_native(self, x, inplace=True):
        """In-place conversion of coordinate array x from native units to base units.

        By default, the input *x* is modified in place and also returned.

        .. versionchanged:: 0.7.5
           Keyword *inplace* can be set to ``False`` so that a
           modified copy is returned *unless* no conversion takes
           place, in which case the reference to the unmodified *x* is
           returned.

        """
        f = units.get_conversion_factor('length', self.units['length'], MDAnalysis.core.flags['length_unit'])
        if f == 1.:
            return x
        if not inplace:
            return f*x
        x *= f
        return x

    def convert_velocities_from_native(self, v, inplace=True):
        """In-place conversion of velocities array *v* from native units to base units.

        By default, the input *v* is modified in place and also returned.

        .. versionadded:: 0.7.5
        """
        f = units.get_conversion_factor('speed', self.units['velocity'], MDAnalysis.core.flags['speed_unit'])
        if f == 1.:
            return v
        if not inplace:
            return f*v
        v *= f
        return v

    def convert_forces_from_native(self, force, inplace=True):
        """In-place conversion of forces array *force* from native units to base units.

        By default, the input *force* is modified in place and also returned.

        .. versionadded:: 0.7.7
        """
        f = units.get_conversion_factor('force', self.units['force'], MDAnalysis.core.flags['force_unit'])
        if f == 1.:
            return force
        if not inplace:
            return f*force
        force *= f
        return force

    def convert_time_from_native(self, t, inplace=True):
        """Convert time *t* from native units to base units.

        By default, the input *t* is modified in place and also
        returned (although note that scalar values *t* are passed by
        value in Python and hence an in-place modification has no
        effect on the caller.)

        .. versionchanged:: 0.7.5
           Keyword *inplace* can be set to ``False`` so that a
           modified copy is returned *unless* no conversion takes
           place, in which case the reference to the unmodified *x* is
           returned.

        """
        f = units.get_conversion_factor('time', self.units['time'], MDAnalysis.core.flags['time_unit'])
        if f == 1.:
            return t
        if not inplace:
            return f*t
        t *= f
        return t

    def convert_pos_to_native(self, x, inplace=True):
        """Conversion of coordinate array x from base units to native units.

        By default, the input *x* is modified in place and also returned.

        .. versionchanged:: 0.7.5
           Keyword *inplace* can be set to ``False`` so that a
           modified copy is returned *unless* no conversion takes
           place, in which case the reference to the unmodified *x* is
           returned.

        """
        f = units.get_conversion_factor('length', MDAnalysis.core.flags['length_unit'], self.units['length'])
        if f == 1.:
            return x
        if not inplace:
            return f*x
        x *= f
        return x

    def convert_velocities_to_native(self, v, inplace=True):
        """In-place conversion of coordinate array *v* from base units to native units.

        By default, the input *v* is modified in place and also returned.

        .. versionadded:: 0.7.5
        """
        f = units.get_conversion_factor('speed', MDAnalysis.core.flags['speed_unit'], self.units['velocity'])
        if f == 1.:
            return v
        if not inplace:
            return f*v
        v *= f
        return v

    def convert_forces_to_native(self, force, inplace=True):
        """In-place conversion of force array *force* from base units to native units.

        By default, the input *force* is modified in place and also returned.

        .. versionadded:: 0.7.7
        """
        f = units.get_conversion_factor('force', MDAnalysis.core.flags['force_unit'], self.units['force'])
        if f == 1.:
            return force
        if not inplace:
            return f*force
        force *= f
        return force

    def convert_time_to_native(self, t, inplace=True):
        """Convert time *t* from base units to native units.

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
        f = units.get_conversion_factor('time', MDAnalysis.core.flags['time_unit'], self.units['time'])
        if f == 1.:
            return t
        if not inplace:
            return f*t
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

    def __getitem__(self, frame):
        """Return the Timestep corresponding to *frame*.

        If *frame* is a integer then the corresponding frame is
        returned. Negative numbers are counted from the end.

        If frame is a :class:`slice` then an iterator is returned that
        allows iteration over that part of the trajectory.

        .. Note:: *frame* is a 0-based frame index.
        """
        if (numpy.dtype(type(frame)) != numpy.dtype(int)) and (type(frame) != slice):
            raise TypeError("The frame index (0-based) must be either an integer or a slice")
        if (numpy.dtype(type(frame)) == numpy.dtype(int)):
            if (frame < 0):
                # Interpret similar to a sequence
                frame = len(self) + frame
            if (frame < 0) or (frame >= len(self)):
                raise IndexError("Index %d exceeds length of trajectory (%d)." % (frame, len(self)))
            return self._read_frame(frame)  # REPLACE WITH APPROPRIATE IMPLEMENTATION
        elif type(frame) == slice: # if frame is a slice object
            if not (((type(frame.start) == int) or (frame.start == None)) and
                    ((type(frame.stop) == int) or (frame.stop == None)) and
                    ((type(frame.step) == int) or (frame.step == None))):
                raise TypeError("Slice indices are not integers")
            return self._sliced_iter(frame)

    def _read_frame(self, frame):
        """Move to *frame* and fill timestep with data."""
        raise TypeError("Reader does not support direct frame indexing.""")
        # Example implementation in the DCDReader:
        #self._jump_to_frame(frame)
        #ts = self.ts
        #ts.frame = self._read_next_frame(ts._x, ts._y, ts._z, ts._unitcell, 1) # XXX required!!
        #return ts

    def _sliced_iter(self, frames):
        """Generator for slicing a trajectory.

        *frames* must be a :class:`slice` object or behave like one.

        If *frames* corresponds to the whole trajectory then the
        standard iterator is used; this should work for any trajectory
        reader. A :exc:`TypeError` is raised if random access to
        frames is not implemented.
        """
        # override with an appropriate implementation e.g. using self[i] might
        # be much slower than skipping steps in a next() loop
        start, stop, step = self._check_slice_indices(frames.start, frames.stop, frames.step)
        if start == 0 and stop == len(self) and step == 1:
            # standard iterator (always implemented)
            _iter = self.__iter__
        else:
            # real slicing
            try:
                # Test if the iterator will work: need to do this here so that code can
                # immediately know if slicing works, even before running through the iterator.
                self[0]  # should raise TypeError if it cannot do random access
            except TypeError:
                raise TypeError("Reader does not support slicing.")
            def _iter(start=start, stop=stop, step=step):
                for i in xrange(start, stop, step):
                    yield self[i]
        return _iter()

    def _check_slice_indices(self, start, stop, step):
        if start == None: start = 0
        if stop == None: stop = len(self)
        if step == None: step = 1
        if (start < 0): start += len(self)
        if (stop < 0): stop += len(self)
        elif (stop > len(self)): stop = len(self)
        if step > 0 and stop <= start:
            raise IndexError("Stop frame is lower than start frame")
        if ((start < 0) or (start >= len(self)) or (stop < 0) or (stop > len(self))):
            raise IndexError("Frame start/stop outside of the range of the trajectory.")
        return (start, stop, step)

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

    - :attr:`time` will not necessarily return the true time but just
      number of frames times a provided time between frames (from the
      keyword *delta*)
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

               Extension: filenames are either single filename or list of file names in either plain file names format or (filename,format) tuple combination

           *skip*
               skip step (also passed on to the individual trajectory
               readers); must be same for all trajectories

           *delta*
               The time between frames in MDAnalysis time units if no
               other information is available. If this is not set then
               any call to :attr:`~ChainReader.time` will raise a
               :exc:`ValueError`.

           *kwargs*
               all other keyword arguments are passed on to each
               trajectory reader unchanged

        .. versionchanged:: 0.8
           The *delta* keyword was added.
        """
        self.filenames = asiterable(filenames)
        self.readers = [core.reader(filename, **kwargs) for filename in self.filenames]
        self.__active_reader_index = 0   # pointer to "active" trajectory index into self.readers

        self.skip = kwargs.get('skip', 1)
        self._default_delta = kwargs.pop('delta', None)
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
        # currently a hack, should really use a list of trajectory lengths and delta * local_frame
        try:
            return self.frame * self._default_delta
        except TypeError:
            raise ValueError("No timestep information available. Set delta to fake a constant time step.")

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
            bad_traj = numpy.array([self.get_flname(fn) for fn in self.filenames])[values != value]
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

    def _read_frame(self, frame):
        """Position trajectory at frame index *frame* and return :class:`Timestep`.

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
        return self.ts

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

    def get_flname(self, filename): #retrieve the actual filename of the list element
        return filename[0] if isinstance(filename,tuple) else filename

    def __repr__(self):
        return "< %s %r with %d frames of %d atoms (%d fixed) >" % \
            (self.__class__.__name__,
             [os.path.basename(self.get_flname(fn)) for fn in self.filenames],
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
        except (TypeError, AttributeError):
            # no trajectory loaded yet or a Writer that does not need e.g. self.numatoms
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
