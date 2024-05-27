# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

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


"""\
Base classes --- :mod:`MDAnalysis.coordinates.base`
===================================================

Derive, FrameIterator, Reader and Writer classes from the classes
in this module. The derived classes must follow the :ref:`Trajectory API`.


.. _FrameIterators:

FrameIterators
--------------

FrameIterators are "sliced trajectories" (a trajectory is a
:ref:`Reader <Readers>`) that can be iterated over. They are typically
created by slicing a trajectory or by fancy-indexing of a trajectory
with an array of frame numbers or a boolean mask of all frames.

Iterator classes used by the by the :class:`ProtoReader`:

.. autoclass:: FrameIteratorBase

.. autoclass:: FrameIteratorSliced

.. autoclass:: FrameIteratorAll

.. autoclass:: FrameIteratorIndices


.. _ReadersBase:

Readers
-------

Readers know how to take trajectory data in a given format and present it in a
common API to the user in MDAnalysis. There are two types of readers:

1. Readers for *multi frame trajectories*, i.e., file formats that typically
   contain many frames. These readers are typically derived from
   :class:`ReaderBase`.

2. Readers for *single frame formats*: These file formats only contain a single
   coordinate set. These readers are derived from
   :class:`SingleFrameReaderBase`.

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



.. _WritersBase:

Writers
-------

Writers know how to write information in a :class:`Timestep` to a trajectory
file.

.. autoclass:: WriterBase
   :members:
   :inherited-members:

Converters
----------
Converters output information to other libraries.

.. deprecated:: 2.7.0
    All converter code has been moved to :mod:`MDAnalysis.converters` and will
    be removed from the :mod:`MDAnalysis.coordinates.base` module in 3.0.0.

.. autoclass:: ConverterBase
   :members:
   :inherited-members:

Helper classes
--------------

The following classes contain basic functionality that all readers and
writers share.

.. autoclass:: IOBase
   :members:

"""
import abc
import numpy as np
import numbers
import warnings
from typing import Any, Union, Optional, List, Dict

from .timestep import Timestep
from . import core
from .. import (
    _READERS, _READER_HINTS,
    _SINGLEFRAME_WRITERS,
    _MULTIFRAME_WRITERS,
    _CONVERTERS,  # remove in 3.0.0 (Issue #3404)
)
from .. import units
from ..auxiliary.base import AuxReader
from ..auxiliary.core import auxreader
from ..auxiliary.core import get_auxreader_for
from ..auxiliary import _AUXREADERS
from ..lib.util import asiterable, Namespace, store_init_arguments
from ..lib.util import NamedStream


class FrameIteratorBase(object):
    """
    Base iterable over the frames of a trajectory.

    A frame iterable has a length that can be accessed with the :func:`len`
    function, and can be indexed similarly to a full trajectory. When indexed,
    indices are resolved relative to the iterable and not relative to the
    trajectory.

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
        return range_length(self.start, self.stop, self.step)

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
            step = (frame.step or 1) * self.step
            if frame.start is None:
                if frame.step is None or frame.step > 0:
                    start = self.start
                else:
                    start = self.start + (len(self) - 1) * self.step
            else:
                start = self.start + (frame.start or 0) * self.step
            if frame.stop is None:
                if frame.step is None or frame.step > 0:
                    last = start + (range_length(start, self.stop, step) - 1) * step
                else:
                    last = self.start
                stop = last + np.sign(step)
            else:
                stop = self.start + (frame.stop or 0) * self.step

            new_slice = slice(start, stop, step)
            frame_iterator = FrameIteratorSliced(self.trajectory, new_slice)
            # The __init__ of FrameIteratorSliced does some conversion between
            # the way indices are handled in slices and the way they are
            # handled by range. We need to overwrite this conversion as we
            # already use the logic for range.
            frame_iterator._start = start
            frame_iterator._stop = stop
            frame_iterator._step = step
            return frame_iterator
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
        self.trajectory.rewind()

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
        f = units.get_conversion_factor('length',
                                        self.units['length'], 'Angstrom')
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
            'speed', self.units['velocity'], 'Angstrom/ps')
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
            'force', self.units['force'], 'kJ/(mol*Angstrom)')
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
            'time', self.units['time'], 'ps')
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
            'length', 'Angstrom', self.units['length'])
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
            'speed', 'Angstrom/ps', self.units['velocity'])
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
            'force', 'kJ/(mol*Angstrom)', self.units['force'])
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
            'time', 'ps', self.units['time'])
        if f == 1.:
            return t
        if not inplace:
            return f * t
        t *= f
        return t

    def close(self):
        """Close the trajectory file."""
        pass # pylint: disable=unnecessary-pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # see http://docs.python.org/2/library/stdtypes.html#typecontextmanager
        self.close()
        return False  # do not suppress exceptions


class _Readermeta(abc.ABCMeta):
    """Automatic Reader registration metaclass

    .. versionchanged:: 1.0.0
       Added _format_hint functionality
    """
    # Auto register upon class creation
    def __init__(cls, name, bases, classdict):
        type.__init__(type, name, bases, classdict)  # pylint: disable=non-parent-init-called
        try:
            fmt = asiterable(classdict['format'])
        except KeyError:
            pass
        else:
            for fmt_name in fmt:
                fmt_name = fmt_name.upper()
                _READERS[fmt_name] = cls

                if '_format_hint' in classdict:
                    # isn't bound yet, so access __func__
                    _READER_HINTS[fmt_name] = classdict['_format_hint'].__func__


class ProtoReader(IOBase, metaclass=_Readermeta):
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
    .. versionchanged:: 2.0.0
       Now supports (un)pickle. Upon unpickling,
       the current timestep is retained by reconstrunction.
    """

    #: The appropriate Timestep class, e.g.
    #: :class:`MDAnalysis.coordinates.xdrfile.XTC.Timestep` for XTC.
    _Timestep = Timestep
    _transformations: list
    _auxs: dict
    _filename: Any
    n_frames: int

    def __init__(self):
        # initialise list to store added auxiliary readers in
        # subclasses should now call super
        self._auxs = {}
        self._transformations=[]

    def __len__(self) -> int:
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

    def next(self) -> Timestep:
        """Forward one step to next frame."""
        try:
            ts = self._read_next_timestep()
        except (EOFError, IOError):
            self.rewind()
            raise StopIteration from None
        else:
            for auxname, reader in self._auxs.items():
                ts = self._auxs[auxname].update_ts(ts)

            ts = self._apply_transformations(ts)

        return ts

    def __next__(self) -> Timestep:
        """Forward one step to next frame when using the `next` builtin."""
        return self.next()

    def rewind(self) -> Timestep:
        """Position at beginning of trajectory"""
        self._reopen()
        self.next()

    @property
    def dt(self) -> float:
        """Time between two trajectory frames in picoseconds."""
        return self.ts.dt

    @property
    def totaltime(self) -> float:
        """Total length of the trajectory

        The time is calculated as ``(n_frames - 1) * dt``, i.e., we assume that
        the first frame no time as elapsed. Thus, a trajectory with two frames will
        be considered to have a length of a single time step `dt` and a
        "trajectory" with a single frame will be reported as length 0.

        """
        return (self.n_frames - 1) * self.dt

    @property
    def frame(self) -> int:
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

    @abc.abstractmethod
    def _read_next_timestep(self, ts=None):
        # Example from DCDReader:
        #     if ts is None:
        #         ts = self.ts
        #     ts.frame = self._read_next_frame(etc)
        #     return ts
        ...

    def __iter__(self):
        """ Iterate over trajectory frames. """
        self._reopen()
        return self

    @abc.abstractmethod
    def _reopen(self):
        """Should position Reader to just before first frame

        Calling next after this should return the first frame
        """
        pass # pylint: disable=unnecessary-pass

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
                frame = np.asarray(frame, dtype=bool)
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
        #                                  ts.dimensions, 1)
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
            errmsg = f"{self.__class__.__name__} does not support slicing."
            raise TypeError(errmsg) from None

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

    def timeseries(self, asel: Optional['AtomGroup']=None,
                   atomgroup: Optional['Atomgroup']=None,
                   start: Optional[int]=None, stop: Optional[int]=None,
                   step: Optional[int]=None,
                   order: Optional[str]='fac') -> np.ndarray:
        """Return a subset of coordinate data for an AtomGroup

        Parameters
        ----------
        asel : AtomGroup (optional)
            The :class:`~MDAnalysis.core.groups.AtomGroup` to read the
            coordinates from. Defaults to ``None``, in which case the full set
            of coordinate data is returned.

            .. deprecated:: 2.7.0
                asel argument will be renamed to atomgroup in 3.0.0

        atomgroup: AtomGroup (optional)
            Same as `asel`, will replace `asel` in 3.0.0
        start :  int (optional)
            Begin reading the trajectory at frame index `start` (where 0 is the
            index of the first frame in the trajectory); the default
            ``None`` starts at the beginning.
        stop : int (optional)
            End reading the trajectory at frame index `stop`-1, i.e, `stop` is
            excluded. The trajectory is read to the end with the default
            ``None``.
        step : int (optional)
            Step size for reading; the default ``None`` is equivalent to 1 and
            means to read every frame.
        order : str (optional)
            the order/shape of the return data array, corresponding
            to (a)tom, (f)rame, (c)oordinates all six combinations
            of 'a', 'f', 'c' are allowed ie "fac" - return array
            where the shape is (frame, number of atoms,
            coordinates)

        See Also
        --------
        :class:`MDAnalysis.coordinates.memory`


        .. versionadded:: 2.4.0
        """
        if asel is not None:
            warnings.warn(
                "asel argument to timeseries will be renamed to"
                "'atomgroup' in 3.0, see #3911",
                category=DeprecationWarning)
            if atomgroup:
                raise ValueError("Cannot provide both asel and atomgroup kwargs")
            atomgroup = asel
        start, stop, step = self.check_slice_indices(start, stop, step)
        nframes = len(range(start, stop, step))

        if atomgroup is not None:
            if len(atomgroup) == 0:
                raise ValueError(
                    "Timeseries requires at least one atom to analyze")
            atom_numbers = atomgroup.indices
            natoms = len(atom_numbers)
        else:
            natoms = self.n_atoms
            atom_numbers = np.arange(natoms)

        # allocate output array in 'fac' order
        coordinates = np.empty((nframes, natoms, 3), dtype=np.float32)
        for i, ts in enumerate(self[start:stop:step]):
            coordinates[i, :] = ts.positions[atom_numbers]

        # switch axes around
        default_order = 'fac'
        if order != default_order:
            try:
                newidx = [default_order.index(i) for i in order]
            except ValueError:
                raise ValueError(f"Unrecognized order key in {order}, "
                                 "must be permutation of 'fac'")

            try:
                coordinates = np.moveaxis(coordinates, newidx, [0, 1, 2])
            except ValueError:
                errmsg = ("Repeated or missing keys passed to argument "
                          f"`order`: {order}, each key must be used once")
                raise ValueError(errmsg)
        return coordinates

# TODO: Change order of aux_spec and auxdata for 3.0 release, cf. Issue #3811
    def add_auxiliary(self,
                      aux_spec: Union[str, Dict[str, str]] = None,
                      auxdata: Union[str, AuxReader] = None,
                      format: str = None,
                      **kwargs) -> None:
        """Add auxiliary data to be read alongside trajectory.

        Auxiliary data may be any data timeseries from the trajectory
        additional to that read in by the trajectory reader. *auxdata* can
        be an :class:`~MDAnalysis.auxiliary.base.AuxReader` instance, or the
        data itself as e.g. a filename; in the latter case an appropriate
        :class:`~MDAnalysis.auxiliary.base.AuxReader` is guessed from the
        data/file format. An appropriate `format` may also be directly provided
        as a key word argument.

        On adding, the AuxReader is initially matched to the current timestep
        of the trajectory, and will be updated when the trajectory timestep
        changes (through a call to :meth:`next()` or jumping timesteps with
        ``trajectory[i]``).

        The representative value(s) of the auxiliary data for each timestep (as
        calculated by the :class:`~MDAnalysis.auxiliary.base.AuxReader`) are
        stored in the current timestep in the ``ts.aux`` namespace under
        *aux_spec*; e.g. to add additional pull force data stored in
        pull-force.xvg::

            u = MDAnalysis.Universe(PDB, XTC)
            u.trajectory.add_auxiliary('pull', 'pull-force.xvg')

        The representative value for the current timestep may then be accessed
        as ``u.trajectory.ts.aux.pull`` or ``u.trajectory.ts.aux['pull']``.


        The following applies to energy readers like the
        :class:`~MDAnalysis.auxiliary.EDR.EDRReader`.

        All data that is present in the (energy) file can be added by omitting
        `aux_spec` like so::

            u.trajectory.add_auxiliary(auxdata="ener.edr")

        *aux_spec* is expected to be a dictionary that maps the desired
        attribute name in the ``ts.aux`` namespace to the precise data to be
        added as identified by a :attr:`data_selector`::

            term_dict = {"temp": "Temperature", "epot": "Potential"}
            u.trajectory.add_auxiliary(term_dict, "ener.edr")

        Adding this data can be useful, for example, to filter trajectory
        frames based on non-coordinate data like the potential energy of each
        time step. Trajectory slicing allows working on a subset of frames::

            selected_frames = np.array([ts.frame for ts in u.trajectory
                                        if ts.aux.epot < some_threshold])
            subset = u.trajectory[selected_frames]


        See Also
        --------
        :meth:`remove_auxiliary`

        Note
        ----
        Auxiliary data is assumed to be time-ordered, with no duplicates. See
        the :ref:`Auxiliary API`.
        """
        if auxdata is None:
            raise ValueError("No input `auxdata` specified, but it needs "
                             "to be provided.")
        if type(auxdata) not in list(_AUXREADERS.values()):
            # i.e. if auxdata is a file, not an instance of an AuxReader
            reader_type = get_auxreader_for(auxdata)
            auxreader = reader_type(auxdata)
        else:
            auxreader = auxdata
        auxreader.attach_auxiliary(self, aux_spec, format, **kwargs)

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
        """Add all transformations to be applied to the trajectory.

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

        The transformations are applied in the order given in the list
        `transformations`, i.e., the first transformation is the first
        or innermost one to be applied to the :class:`Timestep`. The
        example above would be equivalent to

        .. code-block:: python

          for ts in u.trajectory:
             ts = this_transform(another_transform(some_transform(ts)))


        Parameters
        ----------
        transform_list : list
            list of all the transformations that will be applied to the coordinates
            in the order given in the list

        See Also
        --------
        :mod:`MDAnalysis.transformations`

        """

        try:
            self.transformations = transformations
        except ValueError:
            errmsg = ("Can't add transformations again. Please create a new "
                      "Universe object")
            raise ValueError(errmsg) from None
        else:
            self.ts = self._apply_transformations(self.ts)


        # call reader here to apply the newly added transformation on the
        # current loaded frame?

    def _apply_transformations(self, ts):
        """Applies all the transformations given by the user """

        for transform in self.transformations:
            ts = transform(ts)

        return ts

    def __setstate__(self, state):
        self.__dict__ = state
        self[self.ts.frame]


class ReaderBase(ProtoReader):
    """Base class for trajectory readers that extends :class:`ProtoReader` with a
    :meth:`__del__` method.

    New Readers should subclass :class:`ReaderBase` and properly implement a
    :meth:`close` method, to ensure proper release of resources (mainly file
    handles). Readers that are inherently safe in this regard should subclass
    :class:`ProtoReader` instead.

    See the :ref:`Trajectory API` definition in for the required attributes and
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
    .. versionchanged:: 1.0
       Removed deprecated flags functionality, use convert_units kwarg instead

    """
    @store_init_arguments
    def __init__(self, filename, convert_units=True, **kwargs):
        super(ReaderBase, self).__init__()

        if isinstance(filename, NamedStream):
            self.filename = filename
        else:
            self.filename = str(filename)
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

        Will also copy the current state of the Timestep held in the original
        Reader.


        .. versionchanged:: 2.2.0
           Arguments used to construct the reader are correctly captured and
           passed to the creation of the new class. Previously the only
           ``n_atoms`` was passed to class copies, leading to a class created
           with default parameters which may differ from the original class.
        """

        new = self.__class__(**self._kwargs)

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


class WriterBase(IOBase, metaclass=_Writermeta):
    """Base class for trajectory writers.

    See :ref:`Trajectory API` definition in for the required attributes and
    methods.


    .. versionchanged:: 2.0.0
       Deprecated :func:`write_next_timestep` has now been removed, please use
       :func:`write` instead.

    """

    def convert_dimensions_to_unitcell(self, ts, inplace=True):
        """Read dimensions from timestep *ts* and return appropriate unitcell.

        The default is to return ``[A,B,C,alpha,beta,gamma]``; if this
        is not appropriate then this method has to be overriden.
        """
        # override if the native trajectory format does NOT use
        # [A,B,C,alpha,beta,gamma]
        if ts.dimensions is None:
            lengths, angles = np.zeros(3), np.zeros(3)
        else:
            lengths, angles = ts.dimensions[:3], ts.dimensions[3:]
        if not inplace:
            lengths = lengths.copy()
        lengths = self.convert_pos_to_native(lengths)
        return np.concatenate([lengths, angles])

    def write(self, obj):
        """Write current timestep, using the supplied `obj`.

        Parameters
        ----------
        obj : :class:`~MDAnalysis.core.groups.AtomGroup` or :class:`~MDAnalysis.core.universe.Universe`
            write coordinate information associate with `obj`

        Note
        ----
        The size of the `obj` must be the same as the number of atoms provided
        when setting up the trajectory.


        .. versionchanged:: 2.0.0
           Deprecated support for Timestep argument to write has now been
           removed. Use AtomGroup or Universe as an input instead.
        """
        return self._write_next_frame(obj)

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

        Due to rounding, the test is asymmetric (and *min* is supposed to be negative)::

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


class SingleFrameReaderBase(ProtoReader):
    """Base class for Readers that only have one frame.

    To use this base class, define the method :meth:`_read_first_frame` to
    read from file `self.filename`.  This should populate the attribute
    `self.ts` with a :class:`Timestep` object.

    .. versionadded:: 0.10.0
    .. versionchanged:: 0.11.0
       Added attribute "_ts_kwargs" for subclasses
       Keywords "dt" and "time_offset" read into _ts_kwargs
    .. versionchanged:: 2.2.0
       Calling `__iter__` now rewinds the reader before yielding a
       :class:`Timestep` object (fixing behavior that was not
       well defined previously).
    """
    _err = "{0} only contains a single frame"

    @store_init_arguments
    def __init__(self, filename, convert_units=True, n_atoms=None, **kwargs):
        super(SingleFrameReaderBase, self).__init__()

        if isinstance(filename, NamedStream):
            self.filename = filename
        else:
            self.filename = str(filename)

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

        Will also copy the current state of the Timestep held in the original
        Reader.


        .. versionchanged:: 2.2.0
           Arguments used to construct the reader are correctly captured and
           passed to the creation of the new class. Previously the only
           ``n_atoms`` was passed to class copies, leading to a class created
           with default parameters which may differ from the original class.
        """
        new = self.__class__(**self._kwargs)

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

    def _read_next_timestep(self, ts=None):
        raise NotImplementedError(self._err.format(self.__class__.__name__))

    def __iter__(self):
        self.rewind()
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


def range_length(start, stop, step):
    if (step > 0 and start < stop):
        # We go from a lesser number to a larger one.
        return int(1 + (stop - 1 - start) // step)
    elif (step < 0 and start > stop):
        # We count backward from a larger number to a lesser one.
        return int(1 + (start - 1 - stop) // (-step))
    else:
        # The range is empty.
        return 0

# Verbatim copy of code from converters/base.py
# Needed to avoid circular imports before removal in
# MDAnalysis 3.0.0
# Remove in 3.0.0
class _Convertermeta(type):
    # Auto register upon class creation
    def __init__(cls, name, bases, classdict):
        type.__init__(type, name, bases, classdict)
        try:
            fmt = asiterable(classdict['lib'])
        except KeyError:
            pass
        else:
            for f in fmt:
                f = f.upper()
                _CONVERTERS[f] = cls


# Verbatim copy of code from converters/base.py
# Needed to avoid circular imports before removal in
# MDAnalysis 3.0.0
# Remove in 3.0.0
class ConverterBase(IOBase, metaclass=_Convertermeta):
    """Base class for converting to other libraries.

    .. deprecated:: 2.7.0
        This class has been moved to
        :class:`MDAnalysis.converters.base.ConverterBase` and will be removed
        from :mod:`MDAnalysis.coordinates.base` in 3.0.0.
    """

    def __init_subclass__(cls):
        wmsg = ("ConverterBase moved from coordinates.base."
                "ConverterBase to converters.base.ConverterBase "
                "and will be removed from coordinates.base "
                "in MDAnalysis release 3.0.0")
        warnings.warn(wmsg, DeprecationWarning, stacklevel=2)

    def __repr__(self):
        return "<{cls}>".format(cls=self.__class__.__name__)

    def convert(self, obj):
        raise NotImplementedError
