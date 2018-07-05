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
ChainReader --- :mod:`MDAnalysis.coordinates.chain`
===================================================

The :class:`ChainReader` is used by MDAnalysis internally to represent multiple
trajectories as one virtual trajectory. Users typically do not need to use the
:class:`ChainReader` explicitly and the following documentation is primarily of
interest to developers.

.. autoclass:: ChainReader
   :members:

   .. automethod:: _get_local_frame
   .. automethod:: _apply
   .. automethod:: _get
   .. automethod:: _get_same
   .. automethod:: _read_frame
   .. automethod:: _chained_iterator

"""
from __future__ import absolute_import

import warnings

import os.path
import itertools
import bisect
import copy

import numpy as np

from ..lib.util import asiterable
from . import base
from . import core


def multi_level_argsort(l):
    """Return indices to sort a multi value tuple. Sorting is done on the first
    value of the tuple.


    Parameters
    ----------
    l : list

    Returns
    -------
    indices

    Example
    -------
    >>> multi_level_argsort(((0, 2), (4, 9), (0, 4), (7, 9)))
    [0, 2, 1, 3]

    """
    return [el[0] for el in sorted(enumerate(l), key=lambda x: x[1][0])]


def filter_times(times, dt):
    """Given a list of start and end times this function filters out any duplicate
    time steps preferring the last tuple.

    Parameters
    ----------
    times : list
        sorted list of times
    dt : float
        timestep between two frames

    Returns
    -------
    list:
        indices of times to be used with overlaps removed

    Example
    -------
    >>> filter_times(((0, 3), (0, 3)))
    [1, ]
    >>> filter_times(((0, 3), (0, 4)))
    [1, ]
    >>> filter_times(((0, 3), (3, 4)))
    [0, 1]
    >>> filter_times(((0, 3), (2, 5), (4, 9)))
    [1, 2, 3]

    """
    # Special cases
    if len(times) == 1:
        return [0, ]
    elif len(times) == 2:
        if times[0][0] < times[1][0]:
            return [0, 1]
        elif np.allclose(times[0][0], times[1][0]):
            return [1, ]
        else:
            return [0, ]
    if np.unique(times).size == 2:
        return [len(times) - 1, ]

    # more then 2 unique time entries

    used_idx = [0, ]

    for i, (first, middle, last) in enumerate(zip(times[:-2], times[1:-1], times[2:]), start=1):
        if np.allclose(first[0], middle[0]):
            used_idx[-1] = i
        elif not np.allclose(middle[1] - middle[0], dt):
            if (middle[0] <= first[1]) and (last[0] <= middle[1]):
                used_idx.append(i)
        elif (middle[0] <= first[1]):
            used_idx.append(i)

    # take care of first special case
    if times[-2][1] <= times[-1][1]:
        used_idx.append(len(times) - 1)

    return used_idx


class ChainReader(base.ProtoReader):
    """Reader that concatenates multiple trajectories on the fly.

    The :class:`ChainReader` is used by MDAnalysis internally to
    represent multiple trajectories as one virtual trajectory. Users
    typically do not need to use the :class:`ChainReader` explicitly.

    Chainreader can also handle a continuous trajectory split over several
    files. To use this pass the ``continuous == True`` keyword argument.
    Setting ``continuous=True`` will make the reader choose frames from the set
    of trajectories in such a way that the trajectory appears to be as
    continuous in time as possible, i.e. that time is strictly monotonically
    increasing. This means that there will be no duplicate time frames and no
    jumps backwards in time. However, there can be gaps in time (e.g., multiple
    time steps can appear to be missing). Ultimately, it is the user's
    responsibility to ensure that the input trajectories can be virtually
    stitched together in a meaningful manner. As an example take the following
    trajectory that is split into three parts. The column represents the time
    and the trajectory segments overlap. With the continuous chainreader only
    the frames marked with a + will be read.

    ::

        part01:  ++++--
        part02:      ++++++-
        part03:            ++++++++

    .. warning::

        The order in which trajectories are given to the chainreader can change
        what frames are used with the continuous option.

    The default chainreader will read all frames. The continuous option is
    currently only supported for XTC and TRR files.

    Notes
    -----
    The trajectory API attributes exist but most of them only reflect the first
    trajectory in the list; :attr:`ChainReader.n_frames`,
    :attr:`ChainReader.n_atoms`, and :attr:`ChainReader.fixed` are properly
    set, though


    .. versionchanged:: 0.11.0
       Frames now 0-based instead of 1-based
    .. versionchanged:: 0.13.0
       :attr:`time` now reports the time summed over each trajectory's
       frames and individual :attr:`dt`.
    .. versionchanged:: 0.19.0
       added ``continuous`` trajectory option

    """
    format = 'CHAIN'

    def __init__(self, filenames, skip=1, dt=None, continuous=False, **kwargs):
        """Set up the chain reader.

        Parameters
        ----------
        filenames : str or list or sequence
            file name or list of file names; the reader will open all file names
            and provide frames in the order of trajectories from the list. Each
            trajectory must contain the same number of atoms in the same order
            (i.e. they all must belong to the same topology). The trajectory
            format is deduced from the extension of each file name.

            Extension: `filenames` are either a single file name or list of file
            names in either plain file names format or ``(filename, format)``
            tuple combination. This allows explicit setting of the format for
            each individual trajectory file.
        skip : int (optional)
            skip step (also passed on to the individual trajectory readers);
            must be same for all trajectories
        dt : float (optional)
            Passed to individual trajectory readers to enforce a common time
            difference between frames, in MDAnalysis time units. If not set, each
            reader's `dt` will be used (either inferred from the trajectory
            files, or set to the reader's default) when reporting frame times;
            note that this might lead an inconsistent time difference between
            frames.
        continuous : bool (optional)
            treat all trajectories as one single long trajectory. Adds several
            checks; all trajectories have the same dt, they contain at least 2
            frames, and they are all of the same file-type. Not implemented for
            all trajectory formats! This can be used to analyze GROMACS
            simulations without concatenating them prior to analysis.
        **kwargs : dict (optional)
            all other keyword arguments are passed on to each trajectory reader
            unchanged

        """
        super(ChainReader, self).__init__()

        filenames = asiterable(filenames)
        # Override here because single frame readers handle this argument as a
        # kwarg to a timestep which behaves differently if dt is present or not.
        if dt is not None:
            kwargs['dt'] = dt
        self.readers = [core.reader(filename, **kwargs)
                        for filename in filenames]
        self.filenames = np.array([fn[0] if isinstance(fn, tuple) else fn for fn in filenames])
        # pointer to "active" trajectory index into self.readers
        self.__active_reader_index = 0

        self.skip = skip
        self.n_atoms = self._get_same('n_atoms')

        # Translation between virtual frames and frames in individual
        # trajectories. Assumes that individual trajectories i contain frames
        # that can be addressed with an index 0 <= f < n_frames[i]

        # Build a map of frames: ordered list of starting virtual frames; the
        # index i into this list corresponds to the index into self.readers
        #
        # For virtual frame 0 <= k < sum(n_frames) find corresponding
        # trajectory i and local frame f (i.e. readers[i][f] will correspond to
        # ChainReader[k]).
        # build map 'start_frames', which is used by _get_local_frame()
        n_frames = self._get('n_frames')
        # [0]: frames are 0-indexed internally
        # (see Timestep.check_slice_indices())
        self._start_frames = np.cumsum([0] + n_frames)

        self.n_frames = np.sum(n_frames)
        self.dts = np.array(self._get('dt'))
        self.total_times = self.dts * n_frames

        #: source for trajectories frame (fakes trajectory)
        self.__chained_trajectories_iter = None

        # calculate new start_frames to have a time continuous trajectory.
        if continuous:
            filetypes = np.unique([r.format for r in self.readers])
            if not len(filetypes) == 1:
                raise ValueError("ChainReader: continuous=true only supported"
                                 " when all files are using the same format. "
                                 "found {}".format(filetypes))
            if np.any(np.array(n_frames) == 1):
                raise RuntimeError("ChainReader: Need at least two frames in "
                                   "every trajectory with continuous=True")
            if filetypes[0] not in ['XTC', 'TRR']:
                raise NotImplementedError("ChainReader: continuous=True only "
                                          "supported for xtc and trr format")

            # TODO: allow floating point precision in dt check
            dt = self._get_same('dt')
            n_frames = np.asarray(self._get('n_frames'))
            self.dts = np.ones(self.dts.shape) * dt

            # the sorting needs to happen on two levels. The first major level
            # is by start times and the second is by end times.
            # The second level of sorting is needed for cases like:
            # [0 1 2 3 4 5 6 7 8 9] [0 1 2 4]
            # to
            # [0 1 2 4] [0 1 2 3 4 5 6 7 8 9]
            # after that sort the chain reader will work
            times = []
            for r in self.readers:
                r[0]
                start = r.ts.time
                r[-1]
                end = r.ts.time
                times.append((start, end))
            # sort step
            sort_idx = multi_level_argsort(times)
            self.readers = [self.readers[i] for i in sort_idx]
            self.filenames = self.filenames[sort_idx]
            self.total_times = self.dts * n_frames[sort_idx]

            # filter step: remove indices if we have complete overlap
            if len(self.readers) > 1:
                used_idx = filter_times(np.array(times)[sort_idx], dt)

                self.readers = [self.readers[i] for i in used_idx]
                self.filenames = self.filenames[used_idx]
                self.total_times = self.dts[used_idx] * n_frames[used_idx]

            # rebuild lookup table
            sf = [0, ]
            n_frames = 0
            for r1, r2 in zip(self.readers[:-1], self.readers[1:]):
                r2[0], r1[0]
                r1_start_time = r1.time
                start_time = r2.time
                r1[-1]
                if r1.time < start_time:
                    warnings.warn("Missing frame in continuous chain", UserWarning)

                # check for interleaving
                r1[1]
                if r1_start_time < start_time < r1.time:
                    raise RuntimeError("ChainReader: Interleaving not supported with continuous=True.")

                # find end where trajectory was restarted from
                for ts in r1[::-1]:
                    if ts.time < start_time:
                        break
                sf.append(sf[-1] + ts.frame + 1)
                n_frames += ts.frame + 1

            n_frames += self.readers[-1].n_frames

            self._start_frames = sf
            self.n_frames = n_frames
            self._sf = sf

        # make sure that iteration always yields frame 0
        # rewind() also sets self.ts
        self.ts = None
        self.rewind()

    def _get_local_frame(self, k):
        """Find trajectory index and trajectory frame for chained frame `k`.

        Parameters
        ----------

        k : int
           Frame `k` in the chained trajectory can be found in the trajectory at
           index *i* and frame index *f*.

           Frames are internally treated as 0-based indices into the trajectory.

        Returns
        -------
        i : int
            trajectory
        f : int
            frame in trajectory i


        Raises
        ------
        IndexError for `k<0` or `i<0`.


        Note
        ----
        Does not check if `k` is larger than the maximum number of frames in
        the chained trajectory.

        """
        if k < 0:
            raise IndexError("Virtual (chained) frames must be >= 0")
        # trajectory index i
        i = bisect.bisect_right(self._start_frames, k) - 1
        if i < 0:
            raise IndexError("Cannot find trajectory for virtual frame {0:d}".format(k))
        # local frame index f in trajectory i (frame indices are 0-based)
        f = k - self._start_frames[i]
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

    def copy(self):
        new = self.__class__(copy.copy(self.filenames))
        # seek the new reader to the same frame we started with
        new[self.ts.frame]
        # then copy over the current Timestep in case it has
        # been modified since initial load
        new.ts = self.ts.copy()
        return new



    # attributes that can change with the current reader
    @property
    def filename(self):
        """Filename of the currently read trajectory"""
        return self.active_reader.filename

    # TODO: check that skip_timestep is still supported in all readers
    #       or should this be removed?
    @property
    def skip_timestep(self):
        return self.active_reader.skip_timestep

    @property
    def delta(self):
        return self.active_reader.delta

    @property
    def periodic(self):
        """:attr:`periodic` attribute of the currently read trajectory"""
        return self.active_reader.periodic

    @property
    def units(self):
        """:attr:`units` attribute of the currently read trajectory"""
        return self.active_reader.units

    @property
    def compressed(self):
        """:attr:`compressed` attribute of the currently read trajectory"""
        try:
            return self.active_reader.compressed
        except AttributeError:
            return None

    @property
    def frame(self):
        """Cumulative frame number of the current time step."""
        return self.ts.frame

    @property
    def time(self):
        """Cumulative time of the current frame in MDAnalysis time units (typically ps)."""
        # Before 0.13 we had to distinguish between enforcing a common dt or
        #  summing over each reader's times.
        # Now each reader is either already instantiated with a common dt, or
        #  left at its default dt. In any case, we sum over individual times.
        trajindex, subframe = self._get_local_frame(self.frame)
        return self.total_times[:trajindex].sum() + subframe * self.dts[trajindex]

    def _apply(self, method, **kwargs):
        """Execute `method` with `kwargs` for all readers."""
        return [reader.__getattribute__(method)(**kwargs) for reader in self.readers]

    def _get(self, attr):
        """Get value of `attr` for all readers."""
        return [reader.__getattribute__(attr) for reader in self.readers]

    def _get_same(self, attr):
        """Verify that `attr` has the same value for all readers and return value.

        Parameters
        ----------
        attr : str
            attribute name

        Returns
        -------
        value : int or float or str or object
            common value of the attribute

        Raises
        ------
        ValueError if not all readers have the same value

        """
        values = np.array(self._get(attr))
        value = values[0]
        if not np.allclose(values, value):
            bad_traj = np.array(self.filenames)[values != value]
            raise ValueError("The following trajectories do not have the correct {0} "
                             " ({1}):\n{2}".format(attr, value, bad_traj))
        return value

    def __activate_reader(self, i):
        """Make reader `i` the active reader."""
        # private method, not to be used by user to avoid a total mess
        if not (0 <= i < len(self.readers)):
            raise IndexError("Reader index must be 0 <= i < {0:d}".format(len(self.readers)))
        self.__active_reader_index = i

    @property
    def active_reader(self):
        """Reader instance from which frames are currently being read."""
        return self.readers[self.__active_reader_index]

    def _read_frame(self, frame):
        """Position trajectory at frame index `frame` and
        return :class:`~MDAnalysis.coordinates.base.Timestep`.

        The frame is translated to the corresponding reader and local
        frame index and the :class:`Timestep` instance in
        :attr:`ChainReader.ts` is updated.

        Notes
        -----
        `frame` is 0-based, i.e. the first frame in the trajectory is
        accessed with ``frame = 0``.

        See Also
        --------
        :meth:`~ChainReader._get_local_frame`

        """
        i, f = self._get_local_frame(frame)
        # seek to (1) reader i and (2) frame f in trajectory i
        self.__activate_reader(i)
        self.active_reader[f]  # rely on reader to implement __getitem__()
        # update Timestep
        self.ts = self.active_reader.ts
        self.ts.frame = frame  # continuous frames, 0-based
        return self.ts

    def _chained_iterator(self):
        """Iterator that presents itself as a chained trajectory."""
        self._rewind()  # must rewind all readers
        for i in range(self.n_frames):
            j, f = self._get_local_frame(i)
            self.__activate_reader(j)
            self.ts = self.active_reader[f]
            self.ts.frame = i
            yield self.ts

    def _read_next_timestep(self, ts=None):
        self.ts = next(self.__chained_trajectories_iter)
        return self.ts

    def rewind(self):
        """Set current frame to the beginning."""
        self._rewind()
        self.__chained_trajectories_iter = self._chained_iterator()
        # set time step for frame 1
        self.ts = next(self.__chained_trajectories_iter)

    def _rewind(self):
        """Internal method: Rewind trajectories themselves and trj pointer."""
        self._apply('rewind')
        self.__activate_reader(0)

    def close(self):
        self._apply('close')

    def __iter__(self):
        """Generator for all frames, starting at frame 1."""
        self._rewind()
        # start from first frame
        self.__chained_trajectories_iter = self._chained_iterator()
        for ts in self.__chained_trajectories_iter:
            yield ts

    def __repr__(self):
        return ("<{clsname} {fname} with {nframes} frames of {natoms} atoms>"
                "".format(
                    clsname=self.__class__.__name__,
                    fname=[os.path.basename(fn) for fn in self.filenames],
                    nframes=self.n_frames,
                    natoms=self.n_atoms))

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

        super(ChainReader, self).add_transformations(*transformations)
        for r in self.readers:
            r.add_transformations(*transformations)

    def _apply_transformations(self, ts):
        """ Applies the transformations to the timestep."""
        # Overrides :meth:`~MDAnalysis.coordinates.base.ProtoReader.add_transformations`
        # to avoid applying the same transformations multiple times on each frame

        return ts
