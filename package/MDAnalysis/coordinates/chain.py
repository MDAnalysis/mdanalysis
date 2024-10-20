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

"""
import warnings

import os.path
import bisect
from typing import Tuple
import numpy as np

from ..lib import util
from ..lib.util import asiterable, store_init_arguments
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


def check_allowed_filetypes(readers, allowed):
    """
    Make a check that  all readers have the same filetype and are  of the
    allowed files types. Throws Exception on failure.

    Parameters
    ----------
    readers : list of MDA readers
    allowed : list of allowed formats
    """
    classname = type(readers[0])
    only_one_reader = np.all([isinstance(r,  classname) for r in readers])
    if not only_one_reader:
        readernames = [type(r) for r in readers]
        raise ValueError("ChainReader: continuous=true only supported"
                         " when all files are using the same reader. "
                         "Found: {}".format(readernames))
    if readers[0].format not in allowed:
        raise NotImplementedError("ChainReader: continuous=True only "
                                  "supported for formats: {}".format(allowed))


class ChainReader(base.ReaderBase):
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
    currently only supported for XTC, TRR, and LAMMPSDUMP files.

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
    .. versionchanged:: 0.19.0
       limit output of __repr__
    .. versionchanged:: 2.0.0
       Now ChainReader can be (un)pickled. Upon unpickling,
       current timestep is retained.

    """
    format = 'CHAIN'

    @store_init_arguments
    def __init__(self, filenames, skip=1, dt=None, continuous=False,
                 convert_units=True, **kwargs):
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
        super(ChainReader, self).__init__(filename='CHAIN',
                                          skip=skip,
                                          convert_units=convert_units,
                                          dt=dt,
                                          **kwargs)

        filenames = asiterable(filenames)
        # Override here because single frame readers handle this argument as a
        # kwarg to a timestep which behaves differently if dt is present or not.
        if dt is not None:
            kwargs['dt'] = dt
        self.readers = [core.reader(filename, convert_units=convert_units, **kwargs)
                        for filename in filenames]
        # Iterate through all filenames, appending NoneType None for ndarrays
        self.filenames = []
        for fn in filenames:
            if isinstance(fn, np.ndarray):
                self.filenames.append(None)
            elif isinstance(fn, tuple):
                self.filenames.append(fn[0])
            else:
                self.filenames.append(fn)
        self.filenames = np.array(self.filenames)

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

        # calculate new start_frames to have a time continuous trajectory.
        if continuous:
            check_allowed_filetypes(self.readers, ['XTC', 'TRR', 'LAMMPSDUMP',
                                                   'TRC'])
            if np.any(np.array(n_frames) == 1):
                raise RuntimeError("ChainReader: Need at least two frames in "
                                   "every trajectory with continuous=True")
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
                    raise RuntimeError("ChainReader: Interleaving not supported "
                                       "with continuous=True.")

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

    @staticmethod
    def _format_hint(thing):
        """Can ChainReader read the object *thing*

        .. versionadded:: 1.0.0
        """
        return (not isinstance(thing, np.ndarray) and
                util.iterable(thing) and
                not util.isstream(thing))

    def _get_local_frame(self, k) -> Tuple[int, int]:
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

    def __getstate__(self):
        state = self.__dict__.copy()
        #  save ts temporarily otherwise it will be changed during rewinding.
        state['ts'] = self.ts.__deepcopy__()

        #  the ts.frame of each reader is set to the chained frame index during
        #  iteration, thus we need to rewind the readers that have been used.
        #  PR #2723
        for reader in state['readers'][:self.__active_reader_index + 1]:
            reader.rewind()

        #  retrieve the current ts
        self.ts = state['ts']
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self.ts.frame = self.__current_frame

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
        self.filename = self.filenames[i]

    @property
    def active_reader(self):
        """Reader instance from which frames are currently being read."""
        return self.readers[self.__active_reader_index]

    def _read_frame(self, frame):
        """Position trajectory at frame index `frame` and
        return :class:`~MDAnalysis.coordinates.timestep.Timestep`.

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
        self.__current_frame = frame
        return self.ts

    def _read_next_timestep(self, ts=None):
        if ts:
            warnings.warn("ts argument to _read_next_timestep is deprecated as of 2.7.0 and will be removed in 3.0.0, see #3928")

        if ts is None:
            ts = self.ts

        if self.__current_frame < self.n_frames - 1:
            j, f = self._get_local_frame(self.__current_frame + 1)
            self.__activate_reader(j)
            self.ts = self.active_reader[f]
            self.ts.frame = self.__current_frame + 1
            self.__current_frame += 1
            return self.ts
        else:
            raise StopIteration()

    def _reopen(self):
        """Internal method: Rewind trajectories themselves and trj pointer."""
        self.__current_frame = -1
        self._apply('rewind')

    def close(self):
        self._apply('close')

    def __repr__(self):
        if len(self.filenames) > 3:
            fname = (os.path.basename(self.filenames[0])
                     if self.filenames[0] else "numpy.ndarray")
            fnames = "{fname} and {nfnames} more".format(
                    fname=fname,
                    nfnames=len(self.filenames) - 1)
        else:
            fnames = ", ".join([os.path.basename(fn) if fn else "numpy.ndarray"
                                for fn in self.filenames])
        return ("<{clsname} containing {fname} with {nframes} frames of {natoms} atoms>"
                "".format(
                    clsname=self.__class__.__name__,
                    fname=fnames,
                    nframes=self.n_frames,
                    natoms=self.n_atoms))
