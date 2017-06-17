# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.mdanalysis.org
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

import numpy as np

from ..lib.util import asiterable
from . import base
from . import core


class ChainReader(base.ProtoReader):
    """Reader that concatenates multiple trajectories on the fly.

    The :class:`ChainReader` is used by MDAnalysis internally to
    represent multiple trajectories as one virtual trajectory. Users
    typically do not need to use the :class:`ChainReader` explicitly.

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

    """
    format = 'CHAIN'

    def __init__(self, filenames, **kwargs):
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
        **kwargs : dict (optional)
          all other keyword arguments are passed on to each trajectory reader
          unchanged


        .. versionchanged:: 0.8
           The ``delta`` keyword was added.
        .. versionchanged:: 0.13
           The ``delta`` keyword was deprecated in favor of using ``dt``.

        """
        super(ChainReader, self).__init__()

        self.filenames = asiterable(filenames)
        self.readers = [core.reader(filename, **kwargs)
                        for filename in self.filenames]
        # pointer to "active" trajectory index into self.readers
        self.__active_reader_index = 0

        self.skip = kwargs.get('skip', 1)
        self.n_atoms = self._get_same('n_atoms')
        #self.fixed = self._get_same('fixed')

        # Translation between virtual frames and frames in individual
        # trajectories.
        # Assumes that individual trajectories i contain frames that can
        # be addressed with an index 0 <= f < n_frames[i]

        # Build a map of frames: ordered list of starting virtual
        # frames; the index i into this list corresponds to the index
        # into self.readers
        #
        # For virtual frame 0 <= k < sum(n_frames) find corresponding
        # trajectory i and local frame f (i.e. readers[i][f] will
        # correspond to ChainReader[k]).

        # build map 'start_frames', which is used by _get_local_frame()
        n_frames = self._get('n_frames')
        # [0]: frames are 0-indexed internally
        # (see Timestep.check_slice_indices())
        self.__start_frames = np.cumsum([0] + n_frames)

        self.n_frames = np.sum(n_frames)
        self.dts = np.array(self._get('dt'))
        self.total_times = self.dts * n_frames

        #: source for trajectories frame (fakes trajectory)
        self.__chained_trajectories_iter = None

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
        (i, f) : tuple
            **local frame** tuple

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
        i = bisect.bisect_right(self.__start_frames, k) - 1
        if i < 0:
            raise IndexError("Cannot find trajectory for virtual frame {0:d}".format(k))
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
        if not np.all(values == value):
            bad_traj = np.array([self._get_filename(fn) for fn in self.filenames])[values != value]
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
        readers = itertools.chain(*self.readers)
        for frame, ts in enumerate(readers):
            ts.frame = frame  # fake continuous frames, 0-based
            self.ts = ts
            # make sure that the active reader is in sync
            i, f = self._get_local_frame(frame)  # uses 0-based frames!
            self.__activate_reader(i)
            yield ts

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

    @staticmethod
    def _get_filename(filename):
        """retrieve the actual filename of the list element"""
        return filename[0] if isinstance(filename, tuple) else filename

    def __repr__(self):
        return ("<{clsname} {fname} with {nframes} frames of {natoms} atoms>"
                "".format(
                    clsname=self.__class__.__name__,
                    fname=[os.path.basename(self._get_filename(fn))
                           for fn in self.filenames],
                    nframes=self.n_frames,
                    natoms=self.n_atoms))


