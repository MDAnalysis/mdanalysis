"""

Example: Streaming an IMD v3 trajectory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To stream a trajectory from GROMACS or another simulation engine that supports 
IMD v3, ensure that the simulation engine is running and waiting for an IMD connection.

For example, in GROMACS, you can use ``gmx mdrun`` with the ``-imdwait`` flag
to ensure that GROMACS will wait for a client before starting the simulation.
In GROMACS, you will know that the simulation is ready and waiting for the
MDAnalysis IMDReader client when this line is printed to the terminal:

.. code-block:: none

    IMD: Will wait until I have a connection and IMD_GO orders.

Once the simulation is ready for a client connection, setup your :class:`Universe`
like this: ::

    import MDAnalysis as mda
    # Pass host and port of the listening GROMACS simulation
    # server as the trajectory argument
    u = mda.Universe("topology.tpr", "localhost:8888")

Classes
^^^^^^^

.. autoclass:: IMDReader
   :members:
   :inherited-members:

"""

import queue
from MDAnalysis.coordinates.base import (
    ReaderBase,
    FrameIteratorIndices,
    FrameIteratorAll,
    FrameIteratorSliced,
)
from MDAnalysis.coordinates import core
from MDAnalysis.lib.util import store_init_arguments

import imdclient
from imdclient.utils import parse_host_port
import numpy as np
import logging
import warnings
from typing import Optional
import numbers

logger = logging.getLogger("imdclient.IMDClient")


class IMDReader(ReaderBase):
    """
    Reader for IMD protocol packets.
    """

    format = "IMD"

    @store_init_arguments
    def __init__(
        self,
        filename,
        convert_units=True,
        n_atoms=None,
        **kwargs,
    ):
        """
        Parameters
        ----------
        filename : a string of the form "host:port" where host is the hostname
            or IP address of the listening GROMACS server and port
            is the port number.
        n_atoms : int (optional)
            number of atoms in the system. defaults to number of atoms
            in the topology. don't set this unless you know what you're doing.
        """
        self._init_scope = True
        self._reopen_called = False

        super(IMDReader, self).__init__(filename, **kwargs)

        logger.debug("IMDReader initializing")

        if n_atoms is None:
            raise ValueError("IMDReader: n_atoms must be specified")
        self.n_atoms = n_atoms

        host, port = parse_host_port(filename)

        # This starts the simulation
        self._imdclient = imdclient.IMDClient(host, port, n_atoms, **kwargs)

        imdsinfo = self._imdclient.get_imdsessioninfo()
        # NOTE: after testing phase, fail out on IMDv2

        self.ts = self._Timestep(
            self.n_atoms,
            positions=imdsinfo.positions,
            velocities=imdsinfo.velocities,
            forces=imdsinfo.forces,
            **self._ts_kwargs,
        )

        self._frame = -1

        try:
            self._read_next_timestep()
        except StopIteration:
            raise RuntimeError("IMDReader: No data found in stream")

    def _read_next_timestep(self):
        # No rewinding- to both load the first frame after  __init__
        # and access it again during iteration, we need to store first ts in mem
        if not self._init_scope and self._frame == -1:
            self._frame += 1
            # can't simply return the same ts again- transformations would be applied twice
            # instead, return the pre-transformed copy
            return self._first_ts

        return self._read_frame(self._frame + 1)

    def _read_frame(self, frame):

        try:
            imdf = self._imdclient.get_imdframe()
        except EOFError:
            # Not strictly necessary, but for clarity
            raise StopIteration

        self._frame = frame
        self._load_imdframe_into_ts(imdf)

        if self._init_scope:
            self._first_ts = self.ts.copy()
            self._init_scope = False

        logger.debug(f"IMDReader: Loaded frame {self._frame}")
        return self.ts

    def _load_imdframe_into_ts(self, imdf):
        self.ts.frame = self._frame
        if imdf.time is not None:
            self.ts.time = imdf.time
            # NOTE: timestep.pyx "dt" method is suspicious bc it uses "new" keyword for a float
            self.ts.data["dt"] = imdf.dt
        if imdf.energies is not None:
            self.ts.data.update(imdf.energies)
        if imdf.box is not None:
            self.ts.dimensions = core.triclinic_box(*imdf.box)
        if imdf.positions is not None:
            # must call copy because reference is expected to reset
            # see 'test_frame_collect_all_same' in MDAnalysisTests.coordinates.base
            self.ts.positions = imdf.positions
        if imdf.velocities is not None:
            self.ts.velocities = imdf.velocities
        if imdf.forces is not None:
            self.ts.forces = imdf.forces

    @property
    def n_frames(self):
        """Changes as stream is processed unlike other readers"""
        raise RuntimeError("IMDReader: n_frames is unknown")

    def next(self):
        """Don't rewind after iteration. When _reopen() is called,
        an error will be raised
        """
        try:
            ts = self._read_next_timestep()
        except (EOFError, IOError):
            # Don't rewind here like we normally would
            raise StopIteration from None
        else:
            for auxname, reader in self._auxs.items():
                ts = self._auxs[auxname].update_ts(ts)

            ts = self._apply_transformations(ts)

        return ts

    def rewind(self):
        """Raise error on rewind"""
        raise RuntimeError("IMDReader: Stream-based readers can't be rewound")

    @staticmethod
    def _format_hint(thing):
        try:
            parse_host_port(thing)
        except:
            return False
        return True

    def close(self):
        """Gracefully shut down the reader. Stops the producer thread."""
        logger.debug("IMDReader close() called")
        self._imdclient.stop()
        # NOTE: removeme after testing
        logger.debug("IMDReader shut down gracefully.")

    # Incompatible methods
    def copy(self):
        raise NotImplementedError("IMDReader does not support copying")

    def _reopen(self):
        if self._reopen_called:
            raise RuntimeError("IMDReader: Cannot reopen IMD stream")
        self._frame = -1
        self._reopen_called = True

    def __getitem__(self, frame):
        """This method from ProtoReader must be overridden
        to prevent slicing that doesn't make sense in a stream.
        """
        raise RuntimeError("IMDReader: Trajectory can only be read in for loop")

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
        if start is not None:
            raise ValueError(
                "IMDReader: Cannot slice a stream, 'start' must be None"
            )
        if stop is not None:
            raise ValueError(
                "IMDReader: Cannot slice a stream, 'stop' must be None"
            )
        if step is not None:
            if isinstance(step, numbers.Integral):
                if step != 1:
                    raise ValueError(
                        "IMDReader: Cannot slice a stream, 'step' must be None or 1"
                    )

        return start, stop, step

    def __getstate__(self):
        raise NotImplementedError("IMDReader does not support pickling")

    def __setstate__(self, state: object):
        raise NotImplementedError("IMDReader does not support pickling")
