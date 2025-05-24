"""
IMDReader --- :mod:`MDAnalysis.coordinates.IMD`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Read and analyze simulation data interactively using `IMDClient`_.

.. _IMDClient: https://github.com/Becksteinlab/imdclient

Units
-----
The units in IMDv3 are fixed.

.. list-table::
   :widths: 10 10
   :header-rows: 1

   * - Measurement
     - Unit
   * - Length
     - angstrom
   * - Velocity
     - angstrom/picosecond
   * - Force
     - kilojoules/(mol*angstrom)
   * - Time
     - picosecond
   * - Energy
     - kilojoules/mol

Classes
-------

.. autoclass:: IMDReader
   :members:
   :inherited-members:

"""

import numpy as np
import logging

from MDAnalysis.coordinates import core
from MDAnalysis.lib.util import store_init_arguments
from MDAnalysis.coordinates.base import StreamReaderBase

try:
    import imdclient
    from imdclient.IMDClient import IMDClient
except ImportError:
    HAS_IMDCLIENT = False

    # Allow building documentation without imdclient
    import types

    class MockIMDClient:
        pass
    imdclient = types.ModuleType("imdclient")
    imdclient.IMDClient = MockIMDClient

else:
    HAS_IMDCLIENT = True

logger = logging.getLogger("MDAnalysis.coordinates.IMDReader")


class IMDReader(StreamReaderBase):
    """
    Reader for IMD protocol packets.

    Parameters
    ----------
    filename : a string of the form "imd://host:port" where host is the hostname
        or IP address of the listening GROMACS server and port
        is the port number.
    n_atoms : int (optional)
        number of atoms in the system. defaults to number of atoms
        in the topology. Don't set this unless you know what you're doing.
    kwargs : dict (optional)
        keyword arguments passed to the constructed :class:`IMDClient`
    """

    format = "IMD"
    one_pass = True

    @store_init_arguments
    def __init__(
        self,
        filename,
        convert_units=True,
        n_atoms=None,
        **kwargs,
    ):
        if not HAS_IMDCLIENT:
            raise ImportError(
                "IMDReader requires the imdclient package. "
                "Please install it with 'pip install imdclient'."
            )

        super(IMDReader, self).__init__(filename, **kwargs)

        self._imdclient = None
        logger.debug("IMDReader initializing")

        if n_atoms is None:
            raise ValueError("IMDReader: n_atoms must be specified")
        self.n_atoms = n_atoms

        host, port = parse_host_port(filename)

        # This starts the simulation
        self._imdclient = IMDClient(host, port, n_atoms, **kwargs)

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
        except StopIteration as e:
            raise RuntimeError("IMDReader: No data found in stream") from e

    def _read_frame(self, frame):

        try:
            imdf = self._imdclient.get_imdframe()
        except EOFError as e:
            raise e

        self._frame = frame
        self._load_imdframe_into_ts(imdf)

        logger.debug("IMDReader: Loaded frame %d", self._frame)
        return self.ts

    def _load_imdframe_into_ts(self, imdf):
        self.ts.frame = self._frame
        if imdf.time is not None:
            self.ts.time = imdf.time
            # NOTE: timestep.pyx "dt" method is suspicious bc it uses "new" keyword for a float
            self.ts.data["dt"] = imdf.dt
            self.ts.data["step"] = imdf.step
        if imdf.energies is not None:
            self.ts.data.update(
                {k: v for k, v in imdf.energies.items() if k != "step"}
            )
        if imdf.box is not None:
            self.ts.dimensions = core.triclinic_box(*imdf.box)
        if imdf.positions is not None:
            # must call copy because reference is expected to reset
            # see 'test_frame_collect_all_same' in MDAnalysisTests.coordinates.base
            np.copyto(self.ts.positions, imdf.positions)
        if imdf.velocities is not None:
            np.copyto(self.ts.velocities, imdf.velocities)
        if imdf.forces is not None:
            np.copyto(self.ts.forces, imdf.forces)

    @staticmethod
    def _format_hint(thing):
        try:
            parse_host_port(thing)
        except:
            return False
        return HAS_IMDCLIENT and True

    def close(self):
        """Gracefully shut down the reader. Stops the producer thread."""
        logger.debug("IMDReader close() called")
        if self._imdclient is not None:
            self._imdclient.stop()
        # NOTE: removeme after testing
        logger.debug("IMDReader shut down gracefully.")

# NOTE: think of other edge cases as well- should be robust
def parse_host_port(filename):
    if not filename.startswith("imd://"):
        raise ValueError("IMDReader: URL must be in the format 'imd://host:port'")
    # Check if the format is correct
    parts = filename.split("imd://")[1].split(":")
    if len(parts) == 2:
        host = parts[0]
        try:
            port = int(parts[1])
            return (host, port)
        except ValueError as e:
            raise ValueError("IMDReader: Port must be an integer") from e
    else:
        raise ValueError("IMDReader: URL must be in the format 'imd://host:port'")
