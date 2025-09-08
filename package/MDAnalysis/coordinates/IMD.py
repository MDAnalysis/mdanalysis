"""
IMDReader --- :mod:`MDAnalysis.coordinates.IMD`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This module provides support for reading molecular dynamics simulation data via the 
`Interactive Molecular Dynamics (IMD) protocol v3 <https://imdclient.readthedocs.io/en/latest/protocol_v3.html>`_.
The IMD protocol allows two-way communicating molecular simulation data through a socket.
Via IMD, a simulation engine sends data to a receiver (in this case, the IMDClient) and the receiver can send forces and specific control
requests (such as pausing, resuming, or terminating the simulation) back to the simulation engine. 

.. note::
   This reader only supports IMDv3, which is implemented in GROMACS, LAMMPS, and NAMD at varying
   stages of development. See the `imdclient simulation engine docs`_ for more. 
   While IMDv2 is widely available in simulation engines, it was designed primarily for visualization
   and gaps are allowed in the stream (i.e., an inconsistent number of integrator time steps between transmitted coordinate arrays is allowed)

The :class:`IMDReader` connects to a simulation via a socket and receives coordinate,
velocity, force, and energy data as the simulation progresses. This allows for real-time
monitoring and analysis of ongoing simulations. It uses the `imdclient package <https://github.com/Becksteinlab/imdclient>`_ 
(dependency) to implement the IMDv3 protocol and manage the socket connection and data parsing.

.. seealso::
   :class:`IMDReader`
      Technical details and parameter options for the reader class
   
   `imdclient documentation <https://imdclient.readthedocs.io/>`_
      Complete documentation for the IMDClient package
      
   `IMDClient GitHub repository <https://github.com/Becksteinlab/imdclient>`_
      Source code and development resources

.. _`imdclient simulation engine docs`: https://imdclient.readthedocs.io/en/latest/usage.html

Usage Example
-------------

As an example of reading a stream, after configuring GROMACS to run a simulation with IMDv3 enabled
(see the `imdclient simulation engine docs`_ for 
up-to-date resources on configuring each simulation engine), use the following commands:

.. code-block:: bash

    gmx mdrun -v -nt 4 -imdwait -imdport 8889

The :class:`~MDAnalysis.coordinates.IMD.IMDReader` can then connect to the running simulation and stream data in real time:

.. code-block:: python

    import MDAnalysis as mda
    u = mda.Universe("topol.tpr", "imd://localhost:8889", buffer_size=10*1024*1024)

    print("    time [         position         ] [         velocity         ] [           force          ] [            box           ]")
    sel = u.select_atoms("all")  # Select all atoms; adjust selection as needed
    for ts in u.trajectory:
        print(f'{ts.time:8.3f} {sel[0].position} {sel[0].velocity} {sel[0].force} {u.dimensions[0:3]}')

.. important::
   **Jupyter Notebook Users**: When using IMDReader in Jupyter notebooks, be aware that 
   **kernel restarts will not gracefully close active IMD connections**. This can leave 
   socket connections open, potentially preventing new connections to the same stream.

   Always use ``try/except/finally`` blocks to ensure proper cleanup:

   .. code-block:: python

       import MDAnalysis as mda
       
       try:
           u = mda.Universe("topol.tpr", "imd://localhost:8889")
       except Exception as e:
           print(f"Error during connection: {e}")
       else:
           try:
               # Your analysis code here
               for ts in u.trajectory:
                   # Process each frame
                   pass
           finally:
               # Ensure connection is closed
               u.trajectory.close()

   Always explicitly call ``u.trajectory.close()`` when finished with analysis to 
   ensure connection is closed properly.

Important Limitations
---------------------

.. warning::  
   The IMDReader has some important limitations that are inherent in streaming data.  

Since IMD streams data in real-time from a running simulation, it has fundamental 
constraints that differ from traditional trajectory readers:

* **No random access**: Cannot jump to arbitrary frame numbers or seek backwards  
* **Forward-only access**: You can only move forward through frames as they arrive
* **No trajectory length**: The total number of frames is unknown until the simulation ends
* **Single-use iteration**: Cannot restart or rewind once the stream has been consumed
* **No independent copies**: Cannot create separate reader instances for the same stream
* **No stream restart**: Cannot reconnect or reopen once the connection is closed
* **No bulk operations**: Cannot extract all data at once using timeseries methods
* **Limited multiprocessing**: Cannot split reader across processes for parallel analysis
* **Single client connection**: Only one reader can connect to an IMD stream at a time
* **No trajectory Writing**: Complimentary IMD Writer class is not available for streaming data

.. seealso::
   See :class:`~MDAnalysis.coordinates.base.StreamReaderBase` for technical details.

Multiple Client Connections
---------------------------

The ability to establish multiple simultaneous connections to the same IMD port is 
**MD engine implementation dependent**. Some simulation engines may allow multiple 
clients to connect concurrently, while others may reject or fail additional connection 
attempts.

See the `imdclient simulation engine docs`_ for further details. 

.. important::
   Even when multiple connections are supported by the simulation engine, each connection 
   receives its own independent data stream. These streams may contain different data 
   depending on the simulation engine's configuration, so multiple connections should 
   not be assumed to provide identical data streams.

Classes
-------

.. autoclass:: IMDReader
   :members:
   :inherited-members:

"""

import numpy as np
import logging
import warnings

from MDAnalysis.coordinates import core
from MDAnalysis.lib.util import store_init_arguments
from MDAnalysis.coordinates.base import StreamReaderBase


from packaging.version import Version

MIN_IMDCLIENT_VERSION = Version("0.2.2")

try:
    import imdclient
    from imdclient.IMDClient import IMDClient
    from imdclient.utils import parse_host_port
except ImportError:
    HAS_IMDCLIENT = False
    imdclient_version = Version("0.0.0")

    # Allow building documentation without imdclient
    import types

    class MockIMDClient:
        pass

    imdclient = types.ModuleType("imdclient")
    imdclient.IMDClient = MockIMDClient
    imdclient.__version__ = "0.0.0"

else:
    HAS_IMDCLIENT = True
    imdclient_version = Version(imdclient.__version__)

    # Check for compatibility: currently needs to be >=0.2.2
    if imdclient_version < MIN_IMDCLIENT_VERSION:
        warnings.warn(
            f"imdclient version {imdclient_version} is too old; "
            f"need at least {MIN_IMDCLIENT_VERSION}, Your installed version of "
            "imdclient will NOT be used.",
            category=RuntimeWarning,
        )
        HAS_IMDCLIENT = False

logger = logging.getLogger("MDAnalysis.coordinates.IMDReader")


class IMDReader(StreamReaderBase):
    """
    Coordinate reader implementing the IMDv3 protocol for streaming simulation data.

    This class handles the technical aspects of connecting to IMD-enabled simulation
    engines and processing the incoming data stream. For usage examples and protocol
    overview, see the module documentation above.

    The reader manages socket connections, data buffering, and frame parsing according
    to the IMDv3 specification. It automatically handles different data packet types
    (coordinates, velocities, forces, energies, timing) and populates MDAnalysis
    timestep objects accordingly.

    Parameters
    ----------
    filename : a string of the form "imd://host:port" where host is the hostname
        or IP address of the listening simulation engine's IMD server and port
        is the port number.
    n_atoms : int (optional)
        number of atoms in the system. defaults to number of atoms
        in the topology. Don't set this unless you know what you're doing.
    buffer_size: int (optional) default=10*(1024**2)
        number of bytes of memory to allocate to the :class:`~imdclient.IMDClient.IMDClient`'s
        internal buffer. Defaults to 10 megabytes. Larger buffers can improve
        performance for analyses with periodic heavy computation.
    **kwargs : dict (optional)
        keyword arguments passed to the constructed :class:`~imdclient.IMDClient.IMDClient`

    Notes
    -----
    The IMDReader provides access to additional simulation data through the timestep's
    `data` attribute (`ts.data`). The following keys may be available depending on
    what the simulation engine transmits:

    * `dt` : float
        Time step size in picoseconds (from the `IMD_TIME`_ packet of the IMDv3 protocol)
    * `step` : int
        Current simulation step number (from the `IMD_TIME`_ packet of the IMDv3 protocol)
    * Energy terms : float
        Various energy components (e.g., 'potential', 'kinetic', 'total', etc.)
        from the `IMD_ENERGIES`_ packet of the IMDv3 protocol.

    .. _IMD_TIME: https://imdclient.readthedocs.io/en/latest/protocol_v3.html#time
    .. _IMD_ENERGIES: https://imdclient.readthedocs.io/en/latest/protocol_v3.html#energies

    .. note::
       For important limitations inherent to streaming data, see the module documentation above
       and :class:`~MDAnalysis.coordinates.base.StreamReaderBase` for more technical details.

    .. versionadded:: 2.10.0
    """

    format = "IMD"

    @store_init_arguments
    def __init__(
        self,
        filename,
        n_atoms=None,
        buffer_size=10 * (1024**2),
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

        try:
            host, port = parse_host_port(filename)
        except ValueError as e:
            raise ValueError(f"IMDReader: Invalid IMD URL '{filename}': {e}")

        # This starts the simulation
        self._imdclient = IMDClient(
            host, port, n_atoms, buffer_size=buffer_size, **kwargs
        )

        imdsinfo = self._imdclient.get_imdsessioninfo()
        if imdsinfo.version != 3:
            raise ValueError(
                f"IMDReader: Detected IMD version v{imdsinfo.version}, "
                + "but IMDReader is only compatible with v3"
            )

        self.ts = self._Timestep(
            self.n_atoms,
            positions=imdsinfo.positions,
            velocities=imdsinfo.velocities,
            forces=imdsinfo.forces,
            **self._ts_kwargs,
        )

        try:
            self._read_next_timestep()
        except EOFError as e:
            raise RuntimeError(f"IMDReader: Read error: {e}") from e

    def _read_frame(self, frame):

        imdf = self._imdclient.get_imdframe()

        self._frame = frame
        self._load_imdframe_into_ts(imdf)

        logger.debug("IMDReader: Loaded frame %d", self._frame)
        return self.ts

    def _load_imdframe_into_ts(self, imdf):
        self.ts.frame = self._frame
        if imdf.time is not None:
            self.ts.time = imdf.time
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
        if not isinstance(thing, str):
            return False
        # a weaker check for type hint
        if thing.startswith("imd://"):
            return True
        else:
            return False

    def close(self):
        """Gracefully shut down the reader. Stops the producer thread."""
        logger.debug("IMDReader close() called")
        if self._imdclient is not None:
            self._imdclient.stop()
        logger.debug("IMDReader shut down gracefully.")
