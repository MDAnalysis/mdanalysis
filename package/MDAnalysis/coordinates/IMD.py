"""
IMDReader --- :mod:`MDAnalysis.coordinates.IMD`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:class:`MDAnalysis.coordinates.IMD.IMDReader` is a class that implements the
`Interactive Molecular Dynamics (IMD) protocol <https://imdclient.readthedocs.io/en/latest/protocol_v3.html>`_ for reading simulation
data using the IMDClient (see `imdclient <https://github.com/Becksteinlab/imdclient>`_).
The protocol allows two-way communicating molecular simulation data through a socket.
Via IMD, a simulation engine sends data to a receiver (in this case, the IMDClient) and the receiver can send forces and specific control
requests (such as pausing, resuming, or terminating the simulation) back to the simulation engine. 

IMDv3, the newest version of the protocol, is the one supported by this reader class and is implemented in GROMACS, LAMMPS, and NAMD at varying
stages of development. See the `imdclient simulation engine docs <https://imdclient.readthedocs.io/en/latest/usage.html>`_ for more.

IMDv2, the first version to be broadly adopted, is currently available as a part of official releases of GROMACS, LAMMPS, and NAMD. However,
this reader class does not currently provide support for it since it was designed for visualization and gaps are allowed in the stream
(i.e., an inconsistent number of integrator time steps between transmitted coordinate arrays is allowed)

As an example of reading a stream, after configuring GROMACS to run a simulation with IMDv3 enabled
(see the `imdclient simulation engine docs <https://imdclient.readthedocs.io/en/latest/usage.html>`_ for 
up-to-date resources on configuring each simulation engine), use the following commands:

.. code-block:: bash

    gmx grompp -f run-NPT_imd-v3.mdp -c conf.gro -p topol.top -o topol.tpr
    gmx mdrun -v -nt 4 -imdwait -imdport 8889

The :class:`MDAnalysis.coordinates.IMD.IMDReader` can then connect to the running simulation and stream data in real time:

.. code-block:: python

    import MDAnalysis as mda
    u = mda.Universe("topol.tpr", "imd://localhost:8889", buffer_size=10*1024*1024)

    print("    time [         position         ] [         velocity         ] [           force          ] [            box           ]")
    sel = u.select_atoms("all")  # Select all atoms; adjust selection as needed
    for ts in u.trajectory:
        print(f'{ts.time:8.3f} {sel[0].position} {sel[0].velocity} {sel[0].force} {u.dimensions[0:3]}')

Details about the IMD protocol and usage examples can be found in the
`imdclient <https://github.com/Becksteinlab/imdclient>`_ repository.




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
            f"need at least {imdclient_version}, Your installed version of "
            "imdclient will NOT be used.",
            category=RuntimeWarning,
        )
        HAS_IMDCLIENT = False

logger = logging.getLogger("MDAnalysis.coordinates.IMDReader")


class IMDReader(StreamReaderBase):
    """
    Reader that supports the Interactive Molecular Dynamics (IMD) protocol v3 for reading
    simulation data using the IMDClient.

    By using the keyword `buffer_size`, you can change the amount of memory the 
    :class:`~imdclient.IMDClient.IMDClient` allocates to its internal buffer.
    The buffer size determines how many frames can be stored in memory as data is received
    from the socket and awaits reading by the client. For analyses that periodically perform
    heavier computation at fixed intervals, say for example once every 200 received frames,
    increasing this value will decrease the amount of time the simulation engine spends in a
    paused state and potentially decrease total analysis time, but will require more RAM.

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
        internal buffer. Defaults to 10 megabytes.
    kwargs : dict (optional)
        keyword arguments passed to the constructed :class:`~imdclient.IMDClient.IMDClient`

    Notes
    -----
    The IMDReader provides access to additional simulation data through the timestep's
    `data` attribute (`ts.data`). The following keys may be available depending on
    what the simulation engine transmits:

    * `dt` : float
        Time step size in picoseconds (from the `IMD_TIME`_ packet of the IMDv3 protocol)
    * `step` : int
        Current simulation step number (from the `IMD_TIME`_ packet of the IMDv3 prtotocol)
    * Energy terms : float
        Various energy components (e.g., 'potential', 'kinetic', 'total', etc.)
        from the `IMD_ENERGIES`_ packet of the IMDv3 protocol.

    .. _IMD_TIME: https://imdclient.readthedocs.io/en/latest/protocol_v3.html#time
    .. _IMD_ENERGIES: https://imdclient.readthedocs.io/en/latest/protocol_v3.html#energies
        
    
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
