"""Reading of `Gromacs XTC trajectories`_.

.. _Gromacs XTC trajectories: http://www.gromacs.org/Documentation/File_Formats/.xtc_File
.. _Gromacs: http://www.gromacs.org


.. SeeAlso:: :mod:`MDAnalysis.coordinates.xdrfile.libxdrfile` for low-level
   bindings to the Gromacs trajectory file formats

Classes
-------

.. autoclass:: Timestep
   :members:
   :inherited-members:
.. autoclass:: XTCReader
   :members:
   :inherited-members:
.. autoclass:: XTCWriter
   :members:
   :inherited-members:
"""

import core

class Timestep(core.Timestep):
    """Timestep for a Gromacs XTC trajectory."""

class XTCReader(core.TrjReader):
    """Read Gromacs_ XTC trajectory."""
    format = "XTC"
    _Timestep = Timestep

class XTCWriter(core.TrjWriter):
    """Write a Gromacs_ XTC trajectory."""
    format = "XTC"

