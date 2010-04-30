"""Reading of Gromacs xtc trajectories."""

import core

class Timestep(core.Timestep):
    """Timestep for a Gromacs XTC trajectory."""

class XTCReader(core.TrjReader):
    """Read `Gromacs <www.gromacs.org>` XTC trajectory."""
    format = "XTC"
    _Timestep = Timestep

class XTCWriter(core.TrjWriter):
    """Write a `Gromacs <www.gromacs.org>` XTC trajectory."""
    format = "XTC"

