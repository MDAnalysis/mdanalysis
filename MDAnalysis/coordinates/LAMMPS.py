# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2011 Naveen Michaud-Agrawal,
#               Elizabeth J. Denning, Oliver Beckstein,
#               and contributors (see website for details)
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
#     N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
#     O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
#     Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
#     doi:10.1002/jcc.21787
#

"""LAMMPS DCD trajectory I/O  --- :mod:`MDAnalysis.coordinates.LAMMPS`
======================================================================

Classes to read and write LAMMPS_ DCD binary
trajectories. Trajectories can be read regardless of system-endianness
as this is auto-detected.

LAMMPS can write DCD_ trajectories but unlike a `CHARMM trajectory`_
(which is often called a DCD even though CHARMM itself calls them
"trj") the time unit is not fixed to be the AKMA_ time unit (20 AKMA
is 0.978 picoseconds or 1 AKMA = 4.888821e-14 s) but can depend on
settings in LAMMPS. The most common case appears to be that the time
step is actually recorded in picoseconds. Other cases are unit-less
Lennard-Jones time units.

This presents a problem for MDAnalysis because it cannot autodetect
the unit from the file.

.. _LAMMPS: http://lammps.sandia.gov/
.. _DCD: http://lammps.sandia.gov/doc/dump.html
.. _CHARMM trajectory: http://www.charmm.org/documentation/c36b1/dynamc.html#%20Trajectory
.. _AKMA: http://www.charmm.org/documentation/c36b1/usage.html#%20AKMA

Classes
-------

.. autoclass:: Timestep
   :members:
   :inherited-members:
.. autoclass:: TRRReader
   :members:
   :inherited-members:
.. autoclass:: TRRWriter
   :members:
   :inherited-members:

"""

import DCD
from MDAnalysis.core import units

class Timestep(DCD.Timestep):
    """LAMMPS trajectory time step"""

class DCDWriter(DCD.DCDWriter):
    """Write a LAMMPS_ DCD trajectory.

    The units can be set from the constructor with the keyword
    arguments *timeunit* and *lengthunit*. The defaults are "ps" and
    "Angstrom". See :mod:`MDAnalysis.core.units` for other recognized
    values.
    """
    format = "DCD"
    def __init__(self, *args, **kwargs):
        self.units = {'time': 'ps', 'length':'Angstrom'}      # must be instance level
        self.units['time'] = kwargs.pop('timeunit', self.units['time'])
        self.units['length'] = kwargs.pop('lengthunit', self.units['length'])
        for unit_type,unit in self.units.items():
            try:
                if units.unit_types[unit] != unit_type:
                    raise TypeError("LAMMPS DCDWriter: wrong unit %r for unit type %r" % (unit, unit_type))
            except KeyError:
                raise ValueError("LAMMPS DCDWriter: unknown unit %r" % unit)
        super(DCDWriter, self).__init__(*args, **kwargs)

class DCDReader(DCD.DCDReader):
    """Read a LAMMPS_ DCD trajectory.

    The units can be set from the constructor with the keyword
    arguments *timeunit* and *lengthunit*. The defaults are "ps" and
    "Angstrom". See :mod:`MDAnalysis.core.units` for other recognized
    values.
    """
    format = "DCD"

    def __init__(self, dcdfilename, **kwargs):
        self.units = {'time': 'ps', 'length':'Angstrom'}      # must be instance level
        self.units['time'] = kwargs.pop('timeunit', self.units['time'])
        self.units['length'] = kwargs.pop('lengthunit', self.units['length'])
        for unit_type,unit in self.units.items():
            try:
                if units.unit_types[unit] != unit_type:
                    raise TypeError("LAMMPS DCDReader: wrong unit %r for unit type %r" % (unit, unit_type))
            except KeyError:
                raise ValueError("LAMMPS DCDReader: unknown unit %r" % unit)
        super(DCDReader, self).__init__(dcdfilename, **kwargs)
