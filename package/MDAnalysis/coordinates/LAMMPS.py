# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
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
the unit from the file. By default we are assuming that the unit for
length is the ångström and for the time is the picosecond. If this is
not true then the user *should supply the appropriate units* in the as
keywords *timeunit* and/or *lengthunit* to :class:`DCDWriter` and
:class:`DCDReader`.

.. Note::

   Lennard-Jones units are not implemented. See
   :mod:`MDAnalysis.units` for other recognized values.

.. SeeAlso:: For further discussion follow the reports for `Issue 84`_ and `Issue 64`_.

.. _LAMMPS: http://lammps.sandia.gov/
.. _DCD: http://lammps.sandia.gov/doc/dump.html
.. _CHARMM trajectory: http://www.charmm.org/documentation/c36b1/dynamc.html#%20Trajectory
.. _AKMA: http://www.charmm.org/documentation/c36b1/usage.html#%20AKMA
.. _`Issue 64`: https://github.com/MDAnalysis/mdanalysis/issues/64
.. _`Issue 84`: https://github.com/MDAnalysis/mdanalysis/issues/84

Classes
-------

.. autoclass:: DCDReader
   :members:
   :inherited-members:
.. autoclass:: DCDWriter
   :members:
   :inherited-members:
.. autoclass:: DATAReader
   :members:
   :inherited-members:
"""

from . import DCD
from .. import units
from ..topology.LAMMPSParser import DATAParser
from . import base


class DCDWriter(DCD.DCDWriter):
    """Write a LAMMPS_ DCD trajectory.

    The units can be set from the constructor with the keyword
    arguments *timeunit* and *lengthunit*. The defaults are "ps" and
    "Angstrom". See :mod:`MDAnalysis.units` for other recognized
    values.
    """
    format = "DCD"

    def __init__(self, *args, **kwargs):
        self.units = {'time': 'ps', 'length': 'Angstrom'}  # must be instance level
        self.units['time'] = kwargs.pop('timeunit', self.units['time'])
        self.units['length'] = kwargs.pop('lengthunit', self.units['length'])
        for unit_type, unit in self.units.items():
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
    "Angstrom". See :mod:`MDAnalysis.units` for other recognized
    values.
    """
    format = "DCD"

    def __init__(self, dcdfilename, **kwargs):
        self.units = {'time': 'ps', 'length': 'Angstrom'}  # must be instance level
        self.units['time'] = kwargs.pop('timeunit', self.units['time'])
        self.units['length'] = kwargs.pop('lengthunit', self.units['length'])
        for unit_type, unit in self.units.items():
            try:
                if units.unit_types[unit] != unit_type:
                    raise TypeError("LAMMPS DCDReader: wrong unit %r for unit type %r" % (unit, unit_type))
            except KeyError:
                raise ValueError("LAMMPS DCDReader: unknown unit %r" % unit)
        super(DCDReader, self).__init__(dcdfilename, **kwargs)


class DATAReader(base.SingleFrameReader):
    """Reads a single frame of coordinate information from a LAMMPS DATA file.

    .. versionadded:: 0.9.0
    .. versionchanged:: 0.11.0
       Frames now 0-based instead of 1-based
    """
    format = 'DATA'
    units = {'time': None, 'length': 'Angstrom'}

    def __init__(self, filename, **kwargs):
        self.n_atoms = kwargs.pop('n_atoms', None)
        if self.n_atoms is None:  # this should be done by parsing DATA first
            raise ValueError("DATAReader requires n_atoms keyword")
        super(DATAReader, self).__init__(filename, **kwargs)

    def _read_first_frame(self):
        with DATAParser(self.filename) as p:
            self.ts = p.read_DATA_timestep(self.n_atoms, self._Timestep,
                                           self._ts_kwargs)

        self.ts.frame = 0
        if self.convert_units:
            self.convert_pos_from_native(self.ts._pos)  # in-place !
