# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://mdanalysis.googlecode.com
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
   :mod:`MDAnalysis.core.units` for other recognized values.

.. SeeAlso:: For further discussion follow the reports for `Issue 84`_ and `Issue 64`_.

.. _LAMMPS: http://lammps.sandia.gov/
.. _DCD: http://lammps.sandia.gov/doc/dump.html
.. _CHARMM trajectory: http://www.charmm.org/documentation/c36b1/dynamc.html#%20Trajectory
.. _AKMA: http://www.charmm.org/documentation/c36b1/usage.html#%20AKMA
.. _`Issue 64`: http://code.google.com/p/mdanalysis/issues/detail?id=64
.. _`Issue 84`: http://code.google.com/p/mdanalysis/issues/detail?id=84

Classes
-------

.. autoclass:: Timestep
   :members:
   :inherited-members:
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

import DCD
from MDAnalysis.core import units
from MDAnalysis.topology.LAMMPSParser import DATAParser
import MDAnalysis.core
import base


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
    "Angstrom". See :mod:`MDAnalysis.core.units` for other recognized
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


class DATATimestep(base.Timestep):
    """Data file time step"""

    @property
    def dimensions(self):
        """unitcell dimensions (*A*, *B*, *C*, *alpha*, *beta*, *gamma*)

        lengths *a*, *b*, *c* are in the MDAnalysis length unit (Å), and
        angles are in degrees.

        :attr:`dimensions` is read-only because it transforms the
        actual format of the unitcell (which differs between different
        trajectory formats) to the representation described here,
        which is used everywhere in MDAnalysis.
        """
        return self._unitcell

    @dimensions.setter
    def dimensions(self, box):
        """Set unitcell information 

        .. versionadded:: 0.9.0
        """
        self._unitcell[:] = box


class DATAReader(base.Reader):
    """
    Reads a single frame of coordinate information from a LAMMPS DATA file.

    .. versionadded:: 0.9.0
    """
    format = 'DATA'
    units = {'time': None, 'length': 'Angstrom'}

    def __init__(self, pdbfilename, numatoms=None, convert_units=None, **kwargs):
        self.filename = pdbfilename

        if convert_units is None:
            convert_units = MDAnalysis.core.flags['convert_lengths']
        self.convert_units = convert_units  # convert length and time to base units

        if numatoms is None:  # this should be done by parsing DATA first
            raise ValueError("DATAReader requires numatoms keyword")
        self.numatoms = numatoms

        self.ts = DATATimestep(self.numatoms)
        with DATAParser(self.filename) as p:
            p.read_DATA_timestep(self.ts)

        self.numframes = 1
        self.fixed = 0  # parse B field for fixed atoms?
        self.skip = 1
        self.periodic = False
        self.delta = 0
        self.skip_timestep = 1

        self.ts.frame = 1
        if self.convert_units:
            self.convert_pos_from_native(self.ts._pos)  # in-place !

    def __iter__(self):
        yield self.ts  # just a single frame available
        raise StopIteration

    def _read_frame(self, frame):
        """Bounds checking of frame is done in __iter__ """
        return self.ts

    def _read_next_timestep(self):
        # DATA file only contains a single frame
        raise IOError
