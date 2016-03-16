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

LAMMPS can `write DCD`_ trajectories but unlike a `CHARMM trajectory`_
(which is often called a DCD even though CHARMM itself calls them
"trj") the time unit is not fixed to be the AKMA_ time unit (20 AKMA
is 0.978 picoseconds or 1 AKMA = 4.888821e-14 s) but can depend on
settings in LAMMPS. The most common case for biomolecular simulations
appears to be that the time step is recorded in femtoseconds (command
`units real`_ in the input file) and lengths in ångströms. Other cases
are unit-less Lennard-Jones time units.

This presents a problem for MDAnalysis because it cannot autodetect
the unit from the file. By default we are assuming that the unit for
length is the ångström and for the time is the femtosecond. If this is
not true then the user *should supply the appropriate units* in the
keywords *timeunit* and/or *lengthunit* to :class:`DCDWriter` and
:class:`~MDAnalysis.core.universe.Universe` (which then calls
:class:`DCDReader`).

.. Rubric:: Example: Loading a LAMMPS simulation

To load a LAMMPS simulation from a LAMMPS data file (using the
:class:`~MDAnalysis.topology.LAMMPSParser.DATAParser`) together with a
LAMMPS DCD with "*real*" provide the keyword *format="LAMMPS*"::

  u = MDAnalysis.Universe("lammps.data", "lammps_real.dcd", format="LAMMPS")

If the trajectory uses *units nano* then use ::

  u = MDAnalysis.Universe("lammps.data", "lammps_nano.dcd", format="LAMMPS",
                          lengthunit="nm", timeunit="ns")

.. Note::

   Lennard-Jones units are not implemented. See :mod:`MDAnalysis.units` for
   other recognized values and the documentation for the LAMMPS `units
   command`_.

.. SeeAlso::

   For further discussion follow the reports for `Issue 84`_ and `Issue 64`_.

.. _LAMMPS: http://lammps.sandia.gov/
.. _write DCD: http://lammps.sandia.gov/doc/dump.html
.. _CHARMM trajectory: http://www.charmm.org/documentation/c36b1/dynamc.html#%20Trajectory
.. _AKMA: http://www.charmm.org/documentation/c36b1/usage.html#%20AKMA
.. _units real: http://lammps.sandia.gov/doc/units.html
.. _units command: http://lammps.sandia.gov/doc/units.html
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
.. autoclass:: DATAWriter
   :members:
   :inherited-members:
"""

from six.moves import zip, range
import numpy as np

from ..core import flags
from ..lib import util, mdamath
from . import DCD
from .. import units
from ..topology.LAMMPSParser import DATAParser
from ..exceptions import NoDataError
from . import base

btype_sections = {'bond':'Bonds', 'angle':'Angles',
                  'dihedral':'Dihedrals', 'improper':'Impropers'}

class DCDWriter(DCD.DCDWriter):
    """Write a LAMMPS_ DCD trajectory.

    The units can be set from the constructor with the keyword
    arguments *timeunit* and *lengthunit*. The defaults are "fs" and
    "Angstrom". See :mod:`MDAnalysis.units` for other recognized
    values.
    """
    format = 'DCD'
    flavor = 'LAMMPS'

    def __init__(self, *args, **kwargs):
        self.units = {'time': 'fs', 'length': 'Angstrom'}  # must be instance level
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
    arguments *timeunit* and *lengthunit*. The defaults are "fs" and
    "Angstrom", corresponding to LAMMPS `units style`_ "**real**". See
    :mod:`MDAnalysis.units` for other recognized values.

    .. _units style: http://lammps.sandia.gov/doc/units.html
    """
    format = 'DCD'
    flavor = 'LAMMPS'

    def __init__(self, dcdfilename, **kwargs):
        self.units = {'time': 'fs', 'length': 'Angstrom'}  # must be instance level
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
    units = {'time': None, 'length': 'Angstrom', 'velocity': 'Angstrom/fs'}

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
            try:
                self.convert_velocities_from_native(self.ts._velocities)  # in-place !
            except AttributeError:
                pass

class DATAWriter(base.Writer):
    """Write out the current time step as a LAMMPS DATA file.

    This writer supports the sections Atoms, Masses, Velocities, Bonds,
    Angles, Dihedrals, and Impropers. This writer will write the header
    and these sections (if applicable). Atoms section is written in the
    "full" sub-style if charges are available or "molecular" sub-style
    if they are not. Molecule id is set to 0 for all atoms.

    .. Note::

       This writer assumes "conventional" or "real" LAMMPS units where length
       is measured in Angstroms and velocity is measured in Angstroms per
       femtosecond. To write in different units, specify `lengthunit`
    """
    format = 'DATA'

    def __init__(self, filename, convert_units=None, **kwargs):
        """Set up a DATAWriter

        :Arguments:
           *filename*
              output filename
        """
        self.filename = util.filename(filename, ext='data')

        if convert_units is None:
            convert_units = flags['convert_lengths']
        self.convert_units = convert_units

        self.units = {'time': 'fs', 'length': 'Angstrom'}
        self.units['length'] = kwargs.pop('lengthunit', self.units['length'])
        self.units['time'] = kwargs.pop('timeunit', self.units['time'])
        self.units['velocity'] = kwargs.pop('velocityunit', \
                self.units['length']+'/'+self.units['time'])

    def _write_atoms(self, atoms):
        self.file.write('\n')
        self.file.write('Atoms\n')
        self.file.write('\n')

        try:
            charges = atoms.charges
        except (NoDataError, AttributeError):
            has_charges = False
        else:
            has_charges = True

        indices = atoms.indices + 1
        types = atoms.types.astype(np.int32)

        if self.convert_units:
            coordinates = self.convert_pos_to_native(atoms.positions, inplace=False)

        if has_charges:
            for index, atype, charge, coords in zip(indices, types, charges,
                    coordinates):
                self.file.write('%d 0 %d %f %f %f %f\n'%(index, atype,
                        charge, coords[0], coords[1], coords[2]))
        else:
            for index, atype, coords in zip(indices, types, coordinates):
                self.file.write('%d 0 %d %f %f %f\n'%(index, atype,\
                        coords[0], coords[1], coords[2]))

    def _write_velocities(self, atoms):
        self.file.write('\n')
        self.file.write('Velocities\n')
        self.file.write('\n')
        indices = atoms.indices + 1
        velocities = self.convert_velocities_to_native(atoms.velocities, inplace=False)
        for index, vel in zip(indices, velocities):
            self.file.write('%d %f %f %f\n'%(index, vel[0], vel[1], vel[2]))

    def _write_masses(self, atoms):
        self.file.write('\n')
        self.file.write('Masses\n')
        self.file.write('\n')
        mass_dict = {}
        for atype in set(atoms.types.astype(np.int32)):
            masses = set(atoms.select_atoms('type %s'%(atype)).masses)
            mass_dict[atype] = masses.pop()
            if masses:
                raise ValueError('LAMMPS DATAWriter: to write data file, '+\
                        'atoms with same type must have same mass')
        for atype, mass in mass_dict.items():
            self.file.write('%d %f\n'%(atype, mass))

    def _write_bonds(self, bonds):
        num = len(bonds[0])
        self.file.write('\n')
        self.file.write('%s\n'%btype_sections[bonds.btype])
        self.file.write('\n')
        for bond, i in zip(bonds, range(1, len(bonds)+1)):
            self.file.write('%d %s '%(i, bond.type)+\
                    ' '.join((bond.atoms.indices + 1).astype(str))+'\n')

    def _write_dimensions(self, dimensions):
        """
            Convert dimensions to triclinic vectors, convert lengths to native
            units and then write the dimensions section
        """
        if self.convert_units:
            triv = self.convert_pos_to_native(mdamath.triclinic_vectors(\
                    dimensions),inplace=False)
        self.file.write('\n')
        self.file.write('%f %f xlo xhi\n'%(0., triv[0][0]))
        self.file.write('%f %f ylo yhi\n'%(0., triv[1][1]))
        self.file.write('%f %f zlo zhi\n'%(0., triv[2][2]))
        if any(np.array([triv[1][0], triv[2][0], triv[2][1]]) != 0):
            self.file.write('%f %f %f xy xz yz\n'%(triv[1][0], triv[2][0],\
                    triv[2][1]))
        self.file.write('\n')

    def write(self, selection, frame=None):
        """Write selection at current trajectory frame to file, including sections
        Atoms, Masses, Velocities, Bonds, Angles, Dihedrals, and Impropers (if these
        are defined). Atoms section is written in the "full" sub-style if charges
        are available or "molecular" sub-style if they are not.
        Molecule id in atoms section is set to to 0.

        No other sections are written to the DATA file.
        As of this writing, other sections are not parsed into the topology
        by the :class:`DATAReader`.

        .. Note::

           If the selection includes a partial fragment, then the outputted DATA
           file will be invalid, because it will describe bonds between atoms
           which do not exist in the Atoms section.

        :Arguments:
           *selection*
                MDAnalysis AtomGroup (selection or Universe.atoms) or also Universe
        :Keywords:
           *frame*
               optionally move to frame number *frame*
           *renumber_atoms*
                Set to ``True`` to make the writer renumber atoms consecutively from
                1 to N before writing. This may be necessary for some LAMMPS features
                such as `velocity all create` which require consecutively numbered
                atoms. (not yet implemented)
        """
        u = selection.universe
        if frame is not None:
            u.trajectory[frame]
        else:
            frame = u.trajectory.ts.frame

        # make sure to use atoms (Issue 46)
        atoms = selection.atoms

        # in future assert that atoms.fragments.atoms == atoms
        # but right now fragments returns a tuple instead of a Group

        # check that types can be converted to ints if they aren't
        # ints already.
        try:
            u.atoms.types.astype(np.int32)
        except ValueError:
            raise ValueError('LAMMPS.DATAWriter: atom types must be '+\
                    'convertible to integers')

        try:
            velocities = atoms.velocities
        except (NoDataError, AttributeError):
            has_velocities = False
        else:
            has_velocities = True

        features = {}
        with util.openany(self.filename, 'w') as self.file:
            self.file.write('LAMMPS data file via MDAnalysis\n')
            self.file.write('\n')
            self.file.write('%12d  atoms\n'%len(atoms))

            attrs = [('bond', 'bonds'), ('angle', 'angles'),\
                ('dihedral', 'dihedrals'), ('improper', 'impropers')]

            for btype, attr_name in attrs:
                try:
                    features[btype] = atoms.__getattribute__(attr_name)
                    self.file.write('%12d  %s\n'%(len(features[btype]), attr_name))
                except AttributeError:
                    features[btype] = None
                    self.file.write('%12d  %s\n'%(0, attr_name))

            self.file.write('\n')
            self.file.write('%12d  atom types\n'%len(set(atoms.types)))

            for btype, attr in features.items():
                if attr is None:
                    self.file.write('%12d  %s types\n'%(0, btype))
                else:
                    self.file.write('%12d  %s types\n'%(len(attr.types()), btype))

            self._write_dimensions(atoms.dimensions)

            #if renumber_atoms is True:
            self._write_masses(atoms)
            self._write_atoms(atoms)
            for attr in features.values():
                if attr is None or len(attr) == 0:
                    continue
                self._write_bonds(attr)

            if has_velocities:
                self._write_velocities(atoms)
