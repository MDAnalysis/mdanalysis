# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#


"""LAMMPS DCD trajectory and DATA I/O  --- :mod:`MDAnalysis.coordinates.LAMMPS`
===============================================================================

Classes to read and write LAMMPS_ DCD binary trajectories and LAMMPS DATA
files. Trajectories can be read regardless of system-endianness
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

Data file formats
-----------------

By default either the `atomic` or `full` atom styles are expected,
however this can be customised, see :ref:`atom_style_kwarg`.


Example: Loading a LAMMPS simulation
------------------------------------

To load a LAMMPS simulation from a LAMMPS data file (using the
:class:`~MDAnalysis.topology.LAMMPSParser.DATAParser`) together with a
LAMMPS DCD with "*real*" provide the keyword *format="LAMMPS*"::

    >>> u = MDAnalysis.Universe("lammps.data", "lammps_real.dcd", format="LAMMPS")

If the trajectory uses *units nano* then use ::

    >>> u = MDAnalysis.Universe("lammps.data", "lammps_nano.dcd", format="LAMMPS",
    ...                          lengthunit="nm", timeunit="ns")

To scan through a trajectory to find a desirable frame and write to a LAMMPS
data file,

>>> for ts in u.trajectory:
...     # analyze frame
...     if take_this_frame == True:
...         with mda.Writer('frame.data') as W:
...             W.write(u.atoms)
...         break

Note
----
Lennard-Jones units are not implemented. See :mod:`MDAnalysis.units`
for other recognized values and the documentation for the LAMMPS
`units command`_.

See Also
--------

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
from __future__ import absolute_import

from six.moves import zip, range
import numpy as np

from ..core import flags
from ..core.groups import requires
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
    format = 'LAMMPS'
    multiframe = True
    flavor = 'LAMMPS'

    def __init__(self, *args, **kwargs):
        self.units = {'time': 'fs', 'length': 'Angstrom'}  # must be instance level
        self.units['time'] = kwargs.pop('timeunit', self.units['time'])
        self.units['length'] = kwargs.pop('lengthunit', self.units['length'])
        for unit_type, unit in self.units.items():
            try:
                if units.unit_types[unit] != unit_type:
                    raise TypeError("LAMMPS DCDWriter: wrong unit {0!r} for unit type {1!r}".format(unit, unit_type))
            except KeyError:
                raise ValueError("LAMMPS DCDWriter: unknown unit {0!r}".format(unit))
        super(DCDWriter, self).__init__(*args, **kwargs)


class DCDReader(DCD.DCDReader):
    """Read a LAMMPS_ DCD trajectory.

    The units can be set from the constructor with the keyword
    arguments *timeunit* and *lengthunit*. The defaults are "fs" and
    "Angstrom", corresponding to LAMMPS `units style`_ "**real**". See
    :mod:`MDAnalysis.units` for other recognized values.

    .. _units style: http://lammps.sandia.gov/doc/units.html
    """
    format = 'LAMMPS'
    flavor = 'LAMMPS'

    def __init__(self, dcdfilename, **kwargs):
        self.units = {'time': 'fs', 'length': 'Angstrom'}  # must be instance level
        self.units['time'] = kwargs.pop('timeunit', self.units['time'])
        self.units['length'] = kwargs.pop('lengthunit', self.units['length'])
        for unit_type, unit in self.units.items():
            try:
                if units.unit_types[unit] != unit_type:
                    raise TypeError("LAMMPS DCDReader: wrong unit {0!r} for unit type {1!r}".format(unit, unit_type))
            except KeyError:
                raise ValueError("LAMMPS DCDReader: unknown unit {0!r}".format(unit))
        super(DCDReader, self).__init__(dcdfilename, **kwargs)


class DATAReader(base.SingleFrameReaderBase):
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
        self.atom_style = kwargs.pop('atom_style', None)
        super(DATAReader, self).__init__(filename, **kwargs)

    def _read_first_frame(self):
        with DATAParser(self.filename) as p:
            self.ts = p.read_DATA_timestep(self.n_atoms, self._Timestep,
                                           self._ts_kwargs, self.atom_style)

        self.ts.frame = 0
        if self.convert_units:
            self.convert_pos_from_native(self.ts._pos)  # in-place !
            try:
                self.convert_velocities_from_native(self.ts._velocities)  # in-place !
            except AttributeError:
                pass

class DATAWriter(base.WriterBase):
    """Write out the current time step as a LAMMPS DATA file.

    This writer supports the sections Atoms, Masses, Velocities, Bonds,
    Angles, Dihedrals, and Impropers. This writer will write the header
    and these sections (if applicable). Atoms section is written in the
    "full" sub-style if charges are available or "molecular" sub-style
    if they are not. Molecule id is set to 0 for all atoms.

    Note
    ----
    This writer assumes "conventional" or "real" LAMMPS units where length
    is measured in Angstroms and velocity is measured in Angstroms per
    femtosecond. To write in different units, specify `lengthunit`

    If atom types are not already positive integers, the user must set them
    to be positive integers, because the writer will not automatically
    assign new types.

    To preserve numerical atom types when writing a selection, the Masses
    section will have entries for each atom type up to the maximum atom type.
    If the universe does not contain atoms of some type in
    {1, ... max(atom_types)}, then the mass for that type will be set to 1.

    In order to write bonds, each selected bond type must be explicitly set to
    an integer >= 1.

    """
    format = 'DATA'

    def __init__(self, filename, convert_units=None, **kwargs):
        """Set up a DATAWriter

        Parameters
        ----------
        filename : str
            output filename
        """
        self.filename = util.filename(filename, ext='data')

        if convert_units is None:
            convert_units = flags['convert_lengths']
        self.convert_units = convert_units

        self.units = {'time': 'fs', 'length': 'Angstrom'}
        self.units['length'] = kwargs.pop('lengthunit', self.units['length'])
        self.units['time'] = kwargs.pop('timeunit', self.units['time'])
        self.units['velocity'] = kwargs.pop('velocityunit',
                                 self.units['length']+'/'+self.units['time'])

    def _write_atoms(self, atoms):
        self.f.write('\n')
        self.f.write('Atoms\n')
        self.f.write('\n')

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
                self.f.write('{i:d} 0 {t:d} {c:f} {x:f} {y:f} {z:f}\n'.format(
                             i=index, t=atype, c=charge, x=coords[0],
                             y=coords[1], z=coords[2]))
        else:
            for index, atype, coords in zip(indices, types, coordinates):
                self.f.write('{i:d} 0 {t:d} {x:f} {y:f} {z:f}\n'.format(
                             i=index, t=atype, x=coords[0], y=coords[1],
                             z=coords[2]))

    def _write_velocities(self, atoms):
        self.f.write('\n')
        self.f.write('Velocities\n')
        self.f.write('\n')
        indices = atoms.indices + 1
        velocities = self.convert_velocities_to_native(atoms.velocities,
                                                       inplace=False)
        for index, vel in zip(indices, velocities):
            self.f.write('{i:d} {x:f} {y:f} {z:f}\n'.format(i=index, x=vel[0],
                y=vel[1], z=vel[2]))

    def _write_masses(self, atoms):
        self.f.write('\n')
        self.f.write('Masses\n')
        self.f.write('\n')
        mass_dict = {}
        max_type = max(atoms.types.astype(np.int32))
        for atype in range(1, max_type+1):
            # search entire universe for mass info, not just writing selection
            masses = set(atoms.universe.atoms.select_atoms(
                'type {:d}'.format(atype)).masses)
            if len(masses) == 0:
                mass_dict[atype] = 1.0
            else:
                mass_dict[atype] = masses.pop()
            if masses:
                raise ValueError('LAMMPS DATAWriter: to write data file, '+
                        'atoms with same type must have same mass')
        for atype, mass in mass_dict.items():
            self.f.write('{:d} {:f}\n'.format(atype, mass))

    def _write_bonds(self, bonds):
        self.f.write('\n')
        self.f.write('{}\n'.format(btype_sections[bonds.btype]))
        self.f.write('\n')
        for bond, i in zip(bonds, range(1, len(bonds)+1)):
            try:
                self.f.write('{:d} {:d} '.format(i, int(bond.type))+\
                        ' '.join((bond.atoms.indices + 1).astype(str))+'\n')
            except TypeError:
                raise TypeError('LAMMPS DATAWriter: Trying to write bond, '
                                'but bond type {} is not '
                                'numerical.'.format(bond.type))

    def _write_dimensions(self, dimensions):
        """Convert dimensions to triclinic vectors, convert lengths to native
        units and then write the dimensions section
        """
        if self.convert_units:
            triv = self.convert_pos_to_native(mdamath.triclinic_vectors(
                                              dimensions),inplace=False)
        self.f.write('\n')
        self.f.write('{:f} {:f} xlo xhi\n'.format(0., triv[0][0]))
        self.f.write('{:f} {:f} ylo yhi\n'.format(0., triv[1][1]))
        self.f.write('{:f} {:f} zlo zhi\n'.format(0., triv[2][2]))
        if any([triv[1][0], triv[2][0], triv[2][1]]):
            self.f.write('{xy:f} {xz:f} {yz:f} xy xz yz\n'.format(
                xy=triv[1][0], xz=triv[2][0], yz=triv[2][1]))
        self.f.write('\n')

    @requires('types', 'masses')
    def write(self, selection, frame=None):
        """Write selection at current trajectory frame to file.

        The sections for Atoms, Masses, Velocities, Bonds, Angles,
        Dihedrals, and Impropers (if these are defined) are
        written. The Atoms section is written in the "full" sub-style
        if charges are available or "molecular" sub-style if they are
        not. Molecule id in atoms section is set to to 0.

        No other sections are written to the DATA file.
        As of this writing, other sections are not parsed into the topology
        by the :class:`DATAReader`.

        Note
        ----
        If the selection includes a partial fragment, then only the bonds,
        angles, etc. whose atoms are contained within the selection will be
        included.

        Parameters
        ----------
        selection : AtomGroup or Universe
            MDAnalysis AtomGroup (selection or Universe.atoms) or also Universe
        frame : int (optional)
            optionally move to frame number `frame`

        """
        u = selection.universe
        if frame is not None:
            u.trajectory[frame]
        else:
            frame = u.trajectory.ts.frame

        # make sure to use atoms (Issue 46)
        atoms = selection.atoms

        # check that types can be converted to ints if they aren't ints already
        try:
            atoms.types.astype(np.int32)
        except ValueError:
            raise ValueError('LAMMPS.DATAWriter: atom types must be '+
                    'convertible to integers')

        try:
            velocities = atoms.velocities
        except (NoDataError, AttributeError):
            has_velocities = False
        else:
            has_velocities = True

        features = {}
        with util.openany(self.filename, 'w') as self.f:
            self.f.write('LAMMPS data file via MDAnalysis\n')
            self.f.write('\n')
            self.f.write('{:>12d}  atoms\n'.format(len(atoms)))

            attrs = [('bond', 'bonds'), ('angle', 'angles'),
                ('dihedral', 'dihedrals'), ('improper', 'impropers')]

            for btype, attr_name in attrs:
                try:
                    features[btype] = atoms.__getattribute__(attr_name)
                    self.f.write('{:>12d}  {}\n'.format(len(features[btype]),
                                 attr_name))
                    features[btype] = features[btype].atomgroup_intersection(
                                      atoms, strict=True)
                except AttributeError:
                    features[btype] = None
                    self.f.write('{:>12d}  {}\n'.format(0, attr_name))

            self.f.write('\n')
            self.f.write('{:>12d}  atom types\n'.format(max(atoms.types.astype(np.int32))))

            for btype, attr in features.items():
                if attr is None:
                    self.f.write('{:>12d}  {} types\n'.format(0, btype))
                else:
                    self.f.write('{:>12d}  {} types\n'.format(len(attr.types()),
                                                              btype))

            self._write_dimensions(atoms.dimensions)

            self._write_masses(atoms)
            self._write_atoms(atoms)
            for attr in features.values():
                if attr is None or len(attr) == 0:
                    continue
                self._write_bonds(attr)

            if has_velocities:
                self._write_velocities(atoms)
