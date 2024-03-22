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
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#


"""LAMMPS DCD trajectory and DATA I/O  --- :mod:`MDAnalysis.coordinates.LAMMPS`
===============================================================================

Classes to read and write LAMMPS_ DCD binary trajectories, LAMMPS DATA files
and LAMMPS dump files.  Trajectories can be read regardless of system-endianness
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

Dump files
----------

The DumpReader expects ascii dump files written with the default
`LAMMPS dump format`_ of 'atom'


Example: Loading a LAMMPS simulation
------------------------------------
.. testsetup::

To load a LAMMPS simulation from a LAMMPS data file (using the
:class:`~MDAnalysis.topology.LAMMPSParser.DATAParser`) together with a
LAMMPS DCD with "*real*" provide the keyword *format="LAMMPS*"::

>>> import MDAnalysis
>>> from MDAnalysis.tests.datafiles import LAMMPSdata2, LAMMPSdcd2
>>> u = MDAnalysis.Universe(LAMMPSdata2, LAMMPSdcd2, format="LAMMPS")

If the trajectory uses *units nano* then use

>>> import MDAnalysis
>>> from MDAnalysis.tests.datafiles import LAMMPSdata2, LAMMPSdcd2
>>> u = MDAnalysis.Universe(LAMMPSdata2, LAMMPSdcd2, format="LAMMPS",
...                          lengthunit="nm", timeunit="ns")

To scan through a trajectory to find a desirable frame and write to a LAMMPS
data file,

   >>> import MDAnalysis
   >>> from MDAnalysis.tests.datafiles import LAMMPSdata2, LAMMPSdcd2 
   >>> u = MDAnalysis.Universe(LAMMPSdata2, LAMMPSdcd2, format="LAMMPS",
   ...                          lengthunit="nm", timeunit="ns")
   >>> take_this_frame = False
   >>> for ts in u.trajectory:
   ...     # analyze frame
   ...     if ts.frame == 4:
   ...         take_this_frame = True
   ...     if take_this_frame == True:
   ...         with MDAnalysis.Writer('frame.data') as W:
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
.. _`LAMMPS dump format`: http://lammps.sandia.gov/doc/dump.html

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
.. autoclass:: DumpReader
   :members:
   :inherited-members:

"""
import os
import numpy as np

from ..core.groups import requires
from ..lib import util, mdamath, distances
from ..lib.util import cached, store_init_arguments
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
                errmsg = f"LAMMPS DCDWriter: unknown unit {unit}"
                raise ValueError(errmsg) from None
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

    @store_init_arguments
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

    @store_init_arguments
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

    def __init__(self, filename, convert_units=True, **kwargs):
        """Set up a DATAWriter

        Parameters
        ----------
        filename : str
            output filename
        convert_units : bool, optional
            units are converted to the MDAnalysis base format; [``True``]
        """
        self.filename = util.filename(filename, ext='data', keep=True)

        self.convert_units = convert_units

        self.units = {'time': 'fs', 'length': 'Angstrom'}
        self.units['length'] = kwargs.pop('lengthunit', self.units['length'])
        self.units['time'] = kwargs.pop('timeunit', self.units['time'])
        self.units['velocity'] = kwargs.pop('velocityunit',
                                 self.units['length']+'/'+self.units['time'])

    def _write_atoms(self, atoms, data):
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

        moltags = data.get("molecule_tag", np.zeros(len(atoms), dtype=int))

        if self.convert_units:
            coordinates = self.convert_pos_to_native(atoms.positions, inplace=False)

        if has_charges:
            for index, moltag, atype, charge, coords in zip(indices, moltags,
                    types, charges, coordinates):
                x, y, z = coords
                self.f.write(f"{index:d} {moltag:d} {atype:d} {charge:f}"
                             f" {x:f} {y:f} {z:f}\n")
        else:
            for index, moltag, atype, coords in zip(indices, moltags, types,
                    coordinates):
                x, y, z = coords
                self.f.write(f"{index:d} {moltag:d} {atype:d}"
                             f" {x:f} {y:f} {z:f}\n")

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
                errmsg = (f"LAMMPS DATAWriter: Trying to write bond, but bond "
                          f"type {bond.type} is not numerical.")
                raise TypeError(errmsg) from None

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
            errmsg = ("LAMMPS.DATAWriter: atom types must be convertible to "
                      "integers")
            raise ValueError(errmsg) from None

        try:
            velocities = atoms.velocities
        except (NoDataError, AttributeError):
            has_velocities = False
        else:
            has_velocities = True

        features = {}
        with util.openany(self.filename, 'wt') as self.f:
            self.f.write('LAMMPS data file via MDAnalysis\n')
            self.f.write('\n')
            self.f.write('{:>12d}  atoms\n'.format(len(atoms)))

            attrs = [('bond', 'bonds'), ('angle', 'angles'),
                ('dihedral', 'dihedrals'), ('improper', 'impropers')]

            for btype, attr_name in attrs:
                features[btype] = atoms.__getattribute__(attr_name)
                self.f.write('{:>12d}  {}\n'.format(len(features[btype]),
                                                    attr_name))
                features[btype] = features[btype].atomgroup_intersection(
                                    atoms, strict=True)

            self.f.write('\n')
            self.f.write('{:>12d}  atom types\n'.format(max(atoms.types.astype(np.int32))))

            for btype, attr in features.items():
                self.f.write('{:>12d}  {} types\n'.format(len(attr.types()),
                                                          btype))

            self._write_dimensions(atoms.dimensions)

            self._write_masses(atoms)
            self._write_atoms(atoms, u.trajectory.ts.data)
            for attr in features.values():
                if attr is None or len(attr) == 0:
                    continue
                self._write_bonds(attr)

            if has_velocities:
                self._write_velocities(atoms)


class DumpReader(base.ReaderBase):
    """Reads the default `LAMMPS dump format 
    <https://docs.lammps.org/dump.html>`__

    Supports coordinates in various LAMMPS coordinate conventions and both 
    orthogonal and triclinic simulation box dimensions (for more details see 
    `documentation <https://docs.lammps.org/Howto_triclinic.html>`__). In 
    either case, MDAnalysis will always use ``(*A*, *B*, *C*, *alpha*, *beta*,
    *gamma*)`` to represent the unit cell. Lengths *A*, *B*, *C* are in the 
    MDAnalysis length unit (Å), and angles are in degrees.

    Parameters
    ----------
    filename : str
        Filename of LAMMPS dump file
    lammps_coordinate_convention : str (optional) default="auto"
        Convention used in coordinates, can be one of the following according
        to the `LAMMPS documentation <https://docs.lammps.org/dump.html>`__:

         - "auto" - Detect coordinate type from file column header. If auto 
           detection is used, the guessing checks whether the coordinates
           fit each convention in the order "unscaled", "scaled", "unwrapped", 
           "scaled_unwrapped" and whichever set of coordinates is detected 
           first will be used.
         - "scaled" - Coordinates wrapped in box and scaled by box length (see
            note below), i.e., xs, ys, zs
         - "scaled_unwrapped" - Coordinates unwrapped and scaled by box length,
           (see note below) i.e., xsu, ysu, zsu
         - "unscaled" - Coordinates wrapped in box, i.e., x, y, z
         - "unwrapped" - Coordinates unwrapped, i.e., xu, yu, zu

        If coordinates are given in the scaled coordinate convention (xs,ys,zs)
        or scaled unwrapped coordinate convention (xsu,ysu,zsu) they will 
        automatically be converted from their scaled/fractional representation
        to their real values.
    unwrap_images : bool (optional) default=False
        If `True` and the dump file contains image flags, the coordinates 
        will be unwrapped. See `read_data 
        <https://docs.lammps.org/read_data.html>`__  in the lammps 
        documentation for more information.
    **kwargs
       Other keyword arguments used in :class:`~MDAnalysis.coordinates.base.ReaderBase`

    .. versionchanged:: 2.4.0
       Now imports velocities and forces, translates the box to the origin,
       and optionally unwraps trajectories with image flags upon loading.
    .. versionchanged:: 2.2.0
       Triclinic simulation boxes are supported.
       (Issue `#3383 <https://github.com/MDAnalysis/mdanalysis/issues/3383>`__)
    .. versionchanged:: 2.0.0
       Now parses coordinates in multiple lammps conventions (x,xs,xu,xsu)
    .. versionadded:: 0.19.0
    """
    format = 'LAMMPSDUMP'
    _conventions = ["auto", "unscaled", "scaled", "unwrapped",
                    "scaled_unwrapped"]
    _coordtype_column_names = {
        "unscaled": ["x", "y", "z"],
        "scaled": ["xs", "ys", "zs"],
        "unwrapped": ["xu", "yu", "zu"],
        "scaled_unwrapped": ["xsu", "ysu", "zsu"]
    }

    @store_init_arguments
    def __init__(self, filename, 
                 lammps_coordinate_convention="auto",
                 unwrap_images=False,
                 **kwargs):
        super(DumpReader, self).__init__(filename, **kwargs)

        root, ext = os.path.splitext(self.filename)
        if lammps_coordinate_convention in self._conventions:
            self.lammps_coordinate_convention = lammps_coordinate_convention
        else:
            option_string = "'" + "', '".join(self._conventions) + "'"
            raise ValueError("lammps_coordinate_convention="
                             f"'{lammps_coordinate_convention}'"
                             " is not a valid option. "
                             f"Please choose one of {option_string}")

        self._unwrap = unwrap_images

        self._cache = {}

        self._reopen()

        self._read_next_timestep()

    def _reopen(self):
        self.close()
        self._file = util.anyopen(self.filename)
        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
        self.ts.frame = -1

    @property
    @cached('n_atoms')
    def n_atoms(self):
        with util.anyopen(self.filename) as f:
            f.readline()
            f.readline()
            f.readline()
            n_atoms = int(f.readline())
        return n_atoms

    @property
    @cached('n_frames')
    def n_frames(self):
        # 2(timestep) + 2(natoms info) + 4(box info) + 1(atom header) + n_atoms
        lines_per_frame = self.n_atoms + 9
        offsets = []
        counter = 0
        with util.anyopen(self.filename) as f:
            line = True
            while line:
                if not counter % lines_per_frame:
                    offsets.append(f.tell())
                line = f.readline()
                counter += 1
        self._offsets = offsets[:-1]  # last is EOF
        return len(self._offsets)

    def close(self):
        if hasattr(self, '_file'):
            self._file.close()

    def _read_frame(self, frame):
        self._file.seek(self._offsets[frame])
        self.ts.frame = frame - 1  # gets +1'd in next

        return self._read_next_timestep()

    def _read_next_timestep(self):
        f = self._file
        ts = self.ts
        ts.frame += 1
        if ts.frame >= len(self):
            raise EOFError

        f.readline()  # ITEM TIMESTEP
        step_num = int(f.readline())
        ts.data['step'] = step_num
        ts.data['time'] = step_num * ts.dt

        f.readline()  # ITEM NUMBER OF ATOMS
        n_atoms = int(f.readline())
        if n_atoms != self.n_atoms:
            raise ValueError("Number of atoms in trajectory changed "
                             "this is not supported in MDAnalysis")

        triclinic = len(f.readline().split()) == 9  # ITEM BOX BOUNDS
        if triclinic:
            xlo_bound, xhi_bound, xy = map(float, f.readline().split())
            ylo_bound, yhi_bound, xz = map(float, f.readline().split())
            zlo, zhi, yz = map(float, f.readline().split())

            # converts orthogonal bounding box to the conventional format,
            # see https://docs.lammps.org/Howto_triclinic.html
            xlo = xlo_bound - min(0.0, xy, xz, xy + xz)
            xhi = xhi_bound - max(0.0, xy, xz, xy + xz)
            ylo = ylo_bound - min(0.0, yz)
            yhi = yhi_bound - max(0.0, yz)

            box = np.zeros((3, 3), dtype=np.float64)
            box[0] = xhi - xlo, 0.0, 0.0
            box[1] = xy, yhi - ylo, 0.0
            box[2] = xz, yz, zhi - zlo

            xlen, ylen, zlen, alpha, beta, gamma = mdamath.triclinic_box(*box)
        else:
            xlo, xhi = map(float, f.readline().split())
            ylo, yhi = map(float, f.readline().split())
            zlo, zhi = map(float, f.readline().split())
            xlen = xhi - xlo
            ylen = yhi - ylo
            zlen = zhi - zlo
            alpha = beta = gamma = 90.
        ts.dimensions = xlen, ylen, zlen, alpha, beta, gamma

        indices = np.zeros(self.n_atoms, dtype=int)

        atomline = f.readline()  # ITEM ATOMS etc
        attrs = atomline.split()[2:]  # attributes on coordinate line
        attr_to_col_ix = {x: i for i, x in enumerate(attrs)}
        convention_to_col_ix = {}
        for cv_name, cv_col_names in self._coordtype_column_names.items():
            try:
                convention_to_col_ix[cv_name] = [attr_to_col_ix[x] 
                    for x in cv_col_names]
            except KeyError:
                pass

        if self._unwrap:
            try:
                image_cols = [attr_to_col_ix[x] for x in ["ix", "iy", "iz"]]
            except:
                raise ValueError("Trajectory must have image flag in order "
                                 "to unwrap.")

        self._has_vels = all(x in attr_to_col_ix for x in ["vx", "vy", "vz"])
        if self._has_vels:
            ts.has_velocities = True
            vel_cols = [attr_to_col_ix[x] for x in ["vx", "vy", "vz"]]
        self._has_forces = all(x in attr_to_col_ix for x in ["fx", "fy", "fz"])
        if self._has_forces:
            ts.has_forces = True
            force_cols = [attr_to_col_ix[x] for x in ["fx", "fy", "fz"]]

        # this should only trigger on first read of "ATOM" card, after which it
        # is fixed to the guessed value. Auto proceeds unscaled -> scaled
        # -> unwrapped -> scaled_unwrapped
        if self.lammps_coordinate_convention == "auto":
            try:
                # this will automatically select in order of priority
                # unscaled, scaled, unwrapped, scaled_unwrapped
                self.lammps_coordinate_convention = list(convention_to_col_ix)[0]
            except IndexError:
                raise ValueError("No coordinate information detected")
        elif not self.lammps_coordinate_convention in convention_to_col_ix:
            raise ValueError(f"No coordinates following convention "
                             "{self.lammps_coordinate_convention} found in "
                             "timestep")

        coord_cols = convention_to_col_ix[self.lammps_coordinate_convention]
        if self._unwrap:
            coord_cols.extend(image_cols)

        ids = "id" in attr_to_col_ix
        for i in range(self.n_atoms):
            fields = f.readline().split()
            if ids:
                indices[i] = fields[attr_to_col_ix["id"]]
            coords = np.array([fields[dim] for dim in coord_cols], 
                              dtype=np.float32)

            if self._unwrap:
                images = coords[3:]
                coords = coords[:3]
                coords += images * ts.dimensions[:3]
            else:
                coords = coords[:3]
            ts.positions[i] = coords

            if self._has_vels:
                ts.velocities[i] = [fields[dim] for dim in vel_cols]
            if self._has_forces:
                ts.forces[i] = [fields[dim] for dim in force_cols]

        order = np.argsort(indices)
        ts.positions = ts.positions[order]
        if self._has_vels:
            ts.velocities = ts.velocities[order]
        if self._has_forces:
            ts.forces = ts.forces[order]
        if (self.lammps_coordinate_convention.startswith("scaled")):
            # if coordinates are given in scaled format, undo that
            ts.positions = distances.transform_StoR(ts.positions,
                                                    ts.dimensions)
        # Transform to origin after transformation of scaled variables
        ts.positions -= np.array([xlo, ylo, zlo])[None,:]

        return ts
