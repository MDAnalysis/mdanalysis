# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 
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

"""
GRO file format --- :mod:`MDAnalysis.coordinates.GRO`
======================================================

Classes to read and write Gromacs_ GRO_ coordinate files; see the notes on the
`GRO format`_ which includes a conversion routine for the box.

GROWriter format strings
------------------------

The GROWriter class has a .fmt attribute, which is a dictionary of different
strings for writing lines in .gro files.  These are as follows:

n_atoms
  For the first line of the gro file, supply the number of atoms in the system.
  Eg: fmt['n_atoms'].format(42)

xyz
  An atom line without velocities.  Requires that the 'resid', 'resname',
  'name', 'index' and 'pos' keys be supplied.
  Eg: fmt['xyz'].format(resid=1, resname='SOL', name='OW2', index=2,
      pos=(0.0, 1.0, 2.0))

xyz_v
  As above, but with velocities.  Needs an additional keyword 'vel'.

box_orthorhombic
  The final line of the gro file which gives box dimensions.  Requires
  the box keyword to be given, which should be the three cartesian dimensions.
  Eg: fmt['box_orthorhombic'].format(box=(10.0, 10.0, 10.0))

box_triclinic
  As above, but for a non orthorhombic box. Requires the box keyword, but this
  time as a length 9 vector.  This is a flattened version of the (3,3) triclinic
  vector representation of the unit cell.  The rearrangement into the odd
  gromacs order is done automatically.


.. _Gromacs: http://www.gromacs.org
.. _GRO: http://manual.gromacs.org/current/online/gro.html
.. _GRO format: http://chembytes.wikidot.com/g-grofile
"""

from six.moves import range
import warnings
import numpy as np

from ..core import flags
from . import base
from ..lib import util
from .core import triclinic_box, triclinic_vectors
from ..exceptions import NoDataError


class Timestep(base.Timestep):
    _ts_order_x = [0, 3, 4]
    _ts_order_y = [5, 1, 6]
    _ts_order_z = [7, 8, 2]

    def _init_unitcell(self):
        return np.zeros(9, dtype=np.float32)

    @property
    def dimensions(self):
        """unitcell dimensions (A, B, C, alpha, beta, gamma)

        GRO::

          8.00170   8.00170   5.65806   0.00000   0.00000   0.00000   0.00000   4.00085   4.00085

        PDB::

          CRYST1   80.017   80.017   80.017  60.00  60.00  90.00 P 1           1

        XTC: c.trajectory.ts._unitcell::

          array([[ 80.00515747,   0.        ,   0.        ],
                 [  0.        ,  80.00515747,   0.        ],
                 [ 40.00257874,  40.00257874,  56.57218552]], dtype=float32)
        """
        # unit cell line (from http://manual.gromacs.org/current/online/gro.html)
        # v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)
        # 0     1     2      3     4     5     6    7     8
        # This information now stored as _ts_order_x/y/z to keep DRY
        x = self._unitcell[self._ts_order_x]
        y = self._unitcell[self._ts_order_y]
        z = self._unitcell[self._ts_order_z]  # this ordering is correct! (checked it, OB)
        return triclinic_box(x, y, z)

    @dimensions.setter
    def dimensions(self, box):
        x, y, z = triclinic_vectors(box)
        np.put(self._unitcell, self._ts_order_x, x)
        np.put(self._unitcell, self._ts_order_y, y)
        np.put(self._unitcell, self._ts_order_z, z)


class GROReader(base.SingleFrameReader):
    """Reader for the Gromacs GRO structure format.
    
    .. versionchanged:: 0.11.0
       Frames now 0-based instead of 1-based
    """
    format = 'GRO'
    units = {'time': None, 'length': 'nm', 'velocity': 'nm/ps'}
    _Timestep = Timestep

    def _read_first_frame(self):
        with util.openany(self.filename, 'rt') as grofile:
            # Read first two lines to get number of atoms
            grofile.readline()
            self.n_atoms = n_atoms = int(grofile.readline())
            # and the third line to get the spacing between coords (cs)
            # (dependent upon the GRO file precision)
            first_atomline = grofile.readline()
            cs = first_atomline[25:].find('.') + 1

            # Always try, and maybe add them later
            velocities = np.zeros((n_atoms, 3), dtype=np.float32)

            self.ts = ts = self._Timestep(n_atoms,
                                          **self._ts_kwargs)

            missed_vel = False

            grofile.seek(0)
            for pos, line in enumerate(grofile, start=-2):
                # 2 header lines, 1 box line at end
                if pos == n_atoms:
                    unitcell = np.array(list(map(float, line.split())))
                    continue
                if pos < 0:
                    continue

                ts._pos[pos] = [line[20 + cs*i:20 + cs*(i+1)] for i in range(3)]
                try:
                    velocities[pos] = [line[20 + cs*i:20 + cs*(i+1)] for i in range(3, 6)]
                except ValueError:
                    # Remember that we got this error
                    missed_vel = True

        if np.any(velocities):
            ts.velocities = velocities
            if missed_vel:
                warnings.warn("Not all velocities were present.  "
                              "Unset velocities set to zero.")

        self.ts.frame = 0  # 0-based frame number

        if len(unitcell) == 3:
            # special case: a b c --> (a 0 0) (b 0 0) (c 0 0)
            # see Timestep.dimensions() above for format (!)
            self.ts._unitcell[:3] = unitcell
        elif len(unitcell) == 9:
            self.ts._unitcell[:] = unitcell  # fill all
        else:  # or maybe raise an error for wrong format??
            warnings.warn("GRO unitcell has neither 3 nor 9 entries --- might be wrong.")
            self.ts._unitcell[:len(unitcell)] = unitcell  # fill linearly ... not sure about this
        if self.convert_units:
            self.convert_pos_from_native(self.ts._pos)  # in-place !
            self.convert_pos_from_native(self.ts._unitcell)  # in-place ! (all are lengths)
            if self.ts.has_velocities:
                # converts nm/ps to A/ps units
                self.convert_velocities_from_native(self.ts._velocities)

    def Writer(self, filename, **kwargs):
        """Returns a CRDWriter for *filename*.

        :Arguments:
          *filename*
            filename of the output GRO file

        :Returns: :class:`GROWriter`

        """
        return GROWriter(filename, **kwargs)


class GROWriter(base.Writer):
    """GRO Writer that conforms to the Trajectory API.

    .. Note::

       The precision is hard coded to three decimal places and
       velocities are not written (yet).

    .. versionchanged:: 0.11.0
       Frames now 0-based instead of 1-based
    .. versionchanged:: 0.13.0
       Now strictly writes positions with 3dp precision
       and velocities with 4dp
       Removed the `convert_dimensions_to_unitcell` method,
       use `Timestep.triclinic_dimensions` instead
       Now now writes velocities where possible
    """

    format = 'GRO'
    units = {'time': None, 'length': 'nm', 'velocity': 'nm/ps'}
    gro_coor_limits = {'min': -999.9995, 'max': 9999.9995}

    #: format strings for the GRO file (all include newline); precision
    #: of 3 decimal places is hard-coded here.
    fmt = {
        'n_atoms': "{0:5d}\n",  # number of atoms
        # coordinates output format, see http://chembytes.wikidot.com/g-grofile
        'xyz': "{resid:>5d}{resname:<5.5s}{name:>5.5s}{index:>5d}{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}\n",
        # unitcell
        'box_orthorhombic': "{box[0]:10.5f}{box[1]:10.5f}{box[2]:10.5f}\n",
        'box_triclinic': "{box[0]:10.5f}{box[4]:10.5f}{box[8]:10.5f}{box[1]:10.5f}{box[2]:10.5f}{box[3]:10.5f}{box[5]:10.5f}{box[6]:10.5f}{box[7]:10.5f}\n"
    }
    fmt['xyz_v'] = fmt['xyz'][:-1] + "{vel[0]:8.4f}{vel[1]:8.4f}{vel[2]:8.4f}\n"

    def __init__(self, filename, convert_units=None, **kwargs):
        """Set up a GROWriter with a precision of 3 decimal places.

        :Arguments:
           *filename*
              output filename
        """
        self.filename = util.filename(filename, ext='gro')

        if convert_units is None:
            convert_units = flags['convert_lengths']
        self.convert_units = convert_units  # convert length and time to base units

    def write(self, selection, frame=None):
        """Write selection at current trajectory frame to file.

        :Arguments:
          selection
              MDAnalysis AtomGroup (selection or Universe.atoms)
              or also Universe
        :Keywords:
          frame
              optionally move to frame number *frame*

        The GRO format only allows 5 digits for resid and atom
        number. If these number become larger than 99,999 then this
        routine will chop off the leading digits.

        .. versionchanged:: 0.7.6
           resName and atomName are truncated to a maximum of 5 characters
        """
        # write() method that complies with the Trajectory API
        u = selection.universe
        if frame is not None:
            u.trajectory[frame]  # advance to frame
        else:
            frame = u.trajectory.ts.frame

        # make sure to use atoms (Issue 46)
        atoms = selection.atoms
        # can write from selection == Universe (Issue 49)
        coordinates = atoms.positions

        try:
            velocities = atoms.velocities
        except NoDataError:
            has_velocities = False
        else:
            has_velocities = True

        if self.convert_units:
            # Convert back to nm from Angstroms,
            # inplace because coordinates is already a copy
            self.convert_pos_to_native(coordinates)
            if has_velocities:
               self.convert_velocities_to_native(velocities) 
        # check if any coordinates are illegal
        # (checks the coordinates in native nm!)
        if not self.has_valid_coordinates(self.gro_coor_limits, coordinates):
            raise ValueError("GRO files must have coordinate values between "
                             "{0:.3f} and {1:.3f} nm: No file was written."
                             "".format(self.gro_coor_limits["min"],
                                       self.gro_coor_limits["max"]))

        with util.openany(self.filename, 'wt') as output_gro:
            # Header
            output_gro.write('Written by MDAnalysis\n')
            output_gro.write(self.fmt['n_atoms'].format(len(atoms)))
            # Atom descriptions and coords
            for atom_index, atom in enumerate(atoms):
                truncated_atom_index = int(str(atom_index + 1)[-5:])
                if has_velocities:
                    output_gro.write(self.fmt['xyz_v'].format(
                        resid=atom.resid,
                        resname=atom.resname,
                        index=truncated_atom_index,
                        name=atom.name,
                        pos=coordinates[atom_index],
                        vel=velocities[atom_index],
                    ))
                else:
                    output_gro.write(self.fmt['xyz'].format(
                        resid=atom.resid,
                        resname=atom.resname,
                        index=truncated_atom_index,
                        name=atom.name,
                        pos=coordinates[atom_index]
                    ))

            # Footer: box dimensions
            if np.all(u.trajectory.ts.dimensions[3:] == [90., 90., 90.]):
                box = self.convert_pos_to_native(
                    u.coord.dimensions[:3], inplace=False)
                # orthorhombic cell, only lengths along axes needed in gro
                output_gro.write(self.fmt['box_orthorhombic'].format(
                    box=box)
                )
            else:
                # full output
                box = self.convert_pos_to_native(
                    u.coord.triclinic_dimensions.flatten(), inplace=False)
                output_gro.write(self.fmt['box_triclinic'].format(
                    box=box)
                )
