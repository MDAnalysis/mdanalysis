# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
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

"""
GRO file format --- :mod:`MDAnalysis.coordinates.GRO`
======================================================

Classes to read and write Gromacs_ GRO_ coordinate files; see the notes on the
`GRO format`_ which includes a conversion routine for the box.


Writing GRO files
-----------------

By default any written GRO files will renumber the atom ids to move sequentially
from 1.  This can be disabled, and instead the original atom ids kept, by
using the `reindex=False` keyword argument.  This is useful when writing a
subsection of a larger Universe while wanting to preserve the original
identities of atoms.

For example::

   >>> u = mda.Universe()`

   >>> u.atoms.write('out.gro', reindex=False)

   # OR
   >>> with mda.Writer('out.gro', reindex=False) as w:
   ...     w.write(u.atoms)


Classes
-------

.. autoclass:: GROReader
   :members:

.. autoclass:: GROWriter
   :members:


Developer notes: ``GROWriter`` format strings
---------------------------------------------

The :class:`GROWriter` class has a :attr:`GROWriter.fmt` attribute, which is a dictionary of different
strings for writing lines in ``.gro`` files.  These are as follows:

``n_atoms``
  For the first line of the gro file, supply the number of atoms in the system.
  E.g.::

      fmt['n_atoms'].format(42)

``xyz``
  An atom line without velocities.  Requires that the 'resid', 'resname',
  'name', 'index' and 'pos' keys be supplied.
  E.g.::

     fmt['xyz'].format(resid=1, resname='SOL', name='OW2', index=2, pos=(0.0, 1.0, 2.0))

``xyz_v``
  As above, but with velocities.  Needs an additional keyword 'vel'.

``box_orthorhombic``
  The final line of the gro file which gives box dimensions.  Requires
  the box keyword to be given, which should be the three cartesian dimensions.
  E.g.::

     fmt['box_orthorhombic'].format(box=(10.0, 10.0, 10.0))

``box_triclinic``
  As above, but for a non orthorhombic box. Requires the box keyword, but this
  time as a length 9 vector.  This is a flattened version of the (3,3) triclinic
  vector representation of the unit cell.  The rearrangement into the odd
  gromacs order is done automatically.


.. _Gromacs: http://www.gromacs.org
.. _GRO: https://manual.gromacs.org/current/reference-manual/file-formats.html#gro
.. _GRO format: http://chembytes.wikidot.com/g-grofile

"""
import re

import itertools
import warnings

import numpy as np

from . import base
from .core import triclinic_box, triclinic_vectors
from ..exceptions import NoDataError
from ..lib import util
from .timestep import Timestep

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
_TS_ORDER_X = [0, 3, 4]
_TS_ORDER_Y = [5, 1, 6]
_TS_ORDER_Z = [7, 8, 2]


def _gmx_to_dimensions(box):
    # convert gromacs ordered box to [lx, ly, lz, alpha, beta, gamma] form
    x = box[_TS_ORDER_X]
    y = box[_TS_ORDER_Y]
    z = box[_TS_ORDER_Z]  # this ordering is correct! (checked it, OB)
    return triclinic_box(x, y, z)


class GROReader(base.SingleFrameReaderBase):
    """Reader for the Gromacs GRO structure format.

    .. note::
       This Reader will only read the first frame present in a file.


    .. note::
       GRO files with zeroed 3 entry unit cells (i.e. ``0.0   0.0   0.0``)
       are read as missing unit cell information. In these cases ``dimensions``
       will be set to ``None``.

    .. versionchanged:: 0.11.0
       Frames now 0-based instead of 1-based
    .. versionchanged:: 2.0.0
       Reader now only parses boxes defined with 3 or 9 fields.
       Reader now reads a 3 entry zero unit cell (i.e. ``[0, 0, 0]``) as a
       being without dimension information (i.e. will the timestep dimension
       to ``None``).
    """

    format = "GRO"
    units = {"time": None, "length": "nm", "velocity": "nm/ps"}
    _Timestep = Timestep

    def _read_first_frame(self):
        with util.openany(self.filename, "rt") as grofile:
            # Read first two lines to get number of atoms
            grofile.readline()
            self.n_atoms = n_atoms = int(grofile.readline())
            self.ts = ts = self._Timestep(n_atoms, **self._ts_kwargs)
            # Always try, and maybe add them later
            velocities = np.zeros((n_atoms, 3), dtype=np.float32)
            missed_vel = False

            # and the third line to get the spacing between coords (cs)
            # (dependent upon the GRO file precision)
            first_atomline = grofile.readline()
            cs = first_atomline[25:].find(".") + 1
            ts._pos[0] = [
                first_atomline[20 + cs * i : 20 + cs * (i + 1)]
                for i in range(3)
            ]
            try:
                velocities[0] = [
                    first_atomline[20 + cs * i : 20 + cs * (i + 1)]
                    for i in range(3, 6)
                ]
            except ValueError:
                # Remember that we got this error
                missed_vel = True

            # TODO: Handle missing unitcell?
            for pos, line in enumerate(grofile, start=1):
                # 2 header lines, 1 box line at end
                if pos == n_atoms:
                    try:
                        unitcell = np.float32(line.split())
                    except ValueError:
                        # Try to parse floats with 5 digits if no spaces between values...
                        unitcell = np.float32(
                            re.findall(r"(\d+\.\d{5})", line)
                        )
                    break

                ts._pos[pos] = [
                    line[20 + cs * i : 20 + cs * (i + 1)] for i in range(3)
                ]
                try:
                    velocities[pos] = [
                        line[20 + cs * i : 20 + cs * (i + 1)]
                        for i in range(3, 6)
                    ]
                except ValueError:
                    # Remember that we got this error
                    missed_vel = True

        if np.any(velocities):
            ts.velocities = velocities
            if missed_vel:
                warnings.warn(
                    "Not all velocities were present.  "
                    "Unset velocities set to zero."
                )

        self.ts.frame = 0  # 0-based frame number

        if len(unitcell) == 3:
            # special case: a b c --> (a 0 0) (b 0 0) (c 0 0)
            # see docstring above for format (!)
            # Treat empty 3 entry boxes as not having a unit cell
            if np.allclose(unitcell, [0.0, 0.0, 0.0]):
                wmsg = (
                    "Empty box [0., 0., 0.] found - treating as missing "
                    "unit cell. Dimensions set to `None`."
                )
                warnings.warn(wmsg)
                self.ts.dimensions = None
            else:
                self.ts.dimensions = np.r_[unitcell, [90.0, 90.0, 90.0]]
        elif len(unitcell) == 9:
            self.ts.dimensions = _gmx_to_dimensions(unitcell)
        else:  # raise an error for wrong format
            errmsg = "GRO unitcell has neither 3 nor 9 entries."
            raise ValueError(errmsg)

        if self.convert_units:
            self.convert_pos_from_native(self.ts._pos)  # in-place !
            if self.ts.dimensions is not None:
                self.convert_pos_from_native(
                    self.ts.dimensions[:3]
                )  # in-place!
            if self.ts.has_velocities:
                # converts nm/ps to A/ps units
                self.convert_velocities_from_native(self.ts._velocities)

    def Writer(self, filename, n_atoms=None, **kwargs):
        """Returns a CRDWriter for *filename*.

        Parameters
        ----------
        filename: str
            filename of the output GRO file

        Returns
        -------
        :class:`GROWriter`

        """
        if n_atoms is None:
            n_atoms = self.n_atoms
        return GROWriter(filename, n_atoms=n_atoms, **kwargs)


class GROWriter(base.WriterBase):
    """GRO Writer that conforms to the Trajectory API.

    Will attempt to write the following information from the topology:
     - atom name (defaults to 'X')
     - resnames (defaults to 'UNK')
     - resids (defaults to '1')


    .. note::
        The precision is hard coded to three decimal places.


    .. note::
        When dimensions are missing (i.e. set to `None`), a zero width
        unit cell box will be written (i.e. [0.0, 0.0, 0.0]).


    .. versionchanged:: 0.11.0
       Frames now 0-based instead of 1-based
    .. versionchanged:: 0.13.0
       Now strictly writes positions with 3dp precision.
       and velocities with 4dp.
       Removed the `convert_dimensions_to_unitcell` method,
       use `Timestep.triclinic_dimensions` instead.
       Now now writes velocities where possible.
    .. versionchanged:: 0.18.0
       Added `reindex` keyword argument to allow original atom
       ids to be kept.
    .. versionchanged:: 2.0.0
       Raises a warning when writing timestep with missing dimension
       information (i.e. set to ``None``).
    """

    format = "GRO"
    units = {"time": None, "length": "nm", "velocity": "nm/ps"}
    gro_coor_limits = {"min": -999.9995, "max": 9999.9995}

    #: format strings for the GRO file (all include newline); precision
    #: of 3 decimal places is hard-coded here.
    fmt = {
        "n_atoms": "{0:5d}\n",  # number of atoms
        # coordinates output format, see http://chembytes.wikidot.com/g-grofile
        "xyz": "{resid:>5d}{resname:<5.5s}{name:>5.5s}{index:>5d}{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}\n",
        # unitcell
        "box_orthorhombic": "{box[0]:10.5f} {box[1]:9.5f} {box[2]:9.5f}\n",
        "box_triclinic": "{box[0]:10.5f} {box[4]:9.5f} {box[8]:9.5f} {box[1]:9.5f} {box[2]:9.5f} {box[3]:9.5f} {box[5]:9.5f} {box[6]:9.5f} {box[7]:9.5f}\n",
    }
    fmt["xyz_v"] = (
        fmt["xyz"][:-1] + "{vel[0]:8.4f}{vel[1]:8.4f}{vel[2]:8.4f}\n"
    )

    def __init__(self, filename, convert_units=True, n_atoms=None, **kwargs):
        """Set up a GROWriter with a precision of 3 decimal places.

        Parameters
        -----------
        filename : str
            output filename
        n_atoms : int (optional)
            number of atoms
        convert_units : bool (optional)
            units are converted to the MDAnalysis base format; [``True``]
        reindex : bool (optional)
            By default, all the atoms were reindexed to have a atom id starting
            from 1. [``True``] However, this behaviour can be turned off by
            specifying `reindex` ``=False``.

        Note
        ----
        To use the reindex keyword, user can follow the two examples given
        below.::

           u = mda.Universe()

        Usage 1::

           u.atoms.write('out.gro', reindex=False)

        Usage 2::

           with mda.Writer('out.gro', reindex=False) as w:
               w.write(u.atoms)

        """
        self.filename = util.filename(filename, ext="gro", keep=True)
        self.n_atoms = n_atoms
        self.reindex = kwargs.pop("reindex", True)

        self.convert_units = (
            convert_units  # convert length and time to base units
        )

    def write(self, obj):
        """Write selection at current trajectory frame to file.

        Parameters
        -----------
        obj : AtomGroup or Universe

        Note
        ----
        The GRO format only allows 5 digits for *resid* and *atom
        number*. If these numbers become larger than 99,999 then this
        routine will chop off the leading digits.


        .. versionchanged:: 0.7.6
           *resName* and *atomName* are truncated to a maximum of 5 characters
        .. versionchanged:: 0.16.0
           `frame` kwarg has been removed
        .. versionchanged:: 2.0.0
           Deprecated support for calling with Timestep has nwo been removed.
           Use AtomGroup or Universe as an input instead.
        """
        # write() method that complies with the Trajectory API

        try:

            # make sure to use atoms (Issue 46)
            ag = obj.atoms
            # can write from selection == Universe (Issue 49)

        except AttributeError:
            errmsg = "Input obj is neither an AtomGroup or Universe"
            raise TypeError(errmsg) from None

        try:
            velocities = ag.velocities
        except NoDataError:
            has_velocities = False
        else:
            has_velocities = True

        # Check for topology information
        missing_topology = []
        try:
            names = ag.names
        except (AttributeError, NoDataError):
            names = itertools.cycle(("X",))
            missing_topology.append("names")
        try:
            resnames = ag.resnames
        except (AttributeError, NoDataError):
            resnames = itertools.cycle(("UNK",))
            missing_topology.append("resnames")
        try:
            resids = ag.resids
        except (AttributeError, NoDataError):
            resids = itertools.cycle((1,))
            missing_topology.append("resids")

        if not self.reindex:
            try:
                atom_indices = ag.ids
            except (AttributeError, NoDataError):
                atom_indices = range(1, ag.n_atoms + 1)
                missing_topology.append("ids")
        else:
            atom_indices = range(1, ag.n_atoms + 1)
        if missing_topology:
            warnings.warn(
                "Supplied AtomGroup was missing the following attributes: "
                "{miss}. These will be written with default values. "
                "Alternatively these can be supplied as keyword arguments."
                "".format(miss=", ".join(missing_topology))
            )

        positions = ag.positions

        if self.convert_units:
            # Convert back to nm from Angstroms,
            # Not inplace because AtomGroup is not a copy
            positions = self.convert_pos_to_native(positions, inplace=False)
            if has_velocities:
                velocities = self.convert_velocities_to_native(
                    velocities, inplace=False
                )
        # check if any coordinates are illegal
        # (checks the coordinates in native nm!)
        if not self.has_valid_coordinates(self.gro_coor_limits, positions):
            raise ValueError(
                "GRO files must have coordinate values between "
                "{0:.3f} and {1:.3f} nm: No file was written."
                "".format(
                    self.gro_coor_limits["min"], self.gro_coor_limits["max"]
                )
            )

        with util.openany(self.filename, "wt") as output_gro:
            # Header
            output_gro.write("Written by MDAnalysis\n")
            output_gro.write(self.fmt["n_atoms"].format(ag.n_atoms))

            # Atom descriptions and coords
            # Dont use enumerate here,
            # all attributes could be infinite cycles!
            for atom_index, resid, resname, name in zip(
                range(ag.n_atoms), resids, resnames, names
            ):
                truncated_atom_index = util.ltruncate_int(
                    atom_indices[atom_index], 5
                )
                truncated_resid = util.ltruncate_int(resid, 5)
                if has_velocities:
                    output_gro.write(
                        self.fmt["xyz_v"].format(
                            resid=truncated_resid,
                            resname=resname,
                            index=truncated_atom_index,
                            name=name,
                            pos=positions[atom_index],
                            vel=velocities[atom_index],
                        )
                    )
                else:
                    output_gro.write(
                        self.fmt["xyz"].format(
                            resid=truncated_resid,
                            resname=resname,
                            index=truncated_atom_index,
                            name=name,
                            pos=positions[atom_index],
                        )
                    )

            # Footer: box dimensions
            if ag.dimensions is None or np.allclose(
                ag.dimensions[3:], [90.0, 90.0, 90.0]
            ):
                if ag.dimensions is None:
                    wmsg = (
                        "missing dimension - setting unit cell to zeroed "
                        "box [0., 0., 0.]"
                    )
                    warnings.warn(wmsg)
                    box = np.zeros(3)
                else:
                    box = self.convert_pos_to_native(
                        ag.dimensions[:3], inplace=False
                    )
                # orthorhombic cell, only lengths along axes needed in gro
                output_gro.write(self.fmt["box_orthorhombic"].format(box=box))
            else:
                try:  # for AtomGroup/Universe
                    tri_dims = obj.universe.coord.triclinic_dimensions
                except AttributeError:  # for Timestep
                    tri_dims = obj.triclinic_dimensions
                # full output
                box = self.convert_pos_to_native(
                    tri_dims.flatten(), inplace=False
                )
                output_gro.write(self.fmt["box_triclinic"].format(box=box))
