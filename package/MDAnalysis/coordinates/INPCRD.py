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


"""INPCRD structure files in MDAnalysis --- :mod:`MDAnalysis.coordinates.INPCRD`
================================================================================

Read coordinates in Amber_ coordinate/restart file (suffix "inpcrd").

.. _Amber: http://ambermd.org/formats.html#restart


Classes
-------

.. autoclass:: INPReader
   :members:

"""
from __future__ import absolute_import, division, print_function, unicode_literals

from six.moves import range

from . import base
from ..core import flags
import scipy.io.netcdf
import warnings
import logging
import errno
import numpy as np

logger = logging.getLogger("MDAnalysis.coordinates.AMBER")

try:
    import netCDF4
except ImportError:
    netCDF4 = None
    logger.warning("netCDF4 is not available. Writing AMBER NCRST will be slow.")


class INPReader(base.SingleFrameReaderBase):
    """Reader for Amber restart files."""

    format = ['INPCRD', 'RESTRT']
    units = {'length': 'Angstrom'}

    def _read_first_frame(self):
        # Read header
        with open(self.filename, 'r') as inf:
            self.title = inf.readline().strip()
            line = inf.readline().split()
            self.n_atoms = int(line[0])

            self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
            try:
                time = float(line[1])
            except IndexError:
                pass
            else:
                self.ts.time = time
            self.ts.frame = 0

            for p in range(self.n_atoms // 2):
                line = inf.readline()
                # each float is f12.7, 6 floats a line
                for i, dest in enumerate([(2*p, 0), (2*p, 1), (2*p, 2),
                                          (2*p + 1, 0), (2*p + 1, 1), (2*p + 1, 2)]):
                    self.ts._pos[dest] = float(line[i*12:(i+1)*12])
            # Read last coordinate if necessary
            if self.n_atoms % 2:
                line = inf.readline()
                for i in range(3):
                    self.ts._pos[-1, i] = float(line[i*12:(i+1)*12])

    @staticmethod
    def parse_n_atoms(filename, **kwargs):
        with open(filename, 'r') as f:
            f.readline()
            n_atoms = int(f.readline().split()[0])
        return n_atoms


class NCRSTReader(base.SingleFrameReaderBase):
    """Reader for Amber NetCDF restart files.
    To do: detailed information here.
    """

    format = ['NCRST', 'NCRESTRT']
    version = "1.0"
    units = {'time': 'ps',
             'length': 'Angstrom',
             'velocity': 'Angstrom/ps',
             'forces': 'kcal/(mol*Angstrom)'}


    def __init__(self, filename, n_atoms=None, convert_units=None, mmap=None, **kwargs):
        # Add optional mmap here
        self._mmap = mmap
        super(NCRSTReader, self).__init__(filename, convert_units, n_atoms, **kwargs)


    @staticmethod
    def parse_natoms(filename, **kwargs):
        with scipy.io.netcdf.netcdf_file(filename, mmap=None) as f:
            n_atoms = f.dimensions['atom']
        return n_atoms


    def _read_first_frame(self):
        # Open netcdf file as a context manager
        with scipy.io.netcdf.netcdf_file(self.filename, mmap=self._mmap) as rstfile:
            # Check file contents
            if not ('AMBERRESTART' in
                    rstfile.Conventions.decode('utf-8').split(',') or
                    'AMBERRESTART' in
                    rstfile.Conventions.decode('utf-8').split()):
                errmsg = ("NCDF restart file {0} does not conform to AMBER "
                          "specifications, http://ambermd.org/netcdf/nctraj.xhtml "
                          "('AMBERRESTART' must be one of the tokens in attribute "
                          "Conventions)".format(self.filename))
                logger.fatal(errmsg)
                raise TypeError(errmsg)

            if not rstfile.ConventionVersion.decode('utf-8') == self.version:
                wmsg = ("NCDF trajectory format is {0!s} but the reader "
                        "implements format {1!s}".format(
                        rstfile.ConventionVersion, self.version))
                warnings.warn(wmsg)
                logger.warning(wmsg)

            if rstfile.variables['time'].units.decode('utf-8') != "picosecond":
                raise NotImplementedError(
                    "NCRSTReader currently assumes that the trajectory was written "
                    "with a time unit of picoseconds and not {0}.".format(
                        rstfile.variables['time'].units))
            if rstfile.variables['coordinates'].units.decode('utf-8') != "angstrom":
                raise NotImplementedError(
                    "NCRSTReader currently assumes that the trajectory was written "
                    "with a length unit of Angstroem and not {0}.".format(
                        rstfile.variables['coordinates'].units))
            if hasattr(rstfile.variables['coordinates'], 'scale_factor'):
                raise NotImplementedError("scale_factors are not implemented")

            # Get remarks
            try:
                self.remarks = rstfile.title
            except AttributeError:
                self.remarks = ""
            # Get n_atoms - Base class assigns self.n_atom to the passed value
            # of n_atoms which make for a slightly confusing check
            self.n_atoms = rstfile.dimensions['atom']
            if self.n_atom is not None and self.n_atom != self.n_atoms:
                raise ValueError("Supplied n_atoms ({0}) != natom from ncdf ({1}). "
                                 "Note: n_atoms can be None and then the ncdf value "
                                 "is used!".format(n_atoms, self.n_atoms))

            # Check for the presence of optional variables
            self.has_velocities = 'velocities' in rstfile.variables
            self.has_forces = 'forces' in rstfile.variables
            self.periodic = 'cell_lengths' in rstfile.variables
            # Set timestep
            self.ts = self._Timestep(self.n_atom,
                                     velocities=self.has_velocities,
                                     forces=self.has_forces,
                                     reader=self,
                                     **self._ts_kwargs)

            # Fetch file data
            self.ts._pos[:] = rstfile.variables['coordinates'][:]
            self.ts.time = rstfile.variables['time'].getValue()
            if self.has_velocities:
                self.ts._velocities[:] = rstfile.variables['velocities'][:]
            if self.has_forces:
                self.ts._forces[:] = rstfile.variables['forces'][:]
            if self.periodic:
                self.ts._unitcell[:3] = rstfile.variables['cell_lengths'][:]
                self.ts._unitcell[3:] = rstfile.variables['cell_angles'][:]

            # Convert units
            if self.convert_units:
                self.convert_pos_from_native(self.ts._pos)
                self.convert_time_from_native(self.ts.time)
                if self.has_velocities:
                    self.convert_velocities_from_native(self.ts._velocities,
                                                        inplace=True)
            if self.has_forces:
                self.convert_forces_from_native(self.ts._forces, inplace=True)
            if self.periodic:
                self.convert_pos_from_native(self.ts._unitcell[:3])
