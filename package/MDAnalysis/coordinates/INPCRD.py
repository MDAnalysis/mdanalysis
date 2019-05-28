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


"""AMBER restart files in MDAnalysis --- :mod:`MDAnalysis.coordinates.INPCRD`
================================================================================

AMBER_ can write :ref:`ASCII restart<ascii-restart>` ("inpcrd") and
:ref:`binary restart<netcdf-restart>` ("ncrst") coordinate files. MDAnalysis
supports reading of both file formats.

.. rubric:: Units

AMBER restart files are assumed to be in the following units:

* length in Angstrom (Å)
* time in ps
* velocity (NCRST only) in Å / ps
* force (NCRST only) in kcal / (mol * Å)


.. _ascii-restart:

ASCII INPCRD restart files
--------------------------

ASCII AMBER_ INPCRD coordinate files (as defined in `AMBER INPCRD FORMAT`_)
are handled by the :class:`INPReader`.

AMBER ASICC restart files are recognised by the suffix '.inpcrd', '.restrt', or
'.rst7'

.. autoclass:: INPReader
   :members:


.. _netcdf-restart:

Binary NetCDF restart files
---------------------------

The `AMBER netcdf`_ restart format makes use of NetCDF_ (Network Common Data
Form) format. Such binary restart files are recognised in MDAnalysis by the
suffix '.ncrst', '.ncrestrt' or '.ncrst7' and read by the :class:`NCRSTReader`.

Binary restart files can also contain velocity and force information, and can
record the simulation time step. Whilst the `AMBER netcdf`_ format details
default unit values of ångström and picoseconds, these can in theory occupy any
unit type. However, at the moment MDAnalysis only supports the default types
and will raise a :exc:`NotImplementedError` if anything else is detected.

.. autoclass:: NCRSTReader
   :members:



.. Links

.. _AMBER: http://ambermd.org
.. _AMBER INPCRD FORMAT: http://ambermd.org/formats.html#restart
.. _AMBER netcdf: http://ambermd.org/netcdf/nctraj.xhtml
.. _NetCDF: http://www.unidata.ucar.edu/software/netcdf

"""
from __future__ import absolute_import, division, print_function, unicode_literals

from six.moves import range

from . import base
import scipy.io.netcdf
import warnings
import logging

logger = logging.getLogger("MDAnalysis.coordinates.AMBER")


class INPReader(base.SingleFrameReaderBase):
    """Reader for Amber restart files.

    .. rubric:: Limitations

    * Box information is not read (or checked for).
    * Velocities are currently *not supported*.

    """

    format = ['INPCRD', 'RESTRT', 'RST7']
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
    """Reader for `AMBER NETCDF format`_ (version 1.0 rev C) restart files.

    This reader is a :class:`SingleFrameReaderBase` adaptation of the
    :class:`NCDFReader` AMBER NETCDF trajectory reader.

    AMBER binary restart files are automatically recognised by the file
    extensions ".ncrst", ".ncrestrt", and ".ncrst7".

    The number of atoms (`n_atoms`) does not have to be provided as it can
    be read from the input NETCDF file.

    Current simulation time is autodetected and read into the
    :attr:`Timestep.time` attribute. (If this is not available in the restart
    file, then :attr:`Timestep.time` will return 0.0 ps).

    Velocities are autodetected and read into the :attr:`Timestep._velocities`
    attribute.

    Forces are autodetected and read into the :attr:`Timestep._forces`
    attribute.

    Periodic unit cell information is detected and used to populate the
    :attr:`Timestep.dimensions` attribute. (If not unit cell is available in
    the restart file, then :attr:`Timestep.dimensions` will return
    ``[0,0,0,0,0,0]``).

    The NCRST reader uses :mod:`scipy.io.netcdf` and therefore :mod:`scipy`
    must be installed. Support for the *mmap* keyword is available as detailed
    in :class:`NCDFReader` and :mod:`scipy.io.netcdf.netcdf_file`.

    The NCRST reader also uses a custom Timestep object with C-style memory
    mapping in order to match the NCDFReader.

    .. rubric:: Limitations

    * Only NCRST files with time in ps and lengths in Angstroem are processed.
    * scale_factors are not supported (and not checked).
    * Restart files without coordinate information are not supported.
    * Replica exchange variables are not supported.

    .. _AMBER NETCDF format: http://ambermd.org/netcdf/nctraj.xhtml

    See Also
    --------
    :class:`NCDFReader`
    :class:`NCDFWriter`

    .. versionadded: 0.18.1

    """

    format = ['NCRST', 'NCRESTRT', 'NCRST7']
    version = "1.0"
    units = {'time': 'ps',
             'length': 'Angstrom',
             'velocity': 'Angstrom/ps',
             'forces': 'kcal/(mol*Angstrom)'}

    class Timestep(base.Timestep):
        """ Modified Timestep class for AMBER

        Uses C order memory mapping to match the style used AMBER TRAJ

        The Timestep can be initialized with `arg` being an integer
        (the number of atoms) and an optional keyword arguments to allocate
        space for velocities, and forces.

        .. versionchanged:: 0.10.0
           Added ability to contain Forces
        """
        order = 'C'

    _Timestep = Timestep

    def __init__(self, filename, n_atoms=None, convert_units=None, mmap=None,
                 **kwargs):
        # Assign input mmap value
        self._mmap = mmap
        super(NCRSTReader, self).__init__(filename, convert_units, n_atoms,
                                          **kwargs)

    @staticmethod
    def parse_n_atoms(filename, **kwargs):
        with scipy.io.netcdf.netcdf_file(filename, mmap=None) as f:
            n_atoms = f.dimensions['atom']
        return n_atoms

    def _read_first_frame(self):
        """Function to read NetCDF restart file and fill timestep

        Called by: :class:`SingleFrameReaderBase`.__init__
        Overrides :class:`SingleFrameReaderBase` placeholder function
        """
        # Open netcdf file via context manager
        with scipy.io.netcdf.netcdf_file(self.filename, mode='r', mmap=self._mmap) as rstfile:
            # Global attribute checks
            # Conventions should contain the AMBERRESTART string
            if not ('AMBERRESTART' in
                    rstfile.Conventions.decode('utf-8').split(',') or
                    'AMBERRESTART' in
                    rstfile.Conventions.decode('utf-8').split()):
                errmsg = ("NetCDF restart file {0} does not conform to AMBER "
                          "specifications, as detailed in "
                          "http://ambermd.org/netcdf/nctraj.xhtml "
                          "('AMBERRESTART' must be one of the tokens in "
                          "attribute Conventions)".format(self.filename))
                logger.fatal(errmsg)
                raise TypeError(errmsg)

            # The AMBER NetCDF standard enforces 64 bit offsets
            if not rstfile.version_byte == 2:
                errmsg = ("NetCDF restart file {0} does not conform to AMBER "
                          "specifications, as details in "
                          "http://ambermd.org/netcdf/nctraj.xhtml "
                          "(NetCDF file does not use 64 bit offsets "
                          "[version_byte = 2]) ".format(self.filename))
                logger.fatal(errmsg)
                raise TypeError(errmsg)

            # ConventionVersion should exist and be equal to 1.0
            if not rstfile.ConventionVersion.decode('utf-8') == self.version:
                wmsg = ("NCDF restart format is {0!s} but the reader "
                        "implements format {1!s}".format(
                         rstfile.ConventionVersion, self.version))
                warnings.warn(wmsg)
                logger.warning(wmsg)

            # The use of scale_factor is currently unsupported
            if hasattr(rstfile.variables['coordinates'], 'scale_factor'):
                raise NotImplementedError("scale_factors are not implemented")

            # Note: SingleFrameReaderBase class sets parsed n_atoms value to
            # self.n_atom which makes for a confusing check
            try:
                self.n_atoms = rstfile.dimensions['atom']
                if self.n_atom is not None and self.n_atom != self.n_atoms:
                    raise ValueError("Supplied n_atoms ({0}) != n_atoms from "
                                     "NetCDF restart file ({1}). "
                                     "Note: n_atoms can be None and then the "
                                     "restart file value is used.".format(
                                      self.n_atom, self.n_atoms))
            except KeyError:
                errmsg = ("NetCDF restart file {0} does not contain "
                          "atom information ".format(self.filename))
                logger.fatal(errmsg)
                raise KeyError(errmsg)

            # Optional variables
            self.has_velocities = 'velocities' in rstfile.variables
            self.has_forces = 'forces' in rstfile.variables
            self.periodic = 'cell_lengths' in rstfile.variables

            # Set timestep
            self.ts = self._Timestep(self.n_atoms,
                                     velocities=self.has_velocities,
                                     forces=self.has_forces,
                                     reader=self,
                                     **self._ts_kwargs)

            # NetCDF file title is optional
            try:
                self.remarks = rstfile.title
            except AttributeError:
                self.remarks = ""

            # Default to time units of picoseconds
            # Note: unlike trajectories the AMBER NetCDF convention allows
            # restart files to omit time when unecessary (i.e. minimizations)
            try:
                units = rstfile.variables['time'].units.decode('utf-8')
                if units != "picosecond":
                    raise NotImplementedError(
                        "NCRSTReader currently assumes that the restart file "
                        "uses a time unit of picoseconds and not {0}.".format(
                         rstfile.variables['time'].units))
                self.ts.time = rstfile.variables['time'].getValue()
            except KeyError:
                # As of AMBER16 the NetCDF restart files created by
                # minimizations ignore convention and default to 0.0 ps
                # Warn and do the same thing here (IA: need to check this)
                wmsg = ("Restart file {0} does not contain time information "
                        "time will default to 0.0 ps").format(self.filename)
                warnings.warn(wmsg)
                logger.warning(wmsg)
                self.ts.time = 0.0

            # Single frame so we assign it to 0
            self.ts.frame = 0

            # Default to length units of Angstroem
            try:
                units = rstfile.variables['coordinates'].units.decode('utf-8')
                if units != "angstrom":
                    raise NotImplementedError(
                        "NCRSTReader currently assumes that the restart file "
                        "uses a length unit of Angstroem and not {0}.".format(
                         rstfile.variables['coordinates'].units))
                self.ts._pos[:] = rstfile.variables['coordinates'][:]
            except KeyError:
                # Technically coordinate information is not required, however
                # the lack of it in a restart file is highly unlikely
                errmsg = ("NetCDF restart file {0} is missing coordinate "
                          "information ".format(self.filename))
                logger.fatal(errmsg)
                raise KeyError(errmsg)

            if self.has_velocities:
                self.ts._velocities[:] = rstfile.variables['velocities'][:]

            # The presence of forces in a restart file is very rare, but
            # according to AMBER convention, it can exist.
            if self.has_forces:
                self.ts._forces[:] = rstfile.variables['forces'][:]

            # If false u.dimensions is set to [0., 0., 0., 0., 0., 0]
            if self.periodic:
                self.ts._unitcell[:3] = rstfile.variables['cell_lengths'][:]
                self.ts._unitcell[3:] = rstfile.variables['cell_angles'][:]

            # Convert units inplace
            if self.convert_units:
                self.convert_pos_from_native(self.ts._pos)
                self.convert_time_from_native(self.ts.time)
                if self.has_velocities:
                    self.convert_velocities_from_native(self.ts._velocities,
                                                        inplace=True)
                if self.has_forces:
                    self.convert_forces_from_native(self.ts._forces,
                                                    inplace=True)
                if self.periodic:
                    self.convert_pos_from_native(self.ts._unitcell[:3])
