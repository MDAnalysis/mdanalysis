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
default unit values of ångström and picoseconds, these can in theory occupy
any unit type. However, at the moment MDAnalysis only supports the default
types and will raise a :exc:`NotImplementedError` if anything else is detected.

.. autoclass:: NCRSTReader
   :members:


.. Links

.. _AMBER: http://ambermd.org
.. _AMBER INPCRD FORMAT: http://ambermd.org/formats.html#restart
.. _AMBER netcdf: http://ambermd.org/netcdf/nctraj.xhtml
.. _NetCDF: http://www.unidata.ucar.edu/software/netcdf

"""

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

    .. versionchanged: 0.20.0
       Now automatically detects files with .rst7 extension.

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

    Current simulation time is autodetected and if available is read into the
    :attr:`Timestep.time` attribute.

    Velocities are autodetected and read into the :attr:`Timestep._velocities`
    attribute.

    Forces are autodetected and read into the :attr:`Timestep._forces`
    attribute.

    Periodic unit cell information is detected and used to populate the
    :attr:`Timestep.dimensions` attribute. (If no unit cell is available in
    the restart file, then :attr:`Timestep.dimensions` will return
    ``[0,0,0,0,0,0]``).

    The NCRST reader uses :mod:`scipy.io.netcdf` and therefore :mod:`scipy`
    must be installed.

    Support for the *mmap* keyword is available as detailed
    in :class:`NCDFReader` and :mod:`scipy.io.netcdf.netcdf_file`. The use of
    ``mmap=True`` leads to around a 2x read speed improvement in a ~ 1 million
    atom system (AMBER STMV benchmark). As per the :class:`NCDFReader`, the
    default behaviour is ``mmap=None``, which means that the default behaviour
    of :class:`scipy.io.netcdf.netcdf_file` prevails.

    The NCRST reader also uses a custom Timestep object with C-style memory
    mapping in order to match the NCDFReader.

    .. rubric:: Limitations

    * Only NCRST files with time in ps, lengths in Angstroem and angles in
      degree are processed.
    * Restart files without coordinate information are not supported.
    * Replica exchange variables are not supported.

    .. _AMBER NETCDF format: http://ambermd.org/netcdf/nctraj.xhtml

    See Also
    --------
    :class:`NCDFReader`
    :class:`NCDFWriter`

    .. versionadded: 0.20.0

    """

    format = ['NCRST', 'NCRESTRT', 'NCRST7']
    version = "1.0"
    units = {'time': 'ps',
             'length': 'Angstrom',
             'velocity': 'Angstrom/ps',
             'force': 'kcal/(mol*Angstrom)'}

    class Timestep(base.Timestep):
        """ Modified Timestep class for AMBER

        Uses C order memory mapping to match the style used by AMBER TRAJ

        The Timestep can be initialized with `arg` being an integer
        (the number of atoms) and an optional keyword arguments to allocate
        space for velocities, and forces.
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

    @staticmethod
    def _verify_units(eval_units, expected_units):
        if eval_units.decode('utf-8') != expected_units:
            errmsg = ("NCRSTReader currently assumes that the trajectory "
                      "was written in units of {0} instead of {1}".format(
                       eval_units.decode('utf-8'), expected_units))
            raise NotImplementedError(errmsg)

    def _read_first_frame(self):
        """Function to read NetCDF restart file and fill timestep

        Called by: :class:`SingleFrameReaderBase`.__init__
        Overrides :class:`SingleFrameReaderBase` placeholder function
        """
        # Open netcdf file via context manager
        with scipy.io.netcdf.netcdf_file(self.filename, mode='r',
                                         mmap=self._mmap) as rstfile:
            # Conventions should exist and contain the AMBERRESTART string
            try:
                conventions = rstfile.Conventions.decode('utf-8')
                if not ('AMBERRESTART' in conventions.split(',') or
                        'AMBERRESTART' in conventions.split()):
                    errmsg = ("NetCDF restart file {0} does not conform to "
                              "AMBER specifications, as detailed in "
                              "http://ambermd.org/netcdf/nctraj.xhtml "
                              "('AMBERRESTART' must be one of the tokens in "
                              "attribute Conventions)".format(self.filename))
                    logger.fatal(errmsg)
                    raise TypeError(errmsg)
            except AttributeError:
                errmsg = ("NetCDF restart file {0} is missing  a "
                          "Conventions value".format(self.filename))
                raise ValueError(errmsg)

            # ConventionVersion should exist and be equal to 1.0
            try:
                ConventionVersion = rstfile.ConventionVersion.decode('utf-8')
                if not (ConventionVersion == self.version):
                    wmsg = ("NCRST format is {0!s} but the reader implements "
                            "format {1!s}".format(ConventionVersion,
                             self.version))
                    warnings.warn(wmsg)
                    logger.warning(wmsg)
            except AttributeError:
                errmsg = ("NCDF restart file {0} is missing a "
                          "ConventionVersion value".format(self.filename))
                raise ValueError(errmsg)

            # The AMBER NetCDF standard enforces 64 bit offsets
            if not rstfile.version_byte == 2:
                errmsg = ("NetCDF restart file {0} does not conform to AMBER "
                          "specifications, as detailed in "
                          "http://ambermd.org/netcdf/nctraj.xhtml "
                          "(NetCDF file does not use 64 bit offsets "
                          "[version_byte = 2]) ".format(self.filename))
                logger.fatal(errmsg)
                raise TypeError(errmsg)

            # The specs require that dimensions->spatial == 3 exists
            try:
                if not rstfile.dimensions['spatial'] == 3:
                    errmsg = "Incorrect spatial value for NCRST file"
                    raise TypeError(errmsg)
            except KeyError:
                errmsg = ("NCDF restart file {0} does not contain spatial "
                          "dimension".format(self.filename))
                raise ValueError(errmsg)

            # The specs define program and programVersion as required. Here we
            # just warn the users instead of raising an Error.
            if not (hasattr(rstfile, 'program') and
                    hasattr(rstfile, 'programVersion')):
                wmsg = ("This NCRST file {0} may not fully adhere to AMBER "
                        "standards as either the `program` or "
                        "`programVersion` attributes are missing".format(
                        self.filename))
                warnings.warn(wmsg)
                logger.warning(wmsg)

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
                raise ValueError(errmsg)

            # NetCDF file title is optional
            try:
                self.remarks = rstfile.title
            except AttributeError:
                self.remarks = ""

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

            # Check for scale_factor attributes and store them
            scale_factors = {'time': 1.0,
                             'cell_lengths': 1.0,
                             'cell_angles': 1.0,
                             'coordinates': 1.0,
                             'velocities': 1.0,
                             'forces': 1.0}

            for variable in rstfile.variables:
                if hasattr(rstfile.variables[variable], 'scale_factor'):
                    if variable in scale_factors:
                        factor = rstfile.variables[variable].scale_factor
                        scale_factors[variable] = factor
                    else:
                        errmsg = ("scale_factors for variable {0} are not "
                                  "implemented".format(variable))
                        raise NotImplementedError(errmsg)

            # Note: unlike trajectories the AMBER NetCDF convention allows
            # restart files to omit time when unecessary (i.e. minimizations)
            try:
                self._verify_units(rstfile.variables['time'].units,
                                   'picosecond')
                self.ts.time = (rstfile.variables['time'].getValue() *
                                scale_factors['time'])
            except KeyError:
                # Warn the user and move on
                wmsg = ("NCRestart file {0} does not contain time "
                        "information. This should be expected if the file was "
                        "not created from an MD trajectory (e.g. a "
                        "minimization)".format(self.filename))
                warnings.warn(wmsg)
                logger.warning(wmsg)

            # Single frame so we assign it to 0
            self.ts.frame = 0

            # Default to length units of Angstrom
            try:
                self._verify_units(rstfile.variables['coordinates'].units,
                                   'angstrom')
                self.ts._pos[:] = (rstfile.variables['coordinates'][:] *
                                   scale_factors['coordinates'])
            except KeyError:
                # Technically coordinate information is not required, however
                # the lack of it in a restart file is highly unlikely
                errmsg = ("NetCDF restart file {0} is missing coordinate "
                          "information ".format(self.filename))
                logger.fatal(errmsg)
                raise ValueError(errmsg)

            if self.has_velocities:
                self._verify_units(rstfile.variables['velocities'].units,
                                   'angstrom/picosecond')
                self.ts._velocities[:] = (rstfile.variables['velocities'][:] *
                                          scale_factors['velocities'])

            # The presence of forces in an ncrst is rare but possible
            if self.has_forces:
                self._verify_units(rstfile.variables['forces'].units,
                                   'kilocalorie/mole/angstrom')
                self.ts._forces[:] = (rstfile.variables['forces'][:] *
                                      scale_factors['forces'])

            # If false u.dimensions is set to [0., 0., 0., 0., 0., 0]
            # Unlike the NCDFReader, `degrees` is not accepted
            if self.periodic:
                self._verify_units(rstfile.variables['cell_lengths'].units,
                                   'angstrom')
                self._verify_units(rstfile.variables['cell_angles'].units,
                                   'degree')
                self.ts._unitcell[:3] = (rstfile.variables['cell_lengths'][:] *
                                         scale_factors['cell_lengths'])
                self.ts._unitcell[3:] = (rstfile.variables['cell_angles'][:] *
                                         scale_factors['cell_angles'])

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
