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


"""AMBER coordinate/restart files in MDAnalysis --- :mod:`MDAnalysis.coordinates.INPCRD`
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
.. _AMBER INPCRD FORMAT: https://ambermd.org/FileFormats.php#restart
.. _AMBER netcdf: http://ambermd.org/netcdf/nctraj.xhtml
.. _NetCDF: http://www.unidata.ucar.edu/software/netcdf

"""

from . import base
import warnings
import logging
from scipy.io import netcdf_file

from .timestep import Timestep
from ..lib.util import store_init_arguments
from .TRJ import NCDFPicklable, NCDFMixin

logger = logging.getLogger("MDAnalysis.coordinates.AMBER")


class INPReader(base.SingleFrameReaderBase):
    """Reader for Amber restart files.

    .. rubric:: Limitations

    * Box information is not read (or checked for).
    * Velocities are currently *not supported*.

    .. versionchanged: 2.10.0
       Now automatically detects files with .rst7 extension.

    """

    format = ["INPCRD", "RESTRT", "RST7"]
    units = {"length": "Angstrom"}

    def _read_first_frame(self):
        # Read header
        with open(self.filename, "r") as inf:
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
                for i, dest in enumerate(
                    [
                        (2 * p, 0),
                        (2 * p, 1),
                        (2 * p, 2),
                        (2 * p + 1, 0),
                        (2 * p + 1, 1),
                        (2 * p + 1, 2),
                    ]
                ):
                    self.ts._pos[dest] = float(line[i * 12 : (i + 1) * 12])
            # Read last coordinate if necessary
            if self.n_atoms % 2:
                line = inf.readline()
                for i in range(3):
                    self.ts._pos[-1, i] = float(line[i * 12 : (i + 1) * 12])

    @staticmethod
    def parse_n_atoms(filename, **kwargs):
        with open(filename, "r") as f:
            f.readline()
            n_atoms = int(f.readline().split()[0])
        return n_atoms


class NCRSTReader(base.SingleFrameReaderBase, NCDFMixin):
    """Reader for `AMBER NETCDF format`_ (version 1.0 rev C) restart files.

    This reader is a :class:`SingleFrameReaderBase` adaptation of the
    :class:`~MDAnalysis.coordinates.TRJ.NCDFReader` AMBER NETCDF trajectory reader.

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

    Support for the *mmap* keyword is available as detailed
    in :class:`~MDAnalysis.coordinates.TRJ.NCDFReader` and
    :mod:`scipy.io.netcdf.netcdf_file`. The use of
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
    :class:`~MDAnalysis.coordinates.TRJ.NCDFReader`
    :class:`~MDAnalysis.coordinates.TRJ.NCDFWriter`


    .. versionadded: 2.10.0
    """

    format = ["NCRST", "NCRESTRT", "NCRST7"]
    version = "1.0"
    units = {
        "time": "ps",
        "length": "Angstrom",
        "velocity": "Angstrom/ps",
        "force": "kcal/(mol*Angstrom)",
    }

    _Timestep = Timestep

    @store_init_arguments
    def __init__(
        self, filename, n_atoms=None, convert_units=None, mmap=None, **kwargs
    ):
        # Assign input mmap value
        self._mmap = mmap

        super(NCRSTReader, self).__init__(
            filename, convert_units, n_atoms, **kwargs
        )

    @staticmethod
    def parse_n_atoms(filename, **kwargs):
        with scipy.io.netcdf_file(filename, mmap=None) as f:
            n_atoms = f.dimensions["atom"]
        return n_atoms

    @staticmethod
    def _verify_units(eval_units, expected_units):
        if eval_units.decode("utf-8") != expected_units:
            errmsg = (
                "NCRSTReader currently assumes that the trajectory "
                "was written in units of {0} instead of {1}".format(
                    eval_units.decode("utf-8"), expected_units
                )
            )
            raise NotImplementedError(errmsg)

    def _read_first_frame(self):
        """Function to read NetCDF restart file and fill timestep"""
        # Open netcdf file via context manager
        # ensure maskandscale is off so we don't end up double scaling

        with NCDFPicklable(
            self.filename, mode="r", mmap=self._mmap, maskandscale=False
        ) as self.trjfile:

            self._check_conventions(n_atoms=self.n_atoms)

            self.n_frames = 1

            self._read_values(
                frame=()
            )  # AMBERRESTART convention files have dimensionless datasets

            self.ts.frame = 0  # 0-indexed
