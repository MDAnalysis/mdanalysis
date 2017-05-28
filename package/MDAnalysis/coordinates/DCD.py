# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.mdanalysis.org
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


"""DCD trajectory I/O  --- :mod:`MDAnalysis.coordinates.DCD`
============================================================

Classes to read and write DCD binary trajectories, the format used by
CHARMM, NAMD, and also LAMMPS. Trajectories can be read regardless of
system-endianness as this is auto-detected.

Generally, DCD trajectories produced by any code can be read (with the
:class:`DCDReader`) although there can be issues with the unitcell
(simulation box) representation (see
:attr:`Timestep.dimensions`). DCDs can also be written but the
:class:`DCDWriter` follows recent NAMD/VMD convention for the unitcell
but still writes AKMA time. Reading and writing these trajectories
within MDAnalysis will work seamlessly but if you process those
trajectories with other tools you might need to watch out that time
and unitcell dimensions are correctly interpreted.

Note
----
The DCD file format is not well defined. In particular, NAMD and
CHARMM use it differently. Currently, MDAnalysis tries to guess the
correct **format for the unitcell representation** but it can be
wrong. **Check the unitcell dimensions**, especially for triclinic
unitcells (see `Issue 187`_ and :attr:`Timestep.dimensions`). A
second potential issue are the units of time which are AKMA for the
:class:`DCDReader` (following CHARMM) but ps for NAMD. As a
workaround one can employ the configurable
:class:`MDAnalysis.coordinates.LAMMPS.DCDReader` for NAMD
trajectories.

See Also
--------
:mod:`MDAnalysis.coordinates.LAMMPS`
  module provides a more flexible DCD reader/writer.


The classes in this module are the reference implementations for the
Trajectory API.

.. _Issue 187:
   https://github.com/MDAnalysis/mdanalysis/issues/187


Classes
-------

.. autoclass:: Timestep
   :inherited-members:
.. autoclass:: DCDReader
   :inherited-members:
.. autoclass:: DCDWriter
   :inherited-members:

"""
from __future__ import absolute_import, division, print_function, unicode_literals
from six.moves import range

import os
import errno
import numpy as np
from numpy.lib.utils import deprecate
import struct
import types
import warnings

from ..core import flags
from .. import units as mdaunits  # use mdaunits instead of units to avoid a clash
from ..exceptions import NoDataError
from . import base
from . import core


class DCDReader(base.ReaderBase):
    """DCD Reader

    """
    format = 'DCD'
    flavor = 'CHARMM'
    units = {'time': 'AKMA', 'length': 'Angstrom'}

    def __init__(self, filename, convert_units=True, **kwargs):
        """Parameters
        ----------
        filename : str
            trajectory filename
        convert_units : bool (optional)
            convert units to MDAnalysis units
        **kwargs : dict
            General reader arguments.

        """
        super(DCDReader, self).__init__(filename,
                                        convert_units=convert_units,
                                        **kwargs)
        self._file = DCDFile(self.filename)
        self.n_atoms = self._file.n_atoms


        dt = self._file.delta

        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
        frame = self._file.read()
        self._frame = 0
        self._frame_to_ts(frame, self.ts)
        # these should only be initialized once
        self.ts.dt = dt
        self._file.seek(0)
        self.ts.dimensions = frame.unitcell
        if self.convert_units:
            self.convert_pos_from_native(self.ts.dimensions[:3])

    def close(self):
        """close reader"""
        self._file.close()

    @property
    def n_frames(self):
        """number of frames in trajectory"""
        return len(self._file)

    def _reopen(self):
        """reopen trajectory"""
        self.ts.frame = 0
        self._frame = -1
        self._file.close()
        self._file.open(self.filename.encode('utf-8'), 'r')

    def _read_frame(self, i):
        """read frame i"""
        self._frame = i - 1
        try:
            self._file.seek(i)
            timestep = self._read_next_timestep()
        except IOError:
            warnings.warn('seek failed, recalculating offsets and retrying')
            offsets = self._file.calc_offsets()
            self._file.set_offsets(offsets)
            self._read_offsets(store=True)
            self._file.seek(i)
            timestep = self._read_next_timestep()
        return timestep

    def _read_next_timestep(self, ts=None):
        """copy next frame into timestep"""
        if self._frame == self.n_frames - 1:
            raise IOError(errno.EIO, 'trying to go over trajectory limit')
        if ts is None:
            ts = self.ts
        frame = self._file.read()
        self._frame += 1
        self._frame_to_ts(frame, ts)
        return ts

    def Writer(self, filename, n_atoms=None, **kwargs):
        """Return writer for trajectory format"""
        if n_atoms is None:
            n_atoms = self.n_atoms
        return self._writer(filename, n_atoms=n_atoms, **kwargs)

    def _frame_to_ts(self, frame, ts):
        """convert a trr-frame to a mda TimeStep"""
        ts.time = self._file.tell() * self._file.delta
        ts.frame = self._frame
        ts.data['step'] = self._file.tell()

        ts.dimensions = frame.unitcell
        ts.positions = frame.x

        if self.convert_units:
            self.convert_pos_from_native(ts.dimensions[:3])
            self.convert_pos_from_native(ts.positions)

        return ts


class DCDWriter(base.WriterBase):
    """Base class for libmdaxdr file formats xtc and trr"""

    multiframe = True
    flavor = 'CHARMM'
    units = {'time': 'AKMA', 'length': 'Angstrom'}

    def __init__(self, filename, n_atoms, convert_units=True, **kwargs):
        """
        Parameters
        ----------
        filename : str
            filename of trajectory
        n_atoms : int
            number of atoms to be written
        convert_units : bool (optional)
            convert from MDAnalysis units to format specific units
        **kwargs : dict
            General writer arguments
        """
        self.filename = filename
        self._convert_units = convert_units
        self.n_atoms = n_atoms
        self._file = self._file(self.filename, 'w')

    def write_next_timestep(self, ts):
        """Write timestep object into trajectory.

        Parameters
        ----------
        ts: TimeStep

        See Also
        --------
        <FormatWriter>.write(AtomGroup/Universe/TimeStep)
        The normal write() method takes a more general input
        """
        xyz = ts.positions.copy()
        time = ts.time
        step = ts.frame
        dimensions = ts.dimensions

        if self._convert_units:
            xyz = self.convert_pos_to_native(xyz, inplace=False)
            dimensions = self.convert_dimensions_to_unitcell(ts, inplace=False)

        self._file.write(xyz=xyz, box=dimensions, step=step, natoms=xyz.shape[0],
                         charmm=1, time_step=ts.dt * step, ts_between_saves=1, remarks='test')

    def close(self):
        """close trajectory"""
        self._file.close()

    def __del__(self):
        self.close()
