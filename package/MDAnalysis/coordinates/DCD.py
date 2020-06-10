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
"""DCD trajectory I/O  --- :mod:`MDAnalysis.coordinates.DCD`
============================================================

Classes to read and write DCD binary trajectories, the format used by
CHARMM, NAMD, and also LAMMPS. Trajectories can be read regardless of
system-endianness as this is auto-detected.

Generally, DCD trajectories produced by any code can be read (with the
:class:`DCDReader`) although there can be issues with the unitcell (simulation
box) representation (see :attr:`DCDReader.dimensions`). DCDs can also be
written but the :class:`DCDWriter` follows recent NAMD/VMD convention for the
unitcell but still writes AKMA time. Reading and writing these trajectories
within MDAnalysis will work seamlessly but if you process those trajectories
with other tools you might need to watch out that time and unitcell dimensions
are correctly interpreted.


See Also
--------
:mod:`MDAnalysis.coordinates.LAMMPS`
  module provides a more flexible DCD reader/writer.


.. _Issue 187:
   https://github.com/MDAnalysis/mdanalysis/issues/187


Classes
-------

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
import struct
import types
import warnings

from .. import units as mdaunits  # use mdaunits instead of units to avoid a clash
from ..exceptions import NoDataError
from . import base, core
from ..lib.formats.libdcd import DCDFile
from ..lib.mdamath import triclinic_box


class DCDReader(base.ReaderBase):
    """Reader for the DCD format.

    DCD is used by NAMD, CHARMM and LAMMPS as the default trajectory format.
    The DCD file format is not well defined. In particular, NAMD and CHARMM use
    it differently. Currently, MDAnalysis tries to guess the correct **format
    for the unitcell representation** but it can be wrong. **Check the unitcell
    dimensions**, especially for triclinic unitcells (see `Issue 187`_). DCD
    trajectories produced by CHARMM and NAMD( >2.5) record time in AKMA units.
    If other units have been recorded (e.g., ps) then employ the configurable
    :class:MDAnalysis.coordinates.LAMMPS.DCDReader and set the time unit as a
    optional argument. You can find a list of units used in the DCD formats on
    the MDAnalysis `wiki`_.


    MDAnalysis always uses ``(*A*, *B*, *C*, *alpha*, *beta*, *gamma*)`` to
    represent the unit cell. Lengths *A*, *B*, *C* are in the MDAnalysis length
    unit (Å), and angles are in degrees.

    The ordering of the angles in the unitcell is the same as in recent
    versions of VMD's DCDplugin_ (2013), namely the `X-PLOR DCD format`_: The
    original unitcell is read as ``[A, gamma, B, beta, alpha, C]`` from the DCD
    file. If any of these values are < 0 or if any of the angles are > 180
    degrees then it is assumed it is a new-style CHARMM unitcell (at least
    since c36b2) in which box vectors were recorded.

    .. warning::
        The DCD format is not well defined. Check your unit cell
        dimensions carefully, especially when using triclinic boxes.
        Different software packages implement different conventions and
        MDAnalysis is currently implementing the newer NAMD/VMD convention
        and tries to guess the new CHARMM one. Old CHARMM trajectories might
        give wrong unitcell values. For more details see `Issue 187`_.

    .. _`X-PLOR DCD format`: http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html
    .. _Issue 187: https://github.com/MDAnalysis/mdanalysis/issues/187
    .. _DCDplugin: http://www.ks.uiuc.edu/Research/vmd/plugins/doxygen/dcdplugin_8c-source.html#l00947
    .. _wiki: https://github.com/MDAnalysis/mdanalysis/wiki/FileFormats#dcd
    """
    format = 'DCD'
    flavor = 'CHARMM'
    units = {'time': 'AKMA', 'length': 'Angstrom'}

    def __init__(self, filename, convert_units=True, dt=None, **kwargs):
        """
        Parameters
        ----------
        filename : str
            trajectory filename
        convert_units : bool (optional)
            convert units to MDAnalysis units
        dt : float (optional)
            overwrite time delta stored in DCD
        **kwargs : dict
            General reader arguments.


        .. versionchanged:: 0.17.0
           Changed to use libdcd.pyx library and removed the correl function
        """
        super(DCDReader, self).__init__(
            filename, convert_units=convert_units, **kwargs)
        self._file = DCDFile(self.filename)
        self.n_atoms = self._file.header['natoms']

        delta = mdaunits.convert(self._file.header['delta'],
                                 self.units['time'], 'ps')
        if dt is None:
            dt = delta * self._file.header['nsavc']
        self.skip_timestep = self._file.header['nsavc']

        self._ts_kwargs['dt'] = dt
        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
        frame = self._file.read()
        # reset trajectory
        if self._file.n_frames > 1:
            self._file.seek(1)
        else:
            self._file.seek(0)
        self._frame = 0
        self.ts = self._frame_to_ts(frame, self.ts)
        # these should only be initialized once
        self.ts.dt = dt
        if self.convert_units:
            self.convert_pos_from_native(self.ts.dimensions[:3])

    @staticmethod
    def parse_n_atoms(filename, **kwargs):
        with DCDFile(filename) as f:
            n_atoms = f.header['natoms']
        return n_atoms

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
        self._file.open('r')

    def _read_frame(self, i):
        """read frame i"""
        self._frame = i - 1
        self._file.seek(i)
        return self._read_next_timestep()

    def _read_next_timestep(self, ts=None):
        """copy next frame into timestep"""
        if self._frame == self.n_frames - 1:
            raise IOError('trying to go over trajectory limit')
        if ts is None:
            # use a copy to avoid that ts always points to the same reference
            # removing this breaks lammps reader
            ts = self.ts.copy()
        frame = self._file.read()
        self._frame += 1
        ts = self._frame_to_ts(frame, ts)
        self.ts = ts
        return ts

    def Writer(self, filename, n_atoms=None, **kwargs):
        """Return writer for trajectory format"""
        if n_atoms is None:
            n_atoms = self.n_atoms
        return DCDWriter(
            filename,
            n_atoms=n_atoms,
            dt=self.ts.dt,
            convert_units=self.convert_units,
            **kwargs)

    def _frame_to_ts(self, frame, ts):
        """convert a dcd-frame to a :class:`TimeStep`"""
        ts.frame = self._frame
        ts.time = (ts.frame + self._file.header['istart']/self._file.header['nsavc']) * self.ts.dt
        ts.data['step'] = self._file.tell()

        # The original unitcell is read as ``[A, gamma, B, beta, alpha, C]``
        _ts_order = [0, 2, 5, 4, 3, 1]
        uc = np.take(frame.unitcell, _ts_order)

        pi_2 = np.pi / 2
        if (-1.0 <= uc[3] <= 1.0) and (-1.0 <= uc[4] <= 1.0) and (
                -1.0 <= uc[5] <= 1.0):
            # This file was generated by Charmm, or by NAMD > 2.5, with the
            # angle cosines of the periodic cell angles written to the DCD
            # file. This formulation improves rounding behavior for orthogonal
            # cells so that the angles end up at precisely 90 degrees, unlike
            # acos(). (changed in MDAnalysis 0.9.0 to have NAMD ordering of the
            # angles; see Issue 187) */
            uc[3] = 90.0 - np.arcsin(uc[3]) * 90.0 / pi_2
            uc[4] = 90.0 - np.arcsin(uc[4]) * 90.0 / pi_2
            uc[5] = 90.0 - np.arcsin(uc[5]) * 90.0 / pi_2
        # heuristic sanity check: uc = A,B,C,alpha,beta,gamma
        elif np.any(uc < 0.) or np.any(uc[3:] > 180.):
            # might be new CHARMM: box matrix vectors
            H = frame.unitcell.copy()
            e1, e2, e3 = H[[0, 1, 3]], H[[1, 2, 4]], H[[3, 4, 5]]
            uc = triclinic_box(e1, e2, e3)
        else:
            # This file was likely generated by NAMD 2.5 and the periodic cell
            # angles are specified in degrees rather than angle cosines.
            pass

        ts.dimensions = uc
        ts.positions = frame.xyz

        if self.convert_units:
            self.convert_pos_from_native(ts.dimensions[:3])
            self.convert_pos_from_native(ts.positions)

        return ts

    @property
    def dimensions(self):
        """unitcell dimensions (*A*, *B*, *C*, *alpha*, *beta*, *gamma*)
        """
        return self.ts.dimensions

    @property
    def dt(self):
        """timestep between frames"""
        return self.ts.dt

    def timeseries(self,
                   asel=None,
                   start=None,
                   stop=None,
                   step=None,
                   order='afc'):
        """Return a subset of coordinate data for an AtomGroup

        Parameters
        ----------
        asel : :class:`~MDAnalysis.core.groups.AtomGroup`
            The :class:`~MDAnalysis.core.groups.AtomGroup` to read the
            coordinates from. Defaults to None, in which case the full set of
            coordinate data is returned.
        start : int (optional)
            Begin reading the trajectory at frame index `start` (where 0 is the
            index of the first frame in the trajectory); the default ``None``
            starts at the beginning.
        stop : int (optional)
            End reading the trajectory at frame index `stop`-1, i.e, `stop` is
            excluded. The trajectory is read to the end with the default
            ``None``.
        step : int (optional)
            Step size for reading; the default ``None`` is equivalent to 1 and
            means to read every frame.
        order : str (optional)
            the order/shape of the return data array, corresponding
            to (a)tom, (f)rame, (c)oordinates all six combinations
            of 'a', 'f', 'c' are allowed ie "fac" - return array
            where the shape is (frame, number of atoms,
            coordinates)


        .. versionchanged:: 1.0.0
           `skip` and `format` keywords have been removed.
        """

        start, stop, step = self.check_slice_indices(start, stop, step)

        if asel is not None:
            if len(asel) == 0:
                raise NoDataError(
                    "Timeseries requires at least one atom to analyze")
            atom_numbers = list(asel.indices)
        else:
            atom_numbers = list(range(self.n_atoms))

        frames = self._file.readframes(
            start, stop, step, order=order, indices=atom_numbers)
        return frames.xyz


class DCDWriter(base.WriterBase):
    """DCD Writer class

    The writer follows recent NAMD/VMD convention for the unitcell (box lengths
    in Å and angle-cosines, ``[A, cos(gamma), B, cos(beta), cos(alpha), C]``)
    and writes positions in Å and time in AKMA time units.

    """
    format = 'DCD'
    multiframe = True
    flavor = 'NAMD'
    units = {'time': 'AKMA', 'length': 'Angstrom'}

    def __init__(self,
                 filename,
                 n_atoms,
                 convert_units=True,
                 step=1,
                 dt=1,
                 remarks='',
                 nsavc=1,
                 istart=0,
                 **kwargs):
        """Parameters
        ----------
        filename : str
            filename of trajectory
        n_atoms : int
            number of atoms to be written
        convert_units : bool (optional)
            convert from MDAnalysis units to format specific units
        step : int (optional)
            number of steps between frames to be written
        dt : float (optional)
            time between two frames. If ``None`` guess from first written
            :class:`TimeStep`
        remarks : str (optional)
            remarks to be stored in DCD. Shouldn't be more then 240 characters
        nsavc : int (optional)
            DCD files can also store the number of integrator time steps that
            correspond to the interval between two frames as nsavc (i.e., every
            how many MD steps is a frame saved to the DCD). By default, this
            number is just set to one and this should be sufficient for almost
            all cases but if required, `nsavc` can be changed.
        istart : int (optional)
            starting frame number in integrator timesteps. CHARMM defaults to
            `nsavc`, i.e., start at frame number 1 = `istart` / `nsavc`. The value
            ``None`` will set `istart` to `nsavc` (the CHARMM default).
            The MDAnalysis default is 0 so that the frame number and time of the first
            frame is 0.
        **kwargs : dict
            General writer arguments

        """
        self.filename = filename
        self._convert_units = convert_units
        if n_atoms is None:
            raise ValueError("n_atoms argument is required")
        self.n_atoms = n_atoms
        self._file = DCDFile(self.filename, 'w')
        self.step = step
        self.dt = dt
        dt = mdaunits.convert(dt, 'ps', self.units['time'])
        delta = float(dt) / nsavc
        istart = istart if istart is not None else nsavc
        self._file.write_header(
            remarks=remarks,
            natoms=self.n_atoms,
            nsavc=nsavc,
            delta=delta,
            is_periodic=1,
            istart=istart)

    def _write_next_frame(self, ag):
        """Write information associated with ``obj`` at current frame into trajectory

        Parameters
        ----------
        ag : AtomGroup or Universe

        See Also
        --------
        :meth:`DCDWriter.write`  takes a more general input


        .. deprecated:: 1.0.0
           Deprecated use of Timestep as argument. To be removed in version 2.0
        .. versionchanged:: 1.0.0
           Added ability to pass AtomGroup or Universe.
           Renamed from `write_next_timestep` to `_write_next_frame`.
        """
        if isinstance(ag, base.Timestep):
            ts = ag
        else:
            try:
                ts = ag.ts
            except AttributeError:
                try:
                    # Universe?
                    ts = ag.trajectory.ts
                except AttributeError:
                    raise TypeError("No Timestep found in ag argument")
        xyz = ts.positions.copy()
        dimensions = ts.dimensions.copy()

        if self._convert_units:
            xyz = self.convert_pos_to_native(xyz, inplace=True)
            dimensions = self.convert_dimensions_to_unitcell(ts, inplace=True)

        # we only support writing charmm format unit cell info
        # The DCD unitcell is written as ``[A, gamma, B, beta, alpha, C]``
        _ts_order = [0, 5, 1, 4, 3, 2]
        box = np.take(dimensions, _ts_order)

        self._file.write(xyz=xyz, box=box)

    def close(self):
        """close trajectory"""
        self._file.close()

    def __del__(self):
        self.close()
