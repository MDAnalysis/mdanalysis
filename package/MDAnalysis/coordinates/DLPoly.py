# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
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

"""DL_Poly format reader :mod:`MDAnalysis.coordinates.DLPoly`
=============================================================

Read DL Poly_ format coordinate files

.. _Poly: http://www.stfc.ac.uk/SCD/research/app/ccg/software/DL_POLY/44516.aspx
"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import numpy as np

from . import base
from . import core

_DLPOLY_UNITS = {'length': 'Angstrom', 'velocity': 'Angstrom/ps', 'time': 'ps'}


class Timestep(base.Timestep):
    def _init_unitcell(self):
        return np.zeros((3, 3), dtype=np.float32, order='F')

    @property
    def dimensions(self):
        return core.triclinic_box(*self._unitcell)

    @dimensions.setter
    def dimensions(self, new):
        self._unitcell[:] = core.triclinic_vectors(new)


class ConfigReader(base.SingleFrameReader):
    """DLPoly Config file Reader

    .. versionadded:: 0.11.0
    """
    format = 'CONFIG'
    units = _DLPOLY_UNITS
    _Timestep = Timestep

    def _read_first_frame(self):
        with open(self.filename, 'r') as inf:
            self.title = inf.readline().strip()
            levcfg, imcon, megatm = map(int, inf.readline().split()[:3])
            if not imcon == 0:
                cellx = map(float, inf.readline().split())
                celly = map(float, inf.readline().split())
                cellz = map(float, inf.readline().split())

            ids = []
            coords = []
            if levcfg > 0:
                has_vels = True
                velocities = []
            else:
                has_vels = False
            if levcfg == 2:
                has_forces = True
                forces = []
            else:
                has_forces = False

            line = inf.readline().strip()
            # Read records infinitely
            while line:
                try:
                    idx = int(line[8:])
                except ValueError:  # dl_poly classic doesn't have this
                    pass
                else:
                    ids.append(idx)

                xyz = map(float, inf.readline().split())
                coords.append(xyz)
                if has_vels:
                    vxyz = map(float, inf.readline().split())
                    velocities.append(vxyz)
                if has_forces:
                    fxyz = map(float, inf.readline().split())
                    forces.append(fxyz)

                line = inf.readline().strip()

        coords = np.array(coords, dtype=np.float32, order='F')
        if has_vels:
            velocities = np.array(velocities, dtype=np.float32, order='F')
        if has_forces:
            forces = np.array(forces, dtype=np.float32, order='F')
        self.n_atoms = len(coords)

        if ids:
            # If we have indices then sort based on them
            ids = np.array(ids)
            order = np.argsort(ids)

            coords = coords[order]
            if has_vels:
                velocities = velocities[order]
            if has_forces:
                forces = forces[order]

        ts = self.ts = self._Timestep(self.n_atoms,
                                      velocities=has_vels,
                                      forces=has_forces,
                                      **self._ts_kwargs)
        ts._pos = coords
        if has_vels:
            ts._velocities = velocities
        if has_forces:
            ts._forces = forces

        if not imcon == 0:
            ts._unitcell[0][:] = cellx
            ts._unitcell[1][:] = celly
            ts._unitcell[2][:] = cellz

        ts.frame = 0


class HistoryReader(base.Reader):
    """Reads DLPoly format HISTORY files

    .. versionadded:: 0.11.0
    """
    format = 'HISTORY'
    units = _DLPOLY_UNITS
    _Timestep = Timestep

    def __init__(self, filename, **kwargs):
        super(HistoryReader, self).__init__(filename, **kwargs)

        # "private" file handle
        self._file = open(self.filename, 'r')
        self.title = self._file.readline().strip()
        self._levcfg, self._imcon, self.n_atoms = map(int, self._file.readline().split()[:3])
        self._has_vels = True if self._levcfg > 0 else False
        self._has_forces = True if self._levcfg == 2 else False

        self.ts = self._Timestep(self.n_atoms,
                                 velocities=self._has_vels,
                                 forces=self._has_forces,
                                 **self._ts_kwargs)
        self._read_next_timestep()

    def _read_next_timestep(self, ts=None):
        if ts is None:
            ts = self.ts

        line = self._file.readline()  # timestep line
        if not line.startswith('timestep'):
            raise IOError
        if not self._imcon == 0:
            ts._unitcell[0] = map(float, self._file.readline().split())
            ts._unitcell[1] = map(float, self._file.readline().split())
            ts._unitcell[2] = map(float, self._file.readline().split())

        # If ids are given, put them in here
        # and later sort by them
        ids = []

        for i in range(self.n_atoms):
            line = self._file.readline().strip()  # atom info line
            try:
                idx = int(line.split()[1])
            except IndexError:
                pass
            else:
                ids.append(idx)

            # Read in this order for now, then later reorder in place
            ts._pos[i] = map(float, self._file.readline().split())
            if self._has_vels:
                ts._velocities[i] = map(float, self._file.readline().split())
            if self._has_forces:
                ts._forces[i] = map(float, self._file.readline().split())

        if ids:
            ids = np.array(ids)
            # if ids aren't strictly sequential
            if not all(ids == (np.arange(self.n_atoms) + 1)):
                order = np.argsort(ids)
                ts._pos[:] = ts._pos[order]
                if self._has_vels:
                    ts._velocities[:] = ts._velocities[order]
                if self._has_forces:
                    ts._forces[:] = ts._forces[order]

        ts.frame += 1
        return ts

    def _read_frame(self, frame):
        """frame is 0 based, error checking is done in base.getitem"""
        self._file.seek(self._offsets[frame])
        self.ts.frame = frame  # gets +1'd in read_next_frame
        return self._read_next_timestep()

    @property
    def n_frames(self):
        try:
            return self._n_frames
        except AttributeError:
            self._n_frames = self._read_n_frames()
            return self._n_frames

    def _read_n_frames(self):
        """Read the number of frames, and the offset for each frame

        offset[i] - returns the offset in bytes to seek into the file to be
                    just before the frame starts
        """
        offsets = self._offsets = []

        with open(self.filename, 'r') as f:
            n_frames = 0

            f.readline()
            f.readline()
            position = f.tell()
            line = f.readline()
            while line.startswith('timestep'):
                offsets.append(position)
                n_frames += 1
                if not self._imcon == 0:  # box info
                    f.readline()
                    f.readline()
                    f.readline()
                for _ in range(self.n_atoms):
                    f.readline()
                    f.readline()
                    if self._has_vels:
                        f.readline()
                    if self._has_forces:
                        f.readline()
                position = f.tell()
                line = f.readline()

        return n_frames

    def rewind(self):
        self._reopen()
        self.next()

    def _reopen(self):
        self.close()
        self._file = open(self.filename, 'r')
        self._file.readline()  # header is 2 lines
        self._file.readline()
        self.ts.frame = 0

    def close(self):
        self._file.close()
