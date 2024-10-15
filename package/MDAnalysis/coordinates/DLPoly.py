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

"""DL_Poly format reader :mod:`MDAnalysis.coordinates.DLPoly`
=============================================================

Read DL Poly_ format coordinate files

.. _Poly: http://www.stfc.ac.uk/SCD/research/app/ccg/software/DL_POLY/44516.aspx
"""
import numpy as np
import warnings

from . import base
from . import core
from ..lib import util
from ..lib.util import cached, store_init_arguments

_DLPOLY_UNITS = {'length': 'Angstrom', 'velocity': 'Angstrom/ps', 'time': 'ps'}


class ConfigReader(base.SingleFrameReaderBase):
    """DLPoly Config file Reader


    .. versionadded:: 0.11.0
    .. versionchanged:: 2.0.0
       coordinates, velocities, and forces are no longer stored in 'F' memory
       layout, instead now using the numpy default of 'C'.
    """
    format = 'CONFIG'
    units = _DLPOLY_UNITS

    def _read_first_frame(self):
        unitcell = np.zeros((3, 3), dtype=np.float32)

        with open(self.filename, 'r') as inf:
            self.title = inf.readline().strip()
            levcfg, imcon, megatm = np.int64(inf.readline().split()[:3])
            if not imcon == 0:
                unitcell[0][:] = inf.readline().split()
                unitcell[1][:] = inf.readline().split()
                unitcell[2][:] = inf.readline().split()

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

                xyz = np.float32(inf.readline().split())
                coords.append(xyz)
                if has_vels:
                    vxyz = np.float32(inf.readline().split())
                    velocities.append(vxyz)
                if has_forces:
                    fxyz = np.float32(inf.readline().split())
                    forces.append(fxyz)

                line = inf.readline().strip()

        coords = np.array(coords, dtype=np.float32)
        if has_vels:
            velocities = np.array(velocities, dtype=np.float32)
        if has_forces:
            forces = np.array(forces, dtype=np.float32)
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
            ts.dimensions = core.triclinic_box(*unitcell)

        ts.frame = 0


class HistoryReader(base.ReaderBase):
    """Reads DLPoly format HISTORY files

    .. versionadded:: 0.11.0
    """
    format = 'HISTORY'
    units = _DLPOLY_UNITS

    @store_init_arguments
    def __init__(self, filename, **kwargs):
        super(HistoryReader, self).__init__(filename, **kwargs)
        self._cache = {}

        # "private" file handle
        self._file = util.anyopen(self.filename, 'r')
        self.title = self._file.readline().strip()
        header = np.int64(self._file.readline().split()[:3])
        self._levcfg, self._imcon, self.n_atoms = header
        self._has_vels = True if self._levcfg > 0 else False
        self._has_forces = True if self._levcfg == 2 else False

        rwnd = self._file.tell()
        self._file.readline()
        if (len(self._file.readline().split())) == 3:
            self._has_cell = True
        else:
            self._has_cell = False
        self._file.seek(rwnd)

        self.ts = self._Timestep(self.n_atoms,
                                 velocities=self._has_vels,
                                 forces=self._has_forces,
                                 **self._ts_kwargs)
        self._read_next_timestep()

    def _read_next_timestep(self, ts=None):
        if ts:
            warnings.warn("ts argument to _read_next_timestep is deprecated  as of 2.7.0 and will be removed in 3.0.0, see #3928")
        
        if ts is None:
            ts = self.ts

        line = self._file.readline()  # timestep line
        if not line.startswith('timestep'):
            raise IOError

        if self._has_cell:
            unitcell = np.zeros((3, 3))
            unitcell[0] = self._file.readline().split()
            unitcell[1] = self._file.readline().split()
            unitcell[2] = self._file.readline().split()
            ts.dimensions = core.triclinic_box(*unitcell)            

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
            ts._pos[i] = self._file.readline().split()
            if self._has_vels:
                ts._velocities[i] = self._file.readline().split()
            if self._has_forces:
                ts._forces[i] = self._file.readline().split()
            i += 1

        if ids:
            ids = np.array(ids)
            # if ids aren't strictly sequential
            if not np.all(ids == (np.arange(self.n_atoms) + 1)):
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
        self.ts.frame = frame - 1  # gets +1'd in read_next_frame
        return self._read_next_timestep()

    @property
    @cached('n_frames')
    def n_frames(self):
        # Second line is traj_key, imcom, n_atoms, n_frames, n_records
        offsets = []

        with open(self.filename, 'r') as f:
            f.readline()
            f.readline()
            position = f.tell()
            line = f.readline()
            while line.startswith('timestep'):
                offsets.append(position)
                if self._has_cell:
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

        self._offsets = offsets
        return len(self._offsets)

    def _reopen(self):
        self.close()
        self._file = open(self.filename, 'r')
        self._file.readline()  # header is 2 lines
        self._file.readline()
        self.ts.frame = -1

    def close(self):
        self._file.close()
