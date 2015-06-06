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
from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np

from . import base
from . import core

class Timestep(base.Timestep):
    def _init_unitcell(self):
        return np.zeros((3, 3), dtype=np.float32, order='F')

    @property
    def dimensions(self):
        return core.triclinic_box(self._unitcell[0], self._unitcell[1], self._unitcell[2])

    @dimensions.setter
    def dimensions(self, new):
        self._unitcell[:] = core.triclinic_vectors(new)


class ConfigReader(base.SingleFrameReader):
    """DLPoly Config file Reader

    .. versionadded:: 0.10.1
    """
    format = 'DL_Config'
    units = {'length': 'Angstrom', 'velocity': 'Angstrom/ps'}
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
                velocities = []
            if levcfg == 2:
                forces = []

            # Read records infinitely
            while True:
                try:
                    line = inf.readline().strip()
                    if line == "":
                        break

                    try:
                        idx = int(line[8:])
                    except ValueError:  # dl_poly classic doesn't have this
                        pass
                    else:
                        ids.append(idx)

                    xyz = map(float, inf.readline().split())
                    coords.append(xyz)
                    if levcfg > 0:
                        vxyz = map(float, inf.readline().split())
                        velocities.append(vxyz)
                    if levcfg == 2:
                        fxyz = map(float, inf.readline().split())
                        forces.append(fxyz)
                except IOError:
                    break

        coords = np.array(coords, dtype=np.float32, order='F')
        if velocities:
            has_vels = True
            velocities = np.array(velocities, dtype=np.float32, order='F')
        if forces:
            has_forces = True
            forces = np.array(forces, dtype=np.float32, order='F')    
        self.numatoms = len(coords)

        if ids:
            # If we have indices then sort based on them
            ids = np.array(ids)
            order = np.argsort(ids)

            coords = coords[order]
            if has_vels:
                velocities = velocities[order]
            if has_forces:
                forces = forces[order]    

        # TODO: Merge with "new" style Timestep when finished
        ts = self.ts = self._Timestep(self.numatoms)
        ts._pos = coords
        if has_vels:
            ts._velocities = velocities
        if has_forces:
            ts._forces = forces

        if not imcon == 0:
            ts._unitcell[0][:] = cellx
            ts._unitcell[1][:] = celly
            ts._unitcell[2][:] = cellz
