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

"""
DESRES file format --- :mod:`MDAnalysis.coordinates.DMS`
========================================================

Classes to read DESRES_ Molecular Structure file format (DMS_)
coordinate files (as used by the Desmond_ MD package).

.. _DESRES: http://www.deshawresearch.com
.. _Desmond: http://www.deshawresearch.com/resources_desmond.html
.. _DMS: http://www.deshawresearch.com/Desmond_Users_Guide-0.7.pdf
"""

import os
import errno
import warnings
from copy import deepcopy
import numpy
import sqlite3

import MDAnalysis
from . import base
import MDAnalysis.core.util as util
from MDAnalysis.coordinates.core import triclinic_box, triclinic_vectors


class Timestep(base.Timestep):
    def _init_unitcell(self):
        return {'x': numpy.zeros(3),
                'y': numpy.zeros(3),
                'z': numpy.zeros(3)}

    @property
    def dimensions(self):
        """unitcell dimensions (A, B, C, alpha, beta, gamma)"""
        x = self._unitcell['x']
        y = self._unitcell['y']
        z = self._unitcell['z']
        return triclinic_box(x, y, z)

    @dimensions.setter
    def dimensions(self, box):
        x, y, z = triclinic_vectors(box)
        cell = {'x': x, 'y': y, 'z': z}
        self._unitcell = cell


class DMSReader(base.SingleFrameReader):
    """
    Reads both coordinates and velocities.
    """
    format = 'DMS'
    units = {'time': None, 'length': 'A', 'velocity': 'A/ps'}
    _Timestep = Timestep

    def get_coordinates(self, cur):
        cur.execute('SELECT * FROM particle')
        particles = cur.fetchall()
        return [(p['x'], p['y'], p['z']) for p in particles]

    def get_particle_by_columns(self, cur, columns=['x', 'y', 'z']):
        cur.execute('SELECT * FROM particle')
        particles = cur.fetchall()
        return [tuple([p[c] for c in columns]) for p in particles]

    def get_global_cell(self, cur):
        cur.execute('SELECT * FROM global_cell')
        rows = cur.fetchall()
        assert len(rows) == 3
        x = [row["x"] for row in rows]
        y = [row["y"] for row in rows]
        z = [row["z"] for row in rows]
        return {'x': x, 'y': y, 'z': z}

    def _read_first_frame(self):
        coords_list = None
        velocities_list = None
        con = sqlite3.connect(self.filename)

        def dict_factory(cursor, row):
            d = {}
            for idx, col in enumerate(cursor.description):
                d[col[0]] = row[idx]
            return d

        with con:
            # This will return dictionaries instead of tuples, when calling cur.fetch() or fetchall()
            con.row_factory = dict_factory
            cur = con.cursor()
            coords_list = self.get_coordinates(cur)
            velocities_list = self.get_particle_by_columns(cur, columns=['vx', 'vy', 'vz'])
            unitcell = self.get_global_cell(cur)
        assert coords_list
        self.numatoms = len(coords_list)
        coords_list = numpy.array(coords_list)
        self.ts = self._Timestep(coords_list)
        self.ts.frame = 1  # 1-based frame number
        if velocities_list:
            # TODO: use a Timestep that knows about velocities such as TRR.Timestep or better, TRJ.Timestep
            velocities_arr = numpy.array(velocities_list, dtype=numpy.float32)
            if numpy.any(velocities_arr):
                self.ts._velocities = velocities_arr
                self.convert_velocities_from_native(self.ts._velocities)  # converts nm/ps to A/ps units
        # ts._unitcell layout is format dependent; Timestep.dimensions does the conversion
        self.ts._unitcell = unitcell
        if self.convert_units:
            self.convert_pos_from_native(self.ts._pos)  # in-place !
            self.convert_pos_from_native(self.ts._unitcell)  # in-place ! (all are lengths)
