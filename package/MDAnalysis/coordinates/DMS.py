# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2011 Naveen Michaud-Agrawal,
#               Elizabeth J. Denning, Oliver Beckstein,
#               and contributors (see website for details)
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
#     N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
#     O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
#     Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
#     doi:10.1002/jcc.21787
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

from __future__ import with_statement

import os, errno
import warnings

import numpy, sqlite3

import MDAnalysis
import base
import MDAnalysis.core.util as util
from MDAnalysis.coordinates.core import triclinic_box, triclinic_vectors

from copy import deepcopy

class Timestep(base.Timestep):
        @property
        def dimensions(self):
                """unitcell dimensions (A, B, C, alpha, beta, gamma)

                GRO::

                  8.00170   8.00170   5.65806   0.00000   0.00000   0.00000   0.00000   4.00085   4.00085

                PDB::

                  CRYST1   80.017   80.017   80.017  60.00  60.00  90.00 P 1           1

                XTC: c.trajectory.ts._unitcell::

                  array([[ 80.00515747,   0.        ,   0.        ],
                         [  0.        ,  80.00515747,   0.        ],
                         [ 40.00257874,  40.00257874,  56.57218552]], dtype=float32)
                """
                # unit cell line (from http://manual.gromacs.org/current/online/gro.html)
                # v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)
                # 0     1     2      3     4     5     6    7     8
                x = self._unitcell['x']
                y = self._unitcell['y']
                z = self._unitcell['z']  # this ordering is correct! (checked it, OB)
                return triclinic_box(x,y,z)

class DMSReader(base.Reader):
        """
        Reads both coordinates and velocities.
        """
        format = 'DMS'
        units = {'time': None, 'length': 'A', 'velocity': 'A/ps'}
        _Timestep = Timestep


        def get_coordinates(self, cur):
            cur.execute('SELECT * FROM particle')
            particles = cur.fetchall()
            return [ (p['x'], p['y'], p['z'])  for p in particles]

        def get_particle_by_columns(self, cur, columns=['x', 'y', 'z']):
            cur.execute('SELECT * FROM particle')
            particles = cur.fetchall()
            return [ tuple([p[c] for c in columns ])  for p in particles]

        def get_global_cell(self, cur):
            cur.execute('SELECT * FROM global_cell')
            rows = cur.fetchall()
            assert len(rows) == 3
            x = [row["x"] for row in rows]
            y = [row["y"] for row in rows]
            z = [row["z"] for row in rows]
            return {'x': x, 'y': y, 'z': z}

        def __init__(self,filename,convert_units=None,**kwargs):
                if convert_units is None:
                        convert_units = MDAnalysis.core.flags['convert_gromacs_lengths']
                self.convert_units = convert_units  # convert length and time to base units

                coords_list = None
                velocities_list = None
                con = sqlite3.connect(filename)

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
                if velocities_list: #perform this operation only if velocities are present in coord file
                    # TODO: use a Timestep that knows about velocities such as TRR.Timestep or better, TRJ.Timestep
                    velocities_arr = numpy.array(velocities_list, dtype=numpy.float32)
                    if numpy.any(velocities_arr):
                        self.ts._velocities = velocities_arr
                        self.convert_velocities_from_native(self.ts._velocities) #converts nm/ps to A/ps units
                # ts._unitcell layout is format dependent; Timestep.dimensions does the conversion
                self.ts._unitcell = unitcell
                if self.convert_units:
                        self.convert_pos_from_native(self.ts._pos)             # in-place !
                        self.convert_pos_from_native(self.ts._unitcell)        # in-place ! (all are lengths)
                self.numframes = 1
                self.fixed = 0
                self.skip = 1
                self.periodic = False
                self.delta = 0
                self.skip_timestep = 1

        def Writer(self, filename, **kwargs):
                raise NotImplementedError

        def __iter__(self):
                yield self.ts  # Just a single frame
                raise StopIteration

        def _read_frame(self, frame):
                if frame != 0:
                        raise IndexError("DMSReader only handles a single frame at frame index 0")
                return self.ts

        def _read_next_timestep(self):
                # CRD files only contain a single frame
                raise IOError
