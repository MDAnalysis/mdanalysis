# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
# Beckstein and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""MMTF trajectory reader --- :mod:`MDAnalysis.coordinates.MMTF`
================================================================

Reads coordinates data from the `Macromolecular Transmission Format
(MMTF) format`_.  This should generally be a quicker alternative to PDB.

.. versionadded:: 0.16.0

Classes
-------

.. autoclass:: MMTFReader
   :members:


.. _Macromolecular Transmission Format (MMTF) format: https://mmtf.rcsb.org/
"""

import mmtf

from . import base


def _parse_mmtf(fn):
    if fn.endswith('gz'):
        return mmtf.parse_gzip(fn)
    else:
        return mmtf.parse(fn)


class MMTFReader(base.SingleFrameReader):
    format = 'MMTF'

    def _read_first_frame(self):
        # TOOD: Check units?
        top = _parse_mmtf(self.filename)
        self.n_atoms = top.num_atoms

        self.ts = ts = self._Timestep(self.n_atoms,
                                      **self._ts_kwargs)
        ts._pos[:, 0] = top.x_coord_list
        ts._pos[:, 1] = top.y_coord_list
        ts._pos[:, 2] = top.z_coord_list
        ts.dimensions = top.unit_cell

        return ts
