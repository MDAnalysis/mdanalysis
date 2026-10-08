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

"""
PDBx (mmcif) files in MDAnalysis --- :mod:`MDAnalysis.coordinates.PDBx`
=======================================================================

Reads coordinates from a PDBx_ (mmcif) format file.  Will populate the Universe positions from the
``_atom_site.Cartn_x`` field in the PDBx file.  Will populate the unitcell dimensions from the ``_cell`` section.


.. _PDBx:
   https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/beginner’s-guide-to-pdb-structures-and-the-pdbx-mmcif-format
"""
import gemmi
import numpy as np

from . import base


class PDBxReader(base.SingleFrameReaderBase):
    format = ['cif', 'pdbx']
    units = {'time': None, 'length': 'Angstrom'}

    def _read_first_frame(self):
        doc = gemmi.cif.read(self.filename)

        block = doc.sole_block()

        coords = block.find('_atom_site.', ['Cartn_x', 'Cartn_y', 'Cartn_z'])
        self.natoms = len(coords)

        xyz = np.zeros((self.natoms, 3), dtype=np.float32)

        for i, (x, y, z) in enumerate(coords):
            xyz[i, :] = x, y, z

        ts = self.ts = base.Timestep.from_coordinates(xyz, **self._ts_kwargs)
        ts.frame = 0

        box = block.find('_cell.', ['length_a', 'length_b', 'length_c',
                                    'angle_alpha', 'angle_beta', 'angle_gamma'])
        if box:
            unitcell = np.zeros(6, dtype=np.float64)
            unitcell[:] = box[0]

            ts.dimensions = unitcell

        if self.convert_units:
            # in-place !
            self.convert_pos_from_native(self.ts._pos)
            if self.ts.dimensions is not None:
                self.convert_pos_from_native(self.ts.dimensions[:3])

        return ts
