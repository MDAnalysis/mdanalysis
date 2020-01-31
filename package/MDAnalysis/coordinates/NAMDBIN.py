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


"""NAMDBIN files format --- :mod:`MDAnalysis.coordinates.NAMDBIN`
================================================================================

Read/Write coordinates in NAMD_ double-precision binary file (suffix "coor" or "namdbin").

.. _NAMD: https://www.ks.uiuc.edu/Research/namd/2.10/ug/node11.html


Classes
-------

.. autoclass:: NAMDBINReader
   :members:

.. autoclass:: NAMDBINWriter
   :members:

"""
from __future__ import absolute_import, division, print_function, unicode_literals

from struct import unpack, pack
import numpy as np

from . import base
from ..lib import util


class NAMDBINReader(base.SingleFrameReaderBase):
    """Reader for NAMD binary files."""

    format = ['NAMDBIN','COOR']
    units = {'length': 'Angstrom'}

    def _read_first_frame(self):
        # Read header
        with open(self.filename, 'rb') as namdbin:
            self.n_atoms = int(unpack('i', namdbin.read(4))[0])

            self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
            self.ts.frame = 0
            coord_double = unpack('{:d}d'.format(3*self.n_atoms),
                                  namdbin.read(3*self.n_atoms*8))
            self.ts._pos[:] = np.array(
                coord_double, float).reshape(self.n_atoms, 3)

    @staticmethod
    def parse_n_atoms(filename, **kwargs):
        with open(filename, 'rb') as namdbin:
            n_atoms = int(unpack('i', namdbin.read(4))[0])
        return n_atoms

    def Writer(self, filename, **kwargs):
        """Returns a NAMDBINWriter for *filename*.

        Parameters
        ----------
        filename: str
            filename of the output NAMDBIN file

        Returns
        -------
        :class:`NAMDBINWriter`

        """
        return NAMDBINWriter(filename, **kwargs)


class NAMDBINWriter(base.WriterBase):
    """Writer for NAMD binary files."""
    format = ['NAMDBIN','COOR']
    units = {'time': None, 'length': 'Angstrom'}

    def __init__(self, filename, n_atoms=None, **kwargs):
        """
        Parameters
        ----------
        filename : str or :class:`~MDAnalysis.lib.util.NamedStream`
             name of the output file or a stream
        n_atoms  : int
            number of atoms for the output coordinate
        """
        self.n_atoms = n_atoms
        self.filename = util.filename(filename, ext='namdbin')

    def write(self, obj):
        """Write obj at current trajectory frame to file.

        Parameters
        ----------
        obj : :class:`~MDAnalysis.core.groups.AtomGroup` or :class:`~MDAnalysis.core.universe.Universe` or a :class:`Timestep` 
              write coordinate information associate with `obj`
        """

        if isinstance(obj, base.Timestep):
            n_atoms = obj.n_atoms
            coor = obj.positions.reshape(n_atoms*3)
        elif hasattr(obj, 'atoms'):  # AtomGroup or Universe
            atoms = obj.atoms  # make sure to use atoms (Issue 46)
            n_atoms = len(atoms)
            # can write from obj == Universe (Issue 49)
            coor = atoms.positions.reshape(n_atoms*3)
        else:
            raise TypeError

        with util.openany(self.filename, 'wb') as namdbin:
            # Write NUMATOMS
            namdbin.write(pack('i', n_atoms))
            # Write Coordinate
            namdbin.write(pack('{:d}d'.format(len(coor)), *coor))
