# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
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

Read/Write coordinates in `NAMD double-precision binary format`_ (suffix "coor" or "namdbin").

.. _`NAMD double-precision binary format` : https://www.ks.uiuc.edu/Research/namd/2.10/ug/node11.html#SECTION00061500000000000000


Classes
-------

.. autoclass:: NAMDBINReader
   :members:

.. autoclass:: NAMDBINWriter
   :members:

"""
from struct import pack
import numpy as np

from . import base
from ..lib import util


class NAMDBINReader(base.SingleFrameReaderBase):
    """Reader for NAMD binary coordinate files.


    .. versionadded:: 1.0.0
    """

    format = ["COOR", "NAMDBIN"]
    units = {"length": "Angstrom"}

    def _read_first_frame(self):
        # Read header
        with open(self.filename, "rb") as namdbin:
            self.n_atoms = np.fromfile(namdbin, dtype=np.int32, count=1)[0]
            self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
            self.ts.frame = 0
            coord_double = np.fromfile(
                namdbin, dtype=np.float64, count=self.n_atoms * 3
            )
            self.ts._pos[:] = np.array(coord_double, float).reshape(
                self.n_atoms, 3
            )

    @staticmethod
    def parse_n_atoms(filename, **kwargs):
        with open(filename, "rb") as namdbin:
            n_atoms = np.fromfile(namdbin, dtype=np.int32, count=1)[0]
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
    """Writer for NAMD binary coordinate files.


    Note
    ----
    * Does not handle writing to bz2 or gz compressed file types.


    .. versionadded:: 1.0.0
    """

    format = ["COOR", "NAMDBIN"]
    units = {"time": None, "length": "Angstrom"}

    def __init__(self, filename, n_atoms=None, **kwargs):
        """
        Parameters
        ----------
        filename : str or :class:`~MDAnalysis.lib.util.NamedStream`
             name of the output file or a stream
        n_atoms  : int
            number of atoms for the output coordinate
        """
        self.filename = util.filename(filename)

    def _write_next_frame(self, obj):
        """Write information associated with ``obj`` at current frame into
        trajectory


        Parameters
        ----------
        obj : :class:`~MDAnalysis.core.groups.AtomGroup` or
              :class:`~MDAnalysis.core.universe.Universe`
              write coordinate information associated with `obj`


        .. versionchanged:: 1.0.0
           Renamed from `write` to `_write_next_frame`.
        .. versionchanged:: 2.0.0
           Deprecated support for Timestep argument has now been removed.
           Use AtomGroup or Universe as an input instead.
        """
        if hasattr(obj, "atoms"):  # AtomGroup or Universe
            atoms = obj.atoms
            n_atoms = len(atoms)
            coor = atoms.positions.reshape(n_atoms * 3)
        else:
            errmsg = "Input obj is neither an AtomGroup or Universe"
            raise TypeError(errmsg) from None

        with util.openany(self.filename, "wb") as namdbin:
            # Write NUMATOMS
            namdbin.write(pack("i", n_atoms))
            # Write Coordinate
            namdbin.write(pack("{:d}d".format(len(coor)), *coor))
