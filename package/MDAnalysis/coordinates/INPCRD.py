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


"""INPCRD structure files in MDAnalysis --- :mod:`MDAnalysis.coordinates.INPCRD`
================================================================================

Read coordinates in Amber_ coordinate/restart file (suffix "inpcrd").

.. _Amber: https://ambermd.org/FileFormats.php


Classes
-------

.. autoclass:: INPReader
   :members:

"""

from . import base


class INPReader(base.SingleFrameReaderBase):
    """Reader for Amber restart files."""

    format = ["INPCRD", "RESTRT"]
    units = {"length": "Angstrom"}

    def _read_first_frame(self):
        # Read header
        with open(self.filename, "r") as inf:
            self.title = inf.readline().strip()
            line = inf.readline().split()
            self.n_atoms = int(line[0])

            self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
            try:
                time = float(line[1])
            except IndexError:
                pass
            else:
                self.ts.time = time
            self.ts.frame = 0

            for p in range(self.n_atoms // 2):
                line = inf.readline()
                # each float is f12.7, 6 floats a line
                for i, dest in enumerate(
                    [
                        (2 * p, 0),
                        (2 * p, 1),
                        (2 * p, 2),
                        (2 * p + 1, 0),
                        (2 * p + 1, 1),
                        (2 * p + 1, 2),
                    ]
                ):
                    self.ts._pos[dest] = float(line[i * 12 : (i + 1) * 12])
            # Read last coordinate if necessary
            if self.n_atoms % 2:
                line = inf.readline()
                for i in range(3):
                    self.ts._pos[-1, i] = float(line[i * 12 : (i + 1) * 12])

    @staticmethod
    def parse_n_atoms(filename, **kwargs):
        with open(filename, "r") as f:
            f.readline()
            n_atoms = int(f.readline().split()[0])
        return n_atoms
