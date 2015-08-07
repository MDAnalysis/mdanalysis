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


"""INPCRD structure files in MDAnalysis --- :mod:`MDAnalysis.coordinates.INPCRD`
================================================================================

Read and write coordinates in Amber_ coordinate/restart file (suffix
"inpcrd").

.. _Amber: http://ambermd.org/formats.html#restart
"""
from __future__ import absolute_import, division, print_function, unicode_literals

from . import base

class INPReader(base.SingleFrameReader):
    format = 'INPCRD'
    units = {'length': 'Angstrom'}

    def _read_first_frame(self):
        # Read header
        with open(self.filename, 'r') as inf:
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

            for p in xrange(self.n_atoms // 2):
                line = inf.readline()
                # each float is f12.7, 6 floats a line
                for i, dest in enumerate([(2*p, 0), (2*p, 1), (2*p, 2),
                                          (2*p + 1, 0), (2*p + 1, 1), (2*p + 1, 2)]):
                    self.ts._pos[dest] = float(line[i*12:(i+1)*12])
            # Read last coordinate if necessary            
            if self.n_atoms % 2:
                line = inf.readline()
                for i in range(3):
                    self.ts._pos[-1, i] = float(line[i*12:(i+1)*12])
