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
"""Dummy Reader ---:mod:`MDAnalysis.coordinates.dummy`

Doesn't actually read anything! Designed to allow systems
to be built from nothing.

.. versionadded:: 0.13.0
"""
from . import base


class DummyReader(base.SingleFrameReader):
    format = 'Dummy'
    units = {}

    def _read_first_frame(self):
        # Uses self.filename as hack for number of atoms
        self.n_atoms = self.filename

        # Start with not even positions
        self.ts = self._Timestep(self.n_atoms,
                                 positions=False,
                                 **self._ts_kwargs)
