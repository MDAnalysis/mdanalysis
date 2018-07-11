# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""
Dummy coordinate reader
=======================


Classes
-------

.. autoclass:: DummyReader
   :members:


"""
from __future__ import absolute_import


import numpy as np

from .base import SingleFrameReaderBase


class DummyReader(SingleFrameReaderBase):
    """Basic Reader which does not read from any file

    .. versionadded:: 0.17.0
    """
    format = 'dummy'

    def __init__(self, n_atoms=None, velocities=False, forces=False):
        self.n_atoms = n_atoms
        self.filename = 'DummyReader'
        self.n_frames = 1
        self._read_first_frame(velocities, forces)
        self._transformations = []

    def _read_first_frame(self, velocities=False, forces=False):
        ts = self.ts = self._Timestep(self.n_atoms, positions=True,
                                      velocities=velocities, forces=forces)
        return ts
