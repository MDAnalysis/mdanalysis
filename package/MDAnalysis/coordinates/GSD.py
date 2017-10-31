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
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
"""GSD trajectory reader  --- :mod:`MDAnalysis.coordinates.GSD`
============================================================

Classes
-------

.. autoclass:: GSDReader
   :inherited-members:

"""
from __future__ import absolute_import, division

import numpy as np
from . import base
try :
    import gsd.hoomd
except :
    pass

class GSDReader(base.ReaderBase):
    """Reader for the GSD format.

    TODO: write docstring

    """
    format = 'GSD'
    units = {'time': None, 'length': None}
    _Timestep = base.Timestep

    def __init__(self, filename, **kwargs):
        """
        Parameters
        ----------
        filename : str
            trajectory filename
        **kwargs : dict
            General reader arguments.


        .. versionadded:: 0.17.0
        """
        super(GSDReader, self).__init__(filename, **kwargs)
        self.filename = filename
        self.open_trajectory()
        self.ts = self._read_next_timestep()

    def open_trajectory (self) :
        """opens the trajectory file using gsd.hoomd module"""
        self._frame = -1
        self._file = gsd.hoomd.open (self.filename,mode='rb')

    def close(self):
        """close reader"""
        self._file.close()

    @property
    def n_frames(self):
        """number of frames in trajectory"""
        return len(self._file)

    def _reopen(self):
        """reopen trajectory"""
        self.close()
        self.open_trajectory()

    def _read_next_timestep(self):
        """read next frame in trajectory"""
        # again here we assume that the number of particles remains fixed during
        # the trajectory
        self._frame += 1
        self.n_atoms = self._file[self._frame].particles.N
        # instantiate the Timestep
        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
        frame = self._file[self._frame]
        return self._frame_to_ts(frame, self.ts)

    def _frame_to_ts(self, frame, ts):
        """convert a gsd-frame to a :class:`TimeStep`"""
        ts.frame = self._frame

        # set frame box dimensions
        ts.dimensions = frame.configuration.box
        for i in range(3,6) :
            ts.dimensions[i] = np.arccos(ts.dimensions[i]) * 90.0 / np.pi / 2.0

        # set particle positions
        ts.positions = frame.particles.position

        return ts

    @property
    def dimensions(self):
        """unitcell dimensions (*A*, *B*, *C*, *alpha*, *beta*, *gamma*)
        """
        return self.ts.dimensions
