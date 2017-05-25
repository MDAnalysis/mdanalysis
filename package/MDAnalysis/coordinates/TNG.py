# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.mdanalysis.org
# Copyright (c) 2006-2016 The MDAnalysis Development Team and contributors
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
"""\
TNG trajectory files --- :mod:`MDAnalysis.coordinates.TNG`
==========================================================

Read and write GROMACS TNG trajectories.

See Also
--------
MDAnalysis.coordinates.XTC: Read and write GROMACS XTC trajectory files.
MDAnalysis.coordinates.TRR: Read and write GROMACS TRR trajectory files.
"""

from __future__ import absolute_import

import errno
import numpy as np
import pytng

from . import base
from ..lib.mdamath import triclinic_box


class TNGReader(base.ReaderBase):
    """Reader for the Gromacs TNG format.

    The Gromacs TNG trajectory format is a lossless format. The TNG format can
    store *velocoties* and *forces* in addition to the coordinates.

    """

    format = 'TNG'
    units = {'time': 'ps', 'length': 'nm', 'velocity': 'nm/ps',
             'force': 'kJ/(mol*nm)'}
    _writer = TRRWriter
    _file = TRRFile

    def __init__(self, filename, convert_units=True, **kwargs):
        """
        Parameters
        ----------
        filename : str
            trajectory filename
        convert_units : bool (optional)
            convert units to MDAnalysis units
        **kwargs : dict
            General reader arguments.

        """
        super(TNGReader, self).__init__(filename,
                                        convert_units=convert_units,
                                        **kwargs)
        self._tng = pytng.TNGFile(self.filename)

        self.n_atoms = self._tng.n_atoms

        frame = self._tng.read()
        if len(self._tng) > 1:
            frame_2 = self._tng.read()
            dt = frame_2.time - frame.time
        else:
            dt = 0

        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
        self._frame = 0
        self._frame_to_ts(frame, self.ts)
        # these should only be initialized once
        self.ts.dt = dt
        self.ts.dimensions = triclinic_box(*frame.box)
        if self.convert_units:
            self.convert_pos_from_native(self.ts.dimensions[:3])

    def close(self):
        """close reader"""
        self._tng.close()

    @property
    def n_frames(self):
        """number of frames in trajectory"""
        return len(self._tng)

    def _reopen(self):
        """reopen trajectory"""
        self.ts.frame = 0
        self._frame = -1
        self._tng.close()
        self._tng.open(self.filename.encode('utf-8'), 'r')

    def _read_frame(self, i):
        """read frame i"""
        self._frame = i - 1
        self._tng.seek(i)
        return self._read_next_timestep()

    def _read_next_timestep(self, ts=None):
        """copy next frame into timestep"""
        if self._frame == self.n_frames - 1:
            raise IOError(errno.EIO, 'trying to go over trajectory limit')
        if ts is None:
            ts = self.ts
        frame = self._tng.read()
        self._frame += 1
        self._frame_to_ts(frame, ts)
        return ts

    def _frame_to_ts(self, frame, ts):
        """convert a trr-frame to a mda TimeStep"""
        ts.time = frame.time
        ts.frame = self._frame
        ts.data['step'] = frame.step

        ts.positions = frame.xyz
        ts.dimensions = triclinic_box(*frame.box)

        if self.convert_units:
            self.convert_pos_from_native(ts.dimensions[:3])
            self.convert_pos_from_native(ts.positions)

        return ts
