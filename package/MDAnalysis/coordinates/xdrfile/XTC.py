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

"""
Gromacs XTC trajectory I/O --- :mod:`MDAnalysis.coordinates.xdrfile.XTC`
========================================================================

Classes for reading and writing of `Gromacs XTC trajectories`_
together with supporting code.

.. note:: Users should access classes from :mod:`MDAnalysis.coordinates.XTC`.

.. _Gromacs XTC trajectories: http://www.gromacs.org/Documentation/File_Formats/.xtc_File
.. _Gromacs: http://www.gromacs.org


.. SeeAlso:: :mod:`MDAnalysis.coordinates.xdrfile.libxdrfile2` for low-level
   bindings to the Gromacs trajectory file formats

Classes
-------

.. autoclass:: Timestep
   :members:
   :inherited-members:
.. autoclass:: XTCReader
   :members:
   :inherited-members:
.. autoclass:: XTCWriter
   :members:
   :inherited-members:
"""

import numpy as np
import errno

from . import statno
from . import core
from . import libxdrfile2
from .core import Timestep


class XTCWriter(core.TrjWriter):
    """Write a Gromacs_ XTC trajectory.

    .. _Gromacs: http://www.gromacs.org
    """
    format = "XTC"


class XTCReader(core.TrjReader):
    """Read Gromacs_ XTC trajectory.

    .. _Gromacs: http://www.gromacs.org

    .. versionchanged:: 0.11.0
       Timestep attributes status and prec are now stored in the TS.data dictionary
    """
    format = "XTC"
    _Timestep = Timestep
    _Writer = XTCWriter

    def _allocate_sub(self, DIM):
        self._pos_buf = np.zeros((self._trr_n_atoms, DIM), dtype=np.float32, order='C')
        self._velocities_buf = None
        self._forces_buf = None

    def _read_trj_natoms(self, filename):
        return libxdrfile2.read_xtc_natoms(filename)

    def _read_trj_n_frames(self, filename):
        self._n_frames, self._offsets = libxdrfile2.read_xtc_n_frames(filename)
        self._store_offsets()

    def _read_next_timestep(self, ts=None):
        if ts is None:
            ts = self.ts

        if self.xdrfile is None:
            self.open_trajectory()

        if self._sub is None:
            ts.data['status'], ts._frame, ts.time, ts.data['prec'] = libxdrfile2.read_xtc(self.xdrfile, ts._unitcell, ts._pos)
        else:
            ts.data['status'], ts._frame, ts.time, ts.data['prec'] = libxdrfile2.read_xtc(self.xdrfile, ts._unitcell, self._pos_buf)
            ts._pos[:] = self._pos_buf[self._sub]

        if ts.data['status'] == libxdrfile2.exdrENDOFFILE:
            raise IOError(errno.EIO, "End of file reached for {0} file".format(self.format),
                          self.filename)
        elif not ts.data['status'] == libxdrfile2.exdrOK:
            raise IOError(errno.EBADF, ("Problem with {0} file, status {1}"
                                        "".format((self.format, statno.ERRORCODE[ts.data['status']]), self.filename)))

        if self.convert_units:
            self.convert_pos_from_native(ts._pos)  # in-place !
            self.convert_pos_from_native(ts._unitcell)  # in-place ! (note: xtc contain unit vecs!)
            ts.time = self.convert_time_from_native(ts.time)  # in-place does not work with scalars

        ts.frame += 1
        return ts
