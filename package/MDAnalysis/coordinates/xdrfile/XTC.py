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

import numpy
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
    """
    format = "XTC"
    _Timestep = Timestep
    _Writer = XTCWriter

    def _allocate_sub(self, DIM):
        self._pos_buf = numpy.zeros((self._trr_numatoms, DIM), dtype=numpy.float32, order='C')
        self._velocities_buf = None
        self._forces_buf = None

    def _read_trj_natoms(self, filename):
        return libxdrfile2.read_xtc_natoms(filename)

    def _read_trj_numframes(self, filename):
        self._numframes, self._offsets = libxdrfile2.read_xtc_numframes(filename)
        self._store_offsets()

        return

    def _read_next_timestep(self, ts=None):
        if ts is None:
            ts = self.ts

        if self.xdrfile is None:
            self.open_trajectory()

        if self._sub is None:
            ts.status, ts.step, ts.time, ts.prec = libxdrfile2.read_xtc(self.xdrfile, ts._unitcell, ts._pos)
        else:
            ts.status, ts.step, ts.time, ts.prec = libxdrfile2.read_xtc(self.xdrfile, ts._unitcell, self._pos_buf)
            ts._pos[:] = self._pos_buf[self._sub]

        if ts.status == libxdrfile2.exdrENDOFFILE:
            raise IOError(errno.EIO, "End of file reached for %s file" % self.format,
                          self.filename)
        elif not ts.status == libxdrfile2.exdrOK:
            raise IOError(errno.EBADF, "Problem with %s file, status %s" %
                                       (self.format, statno.errorcode[ts.status]), self.filename)

        if self.convert_units:
            self.convert_pos_from_native(ts._pos)  # in-place !
            self.convert_pos_from_native(ts._unitcell)  # in-place ! (note: xtc contain unit vecs!)
            ts.time = self.convert_time_from_native(ts.time)  # in-place does not work with scalars

        ts.frame += 1
        return ts
