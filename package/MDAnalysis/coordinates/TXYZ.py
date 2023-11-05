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
"""
TXYZ file format --- :mod:`MDAnalysis.coordinates.TXYZ`
=======================================================

Coordinate reader for Tinker_ xyz files .txyz and trajectory .arc files.
Differences between Tinker format_ and normal xyz files:

- there is only one header line containing both the number of atoms and a comment
- column 1 contains atom numbers (starting from 1)
- column 6 contains atoms types
- the following columns indicate connectivity (atoms to which that particular atom is
  bonded, according to numbering in column 1)

.. _format: http://chembytes.wikidot.com/tnk-tut00#toc2
.. _Tinker: https://dasher.wustl.edu/tinker/


Classes
-------

.. autoclass:: TXYZReader
   :members:
   :inherited-members:

"""
import numpy as np
import os
import errno

from ..lib import util
from . import base
from ..lib.util import openany, cached, store_init_arguments
from .timestep import Timestep

class TXYZReader(base.ReaderBase):
    """Reads from a TXYZ file"""


    format = ['TXYZ', 'ARC']
    # these are assumed!
    units = {'time': 'ps', 'length': 'Angstrom'}
    _Timestep = Timestep

    @store_init_arguments
    def __init__(self, filename, **kwargs):
        super(TXYZReader, self).__init__(filename, **kwargs)

        # the filename has been parsed to be either be foo.xyz or foo.xyz.bz2 by
        # coordinates::core.py so the last file extension will tell us if it is
        # bzipped or not
        root, ext = os.path.splitext(self.filename)
        self.xyzfile = util.anyopen(self.filename)
        self._cache = dict()
        # Check if file has box information saved
        with util.openany(self.filename) as inp:
            inp.readline()
            line = inp.readline()
            # If second line has float at second position, we have box info
            try:
                float(line.split()[1])
            except ValueError:
                self.periodic = False
            else:
                self.periodic = True
        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)

        self._read_next_timestep()

    @property
    @cached('n_atoms')
    def n_atoms(self):
        """number of atoms in a frame"""
        with util.anyopen(self.filename) as f:
            n = f.readline().split()[0]
        # need to check type of n
        return int(n)

    @property
    @cached('n_frames')
    def n_frames(self):
        try:
            return self._read_xyz_n_frames()
        except IOError:
            return 0

    def _read_xyz_n_frames(self):
        # the number of lines in the XYZ file will be 1 greater than the
        # number of atoms
        linesPerFrame = self.n_atoms + 1
        if self.periodic:
            linesPerFrame += 1
        counter = 0
        offsets = []

        with util.anyopen(self.filename) as f:
            line = True
            while line:
                if not counter % linesPerFrame:
                    offsets.append(f.tell())
                line = f.readline()
                counter += 1

        # need to check this is an integer!
        n_frames = int(counter / linesPerFrame)
        self._offsets = offsets
        return n_frames

    def _read_frame(self, frame):
        self.xyzfile.seek(self._offsets[frame])
        self.ts.frame = frame - 1  # gets +1'd in next
        return self._read_next_timestep()

    def _read_next_timestep(self, ts=None):
        if ts:
            warnings.warn("ts argument to _read_next_timestep is deprecated as of 2.7.0 and will be removed in 3.0.0, see #3928")

        # check that the timestep object exists
        if ts is None:
            ts = self.ts

        f = self.xyzfile

        try:
            # we assume that there is only one header line per frame
            f.readline()
            if self.periodic:
                ts.dimensions = f.readline().split() 
            # convert all entries at the end once for optimal speed
            tmp_buf = []
            for i in range(self.n_atoms):
                tmp_buf.append(f.readline().split()[2:5])
            ts.positions = tmp_buf
            ts.frame += 1
            return ts
        except (ValueError, IndexError) as err:
            raise EOFError(err) from None

    def _reopen(self):
        self.close()
        self.open_trajectory()

    def open_trajectory(self):
        if self.xyzfile is not None:
            raise IOError(
                errno.EALREADY, 'TXYZ file already opened', self.filename)

        self.xyzfile = util.anyopen(self.filename)

        # reset ts
        ts = self.ts
        ts.frame = -1

        return self.xyzfile

    def close(self):
        """Close arc trajectory file if it was open."""
        if self.xyzfile is None:
            return
        self.xyzfile.close()
        self.xyzfile = None
