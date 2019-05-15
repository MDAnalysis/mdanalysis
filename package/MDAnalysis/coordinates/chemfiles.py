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
"""
Reading trajectory with Chemfiles --- :mod:`MDAnalysis.coordinates.chemfiles`
=============================================================================

Classes to read and write files using the `chemfiles`_ library. This library
provides C++ implementation of multiple formats readers and writers, the full
list if available `here <formats>`_.

This module only contains the coordinate reader and writer, which means that the
topology read by chemfiles is ignored when reading a file, and no topology is
written either. This can be surprising when using formats where topology and
coordinates are not separated, such as XYZ or PDB.

.. _chemfiles: https://chemfiles.org
.. _formats: http://chemfiles.org/chemfiles/latest/formats.html

Classes
-------

.. autoclass:: ChemfilesReader

.. autoclass:: ChemfilesWriter

"""
from __future__ import absolute_import, division
from __future__ import print_function, unicode_literals

import numpy as np

from . import base, core

import chemfiles
from chemfiles import Trajectory, Frame, Atom, UnitCell, Residue


assert chemfiles.__version__.startswith("0.9.")


class ChemfilesReader(base.ReaderBase):
    """
    Coordinate reader using chemfiles.

    The file format to used is guessed based on the file extension. If no
    matching format is found, a ``ChemfilesError`` is raised. It is also
    possible to manually specify the format to use for a given file.
    """
    format = 'chemfiles'
    units = {'time': 'fs', 'length': 'Angstrom'}

    def __init__(self, filename, chemfiles_format="", **kwargs):
        """
        Parameters
        ----------
        filename : str
            trajectory filename
        chemfiles_format : str (optional)
            use the given format name instead of guessing from the extension.
            The `list of supported formats <formats>`_ and the associated names
            is available in chemfiles documentation.
        **kwargs : dict
            General reader arguments.


        .. _formats: http://chemfiles.org/chemfiles/latest/formats.html
        """
        super(ChemfilesReader, self).__init__(filename, **kwargs)
        self._format = chemfiles_format
        self._cached_n_atoms = None
        self._open()
        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
        self.next()
        print(kwargs)

    def _open(self):
        self._closed = False
        self._step = 0
        self._frame = -1
        # Open file last to ensure that all other attributes are set
        # in case of exception
        self._file = Trajectory(self.filename, 'r', self._format)

    def close(self):
        """close reader"""
        if not self._closed:
            self._file.close()
            self._closed = True

    @property
    def n_frames(self):
        """number of frames in trajectory"""
        return self._file.nsteps

    @property
    def n_atoms(self):
        """number of atoms in the first frame of the trajectory"""
        if self._cached_n_atoms is None:
            self._cached_n_atoms = len(self._file.read_step(0).atoms)
        return self._cached_n_atoms

    def _reopen(self):
        """reopen trajectory"""
        self.close()
        self._open()

    def _read_frame(self, i):
        """read frame i"""
        self._step = i
        return self._read_next_timestep()

    def _read_next_timestep(self, ts=None):
        """copy next frame into timestep"""
        if self._step >= self.n_frames:
            raise IOError('trying to go over trajectory limit')
        if ts is None:
            ts = self.ts
        self.ts = ts
        frame = self._file.read_step(self._step)
        self._frame_to_ts(frame, ts)
        ts.frame = self._step
        self._step += 1
        return ts

    def _frame_to_ts(self, frame, ts):
        """convert a chemfiles frame to a :class:`TimeStep`"""
        if len(frame.atoms) != self.n_atoms:
            raise IOError(
                "The number of atom changed in the trajectory. This is not "
                "supported by MDAnalysis."
            )

        ts.dimensions[:] = frame.cell.lengths + frame.cell.angles
        ts.positions[:] = frame.positions[:]
        if frame.has_velocities():
            ts.has_velocities = True
            ts.velocities[:] = frame.velocities[:]

    @property
    def dimensions(self):
        """unitcell dimensions (*A*, *B*, *C*, *alpha*, *beta*, *gamma*)
        """
        return self.ts.dimensions

    def Writer(self, filename, n_atoms=None, **kwargs):
        """Return writer for trajectory format"""
        if n_atoms is None:
            n_atoms = self.n_atoms
        return ChemfilesWriter(filename, n_atoms, **kwargs)


class ChemfilesWriter(base.WriterBase):
    """
    Coordinate writer using chemfiles.

    The file format to used is guessed based on the file extension. If no
    matching format is found, a ``ChemfilesError`` is raised. It is also
    possible to manually specify the format to use for a given file.

    Chemfiles support writting to files with varying number of atoms if the
    underlying format support it. This is the case of most of text-based
    formats.
    """

    format = 'chemfiles'
    multiframe = True
    units = {'time': 'fs', 'length': 'Angstrom'}

    def __init__(self, filename, n_atoms=0, mode="w", chemfiles_format="", **kwargs):
        """
        Parameters
        ----------
        filename: str
            filename of trajectory file.
        n_atoms: int
            number of atoms in the trajectory to write. This value is not
            used and can vary during trajectory, if the underlying format
            support it
        mode : str (optional)
            file opening mode: use 'a' to append to an existing file or 'w' to
            create a new file
        chemfiles_format : str (optional)
            use the given format name instead of guessing from the extension.
            The `list of supported formats <formats>`_ and the associated names
            is available in chemfiles documentation.
        **kwargs : dict
            General writer arguments.


        .. _formats: http://chemfiles.org/chemfiles/latest/formats.html
        """
        self.filename = filename
        self.n_atoms = n_atoms
        if mode != "a" and mode != "w":
            raise IOError("Expected 'a' or 'w' as mode in ChemfilesWriter")
        self._file = Trajectory(self.filename, mode, chemfiles_format)
        self._closed = False

    def close(self):
        """Close the trajectory file and finalize the writing"""
        if hasattr(self, "_closed") and not self._closed:
            self._file.close()
            self._closed = True

    def write_next_timestep(self, ts):
        """Write timestep object into trajectory.

        Parameters
        ----------
        ts: TimeStep
        """
        frame = self._ts_to_frame(ts)
        self._file.write(frame)

    def _ts_to_frame(self, ts):
        """
        Convert a Timestep to a chemfiles Frame
        """
        frame = Frame()
        frame.resize(ts.n_atoms)
        if ts.has_positions:
            frame.positions[:] = ts.positions[:]
        if ts.has_velocities:
            frame.add_velocities()
            frame.velocities[:] = ts.velocities[:]
        frame.cell = UnitCell(*ts.dimensions)
        return frame
