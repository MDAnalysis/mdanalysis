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
from six import string_types
from distutils.version import LooseVersion
import warnings

from . import base, core

try:
    import chemfiles
except ImportError:
    HAS_CHEMFILES = False
else:
    HAS_CHEMFILES = True


MIN_CHEMFILES_VERSION = LooseVersion("0.9")
MAX_CHEMFILES_VERSION = LooseVersion("0.10")


def check_chemfiles_version():
    """Check an appropriate Chemfiles is available

    Returns True if a usable chemfiles version is available

    .. versionadded:: 1.0.0
    """
    if not HAS_CHEMFILES:
        warnings.warn(
            "No Chemfiles package found.  "
            "Try installing with 'pip install chemfiles'"
        )
        return False
    version = LooseVersion(chemfiles.__version__)
    wrong = version < MIN_CHEMFILES_VERSION or version >= MAX_CHEMFILES_VERSION
    if wrong:
        warnings.warn(
            "unsupported Chemfiles version {}, we need a version >{} and <{}"
            .format(version, MIN_CHEMFILES_VERSION, MAX_CHEMFILES_VERSION)
        )
    return not wrong


class ChemfilesReader(base.ReaderBase):
    """Coordinate reader using chemfiles.

    The file format to used is guessed based on the file extension. If no
    matching format is found, a ``ChemfilesError`` is raised. It is also
    possible to manually specify the format to use for a given file.

    .. versionadded:: 1.0.0
    """
    format = 'chemfiles'
    units = {'time': 'fs', 'length': 'Angstrom'}

    def __init__(self, filename, chemfiles_format="", **kwargs):
        """
        Parameters
        ----------
        filename : chemfiles.Trajectory or str
            the chemfiles object to read or filename to read
        chemfiles_format : str (optional)
            if *filename* was a string, use the given format name instead of
            guessing from the extension. The `list of supported formats
            <formats>`_ and the associated names is available in chemfiles
            documentation.
        **kwargs : dict
            General reader arguments.


        .. _formats: http://chemfiles.org/chemfiles/latest/formats.html
        """
        if not check_chemfiles_version():
            raise RuntimeError("Please install Chemfiles > {}"
                               "".format(MIN_CHEMFILES_VERSION))
        super(ChemfilesReader, self).__init__(filename, **kwargs)
        self._format = chemfiles_format
        self._cached_n_atoms = None
        self._open()
        self.ts = self._Timestep(self.n_atoms, **self._ts_kwargs)
        self.next()

    @staticmethod
    def _format_hint(thing):
        """Can this Reader read *thing*"""
        # nb, filename strings can still get passed through if
        # format='CHEMFILES' is used
        return HAS_CHEMFILES and isinstance(thing, chemfiles.Trajectory)

    def _open(self):
        self._closed = False
        self._step = 0
        self._frame = -1
        if isinstance(self.filename, chemfiles.Trajectory):
            self._file = self.filename
        else:
            self._file = chemfiles.Trajectory(self.filename, 'r', self._format)

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
        ts.frame = frame.step
        self._step += 1
        return ts

    def _frame_to_ts(self, frame, ts):
        """convert a chemfiles frame to a :class:`TimeStep`"""
        if len(frame.atoms) != self.n_atoms:
            raise IOError(
                "The number of atom changed in the trajectory. "
                "This is not supported by MDAnalysis."
            )

        ts.dimensions[:] = frame.cell.lengths + frame.cell.angles
        ts.positions[:] = frame.positions[:]
        if frame.has_velocities():
            ts.has_velocities = True
            ts.velocities[:] = frame.velocities[:]

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

    .. versionadded:: 1.0.0
    """
    format = 'chemfiles'
    multiframe = True

    # chemfiles mostly[1] uses these units for the in-memory representation,
    # and convert into the file format units when writing.
    #
    # [1] mostly since some format don't have a specified unit
    # (XYZ for example), so then chemfiles just assume they are in A and fs.
    units = {'time': 'fs', 'length': 'Angstrom'}

    def __init__(self, filename, n_atoms=0, mode="w", chemfiles_format="", topology=None, **kwargs):
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
        topology : Universe or AtomGroup (optional)
            use the topology from this :class:`~MDAnalysis.core.groups.AtomGroup`
            or :class:`~MDAnalysis.core.universe.Universe` to write all the
            timesteps to the file
        **kwargs : dict
            General writer arguments.


        .. _formats: http://chemfiles.org/chemfiles/latest/formats.html
        """
        if not check_chemfiles_version():
            raise RuntimeError("Please install Chemfiles > {}"
                               "".format(MIN_CHEMFILES_VERSION))
        self.filename = filename
        self.n_atoms = n_atoms
        if mode != "a" and mode != "w":
            raise IOError("Expected 'a' or 'w' as mode in ChemfilesWriter")
        self._file = chemfiles.Trajectory(self.filename, mode, chemfiles_format)
        self._closed = False
        if topology is not None:
            if isinstance(topology, string_types):
                self._file.set_topology(topology)
            else:
                topology = self._topology_to_chemfiles(topology, n_atoms)
                self._file.set_topology(topology)

    def close(self):
        """Close the trajectory file and finalize the writing"""
        if hasattr(self, "_closed") and not self._closed:
            self._file.close()
            self._closed = True

    def _write_next_frame(self, obj):
        """Write information associated with ``obj`` at current frame into trajectory.

        Topology for the output is taken from the `obj` or default to the value
        of the `topology` keyword supplied to the :class:`ChemfilesWriter`
        constructor.

        If `obj` contains velocities, and the underlying format supports it, the
        velocities are writen to the file. Writing forces is unsupported at the
        moment.

        Parameters
        ----------
        obj : AtomGroup or Universe
            The :class:`~MDAnalysis.core.groups.AtomGroup` or
            :class:`~MDAnalysis.core.universe.Universe` to write.
        """
        # TODO 2.0: Remove timestep logic
        if hasattr(obj, "atoms"):
            if hasattr(obj, 'universe'):
                # For AtomGroup and children (Residue, ResidueGroup, Segment)
                ts_full = obj.universe.trajectory.ts
                if ts_full.n_atoms == obj.atoms.n_atoms:
                    ts = ts_full
                else:
                    # Only populate a time step with the selected atoms.
                    ts = ts_full.copy_slice(atoms.indices)
            elif hasattr(obj, 'trajectory'):
                # For Universe only --- get everything
                ts = obj.trajectory.ts
        else:
            if isinstance(obj, base.Timestep):
                ts = obj
                topology = None
            else:
                raise TypeError("No Timestep found in obj argument")

        frame = self._timestep_to_chemfiles(ts)
        frame.topology = self._topology_to_chemfiles(obj, len(frame.atoms))

        self._file.write(frame)

    def _timestep_to_chemfiles(self, ts):
        """
        Convert a Timestep to a chemfiles Frame
        """
        # TODO: CONVERTERS?
        frame = chemfiles.Frame()
        frame.resize(ts.n_atoms)
        if ts.has_positions:
            frame.positions[:] = ts.positions[:]
        if ts.has_velocities:
            frame.add_velocities()
            frame.velocities[:] = ts.velocities[:]
        frame.cell = chemfiles.UnitCell(*ts.dimensions)
        return frame

    def _topology_to_chemfiles(self, obj, n_atoms):
        """
        Convert an AtomGroup or an Universe to a chemfiles Topology
        """
        topology = chemfiles.Topology()
        if not hasattr(obj, "atoms"):
            # use an empty topology
            topology.resize(n_atoms)
            return topology

        # (1) add all atoms to the topology
        residues = {}
        for atom in obj.atoms:
            name = getattr(atom, 'name', "")
            type = getattr(atom, 'type', name)
            chemfiles_atom = chemfiles.Atom(name, type)

            if hasattr(atom, 'altLoc'):
                chemfiles_atom["altloc"] = str(atom.altLoc)

            if hasattr(atom, 'segid'):
                chemfiles_atom["segid"] = str(atom.segid)

            if hasattr(atom, 'segindex'):
                chemfiles_atom["segindex"] = int(atom.segindex)

            if hasattr(atom, 'resid'):
                resname = getattr(atom, 'resname', "")
                if atom.resid not in residues.keys():
                    residues[atom.resid] = chemfiles.Residue(resname, atom.resid)
                residue = residues[atom.resid]

                atom_idx = len(topology.atoms)
                residue.atoms.append(atom_idx)

                if hasattr(atom, "record_type"):
                    # set corresponding chemfiles residue property
                    if atom.record_type == "ATOM":
                        residue["is_standard_pdb"] = True
                    else:
                        residue["is_standard_pdb"] = False

            topology.atoms.append(chemfiles_atom)

        # (2) add residues to the topology
        for residue in residues.values():
            topology.residues.append(residue)

        # (3) add bonds to the topology
        for bond in getattr(obj, "bonds", []):
            topology.add_bond(bond.atoms[0].ix, bond.atoms[1].ix)

        return topology
