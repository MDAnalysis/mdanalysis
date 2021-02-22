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
Reading trajectories with *chemfiles* --- :mod:`MDAnalysis.coordinates.chemfiles`
=================================================================================

MDAnalysis interoperates with the `chemfiles`_ library. The *chemfiles* C++ library 
supports an expanding set of file formats, some of which are not natively supported by
MDAnalysis. Using the *CHEMFILES* reader you can use  `chemfiles`_ for the low-level 
file reading. Check the list of `chemfile-supported file formats <formats>`_.

.. _chemfiles: https://chemfiles.org
.. _formats: https://chemfiles.org/chemfiles/0.9.3/formats.html#list-of-supported-formats
.. NOTE: MDAnalysis currently restricts chemfiles to 0.9 <= version < 0.10. Update the link
..       above to the latest documentation once this restriction is lifted.
..       https://chemfiles.org/chemfiles/latest/formats.html#list-of-supported-formats

Using the CHEMFILES reader
--------------------------

When reading, set the ``format="CHEMFILES"`` keyword argument and I/O is delegated to 
`chemfiles`_. For example::

   >>> import MDAnalysis as mda
   >>> from MDAnalysis.tests import datafiles as data
   >>> u = mda.Universe(data.TPR, data.TRR, format="CHEMFILES")
   >>> print(u.trajectory)
   <ChemfilesReader ~/anaconda3/envs/mda3/lib/python3.8/site-packages/MDAnalysisTests/data/adk_oplsaa.trr with 10 frames of 47681 atoms>

You can then use the :class:`~MDAnalysis.core.universe.Universe` as usual while chemfiles
is handling the I/O transparently in the background.

`chemfiles`_ can also *write* a number of formats for which there are no Writers in
MDAnalysis. For example, to write a mol2 file::

   >>> u = mda.Universe(data.mol2_ligand)
   >>> with mda.Writer("ligand.mol2", format="CHEMFILES") as W:
   ...     W.write(u.atoms)




Classes
-------

Classes to read and write files using the `chemfiles`_ library. This library
provides C++ implementation of multiple formats readers and writers.

.. autoclass:: ChemfilesReader

.. autoclass:: ChemfilesWriter

.. autoclass:: ChemfilesPicklable

Helper functions
----------------

.. autodata:: MIN_CHEMFILES_VERSION
.. autodata:: MAX_CHEMFILES_VERSION
.. autofunction:: check_chemfiles_version

"""
from distutils.version import LooseVersion
import warnings

from . import base, core

try:
    import chemfiles
except ImportError:
    HAS_CHEMFILES = False

    # Allow building documentation even if chemfiles is not installed
    import types

    class MockTrajectory:
        pass
    chemfiles = types.ModuleType("chemfiles")
    chemfiles.Trajectory = MockTrajectory
else:
    HAS_CHEMFILES = True


#: Lowest version of chemfiles that is supported
#: by MDAnalysis.
MIN_CHEMFILES_VERSION = LooseVersion("0.9")
#: Lowest version of chemfiles that is *not supported*
#: by MDAnalysis.
MAX_CHEMFILES_VERSION = LooseVersion("0.10")


def check_chemfiles_version():
    """Check if an appropriate *chemfiles* is available

    Returns ``True`` if a usable chemfiles version is available,
    with :data:`MIN_CHEMFILES_VERSION` <= version < 
    :data:`MAX_CHEMFILES_VERSION`

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
            <formats>`_ and the associated names is available in the chemfiles
            documentation.
        **kwargs : dict
            General reader arguments.

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
            self._file = ChemfilesPicklable(self.filename, 'r', self._format)

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
            if isinstance(topology, str):
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
        """Write information associated with ``obj`` at current frame into
        trajectory.

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


        .. versionchanged:: 2.0.0
           Deprecated support for Timestep argument has now been removed.
           Use AtomGroup or Universe as an input instead.
        """
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
            errmsg = "Input obj is neither an AtomGroup or Universe"
            raise TypeError(errmsg) from None

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


class ChemfilesPicklable(chemfiles.Trajectory):
    """Chemfiles file object (read-only) that can be pickled.

    This class provides a file-like object (as returned by
    :class:`chemfiles.Trajectory`) that,
    unlike standard Python file objects,
    can be pickled. Only read mode is supported.

    When the file is pickled, path, mode, and format of the file handle
    are saved. On unpickling, the file is opened by path with mode,
    and saved format.
    This means that for a successful unpickle, the original file still has
    to be accessible with its filename.

    Note
    ----
    Can only be used with reading ('r') mode.
    Upon pickling, the current frame is reset. `universe.trajectory[i]` has
    to be used to return to its original frame.

    Parameters
    ----------
    filename : str
        a filename given a text or byte string.
    mode : 'r' , optional
        only 'r' can be used for pickling.
    format : '', optional
        guessed from the file extension if empty.

    Example
    -------
    ::

        f = ChemfilesPicklable(XYZ, 'r', '')
        print(f.read())
        f.close()

    can also be used as context manager::

        with ChemfilesPicklable(XYZ) as f:
            print(f.read())

    See Also
    ---------
    :class:`MDAnalysis.lib.picklable_file_io.FileIOPicklable`
    :class:`MDAnalysis.lib.picklable_file_io.BufferIOPicklable`
    :class:`MDAnalysis.lib.picklable_file_io.TextIOPicklable`
    :class:`MDAnalysis.lib.picklable_file_io.GzipPicklable`
    :class:`MDAnalysis.lib.picklable_file_io.BZ2Picklable`


    .. versionadded:: 2.0.0
    """
    def __init__(self, path, mode="r", format=""):
        if mode != 'r':
            raise ValueError("Only read mode ('r') "
                             "files can be pickled.")
        super().__init__(path=path,
                         mode=mode,
                         format=format)

    def __getstate__(self):
        return self.path, self._Trajectory__mode, self._Trajectory__format

    def __setstate__(self, args):
        self.__init__(*args)
