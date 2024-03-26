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
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""\
=========================================================
Core object: Universe --- :mod:`MDAnalysis.core.universe`
=========================================================

The :class:`~MDAnalysis.core.universe.Universe` class ties a topology
and a trajectory together. Almost all code in MDAnalysis starts with a
``Universe``.

Normally, a ``Universe`` is created from files::

  import MDAnalysis as mda
  u = mda.Universe("topology.psf", "trajectory.dcd")

In order to construct new simulation system it is also convenient to
construct a ``Universe`` from existing
:class:`~MDAnalysis.core.group.AtomGroup` instances with the
:func:`Merge` function.


Classes
=======

.. autoclass:: Universe
   :members:

Functions
=========

.. autofunction:: Merge

"""
import errno
import numpy as np
import logging
import copy
import warnings
import contextlib
import collections

import MDAnalysis
import sys

# When used in an MPI environment with Infiniband, importing MDAnalysis may
# trigger an MPI warning because importing the uuid module triggers a call to
# os.fork(). This happens if MPI_Init() has been called prior to importing
# MDAnalysis. The problem is actually caused by the uuid module and not by
# MDAnalysis itself. Python 3.7 fixes the problem.
import os
import uuid

from .. import _TOPOLOGY_ATTRS, _PARSERS
from ..exceptions import NoDataError
from ..lib import util
from ..lib.log import ProgressBar
from ..lib.util import cached, NamedStream, isstream
from ..lib.mdamath import find_fragments
from . import groups
from ._get_readers import get_reader_for, get_parser_for
from .groups import (ComponentBase, GroupBase,
                     Atom, Residue, Segment,
                     AtomGroup, ResidueGroup, SegmentGroup)
from .topology import Topology
from .topologyattrs import AtomAttr, ResidueAttr, SegmentAttr, BFACTOR_WARNING
from .topologyobjects import TopologyObject
from ..guesser.base import get_guesser

logger = logging.getLogger("MDAnalysis.core.universe")



def _check_file_like(topology):
    if isstream(topology):
        if hasattr(topology, 'name'):
            _name = topology.name
        else:
            _name = None
        return NamedStream(topology, _name)
    return topology

def _topology_from_file_like(topology_file, topology_format=None,
                             **kwargs):
    parser = get_parser_for(topology_file, format=topology_format)

    try:
        with parser(topology_file) as p:
            topology = p.parse(**kwargs)
    except (IOError, OSError) as err:
        # There are 2 kinds of errors that might be raised here:
        # one because the file isn't present
        # or the permissions are bad, second when the parser fails
        if (err.errno is not None and
            errno.errorcode[err.errno] in ['ENOENT', 'EACCES']):
            # Runs if the error is propagated due to no permission / file not found
            raise sys.exc_info()[1] from err
        else:
            # Runs when the parser fails
            raise IOError("Failed to load from the topology file {0}"
                            " with parser {1}.\n"
                            "Error: {2}".format(topology_file, parser, err))
    except (ValueError, NotImplementedError) as err:
        raise ValueError(
            "Failed to construct topology from file {0}"
            " with parser {1}.\n"
            "Error: {2}".format(topology_file, parser, err))
    return topology


def _resolve_formats(*coordinates, format=None, topology_format=None):
    if not coordinates:
        if format is None:
            format = topology_format
        elif topology_format is None:
            topology_format = format
    return format, topology_format


def _resolve_coordinates(filename, *coordinates, format=None,
                         all_coordinates=False):
    if all_coordinates or not coordinates and filename is not None:
        try:
            get_reader_for(filename, format=format)
        except (ValueError, TypeError):
            warnings.warn('No coordinate reader found for {}. Skipping '
                            'this file.'.format(filename))
        else:
            coordinates = (filename,) + coordinates
    return coordinates

def _generate_from_topology(universe):
    # generate Universe version of each class
    # AG, RG, SG, A, R, S
    universe._class_bases, universe._classes = groups.make_classes()

    # Put Group level stuff from topology into class
    for attr in universe._topology.attrs:
        universe._process_attr(attr)

    # Generate atoms, residues and segments.
    # These are the first such groups generated for this universe, so
    #  there are no cached merged classes yet. Otherwise those could be
    #  used directly to get a (very) small speedup. (Only really pays off
    #  the readability loss if instantiating millions of AtomGroups at
    #  once.)
    universe.atoms = AtomGroup(np.arange(universe._topology.n_atoms), universe)

    universe.residues = ResidueGroup(
            np.arange(universe._topology.n_residues), universe)

    universe.segments = SegmentGroup(
            np.arange(universe._topology.n_segments), universe)


class Universe(object):
    """The MDAnalysis Universe contains all the information describing the system.

    The system always requires a *topology file* --- in the simplest case just
    a list of atoms. This can be a CHARMM/NAMD PSF file or a simple coordinate
    file with atom informations such as XYZ, PDB, GROMACS GRO or TPR, or CHARMM
    CRD. See :ref:`Supported topology formats` for what kind of topologies can
    be read.

    A *trajectory file* provides coordinates; the coordinates have to be
    ordered in the same way as the list of atoms in the topology. A trajectory
    can be a single frame such as a PDB, CRD, or GRO file, or it can be a MD
    trajectory (in CHARMM/NAMD/LAMMPS DCD, GROMACS XTC/TRR, AMBER nc, generic
    XYZ format, ...).  See :ref:`Supported coordinate formats` for what can be
    read as a "trajectory".

    As a special case, when the topology is a file that contains atom
    information *and* coordinates (such as XYZ, PDB, GRO or CRD, see
    :ref:`Supported coordinate formats`) then the coordinates are immediately
    loaded from the "topology" file unless a trajectory is supplied.

    Parameters
    ----------
    topology: str, stream, Topology, numpy.ndarray, None
        A CHARMM/XPLOR PSF topology file, PDB file or Gromacs GRO file; used to
        define the list of atoms. If the file includes bond information,
        partial charges, atom masses, ... then these data will be available to
        MDAnalysis. Alternatively, an existing
        :class:`MDAnalysis.core.topology.Topology` instance may be given,
        numpy coordinates, or ``None`` for an empty universe.
    coordinates: str, stream, list of str, list of stream (optional)
        Coordinates can be provided as files of
        a single frame (eg a PDB, CRD, or GRO file); a list of single
        frames; or a trajectory file (in CHARMM/NAMD/LAMMPS DCD, Gromacs
        XTC/TRR, or generic XYZ format). The coordinates must be
        ordered in the same way as the list of atoms in the topology.
        See :ref:`Supported coordinate formats` for what can be read
        as coordinates. Alternatively, streams can be given.
    topology_format: str, ``None``, default ``None``
        Provide the file format of the topology file; ``None`` guesses it from
        the file extension. Can also pass a subclass of
        :class:`MDAnalysis.topology.base.TopologyReaderBase` to define a custom
        reader to be used on the topology file.
    format: str, ``None``, default ``None``
        Provide the file format of the coordinate or trajectory file; ``None``
        guesses it from the file extension. Note that this keyword has no
        effect if a list of file names is supplied because the "chained" reader
        has to guess the file format for each individual list member.
        Can also pass a subclass of :class:`MDAnalysis.coordinates.base.ProtoReader`
        to define a custom reader to be used on the trajectory file.
    all_coordinates: bool, default ``False``
        If set to ``True`` specifies that if more than one filename is passed
        they are all to be used, if possible, as coordinate files (employing a
        :class:`MDAnalysis.coordinates.chain.ChainReader`). The
        default behavior is to take the first file as a topology and the
        remaining as coordinates. The first argument will always always be used
        to infer a topology regardless of *all_coordinates*.
    guess_bonds: bool, default ``False``
        Once Universe has been loaded, attempt to guess the connectivity
        between atoms.  This will populate the .bonds, .angles, and .dihedrals
        attributes of the Universe.
    vdwradii: dict, ``None``, default ``None``
        For use with *guess_bonds*. Supply a dict giving a vdwradii for each
        atom type which are used in guessing bonds.
    context: str or Guesser, default ``default``
        Type of the Guesser to be used in guessing TopologyAttrs
    to_guess: list[str], (optional), default ``['types', 'masses']``
              (in future versions types and masses will be removed)
        TopologyAttrs to be guessed. These TopologyAttrs will be wholly guessed
        if they don't exist in the universe. If they already exist in the Universe,
        only empty or missing values will be guessed.
    force_guess: list[str], (optional)
        TopologyAttrs in this list will be force guessed. If the TopologyAttr
        does not already exist in the Universe, this has no effect. If the TopologyAttr
        does already exist, all values will be overwritten by guessed values.
    fudge_factor: float, default [0.55]
        For use with *guess_bonds*. Supply the factor by which atoms must
        overlap each other to be considered a bond.
    lower_bound: float, default [0.1]
        For use with *guess_bonds*. Supply the minimum bond length.
    transformations: function or list, ``None``, default ``None``
        Provide a list of transformations that you wish to apply to the
        trajectory upon reading. Transformations can be found in
        :mod:`MDAnalysis.transformations`, or can be user-created.
    in_memory: bool, default ``False``
        After reading in the trajectory, transfer it to an in-memory
        representations, which allow for manipulation of coordinates.
    in_memory_step: int, default 1
        Only read every nth frame into in-memory representation.
    continuous: bool, default ``False``
        The `continuous` option is used by the
        :mod:`ChainReader<MDAnalysis.coordinates.chain>`, which contains the
        functionality to treat independent trajectory files as a single virtual
        trajectory.
    **kwargs: extra arguments are passed to the topology parser.

    Attributes
    ----------
    trajectory : base.ReaderBase or base.SingleFrameReaderBase
        currently loaded trajectory reader; readers are described in
        :ref:`Coordinates`
    dimensions : numpy.ndarray
        system dimensions (simulation unit cell, if set in the
        trajectory) at the *current time step*
        (see :attr:`MDAnalysis.coordinates.base.Timestep.dimensions`).
        The unit cell can be set for the current time step (but the change is
        not permanent unless written to a file).
    atoms : AtomGroup
        all particles (:class:`~MDAnalysis.core.groups.Atom`) in the system,
        as read from the `topology` file
    residues : ResidueGroup
        all residues (:class:`~MDAnalysis.core.groups.Residue`) in the system
    segments : SegmentGroup
        all segments (:class:`~MDAnalysis.core.groups.Segment`) in the system
    bonds : topologyattrs.Bonds
        all bonds (if defined in the `topology`) as provided by
        :attr:`Universe.atoms.bonds`
    angles : topologyattrs.Angles
        all angles (if defined in the `topology`), same as
        :attr:`Universe.atoms.angles`
    dihedrals : topologyattrs.Dihedrals
        all dihedral angles (if defined in the `topology`), same as
        :attr:`Universe.atoms.dihedrals`
    impropers : topologyattrs.Impropers
        all improper dihedral angles (if defined in the `topology`), same as
        :attr:`Universe.atoms.impropers`

    Examples
    --------
    Examples for setting up a :class:`Universe`::

       u = Universe(topology, trajectory)          # read system from file(s)
       u = Universe(pdbfile)                       # read atoms and coordinates from PDB or GRO
       u = Universe(topology, [traj1, traj2, ...]) # read from a list of trajectories
       u = Universe(topology, traj1, traj2, ...)   # read from multiple trajectories

    Load new data into a universe (replaces old trajectory and does *not* append)::

       u.load_new(trajectory)                      # read from a new trajectory file

    Selecting atoms with :meth:`~Universe.select_atoms` ::

       ag = u.select_atoms(...)

    returns an :class:`~MDAnalysis.core.groups.AtomGroup`.


    .. versionchanged:: 1.0.0
        Universe() now raises an error. Use Universe(None) or :func:`Universe.empty()` instead.
        Removed instant selectors.

    .. versionchanged:: 2.0.0
        Universe now can be (un)pickled.
        ``topology`` and ``trajectory`` are reserved upon unpickle.
    .. versionchanged:: 2.5.0
        Added fudge_factor and lower_bound parameters for use with
        *guess_bonds*.

    .. versionchanged:: 2.8.0
        Added :meth:`~MDAnalysis.core.universe.Universe.guess_TopologyAttrs` API
        guessing masses and atom types after topology
        is read from a registered parser.

    """
    def __init__(self, topology=None, *coordinates, all_coordinates=False,
                 format=None, topology_format=None, transformations=None,
                 guess_bonds=False, vdwradii=None, fudge_factor=0.55,
                 lower_bound=0.1, in_memory=False, context='default',
                 to_guess=('types', 'masses'), force_guess=(),
                 in_memory_step=1, **kwargs):

        self._trajectory = None  # managed attribute holding Reader
        self._cache = {'_valid': {}}
        self.atoms = None
        self.residues = None
        self.segments = None
        self.filename = None
        self._context = get_guesser(context)
        self._kwargs = {
            'transformations': transformations,
            'guess_bonds': guess_bonds,
            'vdwradii': vdwradii,
            'fudge_factor': fudge_factor,
            'lower_bound': lower_bound,
            'in_memory': in_memory,
            'in_memory_step': in_memory_step,
            'format': format,
            'topology_format': topology_format,
            'all_coordinates': all_coordinates
        }
        self._kwargs.update(kwargs)

        format, topology_format = _resolve_formats(*coordinates, format=format,
                                                   topology_format=topology_format)

        if not isinstance(topology, Topology) and not topology is None:
            self.filename = _check_file_like(topology)
            topology = _topology_from_file_like(self.filename,
                                                topology_format=topology_format,
                                                **kwargs)

        if topology is not None:
            self._topology = topology
        else:
            # point to Universe.empty instead of making empty universe
            raise TypeError('Topology argument required to make Universe. '
                            'Try Universe.empty(n_atoms, ...) to construct '
                            'your own Universe.')

        _generate_from_topology(self)  # make real atoms, res, segments

        coordinates = _resolve_coordinates(self.filename, *coordinates,
                                           format=format,
                                           all_coordinates=all_coordinates)

        if coordinates:
            self.load_new(coordinates, format=format, in_memory=in_memory,
                        in_memory_step=in_memory_step, **kwargs)

        if transformations:
            if callable(transformations):
                transformations = [transformations]
            self._trajectory.add_transformations(*transformations)

        if guess_bonds:
            force_guess = list(force_guess) + ['bonds', 'angles', 'dihedrals']

        self.guess_TopologyAttrs(
            context, to_guess, force_guess, vdwradii=vdwradii, **kwargs)


    def copy(self):
        """Return an independent copy of this Universe"""
        context = self._context.copy()
        new = self.__class__(self._topology.copy(),
                             to_guess=(), context=context)
        new.trajectory = self.trajectory.copy()
        return new

    @classmethod
    def empty(cls, n_atoms, n_residues=1, n_segments=1, n_frames=1,
              atom_resindex=None, residue_segindex=None,
              trajectory=False, velocities=False, forces=False):
        """Create a blank Universe

        Useful for building a Universe without requiring existing files,
        for example for system building.

        If `trajectory` is set to ``True``, a
        :class:`MDAnalysis.coordinates.memory.MemoryReader` will be
        attached to the Universe.

        Parameters
        ----------
        n_atoms: int
          number of Atoms in the Universe
        n_residues: int, default 1
          number of Residues in the Universe, defaults to 1
        n_segments: int, default 1
          number of Segments in the Universe, defaults to 1
        n_frames: int, default 1
          number of Frames in the Universe, defaults to 1
        atom_resindex: array like, optional
          mapping of atoms to residues, e.g. with 6 atoms,
          `atom_resindex=[0, 0, 1, 1, 2, 2]` would put 2 atoms
          into each of 3 residues.
        residue_segindex: array like, optional
          mapping of residues to segments
        trajectory: bool, optional
          if ``True``, attaches a
          :class:`MDAnalysis.coordinates.memory.MemoryReader` allowing
          coordinates to be set and written.
        velocities: bool, optional
          include velocities in the
          :class:`MDAnalysis.coordinates.memory.MemoryReader`
        forces: bool, optional
          include forces in the
          :class:`MDAnalysis.coordinates.memory.MemoryReader`

        Returns
        -------
        Universe
          :class:`~MDAnalysis.core.universe.Universe` instance with dummy
          values for atoms and undefined coordinates/velocities/forces

        Examples
        --------
        For example to create a new Universe with 6 atoms in 2 residues, with
        positions for the atoms and a mass attribute::

          u = mda.Universe.empty(6, 2,
                                 atom_resindex=np.array([0, 0, 0, 1, 1, 1]),
                                 trajectory=True,
                )
          u.add_TopologyAttr('masses')

        .. versionadded:: 0.17.0
        .. versionchanged:: 0.19.0
           The attached Reader when trajectory=True is now a MemoryReader
        .. versionchanged:: 1.0.0
           Universes can now be created with 0 atoms

        """
        if not n_atoms:
            n_residues = 0
            n_segments = 0

        if atom_resindex is None and n_residues > 1:
            warnings.warn(
                'Residues specified but no atom_resindex given.  '
                'All atoms will be placed in first Residue.',
                UserWarning)

        if residue_segindex is None and n_segments > 1:
            warnings.warn(
                'Segments specified but no segment_resindex given.  '
                'All residues will be placed in first Segment',
                UserWarning)

        top = Topology(n_atoms, n_residues, n_segments,
                       atom_resindex=atom_resindex,
                       residue_segindex=residue_segindex,
        )

        u = cls(top, to_guess=())

        if n_frames > 1 or trajectory:
            coords = np.zeros((n_frames, n_atoms, 3), dtype=np.float32)
            vels = np.zeros_like(coords) if velocities else None
            forces = np.zeros_like(coords) if forces else None

            # grab and attach a MemoryReader
            u.trajectory = get_reader_for(coords)(
                coords, order='fac', n_atoms=n_atoms,
                velocities=vels, forces=forces)

        return u

    @property
    def universe(self):
        # for Writer.write(universe), see Issue 49
        # Encapsulation in an accessor prevents the Universe from
        # having to keep a reference to itself,
        #  which might be undesirable if it has a __del__ method.
        # It is also cleaner than a weakref.
        return self

    def load_new(self, filename, format=None, in_memory=False,
                 in_memory_step=1, **kwargs):
        """Load coordinates from `filename`.

        The file format of `filename` is autodetected from the file name suffix
        or can be explicitly set with the `format` keyword. A sequence of files
        can be read as a single virtual trajectory by providing a list of
        filenames.


        Parameters
        ----------
        filename: str or list
            the coordinate file (single frame or trajectory) *or* a list of
            filenames, which are read one after another.
        format: str or list or object (optional)
            provide the file format of the coordinate or trajectory file;
            ``None`` guesses it from the file extension. Note that this
            keyword has no effect if a list of file names is supplied because
            the "chained" reader has to guess the file format for each
            individual list member [``None``]. Can also pass a subclass of
            :class:`MDAnalysis.coordinates.base.ProtoReader` to define a custom
            reader to be used on the trajectory file.
        in_memory: bool (optional)
            Directly load trajectory into memory with the
            :class:`~MDAnalysis.coordinates.memory.MemoryReader`

            .. versionadded:: 0.16.0

        **kwargs: dict
            Other kwargs are passed to the trajectory reader (only for
            advanced use)

        Returns
        -------
        universe: Universe

        Raises
        ------
        TypeError
             if trajectory format can not be determined or no appropriate
             trajectory reader found


        .. versionchanged:: 0.8
           If a list or sequence that is provided for `filename` only contains
           a single entry then it is treated as single coordinate file. This
           has the consequence that it is not read by the
           :class:`~MDAnalysis.coordinates.chain.ChainReader` but directly by
           its specialized file format reader, which typically has more
           features than the
           :class:`~MDAnalysis.coordinates.chain.ChainReader`.

        .. versionchanged:: 0.17.0
           Now returns a :class:`Universe` instead of the tuple of file/array
           and detected file type.
        .. versionchanged:: 2.4.0
           Passes through kwargs if `in_memory=True`.

        """
        # filename==None happens when only a topology is provided
        if filename is None:
            return self

        if len(util.asiterable(filename)) == 1:
            # make sure a single filename is not handed to the ChainReader
            filename = util.asiterable(filename)[0]
        logger.debug("Universe.load_new(): loading {0}...".format(filename))

        try:
            reader = get_reader_for(filename, format=format)
        except ValueError as err:
            raise TypeError(
                "Cannot find an appropriate coordinate reader for file '{0}'.\n"
                "           {1}".format(filename, err))

        # supply number of atoms for readers that cannot do it for themselves
        kwargs['n_atoms'] = self.atoms.n_atoms

        self.trajectory = reader(filename, format=format, **kwargs)
        if self.trajectory.n_atoms != len(self.atoms):
            raise ValueError("The topology and {form} trajectory files don't"
                             " have the same number of atoms!\n"
                             "Topology number of atoms {top_n_atoms}\n"
                             "Trajectory: {fname} Number of atoms {trj_n_atoms}".format(
                                 form=self.trajectory.format,
                                 top_n_atoms=len(self.atoms),
                                 fname=filename,
                                 trj_n_atoms=self.trajectory.n_atoms))

        if in_memory:
            self.transfer_to_memory(step=in_memory_step, **kwargs)

        return self

    def transfer_to_memory(self, start=None, stop=None, step=None,
                           verbose=False, **kwargs):
        """Transfer the trajectory to in memory representation.

        Replaces the current trajectory reader object with one of type
        :class:`MDAnalysis.coordinates.memory.MemoryReader` to support in-place
        editing of coordinates.

        Parameters
        ----------
        start: int, optional
            start reading from the nth frame.
        stop: int, optional
            read upto and excluding the nth frame.
        step: int, optional
            Read in every nth frame. [1]
        verbose: bool, optional
            Will print the progress of loading trajectory to memory, if
            set to True. Default value is False.


        .. versionadded:: 0.16.0
        .. versionchanged:: 2.4.0
           Passes through kwargs to MemoryReader
        """
        from ..coordinates.memory import MemoryReader

        if not isinstance(self.trajectory, MemoryReader):
            n_frames = len(range(
                *self.trajectory.check_slice_indices(start, stop, step)
            ))
            n_atoms = len(self.atoms)
            coordinates = np.zeros((n_frames, n_atoms, 3), dtype=np.float32)
            ts = self.trajectory.ts
            has_vels = ts.has_velocities
            has_fors = ts.has_forces
            has_dims = ts.dimensions is not None

            velocities = np.zeros_like(coordinates) if has_vels else None
            forces = np.zeros_like(coordinates) if has_fors else None
            dimensions = (np.zeros((n_frames, 6), dtype=np.float32)
                          if has_dims else None)

            for i, ts in enumerate(ProgressBar(self.trajectory[start:stop:step],
                                               verbose=verbose,
                                               desc="Loading frames")):
                np.copyto(coordinates[i], ts.positions)
                if has_vels:
                    np.copyto(velocities[i], ts.velocities)
                if has_fors:
                    np.copyto(forces[i], ts.forces)
                if has_dims:
                    np.copyto(dimensions[i], ts.dimensions)

            # Overwrite trajectory in universe with an MemoryReader
            # object, to provide fast access and allow coordinates
            # to be manipulated
            if step is None:
                step = 1
            self.trajectory = MemoryReader(
                coordinates,
                dimensions=dimensions,
                dt=self.trajectory.ts.dt * step,
                filename=self.trajectory.filename,
                velocities=velocities,
                forces=forces, **kwargs)

    # python 2 doesn't allow an efficient splitting of kwargs in function
    # argument signatures.
    # In python3-only we'd be able to explicitly define this function with
    # something like (sel, *othersels, updating=False, **selgroups)
    def select_atoms(self, *args, **kwargs):
        """Select atoms.

        See Also
        --------
        :meth:`MDAnalysis.core.groups.AtomGroup.select_atoms`
        """
        return self.atoms.select_atoms(*args, **kwargs)

    @property
    def bonds(self):
        """Bonds between atoms.

        :meta private:
        """
        return self.atoms.bonds

    @property
    def angles(self):
        """Angles between atoms.

        :meta private:
        """
        return self.atoms.angles

    @property
    def dihedrals(self):
        """Dihedral angles between atoms.

        :meta private:
        """
        return self.atoms.dihedrals

    @property
    def impropers(self):
        """Improper dihedral angles between atoms.

        :meta private:
        """
        return self.atoms.impropers

    def __repr__(self):
        # return "<Universe with {n_atoms} atoms{bonds}>".format(
        #    n_atoms=len(self.atoms),
        #    bonds=" and {0} bonds".format(len(self.bonds)) if self.bonds else "")

        return "<Universe with {n_atoms} atoms>".format(
            n_atoms=len(self.atoms))

    @classmethod
    def _unpickle_U(cls, top, traj, context):
        """Special method used by __reduce__ to deserialise a Universe"""
        #  top is a Topology obj at this point, but Universe can handle that.
        u = cls(top, to_guess=(), context=context)
        u.trajectory = traj

        return u

    def __reduce__(self):
        #  __setstate__/__getstate__ will raise an error when Universe has a
        #  transformation (that has AtomGroup inside). Use __reduce__ instead.
        #  Universe's two "legs" of top and traj both serialise themselves.
        return (self._unpickle_U, (self._topology,
                                   self._trajectory, self._context.copy()))

    # Properties
    @property
    def dimensions(self):
        """Current dimensions of the unitcell.

        :meta private:
        """
        return self.coord.dimensions

    @dimensions.setter
    def dimensions(self, box):
        """Set dimensions if the Timestep allows this

        .. versionadded:: 0.9.0
        """
        # Add fancy error handling here or use Timestep?
        self.coord.dimensions = box

    @property
    def coord(self):
        """Reference to current timestep and coordinates of universe.

        The raw trajectory coordinates are :attr:`Universe.coord.positions`,
        represented as a :class:`numpy.float32` array.

        Because :attr:`coord` is a reference to a
        :class:`~MDAnalysis.coordinates.timestep.Timestep`, it changes its contents
        while one is stepping through the trajectory.

        .. Note::

           In order to access the coordinates it is better to use the
           :meth:`AtomGroup.positions` method; for instance, all coordinates of
           the Universe as a numpy array: :meth:`Universe.atoms.positions`.

        """
        return self.trajectory.ts

    @property
    def kwargs(self):
        """keyword arguments used to initialize this universe"""
        return copy.deepcopy(self._kwargs)

    @property
    def trajectory(self):
        """Reference to trajectory reader object containing trajectory data.

        :meta private:
        """
        if self._trajectory is not None:
            return self._trajectory
        else:
            raise AttributeError("No trajectory loaded into Universe")

    @trajectory.setter
    def trajectory(self, value):
        del self._trajectory  # guarantees that files are closed (?)
        self._trajectory = value

    def add_TopologyAttr(self, topologyattr, values=None):
        """Add a new topology attribute to the Universe

        Adding a TopologyAttribute to the Universe makes it available to
        all AtomGroups etc throughout the Universe.

        Parameters
        ----------
        topologyattr: TopologyAttr or string
          Either a MDAnalysis TopologyAttr object or the name of a possible
          topology attribute.
        values: np.ndarray, optional
          If initiating an attribute from a string, the initial values to
          use.  If not supplied, the new TopologyAttribute will have empty
          or zero values.

        Example
        -------
        For example to add bfactors to a Universe:

        >>> import MDAnalysis as mda
        >>> from MDAnalysis.tests.datafiles import PSF, DCD
        >>> u = mda.Universe(PSF, DCD)
        >>> u.add_TopologyAttr('tempfactors')
        >>> u.atoms.tempfactors
        array([0., 0., 0., ..., 0., 0., 0.])

        .. versionchanged:: 0.17.0
           Can now also add TopologyAttrs with a string of the name of the
           attribute to add (eg 'charges'), can also supply initial values
           using values keyword.

        .. versionchanged:: 1.1.0
            Now warns when adding bfactors to a Universe with
            existing tempfactors, or adding tempfactors to a
            Universe with existing bfactors.
            In version 2.0, MDAnalysis will stop treating
            tempfactors and bfactors as separate attributes. Instead,
            they will be aliases of the same attribute.
        """
        if isinstance(topologyattr, str):
            if topologyattr in ("bfactor", "bfactors"):
                warnings.warn(BFACTOR_WARNING, DeprecationWarning)
            try:
                tcls = _TOPOLOGY_ATTRS[topologyattr]
            except KeyError:
                errmsg = (
                    "Unrecognised topology attribute name: '{}'."
                    "  Possible values: '{}'\n"
                    "To raise an issue go to: "
                    "https://github.com/MDAnalysis/mdanalysis/issues"
                    "".format(
                        topologyattr, ', '.join(
                            sorted(_TOPOLOGY_ATTRS.keys()))))
                raise ValueError(errmsg) from None
            else:
                topologyattr = tcls.from_blank(
                    n_atoms=self._topology.n_atoms,
                    n_residues=self._topology.n_residues,
                    n_segments=self._topology.n_segments,
                    values=values)
        self._topology.add_TopologyAttr(topologyattr)
        self._process_attr(topologyattr)

    def del_TopologyAttr(self, topologyattr):
        """Remove a topology attribute from the Universe

        Removing a TopologyAttribute from the Universe makes it unavailable to
        all AtomGroups etc throughout the Universe.

        Parameters
        ----------
        topologyattr: TopologyAttr or string
          Either a MDAnalysis TopologyAttr object or the name of a possible
          topology attribute.

        Example
        -------
        For example to remove bfactors to a Universe:

        >>> import MDAnalysis as mda
        >>> from MDAnalysis.tests.datafiles import PSF, DCD
        >>> u = mda.Universe(PSF, DCD)
        >>> u.add_TopologyAttr('tempfactors')
        >>> hasattr(u.atoms[:3], 'tempfactors')
        True
        >>>
        >>> u.del_TopologyAttr('tempfactors')
        >>> hasattr(u.atoms[:3], 'tempfactors')
        False


        .. versionadded:: 2.0.0
        """

        if not isinstance(topologyattr, str):
            try:
                topologyattr = topologyattr.attrname
            except AttributeError:
                # either TopologyGroup or not valid
                try:
                    # this may not end well
                    # e.g. matrix -> matrices
                    topologyattr = topologyattr.btype + "s"
                except AttributeError:
                    raise ValueError("Topology attribute must be str or "
                                     "TopologyAttr object or class. "
                                     f"Given: {type(topologyattr)}") from None

        try:
            topologyattr = _TOPOLOGY_ATTRS[topologyattr].attrname
        except KeyError:
            attrs = ', '.join(sorted(_TOPOLOGY_ATTRS))
            errmsg = (f"Unrecognised topology attribute: '{topologyattr}'."
                      f"  Possible values: '{attrs}'\n"
                      "To raise an issue go to: "
                      "https://github.com/MDAnalysis/mdanalysis/issues")
            raise ValueError(errmsg) from None

        try:
            topattr = getattr(self._topology, topologyattr)
        except AttributeError:
            raise ValueError(f"Topology attribute {topologyattr} "
                             "not in Universe.") from None
        self._topology.del_TopologyAttr(topattr)
        self._unprocess_attr(topattr)

    def _process_attr(self, attr):
        """Squeeze a topologyattr for its information

        Grabs:
         - Group properties (attribute access)
         - Component properties
         - Transplant methods
        """
        n_dict = {'atom': self._topology.n_atoms,
                  'residue': self._topology.n_residues,
                  'segment': self._topology.n_segments}
        logger.debug("_process_attr: Adding {0} to topology".format(attr))
        if (attr.per_object is not None and len(attr) != n_dict[attr.per_object]):
            raise ValueError('Length of {attr} does not'
                             ' match number of {obj}s.\n'
                             'Expect: {n:d} Have: {m:d}'.format(
                                 attr=attr.attrname,
                                 obj=attr.per_object,
                                 n=n_dict[attr.per_object],
                                 m=len(attr)))

        for cls in attr.target_classes:
            self._class_bases[cls]._add_prop(attr)

        # TODO: Try and shove this into cls._add_prop
        # Group transplants
        for cls in (Atom, Residue, Segment, GroupBase,
                    AtomGroup, ResidueGroup, SegmentGroup):
            for funcname, meth in attr.transplants[cls]:
                setattr(self._class_bases[cls], funcname, meth)
        # Universe transplants
        for funcname, meth in attr.transplants['Universe']:
            setattr(self.__class__, funcname, meth)

    def _unprocess_attr(self, attr):
        """
        Undo all the stuff in _process_attr.

        If the topology attribute is not present, nothing happens
        (silent fail).
        """
        for cls in attr.target_classes:
            self._class_bases[cls]._del_prop(attr)

        # Universe transplants
        for funcname, _ in attr.transplants.pop("Universe", []):
            delattr(self.__class__, funcname)
        # Group transplants
        for cls, transplants in attr.transplants.items():
            for funcname, _ in transplants:
                delattr(self._class_bases[cls], funcname)

    def add_Residue(self, segment=None, **attrs):
        """Add a new Residue to this Universe

        New Residues will not contain any Atoms, but can be assigned to Atoms
        as per usual.  If the Universe contains multiple segments, this must
        be specified as a keyword.

        Parameters
        ----------
        segment: MDAnalysis.Segment
          If there are multiple segments, then the Segment that the new
          Residue will belong in must be specified.
        attrs: dict
          For each Residue attribute, the value for the new Residue must be
          specified

        Returns
        -------
        A reference to the new Residue

        Raises
        ------
        NoDataError
          If any information was missing.  This happens before any changes have
          been made, ie the change is rolled back.


        Example
        -------

        Adding a new GLY residue, then placing atoms within it:

        >>> import MDAnalysis as mda
        >>> from MDAnalysis.tests.datafiles import PSF, DCD
        >>> u = mda.Universe(PSF, DCD)
        >>> newres = u.add_Residue(segment=u.segments[0], resid=42, resname='GLY', resnum=0)
        >>> u.atoms[[1, 2, 3]].residues = newres
        >>> u.select_atoms('resname GLY and resid 42 and resnum 0')
        <AtomGroup with 3 atoms>

        """
        if len(self.segments) == 1:  # if only one segment, use this
            segment = self.segments[0]
        if segment is None:
            raise NoDataError("")
        # pass this information to the topology
        residx = self._topology.add_Residue(segment, **attrs)
        # resize my residues
        self.residues = ResidueGroup(np.arange(self._topology.n_residues), self)

        # return the new residue
        return self.residues[residx]

    def add_Segment(self, **attrs):
        """Add a new Segment to this Universe

        Parameters
        ----------
        attrs: dict
            For each Segment attribute as a key, give the value in the new
            Segment

        Returns
        -------
        A reference to the new Segment

        Raises
        ------
        NoDataError
            If any attributes were not specified as a keyword.

        """
        # pass this information to the topology
        segidx = self._topology.add_Segment(**attrs)
        # resize my segments
        self.segments = SegmentGroup(np.arange(self._topology.n_segments), self)
        # return the new segment
        return self.segments[segidx]

    def _add_topology_objects(self, object_type, values, types=None, guessed=False,
                           order=None):
        """Add new TopologyObjects to this Universe

        Parameters
        ----------
        object_type : {'bonds', 'angles', 'dihedrals', 'impropers'}
            The type of TopologyObject to add.
        values : TopologyGroup or iterable of tuples, AtomGroups, or TopologyObjects
            An iterable of: tuples of atom indices, or AtomGroups,
            or TopologyObjects. If every value is a TopologyObject, all
            keywords are ignored.
            If AtomGroups or TopologyObjects are passed, they *must* be from the same
            Universe.
        types : iterable (optional, default None)
            None, or an iterable of hashable values with the same length as ``values``
        guessed : bool or iterable (optional, default False)
            bool, or an iterable of hashable values with the same length as ``values``
        order : iterable (optional, default None)
            None, or an iterable of hashable values with the same length as ``values``


        .. versionadded:: 1.0.0
        """
        if all(isinstance(x, TopologyObject) for x in values):
            try:
                types = [t.type for t in values]
            except AttributeError:
                types = None
            guessed = [t.is_guessed for t in values]
            order = [t.order for t in values]

        indices = []
        for x in values:
            if isinstance(x, (AtomGroup, TopologyObject)):
                if x.universe is not self:
                    err_msg = 'Cannot add {} from different Universes.'
                    raise ValueError(err_msg.format(object_type))
                indices.append(x.indices)
            else:
                indices.append(x)

        all_indices = set([i for obj in indices for i in obj])
        nonexistent = all_indices - set(self.atoms.indices)
        if nonexistent:
            istr = ', '.join(map(str, nonexistent))
            err_msg = 'Cannot add {} for nonexistent atom indices: {}'
            raise ValueError(err_msg.format(object_type, istr))

        try:
            attr = getattr(self._topology, object_type)
        except AttributeError:
            self.add_TopologyAttr(object_type, [])
            attr = getattr(self._topology, object_type)


        attr._add_bonds(indices, types=types, guessed=guessed, order=order)

    def add_bonds(self, values, types=None, guessed=False, order=None):
        """Add new Bonds to this Universe.

        Parameters
        ----------
        values : iterable of tuples, AtomGroups, or Bonds; or TopologyGroup
            An iterable of: tuples of 2 atom indices, or AtomGroups with 2 atoms,
            or Bonds. If every value is a Bond, all
            keywords are ignored.
            If AtomGroups, Bonds, or a TopologyGroup are passed,
            they *must* be from the same Universe.
        types : iterable (optional, default None)
            None, or an iterable of hashable values with the same length as ``values``
        guessed : bool or iterable (optional, default False)
            bool, or an iterable of hashable values with the same length as ``values``
        order : iterable (optional, default None)
            None, or an iterable of hashable values with the same length as ``values``


        Example
        -------

        Adding TIP4P water bonds with a list of AtomGroups::

            import MDAnalysis as mda
            from MDAnalysis.tests.datafiles import GRO
            u = mda.Universe(GRO)
            sol = u.select_atoms('resname SOL')
            ow_hw1 = sol.select_atoms('name OW or name HW1').split('residue')
            ow_hw2 = sol.select_atoms('name OW or name HW2').split('residue')
            ow_mw = sol.select_atoms('name OW or name MW').split('residue')
            u.add_bonds(ow_hw1 + ow_hw2 + ow_mw)

        You can only add bonds from the same Universe. If you would like to add
        AtomGroups, Bonds, or a TopologyGroup from a different Universe, convert
        them to indices first. ::

            from MDAnalysis.tests.datafiles import PSF
            u2 = mda.Universe(PSF)

            #  assuming you have already added bonds to u
            u2.add_bonds(u.bonds.to_indices())


        .. versionadded:: 1.0.0
        """
        self._add_topology_objects('bonds', values, types=types,
                                 guessed=guessed, order=order)
        # Invalidate bond-related caches
        self._cache.pop('fragments', None)
        self._cache['_valid'].pop('fragments', None)
        self._cache['_valid'].pop('fragindices', None)

    def add_angles(self, values, types=None, guessed=False):
        """Add new Angles to this Universe.

        Parameters
        ----------
        values : iterable of tuples, AtomGroups, or Angles; or TopologyGroup
            An iterable of: tuples of 3 atom indices, or AtomGroups with 3 atoms,
            or Angles. If every value is a Angle, all
            keywords are ignored.
            If AtomGroups, Angles, or a TopologyGroup are passed,
            they *must* be from the same Universe.
        types : iterable (optional, default None)
            None, or an iterable of hashable values with the same length as ``values``
        guessed : bool or iterable (optional, default False)
            bool, or an iterable of hashable values with the same length as ``values``

        .. versionadded:: 1.0.0
        """
        self._add_topology_objects('angles', values, types=types,
                                 guessed=guessed)

    def add_dihedrals(self, values, types=None, guessed=False):
        """Add new Dihedrals to this Universe.

        Parameters
        ----------
        values : iterable of tuples, AtomGroups, or Dihedrals; or TopologyGroup
            An iterable of: tuples of 4 atom indices, or AtomGroups with 4 atoms,
            or Dihedrals. If every value is a Dihedral, all
            keywords are ignored.
            If AtomGroups, Dihedrals, or a TopologyGroup are passed,
            they *must* be from the same Universe.
        types : iterable (optional, default None)
            None, or an iterable of hashable values with the same length as ``values``
        guessed : bool or iterable (optional, default False)
            bool, or an iterable of hashable values with the same length as ``values``


        .. versionadded:: 1.0.0
        """
        self._add_topology_objects('dihedrals', values, types=types,
                                 guessed=guessed)

    def add_impropers(self, values, types=None, guessed=False):
        """Add new Impropers to this Universe.

        Parameters
        ----------
        values : iterable of tuples, AtomGroups, or Impropers; or TopologyGroup
            An iterable of: tuples of 4 atom indices, or AtomGroups with 4 atoms,
            or Impropers. If every value is an Improper, all
            keywords are ignored.
            If AtomGroups, Impropers, or a TopologyGroup are passed,
            they *must* be from the same Universe.
        types : iterable (optional, default None)
            None, or an iterable of hashable values with the same length as ``values``
        guessed : bool or iterable (optional, default False)
            bool, or an iterable of hashable values with the same length as ``values``


        .. versionadded:: 1.0.0
        """
        self._add_topology_objects('impropers', values, types=types,
                                 guessed=guessed)

    def _delete_topology_objects(self, object_type, values):
        """Delete TopologyObjects from this Universe

        Parameters
        ----------
        object_type : {'bonds', 'angles', 'dihedrals', 'impropers'}
            The type of TopologyObject to add.
        values : iterable of tuples, AtomGroups, or TopologyObjects; or TopologyGroup
            An iterable of: tuples of atom indices, or AtomGroups,
            or TopologyObjects.
            If AtomGroups, TopologyObjects, or a TopologyGroup are passed,
            they *must* be from the same Universe.

        .. versionadded:: 1.0.0
        """
        indices = []
        for x in values:
            if isinstance(x, (AtomGroup, TopologyObject)):
                if x.universe is not self:
                    err_msg = 'Cannot delete {} from different Universes.'
                    raise ValueError(err_msg.format(object_type))
                indices.append(x.indices)
            else:
                indices.append(x)

        try:
            attr = getattr(self._topology, object_type)
        except AttributeError:
            raise ValueError('There are no {} to delete'.format(object_type))

        attr._delete_bonds(indices)

    def delete_bonds(self, values):
        """Delete Bonds from this Universe.

        Parameters
        ----------
        values : iterable of tuples, AtomGroups, or Bonds; or TopologyGroup
            An iterable of: tuples of 2 atom indices, or AtomGroups with 2 atoms,
            or Bonds.
            If AtomGroups, Bonds, or a TopologyGroup are passed,
            they *must* be from the same Universe.


        Example
        -------

        Deleting bonds from a Universe::

            import MDAnalysis as mda
            from MDAnalysis.tests.datafiles import PSF
            u = mda.Universe(PSF)

            #  delete first 5 bonds
            u.delete_bonds(u.bonds[:5])


        If you are deleting bonds in the form of AtomGroups, Bonds, or a
        TopologyGroup, they must come from the same Universe. If you want to
        delete bonds from another Universe, convert them to indices first. ::

            from MDAnalysis.tests.datafiles import PDB
            u2 = mda.Universe(PDB)

            u.delete_bonds(u2.bonds.to_indices())


        .. versionadded:: 1.0.0
        """
        self._delete_topology_objects('bonds', values)
        # Invalidate bond-related caches
        self._cache.pop('fragments', None)
        self._cache['_valid'].pop('fragments', None)
        self._cache['_valid'].pop('fragindices', None)

    def delete_angles(self, values):
        """Delete Angles from this Universe.

        Parameters
        ----------
        values : iterable of tuples, AtomGroups, or Angles; or TopologyGroup
            An iterable of: tuples of 3 atom indices, or AtomGroups with 3 atoms,
            or Angles.
            If AtomGroups, Angles, or a TopologyGroup are passed,
            they *must* be from the same Universe.


        .. versionadded:: 1.0.0
        """
        self._delete_topology_objects('angles', values)

    def delete_dihedrals(self, values):
        """Delete Dihedrals from this Universe.

        Parameters
        ----------
        values : iterable of tuples, AtomGroups, or Dihedrals; or TopologyGroup
            An iterable of: tuples of 4 atom indices, or AtomGroups with 4 atoms,
            or Dihedrals.
            If AtomGroups, Dihedrals, or a TopologyGroup are passed,
            they *must* be from the same Universe.


        .. versionadded:: 1.0.0
        """
        self._delete_topology_objects('dihedrals', values)

    def delete_impropers(self, values):
        """Delete Impropers from this Universe.

        Parameters
        ----------
        values : iterable of tuples, AtomGroups, or Impropers; or TopologyGroup
            An iterable of: tuples of 4 atom indices, or AtomGroups with 4 atoms,
            or Impropers.
            If AtomGroups, Angles, or a TopologyGroup are passed,
            they *must* be from the same Universe.


        .. versionadded:: 1.0.0
        """
        self._delete_topology_objects('impropers', values)

    # TODO: Maybe put this as a Bond attribute transplant
    # Problems: Can we transplant onto Universe?
    # Probably a smarter way to do this too, could generate
    # these on demand *per atom*.
    # Wouldn't then need the Universe linkage here
    #
    # Alternate idea: Bonds Attribute generates a Fragments
    # Attribute (ie, 2 for the price of 1)
    # Fragments then gets its own Class/namespace/jazz.
    @property
    @cached('fragments')
    def _fragdict(self):
        """
        .. versionadded:: 0.9.0
        .. versionchanged:: 0.16.0
           Fragment atoms are sorted by their index, and framgents are sorted
           by their first atom index so their order is predictable.
        .. versionchanged:: 0.19.0
           Uses faster C++ implementation
        .. versionchanged:: 0.20.0
           * _fragdict keys are now atom indices instead of Atoms
           * _fragdict items are now a namedtuple ``fraginfo(ix, fragment)``
             storing the fragindex ``ix`` along with the fragment.
        """
        atoms = self.atoms.ix
        bonds = self.atoms.bonds.to_indices()

        frag_indices = find_fragments(atoms, bonds)
        frags = tuple([AtomGroup(np.sort(ix), self) for ix in frag_indices])

        fragdict = {}
        fraginfo = collections.namedtuple('fraginfo', 'ix, fragment')
        for i, f in enumerate(frags):
            for a in f:
                fragdict[a.ix] = fraginfo(i, f)

        return fragdict

    @classmethod
    def from_smiles(cls, smiles, sanitize=True, addHs=True,
                    generate_coordinates=True, numConfs=1,
                    rdkit_kwargs={}, **kwargs):
        """Create a Universe from a SMILES string with rdkit

        Parameters
        ----------
        smiles : str
            SMILES string
        sanitize : bool (optional, default True)
            Toggle the sanitization of the molecule
        addHs : bool (optional, default True)
            Add all necessary hydrogens to the molecule
        generate_coordinates : bool (optional, default True)
            Generate 3D coordinates using RDKit's
            :func:`AllChem.EmbedMultipleConfs` function. Requires adding
            hydrogens with the `addHs` parameter
        numConfs : int (optional, default 1)
            Number of frames to generate coordinates for. Ignored if
            ``generate_coordinates=False``
        rdkit_kwargs : dict (optional)
            Other arguments passed to the RDKit :func:`EmbedMultipleConfs`
            function
        kwargs : dict
            Parameters passed on Universe creation

        Returns
        -------
        universe : Universe
            contains atom names and topology information (bonds) derived from
            the input SMILES string; coordinates are included if
            `generate_coordinates` was set to ``True``


        Examples
        --------
        To create a Universe with 10 conformers of ethanol:

        >>> from rdkit.Chem import AllChem
        >>> u = mda.Universe.from_smiles('CCO', numConfs=10)
        >>> u
        <Universe with 9 atoms>
        >>> u.trajectory
        <RDKitReader with 10 frames of 9 atoms>

        To use a different conformer generation algorithm, like ETKDGv3:

        >>> u = mda.Universe.from_smiles('CCO', rdkit_kwargs=dict(
        ...      params=AllChem.ETKDGv3()))
        >>> u.trajectory
        <RDKitReader with 1 frames of 9 atoms>


        .. versionadded:: 2.0.0

        """
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
        except ImportError as e:
            raise ImportError(
                "Creating a Universe from a SMILES string requires RDKit but "
                "it does not appear to be installed") from e

        mol = Chem.MolFromSmiles(smiles, sanitize=sanitize)
        if mol is None:
            raise SyntaxError('Error while parsing SMILES {0}'.format(smiles))
        if addHs:
            mol = Chem.AddHs(mol)
        if generate_coordinates:
            if not addHs:
                raise ValueError("Generating coordinates requires adding "
                "hydrogens with `addHs=True`")

            numConfs = rdkit_kwargs.pop("numConfs", numConfs)
            if not (isinstance(numConfs, int) and numConfs > 0):
                raise SyntaxError("numConfs must be a non-zero positive "
                                  "integer instead of {0}".format(numConfs))
            AllChem.EmbedMultipleConfs(mol, numConfs, **rdkit_kwargs)

        return cls(mol, **kwargs)

    def guess_TopologyAttrs(
            self, context=None, to_guess=(), force_guess=(), **kwargs):
        """
        Guess and add attributes through a specific context-aware guesser.

        Parameters
        ----------
        context: str or Guesser class
            For calling a matching guesser class for this specific context
        to_guess: list[str], (optional), default ``['types', 'masses']``
                (in future versions types and masses will be removed)
            TopologyAttrs to be guessed. These TopologyAttrs will be wholly guessed
            if they don't exist in the universe. If they already exist in the Universe,
            only empty or missing values will be guessed.
        force_guess: list[str], (optional)
            TopologyAttrs in this list will be force guessed. If the TopologyAttr
            does not already exist in the Universe, this has no effect. If the TopologyAttr
            does already exist, all values will be overwritten by guessed values.
        **kwargs: extra arguments to be passed to the guesser class

        Examples
        --------
        guess ``masses`` and ``types`` attribute::

            u.guess_TopologyAttrs(context='default', to_guess=['masses', 'types'])

        .. versionadded:: 2.8.0

        """
        if not context:
            context = self._context

        guesser = get_guesser(context, self.universe, **kwargs)
        self._context = guesser

        total_guess = list(force_guess) + list(to_guess)

        # Removing duplicates from the guess list while keeping attributes
        # order as it is more convientent to guess attributes
        # in the same order that the user provided
        total_guess = list(dict.fromkeys(total_guess))

        objects = ['bonds', 'angles', 'dihedrals', 'impropers']

        # Checking if the universe is empty to avoid errors
        # from guesser methods
        if self._topology.n_atoms > 0:

            topology_attrs = [att.attrname for att in
                              self._topology.read_attributes]

            common_attrs = set(to_guess) & set(topology_attrs)
            common_attrs = ", ".join(attr for attr in common_attrs)

            if len(common_attrs) > 0:
                logger.info(
                    f'The attribute(s) {common_attrs} have already been read '
                    'from the topology file. The guesser will '
                    'only guess empty values for this attribute, '
                    'if any exists. To overwrite it by completely '
                    'guessed values, you can pass the attribute to'
                    ' the force_guess parameter instead of '
                    'the to_guess one')

            for attr in total_guess:
                if guesser.is_guessable(attr):
                    fg = True if attr in force_guess else False
                    values = guesser.guess_attr(attr, fg)

                    if values is not None:
                        if attr in objects:
                            self._add_topology_objects(
                                attr, values, guessed=True)
                        else:
                            guessed_attr = _TOPOLOGY_ATTRS[attr](values, True)
                            self.add_TopologyAttr(guessed_attr)
                        logger.info(
                            f'attribute {attr} has been guessed'
                            ' successfully.')

                else:
                    raise ValueError(f'{context} guesser can not guess the'
                                     f' following attribute: {attr}')

        else:
            warnings.warn('Can not guess attributes '
                          'for universe with 0 atoms')


def Merge(*args):
    """Create a new new :class:`Universe` from one or more
    :class:`~MDAnalysis.core.groups.AtomGroup` instances.

    Parameters
    ----------
    *args: :class:`~MDAnalysis.core.groups.AtomGroup`
        One or more AtomGroups.

    Returns
    -------
    universe: :class:`Universe`

    Raises
    ------
    ValueError
        Too few arguments or an AtomGroup is empty and
    TypeError
        Arguments are not :class:`AtomGroup` instances.

    Notes
    -----
    The resulting :class:`Universe` will only inherit the common topology
    attributes that all merged universes share.

    :class:`AtomGroup` instances can come from different Universes, or can come
    directly from a :meth:`~Universe.select_atoms` call.

    :class:`Merge` can also be used with a single :class:`AtomGroup` if the
    user wants to, for example, re-order the atoms in the :class:`Universe`.

    If multiple :class:`AtomGroup` instances from the same :class:`Universe`
    are given, the merge will first simply "add" together the
    :class:`AtomGroup` instances.

    Merging does not create a full trajectory but only a single structure even
    if the input consists of one or more trajectories.  However, one can use
    the :class:`~MDAnalysis.coordinates.memory.MemoryReader` to construct a
    trajectory for the new Universe as described under
    :ref:`creating-in-memory-trajectory-label`.

    Example
    -------
    In this example, protein, ligand, and solvent were externally prepared in
    three different PDB files. They are loaded into separate :class:`Universe`
    objects (where they could be further manipulated, e.g. renumbered,
    relabeled, rotated, ...) The :func:`Merge` command is used to combine all
    of them together::

       u1 = Universe("protein.pdb")
       u2 = Universe("ligand.pdb")
       u3 = Universe("solvent.pdb")
       u = Merge(u1.select_atoms("protein"), u2.atoms, u3.atoms)
       u.atoms.write("system.pdb")

    The complete system is then written out to a new PDB file.


    .. versionchanged:: 0.9.0
       Raises exceptions instead of assertion errors.

    .. versionchanged:: 0.16.0
       The trajectory is now a
       :class:`~MDAnalysis.coordinates.memory.MemoryReader`.

    """
    from ..topology.base import squash_by

    if len(args) == 0:
        raise ValueError("Need at least one AtomGroup for merging")

    for a in args:
        if not isinstance(a, groups.AtomGroup):
            raise TypeError(repr(a) + " is not an AtomGroup")
    for a in args:
        if len(a) == 0:
            raise ValueError("cannot merge empty AtomGroup")

    # Create a new topology using the intersection of topology attributes
    blank_topology_attrs = set(dir(Topology(attrs=[])))
    common_attrs = set.intersection(*[set(dir(ag.universe._topology))
                                      for ag in args])
    tops = set(['bonds', 'angles', 'dihedrals', 'impropers'])

    attrs = []

    # Create set of attributes which are array-valued and can be simply
    # concatenated together
    common_array_attrs = common_attrs - blank_topology_attrs - tops
    # Build up array-valued topology attributes including only attributes
    # that all arguments' universes have
    for attrname in common_array_attrs:
        for ag in args:
            attr = getattr(ag, attrname)
            attr_class = type(getattr(ag.universe._topology, attrname))
            if issubclass(attr_class, AtomAttr):
                pass
            elif issubclass(attr_class, ResidueAttr):
                attr = getattr(ag.residues, attrname)
            elif issubclass(attr_class, SegmentAttr):
                attr = getattr(ag.segments, attrname)
            else:
                raise NotImplementedError("Don't know how to handle"
                                          " TopologyAttr not subclassed"
                                          " from AtomAttr, ResidueAttr,"
                                          " or SegmentAttr.")
            if type(attr) != np.ndarray:
                raise TypeError('Encountered unexpected topology '
                                'attribute of type {}'.format(type(attr)))
            try:
                attr_array.extend(attr)
            except NameError:
                attr_array = list(attr)
        attrs.append(attr_class(np.array(attr_array, dtype=attr.dtype)))
        del attr_array

    # Build up topology groups including only those that all arguments'
    # universes have
    for t in (tops & common_attrs):
        offset = 0
        bondidx = []
        types = []
        for ag in args:
            # create a mapping scheme for this atomgroup
            mapping = {a.index: i for i, a in enumerate(ag, start=offset)}
            offset += len(ag)

            tg = getattr(ag, t)
            bonds_class = type(getattr(ag.universe._topology, t))
            # Create a topology group of only bonds that are within this ag
            # ie we don't want bonds that extend out of the atomgroup
            tg = tg.atomgroup_intersection(ag, strict=True)

            # Map them so they refer to our new indices
            new_idx = [tuple([mapping[x] for x in entry]) for entry in tg.indices]
            bondidx.extend(new_idx)
            if hasattr(tg, '_bondtypes'):
                types.extend(tg._bondtypes)
            else:
                types.extend([None]*len(tg))
        if any(t is None for t in types):
            attrs.append(bonds_class(bondidx))
        else:
            types = np.array(types, dtype='|S8')
            attrs.append(bonds_class(bondidx, types))

    # Renumber residue and segment indices
    n_atoms = sum([len(ag) for ag in args])
    residx = []
    segidx = []
    res_offset = 0
    seg_offset = 0
    for ag in args:
        # create a mapping scheme for this atomgroup's parents
        res_mapping = {r.resindex: i for i, r in enumerate(ag.residues,
                                                           start=res_offset)}
        seg_mapping = {r.segindex: i for i, r in enumerate(ag.segments,
                                                           start=seg_offset)}
        res_offset += len(ag.residues)
        seg_offset += len(ag.segments)

        # Map them so they refer to our new indices
        residx.extend([res_mapping[x] for x in ag.resindices])
        segidx.extend([seg_mapping[x] for x in ag.segindices])

    residx = np.array(residx, dtype=np.int32)
    segidx = np.array(segidx, dtype=np.int32)

    _, _, [segidx] = squash_by(residx, segidx)

    n_residues = len(set(residx))
    n_segments = len(set(segidx))

    top = Topology(n_atoms, n_residues, n_segments,
                   attrs=attrs,
                   atom_resindex=residx,
                   residue_segindex=segidx)

    # Create and populate a universe
    try:
        #Create universe with coordinates if they exists in args
        coords = np.vstack([a.positions for a in args])
        u = Universe(top, coords[None, :, :],
                 format=MDAnalysis.coordinates.memory.MemoryReader)
    except AttributeError:
        #Create universe without coordinates if they dont exists in args
        u = Universe(top)

    return u
