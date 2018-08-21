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


Working with Universes
======================


Quick segid selection
---------------------

.. deprecated:: 0.16.2
   Instant selectors will be removed in the 1.0 release.  See issue `#1377
   <https://github.com/MDAnalysis/mdanalysis/issues/1377>`_ for more details.


If the loaded topology provided segids, then these are made accessible
as attributes of the Universe.  If the segid starts with a number such
as '4AKE', the letter 's' will be prepended to the segid.
For example::

   import MDAnalysis as mda
   from MDAnalysisTests.datafiles import PSF, DCD

   u = mda.Universe(PSF, DCD)
   u.select_atoms('segid 4AKE')  # selects all segments with segid 4AKE

If only a single segment has that segid then a Segment object will
be returned, otherwise a SegmentGroup will be returned.


Classes
=======

.. autoclass:: Universe
   :members:

Functions
=========

.. autofunction:: Merge

"""
from __future__ import absolute_import
from six.moves import range
import six

import errno
import numpy as np
import logging
import copy
import warnings

import MDAnalysis
import sys

# When used in an MPI environment with Infiniband, importing MDAnalysis may
# trigger an MPI warning because importing the uuid module triggers a call to
# os.fork(). This happens if MPI_Init() has been called prior to importing
# MDAnalysis. The problem is actually caused by the uuid module and not by
# MDAnalysis itself. Python 3.7 fixes the problem. However, for Python < 3.7,
# the uuid module works perfectly fine with os.fork() disabled during import.
# A clean solution is therefore simply to disable os.fork() prior to importing
# the uuid module and to re-enable it afterwards.
import os

# Windows doesn't have os.fork so can ignore
# the issue for that platform
if sys.version_info >= (3, 7) or os.name == 'nt':
    import uuid
else:
    _os_dot_fork, os.fork = os.fork, None
    import uuid
    os.fork = _os_dot_fork
    del _os_dot_fork

from .. import _ANCHOR_UNIVERSES, _TOPOLOGY_ATTRS, _PARSERS
from ..exceptions import NoDataError
from ..lib import util
from ..lib.log import ProgressMeter, _set_verbose
from ..lib.util import cached, NamedStream, isstream
from ..lib.mdamath import find_fragments
from . import groups
from ._get_readers import get_reader_for, get_parser_for
from .groups import (ComponentBase, GroupBase,
                     Atom, Residue, Segment,
                     AtomGroup, ResidueGroup, SegmentGroup)
from .topology import Topology
from .topologyattrs import AtomAttr, ResidueAttr, SegmentAttr

logger = logging.getLogger("MDAnalysis.core.universe")


class Universe(object):
    """The MDAnalysis Universe contains all the information describing the system.

    The system always requires a *topology* file --- in the simplest case just
    a list of atoms. This can be a CHARMM/NAMD PSF file or a simple coordinate
    file with atom informations such as XYZ, PDB, Gromacs GRO, or CHARMM
    CRD. See :ref:`Supported topology formats` for what kind of topologies can
    be read.

    A trajectory provides coordinates; the coordinates have to be ordered in
    the same way as the list of atoms in the topology. A trajectory can be a
    single frame such as a PDB, CRD, or GRO file, or it can be a MD trajectory
    (in CHARMM/NAMD/LAMMPS DCD, Gromacs XTC/TRR, or generic XYZ format).  See
    :ref:`Supported coordinate formats` for what can be read as a
    "trajectory".

    As a special case, when the topology is a file that contains atom
    information *and* coordinates (such as XYZ, PDB, GRO or CRD, see
    :ref:`Supported coordinate formats`) then the coordinates are immediately
    loaded from the "topology" file unless a trajectory is supplied.

    Examples for setting up a universe::

       u = Universe(topology, trajectory)          # read system from file(s)
       u = Universe(pdbfile)                       # read atoms and coordinates from PDB or GRO
       u = Universe(topology, [traj1, traj2, ...]) # read from a list of trajectories
       u = Universe(topology, traj1, traj2, ...)   # read from multiple trajectories

    Load new data into a universe (replaces old trajectory and does *not* append)::

       u.load_new(trajectory)                      # read from a new trajectory file

    Select atoms, with syntax similar to CHARMM (see
    :class:`~Universe.select_atoms` for details)::

       u.select_atoms(...)

    Parameters
    ----------
    topology : str, Topology object or stream
        A CHARMM/XPLOR PSF topology file, PDB file or Gromacs GRO file; used to
        define the list of atoms. If the file includes bond information,
        partial charges, atom masses, ... then these data will be available to
        MDAnalysis. A "structure" file (PSF, PDB or GRO, in the sense of a
        topology) is always required. Alternatively, an existing
        :class:`MDAnalysis.core.topology.Topology` instance may also be given.
    topology_format
        Provide the file format of the topology file; ``None`` guesses it from
        the file extension [``None``] Can also pass a subclass of
        :class:`MDAnalysis.topology.base.TopologyReaderBase` to define a custom
        reader to be used on the topology file.
    format
        Provide the file format of the coordinate or trajectory file; ``None``
        guesses it from the file extension. Note that this keyword has no
        effect if a list of file names is supplied because the "chained" reader
        has to guess the file format for each individual list member.
        [``None``] Can also pass a subclass of
        :class:`MDAnalysis.coordinates.base.ProtoReader` to define a custom
        reader to be used on the trajectory file.
    all_coordinates : bool
        If set to ``True`` specifies that if more than one filename is passed
        they are all to be used, if possible, as coordinate files (employing a
        :class:`MDAnalysis.coordinates.chain.ChainReader`). [``False``] The
        default behavior is to take the first file as a topology and the
        remaining as coordinates. The first argument will always always be used
        to infer a topology regardless of *all_coordinates*. This parameter is
        ignored if only one argument is passed.
    guess_bonds : bool, optional
        Once Universe has been loaded, attempt to guess the connectivity
        between atoms.  This will populate the .bonds .angles and .dihedrals
        attributes of the Universe.
    vdwradii : dict, optional
        For use with *guess_bonds*. Supply a dict giving a vdwradii for each
        atom type which are used in guessing bonds.
    is_anchor : bool, optional
        When unpickling instances of
        :class:`MDAnalysis.core.groups.AtomGroup` existing Universes are
        searched for one where to anchor those atoms. Set to ``False`` to
        prevent this Universe from being considered. [``True``]
    anchor_name : str, optional
        Setting to other than ``None`` will cause
        :class:`MDAnalysis.core.groups.AtomGroup` instances pickled from the
        Universe to only unpickle if a compatible Universe with matching
        *anchor_name* is found. Even if *anchor_name* is set *is_anchor* will
        still be honored when unpickling.
    transformations: function or list, optional
        Provide a list of transformations that you wish to apply to the 
        trajectory upon reading. Transformations can be found in 
        :mod:`MDAnalysis.transformations`, or can be user-created.
    in_memory
        After reading in the trajectory, transfer it to an in-memory
        representations, which allow for manipulation of coordinates.
    in_memory_step
        Only read every nth frame into in-memory representation.
    continuous : bool, optional
        The `continuous` option is used by the
        :mod:`ChainReader<MDAnalysis.coordinates.chain>`, which contains the
        functionality to treat independent trajectory files as a single virtual
        trajectory.

    Attributes
    ----------
    trajectory
        currently loaded trajectory reader;
    dimensions
        current system dimensions (simulation unit cell, if set in the
        trajectory)
    atoms, residues, segments
        master Groups for each topology level
    bonds, angles, dihedrals
        master ConnectivityGroups for each connectivity type

    """

    def __init__(self, *args, **kwargs):
        # Store the segments for the deprecated instant selector feature.
        # This attribute has to be defined early to avoid recursion in
        # __getattr__.
        self._instant_selectors = {}
        # hold on to copy of kwargs; used by external libraries that
        # reinitialize universes
        self._kwargs = copy.deepcopy(kwargs)

        # managed attribute holding Reader
        self._trajectory = None
        self._cache = {}

        if not args:
            # create an empty universe
            self._topology = None
            self.atoms = None
        else:
            topology_format = kwargs.pop('topology_format', None)
            if len(args) == 1:
                # special hacks to treat a coordinate file as a coordinate AND
                # topology file
                if kwargs.get('format', None) is None:
                    kwargs['format'] = topology_format
                elif topology_format is None:
                    topology_format = kwargs.get('format', None)

            # if we're given a Topology object, we don't need to parse anything
            if isinstance(args[0], Topology):
                self._topology = args[0]
                self.filename = None
            else:
                if isinstance(args[0], NamedStream):
                    self.filename = args[0]
                elif isstream(args[0]):
                    filename = None
                    if hasattr(args[0], 'name'):
                        filename = args[0].name
                    self.filename = NamedStream(args[0], filename)
                else:
                    self.filename = args[0]
                parser = get_parser_for(self.filename, format=topology_format)
                try:
                    with parser(self.filename) as p:
                        self._topology = p.parse(**kwargs)
                except (IOError, OSError) as err:
                    # There are 2 kinds of errors that might be raised here:
                    # one because the file isn't present
                    # or the permissions are bad, second when the parser fails
                    if (err.errno is not None and
                        errno.errorcode[err.errno] in ['ENOENT', 'EACCES']):
                        # Runs if the error is propagated due to no permission / file not found
                        six.reraise(*sys.exc_info())
                    else:
                        # Runs when the parser fails
                        raise IOError(
                            "Failed to load from the topology file {0}"
                            " with parser {1}.\n"
                            "Error: {2}".format(self.filename, parser, err))
                except (ValueError, NotImplementedError) as err:
                    raise ValueError(
                        "Failed to construct topology from file {0}"
                        " with parser {1}.\n"
                        "Error: {2}".format(self.filename, parser, err))

            # generate and populate Universe version of each class
            self._generate_from_topology()

            # Load coordinates
            if len(args) == 1 or kwargs.get('all_coordinates', False):
                if self.filename is None:
                    # If we got the topology as a Topology object, then we
                    # cannot read coordinates from it.
                    coordinatefile = args[1:]
                else:
                    # Can the topology file also act as coordinate file?
                    try:
                        _ = get_reader_for(self.filename,
                                           format=kwargs.get('format', None))
                    except ValueError:
                        coordinatefile = args[1:]
                    else:
                        coordinatefile = (self.filename,) + args[1:]
            else:
                coordinatefile = args[1:]

            if not coordinatefile:
                coordinatefile = None
                
            self.load_new(coordinatefile, **kwargs)
            # parse transformations
            trans_arg = kwargs.pop('transformations', None)
            if trans_arg:
                transforms =[trans_arg] if callable(trans_arg) else trans_arg
                self.trajectory.add_transformations(*transforms)

        # Check for guess_bonds
        if kwargs.pop('guess_bonds', False):
            self.atoms.guess_bonds(vdwradii=kwargs.pop('vdwradii', None))

        # None causes generic hash to get used.
        # We store the name ieven if is_anchor is False in case the user later
        #  wants to make the universe an anchor.
        self._anchor_name = kwargs.get('anchor_name', None)
        # Universes are anchors by default
        self.is_anchor = kwargs.get('is_anchor', True)
        

    def copy(self):
        """Return an independent copy of this Universe"""
        new = self.__class__(self._topology.copy())
        new.trajectory = self.trajectory.copy()
        return new

    def _generate_from_topology(self):
        # generate Universe version of each class
        # AG, RG, SG, A, R, S
        self._class_bases, self._classes = groups.make_classes()

        # Put Group level stuff from topology into class
        for attr in self._topology.attrs:
            self._process_attr(attr)

        # Generate atoms, residues and segments.
        # These are the first such groups generated for this universe, so
        #  there are no cached merged classes yet. Otherwise those could be
        #  used directly to get a (very) small speedup. (Only really pays off
        #  the readability loss if instantiating millions of AtomGroups at
        #  once.)
        self.atoms = AtomGroup(np.arange(self._topology.n_atoms), self)

        self.residues = ResidueGroup(
                np.arange(self._topology.n_residues), self)

        self.segments = SegmentGroup(
                np.arange(self._topology.n_segments), self)

        # Update Universe namespace with segids
        # Many segments can have same segid, so group together first
        #
        # DEPRECATED in 0.16.2
        # REMOVE in 1.0
        # See https://github.com/MDAnalysis/mdanalysis/issues/1377
        try:
            # returns dict of segid:segment
            segids = self.segments.groupby('segids')
        except AttributeError:
            # no segids, don't do this step
            pass
        else:
            for segid, segment in segids.items():
                if not segid:  # ignore blank segids
                    continue

                # cannot start attribute with number
                if segid[0].isdigit():
                    # prefix 's' if starts with number
                    name = 's' + segid
                else:
                    name = segid
                # if len 1 SegmentGroup, convert to Segment
                if len(segment) == 1:
                    segment = segment[0]
                self._instant_selectors[name] = segment

    @classmethod
    def empty(cls, n_atoms, n_residues=None, n_segments=None,
              atom_resindex=None, residue_segindex=None,
              trajectory=False, velocities=False, forces=False):
        """Create a blank Universe

        Useful for building a Universe without requiring existing files,
        for example for system building.

        Parameters
        ----------
        n_atoms : int
          number of Atoms in the Universe
        n_residues : int, optional
          number of Residues in the Universe, defaults to 1
        n_segments : int, optional
          number of Segments in the Universe, defaults to 1
        atom_resindex : numpy.array, optional
          mapping of atoms to residues
        residue_segindex : numpy.array, optional
          mapping of residues to segments
        trajectory : bool, optional
          if True, attaches a dummy reader to the Universe, therefore
          allowing coordinates to be set and written.  Default is False
        velocities : bool, optional
          include velocities in the dummy Reader
        forces : bool, optional
          include forces in the dummy Reader

        Returns
        -------
        MDAnalysis.Universe object

        Examples
        --------
        For example to create a new Universe with 6 atoms in 2 residues, with
        positions for the atoms and a mass attribute:

        >>> u = mda.Universe.empty(6, 2,
                                   atom_resindex=np.array([0, 0, 0, 1, 1, 1]),
                                   trajectory=True,
                )
        >>> u.add_TopologyAttr('masses')

        .. versionadded:: 0.17.0
        """
        if n_residues is None:
            n_residues = 1
        elif atom_resindex is None:
            warnings.warn(
                'Multiple residues specified but no atom_resindex given.  '
                'All atoms will be placed in first Residue.',
                UserWarning)
        if n_segments is None:
            n_segments = 1
        elif residue_segindex is None:
            warnings.warn(
                'Multiple segments specified but no segment_resindex given.  '
                'All residues will be placed in first Segment',
                UserWarning)

        top = Topology(n_atoms, n_residues, n_segments,
                       atom_resindex=atom_resindex,
                       residue_segindex=residue_segindex,
        )

        u = cls(top)

        if trajectory:
            u.trajectory = get_reader_for('', format='dummy')(
                n_atoms=n_atoms,
                velocities=velocities, forces=forces)

        return u

    def __getattr__(self, key):
        # This implements the instant selector of segments from a Universe.
        # It is implemented as __getattr__ so a deprecation warning can be
        # issued when the feature is used. Instant selectors are deprecated
        # since version 0.16.2 and are tareted to be deleted in version 1.0.
        # self._instant_selectors is populated in self._process_attr and
        # created at the beginning of __init__.
        try:
            segment = self._instant_selectors[key]
        except KeyError:
            raise AttributeError('No attribute "{}".'.format(key))
        else:
            warnings.warn("Instant selector Universe.<segid> "
                          "is deprecated and will be removed in 1.0. "
                          "Use SegmentGroup[SegmentGroup.segids == '<segid>'] "
                          "instead.",
                          DeprecationWarning)
            return segment

    @property
    def universe(self):
        # for Writer.write(universe), see Issue 49
        # Encapsulation in an accessor prevents the Universe from
        # having to keep a reference to itself,
        #  which might be undesirable if it has a __del__ method.
        # It is also cleaner than a weakref.
        return self

    def load_new(self, filename, format=None, in_memory=False, **kwargs):
        """Load coordinates from `filename`.

        The file format of `filename` is autodetected from the file name suffix
        or can be explicitly set with the `format` keyword. A sequence of files
        can be read as a single virtual trajectory by providing a list of
        filenames.


        Parameters
        ----------
        filename : str or list
            the coordinate file (single frame or trajectory) *or* a list of
            filenames, which are read one after another.
        format : str or list or object (optional)
            provide the file format of the coordinate or trajectory file;
            ``None`` guesses it from the file extension. Note that this
            keyword has no effect if a list of file names is supplied because
            the "chained" reader has to guess the file format for each
            individual list member [``None``]. Can also pass a subclass of
            :class:`MDAnalysis.coordinates.base.ProtoReader` to define a custom
            reader to be used on the trajectory file.
        in_memory : bool (optional)
            Directly load trajectory into memory with the
            :class:`~MDAnalysis.coordinates.memory.MemoryReader`

            .. versionadded:: 0.16.0

        **kwargs : dict
            Other kwargs are passed to the trajectory reader (only for
            advanced use)

        Returns
        -------
        universe : Universe

        Raises
        ------
        TypeError if trajectory format can not be
                  determined or no appropriate trajectory reader found


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

        self.trajectory = reader(filename, **kwargs)
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
            self.transfer_to_memory(step=kwargs.get("in_memory_step", 1))

        return self

    def transfer_to_memory(self, start=None, stop=None, step=None,
                           verbose=None, quiet=None):
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
        step : int, optional
            Read in every nth frame. [1]
        verbose : bool, optional
            Will print the progress of loading trajectory to memory, if
            set to True. Default value is False.


        .. versionadded:: 0.16.0
        """
        from ..coordinates.memory import MemoryReader

        verbose = _set_verbose(verbose, quiet, default=False)

        if not isinstance(self.trajectory, MemoryReader):
            # Try to extract coordinates using Timeseries object
            # This is significantly faster, but only implemented for certain
            # trajectory file formats
            try:
                coordinates = self.trajectory.timeseries(
                    self.atoms, start=start, stop=stop, step=step, order='fac')
            # if the Timeseries extraction fails,
            # fall back to a slower approach
            except AttributeError:
                n_frames = len(range(
                    *self.trajectory.check_slice_indices(start, stop, step)
                ))
                pm_format = '{step}/{numsteps} frames copied to memory (frame {frame})'
                pm = ProgressMeter(n_frames, interval=1,
                                   verbose=verbose, format=pm_format)
                coordinates = []  # TODO: use pre-allocated array
                for i, ts in enumerate(self.trajectory[start:stop:step]):
                    coordinates.append(np.copy(ts.positions))
                    pm.echo(i, frame=ts.frame)
                coordinates = np.array(coordinates)

            # Overwrite trajectory in universe with an MemoryReader
            # object, to provide fast access and allow coordinates
            # to be manipulated
            if step is None:
                step = 1
            self.trajectory = MemoryReader(
                coordinates,
                dimensions=self.trajectory.ts.dimensions,
                dt=self.trajectory.ts.dt * step,
                filename=self.trajectory.filename)

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
        """Bonds between atoms"""
        return self.atoms.bonds

    @property
    def angles(self):
        """Angles between atoms"""
        return self.atoms.angles

    @property
    def dihedrals(self):
        """Dihedral angles between atoms"""
        return self.atoms.dihedrals

    @property
    def impropers(self):
        """Improper dihedral angles between atoms"""
        return self.atoms.impropers

    @property
    def anchor_name(self):
        return self._gen_anchor_hash()

    @anchor_name.setter
    def anchor_name(self, name):
        self.remove_anchor()  # clear any old anchor
        self._anchor_name = str(name) if not name is None else name
        self.make_anchor()  # add anchor again

    def _gen_anchor_hash(self):
        # hash used for anchoring.
        # Try and use anchor_name, else use (and store) uuid
        if self._anchor_name is not None:
            return self._anchor_name
        else:
            try:
                return self._anchor_uuid
            except AttributeError:
                # store this so we can later recall it if needed
                self._anchor_uuid = uuid.uuid4()
                return self._anchor_uuid

    @property
    def is_anchor(self):
        """Is this Universe an anchoring for unpickling AtomGroups"""
        return self._gen_anchor_hash() in _ANCHOR_UNIVERSES

    @is_anchor.setter
    def is_anchor(self, new):
        if new:
            self.make_anchor()
        else:
            self.remove_anchor()

    def remove_anchor(self):
        """Remove this Universe from the possible anchor list for unpickling"""
        _ANCHOR_UNIVERSES.pop(self._gen_anchor_hash(), None)

    def make_anchor(self):
        _ANCHOR_UNIVERSES[self._gen_anchor_hash()] = self

    def __repr__(self):
        # return "<Universe with {n_atoms} atoms{bonds}>".format(
        #    n_atoms=len(self.atoms),
        #    bonds=" and {0} bonds".format(len(self.bonds)) if self.bonds else "")

        return "<Universe with {n_atoms} atoms>".format(
            n_atoms=len(self.atoms))

    def __getstate__(self):
        raise NotImplementedError

    def __setstate__(self, state):
        raise NotImplementedError

    # Properties
    @property
    def dimensions(self):
        """Current dimensions of the unitcell"""
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
        :class:`~MDAnalysis.coordinates.base.Timestep`, it changes its contents
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
        """Reference to trajectory reader object containing trajectory data."""
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
        topologyattr : TopologyAttr or string
          Either a MDAnalysis TopologyAttr object or the name of a possible
          topology attribute.
        values : np.ndarray, optional
          If initiating an attribute from a string, the initial values to
          use.  If not supplied, the new TopologyAttribute will have empty
          or zero values.

        Example
        -------
        For example to add bfactors to a Universe:

        >>> u.add_TopologyAttr('bfactors')
        >>> u.atoms.bfactors
        array([ 0.,  0.,  0., ...,  0.,  0.,  0.])

        .. versionchanged:: 0.17.0
           Can now also add TopologyAttrs with a string of the name of the
           attribute to add (eg 'charges'), can also supply initial values
           using values keyword.
        """
        if isinstance(topologyattr, six.string_types):
            try:
                tcls = _TOPOLOGY_ATTRS[topologyattr]
            except KeyError:
                raise ValueError(
                    "Unrecognised topology attribute name: '{}'."
                    "  Possible values: '{}'\n"
                    "To raise an issue go to: http://issues.mdanalysis.org"
                    "".format(
                        topologyattr, ', '.join(sorted(_TOPOLOGY_ATTRS.keys())))
                )
            else:
                topologyattr = tcls.from_blank(
                    n_atoms=self._topology.n_atoms,
                    n_residues=self._topology.n_residues,
                    n_segments=self._topology.n_segments,
                    values=values)
        self._topology.add_TopologyAttr(topologyattr)
        self._process_attr(topologyattr)

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

    def add_Residue(self, segment=None, **attrs):
        """Add a new Residue to this Universe

        New Residues will not contain any Atoms, but can be assigned to Atoms
        as per usual.  If the Universe contains multiple segments, this must
        be specified as a keyword.

        Parameters
        ----------
        segment : MDAnalysis.Segment
          If there are multiple segments, then the Segment that the new
          Residue will belong in must be specified.
        attrs : dict
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

        >>> newres = u.add_Residue(segment=u.segments[0], resid=42, resname='GLY')
        >>> u.atoms[[1, 2, 3]].residues = newres
        >>> u.select_atoms('resname GLY and resid 42')
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
        attrs : dict
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
        """
        atoms = self.atoms.ix
        bonds = self.atoms.bonds.to_indices()

        frag_indices = find_fragments(atoms, bonds)
        frags = tuple([AtomGroup(np.sort(ix), self) for ix in frag_indices])

        fragdict = {}
        for f in frags:
            for a in f:
                fragdict[a] = f

        return fragdict


# TODO: what is the point of this function???
def as_Universe(*args, **kwargs):
    """Return a universe from the input arguments.

    1. If the first argument is a universe, just return it::

         as_Universe(universe) --> universe

    2. Otherwise try to build a universe from the first or the first
       and second argument::

         as_Universe(PDB, **kwargs) --> Universe(PDB, **kwargs)
         as_Universe(PSF, DCD, **kwargs) --> Universe(PSF, DCD, **kwargs)
         as_Universe(*args, **kwargs) --> Universe(*args, **kwargs)

    Returns
    -------
    :class:`~MDAnalysis.core.groups.Universe`
    """
    if len(args) == 0:
        raise TypeError("as_Universe() takes at least one argument (%d given)" % len(args))
    elif len(args) == 1 and isinstance(args[0], Universe):
        return args[0]
    return Universe(*args, **kwargs)


def Merge(*args):
    """Create a new new :class:`Universe` from one or more
    :class:`~MDAnalysis.core.groups.AtomGroup` instances.

    Parameters
    ----------
    *args : :class:`~MDAnalysis.core.groups.AtomGroup`
        One or more AtomGroups.

    Returns
    -------
    universe : :class:`Universe`

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
    coords = np.vstack([a.positions for a in args])
    u = Universe(top, coords[None, :, :],
                 format=MDAnalysis.coordinates.memory.MemoryReader)

    return u
