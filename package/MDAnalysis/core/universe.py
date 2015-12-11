import numpy as np
from numpy.lib.utils import deprecate
import logging

import MDAnalysis
from ..lib import util
from ..lib.util import cached
from . import groups

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
    topologyfile
        A CHARMM/XPLOR PSF topology file, PDB file or Gromacs GRO file; used to
        define the list of atoms. If the file includes bond information,
        partial charges, atom masses, ... then these data will be available to
        MDAnalysis. A "structure" file (PSF, PDB or GRO, in the sense of a
        topology) is always required.
    permissive
        currently only relevant for PDB files: Set to ``True`` in order to
        ignore most errors and read typical MD simulation PDB files; set to
        ``False`` to read with the Bio.PDB reader, which can be useful for real
        Protein Databank PDB files. ``None``  selects the MDAnalysis default
        (which is set in :class:`MDAnalysis.core.flags`) [``None``]
    topology_format
        provide the file format of the topology file; ``None`` guesses it from
        the file extension [``None``] Can also pass a subclass of
        :class:`MDAnalysis.topology.base.TopologyReader` to define a custom
        reader to be used on the topology file.
    format
        provide the file format of the coordinate or trajectory file; ``None``
        guesses it from the file extension. Note that this keyword has no
        effect if a list of file names is supplied because the "chained" reader
        has to guess the file format for each individual list member.
        [``None``] Can also pass a subclass of
        :class:`MDAnalysis.coordinates.base.Reader` to define a custom reader
        to be used on the trajectory file.
    guess_bonds
        Once Universe has been loaded, attempt to guess the connectivity
        between atoms.  This will populate the .bonds .angles and .dihedrals
        attributes of the Universe.
    vdwradii
        For use with *guess_bonds*. Supply a dict giving a vdwradii for each
        atom type which are used in guessing bonds.
    is_anchor
        When unpickling instances of
        :class:`MDAnalysis.core.AtomGroup.AtomGroup` existing Universes are
        searched for one where to anchor those atoms. Set to ``False`` to
        prevent this Universe from being considered. [``True``]
    anchor_name
        Setting to other than ``None`` will cause
        :class:`MDAnalysis.core.AtomGroup.AtomGroup` instances pickled from the
        Universe to only unpickle if a compatible Universe with matching
        *anchor_name* is found. *is_anchor* will be ignored in this case but
        will still be honored when unpickling
        :class:`MDAnalysis.core.AtomGroup.AtomGroup` instances pickled with
        *anchor_name*==``None``. [``None``]

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
        from ..topology.core import get_parser_for
        from ..topology.base import TopologyReader
        from ..coordinates.base import ProtoReader

        # managed attribute holding Reader
        self._trajectory = None

        if len(args) == 0:
            # create an empty universe
            self._topology = None
            self.atoms = None
            return

        self.filename = args[0]
        coordinatefile = args[1:]

        topology_format = kwargs.pop('topology_format', None)

        if len(args) == 1:
            # special hacks to treat a coordinate file as a coordinate AND topology file
            if kwargs.get('format', None) is None:
                kwargs['format'] = topology_format
            elif topology_format is None:
                topology_format = kwargs.get('format', None)

            # if passed a Reader, use that
            fmt = kwargs.get('format', None)
            try:
                if issubclass(fmt, ProtoReader):
                    coordinatefile = self.filename
            except TypeError:
                # or if file is known as a topology & coordinate file, use that
                if fmt is None:
                    fmt = util.guess_format(self.filename)
                if (fmt in MDAnalysis.coordinates._trajectory_readers
                    and fmt in MDAnalysis.topology._topology_parsers):
                    coordinatefile = self.filename
            if len(coordinatefile) == 0:
                coordinatefile = None

        # build the topology (or at least a list of atoms)
        try:  # Try and check if the topology format is a TopologyReader
            if issubclass(topology_format, TopologyReader):
                parser = topology_format
        except TypeError:  # But strings/None raise TypeError in issubclass
            perm = kwargs.get('permissive',
                              MDAnalysis.core.flags['permissive_pdb_reader'])
            parser = get_parser_for(self.filename,
                                    permissive=perm,
                                    format=topology_format)
        try:
            with parser(self.filename, universe=self) as p:
                self._topology = p.parse()
        except IOError as err:
            raise IOError("Failed to load from the topology file {0}"
                          " with parser {1}.\n"
                          "Error: {2}".format(self.filename, parser, err))
        except ValueError as err:
            raise ValueError("Failed to construct topology from file {0}"
                             " with parser {1} \n"
                             "Error: {2}".format(self.filename, parser, err))

        # generate Group classes
        self._make_groupclasses()
        self._make_componentclasses()

        # Generate atoms, residues and segments
        self.atoms = self._groups['atom'](
                np.arange(self._topology.n_atoms), self)

        self.residues = self._groups['residue'](
                np.arange( self._topology.n_residues), self)

        self.segments = self._groups['segment'](np.arange(
            self._topology.n_segments), self)

        # Load coordinates
        self.load_new(coordinatefile, **kwargs)

    def _make_groupclasses(self):
        """Generates Group classes specific to this Universe based on its
        Topology.

        """
        self._groups = {}

        # generate Group class based on Topology
        self._Group = groups.make_group()
        for attr in self._topology.attrs:
            self._Group._add_prop(attr)

        # generate AtomGroup, ResidueGroup, and SegmentGroup classes for this
        # universe
        self._groups['atom'] = groups.make_levelgroup(self._Group,
                                                      level='atom')
        self._groups['residue'] = groups.make_levelgroup(self._Group,
                                                         level='residue')
        self._groups['segment'] = groups.make_levelgroup(self._Group,
                                                         level='segment')

    def _make_componentclasses(self):
        """Generates component classes specific to this Universe based on its
        Topology.

        """
        self._components = {}

        # for each level, generate Component class (Atom, Residue, Segment),
        # and attach attributes appropriate to that level from Topology
        for level in ['atom', 'residue', 'segment']:
            self._components[level] = groups.make_levelcomponent(level=level)
            for attr in self._topology.attrs:
                if attr.level == level:
                    self._components[level]._add_prop(attr)

    @property
    def universe(self):
        # for Writer.write(universe), see Issue 49
        # Encapsulation in an accessor prevents the Universe from having to keep a reference to itself,
        #  which might be undesirable if it has a __del__ method. It is also cleaner than a weakref.
        return self

    def load_new(self, filename, **kwargs):
        """Load coordinates from *filename*, using the suffix to detect file format.

        :Arguments:
             *filename*
                 the coordinate file (single frame or trajectory) *or* a list of
                 filenames, which are read one after another.
             *permissive*
                 currently only relevant for PDB files: Set to ``True`` in order to ignore most errors
                 and read typical MD simulation PDB files; set to ``False`` to read with the Bio.PDB reader,
                 which can be useful for real Protein Databank PDB files. ``None``  selects the
                 MDAnalysis default (which is set in :class:`MDAnalysis.core.flags`) [``None``]
             *format*
                 provide the file format of the coordinate or trajectory file;
                 ``None`` guesses it from the file extension. Note that this
                 keyword has no effect if a list of file names is supplied because
                 the "chained" reader has to guess the file format for each
                 individual list member [``None``]
                 Can also pass a subclass of :class:`MDAnalysis.coordinates.base.Reader`
                 to define a custom reader to be used on the trajectory file.
             *kwargs*
                 Other kwargs are passed to the trajectory reader (only for advanced use)

        :Returns: (filename, trajectory_format) or ``None`` if *filename* == ``None``
        :Raises: :exc:`TypeError` if trajectory format can not be
                  determined or no appropriate trajectory reader found

        .. versionchanged:: 0.8
           If a list or sequence that is provided for *filename*  only contains a single entry
           then it is treated as single coordinate file. This has the consequence that it is
           not read by the :class:`~MDAnalysis.coordinates.base.ChainReader` but directly by
           its specialized file format reader, which typically has more features than the
           :class:`~MDAnalysis.coordinates.base.ChainReader`.
        """
        if filename is None:
            return

        import MDAnalysis.core
        from ..coordinates.core import get_reader_for
        from ..coordinates.base import ProtoReader

        if len(util.asiterable(filename)) == 1:
            # make sure a single filename is not handed to the ChainReader
            filename = util.asiterable(filename)[0]
        logger.debug("Universe.load_new(): loading {0}...".format(filename))

        reader_format = kwargs.pop('format', None)
        perm = kwargs.get('permissive', MDAnalysis.core.flags['permissive_pdb_reader'])
        reader = None

        # Check if we were passed a Reader to use
        try:
            if reader_format is not None and issubclass(reader_format, ProtoReader):
                reader = reader_format
        except TypeError:
            pass

        if not reader:
            # Check if we need to use Chain reader
            if util.iterable(filename):
                # Save the format and pass this to ChainReader
                kwargs.update({'format': reader_format})
                reader_format='CHAIN'
            try:
                reader = get_reader_for(filename,
                                        permissive=perm,
                                        format=reader_format)
            except TypeError as err:
                raise TypeError(
                    "Cannot find an appropriate coordinate reader for file '{0}'.\n"
                    "           {1}".format(filename, err))
        # supply number of atoms for readers that cannot do it for themselves
        kwargs['n_atoms'] = len(self.atoms)

        self.trajectory = reader(filename, **kwargs)    # unified trajectory API
        if self.trajectory.n_atoms != len(self.atoms):
            raise ValueError("The topology and {form} trajectory files don't"
                             " have the same number of atoms!\n"
                             "Topology number of atoms {top_n_atoms}\n"
                             "Trajectory: {fname} Number of atoms {trj_n_atoms}".format(
                                 form=self.trajectory.format,
                                 top_n_atoms=len(self.atoms),
                                 fname=filename,
                                 trj_n_atoms=self.trajectory.n_atoms))

        return filename, self.trajectory.format

    def select_atoms(self, sel, *othersel, **selgroups):
        return self.atoms.select_atoms(sel, *othersel, **selgroups)

    def __repr__(self):
        #return "<Universe with {n_atoms} atoms{bonds}>".format(
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

           In order to access the coordinates it is probably better to use the
           :meth:`AtomGroup.coordinates` method; for instance, all coordinates
           of the Universe as a numpy array:
           :meth:`Universe.atoms.coordinates`.
        """
        return self.trajectory.ts

    @property
    def trajectory(self):
        """Reference to trajectory reader object containing trajectory data."""
        if not self._trajectory is None:
            return self._trajectory
        else:
            raise AttributeError("No trajectory loaded into Universe")

    @trajectory.setter
    def trajectory(self, value):
        del self._trajectory  # guarantees that files are closed (?)
        self._trajectory = value

    def add_TopologyAttr(self, topologyattr):
        self._topology.add_TopologyAttr(topologyattr)
        self._Group._add_prop(topologyattr)


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

    :Returns: an instance of :class:`~MDAnalaysis.AtomGroup.Universe`
    """
    if len(args) == 0:
        raise TypeError("as_Universe() takes at least one argument (%d given)" % len(args))
    elif len(args) == 1 and isinstance(args[0], Universe):
        return args[0]
    return Universe(*args, **kwargs)

asUniverse = deprecate(as_Universe, old_name='asUniverse', new_name='as_Universe')

#TODO: UPDATE ME WITH NEW TOPOLOGY DETAILS
def Merge(*args):
    """Return a :class:`Universe` from two or more :class:`AtomGroup` instances.

    :class:`AtomGroup` instances can come from different Universes, or come
    directly from a :meth:`~Universe.select_atoms` call.

    It can also be used with a single :class:`AtomGroup` if the user wants to,
    for example, re-order the atoms in the Universe.

    :Arguments: One or more :class:`AtomGroup` instances.

    :Returns: an instance of :class:`~MDAnalaysis.AtomGroup.Universe`

    :Raises: :exc:`ValueError` for too few arguments or if an AtomGroup is
             empty and :exc:`TypeError` if arguments are not
             :class:`AtomGroup` instances.

    .. rubric:: Example

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

    .. Note:: Merging does not create a full trajectory but only a single
              structure even if the input consists of one or more trajectories.

    .. versionchanged 0.9.0::
       Raises exceptions instead of assertion errors.

    """
    import MDAnalysis.topology.core

    if len(args) == 0:
        raise ValueError("Need at least one AtomGroup for merging")

    for a in args:
        if not isinstance(a, AtomGroup):
            raise TypeError(repr(a) + " is not an AtomGroup")
    for a in args:
        if len(a) == 0:
            raise ValueError("cannot merge empty AtomGroup")

    coords = np.vstack([a.coordinates() for a in args])
    trajectory = MDAnalysis.coordinates.base.Reader(None)
    ts = MDAnalysis.coordinates.base.Timestep.from_coordinates(coords)
    setattr(trajectory, "ts", ts)
    trajectory.n_frames = 1

    # create an empty Universe object
    u = Universe()
    u.trajectory = trajectory

    # create a list of Atoms, then convert it to an AtomGroup
    atoms = [copy.copy(a) for gr in args for a in gr]
    for a in atoms:
        a.universe = u

    # adjust the atom numbering
    for i, a in enumerate(atoms):
        a.index = i
        a.serial = i + 1
    u.atoms = AtomGroup(atoms)

    # move over the topology
    offset = 0
    tops = ['bonds', 'angles', 'dihedrals', 'impropers']
    idx_lists = {t:[] for t in tops}
    for ag in args:
        # create a mapping scheme for this atomgroup
        mapping = {a.index:i for i, a in enumerate(ag, start=offset)}
        offset += len(ag)

        for t in tops:
            tg = getattr(ag, t)
            # Create a topology group of only bonds that are within this ag
            # ie we don't want bonds that extend out of the atomgroup
            tg = tg.atomgroup_intersection(ag, strict=True)

            # Map them so they refer to our new indices
            new_idx = [tuple(map(lambda x:mapping[x], entry))
                       for entry in tg.to_indices()]
            idx_lists[t].extend(new_idx)

    for t in tops:
        u._topology[t] = idx_lists[t]

    # adjust the residue and segment numbering (removes any remaining references to the old universe)
    MDAnalysis.topology.core.build_residues(u.atoms)
    MDAnalysis.topology.core.build_segments(u.atoms)

    return u
