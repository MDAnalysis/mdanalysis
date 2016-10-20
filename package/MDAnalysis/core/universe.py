import numpy as np
from numpy.lib.utils import deprecate
import logging
import copy

import MDAnalysis
from ..lib import util
from ..lib.util import cached
from ..exceptions import NoDataError
from . import groups
from .groups import (GroupBase, Atom, Residue, Segment,
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
    topology : filename or Topology object
        A CHARMM/XPLOR PSF topology file, PDB file or Gromacs GRO file; used to
        define the list of atoms. If the file includes bond information,
        partial charges, atom masses, ... then these data will be available to
        MDAnalysis. A "structure" file (PSF, PDB or GRO, in the sense of a
        topology) is always required. Alternatively, an existing
        :class:`MDAnalysis.core.topology.Topology` instance may also be given.
    topology_format
        Provide the file format of the topology file; ``None`` guesses it from
        the file extension [``None``] Can also pass a subclass of
        :class:`MDAnalysis.topology.base.TopologyReader` to define a custom
        reader to be used on the topology file.
    format
        Provide the file format of the coordinate or trajectory file; ``None``
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
    in_memory
        After reading in the trajectory, transfer it to an in-memory
        representations, which allow for manipulation of coordinates.
    in_memory_frame_interval
        Only read every nth frame into in-memory representation.

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

        # hold on to copy of kwargs; used by external libraries that
        # reinitialize universes
        self.kwargs = copy.deepcopy(kwargs)

        # managed attribute holding Reader
        self._trajectory = None
        self._cache = {}

        if len(args) == 0:
            # create an empty universe
            self._topology = None
            self.atoms = None
            return

        coordinatefile = args[1:]
        if not coordinatefile:
            coordinatefile = None

        # if we're given a Topology object, we don't need to parse anything
        if isinstance(args[0], Topology):
            self._topology = args[0]
            self.filename = None
        else:
            self.filename = args[0]
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
                    if (fmt in MDAnalysis.coordinates._READERS
                        and fmt in MDAnalysis.topology._PARSERS):
                        coordinatefile = self.filename

            # build the topology (or at least a list of atoms)
            try:  # Try and check if the topology format is a TopologyReader
                if issubclass(topology_format, TopologyReader):
                    parser = topology_format
            except TypeError:  # But strings/None raise TypeError in issubclass
                parser = get_parser_for(self.filename,
                                        format=topology_format)
            try:
                with parser(self.filename) as p:
                    self._topology = p.parse()
            except IOError as err:
                raise IOError("Failed to load from the topology file {0}"
                              " with parser {1}.\n"
                              "Error: {2}".format(self.filename, parser, err))
            except ValueError as err:
                raise ValueError("Failed to construct topology from file {0}"
                                 " with parser {1} \n"
                                 "Error: {2}".format(self.filename, parser, err))

        # generate and populate Universe version of each class
        self._generate_from_topology()

        # Load coordinates
        self.load_new(coordinatefile, **kwargs)

        # Check for guess_bonds
        if kwargs.pop('guess_bonds', False):
            self.atoms.guess_bonds(vdwradii=kwargs.pop('vdwradii', None))

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
        for seg in self.segments:
            if hasattr(seg, 'segid') and seg.segid:
                if seg.segid[0].isdigit():
                    name = 's' + seg.segid
                else:
                    name = seg.segid
                self.__dict__[name] = seg

    @property
    def universe(self):
        # for Writer.write(universe), see Issue 49
        # Encapsulation in an accessor prevents the Universe from
        # having to keep a reference to itself,
        #  which might be undesirable if it has a __del__ method.
        # It is also cleaner than a weakref.
        return self

    def load_new(self, filename, format=None, in_memory=False, **kwargs):
        """Load coordinates from *filename*, using the suffix to detect file
        format.

        Parameters
        ----------
        filename : str
            the coordinate file (single frame or trajectory) *or* a list of
            filenames, which are read one after another.
        format : str (optional)
            provide the file format of the coordinate or trajectory file;
            ``None`` guesses it from the file extension. Note that this
            keyword has no effect if a list of file names is supplied because
            the "chained" reader has to guess the file format for each
            individual list member [``None``] Can also pass a subclass of
            :class:`MDAnalysis.coordinates.base.Reader` to define a custom
            reader to be used on the trajectory file.
        in_memory : bool (optional)
            Directly load trajectory into memory
        kwargs : dict
            Other kwargs are passed to the trajectory reader (only for
            advanced use)


        # TODO: check what happens if filename is None
        Returns
        -------
        filename : str
        trajectory_format : str

        # TODO: look up raises doc formating
        Raises
        ------
        TypeError if trajectory format can not be
                  determined or no appropriate trajectory reader found

        .. versionchanged:: 0.8
           If a list or sequence that is provided for *filename* only contains
           a single entry then it is treated as single coordinate file. This
           has the consequence that it is not read by the
           :class:`~MDAnalysis.coordinates.base.ChainReader` but directly by
           its specialized file format reader, which typically has more
           features than the
           :class:`~MDAnalysis.coordinates.base.ChainReader`.

        """
        # TODO: is this really sensible? Why not require a filename arg?
        if filename is None:
            return

        from ..coordinates.core import get_reader_for
        from ..coordinates.base import ProtoReader

        if len(util.asiterable(filename)) == 1:
            # make sure a single filename is not handed to the ChainReader
            filename = util.asiterable(filename)[0]
        logger.debug("Universe.load_new(): loading {0}...".format(filename))

        reader_format = format
        reader = None

        # Check if we were passed a Reader to use
        try:
            if reader_format is not None and issubclass(reader_format, ProtoReader):
                reader = reader_format
        except TypeError:
            pass

        if reader is None:
            # Check if we need to use Chain reader
            if util.iterable(filename):
                # Save the format and pass this to ChainReader
                reader_format = 'CHAIN'
                kwargs['format'] = reader_format
            try:
                reader = get_reader_for(filename,
                                        format=reader_format)
            except TypeError as err:
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
            self.transfer_to_memory(kwargs.get("in_memory_frame_interval", 1))

        return filename, self.trajectory.format

    def transfer_to_memory(self, frame_interval=1):
        """Transfer the trajectory to in memory representation.

        Replaces the current trajectory reader object with one of type
        :class:`MDAnalysis.coordinates.memory.MemoryReader` to support in-place
        editing of coordinates.

        Parameters
        ----------
        frame_interval : int, optional
            Read in every nth frame. [1]

        .. versionadded:: 0.16.0
        """
        from ..coordinates.memory import MemoryReader

        if not isinstance(self.trajectory, MemoryReader):
            # Try to extract coordinates using Timeseries object
            # This is significantly faster, but only implemented for certain
            # trajectory file formats
            try:
                coordinates = self.trajectory.timeseries(
                    self.atoms, format='afc', step=frame_interval)
            # if the Timeseries extraction fails,
            # fall back to a slower approach
            except AttributeError:
                coordinates = np.array([np.copy(ts.positions[:]) for ts in
                                        self.trajectory[::frame_interval]])
                coordinates = coordinates.swapaxes(0, 1)

          # Overwrite trajectory in universe with an MemoryReader
            # object, to provide fast access and allow coordinates
            # to be manipulated
            self.trajectory = MemoryReader(
                coordinates,
                dimensions=self.trajectory.ts.dimensions,
                dt=self.trajectory.ts.dt)

    def select_atoms(self, sel, *othersel, **selgroups):
        return self.atoms.select_atoms(sel, *othersel, **selgroups)

    @property
    def bonds(self):
        return self.atoms.bonds

    @property
    def angles(self):
        return self.atoms.angles

    @property
    def dihedrals(self):
        return self.atoms.dihedrals

    @property
    def impropers(self):
        return self.atoms.impropers

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
        if (attr.per_object is not None and
            len(attr) != n_dict[attr.per_object]):
            raise ValueError('Length of {attr} does not'
                             ' match number of {obj}s.\n'
                             'Expect: {n:d} Have: {m:d}'.format(
                                 attr=attr.attrname,
                                 obj=attr.per_object,
                                 n=n_dict[attr.per_object],
                                 m=len(attr)))

        self._class_bases[GroupBase]._add_prop(attr)

        for cls in attr.target_classes:
            try:
                self._class_bases[cls]._add_prop(attr)
            except (KeyError, AttributeError):
                pass

        for cls in (Atom, Residue, Segment, GroupBase,
                     AtomGroup, ResidueGroup, SegmentGroup):
            try:
                for funcname, meth in attr.transplants[cls]:
                    setattr(self._class_bases[cls], funcname, meth)
            except AttributeError:
                # not every Attribute will have a transplant dict
                pass

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
        """
        bonds = self.atoms.bonds

        class _fragset(object):
            __slots__ = ['ats']
            """Normal sets aren't hashable, this is"""

            def __init__(self, ats):
                self.ats = set(ats)

            def __iter__(self):
                return iter(self.ats)

            def add(self, other):
                self.ats.add(other)

            def update(self, other):
                self.ats.update(other.ats)

        # each atom starts with its own list
        f = dict.fromkeys(self.atoms, None)

        for a1, a2 in bonds:
            if not (f[a1] or f[a2]):
                # New set made here
                new = _fragset([a1, a2])
                f[a1] = f[a2] = new
            elif f[a1] and not f[a2]:
                # If a2 isn't in a fragment, add it to a1's
                f[a1].add(a2)
                f[a2] = f[a1]
            elif not f[a1] and f[a2]:
                # If a1 isn't in a fragment, add it to a2's
                f[a2].add(a1)
                f[a1] = f[a2]
            elif f[a1] is f[a2]:
                # If they're in the same fragment, do nothing
                continue
            else:
                # If they are both in different fragments, combine fragments
                f[a1].update(f[a2])
                f.update(dict((a, f[a1]) for a in f[a2]))

        # Lone atoms get their own fragment
        f.update(dict((a, _fragset((a,)))
                      for a, val in f.items() if not val))

        # All the unique values in f are the fragments
        frags = tuple(
            [AtomGroup(np.sort(np.array([at.index for at in ag])), self)
             for ag in set(f.values())],
        )

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

    :Returns: an instance of :class:`~MDAnalaysis.AtomGroup.Universe`
    """
    if len(args) == 0:
        raise TypeError("as_Universe() takes at least one argument (%d given)" % len(args))
    elif len(args) == 1 and isinstance(args[0], Universe):
        return args[0]
    return Universe(*args, **kwargs)

asUniverse = deprecate(as_Universe, old_name='asUniverse', new_name='as_Universe')


def Merge(*args):
    """Return a :class:`Universe` from two or more :class:`AtomGroup` instances.

    The resulting universe will only inherit the common topology attributes that
    all merged universes share.

    :class:`AtomGroup` instances can come from different Universes, or come
    directly from a :meth:`~Universe.select_atoms` call.

    It can also be used with a single :class:`AtomGroup` if the user wants to,
    for example, re-order the atoms in the Universe.

    If multiple :class:`AtomGroup` instances from the same Universe are given,
    the merge will first simply "add" together the :class:`AtomGroup` instances.

    Parameters
    ----------
    args : One or more :class:`AtomGroup` instances.

    Returns
    -------
    universe : An instance of :class:`~MDAnalaysis.AtomGroup.Universe`

    :Raises: :exc:`ValueError` for too few arguments or if an AtomGroup is
             empty and :exc:`TypeError` if arguments are not
             :class:`AtomGroup` instances.

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

    Note
    ----
        Merging does not create a full trajectory but only a single
        structure even if the input consists of one or more trajectories.

    .. versionchanged 0.9.0::
       Raises exceptions instead of assertion errors.

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
            mapping = {a.index:i for i, a in enumerate(ag, start=offset)}
            offset += len(ag)

            tg = getattr(ag, t)
            bonds_class = type(getattr(ag.universe._topology, t))
            # Create a topology group of only bonds that are within this ag
            # ie we don't want bonds that extend out of the atomgroup
            tg = tg.atomgroup_intersection(ag, strict=True)

            # Map them so they refer to our new indices
            new_idx = [tuple(map(lambda x:mapping[x], entry))
                       for entry in tg.indices]
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
        residx.extend(map(lambda x:res_mapping[x], ag.resindices))
        segidx.extend(map(lambda x:seg_mapping[x], ag.segindices))

    residx = np.array(residx, dtype=np.int32)
    segidx = np.array(segidx, dtype=np.int32)

    _, _, [segidx] = squash_by(residx, segidx)

    n_residues = len(set(residx))
    n_segments = len(set(segidx))

    top = Topology(n_atoms, n_residues, n_segments,
                   attrs=attrs,
                   atom_resindex=residx,
                   residue_segindex=segidx)

    # Create blank Universe only from topology
    u = Universe(top)

    # Take one frame of coordinates from combined atomgroups
    coords = np.vstack([a.positions for a in args])
    trajectory = MDAnalysis.coordinates.base.Reader(None)
    ts = MDAnalysis.coordinates.base.Timestep.from_coordinates(coords)
    setattr(trajectory, "ts", ts)
    trajectory.n_frames = 1
    u.trajectory = trajectory

    return u
