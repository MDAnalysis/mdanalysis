import numpy as np
from numpy.lib.utils import deprecate
import logging
import itertools

import MDAnalysis
from ..lib import util
from ..lib.util import cached
from . import groups
from .topology import Topology

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
    permissive
        Currently only relevant for PDB files: Set to ``True`` in order to
        ignore most errors and read typical MD simulation PDB files; set to
        ``False`` to read with the Bio.PDB reader, which can be useful for real
        Protein Databank PDB files. ``None``  selects the MDAnalysis default
        (which is set in :class:`MDAnalysis.core.flags`) [``None``]
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
        self._cache = {}

        if len(args) == 0:
            # create an empty universe
            self._topology = None
            self.atoms = None
            return

        coordinatefile = args[1:]

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

        # generate Universe version of each class
        # AG, RG, SG, A, R, S
        self._classes = groups.make_classes()

        # Put Group level stuff from topology into class
        for attr in self._topology.attrs:
            self._process_attr(attr)

        # Generate atoms, residues and segments
        self.atoms = self._classes['atomgroup'](
                np.arange(self._topology.n_atoms), self)

        self.residues = self._classes['residuegroup'](
                np.arange( self._topology.n_residues), self)

        self.segments = self._classes['segmentgroup'](np.arange(
            self._topology.n_segments), self)

        # Update Universe namespace with segids
        for seg in self.segments:
            if hasattr(seg, 'segid'):
                if seg.segid[0].isdigit():
                    name = 's' + seg.segid
                else:
                    name = seg.segid
                self.__dict__[name] = seg

        # Load coordinates
        self.load_new(coordinatefile, **kwargs)

    @property
    def universe(self):
        # for Writer.write(universe), see Issue 49
        # Encapsulation in an accessor prevents the Universe from
        # having to keep a reference to itself,
        #  which might be undesirable if it has a __del__ method.
        # It is also cleaner than a weakref.
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
        self._process_attr(topologyattr)

    def _process_attr(self, attr):
        """Squeeze a topologyattr for its information

        Grabs:
         - Group properties (attribute access)
         - Component properties
         - Transplant methods
        """
        self._classes['group']._add_prop(attr)

        for level in attr.target_levels:
            try:
                self._classes[level]._add_prop(attr)
            except (KeyError, AttributeError):
                pass

        for dest in ['atom', 'residue', 'segment', 'group',
                     'atomgroup', 'residuegroup', 'segmentgroup']:
            try:
                for funcname, meth in attr.transplants[dest]:
                    setattr(self._classes[dest], funcname, meth)
            except AttributeError:
                # not every Attribute will have a transplant dict
                pass

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
        AG = self._classes['atomgroup']
        frags = tuple([AG(np.array([at.index for at in ag]), self)
                       for ag in set(f.values())])

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

    if len(args) == 0:
        raise ValueError("Need at least one AtomGroup for merging")

    for a in args:
        if not isinstance(a, groups.AtomGroup):
            raise TypeError(repr(a) + " is not an AtomGroup")
    for a in args:
        if len(a) == 0:
            raise ValueError("cannot merge empty AtomGroup")

    # If any atom groups come from the same Universe, just add them
    # together first
    ag_dict = {}
    for a in args:
        try:
            ag_dict[a.universe].append(a)
        except KeyError:
            ag_dict[a.universe] = [a]

    disjoint_atom_groups = []
    for ag in ag_dict.values():
        disjoint_atom_groups.append(sum(ag))

    # Create a new topology using the intersection of topology attributes
    blank_topology_attrs = set(dir(Topology(attrs=[])))
    common_attrs = set.intersection(*[set(dir(ag.universe._topology)) 
                                      for ag in disjoint_atom_groups])
    topology_groups = set(['bonds', 'angles', 'dihedrals', 'impropers'])

    # Create set of attributes which are array-valued and can be simply
    # concatenated together
    keep_attrs = common_attrs - blank_topology_attrs - topology_groups

    attrs = []
    dtypes = {}
    for attrname in keep_attrs:
        for ag in disjoint_atom_groups:
            attr = getattr(ag, attrname)
            attr_class = type(getattr(ag.universe._topology, attrname))
            if type(attr) != np.ndarray:
                raise TypeError('Encountered unexpected topology'+
                                'attribute of type {}'.format(
                                type(attr)))
            try:
                attr_array.extend(attr)
            except NameError:
                attr_array = list(attr)
        attrs.append(attr_class(np.array(attr_array,
                                        dtype=attr.dtype)))
        del attr_array

    # Build up topology groups
    for tg in (topology_groups & common_attrs):
        bondidx = []
        types = []
        offset = 0
        for ag in disjoint_atom_groups:
            bonds = getattr(ag, tg)
            bond_class = type(getattr(ag.universe._topology, tg))
            bondidx.extend(bonds.indices + offset)
            if hasattr(bonds, '_bondtypes'):
                types.extend(bonds.types())
            else:
                types.extend([None]*len(bonds))
            offset += len(ag)
        bondidx = np.array(bondidx, dtype=np.int32)
        if any(t is None for t in types):
            attrs.append(bond_class(values))
        else:
            types = np.array(types, dtype='|S8')
            attrs.append(bond_class(bondidx, types))

    # Renumber residue and segment indices
    n_atoms = sum([len(ag) for ag in disjoint_atom_groups])
    residx = []
    segidx = []
    for ag in disjoint_atom_groups:
        res_offset = len(set(residx))
        resdict = {n: i+res_offset for i, n in enumerate(set(ag.resindices))}
        seg_offset = len(set(segidx))
        segdict = {n: i+len(set(segidx)) for i, n in enumerate(set(
                                                     ag.segindices))}
        residx.extend([resdict[n] for n in ag.resindices])
        segidx.extend([segdict[n] for n in ag.segindices])

    residx = np.array(residx, dtype=np.int32)
    segidx = np.array(segidx, dtype=np.int32)

    n_residues = len(set(residx))
    n_segments = len(set(segidx))

    top = Topology(n_atoms, n_residues, n_segments,
                   attrs=attrs,
                   atom_resindex=residx,
                   residue_segindex=segidx)

    # Create blank Universe and put topology in it
    u = Universe()
    u._topology = top

    # generate Universe version of each class
    # AG, RG, SG, A, R, S
    u._classes = groups.make_classes()

    # Put Group level stuff from topology into class
    for attr in u._topology.attrs:
        u._process_attr(attr)

    # Generate atoms, residues and segments
    u.atoms = u._classes['atomgroup'](
            np.arange(u._topology.n_atoms), u)

    u.residues = u._classes['residuegroup'](
            np.arange( u._topology.n_residues), u)

    u.segments = u._classes['segmentgroup'](np.arange(
        u._topology.n_segments), u)

    # Update Universe namespace with segids
    for seg in u.segments:
        if hasattr(seg, 'segid'):
            if seg.segid[0].isdigit():
                name = 's' + seg.segid
            else:
                name = seg.segid
            u.__dict__[name] = seg

    coords = np.vstack([a.positions for a in disjoint_atom_groups])
    trajectory = MDAnalysis.coordinates.base.Reader(None)
    ts = MDAnalysis.coordinates.base.Timestep.from_coordinates(coords)
    setattr(trajectory, "ts", ts)
    trajectory.n_frames = 1
    u.trajectory = trajectory

    return u
