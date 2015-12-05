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

    *Attributes:*

    - :attr:`Universe.trajectory`: currently loaded trajectory reader;
      :attr:`Universe.trajectory.ts` is the current time step
    - :attr:`Universe.dimensions`: current system dimensions (simulation unit cell, if
      set in the trajectory)
    - :attr:`Universe.bonds`: TopologyGroup of bonds in Universe, also
      :attr:`Universe.angles`, :attr:`Universe.dihedrals`, and :attr:`Universe.impropers`
      (low level access through :attr:`Universe._topology`)

    .. Note::

       If atom attributes such as element, mass, or charge are not explicitly
       provided in the topology file then MDAnalysis tries to guess them (see
       :mod:`MDAnalysis.topology.tables`). This does not always work and if you
       require correct values (e.g. because you want to calculate the center of
       mass) then you need to make sure that MDAnalysis gets all the
       information needed.

    .. versionchanged:: 0.7.5
       Can also read multi-frame PDB files with the
       :class:`~MDAnalysis.coordinates.PDB.PrimitivePDBReader`.

    .. versionchanged:: 0.8
       Parse arbitrary number of arguments as a single topology file and a
       sequence of trajectories.

    .. versionchanged:: 0.9.0
       Topology information now loaded lazily, but can be forced with
       :meth:`build_topology`. Changed :attr:`bonds` attribute to be a
       :class:`~MDAnalysis.core.topologyobjects.TopologyGroup`. Added :attr:`angles`
       and :attr:`torsions` attribute as
       :class:`~MDAnalysis.core.topologyobjects.TopologyGroup`. Added fragments to
       Universe cache

    .. versionchanged:: 0.11.0
       :meth:`make_anchor`, :meth:`remove_anchor`, :attr:`is_anchor`, and
       :attr:`anchor_name` were added to support the pickling/unpickling of
       :class:`AtomGroup`.
       Deprecated :meth:`selectAtoms` in favour of :meth:`select_atoms`.
    """

    def __init__(self, *args, **kwargs):
        """Initialize the central MDAnalysis Universe object.

        :Arguments:
          *topologyfile*
             A CHARMM/XPLOR PSF topology file, PDB file or Gromacs GRO file; used to define the
             list of atoms. If the file includes bond information, partial
             charges, atom masses, ... then these data will be available to
             MDAnalysis. A "structure" file (PSF, PDB or GRO, in the sense of a
             topology) is always required.
          *coordinatefile*
             A trajectory (such as CHARMM DCD, Gromacs XTC/TRR/GRO, XYZ, XYZ.bz2) or a PDB that
             will provide coordinates, possibly multiple frames.
             If a **list of filenames** is provided then they are sequentially read and appear
             as one single trajectory to the Universe. The list can contain different file
             formats.

             .. deprecated:: 0.8
                Do not use the *coordinatefile* keyword argument, just provide trajectories as
                positional arguments.

          *permissive*
             currently only relevant for PDB files: Set to ``True`` in order to ignore most errors
             and read typical MD simulation PDB files; set to ``False`` to read with the Bio.PDB reader,
             which can be useful for real Protein Databank PDB files. ``None``  selects the
             MDAnalysis default (which is set in :class:`MDAnalysis.core.flags`) [``None``]
          *topology_format*
             provide the file format of the topology file; ``None`` guesses it from the file
             extension [``None``]
             Can also pass a subclass of :class:`MDAnalysis.topology.base.TopologyReader`
             to define a custom reader to be used on the topology file.
          *format*
             provide the file format of the coordinate or trajectory file;
             ``None`` guesses it from the file extension. Note that this
             keyword has no effect if a list of file names is supplied because
             the "chained" reader has to guess the file format for each
             individual list member. [``None``]
             Can also pass a subclass of :class:`MDAnalysis.coordinates.base.Reader`
             to define a custom reader to be used on the trajectory file.
          *guess_bonds*
              Once Universe has been loaded, attempt to guess the connectivity
              between atoms.  This will populate the .bonds .angles and
              .dihedrals attributes of the Universe.
          *vdwradii*
              For use with *guess_bonds*. Supply a dict giving a vdwradii for each atom type
              which are used in guessing bonds.
          *is_anchor*
              When unpickling instances of :class:`MDAnalysis.core.AtomGroup.AtomGroup`
              existing Universes are searched for one where to anchor those atoms. Set
              to ``False`` to prevent this Universe from being considered. [``True``]
          *anchor_name*
              Setting to other than ``None`` will cause :class:`MDAnalysis.core.AtomGroup.AtomGroup`
              instances pickled from the Universe to only unpickle if a compatible
              Universe with matching *anchor_name* is found. *is_anchor* will be ignored in
              this case but will still be honored when unpickling :class:`MDAnalysis.core.AtomGroup.AtomGroup`
              instances pickled with *anchor_name*==``None``. [``None``]


        This class tries to do the right thing:

        1. If file with topology and coordinate information (such as PDB, GRO,
           CRD, ...) is provided instead of a topology file and no
           *coordinatefile* then the coordinates are taken from the first
           file. Thus you can load a functional universe with ::

              u = Universe('1ake.pdb')

           If you want to specify the coordinate file format yourself you can
           do so using the *format* keyword::

              u = Universe('1ake.ent1', format='pdb')

        2. If only a topology file without coordinate information is provided
           one will have to load coordinates manually using
           :meth:`Universe.load_new`. The file format of the topology file
           can be explicitly set with the *topology_format* keyword.

        .. versionchanged:: 0.7.4
           New *topology_format* and *format* parameters to override the file
           format detection.
        .. versionchanged:: 0.10.0
           Added ``'guess_bonds'`` keyword to cause topology to be guessed on
           Universe creation.
           Deprecated ``'bonds'`` keyword, use ``'guess_bonds'`` instead.
        .. versionchanged:: 0.11.0
           Added the *is_anchor* and *anchor_name* keywords for finer behavior
           control when unpickling instances of :class:`MDAnalysis.core.AtomGroup.AtomGroup`.
        """

        from ..topology.core import get_parser_for
        from ..topology.base import TopologyReader
        from ..coordinates.base import ProtoReader

        # managed attribute holding Reader
        self._trajectory = None

        # Cache is used to store objects which are built lazily into Universe
        # Currently cached objects (managed property name and cache key):
        # - bonds
        # - angles
        # - dihedrals
        # - improper dihedrals
        # - fragments
        # Cached stuff is handled using util.cached decorator
        self._cache = dict()

        if len(args) == 0:
            # create an empty universe
            self._topology = None
            self.atoms = None
            return

        self.filename = args[0]

        # old behaviour (explicit coordfile) overrides new behaviour
        coordinatefile = kwargs.pop('coordinatefile', args[1:])
        topology_format = kwargs.pop('topology_format', None)

        if len(args) == 1 and not coordinatefile:
            # special hacks to treat a coordinate file as a coordinate AND topology file
            # coordinatefile can be None or () (from an empty slice args[1:])
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

        # Generate atoms, residues and segments
        self.atoms = groups.AtomGroup(
            np.arange(self._topology.n_atoms), self)
        self.residues = groups.ResidueGroup(np.arange(
            self._topology.n_residues), self)
        self.segments = groups.SegmentGroup(np.arange(
            self._topology.n_segments), self)
        # Add segids to Universe attribute namespace
        for seg in self.segments:
            if seg.id[0].isdigit():
                name = 's' + seg.id
            else:
                name = seg.id
            self.__dict__[name] = seg

        # Load coordinates
        self.load_new(coordinatefile, **kwargs)

        # Deprecated bonds mode handling here, remove eventually.
        if 'bonds' in kwargs:
            warnings.warn("The 'bonds' keyword has been deprecated"
                          " and will be removed in 0.11.0."
                          " Please use 'guess_bonds' instead.")
            if kwargs.get('bonds') in ['all', True]:
                kwargs['guess_bonds'] = True

        if kwargs.get('guess_bonds', False):
            self.atoms.guess_bonds(vdwradii=kwargs.get('vdwradii',None))

        # For control of AtomGroup unpickling
        if kwargs.get('is_anchor', True):
            self.make_anchor()
        self.anchor_name = kwargs.get('anchor_name')

    def _clear_caches(self, *args):
        """Clear cache for all *args*.

        If not args are provided, all caches are cleared.

        .. versionadded 0.9.0
        """
        if len(args) == 0:
            self._cache = dict()
        else:
            for name in args:
                try:
                    del self._cache[name]
                except KeyError:
                    pass

    def _fill_cache(self, name, value):
        """Populate _cache[name] with value.

        .. versionadded:: 0.9.0
        """
        self._cache[name] = value

    def _init_top(self, cat, Top):
        """Initiate a generic form of topology.

        Arguments:
          *cat*
            The key which will be searched in the _topology dict.
            The key "guessed_" + cat will also be searched.
          *Top*
            Class of the topology object to be created.

        .. versionadded:: 0.10.0
        """
        defined = self._topology.get(cat, set())
        guessed = self._topology.get('guessed_' + cat, set())

        TopSet = top.TopologyGroup.from_indices(defined, self.atoms,
                                                            bondclass=Top, guessed=False,
                                                            remove_duplicates=True)
        TopSet += top.TopologyGroup.from_indices(guessed, self.atoms,
                                                             bondclass=Top, guessed=True,
                                                             remove_duplicates=True)

        return TopSet

    def _init_bonds(self):
        """Set bond information from u._topology['bonds']

        .. versionchanged:: 0.9.0
           Now returns a :class:`~MDAnalysis.core.topologyobjects.TopologyGroup`
        """
        bonds = self._init_top('bonds', top.Bond)

        bondorder = self._topology.get('bondorder', None)
        if bondorder:
            for b in bonds:
                try:
                    b.order = bondorder[b.indices]
                except KeyError:
                    pass

        return bonds

    def _init_angles(self):
        """Builds angle information from u._topology['angles']

        Returns a :class:`~MDAnalysis.core.topologyobjects.TopologyGroup`

        .. versionadded:: 0.9.0
        .. versionchanged:: 0.10.0
           Now reads guessed angles and tags them appropriately.
        """
        return self._init_top('angles', top.Angle)

    def _init_dihedrals(self):
        """Builds dihedral information from u._topology['dihedrals']

        Returns a :class:`~MDAnalysis.core.topologyobjects.TopologyGroup`

        .. versionadded:: 0.9.0
        .. versionchanged:: 0.10.0
           Now reads guessed torsions and tags them appropriately.
        .. versionchanged:: 0.11.0
           Renamed to _init_dihedrals (was _init_torsions)
        """
        return self._init_top('dihedrals', top.Dihedral)

    def _init_impropers(self):
        """Build improper dihedral information from u._topology['impropers']

        Returns a :class:`~MDAnalysis.core.topologyobjects.TopologyGroup`

        .. versionadded:: 0.9.0
        .. versionchanged:: 0.10.0
        """
        return self._init_top('impropers', top.ImproperDihedral)

    def _init_fragments(self):
        """Build all fragments in the Universe

        Generally built on demand by an Atom querying its fragment property.

        .. versionadded:: 0.9.0
        """
        # Check that bond information is present, else inform
        bonds = self.bonds
        if not bonds:
            raise NoDataError("Fragments require that the Universe has Bond information")

        # This current finds all fragments from all Atoms
        # Could redo this to only find fragments for a queried atom (ie. only fill out
        # a single fragment).  This would then make it scale better for large systems.
        # eg:
        # try:
        #    return self._fragDict[a]
        # except KeyError:
        #    self._init_fragments(a)  # builds the fragment a belongs to

        class _fragset(object):
            """Normal sets aren't hashable, this is"""

            def __init__(self, ats):
                self.ats = set(ats)

            def __iter__(self):
                return iter(self.ats)

            def add(self, other):
                self.ats.add(other)

            def update(self, other):
                self.ats.update(other.ats)

        f = dict.fromkeys(self.atoms, None)  # each atom starts with its own list

        for a1, a2 in bonds:  # Iterate through all bonds
            if not (f[a1] or f[a2]):  # New set made here
                new = _fragset([a1, a2])
                f[a1] = f[a2] = new
            elif f[a1] and not f[a2]:  # If a2 isn't in a fragment, add it to a1's
                f[a1].add(a2)
                f[a2] = f[a1]
            elif not f[a1] and f[a2]:  # If a1 isn't in a fragment, add it to a2's
                f[a2].add(a1)
                f[a1] = f[a2]
            elif f[a1] is f[a2]:  # If they're in the same fragment, do nothing
                continue
            else:  # If they are both in different fragments, combine fragments
                f[a1].update(f[a2])
                f.update(dict((a, f[a1]) for a in f[a2]))

                # Lone atoms get their own fragment
        f.update(dict((a, _fragset((a,))) for a, val in f.items() if not val))

        # All the unique values in f are the fragments
        frags = tuple([AtomGroup(list(a.ats)) for a in set(f.values())])

        return frags

    @property
    def universe(self):
        # for Writer.write(universe), see Issue 49
        # Encapsulation in an accessor prevents the Universe from having to keep a reference to itself,
        #  which might be undesirable if it has a __del__ method. It is also cleaner than a weakref.
        return self

    @property
    @cached('fragments')
    def fragments(self):
        """Read only tuple of fragments in the Universe

        .. versionadded 0.9.0
        """
        return self._init_fragments()

    @property
    @cached('bondDict')
    def _bondDict(self):
        """Lazily built dictionary of bonds

        Translates Atom to list of bonds

        .. versionadded:: 0.9.0
        """
        bonds = self.bonds
        bd = defaultdict(list)

        if not bonds:
            pass
        else:
            for b in bonds:
                for a in b:
                    bd[a].append(b)
        return bd

    @property
    @cached('angleDict')
    def _angleDict(self):
        """Lazily built dictionary of angles

        Translates Atom to list of angles

        .. versionadded:: 0.9.0
        """
        bonds = self.angles
        bd = defaultdict(list)

        if not bonds:
            pass
        else:
            for b in bonds:
                for a in b:
                    bd[a].append(b)
        return bd

    @property
    @cached('dihedralDict')
    def _dihedralDict(self):
        """Lazily built dictionary of dihedrals

        Translates Atom to list of dihedrals

        .. versionadded:: 0.9.0
        .. versionchanged:: 0.11.0
           Renamed to _dihedralDict (was _torsionDict)
        """
        bonds = self.dihedrals
        bd = defaultdict(list)

        if bonds is None:
            pass
        else:
            for b in bonds:
                for a in b:
                    bd[a].append(b)
        return bd

    @property
    @cached('improperDict')
    def _improperDict(self):
        """Lazily built dictionary of improper dihedrals

        Translates Atom to list of improper dihedrals

        .. versionadded:: 0.9.0
        """
        bonds = self.impropers
        bd = defaultdict(list)

        if bonds is None:
            pass
        else:
            for b in bonds:
                for a in b:
                    bd[a].append(b)
        return bd

    @property
    @cached('fragDict')
    def _fragmentDict(self):
        """Lazily built dictionary of fragments.

        Translates :class:`Atom` objects into the fragment they belong to.

        The Atom.fragment managed property queries this dictionary.

        .. versionadded 0.9.0
        """
        frags = self.fragments  # will build if not built
        fd = dict()
        for f in frags:
            for a in f:
                fd[a] = f
        return fd

    def build_topology(self):
        """
        Bond angle and dihedral information is lazily constructed into the
        Universe.

        This method forces all this information to be loaded.

        .. versionadded 0.9.0
        """
        if 'bonds' not in self._cache:
            self._cache['bonds'] = self._init_bonds()
        if 'angles' not in self._cache:
            self._cache['angles'] = self._init_angles()
        if 'dihedrals' not in self._cache:
            self._cache['dihedrals'] = self._init_dihedrals()
        if 'impropers' not in self._cache:
            self._cache['impropers'] = self._init_impropers()

    @property
    @cached('bonds')
    def bonds(self):
        """
        Returns a :class:`~MDAnalysis.core.topologyobjects.TopologyGroup` of all
        bonds in the Universe.

        .. versionchanged:: 0.9.0
           Now a lazily built :class:`~MDAnalysis.core.topologyobjects.TopologyGroup`
        .. versionchanged:: 0.9.2
           Now can return empty TopologyGroup
        """
        return self._init_bonds()

    @bonds.setter
    def bonds(self, bondlist):
        """Can set bonds by supplying an iterable of bond tuples.

        Each bond tuple must contain the zero based indices of the two Atoms in
        the bond

        .. versionadded:: 0.9.0
        """
        self._fill_cache('bonds', bondlist)
        self._clear_caches('bondDict')

    @bonds.deleter
    def bonds(self):
        """Delete the bonds from Universe

        This must also remove the per atom record of bonds (bondDict)

        .. versionadded:: 0.9.0
        """
        self._clear_caches('bonds', 'bondDict')

    @property
    @cached('angles')
    def angles(self):
        """
        Returns a :class:`~MDAnalysis.core.topologyobjects.TopologyGroup`
        of all angles in the Universe

        .. versionadded:: 0.9.0
        .. versionchanged:: 0.9.2
           Now can return empty TopologyGroup
        """
        return self._init_angles()

    @angles.setter
    def angles(self, bondlist):
        self._fill_cache('angles', bondlist)
        self._clear_caches('angleDict')

    @angles.deleter
    def angles(self):
        self._clear_caches('angles', 'angleDict')

    @property
    @cached('dihedrals')
    def dihedrals(self):
        """
        Returns a :class:`~MDAnalysis.core.topologyobjects.TopologyGroup`
        of all dihedrals in the Universe

        .. versionadded:: 0.9.0
        .. versionchanged:: 0.9.2
           Now can return empty TopologyGroup
        """
        return self._init_dihedrals()

    @dihedrals.setter
    def dihedrals(self, bondlist):
        self._fill_cache('dihedrals', bondlist)
        self._clear_caches('dihedralDict')

    @dihedrals.deleter
    def dihedrals(self):
        self._clear_caches('dihedrals', 'dihedralDict')

    @property
    @cached('impropers')
    def impropers(self):
        """
        Returns a :class:`~MDAnalysis.core.topologyobjects.TopologyGroup`
        of all improper dihedrals in the Universe

        .. versionadded:: 0.9.0
        .. versionchanged:: 0.9.2
           Now can return empty TopologyGroup
        """
        return self._init_impropers()

    @impropers.setter
    def impropers(self, bondlist):
        self._fill_cache('impropers', bondlist)
        self._clear_caches('improperDict')

    @impropers.deleter
    def impropers(self):
        self._clear_caches('impropers', 'improperDict')

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
        """Select atoms using a CHARMM selection string.

        Returns an :class:`AtomGroup` with atoms sorted according to their
        index in the psf (this is to ensure that there aren't any duplicates,
        which can happen with complicated selections).

        Existing :class:`AtomGroup` objects can be passed as named arguments,
        which will then be available to the selection parser.

        Subselections can be grouped with parentheses.

        Example::
           >>> sel = universe.select_atoms("segid DMPC and not ( name H* or name O* )")
           >>> sel
           <AtomGroup with 3420 atoms>

           >>> universe.select_atoms("around 10 group notHO", notHO=sel)
           <AtomGroup with 1250 atoms>

        .. Note::

           If exact ordering of atoms is required (for instance, for
           :meth:`~AtomGroup.angle` or :meth:`~AtomGroup.dihedral`
           calculations) then one supplies selections *separately* in the
           required order. Also, when multiple :class:`AtomGroup` instances are
           concatenated with the ``+`` operator then the order of :class:`Atom`
           instances is preserved and duplicates are not removed.

        .. SeeAlso:: :ref:`selection-commands-label` for further details and examples.

        The selection parser understands the following CASE SENSITIVE *keywords*:

        **Simple selections**

            protein, backbone, nucleic, nucleicbackbone
                selects all atoms that belong to a standard set of residues;
                a protein is identfied by a hard-coded set of residue names so
                it  may not work for esoteric residues.
            segid *seg-name*
                select by segid (as given in the topology), e.g. ``segid 4AKE``
                or ``segid DMPC``
            resid *residue-number-range*
                resid can take a single residue number or a range of numbers. A
                range consists of two numbers separated by a colon (inclusive)
                such as ``resid 1:5``. A residue number ("resid") is taken
                directly from the topology.
            resnum *resnum-number-range*
                resnum is the canonical residue number; typically it is set to
                the residue id in the original PDB structure.
            resname *residue-name*
                select by residue name, e.g. ``resname LYS``
            name *atom-name*
                select by atom name (as given in the topology). Often, this is
                force field dependent. Example: ``name CA`` (for C&alpha; atoms)
                or ``name OW`` (for SPC water oxygen)
            type *atom-type*
                select by atom type; this is either a string or a number and
                depends on the force field; it is read from the topology file
                (e.g. the CHARMM PSF file contains numeric atom types). It has
                non-sensical values when a PDB or GRO file is used as a topology
            atom *seg-name*  *residue-number*  *atom-name*
                a selector for a single atom consisting of segid resid atomname,
                e.g. ``DMPC 1 C2`` selects the C2 carbon of the first residue of
                the DMPC segment
            altloc *alternative-location*
                a selection for atoms where alternative locations are available,
                which is often the case with high-resolution crystal structures
                e.g. `resid 4 and resname ALA and altloc B` selects only the
                atoms of ALA-4 that have an altloc B record.

        **Boolean**

            not
                all atoms not in the selection, e.g. ``not protein`` selects
                all atoms that aren't part of a protein

            and, or
                combine two selections according to the rules of boolean
                algebra, e.g. ``protein and not (resname ALA or resname LYS)``
                selects all atoms that belong to a protein, but are not in a
                lysine or alanine residue

        **Geometric**

            around *distance*  *selection*
                selects all atoms a certain cutoff away from another selection,
                e.g. ``around 3.5 protein`` selects all atoms not belonging to
                protein that are within 3.5 Angstroms from the protein
            point *x* *y* *z*  *distance*
                selects all atoms within a cutoff of a point in space, make sure
                coordinate is separated by spaces,
                e.g. ``point 5.0 5.0 5.0  3.5`` selects all atoms within 3.5
                Angstroms of the coordinate (5.0, 5.0, 5.0)
            prop [abs] *property*  *operator*  *value*
                selects atoms based on position, using *property*  **x**, **y**,
                or **z** coordinate. Supports the **abs** keyword (for absolute
                value) and the following *operators*: **<, >, <=, >=, ==, !=**.
                For example, ``prop z >= 5.0`` selects all atoms with z
                coordinate greater than 5.0; ``prop abs z <= 5.0`` selects all
                atoms within -5.0 <= z <= 5.0.
            sphzone *radius* *selection*
                Selects all atoms that are within *radius* of the center of
                geometry of *selection*
            sphlayer *inner radius* *outer radius* *selection*
                Similar to sphzone, but also excludes atoms that are within
                *inner radius* of the selection COG

        **Connectivity**

            byres *selection*
                selects all atoms that are in the same segment and residue as
                selection, e.g. specify the subselection after the byres keyword
            bonded *selection*
                selects all atoms that are bonded to selection
                eg: ``select name H bonded name O`` selects only hydrogens
                bonded to oxygens

        **Index**

            bynum *index-range*
                selects all atoms within a range of (1-based) inclusive indices,
                e.g. ``bynum 1`` selects the first atom in the universe;
                ``bynum 5:10`` selects atoms 5 through 10 inclusive. All atoms
                in the :class:`MDAnalysis.Universe` are consecutively numbered,
                and the index runs from 1 up to the total number of atoms.

        **Preexisting selections**

            group *group-name*
                selects the atoms in the :class:`AtomGroup` passed to the
                function as an argument named *group-name*. Only the atoms
                common to *group-name* and the instance :meth:`~select_atoms`
                was called from will be considered. *group-name* will be
                 included in the parsing just by comparison of atom indices.
                This means that it is up to the user to make sure they were
                defined in an appropriate :class:`Universe`.

            fullgroup *group-name*
                just like the ``group`` keyword with the difference that all the
                atoms of *group-name* are included. The resulting selection may
                therefore have atoms that were initially absent from the
                instance :meth:`~select_atoms` was called from.

        .. versionchanged:: 0.7.4
           Added *resnum* selection.
        .. versionchanged:: 0.8.1
           Added *group* and *fullgroup* selections.
        .. versionchanged:: 0.13.0
           Added *bonded* selection
        """
        # can ONLY import in method, otherwise cyclical import!
        from . import Selection

        atomgrp = Selection.Parser.parse(sel, selgroups).apply(self)
        if len(othersel) == 0:
            return atomgrp
        else:
            # Generate a selection for each selection string
            #atomselections = [atomgrp]
            for sel in othersel:
                atomgrp = atomgrp + Selection.Parser.parse(sel, selgroups).apply(self)
                #atomselections.append(Selection.Parser.parse(sel).apply(self))
            #return tuple(atomselections)
            return atomgrp

    selectAtoms = deprecate(select_atoms, old_name='selectAtoms',
                            new_name='select_atoms')

    def __repr__(self):
        return "<Universe with {n_atoms} atoms{bonds}>".format(
            n_atoms=len(self.atoms),
            bonds=" and {0} bonds".format(len(self.bonds)) if self.bonds else "")

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

    def make_anchor(self):
        """Add this Universe to the list where anchors are searched for when unpickling
        :class:`MDAnalysis.core.AtomGroup.AtomGroup` instances. Silently proceeds if it
        is already on the list."""
        MDAnalysis._anchor_universes.add(self)

    def remove_anchor(self):
        """Remove this Universe from the list where anchors are searched for when unpickling
        :class:`MDAnalysis.core.AtomGroup.AtomGroup` instances. Silently proceeds if it
        is already not on the list."""
        MDAnalysis._anchor_universes.discard(self)

    @property
    def is_anchor(self):
        """Whether this Universe will be checked for anchoring when unpickling
        :class:`MDAnalysis.core.AtomGroup.AtomGroup` instances"""
        return self in MDAnalysis._anchor_universes

    @property
    def anchor_name(self):
        return self._anchor_name

    @anchor_name.setter
    def anchor_name(self, name):
        """Setting this attribute to anything other than ``None`` causes this Universe to
        be added to the list where named anchors are searched for when unpickling
        :class:`MDAnalysis.core.AtomGroup.AtomGroup` instances (silently proceeding if
        it already is on the list). Setting to ``None`` causes the removal from said list."""
        self._anchor_name = name
        if name is None:
            MDAnalysis._named_anchor_universes.discard(self)
        else:
            MDAnalysis._named_anchor_universes.add(self)

    def _matches_unpickling(self, anchor_name, n_atoms, fname, trajname):
        if anchor_name is None or anchor_name == self.anchor_name:
            try:
                return len(self.atoms)==n_atoms and self.filename==fname and self.trajectory.filenames==trajname
            except AttributeError: # Only ChainReaders have filenames (plural)
                return len(self.atoms)==n_atoms and self.filename==fname and self.trajectory.filename==trajname
        else:
            return False


