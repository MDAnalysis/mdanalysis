# $Id$
"""
:mod:`MDAnalysis.core.AtomGroup` --- AtomGroup Hierarchy
========================================================

The most important data structure in MDAnalysis is the
:class:`AtomGroup`, which contains :class:`Atom` instances.

A :class:`Universe` is the user-visible entry point and collects all
information needed to analyze a structure or a whole trajectory.

**Class Hierarchy**::
    AtomGroup -> ResidueGroup  -> Segment
              -> Residue
    Atom
    Universe

The important classes and functions:

.. autoclass:: Universe
   :members:
.. autoclass:: AtomGroup
   :members:
.. autoclass:: Atom
   :members:
.. autoclass:: Residue
.. autoclass:: ResidueGroup
.. autoclass:: Segment

.. autofunction:: asUniverse
.. autoexception:: SelectionError
"""
import warnings

import numpy
from MDAnalysis import SelectionError

class Atom(object):
    """A single atom definition

    Data: number, segid, resid, resname, name, type, mass, charge

    Methods:
        a = Atom()
        a.pos     - The current position (as a numpy array) of this atom
    """
    __slots__ = ("number", "id", "name", "type", "resname", "resid", "segid", "mass", "charge", "residue", "segment", "bonds", "__universe", "acceptor", "donor")

    def __init__(self, number, name, type, resname, resid, segid, mass, charge):
        self.number = number
        self.name = name
        self.type = type
        self.resname = resname
        self.resid = resid
        self.segid = segid
        self.mass = mass
        self.charge = charge
        self.donor = None     # H-bond properties (filled in later)
        self.acceptor = None
    def __repr__(self):
        return "< Atom " + repr(self.number+1) + ": name " + repr(self.name) +" of type " + \
               repr(self.type) + " of resname " + repr(self.resname) + ", resid " +repr(self.resid) + " and segid " +repr(self.segid)+'>'
    def __cmp__(self, other):
        return cmp(self.number, other.number)
    def __eq__(self, other):
        return self.number == other.number
    def __hash__(self):
        return hash(self.number)
    def __add__(self, other):
        if not (isinstance(other, Atom) or isinstance(other, AtomGroup)):
            raise TypeError('Can only concatenate Atoms (not "'+repr(other.__class__.__name__)+'") to AtomGroup')
        if isinstance(other, Atom): return AtomGroup([self, other])
        else: return AtomGroup([self]+other.atoms)

    def pos():
        doc = "Current cartesian coordinates of the atom."
        def fget(self):
            return self.universe.coord[self.number] # PDB numbering starts at 0
        return locals()
    pos = property(**pos())

    def bfactor():
        doc = "Crystallographic B-factor (if universe was built from a pdb) or None"
        def fget(self):
            try:
                return self.universe.bfactors[self.number] # PDB numbering starts at 0
            except AttributeError:
                return None
        def fset(self,value):
            try:
                self.universe.bfactors[self.number] = value
            except AttributeError:
                self.universe.bfactors = numpy.zeros(self.universe.atoms.numberOfAtoms())
                self.universe.bfactors[self.number] = value
        return locals()
    bfactor = property(**bfactor())

    def universe():
        doc = "a pointer back to the Universe"
        def fget(self):
            if not self.__universe == None: return self.__universe
            else: raise AttributeError("Atom "+repr(self.number)+" is not assigned to a Universe")
        def fset(self, universe):
            self.__universe = universe
        return locals()
    universe = property(**universe())

class AtomGroup(object):
    """A group of atoms

    Currently contains a list of atoms from the main system that correspond to this group.

    Data: atoms - a list of references to the corresponding atoms in Universe.atoms
          AtomGroups are immutable

    Methods:
        ag = universe.selectAtoms("...")
        ag.numberOfAtoms() - return the number of atoms in group
        ag.indices() - return indices into main atom array
        ag.masses() - array of masses
        ag.totalMass() - total mass
        ag.charges() - array of charges
        ag.totalCharge() - total charge
        ag.centerOfGeometry() - center of geometry
        ag.centerOfMass() - center of mass
        ag.radiusOfGyration() - radius of gyration
        ag.principleAxis() - returns the principle axis of rotation
        ag.bfactors() - returns B-factors (if they were loaded from a PDB)
        c = ag.coordinates() - return array of coordinates

        ag.write() - write all atoms in the group to a file
    """
    def atoms():
        doc = "a list of references to atoms in Universe corresponding to a specifies subset"
        def fget(self):
            return self.__atoms
        return locals()
    atoms = property(**atoms())
    
    def _atoms():
        def fget(self):
            warnings.warn("Usage of '_atoms' is deprecated. Use 'atoms' instead.", category=DeprecationWarning, stacklevel=2)
            return self.__atoms
        return locals()
    _atoms = property(**_atoms())

    # Universe pointer is important for Selections to work on groups
    def universe():
        doc = "The universe to which the atoms belong (read-only)."
        def fget(self):
            try:
                return self.__atoms[0].universe
            except AttributeError:
                return None
        return locals()
    universe = property(**universe())

    def __init__(self, atoms):
        if len(atoms) < 1: raise Exception("No atoms defined for AtomGroup")
        # __atoms property is effectively readonly        
        # check that atoms is indexable:
        try:
            atoms[0]
            self.__atoms = atoms      
        except TypeError:
            self.__atoms = list(atoms)
        # If the number of atoms is very large, create a dictionary cache for lookup
        if len(atoms) > 10000:
            self._atom_cache = dict([(x,None) for x in self.__atoms])
    def __len__(self):
        #import warnings
        #warnings.warn("To prevent confusion with AtomGroup subclasses you should use numberOfAtoms() instead", category=Warning, stacklevel=2)
        return self.numberOfAtoms()
    def __getitem__(self, item):
        if (numpy.dtype(type(item)) == numpy.dtype(int)) or (type(item) == slice):
            return self.atoms[item]
        else: return super(AtomGroup, self).__getitem__(item)
    def __getattr__(self, name):
        # There can be more than one atom with the same name
        atomlist = [atom for atom in self.atoms if name == atom.name]
        if len(atomlist) == 0: raise SelectionError("No atoms with name "+name)
        elif len(atomlist) == 1: return atomlist[0]
        else: return AtomGroup(atomlist)
    def __iter__(self):
        return iter(self.atoms)
    def __contains__(self, other):
        if hasattr(self, "_atom_cache"):
            return other in self._atom_cache
        else: return other in self.atoms
    def __add__(self, other):
        if not (isinstance(other, Atom) or isinstance(other, AtomGroup)):
            raise TypeError('Can only concatenate AtomGroup (not "'+repr(other.__class__.__name__)+'") to AtomGroup')
        if isinstance(other, AtomGroup): return AtomGroup(self.atoms + other.atoms)
        else: return AtomGroup(self.atoms+[other])
    def __repr__(self):
        return '<'+self.__class__.__name__+' with '+repr(self.numberOfAtoms())+' atoms>'
    def numberOfAtoms(self):
        return len(self.atoms)
    def indices(self):
        if not hasattr(self,'_cached_indices'):
            self._cached_indices = numpy.array([atom.number for atom in self.atoms])
        return self._cached_indices
    def masses(self):
        if not hasattr(self, "_cached_masses"):
            self._cached_masses = numpy.array([atom.mass for atom in self.atoms])
        return self._cached_masses
    def totalMass(self):
        return numpy.sum(self.masses(), axis=0)
    def charges(self):
        return numpy.array([atom.charge for atom in self.atoms])
    def totalCharge(self):
        return numpy.sum(self.charges(), axis=0)
    def centerOfGeometry(self):
        return numpy.sum(self.coordinates(), axis=0)/self.numberOfAtoms()
    def centerOfMass(self):
        return numpy.sum(self.coordinates()*self.masses()[:,numpy.newaxis],axis=0)/self.totalMass()
    def radiusOfGyration(self):
        masses = self.masses()
        recenteredpos = self.coordinates() - self.centerOfMass()
        rog_sq = numpy.sum(masses*numpy.sum(numpy.power(recenteredpos, 2), axis=1))/self.totalMass()
        return numpy.sqrt(rog_sq)
    def momentOfInertia(self):
        # Convert to local coordinates
        recenteredpos = self.coordinates() - self.centerOfMass()
        masses = self.masses()
        values = zip(masses, recenteredpos)
        # Create the inertia tensor
        # m_i = mass of atom i
        # (x_i, y_i, z_i) = pos of atom i
        # Ixx = sum(m_i*(y_i^2+z_i^2)); Iyy = sum(m_i*(x_i^2+z_i^2)); Izz = sum(m_i*(x_i^2+y_i^2))
        # Ixy = Iyx = -1*sum(m_i*x_i*y_i)
        # Ixz = Izx = -1*sum(m_i*x_i*z_i)
        # Iyz = Izy = -1*sum(m_i*y_i*z_i)
        Ixx = reduce(lambda t,a: t+a[0]*(a[1][1]*a[1][1]+a[1][2]*a[1][2]), values, 0.)
        Iyy = reduce(lambda t,a: t+a[0]*(a[1][0]*a[1][0]+a[1][2]*a[1][2]), values, 0.)
        Izz = reduce(lambda t,a: t+a[0]*(a[1][0]*a[1][0]+a[1][1]*a[1][1]), values, 0.)
        Ixy = Iyx = -1*reduce(lambda t,a: t+a[0]*a[1][0]*a[1][1], values, 0.)
        Ixz = Izx = -1*reduce(lambda t,a: t+a[0]*a[1][0]*a[1][2], values, 0.)
        Iyz = Izy = -1*reduce(lambda t,a: t+a[0]*a[1][1]*a[1][2], values, 0.)
        return numpy.array([[Ixx, Ixy, Ixz],[Iyx, Iyy, Iyz],[Izx, Izy, Izz]])
    def principleAxes(self):
        from numpy.linalg import eig
        eigenval, eigenvec = eig(self.momentOfInertia())
        # Sort
        indices = numpy.argsort(eigenval)
        return eigenvec[:,indices] 
    def coordinates(self, ts=None, copy=False):
        if ts == None:
            ts = self.universe.coord
        return numpy.array(ts[self.indices()], copy=copy)

    def selectAtoms(self, sel, *othersel):
        import Selection
        atomgrp = Selection.Parser.parse(sel).apply(self)
        if len(othersel) == 0: return atomgrp
        else:
            # Generate a selection for each selection string
            #atomselections = [atomgrp]
            for sel in othersel:
                atomgrp = atomgrp + Selection.Parser.parse(sel).apply(self)
                #atomselections.append(Selection.Parser.parse(sel).apply(self))
            #return tuple(atomselections)
            return atomgrp

    def write(self,filename=None,format="pdb",filenamefmt="%(trjname)s_%(frame)d"):
        """Write AtomGroup to a file.

        AG.write(filename,format='pdb')

        EXPERIMENTAL.
        Only a primitive PDB format and standard CRD is working.

        filename      None: create TRJNAME_FRAME.FORMAT from filenamefmt
        format        pdb, crd; can also be supplied as part of filename
        filenamefmt   format string for default filename; use 'trjname' and 'frame'
        """
        import util
        import os.path

        trj = self.universe.trajectory    # unified trajectory API
        frame = trj.ts.frame

        if filename is None:
            trjname,ext = os.path.splitext(os.path.basename(trj.filename))
            filename = filenamefmt % vars()
        filename = util.filename(filename,ext=format,keep=True)
        format = os.path.splitext(filename)[1][1:]  # strip initial dot!
        try:
            import MDAnalysis.coordinates
            FrameWriter = MDAnalysis.coordinates._frame_writers[format]
        except KeyError:
            raise NotImplementedError("Writing as %r is not implemented; only %r will work." \
                                      % (format, MDAnalysis.coordinates._frame_writers.keys()))
        writer = FrameWriter(filename)
        writer.write(self)         # wants a atomgroup

    # TODO: This is _almost_ the same code as write() --- should unify!
    def write_selection(self,filename=None,format="vmd",filenamefmt="%(trjname)s_%(frame)d",
                        **kwargs):
        """Write AtomGroup selection to a file to be used in another programme.

        :Keywords:
          *filename*
                ``None``: create TRJNAME_FRAME.FORMAT from *filenamefmt*
          *format*
                output file format: VMD (tcl), PyMol (pml), Gromacs (ndx), CHARMM (str);
                can also be supplied as the filename extension. Case insensitive. [vmd]
          *filenamefmt*
                format string for default filename; use '%(trjname)s' and '%(frame)s'
                placeholders; the extension is set according to the *format*
                ["%(trjname)s_%(frame)d"]
          *kwargs*
                additional keywords are passed on to the appropriate
                :class:`~MDAnalysis.selections.base.SelectionWriter`
        """
        import util
        import os.path
        import MDAnalysis.selections

        SelectionWriter = MDAnalysis.selections.get_writer(filename, format)

        trj = self.universe.trajectory    # unified trajectory API
        frame = trj.ts.frame

        # get actual extension from the static class attribute
        extension = SelectionWriter.ext

        if filename is None:
            trjname,ext = os.path.splitext(os.path.basename(trj.filename))
            filename = filenamefmt % vars()
        filename = util.filename(filename,ext=extension,keep=True)

        writer = SelectionWriter(filename, **kwargs)
        writer.write(self)         # wants a atomgroup

    # properties
    @property
    def dimensions(self):
        """Dimensions of the Universe to which the group belongs, at the current time step."""
        if self.universe is not None:
            return self.universe.dimensions
        else:
            raise AttributeError("This AtomGroup does not belong to a Universe with a dimension.")

    @property
    def bfactors(self):
        """B-factors of the AtomGroup"""
        if self.universe is not None:
            return self.universe.bfactors[self.indices()]
        else:
            raise AttributeError("This AtomGroup does not belong to a Universe.")
        

class Residue(AtomGroup):
    """A group of atoms corresponding to a residue

    Data: type, name

    Methods:
        r = Residue()
        r['name'] or r[id] - returns the atom corresponding to that name
        r.name
    """
    __cache = {}
    def __init__(self, name, id, atoms):
        super(Residue, self).__init__(atoms)
        self.name = name
        self.id = id
        self.segment = None
        for i, a in enumerate(atoms):
            a.id = i
            a.residue = self
        # Should I cache the positions of atoms within a residue?
        if not Residue.__cache.has_key(name):
            Residue.__cache[name] = dict([(a.name, i) for i, a in enumerate(self.atoms)])
    def __getitem__(self, item):
        if (type(item) == int) or (type(item) == slice):
            return self.atoms[item]
        else: return self.__getattr__(item)
    def __getattr__(self, name):
        # There can only be one atom with a certain name
        # Use the cache
        #for atom in self.atoms:
        #    if (name == atom.name): return atom
        try:
            index = Residue.__cache[self.name][name]
            return self.atoms[index]
        except KeyError: raise SelectionError("No atom in residue "+self.name+" with name "+name)
    def __repr__(self):
        return '<'+self.__class__.__name__+' '+repr(self.name)+', '+repr(self.id)+'>'

class ResidueGroup(AtomGroup):
    """A group of residues

    Data: residues

    Methods:
       rg = ResidueGroup()
    """
    def __init__(self, residues):
        self._residues = residues
        atoms = []
        for res in residues:
            atoms.extend(res.atoms)
        super(ResidueGroup, self).__init__(atoms)
    def __iter__(self):
        return iter(self._residues)
    def __len__(self):
        return len(self._residues)
    def __getitem__(self, item):
        if (type(item) == int) or (type(item) == slice):
            return self._residues[item]
        else: raise NotImplementedError()
    def __getattr__(self, attr):
        atomlist = [atom for atom in self.atoms if atom.name == attr]
        return AtomGroup(atomlist)
    def __repr__(self):
        return '<'+self.__class__.__name__+' '+repr(self._residues)+'>'

class Segment(ResidueGroup):
    """A group of residues corresponding to one segment of the PSF

    Data: name

    Methods:
       s = Segment()
    """
    def __init__(self, name, residues):
        super(Segment, self).__init__(residues)
        self.name = name
        for res in self._residues:
            res.segment = self
            for atom in res:
                atom.segment = self
    def __getattr__(self, attr):
        if attr[0] == 'r':
            resnum = int(attr[1:])
            return self[resnum]
        else:
            # There can be multiple residues with the same name
            r = []
            for res in self._residues:
                if (res.name == attr): r.append(res)
            if (len(r) == 0): return super(Segment, self).__getattr__(attr)
            elif (len(r) == 1): return r[0]
            else: return ResidueGroup(r)
    def __repr__(self):
        return '<'+self.__class__.__name__+' '+repr(self.name)+'>'

class Universe(object):
    """The MDAnalysis Universe contains all the information describing the system.

    Built from a psf, pdb or gro file.

    Attributes: 
       - bonds, angles, dihedrals, impropers, donors, acceptors [TODO]
       - :attr:`Universe.coord`: current coordinate array for all the atoms in the universe
       - :attr:`Universe.trajectory`: currently loaded trajectory reader; 
          :attr:`Universe.trajectory.ts` is the current time step
       - :attr:`Universe.dimensions`: current system dimensions (simulation unit cell, if 
         set in the trajectory)

    Methods::
       m = Universe(psffile, dcdfile)             # read system from file(s)          
       m = Universe(psffile, pdbfilename=pdbfile) # read coordinates from PDB file
       m = Universe(pdbfile)                      # read atoms and coordinates from PDB
       m = Universe(pdbfile, xtcfile)             # read system from file(s)          
       m.load_new_dcd(dcdfilename)                # read from a new dcd file
       m.selectAtoms(...)                         # select atoms using similar selection string as charmm

    .. Note:: Only a single-frame PDB file is supported; use DCDs or XTC/TRR for
              trajectories. If a PDB is used instead of a PSF then
              neither mass nor charge are correct, and bonds are not available.
    """
    def __init__(self, psffilename, dcdfilename=None, pdbfilename=None):
        """Initialize the central MDAnalysis Universe object.

        :Arguments:
          *psffilename*
             A Charmm/XPLOR PSF topology file or a PDB file; used to define the
             list of atoms. If the file includes bond information, partial
             charges, atom masses, ... then these data will be available to
             MDAnalysis. A "structure" file (PSF or PDB, in the sense of a
             topology) is always required.
          *dcdfilename*
             A CHARMM DCD trajectory or Gromacs XTC/TRR or a PDB; will provide coordinates.
        
        This routine tries to do the right thing: 
          1. If a pdb file is provided instead of a psf and neither a dcd nor a
             pdb structure then the coordinates are taken from the first pdb
             file. Thus you can load a functional universe with ::

                u = Universe('1ake.pdb')
   
          2. If a pdb coordinate file is provided in the *dcdfilename* argument
             then it is silently opened as a pdb file.

          3. If only a psf file is provided one will have to load coordinates
             manually using :meth:`Universe.load_new_dcd` or
             :meth:`Universe.load_new_pdb`.
        """
        import MDAnalysis.topology.core

        # managed attribute holding TRJReader (the Universe.trajectory
        # attribute is also aliase as Universe.<EXT> where <EXT> is the
        # trajectory format type (i.e. the extension))
        self.__trajectory = None

        # for 0.7 we will rename psffile -> topologyfile, and dcdfile -> coordinatefile
        if not pdbfilename is None:
            warnings.warn("Usage of 'pdbfilename=PDB' is deprecated and will be removed. Just use the PDB "
                          "as the second argument of Universe().", category=DeprecationWarning)
        topologyfile = psffilename
        coordinatefile = dcdfilename or pdbfilename

        # build the topology (or at least a list of atoms)
        try:
            parser = MDAnalysis.topology.core.get_parser_for(topologyfile)
            struc = parser(topologyfile)
        except TypeError, err:
            raise ValueError("Failed to build a topology from either a psf, pdb or gro (%s)" % err)

        self.filename = topologyfile
        self._psf = struc
        #for data in struc.keys():
        #    setattr(self, data, struc[data])
        self.atoms = AtomGroup(struc["_atoms"])
        # XXX: add H-bond information here if available from psf (or other sources)
        # 
        segments = MDAnalysis.topology.core.build_segments(self.atoms)
        # Because of weird python rules, attribute names cannot start with a digit
        for seg in segments.keys():
            if seg[0].isdigit():
                newsegname = 's'+seg
                segments[newsegname] = segments[seg]
                del segments[seg]
        self.__dict__.update(segments)

        #MDAnalysis.topology.core.build_bondlists(self.atoms, self._bonds)
        # Let atoms access the universe
        for a in self.atoms: a.universe = self

        # Load coordinates; distinguish file format by extension
        if coordinatefile is None:
            # hack for pdb/gro - only
            coordinatefile = topologyfile
        self.load_new(coordinatefile)

    def load_new(self, filename):
        """Load coordinates from *filename*, using the suffix to detect file format."""
        if filename is None:
            return

        import MDAnalysis.coordinates.core

        try:
            TRJReader = MDAnalysis.coordinates.core.get_reader_for(filename)
        except TypeError, err:
            # only warn because in the past it was possible to build a topology-only
            # Universe and populate it later with coordinates.
            warnings.warn("Universe.load_new() cannot find an appropriate coordinate reader "
                          "for file %r. Universe will be built without coordinates: "
                          "use load_new() later with an appropriate trajectory." % filename,
                          category=UserWarning)
            return None, None
        self.trajectory = TRJReader(filename)    # unified trajectory API
        trjtype = self.trajectory.format.lower() # trjtype is always lower case (see coordinates._frame_readers)
        self.__dict__[trjtype] = self.trajectory # legacy (deprecated)
        # Make sure that they both have the same number of atoms
        if (self.trajectory.numatoms != self.atoms.numberOfAtoms()):
            raise ValueError("The topology and %s trajectory files  don't have the same number of atoms!" % trjtype)
        # hack for PDB
        if trjtype == "pdb":
            # add B-factor to atoms
            for a, pdbatom in zip(self.atoms,self.trajectory.pdb.get_atoms()):
                a.bfactor = pdbatom.get_bfactor()  ## does this work with mmLib??
        return filename, trjtype
            
    def selectAtoms(self, sel, *othersel):
        """Select atoms using a CHARMM selection string. 

        Returns an AtomGroup with atoms sorted according to their index in the
        psf (this is to ensure that there aren't any duplicates, which can
        happen with complicated selections).

        .. Note:: you can group subselections using parentheses, but
                  you need to put spaces around the parentheses due to
                  the way the parser is implemented.

        Example:
           >>> universe.selectAtoms("segid DMPC and not ( name H* or name O* )")
           <AtomGroup with 3420 atoms>

        Here's a list of all keywords; keywords are CASE SENSITIVE:

        **Simple selections**
        .. ------------------

           protein, backbone, nucleic, nucleicbackbone
              selects all atoms that belong to a standard set of residues, may
              not work for esoteric residues

           segid, resid, resname, name, type
              single argument describing selection, resid can take a range of
              numbers separated by a colon (inclusive) ie "segid DMPC", "resname
              LYS", "name CA", "resid 1:5"

           atom 
              a selector for a single atom consisting of segid resid atomname ie
              "DMPC 1 C2" selects the C2 carbon of the first residue of the DMPC
              segment

        **Boolean**
        .. --------

           not
              all atoms not in the selection; ie "not protein" selects all atoms
              that aren't part of a protein

           and, or
              combine two selections according to the rules of boolean algebra ie
              "protein and not ( resname ALA or resname LYS )" selects all atoms
              that belong to a protein, but are not in a lysine or alanine
              residue


        **Geometric**
        .. ----------

           around
              selects all atoms a certain cutoff away from another selection; ie
              "around 3.5 protein" selects all atoms not belonging to protein
              that are within 3.5 Angstroms from the protein point::

                around distance selection

           point
              selects all atoms within a sphere of given radius around point::

                point distance x y z

              Make sure coordinate is separated by spaces, i.e. "point 3.5
              5.0 5.0 5.0" selects all atoms within 3.5 Angstroms of the
              coordinate (5.0, 5.0, 5.0)

           prop
              selects atoms based on position, using x, y, or z coordinate
              supports the abs keyword (for absolute value) and the following
              operators: >, <, >=, <=, ==, != ie "prop z >= 5.0" selects all
              atoms with z coordinate greater than 5.0 "prop abs z <= 5.0"
              selects all atoms within -5.0 <= z <= 5.0

           Periodicity is only taken into account with the 'distance matrix'
           distance functions via a minimum image convention (and this only works
           for rectangular simulation cells). If this behavior is required, set
           these flags::

              MDAnalysis.core.flags['use_periodic_selections'] = True   # default
              MDAnalysis.core.flags['use_KDTree_routines'] = False


        **Connectivity**
        .. -------------

           byres
              selects all atoms that are in the same segment and residue as
              selection ie specify the subselection after the byres keyword


        **Index**
        .. ------

           bynum
              selects all atoms within a range of (1-based) inclusive indices ie
              "bynum 1" selects the first atom in the universe "bynum 5:10"
              selects atoms 5 through 10 inclusive
        """
        import Selection
        atomgrp = Selection.Parser.parse(sel).apply(self)
        if len(othersel) == 0: return atomgrp
        else:
            # Generate a selection for each selection string
            #atomselections = [atomgrp]
            for sel in othersel:
                atomgrp = atomgrp + Selection.Parser.parse(sel).apply(self)
                #atomselections.append(Selection.Parser.parse(sel).apply(self))
            #return tuple(atomselections)
            return atomgrp

    def __repr__(self):
        return '<'+self.__class__.__name__+' with '+repr(len(self.atoms))+' atoms>'

    # Properties
    @property
    def dimensions(self):
        """Current dimensions of the unitcell"""
        return self.coord.dimensions
    
    @property
    def coord(self):
        """Reference to current timestep and coordinates of universe.

        The raw trajectory coordinates are :attr:`Universe.coord._pos`,
        represented as a :attr:`numpy.float32` array.

        Because :attr:`coord` is a reference, it changes its contents while one
        is stepping through the trajectory.

        .. Note:: In order to access the coordinates it is probably better to
                  use the :meth:`AtomGroup.coordinates` method; for instance,
                  all coordinates of the Universe as a numpy array:
                  :meth:`Universe.atoms.coordinates`.

        .. SeeAlso:: :class:`MDAnalysis.coordinates.DCD.Timestep`
        """
        return self.trajectory.ts

    def trajectory():
        doc = "Reference to trajectory reader object containing trajectory data"
        def fget(self):
            if not self.__trajectory == None: return self.__trajectory
            else: raise AttributeError("No trajectory loaded into Universe")
        def fset(self, value):
            del self.__trajectory  # guarantees that files are closed (?)
            self.__trajectory = value
        return locals()
    trajectory = property(**trajectory())

def asUniverse(*args, **kwargs):
    """Return a universe from the input arguments.

    1. If the first argument is a universe, just return it::
 
       as_universe(universe) --> universe

    2. Otherwise try to build a universe from the first or the first
       and second argument::
 
       asUniverse(PDB, **kwargs) --> Universe(PDB, **kwargs)
       asUniverse(PSF, DCD, **kwargs) --> Universe(PSF, DCD, **kwargs)
       asUniverse(*args, **kwargs) --> Universe(*args, **kwargs)

    :Returns: an instance of :class:`~MDAnalaysis.AtomGroup.Universe`
    """
    if len(args) == 0:
        raise TypeError("asUniverse() takes at least one argument (%d given)" % len(args))
    elif len(args) == 1 and type(args[0]) is Universe:
        return args[0]
    return Universe(*args, **kwargs)
        
