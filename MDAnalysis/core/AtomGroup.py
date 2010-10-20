# $Id$
"""
:mod:`MDAnalysis.core.AtomGroup` --- Basic Organization of the library
======================================================================

The most important data structure in MDAnalysis is the
:class:`AtomGroup`, which contains :class:`Atom` instances.

A :class:`Universe` is the user-visible entry point and collects all
information needed to analyze a structure or a whole trajectory.

Segments and Residues are a way to refer to a collection of atoms. By
convention, a :class:`Residue` is a single amino acid, or a water
molecule, ion, or ligand. A :class:`Segment` is a collection of
residues such as a whole protein or a chain in a protein or all the
water in the system. 

Class Hierarchy
---------------

A :class:`Universe` contains Segments, which contain Residues, which
contain Atoms; all containers are derived from :class:`AtomGroup`, and
thus one can always analyze them as a collection of atoms, independent
of the hierarchical level.

    Segment > Residue > Atom

Depending on the use case, it can be more convenient to access data
on, for instance, the basis of residues than atoms, or to write out
individual chains (segments) of a protein. MDAnalysis simply provides
access to these groupings and keeps track of where an atom
belongs. Each object provides three attributes (:attr:`atoms`,
:attr:`residues` or :attr:`~Atom.residue`, :attr:`segments` or
:attr:`~Atom.segment`) that give access to the tiers in the hierarchy
that the object belongs to.


Classes and functions
---------------------

.. autoclass:: Universe
   :members:
.. autoclass:: AtomGroup
   :members:
.. autoclass:: Atom
   :members:
.. autoclass:: Residue
.. autoclass:: ResidueGroup
.. autoclass:: Segment
.. autoclass:: SegmentGroup

.. autofunction:: asUniverse
.. autoexception:: SelectionError
.. autoexception:: SelectionWarning
"""
import warnings

import numpy
from MDAnalysis import SelectionError, NoDataError, SelectionWarning

class Atom(object):
    """A single atom definition

    a = Atom()

    :Data: 
        number
          atom number
        segid
          name of the segment
        resid
           residue number
        resname
           residue name
        residue
           :class:`Residue` object containing the atoms
        name
        type
        mass
        charge
        :attr:`~Atom.pos`
           The current position (as a numpy array) of this atom
    """
    __slots__ = ("number", "id", "name", "type", "resname", "resid", "segid", 
                 "mass", "charge", "residue", "segment", "bonds", "__universe",)

    def __init__(self, number, name, type, resname, resid, segid, mass, charge,
                 residue=None, segment=None):
        self.number = number
        self.name = name
        self.type = type
        self.resname = resname
        self.resid = resid
        self.residue = residue  # typically patched in later
        self.segid = segid
        self.segment = segment  # typically patched in later
        self.mass = mass
        self.charge = charge
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

    @property
    def pos(self):
        """Current cartesian coordinates of the atom."""
        return self.universe.coord[self.number] # PDB numbering starts at 0

    def centroid(self):
        """The centroid of an atom is its position, :attr:`Atom.pos`."""
        # centroid exists for compatibility with AtomGroup
        return self.pos

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
    """A group of atoms.

    ag = universe.selectAtoms(atom-list)

    The AtomGroup contains a list of atoms; typically, a AtomGroup is generated
    from a selection. It is build from any list-like collection of
    :class:`Atom` instances.

    An AtomGroup can be indexed and sliced like a list ::

       ag[0], ag[-1]   

    will return the first and the last :class:`Atom` in the group
    whereas the slice

       ag[0:6:2]

    returns every second element, corresponding to indices 0, 2, and 4.

    It also supports "advanced slicing" when the argument is a
    :class:`numpy.ndarray` or a :class:`list`::

       aslice = [0, 3, -1, 10, 3]
       ag[aslice]

    will return a new :class:`AtomGroup` containing (ag[0], ag[3], ag[-1],
    ag[10], ag[3]).

    .. Note:: AtomGroups originating from a selection are sorted and
       duplicate elements are removed. This is not true for AtomGroups
       produced by slicing. Thus slicing can be used when the order of
       atoms is crucial (for instance, in order to define angles or
       dihedrals).

    Atoms can also be accessed in a Pythonic fashion by using the atom name as
    an attribute. For instance, ::
      
       ag.CA

    will provide a :class:`AtomGroup` of all CA atoms in the group. 

    .. Note:: The name-attribute access to atoms is mainly meant for quick
       interactive work. Thus it either returns a single :class:`Atom` if there
       is only one matching atom, *or* a new :class:`AtomGroup` for multiple
       matches. This makes it difficult to use the feature consistently in
       scripts but is much better for interactive work.

    :Data:
        atoms
            The AtomGroup itself; provided for unified access across the hierarchy. 
            Use :attr:`_atoms` if you really only need a list.
        _atoms
            A list of references to the corresponding atoms in :attr:`Universe.atoms`;
            AtomGroups are immutable, i.e. this list cannot be changed.
        residues
            A :class:`ResidueGroup` of all residues that contain atoms in this group.
        segments
            A :class:`SegmentGroup` of all segments that contain atoms in this group.
        ts
            A :class:`~MDAnalysis.coordinates.base.Timestep` instance, which can
            be passed to a trjectory writer.
        *atom-name*
            Auto-generated attribute for each atom name encountered in the group.
    """
    def __init__(self, atoms):
        if len(atoms) < 1: raise NoDataError("No atoms defined for AtomGroup")
        # __atoms property is effectively readonly        
        # check that atoms is indexable:
        try:
            atoms[0]
            self.__atoms = atoms      
        except TypeError:
            self.__atoms = list(atoms)
        # sanity check
        if not isinstance(self.__atoms[0], Atom):
            raise TypeError("atoms must be a Atom or a list of Atoms.")
        # If the number of atoms is very large, create a dictionary cache for lookup
        if len(atoms) > 10000:
            self._atom_cache = dict(((x,None) for x in self.__atoms))
        # managed timestep object
        self.__ts = None

    # AtomGroup.atoms is guaranteed to be a AtomGroupTo; keeps a consistent API
    # between AtomGroup, Residue, ResidueGroup, Segment; access the list as
    # _atoms (although atoms supports all list-like operations, too).
    @property
    def atoms(self):
        """AtomGroup of all atoms in this group"""
        # Cannot just return self because fails with inheritance from AtomGroup
        if type(self) == AtomGroup:            
            return self
        return AtomGroup(self.__atoms)

    @property
    def _atoms(self):
        """a immutable list of references to the atoms in the group"""
        return self.__atoms
    
    # Universe pointer is important for Selections to work on groups
    @property
    def universe(self):
        """The universe to which the atoms belong (read-only)."""
        try:
            return self._atoms[0].universe
        except AttributeError:
            return None

    def __len__(self):
        return self.numberOfAtoms()
    def __getitem__(self, item):
        """Return Atom (index) or AtomGroup (slicing)"""
        # consistent with the way list indexing/slicing behaves:
        if numpy.dtype(type(item)) == numpy.dtype(int):
            return self._atoms[item]
        elif type(item) == slice:
            return AtomGroup(self._atoms[item])
        elif isinstance(item, (numpy.ndarray, list)):
            # advanced slicing, requires array or list
            return AtomGroup([self._atoms[i] for i in item])
        else:
            return super(AtomGroup, self).__getitem__(item)
    def __getattr__(self, name):
        # There can be more than one atom with the same name
        atomlist = [atom for atom in self._atoms if name == atom.name]
        if len(atomlist) == 0: raise SelectionError("No atoms with name "+name)
        elif len(atomlist) == 1: return atomlist[0]  # XXX: keep this, makes more sense for names
        else: return AtomGroup(atomlist)             # XXX: but inconsistent (see residues and Issue 47)
    def __iter__(self):
        return iter(self._atoms)
    def __contains__(self, other):
        if hasattr(self, "_atom_cache"):
            return other in self._atom_cache
        else: 
            return other in self._atoms
    def __add__(self, other):
        if not (isinstance(other, Atom) or isinstance(other, AtomGroup)):
            raise TypeError('Can only concatenate AtomGroup (not "'+repr(other.__class__.__name__)+'") to AtomGroup')
        if isinstance(other, AtomGroup): 
            return AtomGroup(self._atoms + other._atoms)
        else: 
            return AtomGroup(self._atoms+[other])
    def __repr__(self):
        return '<'+self.__class__.__name__+' with '+repr(self.numberOfAtoms())+' atoms>'
    def numberOfAtoms(self):
        return len(self._atoms)
    def numberOfResidues(self):
        return len(self.residues)
    def numberOfSegments(self):
        return len(self.segments)
    def indices(self):
        if not hasattr(self,'_cached_indices'):
            self._cached_indices = numpy.array([atom.number for atom in self._atoms])
        return self._cached_indices
    @property
    def residues(self):
        """Read-only list of :class:`Residue` objects."""
        if not hasattr(self,'_cached_residues'):
            residues = []
            current_residue = None
            for atom in self._atoms:
                if atom.residue != current_residue:
                    residues.append(atom.residue)
                current_residue = atom.residue
            self._cached_residues = ResidueGroup(residues)
        return self._cached_residues
    def resids(self):
        """Returns a list of residue numbers."""
        return [r.id for r in self.residues]
    def resnames(self):
        """Returns a list of residue names."""
        return [r.name for r in self.residues]
    @property
    def segments(self):
        """Read-only list of :class:`Segment` objects."""
        if not hasattr(self,'_cached_segments'):
            segments = []
            current_segment = None
            for atom in self._atoms:
                if atom.segment != current_segment:
                    segments.append(atom.segment)
                current_segment = atom.segment
            self._cached_segments = SegmentGroup(segments)
        return self._cached_segments
    def segids(self):
        """Returns a list of segment ids (=segment names)."""
        return [s.name for s in self.segments]
    def masses(self):
        """Array of atomic masses (as defined in the topology)"""
        if not hasattr(self, "_cached_masses"):
            self._cached_masses = numpy.array([atom.mass for atom in self._atoms])
        return self._cached_masses
    def totalMass(self):
        """Total mass of the selection (masses are taken from the topology or guessed)."""
        return numpy.sum(self.masses(), axis=0)
    def charges(self):
        """Array of partial charges of the atoms (as defined in the topology)"""
        return numpy.array([atom.charge for atom in self._atoms])
    def totalCharge(self):
        """Sum of all partial charges (must be defined in topology)."""
        return numpy.sum(self.charges(), axis=0)
    def centerOfGeometry(self):
        """Center of geometry (also known as centroid) of the selection."""
        return numpy.sum(self.coordinates(), axis=0)/self.numberOfAtoms()
    centroid = centerOfGeometry
    def centerOfMass(self):
        """Center of mass of the selection."""
        return numpy.sum(self.coordinates()*self.masses()[:,numpy.newaxis],axis=0)/self.totalMass()
    def radiusOfGyration(self):
        """Radius of gyration."""
        masses = self.masses()
        recenteredpos = self.coordinates() - self.centerOfMass()
        rog_sq = numpy.sum(masses*numpy.sum(numpy.power(recenteredpos, 2), axis=1))/self.totalMass()
        return numpy.sqrt(rog_sq)
    def momentOfInertia(self):
        """Tensor of inertia as 3x3 NumPy array."""
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

    def principalAxes(self):
        """Calculate the principal axes from the moment of inertia.

        e1,e2,e3 = AtomGroup.principalAxes()

        The eigenvectors are sorted by eigenvalue, i.e. the first one
        corresponds to the highest eigenvalue and is thus the first principal axes.

        :Returns: numpy.array ``v`` with ``v[0]`` as first, ``v[1]`` as second,
                  and ``v[2]`` as third eigenvector.
        """
        from numpy.linalg import eig
        eigenval, eigenvec = eig(self.momentOfInertia())
        # Sort
        indices = numpy.argsort(eigenval)
        # Return transposed in more logical form. See Issue 33.
        return eigenvec[:,indices].T

    def coordinates(self, ts=None, copy=False, dtype=numpy.float32):
        """NumPy array of the coordinates."""
        if ts == None:
            ts = self.universe.trajectory.ts
        return numpy.array(ts[self.indices()], copy=copy, dtype=dtype)

    def translate(self, t):
        """Apply translation vector *t* to the selection's coordinates.

        AtomGroup.translate(t)

        AtomGroup.translate((A, B))

        The method applies a translation to the AtomGroup from current
        coordinates x to new coordinates x':

            x' = x + t

        The translation can also be given as a tuple of two MDAnalysis objects
        such as two selections (selA, selB), i.e. two :class:`AtomGroup`s, or
        two :class:`Atom`s. The translation vector is computed as the
        difference between the centers of geometry (centroid) of B and A:

        t = B.centroid() - A.centroid()
        """
        try:
            sel1,sel2 = t
            x1,x2 = sel1.centroid(), sel2.centroid()
            vector = x2 - x1
        except (ValueError, AttributeError):
            vector = numpy.asarray(t)
        # changes the coordinates (in place)
        self.universe.trajectory.ts._pos[self.indices()] += vector
        return vector

    def rotate(self, R):
        """Apply a rotation matrix *R* to the selection's coordinates.

        AtomGroup.rotate(R)

        *R* is a 3x3 orthogonal matrix that transforms x --> x':

            x' = R.x
        """
        R = numpy.matrix(R, copy=False, dtype=numpy.float32)
        # changes the coordinates (in place)
        x = self.universe.trajectory.ts._pos
        idx = self.indices()
        x[idx] = x[idx] * R.T     # R.T acts to the left & is broadcasted N times.
        return R

    def rotateby(self, angle, axis, point=None):
        """Apply a rotation to the selection's coordinates.

        AtomGroup.rotateby(angle,axis[,point])

        The transformation from current coordinates x to new coordinates x' is

          x' = R.(x-p) + p

        where R is the rotation by *angle* around the *axis* going through
        *point* p.

        :Arguments:
          *angle*
             rotation angle in degrees
          *axis*
             rotation axis vector, a 3-tuple, list, or array, or a 2-tuple of
             two MDAnalysis objects from which the axis is calculated as the
             vector from the first to the second center of geometry.
          *point*
             point on the rotation axis; by default (``None``) the center of
             geometry of the selection is chosen, or, if *axis* is a tuple of
             selections, it defaults to the first point of the axis. *point*
             can be a 3-tuple, list, or array or a MDAnalysis object (in which
             case its :meth:`centroid` is used).

        :Returns: The 4x4 matrix which consists of the rotation matrix M[:3,:3]
                  and the translation vectort M[:3,3].
        """
        from transformations import rotation_matrix
        alpha = numpy.radians(angle)
        try:
            sel1,sel2 = axis
            x1,x2 = sel1.centroid(), sel2.centroid()
            v = x2 - x1
            n = v/numpy.linalg.norm(v)
            if point is None:
                point = x1
        except (ValueError, AttributeError):
            n = numpy.asarray(axis)
        if point is None:
            p = self.centroid()
        else:
            try:
                p = point.centroid()
            except AttributeError:
                p = numpy.asarray(point)
        M = rotation_matrix(alpha, n, point=p)
        R = M[:3,:3]
        t = M[:3, 3]
        x = self.universe.trajectory.ts._pos
        idx = self.indices() 
        x[idx]  = numpy.dot(x[idx], R.T)
        x[idx] += t
        return M

    def selectAtoms(self, sel, *othersel):
        """Selection of atoms using the MDAnalysis selection syntax.

        AtomGroup.selectAtoms(selection[,selection[,...]])

        .. SeeAlso:: :meth:`Universe.selectAtoms`
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

    def write(self,filename=None,format="PDB",filenamefmt="%(trjname)s_%(frame)d"):
        """Write AtomGroup to a file.

        AtomGroup.write(filename[,format])
        
        :Keywords:
          *filename*
               ``None``: create TRJNAME_FRAME.FORMAT from filenamefmt [``None``]
          *format*
                PDB, CRD, GRO; case-insensitive and can also be supplied as 
                the filename extension [PDB]
          *filenamefmt*   
                format string for default filename; use substitution tokens
                'trjname' and 'frame' ["%(trjname)s_%(frame)d"]
        """
        import util
        import os.path
        import MDAnalysis.coordinates

        trj = self.universe.trajectory    # unified trajectory API
        frame = trj.ts.frame

        if filename is None:
            trjname,ext = os.path.splitext(os.path.basename(trj.filename))
            filename = filenamefmt % vars()
        filename = util.filename(filename,ext=format.lower(),keep=True)
        framewriter = MDAnalysis.coordinates.writer(filename)
        framewriter.write(self)         # wants a atomgroup

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

    @property
    def ts(self):
        """Temporary Timestep that contains the selection coordinates.

        If the :attr:`AtomGroup.ts`
        :class:`~MDAnalysis.coordinates.base.Timestep` object is modified then
        these modifications will be present until the frame number changes
        (which typically happens when the underlying trajectory frame changes).

        It is not possible to assign a new
        :class:`~MDAnalysis.coordinates.base.Timestep` to the
        :attr:`AtomGroup.ts` attribute; change attributes of the object.
        """
        trj_ts = self.universe.trajectory.ts  # original time step
        if self.__ts is None or self.__ts.frame != trj_ts.frame:
            # create a timestep of same type as the underlying trajectory
            ts = trj_ts.__class__(self.coordinates())
            ts._unitcell = trj_ts._unitcell.copy()
            ts.frame = trj_ts.frame
            self.__ts = ts
        return self.__ts

class Residue(AtomGroup):
    """A group of atoms corresponding to a residue.

    Pythonic access to atoms:
      - Using a atom name as attribute returns the matching atom (a
        :class:`Atom` instance), i.e. ``r.name``. Example:
          >>> from MDAnalysis.tests.datafiles import PSF,DCD
          >>> u = Universe(PSF,DCD)
          >>> print(u.s4AKE.r1.CA)  # C-alpha of M1
          < Atom 5: name 'CA' of type '22' of resname 'MET', resid 1 and segid '4AKE'>
      - ``r['name']`` or ``r[id]`` - returns the atom corresponding to that name

    :Data: 
      :attr:`Residue.name`
        Three letter residue name.
      :attr:`Residue.id`
        Numeric (integer) resid, taken from the topology.

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
            Residue.__cache[name] = dict([(a.name, i) for i, a in enumerate(self._atoms)])
    def phi_selection(self):
        """AtomGroup corresponding to the phi protein backbone dihedral C'-N-CA-C.

        :Returns: 4-atom selection in the correct order.  If no C'
                  found in the previous residue (by resid) then this
                  method returns ``None``.
        """
        try:
            return self.universe.selectAtoms('resid %d and name C' % (self.id - 1)) \
                + self.N + self.CA + self.C
        except (SelectionError, NoDataError):
            return None

    def psi_selection(self):
        """AtomGroup corresponding to the psi protein backbone dihedral N-CA-C-N'.

        :Returns: 4-atom selection in the correct order.  If no N'
                  found in the following residue (by resid) then this
                  method returns ``None``.
        """
        try:
            return self.N + self.CA + self.C + \
                self.universe.selectAtoms('resid %d and name N' % (self.id + 1))
        except (SelectionError, NoDataError):
            return None

    def omega_selection(self):
        """AtomGroup corresponding to the omega protein backbone dihedral CA-C-N'-CA'.

        omega described the -C-N- peptide bond. Typically, it is trans
        (180 degrees) although cis-bonds (0 degrees) are also
        occasionally observed (especially near Proline).

        :Returns: 4-atom selection in the correct order.  If no C'
                  found in the previous residue (by resid) then this
                  method returns ``None``.
        """
        nextres = self.id + 1
        try:
            return self.CA + self.C +\
                self.universe.selectAtoms('resid %d and name N' % nextres, 'resid %d and name CA' % nextres)
        except (SelectionError, NoDataError):
            return None

    def __getitem__(self, item):
        if (type(item) == int) or (type(item) == slice):
            return self._atoms[item]
        else: return self.__getattr__(item)
    def __getattr__(self, name):
        # There can only be one atom with a certain name
        # Use the cache
        #for atom in self.atoms:
        #    if (name == atom.name): return atom
        try:
            index = Residue.__cache[self.name][name]
            return self._atoms[index]
        except KeyError: raise SelectionError("No atom in residue "+self.name+" with name "+name)
    def __repr__(self):
        return '<'+self.__class__.__name__+' '+repr(self.name)+', '+repr(self.id)+'>'

class ResidueGroup(AtomGroup):
    """A group of residues.

    Pythonic access to atoms:
      - Using a atom name as attribute returns a list of all atoms (a
        :class:`AtomGroup`) of the same name. Example:
          >>> from MDAnalysis.tests.datafiles import PSF,DCD
          >>> u = Universe(PSF,DCD)
          >>> print(u.s4AKE.MET.CA)  # C-alpha of all Met
          <AtomGroup with 6 atoms>

    :Data: :attr:`ResidueGroup._residues`

    """
    def __init__(self, residues):
        """Initialize the ResidueGroup with a list of :class:`Residue` instances."""
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
    """A group of residues corresponding to one segment of the topology.

    Pythonic access to residues:
      - The attribute rN returns the N-th residue :class:`Residue` of the 
        segment (numbering starts at N=1). Example:
          >>> from MDAnalysis.tests.datafiles import PSF,DCD
          >>> u = Universe(PSF,DCD)
          >>> print(u.s4AKE.r1)
          <Residue 'MET', 1>
      - Using a residue name as attribute returns a list of all residues (a
        :class:`ResidueGroup`) of the same name. Example:
          >>> from MDAnalysis.tests.datafiles import PSF,DCD
          >>> u = Universe(PSF,DCD)
          >>> print(u.s4AKE.CYS)
          <ResidueGroup [<Residue 'CYS', 77>]>
          >>> print(u.s4AKE.MET)
          <ResidueGroup [<Residue 'MET', 1>, <Residue 'MET', 21>, <Residue 'MET', 34>, <Residue 'MET', 53>, <Residue 'MET', 96>, <Residue 'MET', 174>]>

    :Data: :attr:`Segment.name` is the segid from the topology or the
           chain identifier when loaded from a PDB
    """
    def __init__(self, name, residues):
        """Initialize a Segment with segid *name* from a list of :class:`Residue` instances."""
        super(Segment, self).__init__(residues)
        self.name = name
        for res in self._residues:
            res.segment = self
            for atom in res:
                atom.segment = self
    def __getattr__(self, attr):
        if attr[0] == 'r':
            resnum = int(attr[1:]) - 1   # 1-based for the user, 0-based internally
            return self[resnum]
        else:
            # There can be multiple residues with the same name
            r = []
            for res in self._residues:
                if (res.name == attr): r.append(res)
            if (len(r) == 0): return super(Segment, self).__getattr__(attr)
            # elif (len(r) == 1): return r[0]  ## creates unexpected behaviour (Issue 47)
            else: return ResidueGroup(r)
    def __repr__(self):
        return '<'+self.__class__.__name__+' '+repr(self.name)+'>'

class SegmentGroup(ResidueGroup):
    """A group of segments.

    Pythonic access to segments:
      - Using a segid as attribute returns the segment. Because
        of python language rule, any segid starting with a non-letter 
        character is prefixed with 's', thus '4AKE' --> 's4AKE'.
        Example:
          >>> from MDAnalysis.tests.datafiles import PSF,DCD
          >>> u = Universe(PSF,DCD)
          >>> print(u.atoms.segments.s4AKE)  # segment 4AKE
          <AtomGroup with 3314 atoms>
       - Indexing the group returns the appropriate segment.

    :Data: :attr:`SegmentGroup._segments`

    """
    def __init__(self, segments):
        """Initialize the SegmentGroup with a list of :class:`Segment` instances."""
        self._segments = segments
        residues = []
        for s in segments:
            residues.extend(s.residues)
        super(SegmentGroup, self).__init__(residues)
    def __iter__(self):
        return iter(self._segments)
    def __len__(self):
        return len(self._segments)
    def __getitem__(self, item):
        if (type(item) == int) or (type(item) == slice):
            return self._segments[item]
        else: raise NotImplementedError()
    def __getattr__(self, attr):
        if attr.startswith('s') and attr[1].isdigit():
            attr = attr[1:]  # sNxxx only used for python, the name is stored without s-prefix
        seglist = [segment for segment in self._segments if segment.name == attr]
        if len(seglist) == 0:
            return super(SegmentGroup, self).__getattr__(attr)
        if len(seglist) > 1:
            warnings.warn("SegmentGroup: Multiple segments with the same name %r; only the "
                          "FIRST one is returned." % attr, category=SelectionWarning)
        return seglist[0]
    def __repr__(self):
        return '<'+self.__class__.__name__+' '+repr(self._segments)+'>'


class Universe(object):
    """The MDAnalysis Universe contains all the information describing the system.

    Built from a psf, pdb or gro file.

    Attributes: 
       - :attr:`Universe.trajectory`: currently loaded trajectory reader; 
         :attr:`Universe.trajectory.ts` is the current time step
       - :attr:`Universe.dimensions`: current system dimensions (simulation unit cell, if 
         set in the trajectory)
       - bonds, angles, dihedrals, impropers (low level access through :attr:`Universe._psf`)

    Methods:
       - m = Universe(psffile, trajectory)          # read system from file(s)
       - m = Universe(pdbfile)                      # read atoms and coordinates from PDB or GRO
       - m = Universe(psffile, [traj1, traj2, ...]) # read from a list of trajectories
       - m.load_new(trajectory)                     # read from a new trajectory file
       - m.selectAtoms(...)                         # select atoms using similar selection string as CHARMM

    .. Note:: Only a single-frame PDB file is supported; use DCDs or XTC/TRR for
              trajectories. If a PDB is used instead of a PSF then
              charges are not correct, masses are guessed, and bonds are not available.
    """
    def __init__(self, topologyfile, coordinatefile=None, **kwargs):
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
          *permissive*
             Set to ``True`` in order to ignore most errors (currently only relevant
             for PDB files) [``False``]

        This routine tries to do the right thing: 
          1. If a pdb/gro file is provided instead of a psf and no *coordinatefile*
             then the coordinates are taken from the first file. Thus you can load 
             a functional universe with ::

                u = Universe('1ake.pdb')
   
          2. If only a psf file is provided one will have to load coordinates
             manually using :meth:`Universe.load_new_dcd` or
             :meth:`Universe.load_new_pdb`.
        """
        import MDAnalysis.topology.core

        # managed attribute holding TRJReader (the Universe.trajectory
        # attribute is also aliased as Universe.<EXT> where <EXT> is the
        # trajectory format type (i.e. the extension))
        self.__trajectory = None

        # build the topology (or at least a list of atoms)
        try:
            parser = MDAnalysis.topology.core.get_parser_for(topologyfile, permissive=kwargs.get('permissive',False))
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
        # convenience access to residues and segments (these are managed attributes
        # (properties) and are built on the fly or read from a cache) -- does this 
        # create memory problems?
        self.segments = self.atoms.segments
        self.residues = self.atoms.residues

        #MDAnalysis.topology.core.build_bondlists(self.atoms, self._bonds)
        # Let atoms access the universe
        for a in self.atoms:
            a.universe = self

        # Load coordinates
        if coordinatefile is None and \
                MDAnalysis.topology.core.guess_format(topologyfile) in \
                MDAnalysis.coordinates._topology_coordinates_readers:
            coordinatefile = topologyfile         # hack for pdb/gro - only
        self.load_new(coordinatefile, **kwargs)

    def load_new(self, filename, **kwargs):
        """Load coordinates from *filename*, using the suffix to detect file format.

        :Arguments: 
             *filename* 
                 the coordinate file (single frame or trajectory) OR a list of
                 filenames, which are read one after another
             *permissive*
                 If set to ``True``, ignore most errors (at the moment, this is only
                 relevant for PDB) [``False``]
             *kwargs*
                 other kwargs are passed to the trajectory reader (only for advanced use)

        :Returns: (filename, trajectory_format) or ``None`` if *filename* == ``None``
        :Raises: :exc:`TypeError` if trajectory format can not be
                  determined or no appropriate trajectory reader found
        """
        if filename is None:
            return

        import MDAnalysis.coordinates.core
        from itertools import izip

        permissive = kwargs.pop('permissive', False)

        try:
            TRJReader = MDAnalysis.coordinates.core.get_reader_for(filename, permissive=permissive)
        except TypeError, err:
            raise TypeError("Universe.load_new() cannot find an appropriate coordinate reader "
                            "for file %r.\n%r" % (filename, err))
        self.trajectory = TRJReader(filename, **kwargs)    # unified trajectory API
        if self.trajectory.numatoms != self.atoms.numberOfAtoms():
            raise ValueError("The topology and %s trajectory files don't have the same number of atoms!" % self.trajectory.format)
        # hack for PDB
        if self.trajectory.format == "PDB":
            # add B-factor to atoms
            for a, bfactor in izip(self.atoms, self.trajectory.get_bfactors()):
                a.bfactor = bfactor  ## does this work with mmLib??
        return filename, self.trajectory.format
            
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

                point x y z distance

              Make sure coordinate is separated by spaces, i.e. "point 5.0
              5.0 5.0 3.5" selects all atoms within 3.5 Angstroms of the
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
        
