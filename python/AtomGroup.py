"""AtomGroup Hierarchy

Class Hierarchy:
    AtomGroup ->
"""

from Scientific.Geometry import Vector, Tensor

import Numeric

class Atom(object):
    """A single atom definition

Data: number, segid, resid, resname, name, type

Methods:
    a = Atom()
"""
    __slots__ = ("number", "id", "name", "type", "resname", "resid", "segid", "mass", "charge", "_coord", "residue", "segment", "bonds")

    def __init__(self, number, name, type, resname, resid, segid, mass, charge):
        self.number = number
        self.name = name
        self.type = type
        self.resname = resname
        self.resid = resid
        self.segid = segid
        self.mass = mass
        self.charge = charge
    def __repr__(self):
        return "< Atom " + repr(self.number+1) + ": name " + repr(self.name) +" of type " + \
               repr(self.type) + " of resid " + repr(self.resname) + ", " +repr(self.resid) + " and " +repr(self.segid)+'>'
    def __getattr__(self, key):
        if (key == "pos"):
            if hasattr(self, "_coord"):
                return self._coord[self.number] # PDB numbering starts at 0
            else: raise Exception("Atom "+repr(self.number)+" is not assigned to a System")
        else: return super(Atom, self).__getattribute__(key)
    def __cmp__(self, other):
        return cmp(self.number, other.number)
    def __eq__(self, other):
        return self.number == other.number
    def __hash__(self):
        return hash(self.number)

class NameError(Exception):
    pass

class AtomGroup(object):
    """A group of atoms

Currently contains a list of atoms from the main system that correspond to this group.

Data: _atoms - a list of references to the corresponding atoms in System._atoms
      AtomGroups are immutable

Methods:
    ag = AtomGroup()
    ag.indices() - return indices into main atom array
    ag.masses() - array of masses
    ag.totalMass() - total mass
    ag.charges() - array of charges
    ag.totalCharge() - total charge
    ag.centerOfGeom() - center of geometry
    ag.centerOfMass() - center of mass
    ag.radiusOfGyration() - radius of gyration
    ag.principleAxis() - returns the principle axis of rotation
    c = ag.coordinates() - return array of coordinates
"""
    def __init__(self, atoms):
        if len(atoms) < 1: raise Exception("No atoms defined for AtomGroup")
        self._atoms = atoms
        # If the number of atoms is very large, create a dictionary cache for lookup
        if len(atoms) > 1000:
            self.__atom_cache = dict([(x,None) for x in atoms])
    def __len__(self):
        return len(self._atoms)
    def __getitem__(self, item):
        if (type(item) is int) or (type(item) is slice):
            return self._atoms[item]
        else: return super(AtomGroup, self).__getitem__(item)
    def __getattr__(self, name):
        # There can be more than one atom with the same name
        atomlist = [atom for atom in self._atoms if name == atom.name]
        if len(atomlist) == 0: raise NameError("No atoms with name "+name)
        elif len(atomlist) == 1: return atomlist[0]
        else: return AtomGroup(atomlist)
    def __iter__(self):
        return iter(self._atoms)
    def __contains__(self, other):
        if hasattr(self, "__atom_cache"):
            return other in self.__atom_cache
        else: return other in self._atoms
    def __add__(self, other):
        if other.__class__ is not self.__class__:
            raise TypeError('Can only concatenate AtomGroup (not "'+repr(self.__class__.__name__)+'")to AtomGroup')
        return AtomGroup(self._atoms+other._atoms)
    def __repr__(self):
        return '<'+self.__class__.__name__+' with '+repr(len(self._atoms))+' atoms>'
    def indices(self):
        return Numeric.array([atom.number for atom in self._atoms])
    def masses(self):
        if not hasattr(self, "_masses"):
            self._masses = Numeric.array([atom.mass for atom in self._atoms])
        return self._masses
    def totalMass(self):
        return Numeric.sum(self.masses())
    def charges(self):
        return Numeric.array([atom.charge for atom in self._atoms])
    def totalCharge(self):
        return Numeric.sum(self.charges())
    def centerOfGeom(self):
        t = reduce(lambda total,atom: total+atom.pos, self._atoms, Numeric.array([0.,0.,0.]))
        return t/len(self._atoms)
    def centerOfMass(self):
        t = Numeric.sum(self.coordinates()*self.masses()[:,Numeric.NewAxis])
        return t/self.totalMass()
    def radiusOfGyration(self):
        masses = self.masses()
        recenteredpos = self.coordinates() - self.centerOfMass()
        rog_sq = Numeric.sum(masses*Numeric.add.reduce(Numeric.power(recenteredpos, 2), axis=1))/self.totalMass()
        return Numeric.sqrt(rog_sq)
    def principleAxis(self):
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
        from LinearAlgebra import eigenvectors
        eigenval, eigenvec = eigenvectors(Numeric.array([[Ixx, Ixy, Ixz],[Iyx, Iyy, Iyz],[Izx, Izy, Izz]]))
        # Find minimum
        return min(zip(eigenval, eigenvec))[1]
    def coordinates(self, ts=None):
        if ts == None:
            return Numeric.array([atom.pos for atom in self._atoms])
        else:
            return Numeric.array([ts[i] for i in self.indices()])

class Residue(AtomGroup):
    """A group of atoms corresponding to a residue

Data: type, name

Methods:
    r = Residue()
    r['name'] or r[id]
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
            Residue.__cache[name] = dict([(a.name, i) for i, a in enumerate(self._atoms)])
    def __getitem__(self, item):
        if (type(item) is int) or (type(item) is slice):
            return self._atoms[item]
        else: return self.__getattr__(item)
    def __getattr__(self, name):
        # There can only be one atom with a certain name
        # Use the cache
        #for atom in self._atoms:
        #    if (name == atom.name): return atom
        try:
            index = Residue.__cache[self.name][name]
            return self._atoms[index]
        except KeyError: raise NameError("No atom in residue "+self.name+" with name "+name)
    def __repr__(self):
        return '<'+self.__class__.__name__+' '+repr(self.name)+', '+repr(self.id)+'>'

class ResidueGroup(AtomGroup):
    """ A group of residues

Data: residues

Methods:
    rg = ResidueGroup()
"""
    def __init__(self, residues):
        self._residues = residues
        atoms = []
        for res in residues:
            atoms.extend(res._atoms)
        super(ResidueGroup, self).__init__(atoms)
    def __iter__(self):
        return iter(self._residues)
    def __len__(self):
        return len(self._residues)
    def __getitem__(self, item):
        if (type(item) is int) or (type(item) is slice):
            return self._residues[item]
        else: raise NotImplementedError()
    def __getattr__(self, attr):
        atomlist = [atom for atom in self._atoms if atom.name == attr]
        return AtomGroup(atomlist)
    def __repr__(self):
        return '<'+self.__class__.__name__+' '+repr(self._residues)+'>'

class Segment(ResidueGroup):
    """ A group of residues corresponding to one segment of the PSF

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

class Universe(AtomGroup):
    """Contains all the information decribing the system
       Built from a psf file

Data: bonds, angles, dihedrals, impropers, donors, acceptors

Methods:
   m = Universe(psffile, dcdfile) - read system from file(s)

See also:
"""
    def __init__(self, psffilename, dcdfilename=None):
        import PSFParser
        struc = PSFParser.parse(psffilename)
        for data in struc.keys():
            setattr(self, data, struc[data])
        self.__dict__.update(PSFParser._buildstructure_(self._atoms))
        PSFParser._buildbondlists_(self._atoms, self._bonds)
        if dcdfilename is not None:
            self.__init_dcd(dcdfilename)
    def __init_dcd(self, dcdfilename):
        # Read in trajectory file
        from DCD import DCDReader
        self._dcd = DCDReader(dcdfilename)
        # Make sure that they both have the same number of atoms
        if (self._dcd.numatoms != len(self._atoms)):
            raise Exception("The psf and dcd files don't have the same number of atoms!")
        self._coord = self._dcd.ts
        # Let atoms access their current positions
        for a in self._atoms:
            a._coord = self._coord
    def dimensions(self):
        return self._coord.dimensions
    def __getattr__(self, name):
        if name == "_coord":
            raise AttributeError("System was not defined with a dcdfile: no coordinates available")
        elif name == "timesteps":
            return len(self._dcd)
        else: return super(System, self).__getattribute__(name)
    def numberOfAtoms(self):
        return len(self._atoms)
    def getTrajectory(self):
        return self._dcd
    def getConfiguration(self):
        return self._coord
    def selectAtoms(self, sel, *othersel):
        # XXX This uses the outdated selection parser
        # You should be able to do everything else using classes in AtomGroup
        import Selection
        atomgrp = Selection.Parser.parse(sel).apply(self)
        if len(othersel) > 0:
            others = []
            for sel in othersel:
                atomgrp = atomgrp + Selection.Parser.parse(sel).apply(self)
        return atomgrp
    def __repr__(self):
        return '<'+self.__class__.__name__+' with '+repr(len(self._atoms))+' atoms>'
    def recenterAround(self, group):
        com = group.centerOfMass()
        self._coord -= com
