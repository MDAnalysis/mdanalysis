class MoleculeType(object):
    def __init__(self, name, atomkinds, bonds=None, angles=None, 
                 dihe=None, impr=None, donors=None, acceptors=None):
        self.name = name
        self.atomkinds = atomkinds
        self.bonds = bonds
        self.angles = angles
        self.dihe = dihe
        self.impr = impr
        self.donors = donors
        self.acceptors = acceptors

    def number_of_atoms(self):
        return len(self.atomkinds)

    def number_of_residues(self):
        return len(set([a.resid for a in self.atomkinds]))

    # remap_ method returns [blah, blah, ..] or []
    def remap_bonds(self, atom_start_ndx):
        if self.bonds:
            return [[i+atom_start_ndx for i in b] for b in self.bonds]
        else:
            return []

    def remap_angles(self, atom_start_ndx):
        if self.bonds:
            return [[i+atom_start_ndx for i in a] for a in self.angles]
        else:
            return []

    def remap_dihe(self, atom_start_ndx):
        if self.dihe:
            return [[i+atom_start_ndx for i in a] for a in self.dihe]
        else:
            return []

    def remap_impr(self, atom_start_ndx):
        if self.impr:
            return [[i+atom_start_ndx for i in a] for a in self.impr]
        else:
            return []

class AtomKind(object):
    def __init__(self, id, name, type, resid, resname, mass, charge):
        # id is only within the scope of a single molecule, not the whole system
        self.id = id
        self.name = name
        self.type = type
        self.resid = resid
        self.resname = resname
        self.mass = mass
        self.charge = charge

    def __repr__(self):
        return ("< AtomKind: id {0:6d}, name {1:5s}, type {2:2s}, resid {3:6d}, resname {4:3s}, mass {5:8.4f}, charge {6:12.3f} >".format(
                self.id, self.name, self.type, self.resid, self.resname, self.mass, self.charge))

class InteractionType(object):
    def __init__(self, name, long_name, natoms):
        """natoms: number of atoms involved in this type of interaction"""
        self.name = name
        self.long_name = long_name
        self.natoms = natoms
        
    def process(self, atom_ndx):
        while atom_ndx:
            # format for all info: (type, [atom1, atom2, ...])
            # yield atom_ndx.pop(0), [atom_ndx.pop(0) for i in range(self.natoms)]

            # but currently only [atom1, atom2, ...] is interested
            atom_ndx.pop(0)
            yield [atom_ndx.pop(0) for i in range(self.natoms)]
