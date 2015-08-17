import numpy
from Bio.KDTree import KDTree

from MDAnalysis.core.AtomGroup import AtomGroup


class AtomNeighborSearch():

    def __init__(self, atom_group, bucket_size=10):
        self.atom_group = atom_group
        if not hasattr(atom_group, 'coordinates'):
            raise TypeError('atom_group must have a coordinates() method'
                            '(eq a AtomGroup from a selection)')
        self.kdtree = KDTree(dim=3, bucket_size=bucket_size)
        self.kdtree.set_coords(atom_group.coordinates())

    def search(self, atoms, radius, level='A'):
        indices = []
        for atom in atoms.coordinates():
            self.kdtree.search(atom, radius)
            indices.append(self.kdtree.get_indices())
        unique_idx = numpy.unique([i for l in indices for i in l])
        return self._index2level(unique_idx, level)

    def _index2level(self, indices, level):
        n_atom_list = [self.atom_group[i] for i in indices]
        if level == 'A':
            if len(n_atom_list) == 0:
                return []
            else:
                return AtomGroup(n_atom_list)
        elif level == 'R':
            return list(set([a.residue for a in n_atom_list]))
        elif level == 'S':
            return list(set([a.segment for a in n_atom_list]))
        else:
            raise NotImplementedError('{}: level not implemented'.format(level))
