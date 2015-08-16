import numpy
from sklearn.neighbors import KDTree

from MDAnalysis.core.AtomGroup import AtomGroup


class AtomNeighborSearch():

    def __init__(self, atom_group, leaf_size=10):
        self.atom_group = atom_group
        if not hasattr(atom_group, 'coordinates'):
            raise TypeError('atom_group must have a coordinates() method'
                            '(eq a AtomGroup from a selection)')
        self.kdtree = KDTree(atom_group.coordinates(), leaf_size=leaf_size)

    def search(self, atoms, radius, level='A'):
        indices = self.kdtree.query_radius(atoms.coordinates(), radius)
        unique_idx = numpy.unique([i for list in indices for i in list])
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
