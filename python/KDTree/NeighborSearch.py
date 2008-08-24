# $Id$
# Neighbour searching in MDAnalysis. 
# Based on Biopython's Bio/PDB/NeighborSearch.py
# [Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)]
# under the terms of the Biopython license.

"""Fast atom neighbor lookup using a KD tree (implemented in C++)."""

import numpy
from KDTree import KDTree
from MDAnalysis.core.AtomGroup import AtomGroup
import sets

entity_levels = ('A','R','S')

class NeighborSearch:
    """
    This class can be used for two related purposes:

    1. To find all atoms/residues/segments within radius 
    of a given query position. 

    2. To find all atoms/residues/segments that are within 
    a fixed radius of each other.

    NeighborSearch makes use of the KDTree C++ module, so it's fast.
    """
    def __init__(self, atom_list, bucket_size=10):
        """
        o atom_list - list of atoms. This list is used in the queries.
        It can contain atoms from different structures.
        o bucket_size - bucket size of KD tree. You can play around 
        with this to optimize speed if you feel like it.
        """
        self.atom_list=atom_list
        # to Nx3 array of type float
        try:
            self.coords=numpy.asarray(atom_list.coordinates(),dtype=numpy.float32)
        except AttributeError:
            raise TypeError('atom_list must have a coordinates() method '
                            '(eg a AtomGroup resulting from a Selection)')
        assert(bucket_size>1)
        assert(self.coords.shape[1]==3)
        assert(self.coords.dtype == numpy.float32)
        self.kdt=KDTree(3, bucket_size)
        self.kdt.set_coords(self.coords)
    
    # Public

    def search(self, center, radius, level="A"):
        """Neighbor search.

        Return all atoms/residues/segments
        that have at least one atom within radius of center.
        What entitity level is returned (e.g. atoms or residues)
        is determined by level (A=atoms, R=residues, S=Segments).

        o center - numpy array 
        o radius - float
        o level - char (A, R, S)
        """
        if not level in entity_levels:
            raise ValueError("%s: Unknown level" % level)
        self.kdt.search(center, radius)
        indices=self.kdt.get_indices()
        n_atom_list=[]
        atom_list=self.atom_list
        for i in indices:
            a=atom_list[i]
            n_atom_list.append(a)
        if level=="A":
            return AtomGroup(n_atom_list)
        elif level=="R":
            residues = sets.Set([a.residue for a in n_atom_list])
            return list(residues)
        elif level=="S":
            segments = sets.Set([a.segment for a in n_atom_list])
            return list(residues)
        else:
            raise NotImplementedError("level=%s not implemented" % level)
            
    def search_all(self, radius, level="A"):
        """All neighbor search.

        Search all entities that have atoms pairs within
        radius. 

        o radius - float
        o level - char (A, R, S)
        """
        if not level in entity_levels:
            raise ValueError("%s: Unknown level" % level)
        self.kdt.all_search(radius)
        indices=self.kdt.all_get_indices()
        atom_list=self.atom_list
        atom_pair_list=[]
        for i1, i2 in indices:
            a1=atom_list[i1]
            a2=atom_list[i2]
            atom_pair_list.append(AtomGroup([a1, a2]))
        if level=="A":            
            return atom_pair_list  # return atoms as list of AtomGroup pairs
        elif level == "R":
            return self._get_unique_pairs('residue',atom_pair_list)
        elif level == "S":
            return self._get_unique_pairs('segment',atom_pair_list)
        else:
            raise NotImplementedError("level=%s not implemented" % level)

    # Private

    def _get_unique_pairs(self,entity,atom_pair_list):        
        # use sets to remove duplicates
        unique_pairs = sets.Set([
                sets.ImmutableSet((a1.__getattribute__(entity),
                                   a2.__getattribute__(entity)))
                for a1,a2 in atom_pair_list 
                if a1.__getattribute__(entity) != a2.__getattribute__(entity)])
        return [tuple(s) for s in unique_pairs]  # return as list of 2-tuples

if __name__=="__main__":
    import numpy

    class Atom:
        def __init__(self):
            self.coord=100*numpy.random.random(3)
        def coordinates(self):
            return self.coord

    class AtomGroup(list):
        def coordinates(self):
            return numpy.array( [atom.coordinates() for atom in self] )

    for i in range(0, 20):
        al=AtomGroup([Atom() for i in range(0, 1000)])

        ns=NeighborSearch(al)

        print "Found ", len(ns.search_all(5.0))


            

                
        
