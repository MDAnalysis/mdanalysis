"""

optimisation:
-------------
could store 2 tables for each table, so that access is always fast
would mean that updating the tables needs to update two tables.
"""
from scipy import sparse
import numpy as np


class TransTable(object):
    """Who's your daddy?

    Atom
     ^v     AtomToResidue sparse matrix (nres x natoms)
    Residue
     ^v     ResidueToSegment sparse matrix (nseg x nres)
    Segment

    Accepts only Xidx!
    Returns Xidx!
    Ie. row indices for other tables
    Might differ from Ids!
    """
    def __init__(self,
                 n_atoms, n_residues, n_segments,  # Size of tables
                 Ridx=None, Sidx=None,  # Contents of tables
    ):
        self.n_atoms = n_atoms
        self.n_residues = n_residues
        self.n_segments = n_segments

        RA = self.RA = sparse.csr_matrix((n_residues, n_atoms),
                                         dtype=np.bool)
        SR = self.SR = sparse.csr_matrix((n_segments, n_residues),
                                         dtype=np.bool)

        if not Ridx is None:
            # fill in arrays here
            # could optimise by creating different type of SM then converting
            # we want the CSR? format for quick access across rows
            # ie finding all children of a certain parent
            self.AR = Ridx
            for ai, ri in enumerate(Ridx):
                RA[ri, ai] = True   
        if not Sidx is None:
            self.RS = Sidx
            for ri, si in enumerate(Sidx):
                SR[si, ri] = True

    def a2r(self, Aidx):
        """The Residue index for this Atom index"""
        return self.AR[Aidx]

    def r2a_1d(self, Ridx):
        """The Atom indices for this Ridx"""
        # Sort here so atom indices are sorted within each Residue
        return self.RA[Ridx].sorted_indices().indices

    def r2a_2d(self, Ridx):
        return (row.sorted_indices().indices for row in self.RA[Ridx])

    def r2s(self, Ridx):
        """The Sidx for this Ridx"""
        return self.RS[Ridx]

    def s2r_1d(self, Sidx):
        """The Ridxs for this Sidx"""
        return self.SR[Sidx].indices

    def s2r_2d(self, Sidx):
        return (row.sorted_indices().indices for row in self.SR[Sidx])

    # Compound moves, does 2 translations
    def a2s(self, Aidx):
        """The Sidx for this Aidx"""
        Ridx = self.a2r(Aidx)
        return self.r2s(Ridx)

    def s2a_1d(self, Sidx):
        """The Aidxs for this Sidx"""
        Ridxs = self.s2r_2d(Sidx)
        return np.concatenate([self.r2a_1d(Ridx) for Ridx in Ridxs])

    def s2a_2d(self, Sidx):
        Ridxs = self.s2r_2d(Sidx)
        return (self.r2a_1d(Ridx) for Ridx in Ridxs)

    # Move between different groups.
    # In general, delete old address, add new address
    def move_atom(self, Aidx, Ridx):
        """Move Aidx to be in Ridx"""
        pass

    def move_residue(self, Ridx, Sidx):
        """Move Ridx to be in Sidx"""
        pass

    def resize(self):
        # Shit just got real
        pass
