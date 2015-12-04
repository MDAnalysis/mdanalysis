"""

optimisation:
-------------
could store 2 tables for each table, so that access is always fast
would mean that updating the tables needs to update two tables.
"""
from scipy import sparse
import numpy as np


class TransTable(object):
    """Membership tables with methods to translate indices across levels.

    Parameters
    ----------
    n_atoms, n_residues, n_segments : int
        number of atoms, residues, segments in topology
    atom_resindex : 1-D array
        resindex for each atom in the topology; the number of unique values in this
        array must be <= `n_residues`, and the array must be length `n_atoms`;
        giving None defaults to placing all atoms in residue 0
    residue_segindex : 1-D array
        segindex for each atom in the topology; the number of unique values in this
        array must be <= `n_segments`, and the array must be length `n_residues`;
        giving None defaults to placing all residues in segment 0
 

    Attributes
    ----------
    n_atoms, n_residues, n_segments : int
        number of atoms, residues, segments in topology
    AR : 1-D array
        resindex for each atom in the topology; allows fast upward translation
        from atoms to residues 
    RA : sparse matrix
        row ``i`` corresponds to the residue with resindex ``i``, with each
        column giving 1 if the atom with atomindex ``j`` is a member or 0 if it
        is not; this matrix has dimension (nres x natoms); allows fast downward
        translation from residues to atoms
    RS : 1-D array
        segindex for each residue in the topology; allows fast upward
        translation from residues to segments 
    SR : sparse matrix
        row ``i`` corresponds to the segment with segindex ``i``, with each
        column giving 1 if the residue with resindex ``j`` is a member or 0 if
        it is not; this matrix has dimension (nseg x nres); allows fast
        downward translation from segments to residues 

    """
    def __init__(self,
                 n_atoms, n_residues, n_segments,  # Size of tables
                 atom_resindex=None, residue_segindex=None,  # Contents of tables
                ):
        self.n_atoms = n_atoms
        self.n_residues = n_residues
        self.n_segments = n_segments

        self.RA = sparse.csr_matrix((n_residues, n_atoms),
                                         dtype=np.bool)
        self.SR = sparse.csr_matrix((n_segments, n_residues),
                                         dtype=np.bool)

        if not atom_resindex is None:
            # fill in arrays here
            # could optimise by creating different type of SM then converting
            # we want the CSR? format for quick access across rows
            # ie finding all children of a certain parent
            self.AR = atom_resindex
            for ai, ri in enumerate(atom_resindex):
                self.RA[ri, ai] = True   
        if not residue_segindex is None:
            self.RS = residue_segindex
            for ri, si in enumerate(residue_segindex):
                self.SR[si, ri] = True

    def a2r(self, aix):
        """Get residue indices for each atom.

        Parameters
        ----------
        aix : array
            atom indices

        Returns
        -------
        rix : array
            residue index for each atom 

        """
        return self.AR[aix]

    def r2a_1d(self, rix):
        """Get atom indices collectively represented by given residue indices.

        Parameters
        ----------
        rix : array
            residue indices

        Returns
        -------
        aix : array
            indices of atoms present in residues, collectively

        """
        return self.RA[rix].sorted_indices().indices

    def r2a_2d(self, rix):
        """Get atom indices represented by each residue index.

        Parameters
        ----------
        rix : array
            residue indices

        Returns
        -------
        raix : tuple
            each element corresponds to a residue index, in order given in
            `rix`, with each element being an array of the atom indices present
            in that residue

        """
        return (row.sorted_indices().indices for row in self.RA[rix])

    def r2s(self, rix):
        """Get segment indices for each residue.

        Parameters
        ----------
        rix : array
            residue indices 

        Returns
        -------
        six : array
            segment index for each residue

        """
        return self.RS[rix]

    def s2r_1d(self, six):
        """Get residue indices collectively represented by given segment indices.

        Parameters
        ----------
        six : array
            segment indices

        Returns
        -------
        rix : array
            sorted indices of residues present in segments, collectively

        """
        return self.SR[six].sorted_indices().indices

    def s2r_2d(self, six):
        """Get residue indices represented by each segment index.

        Parameters
        ----------
        six : array
            residue indices

        Returns
        -------
        srix : sparse matrix 
            each element corresponds to a segment index, in order given in
            `six`, with each element being an array of the residue indices present
            in that segment

        """
        return (row.sorted_indices().indices for row in self.SR[six])

    # Compound moves, does 2 translations
    def a2s(self, aix):
        """Get segment indices for each atom.

        Parameters
        ----------
        aix : array
            atom indices

        Returns
        -------
        rix : array
            segment index for each atom

        """
        rix = self.a2r(aix)
        return self.r2s(rix)

    def s2a_1d(self, six):
        """Get atom indices collectively represented by given segment indices.

        Parameters
        ----------
        six : array
            segment indices

        Returns
        -------
        aix : array
            sorted indices of atoms present in segments, collectively

        """
        rixs = self.s2r_2d(six)
        return np.concatenate([self.r2a_1d(rix) for rix in rixs])

    def s2a_2d(self, six):
        """Get atom indices represented by each segment index.

        Parameters
        ----------
        six : array
            residue indices

        Returns
        -------
        saix : sparse matrix 
            each element corresponds to a segment index, in order given in
            `six`, with each element being an array of the atom indices present
            in that segment

        """
        rixs = self.s2r_2d(six)
        return (self.r2a_1d(rix) for rix in rixs)

    # Move between different groups.
    # In general, delete old address, add new address
    def move_atom(self, aix, rix):
        """Move aix to be in rix"""
        pass

    def move_residue(self, rix, six):
        """Move rix to be in six"""
        pass

    def resize(self):
        pass
