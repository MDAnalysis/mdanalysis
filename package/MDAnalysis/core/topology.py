
class Topology(object):
    """In-memory, array-based topology database.

    Parameters
    ----------
    topologyattrs : TopologyAttr objects
        components of the topology to be included

    """

    def __init__(self, *topologyattrs):

        # attach the TopologyAttrs
        for topologyattr in topologyattrs:
            self.__setattr__(topologyattr.attrname, topologyattr)


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
        pass

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
        pass

    def r2a(self, rix):
        """Get atom indices collectively represented by given residue indices.

        Parameters
        ----------
        rix : array
            residue indices

        Returns
        -------
        aix : array
            sorted indices of atoms present in residues, collectively

        """
        pass

    def r2ra(self, rix):
        """Get atom indices represented by each residue index.

        Parameters
        ----------
        rix : array
            residue indices

        Returns
        -------
        raix : sparse matrix 
            each row corresponds to a residue index, in order given in `rix`,
            each column corresponds to an atom, with either a 1 or 0 as the
            value indicating membership or not, respectively, in that residue

        """
        pass

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
        pass

    def s2r(self, six):
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
        pass

    def s2sr(self, six):
        """Get residue indices represented by each segment index.

        Parameters
        ----------
        six : array
            residue indices

        Returns
        -------
        srix : sparse matrix 
            each row corresponds to a segment index, in order given in `six`,
            each column corresponds to a residue, with either a 1 or 0 as the
            value indicating membership or not, respectively, in that segment

        """
        pass

    def s2a(self, six):
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
        pass

    def s2sa(self, six):
        """Get atom indices represented by each segment index.

        Parameters
        ----------
        six : array
            residue indices

        Returns
        -------
        saix : sparse matrix 
            each row corresponds to a segment index, in order given in `six`,
            each column corresponds to a atom, with either a 1 or 0 as the
            value indicating membership or not, respectively, in that segment

        """
        pass

