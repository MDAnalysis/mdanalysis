"""The all singing all dancing new Topology system"""

class Topology(object):

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
            residue indices

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
            residue indices

        """
        pass

    def r2a(self, rix):
        pass

    def r2ra(self, rix):
        pass

    def r2s(self, rix):
        pass

    def s2r(self, six):
        pass

    def s2sr(self, six):
        pass

    def s2a(self, six):
        pass

    def s2sa(self, six):
        pass

