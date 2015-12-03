# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""
Topology object --- :mod:`MDAnalysis.core.topology'
===================================================================

"""


class Topology(object):
    """In-memory, array-based topology database.

    The topology model of MDanalysis features atoms, which can each be a member
    of one or zero residues. Each residue, in turn, can be a member of one
    or zero segments. The details of maintaining this heirarchy, and mappings
    of atoms to residues, residues to segments, and vice-versa, are handled
    internally by this object.

    Parameters
    ----------
    n_atoms, n_residues, n_segments : int
        number of atoms, residues, segments in topology
    topologyattrs : TopologyAttr objects
        components of the topology to be included

    """

    def __init__(self, n_atoms, n_residues, n_segments, *topologyattrs):

        # attach the TopologyAttrs
        for topologyattr in topologyattrs:
            self.add_TopologyAttr(topologyattr)


    def add_TopologyAttr(topologyattr):
        """Add a new TopologyAttr to the Topology.

        Parameters
        ----------
        topologyattr : TopologyAttr

        """
        topologyattr.topology = self
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

