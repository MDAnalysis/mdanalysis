# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
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
CRD topology parser
===================

Read a list of atoms from a CHARMM CARD coordinate file (CRD)
to build a basic topology.


Classes
-------

.. autoclass:: CRDParser
   :members:
   :inherited-members:

"""
from __future__ import absolute_import

import numpy as np

from ..lib.util import openany, FORTRANReader
from .base import TopologyReader, change_squash
from ..core.topology import Topology
from ..core.topologyattrs import (
    Atomids,
    Atomnames,
    Tempfactors,
    Resids,
    Resnames,
    Resnums,
    Segids,
)


class CRDParser(TopologyReader):
    """Parse a CHARMM CARD coordinate file for topology information.

    Creates the following Attributes:
     - Atomids
     - Atomnames
     - Tempfactors
     - Resids
     - Resnames
     - Resnums
     - Segids
    """
    format = 'CRD'

    def parse(self):
        """Create the Topology object

        Returns
        -------
        MDAnalysis Topology object

        Todo
        ----
        Could use the resnum and temp factor better
        """
        extformat = FORTRANReader('2I10,2X,A8,2X,A8,3F20.10,2X,A8,2X,A8,F20.10')
        stdformat = FORTRANReader('2I5,1X,A4,1X,A4,3F10.5,1X,A4,1X,A4,F10.5')

        atomids = []
        atomnames = []
        tempfactors = []
        resids = []
        resnames = []
        resnums = []
        segids = []

        with openany(self.filename) as crd:
            for linenum, line in enumerate(crd):
                # reading header
                if line.split()[0] == '*':
                    continue
                elif line.split()[-1] == 'EXT' and int(line.split()[0]):
                    r = extformat
                    continue
                elif line.split()[0] == line.split()[-1] and line.split()[0] != '*':
                    r = stdformat
                    continue
                # anything else should be an atom
                try:
                    (serial, TotRes, resName, name,
                     x, y, z, chainID, resSeq, tempFactor) = r.read(line)
                except:
                    raise ValueError("Check CRD format at line {0}: {1}"
                                     "".format(linenum, line.rstrip()))

                atomids.append(serial)
                atomnames.append(name)
                tempfactors.append(tempFactor)
                resids.append(TotRes)
                resnames.append(resName)
                resnums.append(resSeq)
                segids.append(chainID)

        # Convert to np arrays
        atomids = np.array(atomids, dtype=np.int32)
        atomnames = np.array(atomnames, dtype=object)
        tempfactors = np.array(tempfactors, dtype=np.float32)
        resids = np.array(resids, dtype=np.int32)
        resnames = np.array(resnames, dtype=object)
        resnums = np.array(resnums, dtype=np.int32)
        segids = np.array(segids, dtype=object)

        atom_residx, (res_resids, res_resnames, res_resnums, res_segids) = change_squash(
            (resids, resnames), (resids, resnames, resnums, segids))
        res_segidx, (seg_segids,) = change_squash(
            (res_segids,), (res_segids,))

        top = Topology(len(atomids), len(res_resids), len(seg_segids),
                       attrs=[
                           Atomids(atomids),
                           Atomnames(atomnames),
                           Tempfactors(tempfactors),
                           Resids(res_resids),
                           Resnames(res_resnames),
                           Resnums(res_resnums),
                           Segids(seg_segids),
                       ],
                       atom_resindex=atom_residx,
                       residue_segindex=res_segidx)

        return top
