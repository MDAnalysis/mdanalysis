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
GRO topology parser
===================

Read a list of atoms from a GROMOS/Gromacs GRO coordinate file to
build a basic topology.

Atom types, masses and charges are guessed.

.. SeeAlso:: :mod:`MDAnalysis.coordinates.GRO`

Classes
-------

.. autoclass:: GROParser
   :members:
   :inherited-members:

"""
from __future__ import absolute_import

import numpy as np
from six.moves import range

from ..lib.util import openany
from ..core.topologyattrs import (
    Resids,
    Resnames,
    Atomids,
    Atomnames,
)
from ..core.topology import Topology
from .base import TopologyReader, squash_by


class GROParser(TopologyReader):
    format = 'GRO'

    def parse(self):
        """Return the *Topology* object for this file"""
        # Gro has the following columns
        # resid, resname, name, index, (x,y,z)
        with openany(self.filename, 'rt') as inf:
            inf.readline()
            n_atoms = int(inf.readline())

            # Allocate shizznizz
            resids = np.zeros(n_atoms, dtype=np.int32)
            resnames = np.zeros(n_atoms, dtype=object)
            names = np.zeros(n_atoms, dtype=object)
            indices = np.zeros(n_atoms, dtype=np.int32)

            for i in range(n_atoms):
                line = inf.readline()
                try:
                    resids[i] = int(line[:5])
                    resnames[i] = line[5:10].strip()
                    names[i] = line[10:15].strip()
                    indices[i] = int(line[15:20])
                except (ValueError, TypeError):
                    raise IOError(
                        "Couldn't read the following line of the .gro file:\n"
                        "{0}".format(line))

        residx, new_resids, (new_resnames,) = squash_by(resids, resnames)

        # new_resids is len(residues)
        # so resindex 0 has resid new_resids[0]
        atomnames = Atomnames(names)
        atomids = Atomids(indices)
        residueids = Resids(new_resids)
        residuenames = Resnames(new_resnames)

        top = Topology(n_atoms=n_atoms, n_res=len(new_resids), n_seg=1,
                       attrs=[atomnames, atomids, residueids, residuenames],
                       atom_resindex=residx,
                       residue_segindex=None)

        return top
