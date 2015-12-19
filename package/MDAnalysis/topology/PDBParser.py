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
PDB topology parser
===================

Use a PDB file to build a minimum internal structure representation.

.. SeeAlso:: :mod:`MDAnalysis.coordinates.PDB` and :mod:`Bio.PDB`

.. SeeAlso:: :mod:`MDAnalysis.topology.PrimitivePDBParser` (which
             *can* guess conectivity but does not support all subleties of the full
             PDB format)

Classes
-------

.. autoclass:: PDBParser
   :members:
   :inherited-members:

"""
from __future__ import absolute_import

import numpy as np

try:
    # BioPython is overkill but potentially extensible (altLoc etc)
    import Bio.PDB
except ImportError:
    raise ImportError("Bio.PDB from biopython not found."
                      "Required for PDB topology parser.")

from .base import TopologyReader, squash_by
from ..coordinates.pdb.extensions import get_structure
from ..core.topology import Topology
from ..core.topologyattrs import (
    Atomnames,
    Bfactors,
    Occupancies,
    Resids,
    Resnames,
    Segids,
)


class PDBParser(TopologyReader):
    """Read minimum topology information from a PDB file.

    Creates the following Attributes:
     - names
     - bfactors
     - occupancies
     - resids
     - resnames
     - segids
    """
    def parse(self):
        """Parse atom information from PDB file *pdbfile*.

        This functions uses the :class:`Bio.PDB.PDBParser` as used by
        :func:`MDAnalysis.coordinates.pdb.extensions.get_structure`.
        """
        # use Sloppy PDB parser to cope with big PDBs!
        pdb = get_structure(self.filename, "0UNK")

        names = []
        resids = []
        resnames = []
        segids = []
        bfactors = []
        occupancies = []
        # translate Bio.PDB atom objects to MDAnalysis Atom.
        for atom in pdb.get_atoms():
            residue = atom.parent
            chain_id = residue.parent.id
            atomname = atom.name
            resname = residue.resname
            resid = int(residue.id[1])
            # no empty segids (or Universe throws IndexError)
            segid = residue.get_segid().strip() or chain_id or "SYSTEM"
            bfactor = atom.bfactor
            occupancy = atom.occupancy

            names.append(atomname)
            resids.append(resid)
            resnames.append(resname)
            segids.append(segid)
            bfactors.append(bfactor)
            occupancies.append(occupancy)

        attrs = []
        n_atoms = len(names)
        attrs.append(Atomnames(np.array(names, dtype=object)))
        attrs.append(Bfactors(np.array(bfactors, dtype=np.float32)))
        attrs.append(Occupancies(np.array(occupancies, dtype=np.float32)))

        resids = np.array(resids, dtype=np.int32)
        resnames = np.array(resnames, dtype=object)
        segids = np.array(segids, dtype=object)

        residx, resids, (resnames, segids) = squash_by(
            resids, resnames, segids)
        n_residues = len(resids)
        attrs.append(Resids(resids))
        attrs.append(Resnames(resnames))

        segidx, segids = squash_by(segids)[:2]
        n_segments = len(segids)
        attrs.append(Segids(segids))

        top = Topology(n_atoms, n_residues, n_segments,
                       attrs=attrs,
                       atom_resindex=residx,
                       residue_segindex=segidx)

        return top
