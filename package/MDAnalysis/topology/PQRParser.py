# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
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
PQR topology parser
===================

Read atoms with charges from a PQR_ file (as written by PDB2PQR_). No
connectivity is deduced.

.. SeeAlso:: The file format is described in :mod:`MDAnalysis.coordinates.PQR`.

.. _PQR:     http://www.poissonboltzmann.org/file-formats/biomolecular-structurw/pqr
.. _APBS:    http://www.poissonboltzmann.org/apbs
.. _PDB2PQR: http://www.poissonboltzmann.org/pdb2pqr
.. _PDB:     http://www.rcsb.org/pdb/info.html#File_Formats_and_Standards

Classes
-------

.. autoclass:: PQRParser
   :members:
   :inherited-members:

"""
from __future__ import absolute_import

import numpy as np

from . import guessers
from ..lib.util import openany
from ..core.topologyattrs import (
    Atomids,
    Atomnames,
    Charges,
    Elements,
    Masses,
    Radii,
    Resids,
    Resnums,
    Resnames,
    Segids,
)
from ..core.topology import Topology
from .base import TopologyReader, squash_by


class PQRParser(TopologyReader):
    """Parse atom information from PQR file *filename*.

    Creates a MDAnalysis Topology with the following attributes
     - Atomids
     - Atomnames
     - Charges
     - Radii
     - Resids
     - Resnames
     - Segids

    Guesses the following:
     - elements
     - masses

    .. versionchanged:: 0.9.0
       Read chainID from a PQR file and use it as segid (before we always used
       'SYSTEM' as the new segid).
    """
    format = 'PQR'

    def parse(self):
        """Parse atom information from PQR file *filename*.

        Returns
        -------
        A MDAnalysis Topology object
        """
        serials = []
        names = []
        resnames = []
        chainIDs = []
        resids = []
        charges = []
        radii = []

        with openany(self.filename, 'r') as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    fields = line.split()
                    try:
                        (recordName, serial, name, resName,
                         chainID, resSeq, x, y, z, charge,
                         radius) = fields
                    except ValueError:
                        # files without the chainID
                        (recordName, serial, name, resName,
                         resSeq, x, y, z, charge, radius) = fields
                        chainID = "SYSTEM"
                    serials.append(serial)
                    names.append(name)
                    resnames.append(resName)
                    resids.append(resSeq)
                    charges.append(charge)
                    radii.append(radius)
                    chainIDs.append(chainID)

        n_atoms = len(serials)

        elements = guessers.guess_types(names)
        masses = guessers.guess_masses(elements)

        attrs = []
        attrs.append(Atomids(np.array(serials, dtype=np.int32)))
        attrs.append(Atomnames(np.array(names, dtype=object)))
        attrs.append(Charges(np.array(charges, dtype=np.float32)))
        attrs.append(Elements(elements, guessed=True))
        attrs.append(Masses(masses, guessed=True))
        attrs.append(Radii(np.array(radii, dtype=np.float32)))

        resids = np.array(resids, dtype=np.int32)
        resnames = np.array(resnames, dtype=object)
        chainIDs = np.array(chainIDs, dtype=object)

        residx, resids, (resnames, chainIDs) = squash_by(
            resids, resnames, chainIDs)

        n_residues = len(resids)
        attrs.append(Resids(resids))
        attrs.append(Resnums(resids.copy()))
        attrs.append(Resnames(resnames))

        segidx, chainIDs = squash_by(chainIDs)[:2]

        n_segments = len(chainIDs)
        attrs.append(Segids(chainIDs))

        top = Topology(n_atoms, n_residues, n_segments,
                       attrs=attrs,
                       atom_resindex=residx,
                       residue_segindex=segidx)

        return top
