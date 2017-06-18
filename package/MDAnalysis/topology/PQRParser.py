# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""
PQR topology parser
===================

Read atoms with charges from a PQR_ file (as written by PDB2PQR_). No
connectivity is deduced.

Note
----
The file format is described in :mod:`MDAnalysis.coordinates.PQR`.


Classes
-------

.. autoclass:: PQRParser
   :members:
   :inherited-members:


.. _PQR:     http://www.poissonboltzmann.org/file-formats/biomolecular-structurw/pqr
.. _APBS:    http://www.poissonboltzmann.org/apbs
.. _PDB2PQR: http://www.poissonboltzmann.org/pdb2pqr
.. _PDB:     http://www.rcsb.org/pdb/info.html#File_Formats_and_Standards

"""
from __future__ import absolute_import

import numpy as np

from . import guessers
from ..lib.util import openany
from ..core.topologyattrs import (
    Atomids,
    Atomnames,
    Atomtypes,
    Charges,
    ICodes,
    Masses,
    Radii,
    Resids,
    Resnums,
    Resnames,
    Segids,
)
from ..core.topology import Topology
from .base import TopologyReaderBase, squash_by, change_squash


class PQRParser(TopologyReaderBase):
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
     - atomtypes
     - masses

    .. versionchanged:: 0.9.0
       Read chainID from a PQR file and use it as segid (before we always used
       'SYSTEM' as the new segid).
    .. versionchanged:: 0.16.1
       Now reads insertion codes and splits into new residues around these
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
        icodes = []
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
                    try:
                        resid = int(resSeq)
                    except ValueError:
                        # has icode present
                        resid = int(resSeq[:-1])
                        icode = resSeq[-1]
                    else:
                        icode = ''

                    serials.append(serial)
                    names.append(name)
                    resnames.append(resName)
                    resids.append(resid)
                    icodes.append(icode)
                    charges.append(charge)
                    radii.append(radius)
                    chainIDs.append(chainID)

        n_atoms = len(serials)

        atomtypes = guessers.guess_types(names)
        masses = guessers.guess_masses(atomtypes)

        attrs = []
        attrs.append(Atomids(np.array(serials, dtype=np.int32)))
        attrs.append(Atomnames(np.array(names, dtype=object)))
        attrs.append(Charges(np.array(charges, dtype=np.float32)))
        attrs.append(Atomtypes(atomtypes, guessed=True))
        attrs.append(Masses(masses, guessed=True))
        attrs.append(Radii(np.array(radii, dtype=np.float32)))

        resids = np.array(resids, dtype=np.int32)
        icodes = np.array(icodes, dtype=object)
        resnames = np.array(resnames, dtype=object)
        chainIDs = np.array(chainIDs, dtype=object)

        residx, (resids, resnames, icodes, chainIDs) = change_squash(
            (resids, resnames, icodes, chainIDs),
            (resids, resnames, icodes, chainIDs))

        n_residues = len(resids)
        attrs.append(Resids(resids))
        attrs.append(Resnums(resids.copy()))
        attrs.append(Resnames(resnames))
        attrs.append(ICodes(icodes))

        segidx, chainIDs = squash_by(chainIDs)[:2]

        n_segments = len(chainIDs)
        attrs.append(Segids(chainIDs))

        top = Topology(n_atoms, n_residues, n_segments,
                       attrs=attrs,
                       atom_resindex=residx,
                       residue_segindex=segidx)

        return top
