# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
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


.. _PQR:     https://apbs-pdb2pqr.readthedocs.io/en/latest/formats/pqr.html
.. _APBS:    https://apbs.readthedocs.io/en/latest/
.. _PDB2PQR: https://apbs-pdb2pqr.readthedocs.io/en/latest/pdb2pqr/index.html
.. _PDB:     http://www.wwpdb.org/documentation/file-format

"""
import numpy as np

from ..lib.util import openany
from ..core.topologyattrs import (
    Atomids,
    Atomnames,
    Atomtypes,
    Charges,
    ICodes,
    Radii,
    RecordTypes,
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
     - RecordTypes (ATOM/HETATM)
     - Resids
     - Resnames
     - Segids

     .. note::

        Atomtypes will be read from the input file if they are present
        (e.g. GROMACS PQR files). Otherwise, they will be guessed on Universe
        creation. By default, masses will also be guessed on Universe creation.
        This may change in release 3.0.
        See :ref:`Guessers` for more information.


    .. versionchanged:: 0.9.0
       Read chainID from a PQR file and use it as segid (before we always used
       'SYSTEM' as the new segid).
    .. versionchanged:: 0.16.1
       Now reads insertion codes and splits into new residues around these
    .. versionchanged:: 0.18.0
       Added parsing of Record types
       Can now read PQR files from Gromacs, these provide atom type as last column
       but don't have segids
    .. versionchanged:: 2.8.0
        Removed type and mass guessing (attributes guessing takes place now
        through universe.guess_TopologyAttrs() API).

    """

    format = "PQR"

    @staticmethod
    def guess_flavour(line):
        """Guess which variant of PQR format this line is

        Parameters
        ----------
        line : str
          entire line of PQR file starting with ATOM/HETATM

        Returns
        -------
        flavour : str
          ORIGINAL / GROMACS / NO_CHAINID

        .. versionadded:: 0.18.0
        """
        fields = line.split()
        if len(fields) == 11:
            try:
                float(fields[-1])
            except ValueError:
                flavour = "GROMACS"
            else:
                flavour = "ORIGINAL"
        else:
            flavour = "NO_CHAINID"
        return flavour

    def parse(self, **kwargs):
        """Parse atom information from PQR file *filename*.

        Returns
        -------
        A MDAnalysis Topology object
        """
        record_types = []
        serials = []
        names = []
        resnames = []
        chainIDs = []
        resids = []
        icodes = []
        charges = []
        radii = []
        elements = []

        flavour = None

        with openany(self.filename) as f:
            for line in f:
                if not line.startswith(("ATOM", "HETATM")):
                    continue
                fields = line.split()

                if flavour is None:
                    flavour = self.guess_flavour(line)
                if flavour == "ORIGINAL":
                    (
                        recordName,
                        serial,
                        name,
                        resName,
                        chainID,
                        resSeq,
                        x,
                        y,
                        z,
                        charge,
                        radius,
                    ) = fields
                elif flavour == "GROMACS":
                    (
                        recordName,
                        serial,
                        name,
                        resName,
                        resSeq,
                        x,
                        y,
                        z,
                        charge,
                        radius,
                        element,
                    ) = fields
                    chainID = "SYSTEM"
                    elements.append(element)
                elif flavour == "NO_CHAINID":
                    # files without the chainID
                    (
                        recordName,
                        serial,
                        name,
                        resName,
                        resSeq,
                        x,
                        y,
                        z,
                        charge,
                        radius,
                    ) = fields
                    chainID = "SYSTEM"

                try:
                    resid = int(resSeq)
                except ValueError:
                    # has icode present
                    resid = int(resSeq[:-1])
                    icode = resSeq[-1]
                else:
                    icode = ""

                record_types.append(recordName)
                serials.append(serial)
                names.append(name)
                resnames.append(resName)
                resids.append(resid)
                icodes.append(icode)
                charges.append(charge)
                radii.append(radius)
                chainIDs.append(chainID)

        n_atoms = len(serials)

        attrs = []
        if elements:
            atomtypes = elements
            attrs.append(Atomtypes(atomtypes, False))

        attrs.append(Atomids(np.array(serials, dtype=np.int32)))
        attrs.append(Atomnames(np.array(names, dtype=object)))
        attrs.append(Charges(np.array(charges, dtype=np.float32)))
        attrs.append(RecordTypes(np.array(record_types, dtype=object)))
        attrs.append(Radii(np.array(radii, dtype=np.float32)))

        resids = np.array(resids, dtype=np.int32)
        icodes = np.array(icodes, dtype=object)
        resnames = np.array(resnames, dtype=object)
        chainIDs = np.array(chainIDs, dtype=object)

        residx, (resids, resnames, icodes, chainIDs) = change_squash(
            (resids, resnames, icodes, chainIDs),
            (resids, resnames, icodes, chainIDs),
        )

        n_residues = len(resids)
        attrs.append(Resids(resids))
        attrs.append(Resnums(resids.copy()))
        attrs.append(Resnames(resnames))
        attrs.append(ICodes(icodes))

        segidx, chainIDs = squash_by(chainIDs)[:2]

        n_segments = len(chainIDs)
        attrs.append(Segids(chainIDs))

        top = Topology(
            n_atoms,
            n_residues,
            n_segments,
            attrs=attrs,
            atom_resindex=residx,
            residue_segindex=segidx,
        )

        return top
