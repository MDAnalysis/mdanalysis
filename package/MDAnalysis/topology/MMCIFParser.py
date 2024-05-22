# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
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
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
# TODO

"""
MMCIF Topology Parser # TODO
===================

.. versionadded:: 0.9.1

Reads an xyz file and pulls the atom information from it.  Because
xyz only has atom name information, all information about residues
and segments won't be populated.

Classes
-------

.. autoclass:: XYZParser
   :members:

"""

from .base import TopologyReaderBase
from ..core.topology import Topology
from ..core.topologyattrs import (
    Atomnames,
    Atomids,
    AltLocs,
    ChainIDs,
    Atomtypes,
    Elements,
    # ICodes,
    Masses,
    Occupancies,
    RecordTypes,
    Resids,
    Resnames,
    # Segids,
    Tempfactors,
    FormalCharges,
)


class MMCIFParser(TopologyReaderBase):
    """Parse a list of atoms from an XYZ file.

    Creates the following attributes:
     - Atomnames
     - TODO


    .. versionadded:: 2.8.0
    """

    format = "MMCIF"

    def parse(self, **kwargs):
        """Read the file and return the structure.

        Returns
        -------
        MDAnalysis Topology object
        """
        import gemmi

        structure = gemmi.read_structure(self.filename)

        # here we freaking go

        (
            names,
            atomtypes,
            record_types,
            serials,
            altlocs,
            chainids,
            # icodes,
            tempfactors,
            occupancies,
            resids,
            resnames,
            # segids,
            elements,
            formalcharges,
            weights,
        ) = list(
            zip(
                *[
                    (
                        at.name,  # names
                        at.element.name,  # atomtypes
                        res.het_flag,  # record_types
                        at.serial,  # serials
                        at.altloc,  # altlocs
                        chain.name,  # chainids
                        # res.seqid.icode,  # icodes
                        at.b_iso,  # tempfactores
                        at.occ,  # occupancies
                        res.seqid.num,  # resids
                        res.name,  # resnames
                        # res.segment,  # segids
                        at.element.name,  # elements
                        at.charge,  # formalcharges
                        at.element.weight,  # weights
                    )
                    for model in structure
                    for chain in model
                    for res in chain
                    for at in res
                ]
            )
        )

        attrs = [
            Atomnames(names),
            Atomtypes(atomtypes),
            RecordTypes(record_types),
            Atomids(serials),
            AltLocs(altlocs),
            ChainIDs(chainids),
            # ICodes(icodes),
            Tempfactors(tempfactors),
            Occupancies(occupancies),
            Resids(resids),
            Resnames(resnames),
            # Segids(segids),
            Elements(elements),
            FormalCharges(formalcharges),
            Masses(weights),
        ]

        n_atoms = len(names)
        n_residues = len(resids)
        n_segments = 1
        top = Topology(n_atoms, n_residues, n_segments, attrs=attrs)

        return top
