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

import numpy as np
from .base import TopologyReaderBase, change_squash
from ..core.topology import Topology
from ..core.topologyattrs import (
    Atomnames,
    Atomids,
    AltLocs,
    ChainIDs,
    Atomtypes,
    Elements,
    ICodes,
    Masses,
    Occupancies,
    RecordTypes,
    Resids,
    Resnames,
    Segids,
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
            # atom properties -- no squashing!
            # --
            altlocs,
            atomtypes,
            elements,
            formalcharges,
            names,
            serials,
            tempfactors,
            occupancies,
            weights,
            # --
            # residue properties -- some squashing...
            # --
            icodes,
            record_types,
            resids,
            resnames,
            segids,
            # --
            # chain properties -- lots of squashing...
            # --
            chainids,
        ) = map(
            np.array,
            list(
                zip(
                    *[
                        (
                            # atom properties -- no squashing!
                            # --
                            at.altloc,  # altlocs
                            at.element.name,  # atomtypes
                            at.element.name,  # elements
                            at.charge,  # formalcharges
                            at.name,  # names
                            at.serial,  # serials
                            at.b_iso,  # tempfactores
                            at.occ,  # occupancies
                            at.element.weight,  # weights
                            # --
                            # residue properties -- some squashing...
                            # --
                            res.seqid.icode,  # icodes
                            res.het_flag,  # record_types
                            res.seqid.num,  # resids
                            res.name,  # resnames
                            res.segment,  # segids
                            # --
                            # chain properties -- lots of squashing...
                            # --
                            chain.name,  # chainids
                        )
                        for model in structure
                        for chain in model
                        for res in chain
                        for at in res
                    ]
                )
            ),
        )

        # squash residue-based attributes
        _, (res_icodes, res_record_types, res_resids, res_resnames, res_segids) = (
            change_squash(
                (resids, resnames),
                (icodes, record_types, resids, resnames, segids),
            )
        )
        # squash chain-based attributes
        _, (chain_chainids,) = change_squash((chainids,), (chainids,))

        attrs = [
            # per atom
            AltLocs(altlocs),
            Atomtypes(atomtypes),
            Elements(elements),
            FormalCharges(formalcharges),
            Atomnames(names),
            Atomids(serials),
            Tempfactors(tempfactors),
            Occupancies(occupancies),
            Masses(weights),
            # per residue
            # ICodes(res_icodes),
            # RecordTypes(res_record_types),
            Resids(res_resids),
            Resnames(res_resnames),
            Segids(res_segids),
            # per chain
            ChainIDs(chainids),
        ]

        n_atoms = len(names)
        n_residues = len(res_resids)
        n_segments = len(res_segids)
        top = Topology(n_atoms, n_residues, n_segments, attrs=attrs)

        return top
