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
#
"""
PDBx topology parser
====================


See Also
--------
:class:`MDAnalysis.coordinates.PDBx`

"""
import gemmi
import numpy as np

from .base import TopologyReaderBase, change_squash
from ..core.topology import Topology
from ..core.topologyattrs import (
    Atomnames,
    Atomids,
    AltLocs,
    Elements,
    ICodes,
    RecordTypes,
    Resids,
    Resnames,
    Segids,
)


class PDBxParser(TopologyReaderBase):
    """Read a Topology from a PDBx file

    Creates the following attributes from these "_atom_site" PDBx loop entries
    - "group_PDB" RecordType
    - "id" AtomId
    - "label_alt_id" AltLoc
    - "label_type_symbol" Element
    - "label_atom_id" AtomName
    - "auth_seq_id" Resid
    - "auth_comp_id" Resname
    - "pdbx_PDB_ins_code" ICode
    - "auth_asym_id" ChainID
    """
    format = ['PBDx', 'cif']

    def parse(self, **kwargs) -> Topology:
        doc = gemmi.cif.read(self.filename)
        block = doc.sole_block()

        attrs = []

        def objarr(x):
            return np.array(x, dtype=object)

        # hierarchy correspondence:
        # seq_id -> residues
        # entity_id -> chains
        if recordtypes := block.find('_atom_site.group_PDB'):
            attrs.append(RecordTypes(recordtypes))
        ids = block.find_loop('_atom_site.id')
        n_atoms = len(ids)
        attrs.append(Atomids(ids))
        if altlocs := block.find_loop('_atom_site.label_alt_id'):
            altlocs = np.array(altlocs, dtype=object)
            altlocs[altlocs == '.'] = ''
            attrs.append(AltLocs(altlocs))
        if elements_loop := block.find_loop('_atom_site.type_symbol'):
            attrs.append(Elements(objarr(elements_loop)))
        if names_loop := block.find_loop('_atom_site.label_atom_id'):
            attrs.append(Atomnames(objarr(names_loop)))

        # sort out residues/segments
        # label_seq_id seems to not cover entire model unlike author versions
        resids = block.find_loop('_atom_site.auth_seq_id')
        resnames = block.find_loop('_atom_site.auth_comp_id')
        icodes = block.find_loop('_atom_site.pdbx_PDB_ins_code')
        chainids = block.find_loop('_atom_site.auth_asym_id')

        residx, (resids, icodes, resnames, chainids) = change_squash(
            (resids, icodes), (resids, icodes, resnames, chainids)
        )
        segidx, (chainids,) = change_squash((chainids,), (chainids,))

        attrs.extend((
            Resids(resids),
            Resnames(objarr(resnames)),
            ICodes(objarr(icodes)),
            Segids(chainids),
        ))

        n_residues = len(resids)
        n_segments = len(chainids)

        return Topology(
            n_atoms=n_atoms,
            n_res=n_residues,
            n_seg=n_segments,
            attrs=attrs,
            atom_resindex=residx,
            residue_segindex=segidx,
        )
