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
MMTF Topology Parser
====================

Reads topology data from the `Macromolecular Transmission Format
(MMTF) format`_.  This should generally be a quicker alternative to PDB.

Makes individual models within the MMTF file available via the `models`
attribute on Universe.

.. versionadded:: 0.16.0

Reads the following topology attributes:

Atoms:
 - altLoc
 - atom ID
 - bfactor
 - bonds
 - charge
 - masses (guessed)
 - name
 - occupancy
 - type

Residues:
 - icode
 - resname
 - resid
 - resnum

Segments:
 - segid
 - model

Classes
-------

.. autoclass:: MMTFParser
   :members:

.. _Macromolecular Transmission Format (MMTF) format: https://mmtf.rcsb.org/
"""
from __future__ import absolute_import
from six.moves import zip

from collections import defaultdict
import mmtf
import numpy as np


from . import base
from . import guessers
from ..core.topology import Topology
from ..core.topologyattrs import (
    AltLocs,
    Atomids,
    Atomnames,
    Atomtypes,
    Bfactors,
    Bonds,
    Charges,
    ICodes,
    Masses,
    Occupancies,
    Resids,
    Resnames,
    Resnums,
    Segids,
    SegmentAttr,  # for model
)
from ..core.selection import RangeSelection


def _parse_mmtf(fn):
    if fn.endswith('gz'):
        return mmtf.parse_gzip(fn)
    else:
        return mmtf.parse(fn)


class Models(SegmentAttr):
    attrname = 'models'
    singular = 'model'
    transplants = defaultdict(list)

    def models(self):
        """Models in this Universe.

        The MMTF format can define various models for a given structure. The
        topology (eg residue identity) can change between different models,
        resulting in a different number of atoms in each model.

        Returns
        -------
        A list of AtomGroups, each representing a single model.
        """
        model_ids = np.unique(self.segments.models)

        return [self.select_atoms('model {}'.format(i))
                for i in model_ids]

    transplants['Universe'].append(
        ('models', property(models, None, None, models.__doc__)))


class ModelSelection(RangeSelection):
    token = 'model'
    field = 'models'

    def apply(self, group):
        mask = np.zeros(len(group), dtype=np.bool)
        vals = group.models

        for upper, lower in zip(self.uppers, self.lowers):
            if upper is not None:
                thismask = vals >= lower
                thismask &= vals <= upper
            else:
                thismask = vals == lower

            mask |= thismask
        return group[mask].unique


class MMTFParser(base.TopologyReaderBase):
    format = 'MMTF'

    def parse(self):
        if isinstance(self.filename, mmtf.MMTFDecoder):
            mtop = self.filename
        else:
            mtop = _parse_mmtf(self.filename)

        def iter_atoms(field):
            # iterate through atoms in groups
            for i in mtop.group_type_list:
                g = mtop.group_list[i]
                for val in g[field]:
                    yield val

        natoms = mtop.num_atoms
        nresidues = mtop.num_groups
        nsegments = mtop.num_chains
        attrs = []

        # Atom things
        # required
        altlocs = AltLocs(np.array([val.replace('\x00', '').strip()
                                    for val in mtop.alt_loc_list], dtype=object))
        atomids = Atomids(np.asarray(mtop.atom_id_list, dtype=np.int64))
        bfactors = Bfactors(np.asarray(mtop.b_factor_list, dtype=np.float64))
        charges = Charges(np.array(list(iter_atoms('formalChargeList')), dtype=np.float64))
        names = Atomnames(np.array(list(iter_atoms('atomNameList')), dtype=object))
        occupancies = Occupancies(np.asarray(mtop.occupancy_list, dtype=np.float64))
        types = Atomtypes(np.array(list(iter_atoms('elementList')), dtype=object))

        masses = Masses(guessers.guess_masses(types.values), guessed=True)

        attrs.extend([altlocs, atomids, bfactors, charges, names, occupancies, types, masses])

        # Residue things
        # required
        resids = Resids(np.asarray(mtop.group_id_list, dtype=np.int64))
        resnums = Resnums(resids.values.copy())
        resnames = Resnames(np.array([mtop.group_list[i]['groupName']
                                      for i in mtop.group_type_list], dtype=object))
        # mmtf empty icode is '\x00' rather than ''
        icodes = ICodes(np.array([val.replace('\x00', '').strip()
                                 for val in mtop.ins_code_list], dtype=object))
        attrs.extend([resids, resnums, resnames, icodes])
        # optional

        # Segment things
        segids = Segids(np.asarray(mtop.chain_name_list, dtype=object))
        #chainids = ChainIDs(np.asarray(mtop.chain_id_list, dtype=object))

        attrs.append(segids)

        mods = np.repeat(np.arange(mtop.num_models), mtop.chains_per_model)
        attrs.append(Models(mods))
        #attrs.append(chainids)

        # number of atoms in a given group id
        groupID_2_natoms = {i:len(g['atomNameList'])
                            for i, g in enumerate(mtop.group_list)}
        # mapping of atoms to residues
        resindex = np.repeat(np.arange(nresidues),
                             [groupID_2_natoms[i] for i in mtop.group_type_list])

        # mapping of residues to segments
        segindex = np.repeat(np.arange(nsegments), mtop.groups_per_chain)

        # Bonds
        # bonds are listed as indices within a group,
        # offset pulls out 'global' index
        offset = 0
        bonds = []
        for gtype in mtop.group_type_list:
            g = mtop.group_list[gtype]
            bondlist = g['bondAtomList']

            for x, y in zip(bondlist[1::2], bondlist[::2]):
                if x > y:
                    x, y = y, x  # always have x < y
                bonds.append((x + offset, y + offset))

            offset += groupID_2_natoms[gtype]
        # add inter group bonds
        for x, y in zip(mtop.bond_atom_list[1::2], mtop.bond_atom_list[::2]):
            if x > y:
                x, y = y, x
            bonds.append((x, y))
        attrs.append(Bonds(bonds))

        top = Topology(natoms, nresidues, nsegments,
                       atom_resindex=resindex,
                       residue_segindex=segindex,
                       attrs=attrs)

        return top
