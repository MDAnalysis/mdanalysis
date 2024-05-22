"""
MMCIF Topology Parser # TODO
===================
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
                            at.name,  # atomtypes
                            at.element.name,  # elements
                            at.charge,  # formalcharges
                            at.name,  # names
                            at.serial,  # serials
                            at.b_iso,  # tempfactores
                            at.occ,  # occupancies
                            at.element.weight,  # weights
                            # --
                            # residue properties
                            # --
                            res.seqid.icode,  # icodes
                            res.het_flag,  # record_types
                            res.seqid.num,  # resids
                            res.name,  # resnames
                            res.segment,  # segids
                            # --
                            # chain properties
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
        _, (seg_segids,) = change_squash((res_segids,), (res_segids,))

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
            ICodes(res_icodes),  # for each atom
            RecordTypes(record_types),  # for atom too?
            Resids(res_resids),  # for residue
            Resnames(res_resnames),  # for residue
            #
            Segids(seg_segids),  # for segment (currently for residue)
            # per chain
            ChainIDs(chainids),  # actually for atom
        ]

        n_atoms = len(names)
        n_residues = len(res_resids)
        n_segments = len(seg_segids)
        top = Topology(n_atoms, n_residues, n_segments, attrs=attrs)

        return top
