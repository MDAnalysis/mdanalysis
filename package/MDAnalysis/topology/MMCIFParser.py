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

        # atom properties
        (
            altlocs,
            atomtypes,
            elements,
            formalcharges,
            names,
            serials,
            tempfactors,
            occupancies,
            weights,
        ) = map(
            np.array,
            list(
                zip(
                    *[
                        (
                            at.altloc,  # altlocs
                            at.name,  # atomtypes
                            at.element.name,  # elements
                            at.charge,  # formalcharges
                            at.name,  # names
                            at.serial,  # serials
                            at.b_iso,  # tempfactores
                            at.occ,  # occupancies
                            at.element.weight,  # weights
                        )
                        for model in structure
                        for chain in model
                        for res in chain
                        for at in res
                    ]
                )
            ),
        )
        # per-residue properties
        (
            icodes,
            record_types,
            resids,
            resnames,
            segids,
        ) = map(
            np.array,
            list(
                zip(
                    *[
                        (
                            res.seqid.icode,  # icodes
                            res.het_flag,  # record_types
                            res.seqid.num,  # resids
                            res.name,  # resnames
                            res.segment,  # segids
                        )
                        for model in structure
                        for chain in model
                        for res in chain
                    ]
                )
            ),
        )

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
            ICodes(icodes),  # for each atom
            RecordTypes(record_types),  # for atom too?
            Resids(resids),  # for residue
            Resnames(resnames),  # for residue
            #
            # Segids(segids),  # for segment (currently for residue)
            # per chain
            # ChainIDs(chainids),  # actually for atom
        ]

        n_atoms = len(names)
        n_residues = len(resids)
        n_segments = len(segids)
        top = Topology(n_atoms, n_residues, n_segments, attrs=attrs)

        return top
