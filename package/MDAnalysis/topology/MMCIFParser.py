"""
MMCIF Topology Parser # TODO
===================
"""

import gemmi
import numpy as np
import warnings
import itertools

from ..core.topology import Topology
from ..core.topologyattrs import (
    AltLocs,
    Atomids,
    Atomnames,
    Atomtypes,
    ChainIDs,
    Elements,
    FormalCharges,
    ICodes,
    Masses,
    Occupancies,
    RecordTypes,
    Resids,
    Resnames,
    Resnums,
    Segids,
    Tempfactors,
)
from .base import TopologyReaderBase


def _into_idx(arr: list[int]) -> list[int]:
    return [idx for idx, (_, group) in enumerate(itertools.groupby(arr)) for _ in group]


class MMCIFParser(TopologyReaderBase):
    format = "MMCIF"

    def parse(self, **kwargs):
        """Read the file and return the structure.

        Returns
        -------
        MDAnalysis Topology object
        """
        structure = gemmi.read_structure(self.filename)

        if len(structure) > 1:
            warnings.warn(
                "MMCIF model {self.filename} contains {len(model)=} different models, "
                "but only the first one will be used to assign the topology"
            )
        model = structure[0]

        # atom properties
        (
            altlocs,  # at.altloc
            serials,  # at.serial
            names,  # at.name
            atomtypes,  # at.name
            # ------------------
            chainids,  # chain.name
            elements,  # at.element.name
            formalcharges,  # at.charge
            weights,  # at.element.weight
            # ------------------
            occupancies,  # at.occ
            record_types,  # res.het_flag
            tempfactors,  # at.b_iso
            residx,  # _into_idx(res.seqid.num)
        ) = map(  # this is probably not pretty, but it's efficient -- one loop over the mmcif
            np.array,
            list(
                zip(
                    *[
                        (
                            at.altloc,  # altlocs
                            at.serial,  # serials
                            at.name,  # names
                            at.name,  # atomtypes
                            # ------------------
                            chain.name,  # chainids
                            at.element.name,  # elements
                            at.charge,  # formalcharges
                            at.element.weight,  # weights
                            # ------------------
                            at.occ,  # occupancies
                            res.het_flag,  # record_types
                            at.b_iso,  # tempfactors
                            res.seqid.num,  # residx, later translated into continious repr
                        )
                        for chain in model
                        for res in chain
                        for at in res
                    ]
                )
            ),
        )

        (
            icodes,  # res.seqid.icode
            resids,  # res.seqid.num
            resnames,  # res.name
            segidx,  # chain.name TODO: translate into continious index
            resnums,
        ) = map(
            np.array,
            list(
                zip(
                    *[
                        (
                            res.seqid.icode,
                            res.seqid.num,
                            res.name,
                            chain.name,
                            res.seqid.num,
                        )
                        for chain in model
                        for res in chain
                    ]
                )
            ),
        )

        segids = [chain.name for chain in model]

        # transform *idx into continious numpy arrays
        residx = np.array(_into_idx(residx))
        segidx = np.array(_into_idx(segidx))

        # fill in altlocs
        altlocs = ["A" if not elem else elem for elem in altlocs]
        record_types = [
            "ATOM" if record == "A" else "HETATM" if record == "H" else None
            for record in record_types
        ]
        if any((elem is None for elem in record_types)):
            raise ValueError("Found an atom that is neither ATOM or HETATM")

        attrs = [
            # AtomAttr subclasses
            AltLocs(altlocs),  # at.altloc
            Atomids(serials),  # at.serial
            Atomnames(names),  # at.name
            Atomtypes(atomtypes),  # at.name
            # ---------------------------------------
            ChainIDs(chainids),  # chain.name
            Elements(elements),  # at.element.name
            FormalCharges(formalcharges),  # at.charge
            Masses(weights),  # at.element.weight
            # ---------------------------------------
            Occupancies(occupancies),  # at.occ
            RecordTypes(record_types),  # res.het_flat
            Resnums(resnums),  # res.seqid.num
            Tempfactors(tempfactors),  # at.b_iso
            #
            # ResidueAttr subclasses
            ICodes(icodes),  # res.seqid.icode
            Resids(resids),  # res.seqid.num
            Resnames(resnames),  # res.name
            #
            # SegmentAttr subclasses
            Segids(segids),  # chain.name
        ]

        n_atoms = len(names)
        n_residues = len(resids)
        n_segments = len(segids)

        print(f"{len(residx)=}")
        print(f"{residx=}")
        print(f"{len(segidx)=}")
        print(f"{segidx=}")

        print(f"{n_atoms=}")
        print(f"{n_residues=}")
        print(f"{n_segments=}")

        return Topology(
            n_atoms,
            n_residues,
            n_segments,
            attrs=attrs,
            atom_resindex=residx,
            residue_segindex=segidx,
        )
