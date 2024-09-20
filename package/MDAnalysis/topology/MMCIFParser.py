"""
MMCIF Topology Parser # TODO
===================
"""

import gemmi
import numpy as np
import warnings

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
    Tempfactors,
    Segids,
)
from .base import TopologyReaderBase, change_squash


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
            record_types,  # res.het_flag TODO: match to ATOM/HETATM
            tempfactors,  # at.b_iso
            resids,  # res.seqid.num
            # ------------------
            resnames,  # res.name
            icodes,  # residue.seqid.icode
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
                            res.seqid.num,
                            # ------------------
                            res.name,
                            res.seqid.icode,
                        )
                        for chain in model
                        for res in chain
                        for at in res
                    ]
                )
            ),
        )

        # fill in altlocs
        altlocs = ["A" if not elem else elem for elem in altlocs]

        attrs = [
            # AtomAttr subclasses
            AltLocs(altlocs),  # ✅
            Atomids(serials),  # ✅
            Atomnames(names),  # ✅
            Atomtypes(atomtypes),  # ✅
            # ------------------
            ChainIDs(chainids),  # ✅
            Elements(elements),  # ✅; same as atomtypes
            FormalCharges(formalcharges),  # ✅
            Masses(weights),  # ✅
            # ------------------
            Occupancies(occupancies),  # ✅
            RecordTypes(record_types),  # ✅
            Tempfactors(tempfactors),  # ✅
            # ResidueAttr subclasses
            ICodes(icodes),  # for each atom
            Resids(resids),  # for residue
            Resnames(resnames),  # for residue
            # SegmentAttr subclasses
            Segids(segids),
        ]

        n_atoms = len(names)
        n_residues = len(resids)
        n_segments = len(set(chainids))

        print(resids)
        print(change_squash(resids, resids))

        top = Topology(n_atoms, n_residues, n_segments, attrs=attrs)

        print(f"{n_atoms=}")
        print(f"{n_residues=}")
        print(f"{n_segments=}")

        return top
