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
            altlocs,  # at.altlog
            atomtypes,  # at.name
            elements,  # at.element.name
            formalcharges,  # at.charge
            names,  # at.name
            serials,  # at.serial
            tempfactors,  # at.b_iso
            occupancies,  # at.occ
            weights,  # at.element.weight
            record_types,  # res.het_flag TODO: match to ATOM/HETATM
            chainids,  # chain.name
            resids,  # res.seqid.num
            resnames,  # res.name
            icodes,  # residue.seqid.icode
        ) = map(  # this is probably not pretty, but it's efficient -- one loop over the mmcif
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
                            res.het_flag,
                            chain.name,  # chainids
                            res.seqid.num,
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
            AltLocs(altlocs),
            Atomids(serials),
            Atomnames(names),
            Atomtypes(atomtypes),
            ChainIDs(chainids),  # actually for atom
            Elements(elements),
            FormalCharges(formalcharges),
            Masses(weights),
            Occupancies(occupancies),
            RecordTypes(record_types),  # for atom too?
            Tempfactors(tempfactors),
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
