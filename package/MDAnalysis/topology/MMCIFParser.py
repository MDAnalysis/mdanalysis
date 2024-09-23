"""
MMCIF Topology Parser # 
===================
"""

try:
    import gemmi
except ImportError:
    HAS_GEMMI = False
else:
    HAS_GEMMI = True

import itertools
import warnings

import numpy as np

from ..core.topology import Topology
from ..core.topologyattrs import (
    AltLocs,
    AtomAttr,
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
    ResidueAttr,
    Resnames,
    Resnums,
    Segids,
    SegmentAttr,
    Tempfactors,
)
from .base import TopologyReaderBase


def _into_idx(arr: list[int]) -> list[int]:
    return [idx for idx, (_, group) in enumerate(itertools.groupby(arr)) for _ in group]


def get_Atomattrs(model: gemmi.Model) -> tuple[list[AtomAttr], np.ndarray]:
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

    # transform *idx into continious numpy arrays
    residx = np.array(_into_idx(residx))

    # fill in altlocs
    altlocs = ["A" if not elem else elem for elem in altlocs]
    record_types = [
        "ATOM" if record == "A" else "HETATM" if record == "H" else None
        for record in record_types
    ]
    if any((elem is None for elem in record_types)):
        raise ValueError("Found an atom that is neither ATOM or HETATM")

    attrs = [
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
        Tempfactors(tempfactors),  # at.b_iso
    ]

    return attrs, residx


def get_Residueattrs(model: gemmi.Model) -> tuple[list[ResidueAttr], np.ndarray]:
    (
        icodes,  # res.seqid.icode
        resids,  # res.seqid.num
        resnames,  # res.name
        segidx,  # chain.name
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
    segidx = np.array(_into_idx(segidx))
    attrs = [
        Resnums(resnums),  # res.seqid.num
        ICodes(icodes),  # res.seqid.icode
        Resids(resids),  # res.seqid.num
        Resnames(resnames),  # res.name
    ]
    return attrs, segidx


def get_Segmentattrs(model: gemmi.Model) -> list[SegmentAttr]:
    segids = [chain.name for chain in model]
    return [Segids(segids)]


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

        atomAttrs, residx = get_Atomattrs(model)
        residAttrs, segidx = get_Residueattrs(model)
        segmentAttrs = get_Segmentattrs(model)

        attrs = atomAttrs + residAttrs + segmentAttrs

        n_atoms = len(atomAttrs[0])
        n_residues = len(residAttrs[0])
        n_segments = len(segmentAttrs[0])

        return Topology(
            n_atoms,
            n_residues,
            n_segments,
            attrs=attrs,
            atom_resindex=residx,
            residue_segindex=segidx,
        )
