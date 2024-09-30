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


def _into_idx(arr: list) -> list[int]:
    """Replace consecutive identical elements of an array with their indices.

    Example
    -------
    .. code-block:: python

        arr: list[int] = [1, 1, 5, 5, 7, 3, 3]
        assert _into_idx(arr) == [0, 0, 1, 1, 2, 3, 3]

    Parameters
    ----------
    arr
        input array of elements that can be compared with `__eq__`

    Returns
    -------
        list[int] -- array where these elements are replaced with their unique indices, in order of appearance.
    """
    return [idx for idx, (_, group) in enumerate(itertools.groupby(arr)) for _ in group]


def get_Atomattrs(model: gemmi.Model) -> tuple[list[AtomAttr], np.ndarray]:
    """Extract all attributes that are subclasses of :class:`..core.topologyattrs.AtomAttr` from a  ``gemmi.Model`` object,
    and a `segidx` index.

    Parameters
    ----------
    model
        input `gemmi.Model`, e.g. `gemmi.read_structure('file.cif')[0]`

    Returns
    -------
        tuple[list[AtomAttr], np.ndarray] -- first element is list of all extracted attributes, second element is `segidx`

    Raises
    ------
    ValueError
        if any of the records is neither 'ATOM' nor 'HETATM'
    """
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
    ) = map(  # this construct takes np.ndarray of all lists of attributes, extracted from the `gemmi.Model`
        np.array,
        list(
            zip(
                *[
                    (
                        # tuple of attributes
                        # extracted from residue, atom or chain in the structure
                        # ------------------
                        atom.altloc,  # altlocs
                        atom.serial,  # serials
                        atom.name,  # names
                        atom.name,  # atomtypes
                        # ------------------
                        chain.name,  # chainids
                        atom.element.name,  # elements
                        atom.charge,  # formalcharges
                        atom.element.weight,  # weights
                        # ------------------
                        atom.occ,  # occupancies
                        residue.het_flag,  # record_types
                        atom.b_iso,  # tempfactors
                        residue.seqid.num,  # residx, later translated into continious repr
                    )
                    # the main loop over the `gemmi.Model` object
                    for chain in model
                    for residue in chain
                    for atom in residue
                ]
            )
        ),
    )

    # transform *idx into continious numpy arrays
    residx = np.array(_into_idx(residx))

    # fill in altlocs, since gemmi has '' as default
    altlocs = ["A" if not elem else elem for elem in altlocs]

    # convert default gemmi record types to default MDAnalysis record types
    record_types = [
        "ATOM" if record == "A" else "HETATM" if record == "H" else None
        for record in record_types
    ]
    if any((elem is None for elem in record_types)):
        raise ValueError("Found an atom that is neither ATOM or HETATM")

    attrs = [
        AltLocs(altlocs),
        Atomids(serials),
        Atomnames(names),
        Atomtypes(atomtypes),
        # ----------------------------
        ChainIDs(chainids),
        Elements(elements),
        FormalCharges(formalcharges),
        Masses(weights),
        # ----------------------------
        Occupancies(occupancies),
        RecordTypes(record_types),
        Tempfactors(tempfactors),
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
