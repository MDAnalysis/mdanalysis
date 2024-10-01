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

    .. versionadded:: 2.8.0
    """
    return [idx for idx, (_, group) in enumerate(itertools.groupby(arr)) for _ in group]


def get_Atomattrs(model: gemmi.Model) -> tuple[list[AtomAttr], np.ndarray]:
    """Extract all attributes that are subclasses of :class:`..core.topologyattrs.AtomAttr` from a  ``gemmi.Model`` object,
    and a `residx` index with indices of all atoms in residues.

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

    .. versionadded:: 2.8.0
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
    """Extract all attributes that are subclasses of :class:`..core.topologyattrs.ResidueAttr` from a  ``gemmi.Model`` object,
    and a `segidx` index witn indices of all residues in segments.

    Parameters
    ----------
    model
        input `gemmi.Model`, e.g. `gemmi.read_structure('file.cif')[0]`

    Returns
    -------
        tuple[list[ResidueAttr], np.ndarray] -- first element is list of all extracted attributes, second element is `segidx`

    .. versionadded:: 2.8.0
    """
    (
        icodes,  # residue.seqid.icode
        resids,  # residue.seqid.num
        resnames,  # residue.name
        segidx,  # chain.name
        resnums,
    ) = map(
        np.array,
        list(
            zip(
                *[
                    (
                        residue.seqid.icode,
                        residue.seqid.num,
                        residue.name,
                        chain.name,
                        residue.seqid.num,
                    )
                    for chain in model
                    for residue in chain
                ]
            )
        ),
    )
    segidx = np.array(_into_idx(segidx))
    attrs = [
        Resnums(resnums),
        ICodes(icodes),
        Resids(resids),
        Resnames(resnames),
    ]
    return attrs, segidx


def get_Segmentattrs(model: gemmi.Model) -> SegmentAttr:
    """Extract all attributes that are subclasses of :class:`..core.topologyattrs.SegmentAttr` from a  ``gemmi.Model`` object.

    Parameters
    ----------
    model
        input `gemmi.Model`, e.g. `gemmi.read_structure('file.cif')[0]`

    Returns
    -------
        list[SegmentAttr] -- list of all extracted attributes

    .. versionadded:: 2.8.0
    """
    segids = [chain.name for chain in model]
    return [Segids(segids)]


class MMCIFParser(TopologyReaderBase):
    """Parser that obtains a list of atoms from a standard MMCIF/PDBx file using ``gemmi`` library (https://github.com/project-gemmi/gemmi).

    Creates the following Attributes (if present):
        - :class:`..core.topologyattrs.AtomAttr` subclasses:
            - :class:`..core.topologyattrs.AltLocs`
            - :class:`..core.topologyattrs.Atomids`
            - :class:`..core.topologyattrs.Atomnames`
            - :class:`..core.topologyattrs.Atomtypes`
            - :class:`..core.topologyattrs.ChainIDs`
            - :class:`..core.topologyattrs.Elements`
            - :class:`..core.topologyattrs.FormalCharges`
            - :class:`..core.topologyattrs.Masses`
            - :class:`..core.topologyattrs.Occupancies`
            - :class:`..core.topologyattrs.RecordTypes`
            - :class:`..core.topologyattrs.Tempfactors`
        - :class:`..core.topologyattrs.ResidueAttr` subclasses:
            - :class:`..core.topologyattrs.Resnums`
            - :class:`..core.topologyattrs.ICodes`
            - :class:`..core.topologyattrs.Resids`
            - :class:`..core.topologyattrs.Resnames`
        - :class:`..core.topologyattrs.SegmentAttr` subclasses:
            - :class:`..core.topologyattrs.Segids`

    .. versionadded:: 2.8.0
    """

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

        # due to the list(map(...)) construction, all elements in array have equal lengths
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
