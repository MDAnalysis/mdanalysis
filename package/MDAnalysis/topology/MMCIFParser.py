# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
"""
MMCIF Topology Parser
=====================

.. versionadded:: 2.11.0

Read topology information from mmCIF/PDBx coordinate files using the
`Gemmi library <https://gemmi.readthedocs.io>`_. PDB files can be read the
same way with ``format="pdb_gemmi"``; see :ref:`PDB-gemmi-backend` for how
this differs from the default :mod:`~MDAnalysis.topology.PDBParser`.



mmCIF files contain topology information about the molecules in the
structure. For each atom the following attributes are read and stored
in the relevant topology attributes:

    - :class:`MDAnalysis.core.topologyattrs.AtomAttr` subclasses:
        - :class:`MDAnalysis.core.topologyattrs.AltLocs`
        - :class:`MDAnalysis.core.topologyattrs.Atomids`
        - :class:`MDAnalysis.core.topologyattrs.Atomnames`
        - :class:`MDAnalysis.core.topologyattrs.Atomtypes`
        - :class:`MDAnalysis.core.topologyattrs.ChainIDs`
        - :class:`MDAnalysis.core.topologyattrs.Elements`
        - :class:`MDAnalysis.core.topologyattrs.FormalCharges`
        - :class:`MDAnalysis.core.topologyattrs.Masses`
        - :class:`MDAnalysis.core.topologyattrs.Occupancies`
        - :class:`MDAnalysis.core.topologyattrs.RecordTypes`
        - :class:`MDAnalysis.core.topologyattrs.Tempfactors`
    - :class:`MDAnalysis.core.topologyattrs.ResidueAttr` subclasses:
        - :class:`MDAnalysis.core.topologyattrs.Resnums`
        - :class:`MDAnalysis.core.topologyattrs.ICodes`
        - :class:`MDAnalysis.core.topologyattrs.Resids`
        - :class:`MDAnalysis.core.topologyattrs.Resnames`
    - :class:`MDAnalysis.core.topologyattrs.SegmentAttr` subclasses:
        - :class:`MDAnalysis.core.topologyattrs.Segids`

Classes
-------

.. autoclass:: MMCIFParser
   :members:
   :inherited-members:

"""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from gemmi import Structure

import logging
import warnings

import numpy as np

from ..coordinates.MMCIF import HAS_GEMMI, _read_gemmi_structure
from ..core.topology import Topology
from ..guesser.tables import SYMB2Z
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
from .base import TopologyReaderBase, change_squash

logger = logging.getLogger("MDAnalysis.topology.MMCIFParser")


def _factorize(column: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Split a fixed-width bytes column into an index per entry and the
    distinct decoded strings, sorting the bytes as unsigned integers"""
    column = np.ascontiguousarray(column)
    codes, nmidx = np.unique(
        column.view(f"u{column.dtype.itemsize}"), return_inverse=True
    )
    return nmidx, codes.view(column.dtype).astype(str).astype(object)


def _dedupe(nmidx: np.ndarray, values) -> tuple[np.ndarray, np.ndarray]:
    """Merge distinct values that became equal after being transformed"""
    values, remap = np.unique(
        np.asarray(values, dtype=object), return_inverse=True
    )
    return remap[nmidx], values


class MMCIFParser(TopologyReaderBase):
    """Parser that obtains a list of atoms from a standard MMCIF/PDBx file using
    the `gemmi library <https://gemmi.readthedocs.io>`_.

    The *filename* argument accepts a file path, a compressed ``.cif.gz`` file,
    or a stream/file-like object.

    Creates the following Attributes (if present):
        - :class:`MDAnalysis.core.topologyattrs.AtomAttr` subclasses:
            - :class:`MDAnalysis.core.topologyattrs.AltLocs`
            - :class:`MDAnalysis.core.topologyattrs.Atomids`
            - :class:`MDAnalysis.core.topologyattrs.Atomnames`
            - :class:`MDAnalysis.core.topologyattrs.Atomtypes`
            - :class:`MDAnalysis.core.topologyattrs.ChainIDs`
            - :class:`MDAnalysis.core.topologyattrs.Elements`
            - :class:`MDAnalysis.core.topologyattrs.FormalCharges`
            - :class:`MDAnalysis.core.topologyattrs.Masses`
            - :class:`MDAnalysis.core.topologyattrs.Occupancies`
            - :class:`MDAnalysis.core.topologyattrs.RecordTypes`
            - :class:`MDAnalysis.core.topologyattrs.Tempfactors`
        - :class:`MDAnalysis.core.topologyattrs.ResidueAttr` subclasses:
            - :class:`MDAnalysis.core.topologyattrs.Resnums`
            - :class:`MDAnalysis.core.topologyattrs.ICodes`
            - :class:`MDAnalysis.core.topologyattrs.Resids`
            - :class:`MDAnalysis.core.topologyattrs.Resnames`
        - :class:`MDAnalysis.core.topologyattrs.SegmentAttr` subclasses:
            - :class:`MDAnalysis.core.topologyattrs.Segids`
    """

    format = ["cif", "cif.gz", "mmcif", "mmcif.gz", "pdb_gemmi"]

    def __init__(self, filename):
        if not HAS_GEMMI:
            errmsg = (
                "MMCIFParser: To read a Topology from an mmCIF file, "
                "please install gemmi"
            )
            raise ImportError(errmsg)
        super(MMCIFParser, self).__init__(filename)

    def parse(self, **kwargs) -> Topology:
        """Read the file and return the structure.

        Returns
        -------
        MDAnalysis Topology object
        """
        structure = self._get_structure()

        if len(structure) > 1:
            wmsg = (
                f"MMCIF model {self.filename} contains {len(structure)} different models, "
                "but only the first one will be used to assign the topology"
            )
            warnings.warn(wmsg)
            logger.warning(wmsg)
            del structure[1:]

        import gemmi

        flat = gemmi.FlatStructure(structure)
        flat.strings_as_numbers = False

        # String columns are factorized in bulk into (per-atom index, distinct
        # values) pairs, the layout MDAnalysis interns strings in anyway
        het_idx, het_flags = _factorize(np.asarray(flat.het_flags).view("S1"))
        # Intentionally strict: this also rejects the "\0" (unspecified)
        # het_flag, so an mmCIF lacking the ``group_PDB`` field will fail
        # to parse.
        if not set(het_flags) <= {"A", "H"}:
            raise ValueError("Found an atom that is neither ATOM nor HETATM")
        record_types = np.where(het_flags == "A", "ATOM", "HETATM")

        # gemmi stores "no altloc" as "\0", which numpy decodes to ""
        altlocs = _factorize(np.asarray(flat.altlocs).view("S1"))
        names = _factorize(np.asarray(flat.atom_names))
        chain_idx, chainids = _factorize(np.asarray(flat.chain_ids))
        elem_idx, elements = _factorize(np.asarray(flat.element_names))
        resname_idx, resnames = _factorize(np.asarray(flat.residue_names))
        icode_idx, icodes = _factorize(np.asarray(flat.icodes).view("S1"))
        icode_idx, icodes = _dedupe(icode_idx, [c.strip() for c in icodes])
        resids = np.asarray(flat.resnums)

        # As in PDBParser: feed atomtypes the raw element column, but
        # validate elements against known symbols before storing them
        validated_elements = []
        for elem in elements:
            if elem.capitalize() in SYMB2Z:
                validated_elements.append(elem.capitalize())
            else:
                wmsg = (
                    f"Unknown element {elem} found for some atoms. "
                    f"These have been given an empty element record. "
                    f"If needed they can be guessed using "
                    f"universe.guess_TopologyAttrs(context='default',"
                    " to_guess=['elements'])."
                )
                warnings.warn(wmsg)
                validated_elements.append("")

        # Atom Attributes
        attrs = [
            AltLocs._from_lookup(*altlocs),
            Atomids(np.asarray(flat.serials)),
            Atomnames._from_lookup(*names),
            Atomtypes._from_lookup(elem_idx, elements),
            # ----------------------------
            ChainIDs._from_lookup(chain_idx, chainids),
            Elements._from_lookup(*_dedupe(elem_idx, validated_elements)),
            FormalCharges(np.asarray(flat.charge)),
            Masses(np.asarray(flat.element_weights)),
            # ----------------------------
            Occupancies(np.asarray(flat.occ)),
            RecordTypes._from_lookup(het_idx, record_types),
            Tempfactors(np.asarray(flat.b_iso)),
        ]
        n_atoms = len(resids)

        # Residue Attributes
        residx, (resids, resname_idx, icode_idx, chain_idx) = change_squash(
            (resids, resname_idx, icode_idx, chain_idx),
            (resids, resname_idx, icode_idx, chain_idx),
        )
        attrs.append(Resids(resids))
        attrs.append(Resnames._from_lookup(resname_idx, resnames))
        attrs.append(Resnums(resids.copy()))
        attrs.append(ICodes._from_lookup(icode_idx, icodes))
        n_residues = len(resids)

        # Segment Attributes
        segidx, (chain_idx,) = change_squash((chain_idx,), (chain_idx,))
        attrs.append(Segids._from_lookup(chain_idx, chainids))
        n_segments = len(chain_idx)

        return Topology(
            n_atoms,
            n_residues,
            n_segments,
            attrs=attrs,
            atom_resindex=residx,
            residue_segindex=segidx,
        )

    def _get_structure(self) -> "Structure":
        return _read_gemmi_structure(self.filename)
