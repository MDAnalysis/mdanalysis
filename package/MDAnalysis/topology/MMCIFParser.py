# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
"""
MMCIF Topology Parser
=====================

.. versionadded:: 2.11.0

Read topology information from mmCIF/PDBx coordinate files using the
`Gemmi library <https://gemmi.readthedocs.io>`_.



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

    format = ["cif", "cif.gz", "mmcif", "mmcif.gz"]

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
        model = structure[0]

        # TODO: gemmi.FlatStructure provides vectorised column access to all atom
        # fields and could replace the per-atom Python loop below for a speed-up
        # on large structures. gemmi API is still not 100% there yes so worth revisiting
        # later once it matures more
        altlocs = []
        serials = []
        names = []
        chainids = []
        elements = []
        formalcharges = []
        weights = []
        occupancies = []
        record_types = []
        tempfactors = []
        icodes = []
        resids = []
        resnames = []

        for chain in model:
            for residue in chain:
                match residue.het_flag:
                    case "A":
                        rec = "ATOM"
                    case "H":
                        rec = "HETATM"
                    case _:
                        raise ValueError(
                            "Found an atom that is neither ATOM nor HETATM"
                        )
                for atom in residue:
                    altlocs.append(atom.altloc if atom.has_altloc() else "")
                    serials.append(atom.serial)
                    names.append(atom.name)
                    chainids.append(chain.name)
                    elements.append(atom.element.name)
                    formalcharges.append(atom.charge)
                    weights.append(atom.element.weight)
                    occupancies.append(atom.occ)
                    record_types.append(rec)
                    tempfactors.append(atom.b_iso)
                    icodes.append(residue.seqid.icode.strip())
                    resids.append(residue.seqid.num)
                    resnames.append(residue.name)

        # Atom Attributes
        attrs = [
            AltLocs(altlocs),
            Atomids(serials),
            Atomnames(names),
            Atomtypes(names),
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
        n_atoms = len(altlocs)

        # Residue Attributes
        resids = np.array(resids)
        resnames = np.array(resnames)
        icodes = np.array(icodes)
        chainids = np.array(chainids)
        residx, (resids, resnames, icodes, chainids) = change_squash(
            (resids, resnames, icodes, chainids),
            (resids, resnames, icodes, chainids),
        )
        attrs.append(Resids(resids))
        attrs.append(Resnames(resnames))
        attrs.append(Resnums(resids.copy()))
        attrs.append(ICodes(icodes))
        n_residues = len(resids)

        # Segment Attributes
        segidx, (segids,) = change_squash((chainids,), (chainids,))
        attrs.append(Segids(segids))
        n_segments = len(segids)

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
