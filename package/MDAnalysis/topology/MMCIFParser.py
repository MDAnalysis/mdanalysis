# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
"""
MMCIF Topology Parser
====================

Read topology information from mmCIF/PDBx coordinate files using the `Gemmi library <https://github.com/project-gemmi/gemmi>`_

mmCIF files contain topology information about the molecules in the structure. For each atom the following attributes are read
and stored in the relevant topology attributes:
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

.. versionadded:: 2.11.0
"""

import warnings

import numpy as np

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
from ..coordinates.MMCIF import _read_gemmi_structure


class MMCIFParser(TopologyReaderBase):
    """Parser that obtains a list of atoms from a standard MMCIF/PDBx file using the `gemmi library <https://github.com/project-gemmi/gemmi>`_.

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

    .. versionadded:: 2.11.0
    """

    format = ["cif", "cif.gz", "mmcif", "mmcif.gz"]

    def parse(self, **kwargs) -> Topology:
        """Read the file and return the structure.

        Returns
        -------
        MDAnalysis Topology object
        """
        structure = self._get_structure()

        if len(structure) > 1:
            warnings.warn(
                f"MMCIF model {self.filename} contains {len(structure)} different models, "
                "but only the first one will be used to assign the topology"
            )
        model = structure[0]

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
            # ------------------
            icodes,  # residue.seqid.icode
            resids,  # residue.seqid.num
            resnames,  # residue.name
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
                            # ------------------
                            residue.seqid.icode,  # icodes
                            residue.seqid.num,  # resids
                            residue.name,  # resnames
                        )
                        # the main loop over the `gemmi.Model` object
                        for chain in model
                        for residue in chain
                        for atom in residue
                    ]
                )
            ),
        )

        # fill in altlocs, since gemmi has '' as default
        altlocs = ["A" if not elem else elem for elem in altlocs]

        # convert default gemmi record types to default MDAnalysis record types
        record_types = [
            "ATOM" if record == "A" else "HETATM" if record == "H" else None
            for record in record_types
        ]
        if any((elem is None for elem in record_types)):
            raise ValueError("Found an atom that is neither ATOM or HETATM")

        # Atom Attr's
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
        n_atoms = len(altlocs)

        # Residue Attr's
        residx, (resids, resnames, icodes, chainids) = change_squash(
            (resids, resnames, icodes, chainids),
            (resids, resnames, icodes, chainids),
        )
        attrs.append(Resids(resids))
        attrs.append(Resnames(resnames))
        attrs.append(Resnums(resids.copy()))
        attrs.append(ICodes([icode.strip() for icode in icodes]))
        n_residues = len(resids)

        # Segment Attr's
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

    def _get_structure(self):
        return _read_gemmi_structure(self.filename)
