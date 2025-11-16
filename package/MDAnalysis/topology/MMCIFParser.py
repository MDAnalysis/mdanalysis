# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
"""
MMCIF Topology Parser
====================

Read topology information from mmCIF/PDBx coordinate files using the `Gemmi library <https://github.com/project-gemmi/gemmi>`_

mmCIF files contain topology information about the molecules in the structure. For each atom the following attributes are read
and stored in the relevant topology attributes:
    - :class:`..core.topologyattrs.AtomAttr` subclasses:
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
    - :class:`..core.topologyattrs.ResidueAttr` subclasses:
        - :class:`MDAnalysis.core.topologyattrs.Resnums`
        - :class:`MDAnalysis.core.topologyattrs.ICodes`
        - :class:`MDAnalysis.core.topologyattrs.Resids`
        - :class:`MDAnalysis.core.topologyattrs.Resnames`
    - :class:`..core.topologyattrs.SegmentAttr` subclasses:
        - :class:`MDAnalysis.core.topologyattrs.Segids`

Classes
-------

.. autoclass:: MMCIFParser
   :members:
   :inherited-members:

.. versionadded:: 2.9.0
"""

try:
    import gemmi
except ImportError:
    HAS_GEMMI = False
else:
    HAS_GEMMI = True

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
from ..lib import util


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

    .. versionadded:: 2.9.0
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
                f"MMCIF model {self.filename} contains {len(structure)=} different models, "
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
        # This method exists because of some lacking methods in the gemmi Python API.
        # within gemmi in C++, one can call `read_structure` and in-memory, string, and filepath
        # arguments will all be accepted:
        # https://github.com/project-gemmi/gemmi/blob/4416e298f204b7b57bf5b3051d7efd4fe02957cf/include/gemmi/mmread.hpp#L86

        # However, for MDA to similarly accept common input types like streams (open File-like objs and StringIO objs)
        # as well as pathlib.Path() objects, we have to use the Python API methods available currently (as of 0.7.3)
        # with a string as a common target for all input types
        # For this, we call gemmi.cif.read_string (https://gemmi.readthedocs.io/en/latest/cif.html#reading) to handle CIF
        # strings and gemmi.read_pdb to handle PDB strings (no one method can handle both formats currently Py-side)

        # openany() is called instead of passing file paths (when available) differently from streams
        # even though reading the file into a string is less efficient, this is easier to maintain

        # if the gemmi Python API is extended, this method can be simplified/removed and replaced with something like
        # gemmi.read_structure

        with util.openany(self.filename) as f:
            content_as_str = f.read()
            try:
                # String -> Doc -> Block -> Structure
                # making Structure from first Block in Document as is done internally in gemmi:
                # https://github.com/project-gemmi/gemmi/blob/4416e298f204b7b57bf5b3051d7efd4fe02957cf/include/gemmi/mmcif.hpp#L32
                return gemmi.make_structure_from_block(
                    gemmi.cif.read_string(content_as_str)[0]
                )
            except ValueError:
                return gemmi.read_pdb_string(content_as_str)
