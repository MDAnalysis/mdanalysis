""" Standard residue templates for biomolecular residues """

import json
from collections import Counter
from io import TextIOWrapper
try:
    from functools import cache
except ImportError:
    # Patch to support a `cache` decorator in Python 3.8 equivalent to functools.cache in Python 3.9+
    from functools import lru_cache
    cache = lru_cache(maxsize=None)
from typing import Mapping, Dict, Any, List, Set, TextIO
try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal
from pathlib import Path
import os
from ..amber.offlib import AmberOFFLibrary
from .residue import ResidueTemplate
from ..topologyobjects import Atom, QualitativeBondType, Hybridization
from ..periodic_table import AtomicNum

__all__ = ['StandardBiomolecularResidues']

StandardBiomolecularResidues = AmberOFFLibrary.parse(
    os.path.join(os.path.split(__file__)[0], 'data', 'standard_residues.lib')
)

def load_ccd_residue_templates(fileobj: TextIO) -> Mapping[str, ResidueTemplate]:
    """Loads a biopolymer residue database from a dictionary prepared from the chemical component dictionary

    Parameters
    ----------
    fileobj: file-like
        The file object to read the CCD residue templates from

    Returns
    -------
    The dictionary mapping residue template name with the ResidueTemplate object
    """
    residues = dict()
    all_residue_data = json.load(fileobj)

    for single_residue in all_residue_data:
        residues[single_residue["name"]] = _process_single_residue(single_residue)

    # Augment common residues from the Amber-standard residues
    extra_residues = set(StandardBiomolecularResidues.keys()) - set(residues.keys())
    residues.update({name: StandardBiomolecularResidues[name] for name in extra_residues})
    return residues


def _process_single_residue(data: Dict[str, Any]) -> ResidueTemplate:
    """Processes data from a single residue and builds a template from it."""
    # In the CCD, a "leaving" atom is one that is "removed" when a residue forms a polymeric chain
    # with an adjacent residue. We adopt the convention here that the 'head' and 'tail' of the polymeric
    # residues are the first and last atoms attached to those that are "leaving"
    #
    # It is harder to determine whether a single linkage point is a head or tail. That is done by
    # <insert heuristic algorithm here>
    template = ResidueTemplate(data["name"])
    _add_atoms(data["atoms"], template)
    _add_bonds(data["bonds"], template)
    _assign_head_tail(template, {atom["name"] for atom in data["atoms"] if atom["leaving"]}, data["type"])
    return template

def _add_atoms(atoms: List[Dict[str, Any]], template: ResidueTemplate) -> None:
    """Adds atoms to a residue template"""
    for atom in atoms:
        name = atom["name"]
        symbol = atom["symbol"]
        formal_charge = atom["formal_charge"]
        aromatic = atom["aromatic"]

        hybridization_value = atom.get("hybridization", None)
        hybridization = Hybridization(hybridization_value) if hybridization_value is not None else None

        atom = Atom(
            atomic_number=AtomicNum[symbol.title()], type=symbol, name=name, formal_charge=formal_charge,
            hybridization=hybridization, aromatic=aromatic
        )
        template.add_atom(atom)

def _add_bonds(bonds: List[Dict[str, Any]], template: ResidueTemplate) -> None:
    """Adds bonds to a template"""
    atom_map = {atom.name: atom for atom in template.atoms}
    _order_map = {
        ("SING", False): 1.0, ("SING", True): 1.5, ("DOUB", False): 2.0, ("DOUB", True): 1.5, ("TRIP", False): 3.0,
        ("TRIP", True): 2.5
    }
    _qualitative_map = {
        ("SING", False): QualitativeBondType.SINGLE, ("DOUB", False): QualitativeBondType.DOUBLE,
        ("SING", True): QualitativeBondType.AROMATIC, ("DOUB", True): QualitativeBondType.AROMATIC,
        ("TRIP", False): QualitativeBondType.TRIPLE, ("TRIP", True): QualitativeBondType.AROMATIC,
    }
    for bond in bonds:
        key = (bond["order"], bond["aromatic"])
        template.add_bond(
            atom_map[bond["name_a"]], atom_map[bond["name_b"]],
            order=_order_map[key], qualitative_type=_qualitative_map[key]
        )


def _head_or_tail(atom: Atom, residue_type: str) -> Literal["head", "tail", "neither", "maybe head", "maybe tail"]:
    if 'peptide' in residue_type.lower():
        if atom.atomic_number == AtomicNum["N"]:
            return "head" if atom.name == "N" else "maybe head"
        if atom.atomic_number == AtomicNum["C"]:
            return "tail" if atom.name == "C" else "maybe tail"
        return "neither"

    if 'rna' in residue_type.lower() or 'dna' in residue_type.lower():
        if atom.atomic_number == AtomicNum["P"]:
            return "head" if atom.name == "P" else "maybe head"
        elif atom.atomic_number == AtomicNum["O"]:
            return "tail" if atom.name == "O3'" else "maybe tail"
        return "neither"

    return "neither"

def _assign_head_tail(template: ResidueTemplate, removable_atoms: Set[str], residue_type: str) -> None:
    """Assigns the head and tail of a residue based on the atoms that can be removed"""
    atoms_joined_to_removable = set()
    for atom in template.atoms:
        if atom.name not in removable_atoms:
            continue
        for partner in atom.bond_partners:
            if partner.name in removable_atoms:
                continue
            atoms_joined_to_removable.add(partner)

    sorted_connectivity_atoms = sorted(atoms_joined_to_removable)
    head_tail_assignments = [_head_or_tail(atom, residue_type) for atom in sorted_connectivity_atoms]
    head_tail_counts = Counter(head_tail_assignments)
    assert head_tail_counts.get("head", 0) <= 1, "Too many heads. Check residue template file"
    assert head_tail_counts.get("tail", 0) <= 1, "Too many tails. Check residue template file"

    # Assign the head. In a spot-check of the residues where there was no atom that should clearly be assigned as
    # the head, the "first" atom in the sequence that was tagged as "maybe head" was always the correct choice.
    # Use that strategy
    if head_tail_counts.get("head", 0) == 1:
        for atom, label in zip(sorted_connectivity_atoms, head_tail_assignments):
            if label == "head":
                template.head = atom
                break
    else:
        for atom, label in zip(sorted_connectivity_atoms, head_tail_assignments):
            if label == "maybe head":
                # Stop at the first
                template.head = atom
                break

    # Assign the tail. In a spot-check of the residues where there was no atom that should clearly be assigned as
    # the head, the "last" atom in the sequence that was tagged as "maybe tail" was always the correct choice.
    # Use that strategy
    if head_tail_counts.get("tail", 0) == 1:
        for atom, label in zip(sorted_connectivity_atoms, head_tail_assignments):
            if label == "tail":
                template.tail = atom
                break
    else:
        for atom, label in zip(sorted_connectivity_atoms, head_tail_assignments):
            if label == "maybe tail":
                template.tail = atom
                # Do not stop at the first in case there's another one later.
                # Reversing would be faster, but optimization is not needed

    # Assign the remaining connectivity atoms that are neither head nor tail
    for atom in sorted_connectivity_atoms:
        if atom is not template.head and atom is not template.tail:
            template.connections.append(atom)


_CACHED_RESIDUE_TEMPLATE_LIBRARY = None

@cache
def get_standard_residue_template_library() -> Dict[str, ResidueTemplate]:
    path = Path(__file__).parent / "data" / "ccd_residue_templates.json"
    with path.open("r") as ifh:
        return load_ccd_residue_templates(ifh)


@cache
def get_nonstandard_ccd_residues() -> Dict[str, ResidueTemplate]:
    import gzip
    path = Path(__file__).parent / "data" / "nonstandard_ccd_residue_templates.json.gz"
    with TextIOWrapper(gzip.open(path, "rb")) as ifh:
        return load_ccd_residue_templates(ifh)
