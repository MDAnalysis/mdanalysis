""" Various utilities used by ParmEd that don't really fit elsewhere """
from ..exceptions import MoleculeError as _MoleculeError
from .pairlist import find_atom_pairs
from ..topologyobjects import Atom
from shutil import which
from typing import List, Set, Tuple
import sys

__all__ = ['io', 'timer', 'which', 'tag_molecules', 'PYPY', 'find_atom_pairs']

PYPY = '__pypy__' in sys.builtin_module_names

def tag_molecules(struct) -> List[Set[int]]:
    """
    Sets the ``marked`` attribute of every Atom in struct to the molecule number
    it is a part of. If no bonds are present, every atom is its own molecule.

    Parameters
    ----------
    struct : :class:`parmed.Structure`
        Input structure to tag the molecules for
    """
    # Make sure our recursion limit is large enough, but never shrink it
    from sys import setrecursionlimit, getrecursionlimit
    setrecursionlimit(max(int(1.2 * len(struct.atoms)), getrecursionlimit()))

    if not struct.bonds:
        for i, atom in enumerate(struct.atoms):
            atom.marked = i + 1
        return [{i + 1} for i in range(len(struct.atoms))]
    owner = []
    # We do have bonds, this is the interesting part
    struct.atoms.unmark()
    mol_id = 1
    for atom in struct.atoms:
        if atom.marked:
            continue
        atom.marked = mol_id
        owner.append({atom.idx})
        _set_owner(atom, owner[-1], mol_id)
        mol_id += 1

    return owner

def _set_owner(atm: Atom, owner: Set[int], mol_id: int) -> None:
    """ Recursively sets ownership of given atom and all bonded partners """
    for partner in atm.bond_partners:
        if not partner.marked:
            owner.add(partner.idx)
            partner.marked = mol_id
            _set_owner(partner, owner, mol_id)
        elif partner.marked != mol_id:
            raise _MoleculeError(f'Atom {partner.idx} in multiple molecules')
