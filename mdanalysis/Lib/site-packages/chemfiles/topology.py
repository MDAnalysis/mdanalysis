from ctypes import c_bool, c_uint64
from enum import IntEnum

import numpy as np

from ._c_api import chfl_bond_order
from .atom import Atom
from .residue import Residue
from .utils import CxxPointer


class BondOrder(IntEnum):
    """
    Possible bond orders are:

    - ``BondOrder.Unknown``: the bond order is not specified
    - ``BondOrder.Single``: bond order for single bond
    - ``BondOrder.Double``: bond order for double bond
    - ``BondOrder.Triple``: bond order for triple bond
    - ``BondOrder.Quadruple``: bond order for quadruple bond (present in some metals)
    - ``BondOrder.Quintuplet``: bond order for quintuplet bond (present in some metals)
    - ``BondOrder.Amide``: bond order for amide bond
    - ``BondOrder.Aromatic``: bond order for aromatic bond
    """

    Unknown = chfl_bond_order.CHFL_BOND_UNKNOWN
    Single = chfl_bond_order.CHFL_BOND_SINGLE
    Double = chfl_bond_order.CHFL_BOND_DOUBLE
    Triple = chfl_bond_order.CHFL_BOND_TRIPLE
    Quadruple = chfl_bond_order.CHFL_BOND_QUADRUPLE
    Quintuplet = chfl_bond_order.CHFL_BOND_QUINTUPLET
    Amide = chfl_bond_order.CHFL_BOND_AMIDE
    Aromatic = chfl_bond_order.CHFL_BOND_AROMATIC


class TopologyAtoms(object):
    """Proxy object to get the atoms in a topology"""

    def __init__(self, topology):
        self.topology = topology

    def __len__(self):
        """Get the current number of atoms in this :py:class:`Topology`."""
        count = c_uint64()
        self.topology.ffi.chfl_topology_atoms_count(self.topology.ptr, count)
        return count.value

    def __getitem__(self, index):
        """
        Get a reference to the :py:class:`Atom` at the given ``index`` in the
        associated :py:class:`Topology`.
        """
        if index >= len(self):
            raise IndexError(f"atom index ({index}) out of range for this topology")
        else:
            ptr = self.topology.ffi.chfl_atom_from_topology(
                self.topology.mut_ptr, c_uint64(index)
            )
            return Atom.from_mutable_ptr(self, ptr)

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    def __delitem__(self, index):
        self.remove(index)

    def __repr__(self):
        return "[" + ", ".join([atom.__repr__() for atom in self]) + "]"

    def remove(self, index):
        """
        Remove the :py:class:`Atom` at the given ``index`` from the associated
        :py:class:`Topology`.

        This shifts all the atoms indexes larger than ``index`` by 1  (``n``
        becomes ``n - 1``);
        """
        self.topology.ffi.chfl_topology_remove(self.topology.mut_ptr, c_uint64(index))

    def append(self, atom):
        """
        Add a copy of the :py:class:`Atom` ``atom`` at the end of this
        :py:class:`Topology`.
        """
        self.topology.ffi.chfl_topology_add_atom(self.topology.mut_ptr, atom.ptr)


class TopologyResidues(object):
    """Proxy object to get the residues in a topology"""

    def __init__(self, topology):
        self.topology = topology

    def __len__(self):
        """Get the current number of residues in this :py:class:`Topology`."""
        count = c_uint64()
        self.topology.ffi.chfl_topology_residues_count(self.topology.ptr, count)
        return count.value

    def __getitem__(self, index):
        """
        Get read-only access to the :py:class:`Residue` at the given ``index``
        from the associated :py:class:`Topology`. The residue index in the
        topology does not necessarily match the residue id.
        """
        if index >= len(self):
            raise IndexError(f"residue index ({index}) out of range for this topology")
        else:
            ptr = self.topology.ffi.chfl_residue_from_topology(
                self.topology.ptr, c_uint64(index)
            )
            return Residue.from_const_ptr(self, ptr)

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    def __repr__(self):
        return "[" + ", ".join([residue.__repr__() for residue in self]) + "]"

    def append(self, residue):
        """
        Add the :py:class:`Residue` ``residue`` to this :py:class:`Topology`.

        The residue ``id`` must not already be in the topology, and the residue
        must contain only atoms that are not already in another residue.
        """
        self.topology.ffi.chfl_topology_add_residue(self.topology.mut_ptr, residue.ptr)


class Topology(CxxPointer):
    """
    A :py:class:`Topology` contains the definition of all the atoms in the
    system, and the liaisons between the atoms (bonds, angles, dihedrals, ...).
    It will also contain all the residues information if it is available.
    """

    def __init__(self):
        """Create a new empty :py:class:`Topology`."""
        super(Topology, self).__init__(self.ffi.chfl_topology(), is_const=False)

    def __copy__(self):
        return Topology.from_mutable_ptr(None, self.ffi.chfl_topology_copy(self.ptr))

    def __repr__(self):
        return f"Topology with {len(self.atoms)} atoms"

    @property
    def atoms(self):
        # late import to break circular dependency
        from .frame import FrameAtoms

        # if the topology comes from a frame, allow accessing atoms as
        # frame.topology.atoms anyway
        if self._CxxPointer__is_const:
            return FrameAtoms(self._CxxPointer__origin)
        else:
            return TopologyAtoms(self)

    def resize(self, count):
        """
        Resize this :py:class:`Topology` to contain ``count`` atoms. If the
        new number of atoms is bigger than the current number, new atoms will
        be created with an empty name and type. If it is lower than the current
        number of atoms, the last atoms will be removed, together with the
        associated bonds, angles and dihedrals.
        """
        self.ffi.chfl_topology_resize(self.mut_ptr, count)

    @property
    def residues(self):
        return TopologyResidues(self)

    def residue_for_atom(self, index):
        """
        Get read-only access to the :py:class:`Residue` containing the atom at
        the given ``index`` from this :py:class:`Topology`; or ``None`` if the
        atom is not part of a residue.
        """
        if index >= len(self.atoms):
            raise IndexError(f"residue index ({index}) out of range for this topology")

        ptr = self.ffi.chfl_residue_for_atom(self.ptr, c_uint64(index))
        if ptr:
            return Residue.from_const_ptr(self, ptr)
        else:
            return None

    def residues_linked(self, first, second):
        """
        Check if the two :py:class:`Residue` ``first`` and ``second`` from this
        :py:class:`Topology` are linked together, *i.e.* if there is a bond
        between one atom in the first residue and one atom in the second one.
        """
        linked = c_bool()
        self.ffi.chfl_topology_residues_linked(self.ptr, first.ptr, second.ptr, linked)
        return linked.value

    def bonds_count(self):
        """Get the number of bonds in this :py:class:`Topology`."""
        bonds = c_uint64()
        self.ffi.chfl_topology_bonds_count(self.ptr, bonds)
        return bonds.value

    def angles_count(self):
        """Get the number of angles in this :py:class:`Topology`."""
        angles = c_uint64()
        self.ffi.chfl_topology_angles_count(self.ptr, angles)
        return angles.value

    def dihedrals_count(self):
        """Get the number of dihedral angles in this :py:class:`Topology`."""
        dihedrals = c_uint64()
        self.ffi.chfl_topology_dihedrals_count(self.ptr, dihedrals)
        return dihedrals.value

    def impropers_count(self):
        """Get the number of improper angles in this :py:class:`Topology`."""
        impropers = c_uint64()
        self.ffi.chfl_topology_impropers_count(self.ptr, impropers)
        return impropers.value

    @property
    def bonds(self):
        """Get the list of bonds in this :py:class:`Topology`."""
        count = self.bonds_count()
        bonds = np.zeros((count, 2), np.uint64)
        self.ffi.chfl_topology_bonds(self.ptr, bonds, c_uint64(count))
        return bonds

    def bonds_order(self, i, j):
        """
        Get the bonds order corresponding to the bond between atoms i and j
        """
        order = chfl_bond_order()
        self.ffi.chfl_topology_bond_order(self.ptr, c_uint64(i), c_uint64(j), order)
        return BondOrder(order.value)

    @property
    def bonds_orders(self):
        """
        Get the list of bonds order for each bond in this :py:class:`Topology`.
        """
        count = self.bonds_count()
        orders = np.zeros(count, chfl_bond_order)
        self.ffi.chfl_topology_bond_orders(self.ptr, orders, c_uint64(count))
        return list(map(BondOrder, orders))

    @property
    def angles(self):
        """Get the list of angles in this :py:class:`Topology`."""
        count = self.angles_count()
        angles = np.zeros((count, 3), np.uint64)
        self.ffi.chfl_topology_angles(self.ptr, angles, c_uint64(count))
        return angles

    @property
    def dihedrals(self):
        """Get the list of dihedral angles in this :py:class:`Topology`."""
        count = self.dihedrals_count()
        dihedrals = np.zeros((count, 4), np.uint64)
        self.ffi.chfl_topology_dihedrals(self.ptr, dihedrals, c_uint64(count))
        return dihedrals

    @property
    def impropers(self):
        """Get the list of improper angles in this :py:class:`Topology`."""
        count = self.impropers_count()
        impropers = np.zeros((count, 4), np.uint64)
        self.ffi.chfl_topology_impropers(self.ptr, impropers, c_uint64(count))
        return impropers

    def add_bond(self, i, j, order=None):
        """
        Add a bond between the atoms at indexes ``i`` and ``j`` in this
        :py:class:`Topology`, optionally setting the bond ``order``.
        """
        if order is None:
            self.ffi.chfl_topology_add_bond(self.mut_ptr, c_uint64(i), c_uint64(j))
        else:
            self.ffi.chfl_topology_bond_with_order(
                self.mut_ptr, c_uint64(i), c_uint64(j), chfl_bond_order(order)
            )

    def remove_bond(self, i, j):
        """
        Remove any existing bond between the atoms at indexes ``i`` and ``j``
        in this :py:class:`Topology`.

        This function does nothing if there is no bond between ``i`` and ``j``.
        """
        self.ffi.chfl_topology_remove_bond(self.mut_ptr, c_uint64(i), c_uint64(j))

    def clear_bonds(self):
        """
        Remove all existing bonds, angles, dihedral angles and improper dihedral
        angles in this topology.
        """
        self.ffi.chfl_topology_clear_bonds(self.mut_ptr)
