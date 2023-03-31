from ctypes import POINTER, c_bool, c_char_p, c_double, c_uint64

import numpy as np

from ._c_api import chfl_bond_order, chfl_vector3d
from .atom import Atom
from .cell import UnitCell
from .misc import ChemfilesError
from .property import Property
from .topology import Topology
from .utils import CxxPointer


class FrameAtoms(object):
    """Proxy object to get the atoms in a frame"""

    def __init__(self, frame):
        self.frame = frame

    def __len__(self):
        """Get the current number of atoms in this :py:class:`Frame`."""
        count = c_uint64()
        self.frame.ffi.chfl_frame_atoms_count(self.frame.ptr, count)
        return count.value

    def __getitem__(self, index):
        """
        Get a reference to the :py:class:`Atom` at the given ``index`` in the
        associated :py:class:`Frame`.
        """
        if index >= len(self):
            raise IndexError(f"atom index ({index}) out of range for this frame")
        else:
            ptr = self.frame.ffi.chfl_atom_from_frame(
                self.frame.mut_ptr, c_uint64(index)
            )
            return Atom.from_mutable_ptr(self, ptr)

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    def __repr__(self):
        return "[" + ", ".join([atom.__repr__() for atom in self]) + "]"


class Frame(CxxPointer):
    """
    A :py:class:`Frame` contains data from one simulation step: the current
    unit cell, the topology, the positions, and the velocities of the particles
    in the system. If some information is missing (topology or velocity or unit
    cell), the corresponding data is filled with a default value.
    """

    def __init__(self):
        """
        Create an empty :py:class:`Frame` that will be resized by the runtime
        as needed.
        """
        super(Frame, self).__init__(self.ffi.chfl_frame(), is_const=False)

    def __copy__(self):
        return Frame.from_mutable_ptr(None, self.ffi.chfl_frame_copy(self.ptr))

    def __repr__(self):
        return f"Frame with {len(self.atoms)} atoms"

    @property
    def atoms(self):
        return FrameAtoms(self)

    def resize(self, count):
        """
        Resize the positions, velocities and topology in this
        :py:class:`Frame`, to have space for `count` atoms.

        This function may invalidate any array of the positions or the
        velocities if the new size is bigger than the old one. In all cases,
        previous data is conserved. This function conserve the presence or
        absence of velocities.
        """
        self.ffi.chfl_frame_resize(self.mut_ptr, c_uint64(count))

    def add_atom(self, atom, position, velocity=None):
        """
        Add a copy of the :py:class:`Atom` ``atom`` and the corresponding
        ``position`` and ``velocity`` to this :py:class:`Frame`.

        ``velocity`` can be ``None`` if no velocity is associated with the
        atom.
        """
        position = chfl_vector3d(position[0], position[1], position[2])
        if velocity:
            velocity = chfl_vector3d(velocity[0], velocity[1], velocity[2])
        self.ffi.chfl_frame_add_atom(self.mut_ptr, atom.ptr, position, velocity)

    def remove(self, index):
        """
        Remove the atom at the given ``index`` in this :py:class:`Frame`.

        This shifts all the atoms indexes larger than ``index`` by 1  (``n``
        becomes ``n - 1``); and invalidate any array obtained using
        :py:func:`Frame.positions` or :py:func:`Frame.velocities`.
        """
        self.ffi.chfl_frame_remove(self.mut_ptr, c_uint64(index))

    def add_bond(self, i, j, order=None):
        """
        Add a bond between the atoms at indexes ``i`` and ``j`` in this
        :py:class:`Frame`'s topology, optionally setting the bond ``order``.
        """
        if order is None:
            self.ffi.chfl_frame_add_bond(self.mut_ptr, c_uint64(i), c_uint64(j))
        else:
            self.ffi.chfl_frame_bond_with_order(
                self.mut_ptr, c_uint64(i), c_uint64(j), chfl_bond_order(order)
            )

    def remove_bond(self, i, j):
        """
        Remove any existing bond between the atoms at indexes ``i`` and ``j``
        in this :py:class:`Frame`'s topology.

        This function does nothing if there is no bond between ``i`` and ``j``.
        """
        self.ffi.chfl_frame_remove_bond(self.mut_ptr, c_uint64(i), c_uint64(j))

    def clear_bonds(self):
        """
        Remove all existing bonds, angles, dihedral angles and improper dihedral
        angles in this frame.
        """
        self.ffi.chfl_frame_clear_bonds(self.mut_ptr)

    def add_residue(self, residue):
        """
        Add the :py:class:`Residue` ``residue`` to this :py:class:`Frame`'s
        topology.

        The residue ``id`` must not already be in the topology, and the residue
        must contain only atoms that are not already in another residue.
        """
        self.ffi.chfl_frame_add_residue(self.mut_ptr, residue.ptr)

    @property
    def positions(self):
        """
        Get a view into the positions of this :py:class:`Frame`.

        This function gives direct access to the positions as a numpy array.
        Modifying the array will change the positions in the frame.

        If the frame is resized (by writing to it, calling
        :py:func:`Frame.resize`, :py:func:`Frame.add_atom`,
        :py:func:`Frame.remove`), the array is invalidated. Accessing it can
        cause a segfault.
        """
        count = c_uint64()
        data = POINTER(chfl_vector3d)()
        self.ffi.chfl_frame_positions(self.mut_ptr, data, count)
        count = count.value
        if count != 0:
            positions = np.ctypeslib.as_array(data, shape=(count,))
            return positions.view(np.float64).reshape((count, 3))
        else:
            return np.array([[], [], []], dtype=np.float64)

    @property
    def velocities(self):
        """
        Get a view into the velocities of this :py:class:`Frame`.

        This function gives direct access to the velocities as a numpy array.
        Modifying the array will change the velocities in the frame.

        If the frame is resized (by writing to it, calling
        :py:func:`Frame.resize`, :py:func:`Frame.add_atom`,
        :py:func:`Frame.remove`), the array is invalidated. Accessing it can
        cause a segfault.
        """
        count = c_uint64()
        data = POINTER(chfl_vector3d)()
        self.ffi.chfl_frame_velocities(self.mut_ptr, data, count)
        count = count.value
        if count != 0:
            velocities = np.ctypeslib.as_array(data, shape=(count,))
            return velocities.view(np.float64).reshape((count, 3))
        else:
            return np.array([[], [], []], dtype=np.float64)

    def add_velocities(self):
        """
        Add velocity data to this :py:class:`Frame`.

        The velocities are initialized to zero. If the frame already contains
        velocities, this function does nothing.
        """
        self.ffi.chfl_frame_add_velocities(self.mut_ptr)

    def has_velocities(self):
        """Check if this :py:class:`Frame` contains velocity."""
        velocities = c_bool()
        self.ffi.chfl_frame_has_velocities(self.ptr, velocities)
        return velocities.value

    @property
    def cell(self):
        """
        Get a mutable reference to the :py:class:`UnitCell` of this
        :py:class:`Frame`. Any modification to the cell will be reflected in
        the frame.
        """
        return UnitCell.from_mutable_ptr(
            self, self.ffi.chfl_cell_from_frame(self.mut_ptr)
        )

    @cell.setter
    def cell(self, cell):
        """
        Set the :py:class:`UnitCell` of this :py:class:`Frame` to ``cell``.
        """
        self.ffi.chfl_frame_set_cell(self.mut_ptr, cell.ptr)

    @property
    def topology(self):
        """
        Get read-only access to the :py:class:`Topology` of this
        :py:class:`Frame`.
        """
        return Topology.from_const_ptr(
            self, self.ffi.chfl_topology_from_frame(self.ptr)
        )

    @topology.setter
    def topology(self, topology):
        """
        Set the :py:class:`Topology` of this :py:class:`Frame` to ``topology``.
        """
        self.ffi.chfl_frame_set_topology(self.mut_ptr, topology.ptr)

    @property
    def step(self):
        """
        Get the step of this :py:class:`Frame`, i.e. the frame number in the
        trajectory.
        """
        step = c_uint64()
        self.ffi.chfl_frame_step(self.ptr, step)
        return step.value

    @step.setter
    def step(self, value):
        """Set the step for this :py:class:`Frame` to the given ``value``."""
        self.ffi.chfl_frame_set_step(self.mut_ptr, c_uint64(value))

    def guess_bonds(self):
        """
        Guess the bonds, angles and dihedrals in this :py:class:`Frame`.

        The bonds are guessed using a distance-based algorithm, and then angles
        and dihedrals are guessed from the bonds.
        """
        self.ffi.chfl_frame_guess_bonds(self.mut_ptr)

    def distance(self, i, j):
        """
        Get the distance (in Ångströms) between the atoms at indexes ``i`` and
        ``j`` in this :py:class:`Frame`, taking periodic boundary conditions
        into account.
        """
        distance = c_double()
        self.ffi.chfl_frame_distance(self.ptr, c_uint64(i), c_uint64(j), distance)
        return distance.value

    def angle(self, i, j, k):
        """
        Get the angle (in radians) formed by the atoms at indexes ``i``, ``j``
        and ``k`` in this :py:class:`Frame`, taking periodic boundary conditions
        into account.
        """
        angle = c_double()
        self.ffi.chfl_frame_angle(
            self.ptr, c_uint64(i), c_uint64(j), c_uint64(k), angle
        )
        return angle.value

    def dihedral(self, i, j, k, m):
        """
        Get the dihedral angle (in radians) formed by the atoms at indexes
        ``i``, ``j``, ``k`` and ``m`` in this :py:class:`Frame`, taking periodic
        boundary conditions into account.
        """
        dihedral = c_double()
        self.ffi.chfl_frame_dihedral(
            self.ptr, c_uint64(i), c_uint64(j), c_uint64(k), c_uint64(m), dihedral
        )
        return dihedral.value

    def out_of_plane(self, i, j, k, m):
        """
        Get the out of plane distance (in Ångströms) formed by the atoms at
        indexes ``i``, ``j``, ``k`` and ``m`` in this :py:class:`Frame`, taking
        periodic boundary conditions into account.

        This is the distance betweent the atom j and the ikm plane. The j atom
        is the center of the improper dihedral angle formed by i, j, k and m.
        """
        distance = c_double()
        self.ffi.chfl_frame_out_of_plane(
            self.ptr, c_uint64(i), c_uint64(j), c_uint64(k), c_uint64(m), distance
        )
        return distance.value

    def __iter__(self):
        # Disable automatic iteration from __getitem__
        raise TypeError("use Frame.atoms to iterate over a frame")

    def __getitem__(self, name):
        """
        Get a property of this frame with the given ``name``, or raise an error
        if the property does not exists.
        """
        if not isinstance(name, str):
            raise ChemfilesError(f"Invalid type {type(name)} for a frame property name")

        ptr = self.ffi.chfl_frame_get_property(self.ptr, name.encode("utf8"))
        return Property.from_mutable_ptr(self, ptr).get()

    def __setitem__(self, name, value):
        """
        Set a property of this frame, with the given ``name`` and ``value``.
        The new value overwrite any pre-existing property with the same name.
        """
        if not isinstance(name, str):
            raise ChemfilesError(f"Invalid type {type(name)} for a frame property name")

        property = Property(value)
        self.ffi.chfl_frame_set_property(
            self.mut_ptr, name.encode("utf8"), property.ptr
        )

    def properties_count(self):
        """Get the number of properties in this frame."""
        count = c_uint64()
        self.ffi.chfl_frame_properties_count(self.ptr, count)
        return count.value

    def list_properties(self):
        """Get the name of all properties in this frame."""
        count = self.properties_count()
        StringArray = c_char_p * count
        names = StringArray()
        self.ffi.chfl_frame_list_properties(self.ptr, names, count)
        return list(map(lambda n: n.decode("utf8"), names))
