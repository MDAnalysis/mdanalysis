from ctypes import c_bool, c_char_p, c_int64, c_uint64

import numpy as np

from .misc import ChemfilesError
from .property import Property
from .utils import CxxPointer, _call_with_growing_buffer


class ResidueAtoms(object):
    """Proxy object to get the atomic indexes in a residue"""

    def __init__(self, residue):
        self.residue = residue
        # Use a cache for indexes
        self.indexes = None

    def __len__(self):
        """Get the current number of atoms in this :py:class:`Residue`."""
        if self.indexes is not None:
            return len(self.indexes)
        else:
            count = c_uint64()
            self.residue.ffi.chfl_residue_atoms_count(self.residue.ptr, count)
            return count.value

    def __contains__(self, index):
        """
        Check if the :py:class:`Residue` contains the atom at index ``atom``.
        """
        if self.indexes is not None:
            return index in self.indexes
        else:
            result = c_bool()
            self.residue.ffi.chfl_residue_contains(
                self.residue.ptr, c_uint64(index), result
            )
            return result.value

    def __getitem__(self, i):
        """
        Get the atomic index number ``i`` in the associated :py:class:`Residue`.
        """
        if self.indexes is None:
            count = len(self)
            self.indexes = np.zeros(count, np.uint64)
            self.residue.ffi.chfl_residue_atoms(
                self.residue.ptr, self.indexes, c_uint64(count)
            )

        return self.indexes[i]

    def __iter__(self):
        for i in range(len(self)):
            yield self[i]

    def __repr__(self):
        return "[" + ", ".join([i.__repr__() for i in self]) + "]"

    def append(self, atom):
        """Add the atom index ``atom`` in the :py:class:`Residue`."""
        self.residue.ffi.chfl_residue_add_atom(self.residue.ptr, c_uint64(atom))
        # reset the cache for indexes
        self.indexes = None


class Residue(CxxPointer):
    """
    A :py:class:`Residue` is a group of atoms belonging to the same logical
    unit. They can be small molecules, amino-acids in a protein, monomers in
    polymers, etc.
    """

    def __init__(self, name, resid=None):
        """
        Create a new :py:class:`Residue` from a ``name`` and optionally a
        residue id ``resid``.
        """

        if resid:
            ptr = self.ffi.chfl_residue_with_id(name.encode("utf8"), c_int64(resid))
        else:
            ptr = self.ffi.chfl_residue(name.encode("utf8"))
        super(Residue, self).__init__(ptr, is_const=False)

    def __copy__(self):
        return Residue.from_mutable_ptr(None, self.ffi.chfl_residue_copy(self.ptr))

    def __repr__(self):
        return f"Residue('{self.name}') with {len(self.atoms)} atoms"

    @property
    def name(self):
        """Get the name of this :py:class:`Residue`."""
        return _call_with_growing_buffer(
            lambda buffer, size: self.ffi.chfl_residue_name(self.ptr, buffer, size),
            initial=32,
        )

    @property
    def id(self):
        """Get the :py:class:`Residue` index in the initial topology."""
        id = c_int64()
        self.ffi.chfl_residue_id(self.ptr, id)
        return id.value

    @property
    def atoms(self):
        return ResidueAtoms(self)

    def __iter__(self):
        # Disable automatic iteration from __getitem__
        raise TypeError("can not iterate over a residue")

    def __getitem__(self, name):
        """
        Get a property of this residude with the given ``name``, or raise an
        error if the property does not exists.
        """
        if not isinstance(name, str):
            raise ChemfilesError(
                f"Invalid type {type(name)} for a residue property name"
            )

        ptr = self.ffi.chfl_residue_get_property(self.ptr, name.encode("utf8"))
        return Property.from_mutable_ptr(self, ptr).get()

    def __setitem__(self, name, value):
        """
        Set a property of this residue, with the given ``name`` and ``value``.
        The new value overwrite any pre-existing property with the same name.
        """
        if not isinstance(name, str):
            raise ChemfilesError(
                f"Invalid type {type(name)} for a residue property name"
            )

        property = Property(value)
        self.ffi.chfl_residue_set_property(
            self.mut_ptr, name.encode("utf8"), property.ptr
        )

    def properties_count(self):
        """Get the number of properties in this residue."""
        count = c_uint64()
        self.ffi.chfl_residue_properties_count(self.ptr, count)
        return count.value

    def list_properties(self):
        """Get the name of all properties in this residue."""
        count = self.properties_count()
        StringArray = c_char_p * count
        names = StringArray()
        self.ffi.chfl_residue_list_properties(self.ptr, names, count)
        return list(map(lambda n: n.decode("utf8"), names))
