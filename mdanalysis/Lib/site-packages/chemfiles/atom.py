from ctypes import c_char_p, c_double, c_uint64

from .misc import ChemfilesError
from .property import Property
from .utils import CxxPointer, _call_with_growing_buffer


class Atom(CxxPointer):
    """
    An :py:class:`Atom` is a particle in the current :py:class:`Frame`. It
    stores the following atomic properties:

    - atom name;
    - atom type;
    - atom mass;
    - atom charge.

    The atom name is usually an unique identifier (``"H1"``, ``"C_a"``) while
    the atom type will be shared between all particles of the same type:
    ``"H"``, ``"Ow"``, ``"CH3"``.
    """

    def __init__(self, name, type=None):
        """
        Create a new :py:class:`Atom` with the given ``name``. If ``type`` is
        present, use it as the atom type. Else the atom name is used as atom
        type.
        """
        ptr = self.ffi.chfl_atom(name.encode("utf8"))
        super(Atom, self).__init__(ptr, is_const=False)
        if type:
            self.type = type

    def __copy__(self):
        return Atom.from_mutable_ptr(None, self.ffi.chfl_atom_copy(self.ptr))

    def __repr__(self):
        if self.type == self.name:
            return f"Atom('{self.name}')"
        else:
            return f"Atom('{self.name}', '{self.type}')"

    @property
    def mass(self):
        """Get this :py:class:`Atom` mass, in atomic mass units."""
        mass = c_double()
        self.ffi.chfl_atom_mass(self.ptr, mass)
        return mass.value

    @mass.setter
    def mass(self, mass):
        """Set this :py:class:`Atom` mass, in atomic mass units."""
        self.ffi.chfl_atom_set_mass(self.mut_ptr, c_double(mass))

    @property
    def charge(self):
        """
        Get this :py:class:`Atom` charge, in number of the electron charge *e*.
        """
        charge = c_double()
        self.ffi.chfl_atom_charge(self.ptr, charge)
        return charge.value

    @charge.setter
    def charge(self, charge):
        """
        Set this :py:class:`Atom` charge, in number of the electron charge *e*.
        """
        self.ffi.chfl_atom_set_charge(self.mut_ptr, c_double(charge))

    @property
    def name(self):
        """Get this :py:class:`Atom` name."""
        return _call_with_growing_buffer(
            lambda buffer, size: self.ffi.chfl_atom_name(self.ptr, buffer, size),
            initial=32,
        )

    @name.setter
    def name(self, name):
        """Set this :py:class:`Atom` name to ``name``."""
        self.ffi.chfl_atom_set_name(self.mut_ptr, name.encode("utf8"))

    @property
    def type(self):
        """Get this :py:class:`Atom` type."""
        return _call_with_growing_buffer(
            lambda buffer, size: self.ffi.chfl_atom_type(self.ptr, buffer, size),
            initial=32,
        )

    @type.setter
    def type(self, type):
        """Set this :py:class:`Atom` type to ``type``."""
        self.ffi.chfl_atom_set_type(self.mut_ptr, type.encode("utf8"))

    @property
    def full_name(self):
        """
        Try to get the full name of this :py:class:`Atom` from its type. For
        example, the full name of "He" is "Helium". If the name can not be
        found, returns the empty string.
        """
        return _call_with_growing_buffer(
            lambda buffer, size: self.ffi.chfl_atom_full_name(self.ptr, buffer, size),
            initial=64,
        )

    @property
    def vdw_radius(self):
        """
        Try to get the Van der Waals radius of this :py:class:`Atom` from its
        type. If the radius can not be found, returns 0.
        """
        radius = c_double()
        self.ffi.chfl_atom_vdw_radius(self.ptr, radius)
        return radius.value

    @property
    def covalent_radius(self):
        """
        Try to get the covalent radius of this :py:class:`Atom` from its type.
        If the radius can not be found, returns 0.
        """
        radius = c_double()
        self.ffi.chfl_atom_covalent_radius(self.ptr, radius)
        return radius.value

    @property
    def atomic_number(self):
        """
        Try to get the atomic number of this :py:class:`Atom` from its type. If
        the atomic number can not be found, returns 0.
        """
        number = c_uint64()
        self.ffi.chfl_atom_atomic_number(self.ptr, number)
        return number.value

    def __iter__(self):
        # Disable automatic iteration from __getitem__
        raise TypeError("can not iterate over an atom")

    def __getitem__(self, name):
        """
        Get a property of this atom with the given ``name``, or raise an error
        if the property does not exists.
        """
        if not isinstance(name, str):
            raise ChemfilesError(
                f"Invalid type {type(name)} for an atomic property name"
            )

        ptr = self.ffi.chfl_atom_get_property(self.ptr, name.encode("utf8"))
        return Property.from_mutable_ptr(self, ptr).get()

    def __setitem__(self, name, value):
        """
        Set a property of this atom, with the given ``name`` and ``value``.
        The new value overwrite any pre-existing property with the same name.
        """
        if not isinstance(name, str):
            raise ChemfilesError(f"invalid type {type(name)} for a property name")

        property = Property(value)
        self.ffi.chfl_atom_set_property(self.mut_ptr, name.encode("utf8"), property.ptr)

    def properties_count(self):
        """Get the number of properties in this atom."""
        count = c_uint64()
        self.ffi.chfl_atom_properties_count(self.ptr, count)
        return count.value

    def list_properties(self):
        """Get the name of all properties in this atom."""
        count = self.properties_count()
        StringArray = c_char_p * count
        names = StringArray()
        self.ffi.chfl_atom_list_properties(self.ptr, names, count)
        return list(map(lambda n: n.decode("utf8"), names))
