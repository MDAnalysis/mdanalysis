from ctypes import c_bool, c_double

import numpy as np

from ._c_api import chfl_property_kind, chfl_vector3d
from .misc import ChemfilesError
from .utils import CxxPointer, _call_with_growing_buffer


class Property(CxxPointer):
    """
    A :py:class:`Property` holds the data used in properties in
    :py:class:`Frame` and :py:class:`Atom`. A property can have various types:
    bool, double, string or 3D vectors.

    This class is not meant for direct use, but is an internal class.
    """

    def __init__(self, value):
        """Create a new property containing the given value"""
        if isinstance(value, bool):
            ptr = self.ffi.chfl_property_bool(c_bool(value))
        elif isinstance(value, (float, int)):
            ptr = self.ffi.chfl_property_double(c_double(value))
        elif isinstance(value, str):
            ptr = self.ffi.chfl_property_string(value.encode("utf8"))
        elif _is_vector3d(value):
            value = chfl_vector3d(value[0], value[1], value[2])
            ptr = self.ffi.chfl_property_vector3d(value)
        else:
            raise ChemfilesError(
                f"can not create a Property with a value of type '{type(value)}'"
            )

        super(Property, self).__init__(ptr, is_const=False)

    def get(self):
        kind = chfl_property_kind()
        self.ffi.chfl_property_get_kind(self.ptr, kind)
        if kind.value == chfl_property_kind.CHFL_PROPERTY_BOOL:
            value = c_bool()
            self.ffi.chfl_property_get_bool(self.ptr, value)
            return value.value
        if kind.value == chfl_property_kind.CHFL_PROPERTY_DOUBLE:
            value = c_double()
            self.ffi.chfl_property_get_double(self.ptr, value)
            return value.value
        if kind.value == chfl_property_kind.CHFL_PROPERTY_STRING:

            def callback(buffer, size):
                self.ffi.chfl_property_get_string(self.ptr, buffer, size)

            return _call_with_growing_buffer(callback, initial=32)
        if kind.value == chfl_property_kind.CHFL_PROPERTY_VECTOR3D:
            value = chfl_vector3d()
            self.ffi.chfl_property_get_vector3d(self.ptr, value)
            return value[0], value[1], value[2]
        else:
            raise ChemfilesError("unknown property kind, this is a bug")


def _is_vector3d(value):
    try:
        a = np.array(value, dtype="double")
        return len(a) >= 3
    except Exception:
        return False
