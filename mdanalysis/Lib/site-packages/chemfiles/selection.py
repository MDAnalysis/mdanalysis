from ctypes import c_uint64

import numpy as np

from ._c_api import chfl_match
from .utils import CxxPointer, _call_with_growing_buffer


class Selection(CxxPointer):
    """
    Select atoms in a :py:class:`Frame` with a selection language.

    The selection language is built by combining basic operations. Each basic
    operation follows the ``<selector>[(<variable>)] <operator> <value>``
    structure, where ``<operator>`` is a comparison operator in
    ``== != < <= > >=``. Refer to the `full documentation
    <selections-doc>`_ to know the allowed selectors and how to use them.

    .. selections-doc: https://chemfiles.org/chemfiles/latest/selections.html
    """

    def __init__(self, selection):
        """
        Create a new :py:class:`Selection` from the given ``selection`` string.
        """
        ptr = self.ffi.chfl_selection(selection.encode("utf8"))
        super(Selection, self).__init__(ptr, is_const=False)

    def __copy__(self):
        return Selection.from_mutable_ptr(None, self.ffi.chfl_selection_copy(self.ptr))

    def __repr__(self):
        return f"Selection('{self.string}')"

    @property
    def size(self):
        """
        Get the size of this :py:class:`Selection`.

        The size of a selection is the number of atoms we are selecting
        together. This value is 1 for the 'atom' context, 2 for the 'pair' and
        'bond' context, 3 for the 'three' and 'angles' contextes and 4 for the
        'four' and 'dihedral' contextes.
        """
        size = c_uint64()
        self.ffi.chfl_selection_size(self.ptr, size)
        return size.value

    @property
    def string(self):
        """
        Get the selection string used to create this :py:class:`Selection`.
        """
        return _call_with_growing_buffer(
            lambda buffer, size: self.ffi.chfl_selection_string(self.ptr, buffer, size),
            initial=128,
        )

    def evaluate(self, frame):
        """
        Evaluate a :py:class:`Selection` for a given :py:class:`Frame`, and
        return a list of matching atoms, either as a list of index or a list
        of tuples of indexes.
        """
        matching = c_uint64()
        self.ffi.chfl_selection_evaluate(self.mut_ptr, frame.ptr, matching)

        matches = np.zeros(matching.value, chfl_match)
        self.ffi.chfl_selection_matches(self.mut_ptr, matches, matching)

        size = self.size
        result = []
        for match in matches:
            assert match[0] == size
            atoms = match[1]
            if size == 1:
                result.append(atoms[0])
            elif size == 2:
                result.append((atoms[0], atoms[1]))
            elif size == 3:
                result.append((atoms[0], atoms[1], atoms[2]))
            elif size == 4:
                result.append((atoms[0], atoms[1], atoms[2], atoms[3]))
        return result
