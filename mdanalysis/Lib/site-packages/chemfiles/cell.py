from ctypes import ARRAY, c_double
from enum import IntEnum

import numpy as np

from ._c_api import chfl_cellshape, chfl_vector3d
from .utils import ChemfilesError, CxxPointer


class CellShape(IntEnum):
    """
    Available cell shapes in Chemfiles:

    - ``CellType.Orthorhombic``: for cells where the three angles are 90°;
    - ``CellType.Triclinic``: for cells where the three angles may not be 90°;
    - ``CellType.Infinite``: for cells without periodic boundary conditions;
    """

    Orthorhombic = chfl_cellshape.CHFL_CELL_ORTHORHOMBIC
    Triclinic = chfl_cellshape.CHFL_CELL_TRICLINIC
    Infinite = chfl_cellshape.CHFL_CELL_INFINITE


class UnitCell(CxxPointer):
    """
    An :py:class:`UnitCell` represent the box containing the atoms, and its
    periodicity.

    An unit cell is fully represented by three lengths (a, b, c); and three
    angles (alpha, beta, gamma). The angles are stored in degrees, and the
    lengths in Angstroms. The cell angles are defined as follow: alpha is the
    angles between the cell vectors `b` and `c`; beta as the angle between `a`
    and `c`; and gamma as the angle between `a` and `b`.

    A cell also has a matricial representation, by projecting the three base
    vector into an orthonormal base. We choose to represent such matrix as an
    upper triangular matrix:

    .. code-block:: sh

        | a_x   b_x   c_x |
        |  0    b_y   c_y |
        |  0     0    c_z |
    """

    def __init__(self, lengths, angles=(90.0, 90.0, 90.0)):
        """
        Create a new :py:class:`UnitCell` with the given cell ``lengths`` and
        cell ``angles``. If ``lengths`` is a 3x3 matrix, it is taken to be the
        unit cell matrix, and ``angles`` is ignored.

        If the three angles are equal to 90.0, the new unit cell shape is
        ``CellShape.Orthorhombic``. Else it is ``CellShape.Infinite``.
        """
        lengths = np.array(lengths)
        if len(lengths.shape) == 1:
            lengths = chfl_vector3d(*lengths)
            angles = chfl_vector3d(*angles)
            ptr = self.ffi.chfl_cell(lengths, angles)
        else:
            if lengths.shape != (3, 3):
                raise ChemfilesError(
                    f"expected the cell matrix to have 3x3 shape, got {lengths.shape}"
                )
            matrix = ARRAY(chfl_vector3d, (3))()
            matrix[0][0] = lengths[0, 0]
            matrix[0][1] = lengths[0, 1]
            matrix[0][2] = lengths[0, 2]

            matrix[1][0] = lengths[1, 0]
            matrix[1][1] = lengths[1, 1]
            matrix[1][2] = lengths[1, 2]

            matrix[2][0] = lengths[2, 0]
            matrix[2][1] = lengths[2, 1]
            matrix[2][2] = lengths[2, 2]
            ptr = self.ffi.chfl_cell_from_matrix(matrix)

        super(UnitCell, self).__init__(ptr, is_const=False)

    def __copy__(self):
        return UnitCell.from_mutable_ptr(None, self.ffi.chfl_cell_copy(self.ptr))

    def __repr__(self):
        return """UnitCell(
    lengths=({:.9g}, {:.9g}, {:.9g}),
    angles=({:.7g}, {:.7g}, {:.7g})
)""".format(
            *(self.lengths + self.angles)
        )

    @property
    def lengths(self):
        """Get the three lengths of this :py:class:`UnitCell`, in Angstroms."""
        lengths = chfl_vector3d(0, 0, 0)
        self.ffi.chfl_cell_lengths(self.ptr, lengths)
        return lengths[0], lengths[1], lengths[2]

    @lengths.setter
    def lengths(self, lengths):
        """
        Set the three lengths of this :py:class:`UnitCell` to ``lengths``. The
        values should be in Angstroms.
        """
        a, b, c = lengths
        self.ffi.chfl_cell_set_lengths(self.mut_ptr, chfl_vector3d(a, b, c))

    @property
    def angles(self):
        """Get the three angles of this :py:class:`UnitCell`, in degrees."""
        angles = chfl_vector3d(0, 0, 0)
        self.ffi.chfl_cell_angles(self.ptr, angles)
        return angles[0], angles[1], angles[2]

    @angles.setter
    def angles(self, angles):
        """
        Set the three angles of this :py:class:`UnitCell` to ``alpha``,
        ``beta`` and ``gamma``. These values should be in degrees. Setting
        angles is only possible for ``CellShape.Triclinic`` cells.
        """
        alpha, beta, gamma = angles
        self.ffi.chfl_cell_set_angles(self.mut_ptr, chfl_vector3d(alpha, beta, gamma))

    @property
    def matrix(self):
        """
        Get this :py:class:`UnitCell` matricial representation.

        The matricial representation is obtained by aligning the a vector along
        the *x* axis and putting the b vector in the *xy* plane. This make the
        matrix an upper triangular matrix:

        .. code-block:: sh

            | a_x   b_x   c_x |
            |  0    b_y   c_y |
            |  0     0    c_z |
        """
        m = ARRAY(chfl_vector3d, 3)()
        self.ffi.chfl_cell_matrix(self.ptr, m)
        return np.array(
            (
                (m[0][0], m[0][1], m[0][2]),
                (m[1][0], m[1][1], m[1][2]),
                (m[2][0], m[2][1], m[2][2]),
            )
        )

    @property
    def shape(self):
        """Get the shape of this :py:class:`UnitCell`."""
        shape = chfl_cellshape()
        self.ffi.chfl_cell_shape(self.ptr, shape)
        return CellShape(shape.value)

    @shape.setter
    def shape(self, shape):
        """Set the shape of this :py:class:`UnitCell` to ``shape``."""
        self.ffi.chfl_cell_set_shape(self.mut_ptr, chfl_cellshape(shape))

    @property
    def volume(self):
        """Get the volume of this :py:class:`UnitCell`."""
        volume = c_double()
        self.ffi.chfl_cell_volume(self.ptr, volume)
        return volume.value

    def wrap(self, vector):
        """
        Wrap a ``vector`` in this :py:class:`UnitCell`, and return the wrapped
        vector.
        """
        vector = chfl_vector3d(vector[0], vector[1], vector[2])
        self.ffi.chfl_cell_wrap(self.ptr, vector)
        return (vector[0], vector[1], vector[2])
