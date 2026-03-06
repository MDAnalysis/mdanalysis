"""
This module contains classes for reading and writing (single frame) NAMD binary
files.  The main functionality is currently aimed at manipulation rather than
analysis.
"""
from .. import unit as u
from struct import unpack, pack

import numpy as np

class NamdBinFile:
    """From the NAMD manual:

        NAMD uses a trivial double-precision binary file format for
        coordinates, velocities, and forces ... The file consists of the atom
        count as a 32-bit integer followed by all three position or velocity
        components for each atom as 64-bit double-precision floating point ...

    The main attributes are the number of atom entries (natom) and a (flat)
    numpy array of size 3*natom (values).  The meaning of "values" is
    effectively arbitrary, but for convenience derived classes are provided
    which alias the values to more descriptive names (e.g. "coordinates").

    See Also
    --------
    :class:`NamdBinCoor`
    :class:`NamdBinVel`
    """
    _SCALE_FACTOR = 1.0
    def __init__(self, values=[]):
        self._values = np.asarray(values,np.float64) * self._SCALE_FACTOR

    @property
    def natom(self):
        """The current number of atom entries."""
        return int(self._values.size / 3)

    @classmethod
    def read(cls, fname):
        """Return an object from the values in a NAMD binary file."""
        infile = open(fname,'rb')
        natoms = int(unpack('i',infile.read(4))[0])
        values = [unpack('d',infile.read(8))[0] for n in range(3*natoms)]
        infile.close()
        return cls(values)

    def write(self, fname):
        """Write the current attributes to a file."""
        outfile = open(fname,'wb')
        outfile.write(pack('i',self.natom))
        for x in self._values / self._SCALE_FACTOR:
            outfile.write(pack('d',x))
        outfile.close()

    def delatoms(self, indices):
        """Delete entries corresponding to the given atom indices."""
        del_indices = np.atleast_1d(np.asarray(indices,np.int32))
        mask = np.ones(self.natom,np.bool)
        mask[del_indices] = False
        newvalues = np.zeros(3*(self.natom - del_indices.size),np.float64)
        newvalues += self._values.reshape((self.natom,3))[mask].flatten()
        self._values = newvalues

    def insertatoms(self, start_index, natoms, values=None):
        """Insert space for natom entries beginning at start_index. If
        specified, give them the provided values, otherwise set them to zero.
        """
        if values is not None:
            values = np.asarray(values)
            assert values.size == 3*natoms
        else:
            values = np.zeros(3*natoms)
        hinge = 3*start_index
        newvalues = np.concatenate((self._values[:hinge],values,self._values[hinge:]))
        self._values = newvalues

    def copyatoms(self, start_index, natoms):
        """Convenience function, same as insertatoms() but set 'values' to
        be the same as the previous natoms' values (i.e. make a copy of them).
        """
        values = self._values[3*(start_index-natoms):3*start_index]
        self.insertatoms(start_index,natoms,values)


class NamdBinCoor(NamdBinFile):
    """ Class to read or write NAMD "bincoordinates" files. """
    @property
    def coordinates(self):
        return self._values.reshape((-1, self.natom, 3))

    @coordinates.setter
    def coordinates(self, value):
        if u.is_quantity(value):
            value = value.value_in_unit(u.angstroms)
        self._values = np.array(value).flatten()

class NamdBinVel(NamdBinFile):
    """ Class to read or write NAMD "binvelocities" files. """

    _SCALE_FACTOR = 20.45482706

    @property
    def velocities(self):
        return self._values.reshape((-1, self.natom, 3))

    @velocities.setter
    def velocities(self, value):
        if u.is_quantity(value):
            value = value.value_in_unit(u.angstroms/u.picosecond)
        self._values = np.array(value).flatten()
