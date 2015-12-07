"""Containers for objects in MDA


"""

class Group(object):
    level = ''

    def __init__(self, ix, u):
        # indices for the objects I hold
        self._ix = ix
        self._u = u
        self._cache = dict()

    def _addprop(self, attr):
        getter = lambda self: attr.__getitem__(self)
        setter = lambda self, values: attr.__setitem__(self, values)

        setattr(self.__class__, attr.attrname,
                property(getter, setter))

    def __len__(self):
        return len(self._ix)

    def __getitem__(self, item):
        # supports
        # - integer access
        # - boolean slicing
        # - fancy indexing
        # because our _ix attribute is a numpy array
        # it can be sliced by all of these already,
        # so just return ourselves sliced by the item
        return self.__class__(self._ix[item], self._u)

    def __repr__(self):
        return ("<{}Group with {} {}s>"
                "".format(self.level.capitalize(), len(self), self.level))

    def __add__(self, other):
        if self.level != other.level:
            raise TypeError("Can't add different level objects")
        if not self._u is other._u:
            raise ValueError("Can't add objects from different Universe")
        return self.__class__(np.concatenate([self._ix, other._ix]), self._u)

    @property
    def universe(self):
        return self._u


class AtomGroup(Group):
    level = 'atom'

    @property
    def names(self):
        return self._u._topology.atomnames.get_atoms(self._ix)

    @names.setter
    def names(self, values):
        return self._u._topology.atomnames.set_atoms(self._ix, values)

    @property
    def ids(self):
        return self._u._topology.atomids.get_atoms(self._ix)

    @ids.setter
    def ids(self, values):
        return self._u._topology.atomids.set_atoms(self._ix, values)

    @property
    def masses(self):
        return self._u._topology.masses.get_atoms(self._ix)

    @masses.setter
    def masses(self, values):
        return self._u._topology.masses.set_atoms(self._ix, values)


class ResidueGroup(Group):
    level = 'residue'

    @property
    def names(self):
        return self._u._topology.atomnames.get_residues(self._ix)

    @property
    def ids(self):
        return self._u._topology.atomids.get_residues(self._ix)

    @property
    def resids(self):
        return self._u._topology.resids.get_residues(self._ix)

    @resids.setter
    def resids(self, values):
        return self._u._topology.resids.set_residues(self._ix, values)

    @property
    def resnames(self):
        return self._u._topology.resnames.get_residues(self._ix)

    @resnames.setter
    def resnames(self, values):
        return self._u._topology.resnames.set_residues(self._ix, values)

    @property
    def masses(self):
        return self._u._topology.masses.get_residues(self._ix)

class SegmentGroup(Group):
    level = 'segment'

    @property
    def id(self):
        return self._u._topology.segids.get_segments(self._ix)
