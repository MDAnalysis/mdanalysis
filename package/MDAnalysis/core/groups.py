"""Containers for objects in MDA


"""

def make_group(top):
    """Generate the Group class with attributes according to the topology.

    """
    return type('Group', (GroupBase,), {})


def make_levelgroup(top, Groupclass, level):
    """Generate the AtomGroup class with attributes according to the topology.

    """
    if level == 'atom':
        levelgroup = 'AtomGroup'
        baseclass = AtomGroupBase
    elif level == 'residue':
        levelgroup = 'ResidueGroup'
        baseclass = ResidueGroupBase
    elif level == 'segment':
        levelgroup = 'SegmentGroup'
        baseclass = SegmentGroupBase 

    return type(levelgroup, (Groupclass, baseclass), {})


class GroupBase(object):
    def __init__(self, ix, u):
        # indices for the objects I hold
        self._ix = ix
        self._u = u
        self._cache = dict()

    @classmethod
    def _add_prop(cls, attr):
        """Add attr into the namespace for this class

        Arguments
        ---------
        attr - A TopologyAttr object
        """
        getter = lambda self: attr.__getitem__(self)
        setter = lambda self, values: attr.__setitem__(self, values)

        setattr(cls, attr.attrname,
                property(getter, setter, None, attr.__doc__))

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
    
    @property
    def atoms(self):
        return self._u.atoms[np.unique(self.indices)]

    @property
    def residues(self):
        return self._u.residues[np.unique(self.resindices)]

    @property
    def segments(self):
        return self._u.segments[np.unique(self.segindices)]


class AtomGroupBase(object):
    """AtomGroup base class.

    This class is used by a Universe for generating its Topology-specific
    AtomGroup class. All the TopologyAttr components are obtained from
    GroupBase, so this class only includes ad-hoc methods specific to
    AtomGroups.

    """
    level = 'atom'

    @property
    def position(self):
        """coordinates of the atom

        Get the current cartesian coordinates of the atom.

        :Returns: a (3,) shape numpy array
        """
        return self.universe.coord.positions[self.index]  # internal numbering starts at 0

    @position.setter
    def position(self, coords):
        """
        Set the current cartesian coordinates of the atom.
        @param coords: a 1x3 numpy array of {x,y,z} coordinates, or optionally
            a single scalar if you should want to set all coordinates to the same value.
        """
        self.universe.coord.positions[self.index, :] = coords  # internal numbering starts at 0

    @property
    def velocity(self):
        """Current velocity of the atom.

        :Returns: a (3,) shape numpy array

        A :exc:`~MDAnalysis.NoDataError` is raised if the trajectory
        does not contain velocities.

        .. versionadded:: 0.7.5
        """
        # TODO: Remove error checking here (and all similar below)
        # and add to Timestep
        try:
            return self.universe.coord.velocities[self.index]
        except (AttributeError, NoDataError):
            raise NoDataError("Timestep does not contain velocities")

    @velocity.setter
    def velocity(self, vals):
        """Set the current velocity of the atom.

        A :exc:`~MDAnalysis.NoDataError` is raised if the trajectory
        does not contain velocities.

        .. versionadded:: 0.9.2
        """
        try:
            self.universe.coord.velocities[self.index] = vals
        except (AttributeError, NoDataError):
            raise NoDataError("Timestep does not contain velocities")


class ResidueGroupBase(object):
    level = 'residue'


class SegmentGroupBase(object):
    level = 'segment'
