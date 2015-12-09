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

    # TODO: need to generalize for checking Levelobject isin Group
    def __contains__(self, other):
        # If the number of atoms is very large, create a dictionary cache for lookup
        if len(self) > self._atomcache_size and not 'atoms' in self._cache:
            self._cache['atoms'] = dict(((x, None) for x in self.__atoms))
        try:
            return other in self._cache['atoms']
        except KeyError:
            return other in self._atoms

    @property
    def universe(self):
        return self._u
    
    # rethink: perhaps a set only uphill? Need to place in separate *GroupBase
    # classes if so

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


    # TODO: CHECK THAT THIS IS APPROPRIATE IN NEW SCHEME
    def __getattr__(self, name):
        try:
            return self._get_named_atom(name)
        except SelectionError:
            raise AttributeError("'{0}' object has no attribute '{1}'".format(
                    self.__class__.__name__, name))

    # TODO: fix me; we need Atom objects!
    def _get_named_atom(self, name):
        """Get all atoms with name *name* in the current AtomGroup.

        For more than one atom it returns a list of :class:`Atom`
        instance. A single :class:`Atom` is returned just as such. If
        no atoms are found, a :exc:`SelectionError` is raised.

        .. versionadded:: 0.9.2
        """
        # There can be more than one atom with the same name
        atomlist = [atom for atom in self._atoms if name == atom.name]
        if len(atomlist) == 0:
            raise SelectionError("No atoms with name '{0}'".format(name))
        elif len(atomlist) == 1:
            return atomlist[0]  # XXX: keep this, makes more sense for names
        else:
            return AtomGroup(atomlist)  # XXX: but inconsistent (see residues and Issue 47)

    @property
    def atoms(self):
        """Get a unique (non-repeating) AtomGroup of the atoms in the group.
        """
        return self._u.atoms[np.unique(self.indices)]

    @property
    def n_atoms(self):
        """Total number of atoms in the group.
        """
        return len(self)


class ResidueGroupBase(object):
    level = 'residue'

    @property
    def atoms(self):
        """Get an AtomGroup of the atoms in each residue.

        Atoms are ordered according the order of residues in the group.
        """
        return self._u.atoms[self.indices]

    @property
    def n_atoms(self):
        """Total number of atoms in the group, including repeats.
        """
        return len(self.atoms)


class SegmentGroupBase(object):
    level = 'segment'

    @property
    def atoms(self):
        """Get an AtomGroup of the atoms in each segment.

        Atoms are ordered according the order of segments in the group.
        """
        return self._u.atoms[self.indices]
