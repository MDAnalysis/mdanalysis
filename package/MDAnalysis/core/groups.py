"""Containers for objects in MDA


"""
import numpy as np

from . import selection
from . import flags

def make_group():
    """Generate the Group class with attributes according to the topology.

    """
    return type('Group', (GroupBase,), {})


def make_levelgroup(Groupclass, level):
    """Generate the *Group class at `level` with attributes according to the
    topology.

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


def make_levelcomponent(level):
    """Generate a copy of the Component class specified by `level`.

    """
    if level == 'atom':
        levelgroup = 'Atom'
        baseclass = AtomBase
    elif level == 'residue':
        levelgroup = 'Residue'
        baseclass = ResidueBase
    elif level == 'segment':
        levelgroup = 'Segment'
        baseclass = SegmentBase 

    return type(levelgroup, (baseclass,), {})


class GroupBase(object):
    """Base class from which a Universe's Group class is built.

    """
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
        if not isinstance(item, (int, np.int_)):
            return self.__class__(self._ix[item], self._u)
        else:
            return self._u._components[self.level](self._ix[item], self._u)

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
    
    @property
    def atoms(self):
        """Get a unique (non-repeating) AtomGroup of the atoms in the group.
        """
        return self._u.atoms[np.unique(self.indices)]

    @property
    def n_atoms(self):
        """Total number of unique atoms represented in the group.
        """
        return len(self.atoms)

    @property
    def residues(self):
        """Get a unique (non-repeating) ResidueGroup of the residues
        represented in the group.
        """
        return self._u.residues[np.unique(self.resindices)]

    @property
    def n_residues(self):
        """Total number of unique residues represented in the group.

        """
        return len(self.residues)

    @property
    def segments(self):
        """Get a unique (non-repeating) SegmentGroup of the segments 
        represented in the group.
        """
        return self._u.segments[np.unique(self.segindices)]
                                                
    @property
    def n_segments(self):
        """Total number of unique segments represented in the group.

        """
        return len(self.segments)


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
    def positions(self):
        """Coordinates of the atoms in the AtomGroup.

        The positions can be changed by assigning an array of the appropriate
        shape, i.e. either Nx3 to assign individual coordinates or 3, to assign
        the *same* coordinate to all atoms (e.g. ``ag.positions = array([0,0,0])``
        will move all particles to the origin).

        """
        return self._u.trajectory.ts.positions[self._ix]
    
    @positions.setter
    def positions(self, values):
        ts = self._u.trajectory.ts
        ts.positions[self._ix, :] = values

    @property
    def velocities(self):
        """Velocities of the atoms in the AtomGroup.

        The velocities can be changed by assigning an array of the appropriate
        shape, i.e. either Nx3 to assign individual velocities or 3 to assign
        the *same* velocity to all atoms (e.g. ``ag.velocity = array([0,0,0])``
        will give all particles zero velocity).

        Raises a :exc:`NoDataError` if the underlying
        :class:`~MDAnalysis.coordinates.base.Timestep` does not contain
        :attr:`~MDAnalysis.coordinates.base.Timestep.velocities`.

        """
        ts = self._u.trajectory.ts
        try:
            return np.array(ts.velocities[self._ix])
        except (AttributeError, NoDataError):
            raise NoDataError("Timestep does not contain velocities")

    @velocities.setter
    def velocities(self, values):
        ts = self._u.trajectory.ts
        try:
            ts.velocities[self._ix, :] = values
        except AttributeError:
            raise NoDataError("Timestep does not contain velocities")

    @property
    def forces(self):
        """Forces on each atom in the AtomGroup.

        The velocities can be changed by assigning an array of the appropriate
        shape, i.e. either Nx3 to assign individual velocities or 3 to assign
        the *same* velocity to all atoms (e.g. ``ag.velocity = array([0,0,0])``
        will give all particles zero velocity).


        """
        ts = self._u.trajectory.ts
        try:
            return ts.forces[self._ix]
        except (AttributeError, NoDataError):
            raise NoDataError("Timestep does not contain forces")

    @forces.setter
    def forces(self, values):
        ts = self._u.trajectory.ts
        try:
            ts.forces[self._ix, :] = forces
        except AttributeError:
            raise NoDataError("Timestep does not contain forces")

    def center_of_geometry(self, **kwargs):
        """Center of geometry (also known as centroid) of the selection.

        Keywords
        --------
          *pbc*
            ``True``: Move all atoms within the primary unit cell
                      before calculation [``False``]

        Notes
        -----
        The :class:`MDAnalysis.core.flags` flag *use_pbc* when set to
        ``True`` allows the *pbc* flag to be used by default.

        .. versionchanged:: 0.8 Added *pbc* keyword
        """
        pbc = kwargs.pop('pbc', flags['use_pbc'])
        if pbc:
            return np.sum(self.pack_into_box(inplace=False), axis=0) / len(self)
        else:
            return np.sum(self.positions, axis=0) / len(self)

    def select_atoms(self, sel, *othersel, **selgroups):
        """Select atoms using a CHARMM selection string.

        Returns an :class:`AtomGroup` with atoms sorted according to their
        index in the psf (this is to ensure that there aren't any duplicates,
        which can happen with complicated selections).

        Existing :class:`AtomGroup` objects can be passed as named arguments,
        which will then be available to the selection parser.

        Subselections can be grouped with parentheses.

        Example::
           >>> sel = universe.select_atoms("segid DMPC and not ( name H* or name O* )")
           >>> sel
           <AtomGroup with 3420 atoms>

           >>> universe.select_atoms("around 10 group notHO", notHO=sel)
           <AtomGroup with 1250 atoms>

        .. Note::

           If exact ordering of atoms is required (for instance, for
           :meth:`~AtomGroup.angle` or :meth:`~AtomGroup.dihedral`
           calculations) then one supplies selections *separately* in the
           required order. Also, when multiple :class:`AtomGroup` instances are
           concatenated with the ``+`` operator then the order of :class:`Atom`
           instances is preserved and duplicates are not removed.

        .. SeeAlso:: :ref:`selection-commands-label` for further details and examples.

        The selection parser understands the following CASE SENSITIVE *keywords*:

        **Simple selections**

            protein, backbone, nucleic, nucleicbackbone
                selects all atoms that belong to a standard set of residues;
                a protein is identfied by a hard-coded set of residue names so
                it  may not work for esoteric residues.
            segid *seg-name*
                select by segid (as given in the topology), e.g. ``segid 4AKE``
                or ``segid DMPC``
            resid *residue-number-range*
                resid can take a single residue number or a range of numbers. A
                range consists of two numbers separated by a colon (inclusive)
                such as ``resid 1:5``. A residue number ("resid") is taken
                directly from the topology.
            resnum *resnum-number-range*
                resnum is the canonical residue number; typically it is set to
                the residue id in the original PDB structure.
            resname *residue-name*
                select by residue name, e.g. ``resname LYS``
            name *atom-name*
                select by atom name (as given in the topology). Often, this is
                force field dependent. Example: ``name CA`` (for C&alpha; atoms)
                or ``name OW`` (for SPC water oxygen)
            type *atom-type*
                select by atom type; this is either a string or a number and
                depends on the force field; it is read from the topology file
                (e.g. the CHARMM PSF file contains numeric atom types). It has
                non-sensical values when a PDB or GRO file is used as a topology
            atom *seg-name*  *residue-number*  *atom-name*
                a selector for a single atom consisting of segid resid atomname,
                e.g. ``DMPC 1 C2`` selects the C2 carbon of the first residue of
                the DMPC segment
            altloc *alternative-location*
                a selection for atoms where alternative locations are available,
                which is often the case with high-resolution crystal structures
                e.g. `resid 4 and resname ALA and altloc B` selects only the
                atoms of ALA-4 that have an altloc B record.

        **Boolean**

            not
                all atoms not in the selection, e.g. ``not protein`` selects
                all atoms that aren't part of a protein

            and, or
                combine two selections according to the rules of boolean
                algebra, e.g. ``protein and not (resname ALA or resname LYS)``
                selects all atoms that belong to a protein, but are not in a
                lysine or alanine residue

        **Geometric**

            around *distance*  *selection*
                selects all atoms a certain cutoff away from another selection,
                e.g. ``around 3.5 protein`` selects all atoms not belonging to
                protein that are within 3.5 Angstroms from the protein
            point *x* *y* *z*  *distance*
                selects all atoms within a cutoff of a point in space, make sure
                coordinate is separated by spaces,
                e.g. ``point 5.0 5.0 5.0  3.5`` selects all atoms within 3.5
                Angstroms of the coordinate (5.0, 5.0, 5.0)
            prop [abs] *property*  *operator*  *value*
                selects atoms based on position, using *property*  **x**, **y**,
                or **z** coordinate. Supports the **abs** keyword (for absolute
                value) and the following *operators*: **<, >, <=, >=, ==, !=**.
                For example, ``prop z >= 5.0`` selects all atoms with z
                coordinate greater than 5.0; ``prop abs z <= 5.0`` selects all
                atoms within -5.0 <= z <= 5.0.
            sphzone *radius* *selection*
                Selects all atoms that are within *radius* of the center of
                geometry of *selection*
            sphlayer *inner radius* *outer radius* *selection*
                Similar to sphzone, but also excludes atoms that are within
                *inner radius* of the selection COG

        **Connectivity**

            byres *selection*
                selects all atoms that are in the same segment and residue as
                selection, e.g. specify the subselection after the byres keyword
            bonded *selection*
                selects all atoms that are bonded to selection
                eg: ``select name H bonded name O`` selects only hydrogens
                bonded to oxygens

        **Index**

            bynum *index-range*
                selects all atoms within a range of (1-based) inclusive indices,
                e.g. ``bynum 1`` selects the first atom in the universe;
                ``bynum 5:10`` selects atoms 5 through 10 inclusive. All atoms
                in the :class:`MDAnalysis.Universe` are consecutively numbered,
                and the index runs from 1 up to the total number of atoms.

        **Preexisting selections**

            group *group-name*
                selects the atoms in the :class:`AtomGroup` passed to the
                function as an argument named *group-name*. Only the atoms
                common to *group-name* and the instance :meth:`~select_atoms`
                was called from will be considered. *group-name* will be
                 included in the parsing just by comparison of atom indices.
                This means that it is up to the user to make sure they were
                defined in an appropriate :class:`Universe`.

            fullgroup *group-name*
                just like the ``group`` keyword with the difference that all the
                atoms of *group-name* are included. The resulting selection may
                therefore have atoms that were initially absent from the
                instance :meth:`~select_atoms` was called from.

        .. versionchanged:: 0.7.4
           Added *resnum* selection.
        .. versionchanged:: 0.8.1
           Added *group* and *fullgroup* selections.
        .. versionchanged:: 0.13.0
           Added *bonded* selection
        """
        atomgrp = selection.Parser.parse(sel, selgroups).apply(self)
        if othersel:
            # Generate a selection for each selection string
            for sel in othersel:
                atomgrp += selection.Parser.parse(sel, selgroups).apply(self)
        return atomgrp


class ResidueGroupBase(object):
    """ResidueGroup base class.

    This class is used by a Universe for generating its Topology-specific
    ResidueGroup class. All the TopologyAttr components are obtained from
    GroupBase, so this class only includes ad-hoc methods specific to
    ResidueGroups.

    """
    level = 'residue'

    def sequence(self, **kwargs):
        """Returns the amino acid sequence.

        The format of the sequence is selected with the keyword *format*:

        ============== ============================================
        *format*       description
        ============== ============================================
        'SeqRecord'    :class:`Bio.SeqRecord.SeqRecord` (default)
        'Seq'          :class:`Bio.Seq.Seq`
        'string'       string
        ============== ============================================

        The sequence is returned by default (keyword ``format = 'SeqRecord'``)
        as a :class:`Bio.SeqRecord.SeqRecord` instance, which can then be
        further processed. In this case, all keyword arguments (such as the
        *id* string or the *name* or the *description*) are directly passed to
        :class:`Bio.SeqRecord.SeqRecord`.

        If the keyword *format* is set to ``'Seq'``, all *kwargs* are ignored
        and a :class:`Bio.Seq.Seq` instance is returned. The difference to the
        record is that the record also contains metadata and can be directly
        used as an input for other functions in :mod:`Bio`.

        If the keyword *format* is set to ``'string'``, all *kwargs* are
        ignored and a Python string is returned.

        .. rubric:: Example: Write FASTA file

        Use :func:`Bio.SeqIO.write`, which takes sequence records::

           import Bio.SeqIO

           # get the sequence record of a protein component of a Universe
           protein = u.select_atoms("protein")
           record = protein.sequence(id="myseq1", name="myprotein")

           Bio.SeqIO.write(record, "single.fasta", "fasta")

        A FASTA file with multiple entries can be written with ::

           Bio.SeqIO.write([record1, record2, ...], "multi.fasta", "fasta")

        :Keywords:
            *format*

                - ``"string"``: return sequence as a string of 1-letter codes
                - ``"Seq"``: return a :class:`Bio.Seq.Seq` instance
                - ``"SeqRecord"``: return a :class:`Bio.SeqRecord.SeqRecord`
                  instance

                Default is ``"SeqRecord"``

             *id*
                Sequence ID for SeqRecord (should be different for different
                sequences)
             *name*
                Name of the protein.
             *description*
                Short description of the sequence.
             *kwargs*
                Any other keyword arguments that are understood by
                :class:`Bio.SeqRecord.SeqRecord`.

        :Raises: :exc:`ValueError` if a residue name cannot be converted to a
                 1-letter IUPAC protein amino acid code; make sure to only
                 select protein residues. Raises :exc:`TypeError` if an unknown
                 *format* is selected.

        .. versionadded:: 0.9.0
        """
        import Bio.Seq
        import Bio.SeqRecord
        import Bio.Alphabet
        formats = ('string', 'Seq', 'SeqRecord')

        format = kwargs.pop("format", "SeqRecord")
        if format not in formats:
            raise TypeError("Unknown format='{0}': must be one of: {1}".format(
                    format, ", ".join(formats)))
        try:
            sequence = "".join([util.convert_aa_code(r) for r in self.residues.resnames])
        except KeyError as err:
            raise ValueError("AtomGroup contains a residue name '{0}' that "
                             "does not have a IUPAC protein 1-letter "
                             "character".format(err.message))
        if format == "string":
            return sequence
        seq = Bio.Seq.Seq(sequence, alphabet=Bio.Alphabet.IUPAC.protein)
        if format == "Seq":
            return seq
        return Bio.SeqRecord.SeqRecord(seq, **kwargs)


class SegmentGroupBase(object):
    """SegmentGroup base class.

    This class is used by a Universe for generating its Topology-specific
    SegmentGroup class. All the TopologyAttr components are obtained from
    GroupBase, so this class only includes ad-hoc methods specific to
    SegmentGroups.

    """
    level = 'segment'


class ComponentBase(object):
    """Base class from which a Universe's Component class is built.

    Components are the individual objects that are found in Groups.
    """
    def __init__(self, ix, u):
        # index of component
        self._ix = ix
        self._u = u

    def __repr__(self):
        return ("<{} {}>"
                "".format(self.level.capitalize(), self._ix))

    # TODO: put in mixin with GroupBase method of same name
    @classmethod
    def _add_prop(cls, attr):
        """Add attr into the namespace for this class

        Arguments
        ---------
        attr 
            TopologyAttr object to add
        """
        getter = lambda self: attr.__getitem__(self)
        setter = lambda self, values: attr.__setitem__(self, values)

        setattr(cls, attr.singular,
                property(getter, setter, None, attr.__doc__))


class AtomBase(ComponentBase):
    """Atom base class.

    This class is used by a Universe for generating its Topology-specific Atom
    class. All the TopologyAttr components are obtained from ComponentBase, so
    this class only includes ad-hoc methods specific to Atoms.

    """
    level = 'atom'

    @property
    def residue(self):
        residueclass = self._u._components['residue']
        return residueclass(self._u._topology.resindices[self],
                            self._u)

    @property
    def segment(self):
        segmentclass = self._u._components['segment']
        return segmentclass(self._u._topology.segindices[self],
                            self._u)


class ResidueBase(ComponentBase):
    """Residue base class.

    This class is used by a Universe for generating its Topology-specific
    Residue class. All the TopologyAttr components are obtained from
    ComponentBase, so this class only includes ad-hoc methods specific to
    Residues.

    """
    level = 'residue'

    @property
    def atoms(self):
        atomsclass = self._u._groups['atom']

        # we need to pass a fake residue with a numpy array as its self._ix for
        # downward translation tables to work; this is because they accept
        # arrays only for performance
        r_proxy = self.__class__(np.array([self._ix]), self._u)
        return atomsclass(self._u._topology.indices[r_proxy],
                          self._u)

    @property
    def segment(self):
        segmentclass = self._u._components['segment']
        return segmentclass(self._u._topology.segindices[self],
                            self._u)


class SegmentBase(ComponentBase):
    """Segment base class.

    This class is used by a Universe for generating its Topology-specific Segment
    class. All the TopologyAttr components are obtained from ComponentBase, so
    this class only includes ad-hoc methods specific to Segments.

    """
    level = 'segment'

    @property
    def atoms(self):
        atomsclass = self._u._groups['atom']

        # we need to pass a fake segment with a numpy array as its self._ix for
        # downward translation tables to work; this is because they accept
        # arrays only for performance
        s_proxy = self.__class__(np.array([self._ix]), self._u)
        return atomsclass(self._u._topology.indices[s_proxy],
                          self._u)

    @property
    def residues(self):
        residuesclass = self._u._groups['residue']

        # we need to pass a fake residue with a numpy array as its self._ix for
        # downward translation tables to work; this is because they accept arrays only
        # for performance
        s_proxy = self.__class__(np.array([self._ix]), self._u)
        return residuesclass(self._u._topology.resindices[s_proxy],
                             self._u)
