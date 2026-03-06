"""
This module contains a chamber prmtop class that will read in all
parameters and allow users to manipulate that data and write a new
prmtop object.
"""
import copy as _copy
import warnings
from math import pi, sqrt

from ..constants import DEG_TO_RAD, PrmtopPointers, RAD_TO_DEG, SMALL, TINY
from ..exceptions import AmberError, AmberWarning
from ..topologyobjects import BondType, ExtraPoint, Improper, ImproperType, UreyBradley
from ._amberparm import AmberParm

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class ChamberParm(AmberParm):
    """Chamber Topology (parm7 format) class.

    Gives low, and some high, level access to topology data or
    interact with some of the high-level classes comprising the system
    topology and parameters.  The ChamberParm class uses the same
    attributes that the AmberParm class uses, and only the ones unique
    to ChamberParm will be shown below.

    Parameters
    ----------
    prm_name : str, optional
        If provided, this file is parsed and the data structures will be loaded
        from the data in this file
    xyz : str or array, optional
        If provided, the coordinates and unit cell dimensions from the provided
        Amber inpcrd/restart file will be loaded into the molecule, or the
        coordinates will be loaded from the coordinate array
    box : array, optional
        If provided, the unit cell information will be set from the provided
        unit cell dimensions (a, b, c, alpha, beta, and gamma, respectively)

    Attributes
    ----------
    LJ_14_radius : list(float)
        The same as LJ_radius, except specific for 1-4 nonbonded parameters,
        which may differ in the CHARMM force field
    LJ_14_depth : list(float)
        The same as LJ_depth, except specific for 1-4 nonbonded parameters,
        which may differ in the CHARMM force field
    urey_bradleys : TrackedList(UreyBradley)
        List of Urey-Bradley terms between two atoms in a valence angle
    impropers : TrackedList(Improper)
        List of CHARMM-style improper torsions
    cmaps : TrackedList(Cmap)
        List of coupled-torsion correction map parameters
    urey_bradley_types : TrackedList(UreyBradleyType)
        List of parameters defining the Urey-Bradley terms
    improper_types : TrackedList(Improper)
        List of parameters defining the Improper terms
    cmap_types : TrackedList(CmapType)
        List of parameters defining the CMAP terms
    chamber : bool=True
        On ChamberParm instances, this is always True to indicate that it is a
        CHAMBER-style topology file
    amoeba : bool=False
        On ChamberParm instances, this is always False to indicate that it is
        not an AMOEBA-style topology file
    has_cmap : bool
        True if CMAP parameters are present in this system; False otherwise

    See Also
    --------
    :class:`AmberParm`

    """

    _cmap_prefix = "CHARMM_"

    #===================================================

    def initialize_topology(self, xyz=None, box=None):
        """
        Initializes topology data structures, like the list of atoms, bonds,
        etc., after the topology file has been read. The following methods are
        called:
        """
        self.LJ_14_radius = []
        self.LJ_14_depth = []
        AmberParm.initialize_topology(self, xyz, box)

    #===================================================

    def _copy_lj_data(self, other):
        """ Copies Lennard-Jones lists and dicts from myself to a copy """
        super(ChamberParm, self)._copy_lj_data(other)
        other.LJ_14_radius = _copy.copy(self.LJ_14_radius)
        other.LJ_14_depth = _copy.copy(self.LJ_14_depth)
        for atom in other.atoms:
            other.LJ_14_radius[atom.nb_idx-1] = atom.atom_type.rmin_14
            other.LJ_14_depth[atom.nb_idx-1] = atom.atom_type.epsilon_14

    #===================================================

    def load_pointers(self):
        """
        Loads the data in POINTERS section into a pointers dictionary with each
        key being the pointer name according to http://ambermd.org/formats.html
        """
        AmberParm.load_pointers(self)
        # Other pointers
        nub, nubtypes = self.parm_data['CHARMM_UREY_BRADLEY_COUNT'][:2]
        self.pointers['NUB'] = nub
        self.pointers['NUBTYPES'] = nubtypes
        self.pointers['NIMPHI'] = self.parm_data['CHARMM_NUM_IMPROPERS'][0]
        self.pointers['NIMPRTYPES'] = self.parm_data['CHARMM_NUM_IMPR_TYPES'][0]

    #===================================================

    def load_structure(self):
        """
        Loads all of the topology instance variables. This is necessary if we
        actually want to modify the topological layout of our system
        (like deleting atoms)
        """
        super(ChamberParm, self).load_structure()
        self._load_urey_brad_info()
        self._load_improper_info()
        # All of our angles are urey-bradley types
        for angle in self.angles: angle.funct = 5
        super(ChamberParm, self).unchange()

    #===================================================

    @classmethod
    def from_structure(cls, struct, copy=False):
        """
        Take a Structure instance and initialize a ChamberParm instance from
        that data.

        Parameters
        ----------
        struct : Structure
            The input structure from which to construct a ChamberParm instance
        copy : bool
            If True, the input struct is deep-copied to make sure it does not
            share any objects with the original ``struct``. Default is False

        Returns
        -------
        inst : :class:`ChamberParm`
            The ChamberParm instance derived from the input structure

        Notes
        -----
        Due to the nature of the prmtop file, struct almost *always* returns a
        deep copy. The one exception is when struct is already of type
        :class:`ChamberParm`, in which case the original object is returned
        unless ``copy`` is ``True``.
        """
        if isinstance(struct, cls):
            if copy:
                return _copy.copy(struct)
            return struct
        if (struct.rb_torsions or struct.trigonal_angles or struct.pi_torsions
                or struct.out_of_plane_bends or struct.stretch_bends
                or struct.torsion_torsions or struct.multipole_frames):
            raise TypeError('ChamberParm does not support all potential terms '
                             'defined in the input Structure')
        inst = struct.copy(cls, split_dihedrals=True)
        inst.update_dihedral_exclusions()
        inst._add_missing_13_14(ignore_inconsistent_vdw=True)
        inst.pointers = {}
        inst.LJ_types = {}
        nbfixes = inst.atoms.assign_nbidx_from_types()
        # Give virtual sites a name that Amber understands
        for atom in inst.atoms:
            if isinstance(atom, ExtraPoint): atom.type = 'EP'
        # Fill the Lennard-Jones arrays/dicts
        ntyp = 0
        for atom in inst.atoms:
            inst.LJ_types[atom.type] = atom.nb_idx
            ntyp = max(ntyp, atom.nb_idx)
        inst.LJ_radius = [0 for i in range(ntyp)]
        inst.LJ_depth = [0 for i in range(ntyp)]
        inst.LJ_14_radius = [0 for i in range(ntyp)]
        inst.LJ_14_depth = [0 for i in range(ntyp)]
        for atom in inst.atoms:
            inst.LJ_radius[atom.nb_idx-1] = atom.rmin
            inst.LJ_depth[atom.nb_idx-1] = atom.epsilon
            inst.LJ_14_radius[atom.nb_idx-1] = atom.rmin_14
            inst.LJ_14_depth[atom.nb_idx-1] = atom.epsilon_14
        inst._add_standard_flags()
        inst.pointers['NATOM'] = len(inst.atoms)
        inst.parm_data['POINTERS'][PrmtopPointers.NATOM] = len(inst.atoms)
        inst.box = _copy.copy(struct.box)
        if struct.box is None:
            inst.parm_data['POINTERS'][PrmtopPointers.IFBOX] = 0
            inst.pointers['IFBOX'] = 0
        elif (abs(struct.box[3] - 90) > TINY or abs(struct.box[4] - 90) > TINY
                or abs(struct.box[5] - 90) > TINY):
            inst.parm_data['POINTERS'][PrmtopPointers.IFBOX] = 2
            inst.pointers['IFBOX'] = 2
            inst.parm_data['BOX_DIMENSIONS'] = [struct.box[3]] + list(struct.box[:3])
        else:
            inst.parm_data['POINTERS'][PrmtopPointers.IFBOX] = 1
            inst.pointers['IFBOX'] = 1
            inst.parm_data['BOX_DIMENSIONS'] = [90] + list(struct.box[:3])
        # pmemd likes to skip torsions with periodicities of 0, which may be
        # present as a way to hack entries into the 1-4 pairlist. See
        # https://github.com/ParmEd/ParmEd/pull/145 for discussion. The solution
        # here is to simply set that periodicity to 1.
        for dt in inst.dihedral_types:
            if dt.phi_k == 0 and dt.per == 0:
                dt.per = 1.0
            elif dt.per == 0:
                warnings.warn('Periodicity of 0 detected with non-zero force constant. Changing '
                              'periodicity to 1 and force constant to 0 to ensure 1-4 nonbonded '
                              'pairs are properly identified. This might cause a shift in the '
                              'energy, but will leave forces unaffected', AmberWarning)
                dt.phi_k = 0.0
                dt.per = 1.0
        inst.remake_parm()
        inst._set_nonbonded_tables(nbfixes)
        del inst.adjusts[:]
        inst.parm_data['FORCE_FIELD_TYPE'] = fftype = []
        fftype.extend([1, 'CHARMM force field: No FF information parsed...'])

        return inst

    #===================================================

    def remake_parm(self):
        """
        Fills :attr:`parm_data` from the data in the parameter and topology
        arrays (e.g., :attr:`atoms`, :attr:`bonds`, :attr:`bond_types`, ...)
        """
        # Get rid of terms containing deleted atoms and empty residues
        self.prune_empty_terms()
        self.residues.prune()
        self.rediscover_molecules()

        # Transfer information from the topology lists
        self._xfer_atom_info()
        self._xfer_residue_info()
        self._xfer_bond_info()
        self._xfer_angle_info()
        self._xfer_dihedral_info()
        self._xfer_urey_bradley_properties()
        self._xfer_improper_properties()
        self._xfer_cmap_properties()
        # Load the pointers dict
        self.load_pointers()
        # Mark atom list as unchanged
        super(ChamberParm, self).unchange()

    #===================================================

    def fill_LJ(self):
        """
        Fills the LJ_radius, LJ_depth arrays and LJ_types dictionary with data
        from LENNARD_JONES_ACOEF and LENNARD_JONES_BCOEF sections of the prmtop
        files, by undoing the canonical combining rules.
        """
        AmberParm.fill_LJ(self)

        acoef = self.parm_data['LENNARD_JONES_14_ACOEF']
        bcoef = self.parm_data['LENNARD_JONES_14_BCOEF']
        ntypes = self.pointers['NTYPES']

        self.LJ_14_radius = []  # empty LJ_radii so it can be re-filled
        self.LJ_14_depth = []   # empty LJ_depths so it can be re-filled
        one_sixth = 1.0 / 6.0 # we need to raise some numbers to the 1/6th power

        for i in range(ntypes):
            lj_index = self.parm_data["NONBONDED_PARM_INDEX"][ntypes*i+i] - 1
            if acoef[lj_index] < 1.0e-6:
                self.LJ_14_radius.append(0)
                self.LJ_14_depth.append(0)
            else:
                factor = 2 * acoef[lj_index] / bcoef[lj_index]
                self.LJ_14_radius.append(pow(factor, one_sixth) * 0.5)
                self.LJ_14_depth.append(bcoef[lj_index] / 2 / factor)

    #===================================================

    def recalculate_LJ(self):
        """
        Takes the values of the LJ_radius and LJ_depth arrays and recalculates
        the LENNARD_JONES_A/BCOEF topology sections from the canonical
        combining rules.
        """
        AmberParm.recalculate_LJ(self)
        ntypes = self.pointers['NTYPES']
        acoef = self.parm_data['LENNARD_JONES_14_ACOEF']
        bcoef = self.parm_data['LENNARD_JONES_14_BCOEF']
        for i in range(ntypes):
            for j in range(i, ntypes):
                index = self.parm_data['NONBONDED_PARM_INDEX'][ntypes*i+j] - 1
                rij = self.LJ_14_radius[i] + self.LJ_14_radius[j]
                wdij = sqrt(self.LJ_14_depth[i] * self.LJ_14_depth[j])
                acoef[index] = wdij * rij ** 12
                bcoef[index] = 2 * wdij * rij**6

    #===================================================

    @property
    def chamber(self):
        return True

    @property
    def amoeba(self):
        return False

    #===========  PRIVATE INSTANCE METHODS  ============

    def _load_urey_brad_info(self):
        """ Loads the Urey-Bradley types and array """
        del self.urey_bradleys[:]
        del self.urey_bradley_types[:]
        for k, req in zip(self.parm_data['CHARMM_UREY_BRADLEY_FORCE_CONSTANT'],
                          self.parm_data['CHARMM_UREY_BRADLEY_EQUIL_VALUE']):
            self.urey_bradley_types.append(
                    BondType(k, req, self.urey_bradley_types)
            )
        it = iter(self.parm_data['CHARMM_UREY_BRADLEY'])
        for i, j, k in zip(it, it, it):
            self.urey_bradleys.append(
                    UreyBradley(self.atoms[i-1], self.atoms[j-1],
                                self.urey_bradley_types[k-1])
            )

    #===================================================

    def _load_improper_info(self):
        """ Loads the CHARMM Improper types and array """
        del self.impropers[:]
        del self.improper_types[:]
        for k, eq in zip(self.parm_data['CHARMM_IMPROPER_FORCE_CONSTANT'],
                         self.parm_data['CHARMM_IMPROPER_PHASE']):
            # Previous versions of ParmEd stored improper phases as degrees,
            # whereas it should really be stored in radians. So do a simple
            # heuristic check to see if a conversion is necessary so we support
            # all versions.
            eq = eq * RAD_TO_DEG if abs(eq) <= 2*pi else eq
            self.improper_types.append(
                    ImproperType(k, eq, self.improper_types)
            )
        it = iter(self.parm_data['CHARMM_IMPROPERS'])
        for i, j, k, l, m in zip(it, it, it, it, it):
            self.impropers.append(
                    Improper(self.atoms[i-1], self.atoms[j-1], self.atoms[k-1],
                             self.atoms[l-1], self.improper_types[m-1])
            )
        # Make sure that if we have a comment in the CHARMM impropers, we fix it
        # to say the units are in radians
        for i in range(len(self.parm_comments.get('CHARMM_IMPROPER_PHASE', []))):
            comment = self.parm_comments['CHARMM_IMPROPER_PHASE'][i]
            if 'degrees' in comment:
                self.parm_comments['CHARMM_IMPROPER_PHASE'][i] = \
                        comment.replace('degrees', 'radians')

    #===================================================

    def _check_section_lengths(self):
        """ Make sure each section has the necessary number of entries """
        super(ChamberParm, self)._check_section_lengths()

        def check_length(key, length, required=True):
            if not required and not key in self.parm_data: return
            if len(self.parm_data[key]) != length:
                raise AmberError('FLAG %s has %d elements; expected %d' %
                                 (key, len(self.parm_data[key]), length))

        check_length('CHARMM_UREY_BRADLEY_COUNT', 2)
        check_length('CHARMM_UREY_BRADLEY', self.pointers['NUB']*3)
        check_length('CHARMM_UREY_BRADLEY_FORCE_CONSTANT',
                     self.pointers['NUBTYPES'])
        check_length('CHARMM_UREY_BRADLEY_EQUIL_VALUE',
                     self.pointers['NUBTYPES'])
        check_length('CHARMM_NUM_IMPROPERS', 1)
        check_length('CHARMM_IMPROPERS', self.pointers['NIMPHI']*5)
        check_length('CHARMM_NUM_IMPR_TYPES', 1)
        check_length('CHARMM_IMPROPER_FORCE_CONSTANT',
                     self.pointers['NIMPRTYPES'])
        check_length('CHARMM_IMPROPER_PHASE', self.pointers['NIMPRTYPES'])

        ntypes = self.pointers['NTYPES']
        check_length('LENNARD_JONES_14_ACOEF', ntypes*(ntypes+1)//2)
        check_length('LENNARD_JONES_14_BCOEF', ntypes*(ntypes+1)//2)

    #===================================================

    def _xfer_urey_bradley_properties(self):
        """
        Sets the various topology file section data from the Urey-Bradley arrays
        """
        data = self.parm_data
        for urey_type in self.urey_bradley_types:
            urey_type.used = False
        for urey in self.urey_bradleys:
            urey.type.used = True
        self.urey_bradley_types.prune_unused()
        data['CHARMM_UREY_BRADLEY_FORCE_CONSTANT'] = \
                [type.k for type in self.urey_bradley_types]
        data['CHARMM_UREY_BRADLEY_EQUIL_VALUE'] = \
                [type.req for type in self.urey_bradley_types]
        data['CHARMM_UREY_BRADLEY'] = bond_array = []
        for urey in self.urey_bradleys:
            bond_array.extend([urey.atom1.idx+1, urey.atom2.idx+1,
                               urey.type.idx+1])
        nub = len(self.urey_bradleys)
        nubt = len(self.urey_bradley_types)
        data['CHARMM_UREY_BRADLEY_COUNT'] = [nub, nubt]
        self.pointers['NUB'] = nub
        self.pointers['NUBTYPES'] = nubt

    #===================================================

    def _xfer_improper_properties(self):
        """ Sets the topology file section data from the improper arrays """
        data = self.parm_data
        for improper_type in self.improper_types:
            improper_type.used = False
        for improper in self.impropers:
            improper.type.used = True
        self.improper_types.prune_unused()
        data['CHARMM_IMPROPER_FORCE_CONSTANT'] = \
                [type.psi_k for type in self.improper_types]
        data['CHARMM_IMPROPER_PHASE'] = \
                [type.psi_eq*DEG_TO_RAD for type in self.improper_types]
        data['CHARMM_IMPROPERS'] = improper_array = []
        for imp in self.impropers:
            improper_array.extend([imp.atom1.idx+1, imp.atom2.idx+1,
                                   imp.atom3.idx+1, imp.atom4.idx+1,
                                   imp.type.idx+1])
        data['CHARMM_NUM_IMPROPERS'] = [len(self.impropers)]
        data['CHARMM_NUM_IMPR_TYPES'] = [len(self.improper_types)]
        self.pointers['NIMPHI'] = len(improper_array)
        self.pointers['NIMPRTYPES'] = len(self.improper_types)

    #===================================================

    def _add_standard_flags(self):
        """ Adds all of the standard flags to the parm_data array """
        self.set_version()
        self.add_flag('CTITLE', '20a4', num_items=0)
        self.add_flag('POINTERS', '10I8', num_items=31)
        self.add_flag('FORCE_FIELD_TYPE', 'i2,a78', num_items=0)
        self.add_flag('ATOM_NAME', '20a4', num_items=0)
        self.add_flag('CHARGE', '3E24.16', num_items=0,
                      comments=['Atomic charge multiplied by sqrt(332.0716D0) (CCELEC)'])
        self.add_flag('ATOMIC_NUMBER', '10I8', num_items=0)
        self.add_flag('MASS', '5E16.8', num_items=0)
        self.add_flag('ATOM_TYPE_INDEX', '10I8', num_items=0)
        self.add_flag('NUMBER_EXCLUDED_ATOMS', '10I8', num_items=0)
        self.add_flag('NONBONDED_PARM_INDEX', '10I8', num_items=0)
        self.add_flag('RESIDUE_LABEL', '20a4', num_items=0)
        self.add_flag('RESIDUE_POINTER', '10I8', num_items=0)
        self.add_flag('BOND_FORCE_CONSTANT', '5E16.8', num_items=0)
        self.add_flag('BOND_EQUIL_VALUE', '5E16.8', num_items=0)
        self.add_flag('ANGLE_FORCE_CONSTANT', '5E16.8', num_items=0)
        self.add_flag('ANGLE_EQUIL_VALUE', '3E25.17', num_items=0)
        self.add_flag('CHARMM_UREY_BRADLEY_COUNT', '2I8', num_items=2,
                      comments=['V(ub) = K_ub(r_ik - R_ub)**2',
                                'Number of Urey Bradley terms and types'])
        self.add_flag('CHARMM_UREY_BRADLEY', '10I8', num_items=0,
                      comments=['List of the two atoms and its parameter index',
                                'in each UB term: i,k,index'])
        self.add_flag('CHARMM_UREY_BRADLEY_FORCE_CONSTANT', '5E16.8', num_items=0,
                      comments=['K_ub: kcal/mol/A**2'])
        self.add_flag('CHARMM_UREY_BRADLEY_EQUIL_VALUE', '5E16.8', num_items=0,
                      comments=['r_ub: A'])
        self.add_flag('DIHEDRAL_FORCE_CONSTANT', '5E16.8', num_items=0)
        self.add_flag('DIHEDRAL_PERIODICITY', '5E16.8', num_items=0)
        self.add_flag('DIHEDRAL_PHASE', '5E16.8', num_items=0)
        self.add_flag('SCEE_SCALE_FACTOR', '5E16.8', num_items=0)
        self.add_flag('SCNB_SCALE_FACTOR', '5E16.8', num_items=0)
        self.add_flag('CHARMM_NUM_IMPROPERS', '10I8', num_items=0,
                      comments=['Number of terms contributing to the',
                                'quadratic four atom improper energy term:',
                                'V(improper) = K_psi(psi - psi_0)**2'])
        self.add_flag('CHARMM_IMPROPERS', '10I8', num_items=0,
                      comments=['List of the four atoms in each improper term',
                                'i,j,k,l,index  i,j,k,l,index',
                                'where index is into the following two lists:',
                                'CHARMM_IMPROPER_{FORCE_CONSTANT,IMPROPER_PHASE}'])
        self.add_flag('CHARMM_NUM_IMPR_TYPES', '1I8', num_items=1,
                      comments=['Number of unique parameters contributing to the',
                                'quadratic four atom improper energy term'])
        self.add_flag('CHARMM_IMPROPER_FORCE_CONSTANT', '5E16.8', num_items=0,
                      comments=['K_psi: kcal/mole/rad**2'])
        self.add_flag('CHARMM_IMPROPER_PHASE', '5E16.8', num_items=0, comments=['psi: radians'])
        self.pointers['NATYP'] = 1
        self.parm_data['POINTERS'][PrmtopPointers.NATYP] = 1
        self.add_flag('SOLTY', '5E16.8', num_items=1)
        self.add_flag('LENNARD_JONES_ACOEF', '3E24.16', num_items=0)
        self.add_flag('LENNARD_JONES_BCOEF', '3E24.16', num_items=0)
        self.add_flag('LENNARD_JONES_14_ACOEF', '3E24.16', num_items=0)
        self.add_flag('LENNARD_JONES_14_BCOEF', '3E24.16', num_items=0)
        self.add_flag('BONDS_INC_HYDROGEN', '10I8', num_items=0)
        self.add_flag('BONDS_WITHOUT_HYDROGEN', '10I8', num_items=0)
        self.add_flag('ANGLES_INC_HYDROGEN', '10I8', num_items=0)
        self.add_flag('ANGLES_WITHOUT_HYDROGEN', '10I8', num_items=0)
        self.add_flag('DIHEDRALS_INC_HYDROGEN', '10I8', num_items=0)
        self.add_flag('DIHEDRALS_WITHOUT_HYDROGEN', '10I8', num_items=0)
        self.add_flag('EXCLUDED_ATOMS_LIST', '10I8', num_items=0)
        self.add_flag('HBOND_ACOEF', '5E16.8', num_items=0)
        self.add_flag('HBOND_BCOEF', '5E16.8', num_items=0)
        self.add_flag('HBCUT', '5E16.8', num_items=0)
        self.add_flag('AMBER_ATOM_TYPE', '20a4', num_items=0)
        self.add_flag('TREE_CHAIN_CLASSIFICATION', '20a4', num_items=0)
        self.add_flag('JOIN_ARRAY', '10I8', num_items=0)
        self.add_flag('IROTAT', '10I8', num_items=0)
        if self.has_cmap:
            self.add_flag(self._cmap_prefix + 'CMAP_COUNT', '2I8', num_items=2,
                          comments=['Number of CMAP terms, number of unique CMAP parameters'])
            self.add_flag(self._cmap_prefix + 'CMAP_RESOLUTION', '20I4', num_items=0,
                          comments=['Number of steps along each phi/psi CMAP axis',
                                    'for each CMAP_PARAMETER grid'])
            self.add_flag(self._cmap_prefix + 'CMAP_INDEX', '6I8', num_items=0,
                          comments=['Atom index i,j,k,l,m of the cross term',
                                    'and then pointer to CMAP_PARAMETER_n'])
        if self.box is not None:
            self.add_flag('SOLVENT_POINTERS', '3I8', num_items=3)
            self.add_flag('ATOMS_PER_MOLECULE', '10I8', num_items=0)
            self.add_flag('BOX_DIMENSIONS', '5E16.8', num_items=4)
        self.add_flag('RADIUS_SET', '1a80', num_items=1)
        self.add_flag('RADII', '5E16.8', num_items=0)
        self.add_flag('SCREEN', '5E16.8', num_items=0)
        self.add_flag('IPOL', '1I8', num_items=1)

    #===================================================

    def _set_nonbonded_tables(self, nbfixes=None):
        """ Sets the tables of Lennard-Jones nonbonded interaction pairs """
        from ..tools.actions import addLJType
        data = self.parm_data
        ntypes = data['POINTERS'][PrmtopPointers.NTYPES]
        ntypes2 = ntypes * ntypes
        # Set up the index lookup tables (not a unique solution)
        data['NONBONDED_PARM_INDEX'] = [0 for i in range(ntypes2)]
        holder = [0 for i in range(ntypes2)]
        idx = 0
        for i in range(ntypes):
            for j in range(i+1):
                idx += 1
                holder[ntypes*i+j] = holder[ntypes*j+i] = idx
        idx = 0
        for i in range(ntypes):
            for j in range(ntypes):
                data['NONBONDED_PARM_INDEX'][idx] = holder[ntypes*i+j]
                idx += 1
        nttyp = ntypes * (ntypes + 1) // 2
        # Now build the Lennard-Jones arrays
        data['LENNARD_JONES_14_ACOEF'] = [0 for i in range(nttyp)]
        data['LENNARD_JONES_14_BCOEF'] = [0 for i in range(nttyp)]
        data['LENNARD_JONES_ACOEF'] = [0 for i in range(nttyp)]
        data['LENNARD_JONES_BCOEF'] = [0 for i in range(nttyp)]
        self.recalculate_LJ()
        # Now make any NBFIX modifications we had
        if nbfixes is not None:
            for i, fix in enumerate(nbfixes):
                for terms in fix:
                    j, rmin, eps, rmin14, eps14 = terms
                    k, l = min(i, j-1), max(i, j-1)
                    eps = abs(eps)
                    eps14 = abs(eps14)
                    idx = data['NONBONDED_PARM_INDEX'][ntypes*k+l] - 1
                    data['LENNARD_JONES_ACOEF'][idx] = eps * rmin**12
                    data['LENNARD_JONES_BCOEF'][idx] = 2 * eps * rmin**6
                    data['LENNARD_JONES_14_ACOEF'][idx] = eps14 * rmin14**12
                    data['LENNARD_JONES_14_BCOEF'][idx] = 2 * eps14 * rmin14**6
        # If we had an explicit set of exceptions, we need to implement all of
        # those exclusions. The electrostatic component of that exception is
        # already handled (since a scaling factor *can* represent the full
        # flexibility of electrostatic exceptions). The vdW component must be
        # handled as 1-4 A- and B-coefficients. If vdW types are too compressed
        # for this to be handled correctly, use "addLJType" to expand the types
        # by 1 (so the tables only get as big as they *need* to get, to make
        # sure they continue to fit in CUDA shared memory).
        #
        # The way we do this is to fill all of the elements with None, then fill
        # them in as we walk through the 1-4 exceptions. If we hit a case where
        # the target element is not None and *doesn't* equal the computed A- and
        # B-coefficients, we have to expand our type list (assume all type
        # names will have the same L-J parameters, which has been a fair
        # assumption in my experience).
        if not self.adjusts: return
        for i in range(len(data['LENNARD_JONES_14_ACOEF'])):
            data['LENNARD_JONES_14_ACOEF'][i] = None
            data['LENNARD_JONES_14_BCOEF'][i] = None
        ii = 0
        replaced_atoms = set()
        while True:
            needed_split = False
            for pair in self.adjusts:
                a1, a2 = pair.atom1, pair.atom2
                i, j = sorted([a1.nb_idx - 1, a2.nb_idx - 1])
                idx = data['NONBONDED_PARM_INDEX'][ntypes*i+j] - 1
                eps = pair.type.epsilon
                rmin = pair.type.rmin
                rmin6 = rmin * rmin * rmin * rmin * rmin * rmin
                acoef = eps * rmin6*rmin6
                bcoef = 2 * eps * rmin6
                if data['LENNARD_JONES_14_ACOEF'][idx] is not None:
                    if abs(data['LENNARD_JONES_14_ACOEF'][idx] - acoef) > SMALL:
                        # Need to split out another type
                        needed_split = True
                        assert a1 not in replaced_atoms or a2 not in replaced_atoms
                        # Only add each atom as a new type ONCE
                        if a1 in replaced_atoms:
                            mask = '@%d' % (a2.idx+1)
                            replaced_atoms.add(a2)
                        else:
                            mask = '@%d' % (a1.idx+1)
                            replaced_atoms.add(a1)
                        addLJType(self, mask, radius_14=0, epsilon_14=0).execute()
                        ntypes += 1
                        # None-out all of the added terms
                        j = ntypes - 1
                        for i in range(j):
                            idx2 = data['NONBONDED_PARM_INDEX'][ntypes*i+j] - 1
                            data['LENNARD_JONES_14_ACOEF'][idx2] = None
                            data['LENNARD_JONES_14_BCOEF'][idx2] = None
                        # We can stop here, since the next loop through the
                        # explicit exclusions will fill this in
                else:
                    data['LENNARD_JONES_14_ACOEF'][idx] = acoef
                    data['LENNARD_JONES_14_BCOEF'][idx] = bcoef
            ii += 1
            if not needed_split:
                break
            # The following should never happen
            assert ii <= len(self.atoms)+1, 'Could not resolve all exceptions. This is a bug'
        # Now go through and change all None's to 0s, as these terms won't be
        # used for any exceptions, anyway
        for i, item in enumerate(data['LENNARD_JONES_14_ACOEF']):
            if item is None:
                assert data['LENNARD_JONES_14_BCOEF'][i] is None, \
                       'A- and B- coefficients must be in lock-step!'
                data['LENNARD_JONES_14_ACOEF'][i] = 0.0
                data['LENNARD_JONES_14_BCOEF'][i] = 0.0

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def ConvertFromPSF(struct, params, title=''):
    """
    This function instantiates a ChamberParm instance from a data structure
    instantiated by a CHARMM PSF.

    Parameters
    ----------
    struct : Structure
        The structure object (typically loaded from a PSF file)
    params : CharmmParameterSet
        The parameter set describing the parameters of the input struct
    title : str=''
        The title to assign to the topology file

    Returns
    -------
    ChamberParm
        ChamberParm instance with all parameters loaded
    """
    # Make sure all atom types are strings, not integers
    int_starting = str(struct.atoms[0].atom_type) != struct.atoms[0].type
    if int_starting:
        for atom in struct.atoms:
            atom.type = str(atom.atom_type)
    parm = ChamberParm.from_structure(struct)
    parm.parm_data['FORCE_FIELD_TYPE'] = fftype = []
    if params.parametersets == []:
        params.parametersets.append('')
    for pset in params.parametersets:
        if 'CHARMM' not in pset: # needed to trigger "charmm_active"...
            pset = 'CHARMM: %s' % pset
        fftype.extend([len(params.parametersets), pset])

    # Convert atom types back to integers if that's how they started
    if int_starting:
        for atom in struct.atoms:
            atom.type = int(atom.atom_type)
    return parm
