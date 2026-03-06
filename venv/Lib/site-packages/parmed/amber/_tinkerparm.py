"""
This module contains an amber prmtop class that will read in all
parameters and allow users to manipulate that data and write a new
prmtop object.
"""
import numpy as np

from ..constants import DEG_TO_RAD, PrmtopPointers, RAD_TO_DEG
from ..exceptions import AmberError
from ..formats.registry import load_file
from ..topologyobjects import AmoebaNonbondedExceptionType as NonbondedExceptionType
from ..topologyobjects import (Angle, AngleType, Bond, BondType, ChiralFrame,
                               Dihedral, DihedralType, MultipoleFrame,
                               NonbondedException, OutOfPlaneBend,
                               OutOfPlaneBendType, PiTorsion, StretchBend,
                               StretchBendType, TorsionTorsion,
                               TorsionTorsionType, TrigonalAngle, UreyBradley)
from ._amberparm import AmberParm
from .amberformat import AmberFormat


class AmoebaParm(AmberParm):
    """
    Tinker Topology (parm7 format) class defining the AMOEBA force field. Gives
    low, and some high, level access to topology data.  You can interact with
    the raw data in the topology file directly or interact with some of the
    high-level classes comprising the system topology and parameters.

    Parameters
    ----------
    prm_name : str=None
        If provided, this file is parsed and the data structures will be loaded
        from the data in this file
    rst7_name : str=None
        If provided, the coordinates and unit cell dimensions from the provided
        Amber inpcrd/restart file will be loaded into the molecule

    Attributes
    ----------
    parm_data : dict {str : list}
        A dictionary that maps FLAG names to all of the data contained in that
        section of the Amber file.
    formats : dict {str : FortranFormat}
        A dictionary that maps FLAG names to the FortranFormat instance in which
        the data is stored in that section
    parm_comments : dict {str : list}
        A dictionary that maps FLAG names to the list of COMMENT lines that were
        stored in the original file
    flag_list : list
        An ordered list of all FLAG names. This must be kept synchronized with
        `parm_data`, `formats`, and `parm_comments` such that every item in
        `flag_list` is a key to those 3 dicts and no other keys exist
    version : str
        The VERSION string from the Amber file
    name : str
        The file name of the originally parsed file (set to the fname parameter)
    atoms : AtomList(Atom)
        List of all atoms in the system
    residues : ResidueList(Residue)
        List of all residues in the system
    bonds : TrackedList(Bond)
        List of bonds between two atoms in the system
    angles : TrackedList(Angle)
        List of regular angles between three atoms in the system
    dihedrals : TrackedList(Angle)
        List of all proper torsions between 4 atoms in the system
    urey_bradleys : TrackedList(UreyBradley)
        List of all Urey-Bradley terms between 2 atoms connected by an angle
    trigonal_angles : TrackedList(TrigonalAngle)
        List of all trigonal angle terms
    out_of_plane_bends : TrackedList(OutOfPlaneBend)
        List of all out-of-plane bending terms
    pi_torsions : TrackedList(PiTorsion)
        List of all pi-torsion terms
    stretch_bends : TrackedList(StretchBend)
        List of all stretch-bending terms
    torsion_torsions : TrackedList(TorsionTorsion)
        List of all coupled torsion-torsion terms
    chiral_frames : TrackedList(ChiralFrame)
        List of all chiral centers
    multipole_frames : TrackedList(MultipoleFrame)
        List of all multipole frames of reference
    adjusts : TrackedList(NonbondedException)
        List of all nonbonded exception parameters used for adjusting nonbonded
        interactions between particular pairs of atoms
    box : list of 6 floats
        Periodic boundary unit cell dimensions and angles
    bond_types : TrackedList(BondType)
        The bond types containing the parameters for each bond stretching term
    angle_types : TrackedList(AngleType)
        The angle types containing the parameters for each angle bending term
    dihedral_types : TrackedList(DihedralType)
        The dihedral types containing the parameters for each torsional term
    urey_bradley_types : TrackedList(BondType)
        The Urey-Bradley types containing the parameters for each term
    trigonal_angle_types : TrackedList(AngleType)
        The trigonal angle types containing the parameters for each term
    out_of_plane_bend_types : TrackedList(AngleType)
        The out-of-plane bending angle type containing parameters for each term
    pi_torsion_types : TrackedList(DihedralType)
        The pi-torsion type containing parameters for each torsional term
    stretch_bend_types : TrackedList(StretchBendType)
        The stretch-bend type containing parameters for each term
    torsion_torsion_types : TrackedList(TorsionTorsionType)
        The coupled torsion-torsion type containing parameters for coupled
        torsions
    adjust_types : TrackedList(NonbondedExceptionType)
        The nonbonded exception scaling factors for pairs of particles
    bonds_inc_h : iterator(Bond)
        Read-only generator that loops through all bonds that contain Hydrogen
    bonds_without_h : iterator(Bond)
        Read-only generator that loops through all bonds that do not contain
        Hydrogen
    angles_inc_h : iterator(Angle)
        Read-only generator that loops through all angles that contain Hydrogen
    angles_without_h : iterator(Angle)
        Read-only generator that loops through all angles that do not contain
        Hydrogen
    dihedrals_inc_h : iterator(Dihedral)
        Read-only generator that loops through all dihedrals that contain
        Hydrogen
    dihedrals_without_h : iterator(Dihedral)
        Read-only generator that loops through all dihedrals that do not contain
        Hydrogen
    chamber : bool=False
        On AmberParm instances, this is always False to indicate that it is not
        a CHAMBER-style topology file
    amoeba : bool=False
        On AmberParm instances, this is always False to indicate that it is not
        an AMOEBA-style topology file
    has_cmap : bool=False
        On AmberParm instances, this is always False to indicate that it does
        not have correction maps (unique to CHARMM force field and chamber
        topologies)
    """

    #=============================================

    def initialize_topology(self, xyz=None, box=None):
        """
        Initializes topology data structures, like the list of atoms, bonds,
        etc., after the topology file has been read.

        Raises
        ------
        AmberError if it is not an Amoeba-styled topology file
        """
        try:
            if self.parm_data['AMOEBA_FORCEFIELD'][0] != 1:
                raise AmberError('Bad AMOEBA-format topology')
        except KeyError:
            raise AmberError('Bad AMOEBA-format topology')

        # We need to handle RESIDUE_ICODE properly since it may have picked up
        # some extra values
        if 'RESIDUE_ICODE' in self.parm_data:
            self._truncate_array('RESIDUE_ICODE', self.parm_data['POINTERS'][PrmtopPointers.NRES])

        self.LJ_types = {}
        self.LJ_radius = []
        self.LJ_depth = []
        # If we were given a prmtop, read it in
        self.pointers = {}
        self.load_pointers()

        # Load the structure arrays
        self.load_structure()

        if isinstance(xyz, str):
            f = load_file(xyz)
            if not hasattr(f, 'coordinates') or f.coordinates is None:
                raise TypeError('%s does not have coordinates' % xyz)
            self.coordinates = f.coordinates
            if hasattr(f, 'box') and f.box is not None and box is None:
                self.box = f.box
        else:
            self.coordinates = xyz
        if box is not None:
            self.box = box

        # If all else fails, set the box from the prmtop file
        if self.parm_data['POINTERS'][PrmtopPointers.IFBOX] > 0 and self.box is None:
            box = self.parm_data['BOX_DIMENSIONS']
            self.box = list(box[1:]) + [box[0], box[0], box[0]]

        self.hasbox = self.box is not None

    #=============================================

    def load_structure(self):
        """
        Loads all of the topology instance variables. This is necessary if we
        actually want to modify the topological layout of our system
        """
        self._load_atoms_and_residues()
        self.load_atom_info()
        self._load_bond_info()
        self._load_angle_info()
        self._load_urey_bradley_info()
        self._load_trigonal_angle_info()
        self._load_oopbend_info()
        self._load_dihedral_info()
        self._load_pitorsion_info()
        self._load_stretch_bend_info()
        self._load_torsion_torsion_info()
        self._load_frame_info()
        self._load_exception_info()
        super(AmoebaParm, self).unchange()

    #=============================================

    def load_atom_info(self):
        """
        Loads atom properties into the atoms that have been loaded. If any
        arrays are too short or too long, an IndexError will be raised
        """
        anam = self.parm_data['ATOM_NAME']
        mass = self.parm_data['MASS']
        atyp = self.parm_data['AMBER_ATOM_TYPE']
        tree = self.parm_data['TREE_CHAIN_CLASSIFICATION']
        join = self.parm_data['JOIN_ARRAY']
        irot = self.parm_data['IROTAT']
        typi = self.parm_data['AMOEBA_ATOM_TYPE_INDEX']
        atnum = self.parm_data['AMOEBA_ATOMIC_NUMBER']
        clsi = self.parm_data['AMOEBA_ATOM_CLASS_INDEX']
        nbtyp = self.parm_data['AMOEBA_VDW_ATOM_TYPES_LIST']
        vdwp = self.parm_data['AMOEBA_VDW_ATOM_PARENT_LIST']
        vdww = self.parm_data['AMOEBA_VDW_PARENT_COORD_WEIGHT_LIST']
        mpole = self.parm_data['AMOEBA_LOCAL_FRAME_MULTIPOLES_LIST']
        pol = self.parm_data['AMOEBA_POLARIZABILITY_LIST']
        for i, atom in enumerate(self.atoms):
            i10 = i * 10
            multipoles = mpole[i10:i10+10]
            atom.name = anam[i]
            atom.type = atyp[i]
            atom.mass = mass[i]
            atom.tree = tree[i]
            atom.join = join[i]
            atom.irotat = irot[i]
            atom.type_idx = typi[i]
            atom.atomic_number = atnum[i]
            atom.class_idx = clsi[i]
            atom.vdw_parent = self.atoms[vdwp[i]-1]
            atom.vdw_weight = vdww[i]
            atom.charge = multipoles[0] # monopole is charge
            atom.multipoles = multipoles
            atom.polarizability = pol[i]
            atom.nb_idx = nbtyp[i]

    #=============================================

    def remake_parm(self):
        """ Recomputes the topology file parameters and fills parm_data """
        # Get rid of terms containing deleted atoms and empty residues
        self.prune_empty_terms()
        self.residues.prune()
        self.rediscover_molecules()

        # Transfer information from the topology lists
        self._xfer_atom_info()
        self._xfer_residue_info()
        self._xfer_bond_info()
        self._xfer_angle_info()
        self._xfer_urey_bradley_info()
        self._xfer_trigonal_angle_info()
        self._xfer_oopbend_info()
        self._xfer_dihedral_info()
        self._xfer_pitorsion_info()
        self._xfer_stretch_bend_info()
        self._xfer_torsion_torsion_info()
        self._xfer_frame_info()
        self._xfer_exception_info()

    #=============================================

    def mdin_skeleton(self):
        """
        Returns the skeleton of an mdin file with the &amoeba namelist set up
        correctly for the potential terms that are present in this topology
        file.

        Returns
        -------
        str
            A skeleton MDIN file with all of the do_* variables in the &amoeba
            section set correctly. It is commented for easy editing
        """
        return ('Input file for AMOEBA simulations.\n'
                ' &cntrl\n'
                '     ! Add whatever variables you need here\n'
                '     ntb=1, ntt=1, ntp=0, ! PBC, thermostat, barostat\n'
                '     irest=0, ntx=1,      ! restart flags\n'
                ' /\n'
                ' &amoeba\n'
                '     ! Some basic potential parameters. For better\n'
                '     ! energy conservation you need to adjust these\n'
                '     ! defaults\n'
                '     beeman_integrator=1,   ! Use Beeman integrator\n'
                '     dipole_scf_tol=0.01,   ! 10e-6 gives good NVE\n'
                '\n'
                '     ! You should not generally modify these variables:\n'
                '     do_valence=1, do_bond=%d, do_ureyb=%d,\n'
                '     do_reg_angle=%d, do_trig_angle=%d, do_opbend=%d,\n'
                '     do_torsion=%d, do_pi_torsion=%d, do_strbend=%d,\n'
                '     do_torsion_torsion=%d,\n'
                ' /\n' % (bool(self.bonds), bool(self.urey_bradleys),
                bool(self.angles), bool(self.trigonal_angles),
                bool(self.out_of_plane_bends), bool(self.dihedrals),
                bool(self.pi_torsions), bool(self.stretch_bends),
                bool(self.torsion_torsions))
        )

    #=============================================

    @property
    def chamber(self):
        return False
   
    @property
    def amoeba(self):
        return True

    #===========  PRIVATE INSTANCE METHODS  ============

    def _load_bond_info(self):
        """ Load the regular AMOEBA bonds, if they exist """
        if not 'AMOEBA_REGULAR_BOND_LIST' in self.parm_data: return
        data = self.parm_data
        del self.bonds[:]
        del self.bond_types[:]
        for k, req in zip(data['AMOEBA_REGULAR_BOND_FORCE_CONSTANT'],
                          data['AMOEBA_REGULAR_BOND_EQUIL_VALUE']):
            self.bond_types.append(BondType(k, req, self.bond_types))
        self.bond_types.degree = data['AMOEBA_REGULAR_BOND_FTAB_DEGREE'][0]
        self.bond_types.coeffs = data['AMOEBA_REGULAR_BOND_FTAB_COEFFS'][:]
        if len(self.bond_types.coeffs) != self.bond_types.degree + 1:
            raise AmberError('Bond degree (%d) does not make sense with %d '
                             'coefficients' % (self.bond_types.degree,
                             len(self.bond_types.coeffs)))
        it = iter(self.parm_data['AMOEBA_REGULAR_BOND_LIST'])
        for i, j, k in zip(it, it, it):
            self.bonds.append(
                    Bond(self.atoms[i-1], self.atoms[j-1],
                         self.bond_types[k-1])
            )

    #=============================================

    def _load_angle_info(self):
        """ Load the regular AMOEBA angles, if they exist """
        if not 'AMOEBA_REGULAR_ANGLE_LIST' in self.parm_data: return
        data = self.parm_data
        del self.angles[:]
        del self.angle_types[:]
        for k, eq in zip(data['AMOEBA_REGULAR_ANGLE_FORCE_CONSTANT'],
                         data['AMOEBA_REGULAR_ANGLE_EQUIL_VALUE']):
            self.angle_types.append(AngleType(k, eq, self.angle_types))
        self.angle_types.degree = data['AMOEBA_REGULAR_ANGLE_FTAB_DEGREE'][0]
        self.angle_types.coeffs = data['AMOEBA_REGULAR_ANGLE_FTAB_COEFFS'][:]
        if len(self.angle_types.coeffs) != self.angle_types.degree + 1:
            raise AmberError('Angle degree (%d) does not make sense with %d '
                              'coefficients' % (self.angle_types.degree,
                              len(self.angle_types.coeffs)))
        it = iter(data['AMOEBA_REGULAR_ANGLE_LIST'])
        for i, j, k, l in zip(it, it, it, it):
            self.angles.append(
                    Angle(self.atoms[i-1], self.atoms[j-1], self.atoms[k-1],
                          self.angle_types[l-1])
            )

    #=============================================

    def _load_urey_bradley_info(self):
        """ Loads the AMOEBA Urey-Bradley terms, if they exist """
        if not 'AMOEBA_UREY_BRADLEY_BOND_LIST' in self.parm_data: return
        data = self.parm_data
        del self.urey_bradleys[:]
        del self.urey_bradley_types[:]
        for k, eq in zip(data['AMOEBA_UREY_BRADLEY_BOND_FORCE_CONSTANT'],
                         data['AMOEBA_UREY_BRADLEY_BOND_EQUIL_VALUE']):
            self.urey_bradley_types.append(
                    BondType(k, eq, self.urey_bradley_types)
            )
        self.urey_bradley_types.degree = degree = \
                data['AMOEBA_UREY_BRADLEY_BOND_FTAB_DEGREE'][0]
        self.urey_bradley_types.coeffs = coeffs = \
                data['AMOEBA_UREY_BRADLEY_BOND_FTAB_COEFFS'][:]
        if len(coeffs) != degree + 1:
            raise AmberError('Urey-Bradley degree (%d) does not make sense '
                              'with %d coefficients' % (degree, len(coeffs)))
        it = iter(data['AMOEBA_UREY_BRADLEY_BOND_LIST'])
        for i, j, k in zip(it, it, it):
            self.urey_bradleys.append(
                    UreyBradley(self.atoms[i-1], self.atoms[j-1],
                                self.urey_bradley_types[k-1])
            )

    #=============================================

    def _load_trigonal_angle_info(self):
        """ Loads the AMOEBA trigonal angle terms, if they exist """
        if not 'AMOEBA_TRIGONAL_ANGLE_LIST' in self.parm_data: return
        data = self.parm_data
        del self.trigonal_angles[:]
        del self.trigonal_angle_types[:]
        for k, eq in zip(data['AMOEBA_TRIGONAL_ANGLE_FORCE_CONSTANT'],
                         data['AMOEBA_TRIGONAL_ANGLE_EQUIL_VALUE']):
            self.trigonal_angle_types.append(
                    AngleType(k, eq, self.trigonal_angle_types)
            )
        self.trigonal_angle_types.degree = degree = \
                data['AMOEBA_TRIGONAL_ANGLE_FTAB_DEGREE'][0]
        self.trigonal_angle_types.coeffs = coeffs = \
                data['AMOEBA_TRIGONAL_ANGLE_FTAB_COEFFS'][:]
        if len(coeffs) != degree + 1:
            raise AmberError('Trigonal Angle degree (%d) does not make sense '
                              'with %d coefficients' % (degree, len(coeffs)))
        it = iter(data['AMOEBA_TRIGONAL_ANGLE_LIST'])
        for i, j, k, l, m in zip(it, it, it, it, it):
            self.trigonal_angles.append(
                    TrigonalAngle(self.atoms[i-1], self.atoms[j-1],
                                  self.atoms[k-1], self.atoms[l-1],
                                  self.trigonal_angle_types[m-1])
            )

    #=============================================

    def _load_oopbend_info(self):
        """ Loads the AMOEBA out-of-plane bending terms, if they exist """
        if not 'AMOEBA_OPBEND_ANGLE_LIST' in self.parm_data: return
        data = self.parm_data
        del self.out_of_plane_bends[:]
        del self.out_of_plane_bend_types[:]
        for k in data['AMOEBA_OPBEND_ANGLE_FORCE_CONSTANT']:
            self.out_of_plane_bend_types.append(
                    OutOfPlaneBendType(k, self.out_of_plane_bend_types)
            )
        self.out_of_plane_bend_types.degree = degree = \
                data['AMOEBA_OPBEND_ANGLE_FTAB_DEGREE'][0]
        self.out_of_plane_bend_types.coeffs = coeffs = \
                data['AMOEBA_OPBEND_ANGLE_FTAB_COEFFS'][:]
        if len(coeffs) != degree + 1:
            raise AmberError('OOP-bend angle degree (%d) does not make sense '
                              'with %d coefficients.' % (degree, len(coeffs)))
        it = iter(data['AMOEBA_OPBEND_ANGLE_LIST'])
        for i, j, k, l, m in zip(it, it, it, it, it):
            self.out_of_plane_bends.append(
                    OutOfPlaneBend(self.atoms[i-1], self.atoms[j-1],
                                   self.atoms[k-1], self.atoms[l-1],
                                   self.out_of_plane_bend_types[m-1])
            )

    #=============================================

    def _load_dihedral_info(self):
        """ Loads the AMOEBA regular torsion terms, if they exist """
        if not 'AMOEBA_TORSION_LIST' in self.parm_data: return
        data = self.parm_data
        del self.dihedrals[:]
        del self.dihedral_types[:]
        for k, per, phase in zip(data['AMOEBA_TORSION_FORCE_CONSTANT'],
                                 data['AMOEBA_TORSION_PERIODICITY'],
                                 data['AMOEBA_TORSION_PHASE']):
            self.dihedral_types.append(
                    DihedralType(k, per, phase*RAD_TO_DEG,
                                 list=self.dihedral_types)
            )
        it = iter(data['AMOEBA_TORSION_LIST'])
        for i, j, k, l, m in zip(it, it, it, it, it):
            self.dihedrals.append(
                    Dihedral(self.atoms[i-1], self.atoms[j-1], self.atoms[k-1],
                             self.atoms[l-1], type=self.dihedral_types[m-1])
            )

    #=============================================

    def _load_pitorsion_info(self):
        """ Loads the AMOEBA pi-torsion terms, if they exist """
        if not 'AMOEBA_PI_TORSION_LIST' in self.parm_data: return
        data = self.parm_data
        del self.pi_torsions[:]
        del self.pi_torsion_types[:]
        for k, per, phase in zip(data['AMOEBA_PI_TORSION_FORCE_CONSTANT'],
                                 data['AMOEBA_PI_TORSION_PERIODICITY'],
                                 data['AMOEBA_PI_TORSION_PHASE']):
            self.pi_torsion_types.append(
                    DihedralType(k, per, phase*RAD_TO_DEG,
                                 list=self.pi_torsion_types)
            )
        it = iter(data['AMOEBA_PI_TORSION_LIST'])
        for i, j, k, l, m, n, o in zip(it, it, it, it, it, it, it):
            self.pi_torsions.append(
                    PiTorsion(self.atoms[i-1], self.atoms[j-1], self.atoms[k-1],
                              self.atoms[l-1], self.atoms[m-1], self.atoms[n-1],
                              self.pi_torsion_types[o-1])
            )

    #=============================================

    def _load_stretch_bend_info(self):
        """ Loads the AMOEBA stretch-bend terms, if they exist """
        if not 'AMOEBA_STRETCH_BEND_LIST' in self.parm_data: return
        data = self.parm_data
        del self.stretch_bends[:]
        del self.stretch_bend_types[:]
        if 'AMOEBA_STRETCH_BEND_FORCE_CONSTANT' in data:
            for a,b,c,d in zip(data['AMOEBA_STRETCH_BEND_FORCE_CONSTANT'],
                               data['AMOEBA_STRETCH_BEND_BOND1_EQUIL_VALUE'],
                               data['AMOEBA_STRETCH_BEND_BOND2_EQUIL_VALUE'],
                               data['AMOEBA_STRETCH_BEND_ANGLE_EQUIL_VALUE']):
                self.stretch_bend_types.append(
                        StretchBendType(a, a, b, c, d,
                                list=self.stretch_bend_types)
                )
        elif 'AMOEBA_STRETCH_BEND_FORCE_CONSTANT_1' in data:
            for a,b,c,d,e in zip(data['AMOEBA_STRETCH_BEND_FORCE_CONSTANT_1'],
                                 data['AMOEBA_STRETCH_BEND_FORCE_CONSTANT_2'],
                                 data['AMOEBA_STRETCH_BEND_BOND1_EQUIL_VALUE'],
                                 data['AMOEBA_STRETCH_BEND_BOND2_EQUIL_VALUE'],
                                 data['AMOEBA_STRETCH_BEND_ANGLE_EQUIL_VALUE']):
                self.stretch_bend_types.append(
                        StretchBendType(a, b, c, d, e,
                                list=self.stretch_bend_types)
                )
        it = iter(data['AMOEBA_STRETCH_BEND_LIST'])
        for i, j, k, l in zip(it, it, it, it):
            self.stretch_bends.append(
                    StretchBend(self.atoms[i-1], self.atoms[j-1],
                                self.atoms[k-1], self.stretch_bend_types[l-1])
            )

    #=============================================

    def _load_torsion_torsion_info(self):
        """ Loads the AMOEBA coupled torsion-torsion terms, if they exist """
        if not 'AMOEBA_TORSION_TORSION_LIST' in self.parm_data: return
        del self.torsion_torsion_types[:]
        del self.torsion_torsions[:]
        data = self.parm_data
        ntypes = data['AMOEBA_TORSION_TORSION_NUM_PARAMS'][0]
        for i in range(ntypes):
            prefix = 'AMOEBA_TORSION_TORSION_TORTOR_TABLE_%02d_' % (i + 1)
            dims = tuple(data[prefix + 'DIMS'])
            ang1 = data[prefix + 'ANGLE1']
            ang2 = data[prefix + 'ANGLE2']
            f = data[prefix + 'FUNC']
            dfda1 = data[prefix + 'DFUNC_DANGLE1']
            dfda2 = data[prefix + 'DFUNC_DANGLE2']
            d2fda1da2 = data[prefix + 'D2FUNC_DANGLE1_DANGLE2']
            self.torsion_torsion_types.append(
                    TorsionTorsionType(dims, ang1, ang2, f, dfda1,
                                       dfda2, d2fda1da2,
                                       list=self.torsion_torsion_types)
            )
        it = iter(data['AMOEBA_TORSION_TORSION_LIST'])
        for i, j, k, l, m, n in zip(it, it, it, it, it, it):
            self.torsion_torsions.append(
                TorsionTorsion(
                    self.atoms[i-1],
                    self.atoms[j-1],
                    self.atoms[k-1],
                    self.atoms[l-1],
                    self.atoms[m-1],
                    self.torsion_torsion_types[n-1],
                )
            )

    #=============================================

    def _load_frame_info(self):
        """ Loads the AMOEBA chiral and multipole frames """
        data = self.parm_data
        del self.chiral_frames[:]
        if 'AMOEBA_CHIRAL_FRAME_LIST' in data:
            it = iter(data['AMOEBA_CHIRAL_FRAME_LIST'])
            for i, j, k in zip(it, it, it):
                self.chiral_frames.append(
                    ChiralFrame(self.atoms[i-1], self.atoms[j-1], k)
                )
        del self.multipole_frames[:]
        if 'AMOEBA_FRAME_DEF_LIST' in data:
            it = iter(data['AMOEBA_FRAME_DEF_LIST'])
            for i, j, k, l, m in zip(it, it, it, it, it):
                self.multipole_frames.append(
                    MultipoleFrame(self.atoms[i-1], j, k, l, m)
                )

    #=============================================

    def _load_exception_info(self):
        """ Loads all of the pairwise nonbonded exception rules """
        del self.adjust_types[:]
        del self.adjusts[:]
        data = self.parm_data
        # This section should always be present
        for a,b,c,d,e in zip(data['AMOEBA_ADJUST_VDW_WEIGHTS_LIST'],
                             data['AMOEBA_ADJUST_MPOLE_WEIGHTS_LIST'],
                             data['AMOEBA_ADJUST_DIRECT_WEIGHTS_LIST'],
                             data['AMOEBA_ADJUST_POLAR_WEIGHTS_LIST'],
                             data['AMOEBA_ADJUST_MUTUAL_WEIGHTS_LIST']):
            self.adjust_types.append(
                NonbondedExceptionType(a,b,c,d,e,list=self.adjust_types)
            )
        it = iter(data['AMOEBA_ADJUST_LIST'])
        for i, j, k in zip(it, it, it):
            self.adjusts.append(
                NonbondedException(self.atoms[i-1], self.atoms[j-1], self.adjust_types[k-1])
            )

    #=============================================

    def _xfer_atom_info(self):
        """ Transfers atom info to the topology file data arrays """
        data = self.parm_data
        data['POINTERS'][PrmtopPointers.NATOM] = len(self.atoms)
        self.pointers['NATOM'] = len(self.atoms)
        data['ATOM_NAME'] = [a.name for a in self.atoms]
        data['MASS'] = [a.mass for a in self.atoms]
        data['CHARGE'] = [0.0 for a in self.atoms] # charge is in multipoles
        data['AMBER_ATOM_TYPE'] = [a.type for a in self.atoms]
        data['TREE_CHAIN_CLASSIFICATION'] = [a.tree for a in self.atoms]
        data['JOIN_ARRAY'] = [a.join for a in self.atoms]
        data['IROTAT'] = [a.irotat for a in self.atoms]
        data['AMOEBA_ATOM_TYPE_INDEX'] = [a.type_idx for a in self.atoms]
        data['AMOEBA_ATOMIC_NUMBER'] = [a.atomic_number for a in self.atoms]
        data['AMOEBA_ATOM_CLASS_INDEX'] = [a.class_idx for a in self.atoms]
        data['AMOEBA_VDW_ATOM_TYPES_LIST'] = [a.nb_idx for a in self.atoms]
        data['AMOEBA_VDW_ATOM_PARENT_LIST'] = [a.vdw_parent.idx+1 for a in self.atoms]
        data['AMOEBA_VDW_PARENT_COORD_WEIGHT_LIST'] = [a.vdw_weight for a in self.atoms]
        data['AMOEBA_POLARIZABILITY_LIST'] = [a.polarizability for a in self.atoms]
        data['AMOEBA_LOCAL_FRAME_MULTIPOLES_LIST'] = mpoles = []
        for atom in self.atoms:
            mpoles.extend(atom.multipoles)

    #=============================================

    def _xfer_bond_info(self):
        """
        Transfers the bond information from the bond arrays to the raw data
        arrays
        """
        if len(self.bonds) == 0:
            self.delete_flag('AMOEBA_REGULAR_BOND_NUM_PARAMS')
            self.delete_flag('AMOEBA_REGULAR_BOND_FORCE_CONSTANT')
            self.delete_flag('AMOEBA_REGULAR_BOND_EQUIL_VALUE')
            self.delete_flag('AMOEBA_REGULAR_BOND_FTAB_DEGREE')
            self.delete_flag('AMOEBA_REGULAR_BOND_FTAB_COEFFS')
            self.delete_flag('AMOEBA_REGULAR_BOND_NUM_LIST')
            self.delete_flag('AMOEBA_REGULAR_BOND_LIST')
            return
        data = self.parm_data
        for bond_type in self.bond_types:
            bond_type.used = False
        for bond in self.bonds:
            bond.type.used = True
        self.bond_types.prune_unused()
        data['AMOEBA_REGULAR_BOND_NUM_PARAMS'] = [len(self.bond_types)]
        data['AMOEBA_REGULAR_BOND_FORCE_CONSTANT'] = [bt.k for bt in self.bond_types]
        data['AMOEBA_REGULAR_BOND_EQUIL_VALUE'] = [bt.req for bt in self.bond_types]
        data['AMOEBA_REGULAR_BOND_FTAB_DEGREE'] = [self.bond_types.degree]
        data['AMOEBA_REGULAR_BOND_FTAB_COEFFS'] = self.bond_types.coeffs[:]
        data['AMOEBA_REGULAR_BOND_NUM_LIST'] = [len(self.bonds)]
        data['AMOEBA_REGULAR_BOND_LIST'] = bond_array = []
        for bond in self.bonds:
            bond_array.extend([bond.atom1.idx+1, bond.atom2.idx+1, bond.type.idx+1])

    #=============================================

    def _xfer_angle_info(self):
        """
        Transfers the regular AMOEBA angle information from the topology arrays
        to the raw data arrays
        """
        if len(self.angles) == 0:
            self.delete_flag('AMOEBA_REGULAR_ANGLE_NUM_PARAMS')
            self.delete_flag('AMOEBA_REGULAR_ANGLE_FORCE_CONSTANT')
            self.delete_flag('AMOEBA_REGULAR_ANGLE_EQUIL_VALUE')
            self.delete_flag('AMOEBA_REGULAR_ANGLE_FTAB_DEGREE')
            self.delete_flag('AMOEBA_REGULAR_ANGLE_FTAB_COEFFS')
            self.delete_flag('AMOEBA_REGULAR_ANGLE_NUM_LIST')
            self.delete_flag('AMOEBA_REGULAR_ANGLE_LIST')
            return
        data = self.parm_data
        for angle_type in self.angle_types:
            angle_type.used = False
        for angle in self.angles:
            angle.type.used = True
        self.angle_types.prune_unused()
        data['AMOEBA_REGULAR_ANGLE_NUM_PARAMS'] = [len(self.angle_types)]
        data['AMOEBA_REGULAR_ANGLE_FORCE_CONSTANT'] = [at.k for at in self.angle_types]
        data['AMOEBA_REGULAR_ANGLE_EQUIL_VALUE'] = [at.theteq for at in self.angle_types]
        data['AMOEBA_REGULAR_ANGLE_FTAB_DEGREE'] = [self.angle_types.degree]
        data['AMOEBA_REGULAR_ANGLE_FTAB_COEFFS'] = self.angle_types.coeffs[:]
        data['AMOEBA_REGULAR_ANGLE_NUM_LIST'] = [len(self.angles)]
        data['AMOEBA_REGULAR_ANGLE_LIST'] = angle_array = []
        for angle in self.angles:
            angle_array.extend([angle.atom1.idx+1, angle.atom2.idx+1,
                                angle.atom3.idx+1, angle.type.idx+1])

    #=============================================

    def _xfer_urey_bradley_info(self):
        """
        Transfers the AMOEBA Urey-Bradley bond information from the topology
        arrays to the raw data arrays
        """
        if len(self.urey_bradleys) == 0:
            self.delete_flag('AMOEBA_UREY_BRADLEY_BOND_NUM_PARAMS')
            self.delete_flag('AMOEBA_UREY_BRADLEY_BOND_FORCE_CONSTANT')
            self.delete_flag('AMOEBA_UREY_BRADLEY_BOND_EQUIL_VALUE')
            self.delete_flag('AMOEBA_UREY_BRADLEY_BOND_FTAB_DEGREE')
            self.delete_flag('AMOEBA_UREY_BRADLEY_BOND_FTAB_COEFFS')
            self.delete_flag('AMOEBA_UREY_BRADLEY_BOND_NUM_LIST')
            self.delete_flag('AMOEBA_UREY_BRADLEY_BOND_LIST')
            return
        data = self.parm_data
        for urey_bradley_type in self.urey_bradley_types:
            urey_bradley_type.used = False
        for urey_bradley in self.urey_bradleys:
            urey_bradley.type.used = True
        self.urey_bradley_types.prune_unused()
        data['AMOEBA_UREY_BRADLEY_BOND_NUM_PARAMS'] = [len(self.urey_bradley_types)]
        data['AMOEBA_UREY_BRADLEY_BOND_FORCE_CONSTANT'] = [ut.k for ut in self.urey_bradley_types]
        data['AMOEBA_UREY_BRADLEY_BOND_EQUIL_VALUE'] = [ut.req for ut in self.urey_bradley_types]
        data['AMOEBA_UREY_BRADLEY_BOND_FTAB_DEGREE'] = [self.urey_bradley_types.degree]
        data['AMOEBA_UREY_BRADLEY_BOND_FTAB_COEFFS'] = self.urey_bradley_types.coeffs[:]
        data['AMOEBA_UREY_BRADLEY_BOND_NUM_LIST'] = [len(self.urey_bradleys)]
        data['AMOEBA_UREY_BRADLEY_BOND_LIST'] = urey_array = []
        for urey in self.urey_bradleys:
            urey_array.extend([urey.atom1.idx+1, urey.atom2.idx+1,
                               urey.type.idx+1])

    #=============================================

    def _xfer_trigonal_angle_info(self):
        """
        Transfers the AMOEBA trigonal angle information from the topology arrays
        to the raw data arrays
        """
        if len(self.trigonal_angles) == 0:
            self.delete_flag('AMOEBA_TRIGONAL_ANGLE_NUM_PARAMS')
            self.delete_flag('AMOEBA_TRIGONAL_ANGLE_FORCE_CONSTANT')
            self.delete_flag('AMOEBA_TRIGONAL_ANGLE_EQUIL_VALUE')
            self.delete_flag('AMOEBA_TRIGONAL_ANGLE_FTAB_DEGREE')
            self.delete_flag('AMOEBA_TRIGONAL_ANGLE_FTAB_COEFFS')
            self.delete_flag('AMOEBA_TRIGONAL_ANGLE_NUM_LIST')
            self.delete_flag('AMOEBA_TRIGONAL_ANGLE_LIST')
            return
        data = self.parm_data
        for trigonal_angle_type in self.trigonal_angle_types:
            trigonal_angle_type.used = False
        for trigonal_angle in self.trigonal_angles:
            trigonal_angle.type.used = True
        self.trigonal_angle_types.prune_unused()
        data['AMOEBA_TRIGONAL_ANGLE_NUM_PARAMS'] = [len(self.trigonal_angle_types)]
        data['AMOEBA_TRIGONAL_ANGLE_FORCE_CONSTANT'] = [at.k for at in self.trigonal_angle_types]
        data['AMOEBA_TRIGONAL_ANGLE_EQUIL_VALUE'] = [at.theteq for at in self.trigonal_angle_types]
        data['AMOEBA_TRIGONAL_ANGLE_FTAB_DEGREE'] = [self.trigonal_angle_types.degree]
        data['AMOEBA_TRIGONAL_ANGLE_FTAB_COEFFS'] = self.trigonal_angle_types.coeffs[:]
        data['AMOEBA_TRIGONAL_ANGLE_NUM_LIST'] = [len(self.trigonal_angles)]
        data['AMOEBA_TRIGONAL_ANGLE_LIST'] = angle_array = []
        for angle in self.trigonal_angles:
            angle_array.extend([angle.atom1.idx+1, angle.atom2.idx+1,
                                angle.atom3.idx+1, angle.atom4.idx+1,
                                angle.type.idx+1])

    #=============================================

    def _xfer_oopbend_info(self):
        """
        Transfers the AMOEBA out-of-plane bending angle information from the
        topology arrays to the raw data arrays
        """
        if len(self.out_of_plane_bends) == 0:
            self.delete_flag('AMOEBA_OPBEND_ANGLE_NUM_PARAMS')
            self.delete_flag('AMOEBA_OPBEND_ANGLE_FORCE_CONSTANT')
            self.delete_flag('AMOEBA_OPBEND_ANGLE_FTAB_DEGREE')
            self.delete_flag('AMOEBA_OPBEND_ANGLE_FTAB_COEFFS')
            self.delete_flag('AMOEBA_OPBEND_ANGLE_NUM_LIST')
            self.delete_flag('AMOEBA_OPBEND_ANGLE_LIST')
            return
        data = self.parm_data
        for out_of_plane_bend_type in self.out_of_plane_bend_types:
            out_of_plane_bend_type.used = False
        for out_of_plane_bend in self.out_of_plane_bends:
            out_of_plane_bend.type.used = True
        self.out_of_plane_bend_types.prune_unused()
        data['AMOEBA_OPBEND_ANGLE_NUM_PARAMS'] = [len(self.out_of_plane_bend_types)]
        data['AMOEBA_OPBEND_ANGLE_FORCE_CONSTANT'] = [at.k for at in self.out_of_plane_bend_types]
        data['AMOEBA_OPBEND_ANGLE_FTAB_DEGREE'] = [self.out_of_plane_bend_types.degree]
        data['AMOEBA_OPBEND_ANGLE_FTAB_COEFFS'] = self.out_of_plane_bend_types.coeffs[:]
        data['AMOEBA_OPBEND_ANGLE_NUM_LIST'] = [len(self.out_of_plane_bends)]
        data['AMOEBA_OPBEND_ANGLE_LIST'] = angle_array = []
        for angle in self.out_of_plane_bends:
            angle_array.extend([angle.atom1.idx+1, angle.atom2.idx+1,
                                angle.atom3.idx+1, angle.atom4.idx+1,
                                angle.type.idx+1])

    #=============================================

    def _xfer_dihedral_info(self):
        """
        Transfers the AMOEBA regular torsion angle information from the topology
        arrays to the raw data arrays
        """
        if len(self.dihedrals) == 0:
            self.delete_flag('AMOEBA_TORSION_NUM_PARAMS')
            self.delete_flag('AMOEBA_TORSION_FORCE_CONSTANT')
            self.delete_flag('AMOEBA_TORSION_PERIODICITY')
            self.delete_flag('AMOEBA_TORSION_PHASE')
            self.delete_flag('AMOEBA_TORSION_NUM_LIST')
            self.delete_flag('AMOEBA_TORSION_LIST')
            return
        data = self.parm_data
        for dihedral_type in self.dihedral_types:
            dihedral_type.used = False
        for dihedral in self.dihedrals:
            dihedral.type.used = True
        self.dihedral_types.prune_unused()
        data['AMOEBA_TORSION_NUM_PARAMS'] = [len(self.dihedral_types)]
        data['AMOEBA_TORSION_FORCE_CONSTANT'] = [dt.phi_k for dt in self.dihedral_types]
        data['AMOEBA_TORSION_PEROIDICITY'] = [dt.per for dt in self.dihedral_types]
        data['AMOEBA_TORSION_PHASE'] = [dt.phase*DEG_TO_RAD for dt in self.dihedral_types]
        data['AMOEBA_TORSION_NUM_LIST'] = [len(self.dihedrals)]
        data['AMOEBA_TORSION_LIST'] = dlist = []
        for dih in self.dihedrals:
            dlist.extend([dih.atom1.idx+1, dih.atom2.idx+1, dih.atom3.idx+1,
                          dih.atom4.idx+1, dih.type.idx+1])

    #=============================================

    def _xfer_pitorsion_info(self):
        """
        Transfers the AMOEBA pi-torsion information from the topology arrays to
        the raw data arrays
        """
        if len(self.pi_torsions) == 0:
            self.delete_flag('AMOEBA_PI_TORSION_NUM_PARAMS')
            self.delete_flag('AMOEBA_PI_TORSION_FORCE_CONSTANT')
            self.delete_flag('AMOEBA_PI_TORSION_PERIODICITY')
            self.delete_flag('AMOEBA_PI_TORSION_PHASE')
            self.delete_flag('AMOEBA_PI_TORSION_NUM_LIST')
            self.delete_flag('AMOEBA_PI_TORSION_LIST')
            return
        data = self.parm_data
        for pitor_type in self.pi_torsion_types:
            pitor_type.used = False
        for pi_torsion in self.pi_torsions:
            pi_torsion.type.used = True
        self.pi_torsion_types.prune_unused()
        data['AMOEBA_PI_TORSION_NUM_PARAMS'] = [len(self.pi_torsion_types)]
        data['AMOEBA_PI_TORSION_FORCE_CONSTANT'] = [dt.phi_k for dt in self.pi_torsion_types]
        data['AMOEBA_PI_TORSION_PEROIDICITY'] = [dt.per for dt in self.pi_torsion_types]
        data['AMOEBA_PI_TORSION_PHASE'] = [dt.phase*DEG_TO_RAD for dt in self.pi_torsion_types]
        data['AMOEBA_PI_TORSION_NUM_LIST'] = [len(self.pi_torsions)]
        data['AMOEBA_PI_TORSION_LIST'] = dlist = []
        for pit in self.pi_torsions:
            dlist.extend([pit.atom1.idx+1, pit.atom2.idx+1, pit.atom3.idx+1,
                          pit.atom4.idx+1, pit.atom5.idx+1, pit.atom6.idx+1,
                          pit.type.idx+1])

    #=============================================

    def _xfer_stretch_bend_info(self):
        """
        Transfers the AMOEBA stretch-bend information from the topology arrays
        to the raw data arrays
        """
        if len(self.stretch_bends) == 0:
            self.delete_flag('AMOEBA_STRETCH_BEND_FORCE_CONSTANT')
            self.delete_flag('AMOEBA_STRETCH_BEND_FORCE_CONSTANT_1')
            self.delete_flag('AMOEBA_STRETCH_BEND_FORCE_CONSTANT_2')
            self.delete_flag('AMOEBA_STRETCH_BEND_BOND1_EQUIL_VALUE')
            self.delete_flag('AMOEBA_STRETCH_BEND_BOND2_EQUIL_VALUE')
            self.delete_flag('AMOEBA_STRETCH_BEND_ANGLE_EQUIL_VALUE')
            self.delete_flag('AMOEBA_STRETCH_BEND_NUM_LIST')
            self.delete_flag('AMOEBA_STRETCH_BEND_LIST')
            return
        # This flag is deprecated... get rid of it and replace it with the 2
        # force constant flags instead
        # TODO: Deprecate the AMOEBA_STRETCH_BEND_FORCE_CONSTANT flag in
        # tinker_to_amber and then do it here as well.
#       self.delete_flag('AMOEBA_STRETCH_BEND_FORCE_CONSTANT')
        data = self.parm_data
        for strbnd_type in self.stretch_bend_types:
            strbnd_type.used = False
        for strbnd in self.stretch_bends:
            strbnd.type.used = True
        self.stretch_bend_types.prune_unused()
        data['AMOEBA_STRETCH_BEND_NUM_PARAMS'] = [len(self.stretch_bend_types)]
#       if not 'AMOEBA_STRETCH_BEND_FORCE_CONSTANT_1' in self.flag_list:
#           self.add_flag('AMOEBA_STRETCH_BEND_FORCE_CONSTANT_1', '5E16.8',
#                   data=[strbnd.k1 for strbnd in self.stretch_bend_types])
#       else:
#           data['AMOEBA_STRETCH_BEND_FORCE_CONSTANT_1'] = \
#                       [strbnd.k1 for strbnd in self.stretch_bend_types]
#       if not 'AMOEBA_STRETCH_BEND_FORCE_CONSTANT_2' in self.flag_list:
#           self.add_flag('AMOEBA_STRETCH_BEND_FORCE_CONSTANT_2', '5E16.8',
#                   data=[strbnd.k2 for strbnd in self.stretch_bend_types])
#       else:
#           data['AMOEBA_STRETCH_BEND_FORCE_CONSTANT_2'] = \
#                       [strbnd.k2 for strbnd in self.stretch_bend_types]
        data['AMOEBA_STRETCH_BEND_FORCE_CONSTANT'] = \
                    [strbnd.k1 for strbnd in self.stretch_bend_types]
        data['AMOEBA_STRETCH_BEND_BOND1_EQUIL_VALUE'] = \
                    [strbnd.req1 for strbnd in self.stretch_bend_types]
        data['AMOEBA_STRETCH_BEND_BOND2_EQUIL_VALUE'] = \
                    [strbnd.req2 for strbnd in self.stretch_bend_types]
        data['AMOEBA_STRETCH_BEND_ANGLE_EQUIL_VALUE'] = \
                    [strbnd.theteq for strbnd in self.stretch_bend_types]
        data['AMOEBA_STRETCH_BEND_NUM_LIST'] = [len(self.stretch_bends)]
        data['AMOEBA_STRETCH_BEND_LIST'] = slist = []
        for strbnd in self.stretch_bends:
            slist.extend([strbnd.atom1.idx+1, strbnd.atom2.idx+1,
                          strbnd.atom3.idx+1, strbnd.type.idx+1])

    #=============================================

    def _xfer_torsion_torsion_info(self):
        """
        Transfers the AMOEBA coupled torsion-torsion information from the
        topology arrays to the raw data arrays
        """
        if len(self.torsion_torsions) == 0:
            delete_flags = set(
                flag for flag in self.flag_list if flag.startswith('AMOEBA_TORSION_TORSION')
            )
            for flag in delete_flags:
                self.delete_flag(flag)
            return
        data = self.parm_data
        for tortor_type in self.torsion_torsion_types:
            tortor_type.used = False
        for tortor in self.torsion_torsions:
            tortor.type.used = True
        self.torsion_torsion_types.prune_unused()
        after = 'AMOEBA_TORSION_TORSION_NUM_PARAMS'
        data[after] = [len(self.torsion_torsion_types)]
        delete_flags = set(flag for flag in self.flag_list
                    if flag.startswith('AMOEBA_TORSION_TORSION_TORTOR_TABLE'))
        for flag in delete_flags:
            self.delete_flag(flag)
        for i, tt in enumerate(self.torsion_torsion_types):
            tblsize = ['dimension = (%d,%d)' % tt.dims]
            prefix = 'AMOEBA_TORSION_TORSION_TORTOR_TABLE_%02d_' % (i+1)
            self.add_flag(prefix+'DIMS', '2I8', data=list(tt.dims),
                          comments=['dimension = (2)'], after=after)
            self.add_flag(prefix+'ANGLE1', '5E16.8', data=tt.ang1[:],
                          comments=['dimension = (%d)' % tt.dims[0]], after=prefix+'DIMS')
            self.add_flag(prefix+'ANGLE2', '5E16.8', data=tt.ang2[:],
                          comments=['dimension = (%d)' % tt.dims[1]], after=prefix+'ANGLE1')
            self.add_flag(prefix+'FUNC', '5E16.8', data=tt.f.data[:],
                          comments=tblsize[:], after=prefix+'ANGLE2')
            self.add_flag(prefix+'DFUNC_DANGLE1', '5E16.8', data=tt.dfda1.data[:],
                          comments=tblsize[:], after=prefix+'FUNC')
            self.add_flag(prefix+'DFUNC_DANGLE2', '5E16.8', data=tt.dfda2.data[:],
                          comments=tblsize[:], after=prefix+'DFUNC_DANGLE1')
            self.add_flag(prefix+'D2FUNC_DANGLE1_DANGLE2', '5E16.8', data=tt.d2fda1da2.data[:],
                          comments=tblsize[:], after=prefix+'DFUNC_DANGLE2')
            after = prefix + 'D2FUNC_DANGLE1_DANGLE2'
        data['AMOEBA_TORSION_TORSION_NUM_LIST'] = [len(self.torsion_torsions)]
        data['AMOEBA_TORSION_TORSION_LIST'] = tlist = []
        for tortor in self.torsion_torsions:
            tlist.extend([tortor.atom1.idx+1, tortor.atom2.idx+1,
                          tortor.atom3.idx+1, tortor.atom4.idx+1,
                          tortor.atom5.idx+1, tortor.type.idx+1])

    #=============================================

    def _xfer_frame_info(self):
        """
        Transfers the chiral and multipole frame data from the topology arrays
        to the raw data arrays
        """
        data = self.parm_data
        if len(self.chiral_frames) > 0:
            data['AMOEBA_CHIRAL_FRAME_NUM_LIST'] = [len(self.chiral_frames)]
            data['AMOEBA_CHIRAL_FRAME_LIST'] = clist = []
            for cf in self.chiral_frames:
                clist.extend([cf.atom1.idx+1, cf.atom2.idx+1, cf.chirality])
        else:
            self.delete_flag('AMOEBA_CHIRAL_FRAME_NUM_LIST')
            self.delete_flag('AMOEBA_CHIRAL_FRAME_LIST')

        if len(self.multipole_frames) > 0:
            data['AMOEBA_FRAME_DEF_NUM_LIST'] = [len(self.multipole_frames)]
            data['AMOEBA_FRAME_DEF_LIST'] = flist = []
            for mf in self.multipole_frames:
                flist.extend([mf.atom.idx+1, mf.frame_pt_num, mf.vectail, mf.vechead, mf.nvec])
        else:
            self.delete_flag('AMOEBA_FRAME_DEF_NUM_LIST')
            self.delete_flag('AMOEBA_FRAME_DEF_LIST')

    #=============================================

    def _xfer_exception_info(self):
        """
        Transfers the nonboned exception (adjust) info from the topology arrays
        to the raw data arrays
        """
        # adjust type arrays hard-coded in length... do not purge unused.
        data = self.parm_data
        data['AMOEBA_ADJUST_VDW_WEIGHTS_LIST'] = [at.vdw_weight for at in self.adjust_types]
        data['AMOEBA_ADJUST_MPOLE_WEIGHTS_LIST'] = [at.multipole_weight for at in self.adjust_types]
        data['AMOEBA_ADJUST_DIRECT_WEIGHTS_LIST'] = [at.direct_weight for at in self.adjust_types]
        data['AMOEBA_ADJUST_POLAR_WEIGHTS_LIST'] = [at.polar_weight for at in self.adjust_types]
        data['AMOEBA_ADJUST_MUTUAL_WEIGHTS_LIST'] = [at.mutual_weight for at in self.adjust_types]
        data['AMOEBA_ADJUST_NUM_LIST'] = [len(self.adjusts)]
        data['AMOEBA_ADJUST_LIST'] = alist = []
        for adj in self.adjusts:
            alist.extend([adj.atom1.idx + 1, adj.atom2.idx + 1, adj.type.idx + 1])

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class BeemanRestart(AmberFormat):
    """
    The restart files written for/by the Beeman integrator has the same type of
    format as the topology file
    """

    @classmethod
    def from_rawdata(cls, rawdata):
        """
        Take the raw data from a AmberFormat object and initialize a
        BeemanRestart from that data.

        Parameters
        ----------
        rawdata : :class:`AmberFormat`
            An AmberFormat instance that has already been instantiated

        Returns
        -------
        inst : :class:`BeemanRestart`
            An instance of this type from the data in rawdata
        """
        inst = cls()
        inst.name = rawdata.name
        inst.version = rawdata.version
        inst.formats = rawdata.formats
        inst.parm_data = rawdata.parm_data
        inst.parm_comments = rawdata.parm_comments
        inst.flag_list = rawdata.flag_list
        return inst

    @property
    def natom(self):
        return self.parm_data['ATOMIC_COORDS_NUM_LIST'][0]

    @natom.setter
    def natom(self, value):
        value = int(value) * 3
        old_natom3 = self.natom * 3
        if value > old_natom3:
            # Getting bigger
            zeros = [0 for i in range(value)]
            new_data = zeros[:]
            new_data[:old_natom3] = self.parm_data['ATOMIC_COORDS_LIST']
            self.parm_data['ATOMIC_COORDS_LIST'] = new_data
            if 'ATOMIC_VELOCITIES_LIST' in self.parm_data:
                new_data = zeros[:]
                new_data[:old_natom3] = self.parm_data['ATOMIC_VELOCITIES_LIST']
                self.parm_data['ATOMIC_VELOCITIES_LIST'] = new_data
            if 'ATOMIC_ACCELERATIONS_LIST' in self.parm_data:
                new_data = zeros[:]
                new_data[:old_natom3] = self.parm_data['ATOMIC_ACCELERATIONS_LIST']
                self.parm_data['ATOMIC_ACCELERATIONS_LIST'] = new_data
            if 'OLD_ATOMIC_ACCELERATIONS_LIST' in self.parm_data:
                new_data = zeros[:]
                new_data[:old_natom3] = self.parm_data['OLD_ATOMIC_ACCELERATIONS_LIST']
                self.parm_data['OLD_ATOMIC_ACCELERATIONS_LIST'] = new_data
        else:
            # Getting smaller
            del self.parm_data['ATOMIC_COORDS_LIST'][value:]
            if 'ATOMIC_VELOCITIES_LIST' in self.parm_data:
                del self.parm_data['ATOMIC_VELOCITIES_LIST'][value:]
            if 'ATOMIC_ACCELERATIONS_LIST' in self.parm_data:
                del self.parm_data['ATOMIC_ACCELERATIONS_LIST'][value:]
            if 'OLD_ATOMIC_ACCELERATIONS_LIST' in self.parm_data:
                del self.parm_data['OLD_ATOMIC_ACCELERATIONS_LIST'][value:]
        value //= 3
        self.parm_data['ATOMIC_COORDS_NUM_LIST'][0] = value
        if 'ATOMIC_VELOCITIES_LIST' in self.parm_data:
            self.parm_data['ATOMIC_VELOCITIES_NUM_LIST'][0] = value
        if 'ATOMIC_ACCELERATIONS_LIST' in self.parm_data:
            self.parm_data['ATOMIC_ACCELERATIONS_NUM_LIST'][0] = value
        if 'OLD_ATOMIC_ACCELERATIONS_NUM_LIST' in self.parm_data:
            self.parm_data['OLD_ATOMIC_ACCELERATIONS_NUM_LIST'][0] = value

    @property
    def coordinates(self):
        return np.array(self.parm_data['ATOMIC_COORDS_LIST']).reshape((1, self.natom, 3))

    @coordinates.setter
    def coordinates(self, value):
        value = np.asarray(value).flatten()
        if value.shape != (3*self.natom,):
            raise ValueError(f'Require {3 * self.natom}-length sequence for coordinates')
        self.parm_data['ATOMIC_COORDS_LIST'] = value.tolist()

    @property
    def velocities(self):
        try:
            return np.array(self.parm_data['ATOMIC_VELOCITIES_LIST']).reshape((1, self.natom, 3))
        except KeyError:
            raise AttributeError('Beeman restart does not have velocities')

    @velocities.setter
    def velocities(self, value):
        value = np.asarray(value).flatten()
        if value.shape != (3*self.natom,):
            raise ValueError(f'Require {3 * self.natom}-length sequence for velocities')
        if not 'ATOMIC_VELOCITIES_LIST' in self.flag_list:
            self.add_flag('ATOMIC_VELOCITIES_NUM_LIST', 'i8', data=[self.natom])
            self.add_flag('ATOMIC_VELOCITIES_LIST', '3e20.12', data=value.tolist())
        else:
            self.parm_data['ATOMIC_VELOCITIES_LIST'] = value.tolist()

    @property
    def accelerations(self):
        try:
            return np.array(
                self.parm_data['ATOMIC_ACCELERATIONS_LIST']).reshape((1, self.natom, 3)
            )
        except KeyError:
            raise AttributeError('Accelerations not present in Beeman restart')

    @accelerations.setter
    def accelerations(self, value):
        value = np.asarray(value).flatten()
        if value.shape != (3*self.natom,):
            raise ValueError(f'Require {3 * self.natom}-length sequence for accelerations')
        if not 'ATOMIC_ACCELERATIONS_LIST' in self.flag_list:
            self.add_flag('ATOMIC_ACCELERATIONS_NUM_LIST', 'i8', data=[self.natom])
            self.add_flag('ATOMIC_ACCELERATIONS_LIST', '3e20.12', data=value.tolist())
        else:
            self.parm_data['ATOMIC_ACCELERATIONS_LIST'] = value.tolist()

    @property
    def old_accelerations(self):
        try:
            return np.array(
                self.parm_data['OLD_ATOMIC_ACCELERATIONS_LIST']).reshape((1, self.natom, 3)
            )
        except KeyError:
            raise AttributeError('Old accelerations not present in Beeman restart')

    @old_accelerations.setter
    def old_accelerations(self, value):
        value = np.asarray(value).flatten()
        if value.shape != (3*self.natom,):
            raise ValueError(f'Require {3 * self.natom}-length sequence for accelerations')
        if not 'OLD_ATOMIC_ACCELERATIONS_LIST' in self.flag_list:
            self.add_flag('OLD_ATOMIC_ACCELERATIONS_NUM_LIST', 'i8', data=[self.natom])
            self.add_flag('OLD_ATOMIC_ACCELERATIONS_LIST', '3e20.12', data=value.tolist())
        else:
            self.parm_data['OLD_ATOMIC_ACCELERATIONS_LIST'] = value.tolist()

    @property
    def box(self):
        return np.array(self.parm_data['UNIT_CELL_PARAMETERS'])

    @box.setter
    def box(self, stuff):
        if len(stuff) != 6:
            raise ValueError('Expected 3 box lengths and 3 box angles')
        self.parm_data['UNIT_CELL_PARAMETERS'] = list(stuff)
