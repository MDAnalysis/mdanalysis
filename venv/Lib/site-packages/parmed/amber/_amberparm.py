"""
This module contains an amber prmtop class that will read in all
parameters and allow users to manipulate that data and write a new
prmtop object.
"""
import copy as _copy
from collections import defaultdict
from math import sqrt
from warnings import warn

import numpy as np

from .. import unit as u
from ..constants import PrmtopPointers, TRUNCATED_OCTAHEDRON_ANGLE, RAD_TO_DEG, SMALL, DEG_TO_RAD
from ..exceptions import AmberError, AmberWarning, MoleculeError
from ..geometry import box_lengths_and_angles_to_vectors
from ..periodic_table import AtomicNum, element_by_mass
from ..residue import ALLION_NAMES, SOLVENT_NAMES
from ..structure import Structure, needs_openmm
from ..topologyobjects import (Angle, AngleType, Atom, AtomList, AtomType,
                               Bond, BondType, Dihedral, DihedralType,
                               DihedralTypeList, ExtraPoint, Cmap, CmapType)
from ..vec3 import Vec3
from .amberformat import AmberFormat
from .asciicrd import AmberAsciiRestart
from .netcdffiles import NetCDFRestart

try:
    import openmm as mm
    from openmm import app
except ImportError:
    mm = app = None

class AmberParm(AmberFormat, Structure):
    """
    Amber Topology (parm7 format) class. Gives low, and some high, level access
    to topology data. You can interact with the raw data in the topology file
    directly or interact with some of the high-level classes comprising the
    system topology and parameters.

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
    charge_flag : str='CHARGE'
        The name of the name of the FLAG that describes partial atomic charge
        data. If this flag is found, then its data are multiplied by the
        ELECTROSTATIC_CONSTANT to convert back to fractions of electrons
    version : str
        The VERSION string from the Amber file
    name : str
        The file name of the originally parsed file (set to the fname parameter)
    LJ_types : dict {str : int}
        A mapping whose keys are atom types paired with the nonbonded index of
        that type
    LJ_radius : list(float)
        A list of floats corresponding to the Rmin/2 parameter for every
        Lennard-Jones type. The indexes are the nonbonded index (`nb_idx`
        attribute of the `Atom` class) minus 1 to account for indexing from 0 in
        Python. The values are in Angstroms. To get the radius for a particular
        atom type, you can use `LJ_radius[LJ_types["type"]-1]`
    LJ_depth : list(float)
        A list of Lennard-Jones epsilon parameters laid out the same way as
        LJ_radius, described above.
    atoms : AtomList(Atom)
        List of all atoms in the system
    residues : ResidueList(Residue)
        List of all residues in the system
    bonds : TrackedList(Bond)
        List of bonds between two atoms in the system
    angles : TrackedList(Angle)
        List of angles between three atoms in the system
    dihedrals : TrackedList(Angle)
        List of all proper and improper torsions between 4 atoms in the system
    box : list of 6 floats
        Periodic boundary unit cell dimensions and angles
    bond_types : TrackedList(BondType)
        The bond types containing the parameters for each bond stretching term
    angle_types : TrackedList(AngleType)
        The angle types containing the parameters for each angle bending term
    dihedral_types : TrackedList(DihedralType)
        The dihedral types containing the parameters for each torsional term
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

    _cmap_prefix = ""

    #===================================================

    def __init__(self, prm_name=None, xyz=None, box=None):
        """
        Instantiates an AmberParm object from data in prm_name and establishes
        validity based on presence of POINTERS and CHARGE sections. In general,
        you should use LoadParm from the readparm module instead. LoadParm will
        correctly dispatch the object to the 'correct' flavor of AmberParm
        """
        AmberFormat.__init__(self, prm_name)
        Structure.__init__(self)
        self.hasvels = False
        self.hasbox = False
        self._box = None
        self.crdname = None
        if isinstance(xyz, str):
            self.crdname = xyz
        if prm_name is not None:
            self.initialize_topology(xyz, box)

    #===================================================

    def initialize_topology(self, xyz=None, box=None):
        """
        Initializes topology data structures, like the list of atoms, bonds,
        etc., after the topology file has been read.
        """
        from .. import load_file
        # We need to handle RESIDUE_ICODE properly since it may have picked up
        # some extra values
        if 'RESIDUE_ICODE' in self.flag_list:
            self._truncate_array('RESIDUE_ICODE', self.parm_data['POINTERS'][PrmtopPointers.NRES])
        # Set up some of the attributes provided
        self.pointers = {}
        self.LJ_types = {}
        self.LJ_radius = []
        self.LJ_depth = []
        self.load_pointers()
        self.fill_LJ()
        # The way we check if the combining rules are geometric vs.
        # Lorentz-Berthelot is to try and change the combining rules to
        # geometric *if we detect NBFIXes* and see if we still have NBFIXes. If
        # not, then our combining rules are clearly geometric.
        if not self.has_1012() and self.has_NBFIX():
            self.combining_rule = 'geometric'
            if self.has_NBFIX():
                self.combining_rule = 'lorentz'
        # Instantiate the Structure data structures
        self.load_structure()

        if isinstance(xyz, str):
            f = load_file(xyz, skip_bonds=True)
            if not hasattr(f, 'coordinates') or f.coordinates is None:
                raise TypeError('%s does not have coordinates' % xyz)
            self.coordinates = f.coordinates
            if hasattr(f, 'velocities') and f.velocities is not None:
                self.velocities = f.velocities
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

    #===================================================

    @classmethod
    def from_rawdata(cls, rawdata):
        """
        Take the raw data from a AmberFormat object and initialize an AmberParm
        from that data.

        Parameters
        ----------
        rawdata : :class:`AmberFormat`
            An AmberFormat instance that has already been instantiated

        Returns
        -------
        parm : :class:`AmberParm`
            An instance of this type from the data in rawdata
        """
        inst = cls()
        inst.name = rawdata.name
        inst.version = rawdata.version
        inst.formats = rawdata.formats
        inst.parm_data = rawdata.parm_data
        inst.parm_comments = rawdata.parm_comments
        inst.flag_list = rawdata.flag_list
        inst.initialize_topology()
        # See if the rawdata has any kind of structural attributes, like
        # coordinates and an atom list with positions and/or velocities
        if hasattr(rawdata, 'coordinates'):
            inst.coordinates = _copy.copy(rawdata.coordinates)
        if hasattr(rawdata, 'velocities'):
            inst.velocities = _copy.copy(rawdata.velocities)
        if hasattr(rawdata, 'box'):
            inst.box = _copy.copy(rawdata.box)
        inst.hasbox = inst.box is not None
        inst.hasvels = inst.velocities is not None
        n_copy = inst.pointers.get('NCOPY', 1)
        if n_copy >= 2:
            inst._label_alternates()
        return inst

    #===================================================

    @classmethod
    def from_structure(cls, struct, copy=False):
        """
        Take a Structure instance and initialize an AmberParm instance from that
        data.

        Parameters
        ----------
        struct : :class:`Structure`
            The input structure from which to construct an AmberParm instance
        copy : bool
            If True, the input struct is deep-copied to make sure it does not
            share any objects with the original ``struct``. Default is False

        Returns
        -------
        inst : :class:`AmberParm`
            The AmberParm instance derived from the input structure

        Raises
        ------
        TypeError
            If the structure has parameters not supported by the standard Amber
            force field (i.e., standard bond, angle, and dihedral types)

        Notes
        -----
        Due to the nature of the prmtop file, struct almost *always* returns a
        deep copy. The one exception is when struct is already of type
        :class:`AmberParm`, in which case the original object is returned unless
        ``copy`` is ``True``.
        """
        if type(struct) is cls:
            if copy:
                return _copy.copy(struct)
            return struct
        if struct.unknown_functional:
            raise TypeError('Cannot instantiate an AmberParm from unknown functional')
        if (struct.urey_bradleys or struct.impropers or
                struct.trigonal_angles or struct.pi_torsions or
                struct.out_of_plane_bends or struct.stretch_bends or
                struct.torsion_torsions or struct.multipole_frames):
            if (struct.trigonal_angles or struct.pi_torsions or
                    struct.out_of_plane_bends or struct.torsion_torsions or
                    struct.multipole_frames or struct.stretch_bends):
                raise TypeError('AmberParm does not support all of the '
                                'parameters defined in the input Structure')
            # Maybe it just has CHARMM parameters?
            raise TypeError('AmberParm does not support all of the parameters '
                            'defined in the input Structure. Try ChamberParm')

        # Convert all RB torsions to propers and delete the originals.
        for dihedral in struct.rb_torsions:
            proper_types = DihedralTypeList.from_rbtorsion(dihedral.type)
            for proper_type in proper_types:
                proper = _copy.copy(dihedral)
                proper.type = proper_type
                struct.dihedrals.append(proper)
                struct.dihedral_types.append(proper.type)

        del struct.rb_torsions[:]
        del struct.rb_torsion_types[:]
        struct.dihedrals.claim()
        struct.dihedral_types.claim()

        inst = struct.copy(cls, split_dihedrals=True)
        inst.update_dihedral_exclusions()
        inst._add_missing_13_14()
        del inst.adjusts[:]
        inst.pointers = {}
        inst.LJ_types = {}
        nbfixes = inst._set_nbidx_lj_params()
        inst._add_standard_flags()
        inst.pointers['NATOM'] = len(inst.atoms)
        inst.parm_data['POINTERS'][PrmtopPointers.NATOM] = len(inst.atoms)
        # pmemd likes to skip torsions with periodicities of 0, which may be
        # present as a way to hack entries into the 1-4 pairlist. See
        # https://github.com/ParmEd/ParmEd/pull/145 for discussion. The solution
        # here is to simply set that periodicity to 1.
        inst._cleanup_dihedrals_with_periodicity_zero()
        inst.remake_parm()
        inst._set_nonbonded_tables(nbfixes)
        n_copy = inst.pointers.get('NCOPY', 1)
        if n_copy >= 2:
            inst._label_alternates()

        # Setting box sets some of the flag data appropriately
        inst.box = inst.get_box()
        return inst

    #===================================================

    def _set_nbidx_lj_params(self):
        nbfixes = self.atoms.assign_nbidx_from_types()
        # Give virtual sites a name that Amber understands
        for atom in self.atoms:
            atom.type = atom.type if not isinstance(atom, ExtraPoint) else 'EP'
        # Fill the Lennard-Jones arrays/dicts
        ntyp = 0
        for atom in self.atoms:
            self.LJ_types[atom.type] = atom.nb_idx
            ntyp = max(ntyp, atom.nb_idx)
        self.LJ_radius = [0 for i in range(ntyp)]
        self.LJ_depth = [0 for i in range(ntyp)]
        for atom in self.atoms:
            self.LJ_radius[atom.nb_idx-1] = atom.atom_type.rmin
            self.LJ_depth[atom.nb_idx-1] = atom.atom_type.epsilon
        return nbfixes

    #===================================================

    def __copy__(self):
        """ Needs to copy a few additional data structures """
        other = super(AmberParm, self).__copy__()
        other.initialize_topology()
        # Copy coordinates, box, etc. if applicable
        other.coordinates = _copy.copy(self.get_coordinates('all'))
        other._box = _copy.copy(self._box)
        # Now we should have a full copy
        return other

    #===================================================

    def __getitem__(self, selection):
        other = super(AmberParm, self).__getitem__(selection)
        if isinstance(other, Atom):
            return other
        other.pointers = {}
        self._copy_lj_data(other)
        other._add_standard_flags()
        other.remake_parm()
        other._set_nonbonded_tables()
        # Add back the original L-J tables to restore any NBFIXes that may have
        # been there
        other.parm_data['LENNARD_JONES_ACOEF'] = self.parm_data['LENNARD_JONES_ACOEF'][:]
        other.parm_data['LENNARD_JONES_BCOEF'] = self.parm_data['LENNARD_JONES_BCOEF'][:]
        if 'LENNARD_JONES_CCOEF' in self.parm_data:
            other.parm_data['LENNARD_JONES_CCOEF'] = self.parm_data['LENNARD_JONES_CCOEF'][:]
        return other

    def _copy_lj_data(self, other):
        """ Copies Lennard-Jones lists and dicts from myself to a copy """
        other.LJ_types = self.LJ_types.copy()
        other.LJ_radius = _copy.copy(self.LJ_radius)
        other.LJ_depth = _copy.copy(self.LJ_depth)
        for atom in other.atoms:
            other.LJ_radius[atom.nb_idx-1] = atom.atom_type.rmin
            other.LJ_depth[atom.nb_idx-1] = atom.atom_type.epsilon

    #===================================================

    def __imul__(self, ncopies, other=None):
        super(AmberParm, self).__imul__(ncopies, other)
        self.remake_parm()
        return self

    def __iadd__(self, other):
        if self.has_NBFIX() or self.has_1012():
            raise ValueError('Cannot combine Amber systems with NBFIX or 10-12 parameters')
        super(AmberParm, self).__iadd__(other)
        # Make sure we properly recompute the LJ tables
        nbfixes = self._set_nbidx_lj_params()
        self.remake_parm()
        self._set_nonbonded_tables(nbfixes)
        return self

    #===================================================

    def load_pointers(self):
        """
        Loads the data in POINTERS section into a pointers dictionary with each
        key being the pointer name according to http://ambermd.org/formats.html
        """
        self.pointers["NATOM"] = self.parm_data["POINTERS"][PrmtopPointers.NATOM]
        self.pointers["NTYPES"] = self.parm_data["POINTERS"][PrmtopPointers.NTYPES]
        self.pointers["NBONH"] = self.parm_data["POINTERS"][PrmtopPointers.NBONH]
        self.pointers["MBONA"] = self.parm_data["POINTERS"][PrmtopPointers.MBONA]
        self.pointers["NTHETH"] = self.parm_data["POINTERS"][PrmtopPointers.NTHETH]
        self.pointers["MTHETA"] = self.parm_data["POINTERS"][PrmtopPointers.MTHETA]
        self.pointers["NPHIH"] = self.parm_data["POINTERS"][PrmtopPointers.NPHIH]
        self.pointers["MPHIA"] = self.parm_data["POINTERS"][PrmtopPointers.MPHIA]
        self.pointers["NHPARM"] = self.parm_data["POINTERS"][PrmtopPointers.NHPARM]
        self.pointers["NPARM"] = self.parm_data["POINTERS"][PrmtopPointers.NPARM]
        self.pointers["NEXT"] = self.parm_data["POINTERS"][PrmtopPointers.NEXT]
        self.pointers["NNB"] = self.parm_data["POINTERS"][PrmtopPointers.NNB] # alias for above
        self.pointers["NRES"] = self.parm_data["POINTERS"][PrmtopPointers.NRES]
        self.pointers["NBONA"] = self.parm_data["POINTERS"][PrmtopPointers.NBONA]
        self.pointers["NTHETA"] = self.parm_data["POINTERS"][PrmtopPointers.NTHETA]
        self.pointers["NPHIA"] = self.parm_data["POINTERS"][PrmtopPointers.NPHIA]
        self.pointers["NUMBND"] = self.parm_data["POINTERS"][PrmtopPointers.NUMBND]
        self.pointers["NUMANG"] = self.parm_data["POINTERS"][PrmtopPointers.NUMANG]
        self.pointers["NPTRA"] = self.parm_data["POINTERS"][PrmtopPointers.NPTRA]
        self.pointers["NATYP"] = self.parm_data["POINTERS"][PrmtopPointers.NATYP]
        self.pointers["NPHB"] = self.parm_data["POINTERS"][PrmtopPointers.NPHB]
        self.pointers["IFPERT"] = self.parm_data["POINTERS"][PrmtopPointers.IFPERT]
        self.pointers["NBPER"] = self.parm_data["POINTERS"][PrmtopPointers.NBPER]
        self.pointers["NGPER"] = self.parm_data["POINTERS"][PrmtopPointers.NGPER]
        self.pointers["NDPER"] = self.parm_data["POINTERS"][PrmtopPointers.NDPER]
        self.pointers["MBPER"] = self.parm_data["POINTERS"][PrmtopPointers.MBPER]
        self.pointers["MGPER"] = self.parm_data["POINTERS"][PrmtopPointers.MGPER]
        self.pointers["MDPER"] = self.parm_data["POINTERS"][PrmtopPointers.MDPER]
        self.pointers["IFBOX"] = self.parm_data["POINTERS"][PrmtopPointers.IFBOX]
        self.pointers["NMXRS"] = self.parm_data["POINTERS"][PrmtopPointers.NMXRS]
        self.pointers["IFCAP"] = self.parm_data["POINTERS"][PrmtopPointers.IFCAP]
        self.pointers["NUMEXTRA"] = self.parm_data["POINTERS"][PrmtopPointers.NUMEXTRA]
        if self.parm_data['POINTERS'][PrmtopPointers.IFBOX] > 0:
            self.pointers['IPTRES'] = self.parm_data['SOLVENT_POINTERS'][0]
            self.pointers['NSPM'] = self.parm_data['SOLVENT_POINTERS'][1]
            self.pointers['NSPSOL'] = self.parm_data['SOLVENT_POINTERS'][2]
        # The next is probably only there for LES-prmtops
        try:
            self.pointers["NCOPY"] = self.parm_data["POINTERS"][PrmtopPointers.NCOPY]
        except (KeyError, IndexError):
            pass
        # If CMAP is not present, don't load the pointers
        if self.has_cmap:
            self.pointers['CMAP'] = self.parm_data[f"{self._cmap_prefix}CMAP_COUNT"][0]
            self.pointers['CMAP_TYPES'] = self.parm_data[f"{self._cmap_prefix}CMAP_COUNT"][1]

    #===================================================

    def load_structure(self):
        """
        Loads all of the topology instance variables. This is necessary if we
        actually want to modify the topological layout of our system
        (like deleting atoms)
        """
        self._check_section_lengths()
        self._load_atoms_and_residues()
        self.load_atom_info()
        self._load_bond_info()
        self._load_angle_info()
        self._load_dihedral_info()
        self._load_cmap_info()
        self._load_extra_exclusions()
        super().unchange()

    #===================================================

    def load_atom_info(self):
        """
        Loads atom properties into the atoms that have been loaded. If any
        arrays are too short or too long, an IndexError will be raised
        """
        # Collect all of the atom properties present in our topology file
        zeros = _zeros(len(self.atoms))
        anam = self.parm_data['ATOM_NAME']
        chg = self.parm_data['CHARGE']
        mass = self.parm_data['MASS']
        nbtyp = self.parm_data['ATOM_TYPE_INDEX']
        atyp = self.parm_data['AMBER_ATOM_TYPE']
        join = self.parm_data['JOIN_ARRAY']
        irot = self.parm_data['IROTAT']
        tree = self.parm_data['TREE_CHAIN_CLASSIFICATION']
        try:
            radii = self.parm_data['RADII']
        except KeyError:
            radii = zeros
        try:
            screen = self.parm_data['SCREEN']
        except KeyError:
            screen = zeros
        try:
            atnum = self.parm_data['ATOMIC_NUMBER']
            replace_atnum = True
        except KeyError:
            atnum = [AtomicNum[element_by_mass(m)] for m in mass]
            replace_atnum = False
        try:
            occu = self.parm_data['ATOM_OCCUPANCY']
        except KeyError:
            occu = zeros
        try:
            bfac = self.parm_data['ATOM_BFACTOR']
        except KeyError:
            bfac = zeros
        try:
            anum = self.parm_data['ATOM_NUMBER']
        except KeyError:
            anum = [-1 for atom in self.atoms]
        for i, atom in enumerate(self.atoms):
            atom.name = anam[i]
            atom.charge = chg[i]
            atom.mass = mass[i]
            atom.nb_idx = nbtyp[i]
            atom.type = atyp[i]
            atom.join = join[i]
            atom.irotat = irot[i]
            atom.tree = tree[i]
            atom.solvent_radius = radii[i]
            atom.screen = screen[i]
            if replace_atnum or atom.atomic_number == 0:
                if atnum[i] == -1:
                    atom.atomic_number = AtomicNum[element_by_mass(mass[i])]
                else:
                    atom.atomic_number = atnum[i]
            atom.atom_type = AtomType(atyp[i], None, mass[i], atnum[i])
            atom.occupancy = occu[i]
            atom.bfactor = bfac[i]
            atom.number = anum[i]
            depth = self.LJ_depth[atom.nb_idx-1]
            radius = self.LJ_radius[atom.nb_idx-1]
            try:
                depth14 = self.LJ_14_depth[atom.nb_idx-1]
                radius14 = self.LJ_14_radius[atom.nb_idx-1]
            except AttributeError:
                depth14 = radius14 = None
            atom.atom_type.set_lj_params(depth, radius, depth14, radius14)

    #===================================================

    def ptr(self, pointer: str) -> int:
        """
        Returns the value of the given pointer, and converts to upper-case so
        it's case-insensitive. A non-existent pointer meets with a KeyError

        Parameters
        ----------
        pointer : str
            The AMBER pointer for which to extract the value

        Returns
        -------
        int
            The returned integer is the value of that pointer
        """
        return self.pointers[pointer.upper()]

    #===================================================

    def write_rst7(self, name, netcdf=None):
        """
        Writes a restart file with the current coordinates and velocities and
        box info if it's present

        Parameters
        ----------
        name : str
            Name of the file to write the restart file to
        netcdf : bool=False
            If True, write a NetCDF-format restart file (requires a NetCDF
            backend; scipy, netCDF4, or ScientificPython; to be installed)

        Notes
        -----
        If `netcdf` is not specified and the filename extension given by `name`
        is `.ncrst`, the a NetCDF restart file will be written. However, an
        explicit value for `netcdf` will override any filename extensions.
        """
        # By default, determine file type by extension (.ncrst is NetCDF)
        netcdf = netcdf or (netcdf is None and name.endswith('.ncrst'))

        # Check that we have a rst7 loaded, then overwrite it with a new one if
        # necessary
        rst7 = Rst7(natom=len(self.atoms))

        # Now fill in the rst7 coordinates
        rst7.coordinates = [0.0 for i in range(len(self.atoms)*3)]
        if self.velocities is not None:
            rst7.vels = [0.0 for i in range(len(self.atoms)*3)]

        for i, at in enumerate(self.atoms):
            i3 = i * 3
            rst7.coordinates[i3  ] = at.xx
            rst7.coordinates[i3+1] = at.xy
            rst7.coordinates[i3+2] = at.xz
            if rst7.hasvels:
                rst7.vels[i3  ] = at.vx
                rst7.vels[i3+1] = at.vy
                rst7.vels[i3+2] = at.vz

        rst7.box = _copy.copy(self.box)
        # Now write the restart file
        rst7.write(name, netcdf)

    #===================================================

    def write_parm(self, name):
        """
        Writes the current data in parm_data into a new topology file with a given name.

        Parameters
        ----------
        name : str or file-like
            The name of the file to write the prmtop to or the file object to write to
        """
        self.remake_parm()
        AmberFormat.write_parm(self, name)

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
        self._xfer_cmap_properties()
        self._set_ifbox()
        # Load the pointers dict
        self.load_pointers()
        # Mark atom list as unchanged
        super().unchange()

    #===================================================

    def is_changed(self):
        """
        Determines if any of the topological arrays have changed since the
        last upload
        """
        is_changed = super(AmberParm, self).is_changed()
        if is_changed and hasattr(self, '_topology'):
            del self._topology
        return is_changed

    #===================================================

    def strip(self, selection):
        """
        Deletes a subset of the atoms corresponding to an atom-based selection.

        Parameters
        ----------
        selection : AmberMask, str, or iterable of bool
            This is the selection of atoms that will be deleted from this
            structure. If it is a string, it will be interpreted as an
            AmberMask. If it is an AmberMask, it will be converted to a
            selection of atoms. If it is an iterable, it must be the same length
            as the `atoms` list.
        """
        super(AmberParm, self).strip(selection)
        self.remake_parm()

    #===================================================

    def rediscover_molecules(self, solute_ions=True, fix_broken=True):
        """
        This determines the molecularity and sets the ATOMS_PER_MOLECULE and
        SOLVENT_POINTERS sections of the prmtops. Returns the new atom sequence
        in terms of the 'old' atom indexes if re-ordering was necessary to fix
        the tleap bug. Returns None otherwise.
        """
        from ..utils import tag_molecules
        # Bail out of we are not doing a solvated prmtop
        if self.parm_data['POINTERS'][PrmtopPointers.IFBOX] == 0 or self.box is None:
            return None

        owner = tag_molecules(self)
        all_solvent = SOLVENT_NAMES
        if not solute_ions:
            all_solvent = all_solvent | ALLION_NAMES
        first_solvent = None
        for i, res in enumerate(self.residues):
            if res.name in all_solvent:
                first_solvent = i
                break
        # Now remake our SOLVENT_POINTERS and ATOMS_PER_MOLECULE section
        self.parm_data['SOLVENT_POINTERS'] = [first_solvent, len(owner), 0]
        # Find the first solvent molecule if there are any
        if first_solvent is None:
            self.parm_data['SOLVENT_POINTERS'][0] = len(self.residues)
            self.parm_data['SOLVENT_POINTERS'][2] = len(owner) + 1
        else:
            self.parm_data['SOLVENT_POINTERS'][2] = self.residues[first_solvent].atoms[0].marked

        # Now set up ATOMS_PER_MOLECULE and catch any errors
        self.parm_data['ATOMS_PER_MOLECULE'] = [len(mol) for mol in owner]

        # Check that all of our molecules are contiguous, because we have to
        # re-order atoms if they're not
        try:
            for mol in owner:
                sortedmol = sorted(list(mol))
                for i in range(1, len(sortedmol)):
                    if sortedmol[i] != sortedmol[i-1] + 1:
                        raise StopIteration()
        except StopIteration:
            if not fix_broken:
                raise MoleculeError('Molecule atoms are not contiguous!')
            # Non-contiguous molecules detected... time to fix (ugh!)
            warn('Molecule atoms are not contiguous! Attempting to reorder to fix.', AmberWarning)
            # Make sure that no residues are split up by this
            for res in self.residues:
                molid = res.atoms[0].marked
                for atom in res:
                    if molid != atom.marked:
                        warn('Residues cannot be part of 2 molecules! Molecule '
                             'section will not be correctly set. [Offending '
                             f'residue is {repr(res)}', AmberWarning)
                        return None
            new_atoms = AtomList()
            for mol in owner:
                for idx in sorted(mol):
                    new_atoms.append(self.atoms[idx])
            self.atoms = new_atoms
            # Re-sort our residues and residue list for new atom ordering
            for residue in self.residues:
                residue.sort()
            self.residues.sort()
            return owner

        return None

    #===================================================

    def fill_LJ(self):
        """
        Calculates the Lennard-Jones parameters (Rmin/2 and epsilon) for each
        atom type by computing their values from the A and B coefficients of
        each atom interacting with itself.

        This fills the :attr:`LJ_radius`, :attr:`LJ_depth`, and :attr:`LJ_types`
        data structures.
        """
        self.LJ_radius = []  # empty LJ_radii so it can be re-filled
        self.LJ_depth = []   # empty LJ_depths so it can be re-filled
        self.LJ_types = {}   # empty LJ_types so it can be re-filled
        one_sixth = 1 / 6    # we need to raise some numbers to the 1/6th power

        pd = self.parm_data
        acoef = pd['LENNARD_JONES_ACOEF']
        bcoef = pd['LENNARD_JONES_BCOEF']
        natom = self.pointers['NATOM']
        ntypes = self.pointers['NTYPES']
        for i in range(natom): # fill the LJ_types array
            self.LJ_types[pd["AMBER_ATOM_TYPE"][i]] = pd["ATOM_TYPE_INDEX"][i]

        for i in range(ntypes):
            lj_index = pd["NONBONDED_PARM_INDEX"][ntypes*i+i] - 1
            if lj_index < 0 or acoef[lj_index] < 1.0e-10 or bcoef[lj_index] < 1.0e-10:
                self.LJ_radius.append(0)
                self.LJ_depth.append(0)
            else:
                factor = 2 * acoef[lj_index] / bcoef[lj_index]
                self.LJ_radius.append(pow(factor, one_sixth) * 0.5)
                self.LJ_depth.append(bcoef[lj_index] / 2 / factor)

    #===================================================

    def recalculate_LJ(self):
        """
        Fills the ``LENNARD_JONES_ACOEF`` and ``LENNARD_JONES_BCOEF`` arrays in
        the :attr:`parm_data` raw data dictionary by applying the canonical
        Lorentz-Berthelot combining rules to the values in :attr:`LJ_radius` and
        :attr:`LJ_depth`.

        Notes
        -----
        This will undo any off-diagonal L-J modifications you may have made, so
        call this function with care.
        """
        assert self.combining_rule in ('lorentz', 'geometric'), "Unrecognized combining rule"
        if self.combining_rule == 'lorentz':
            comb_sig = lambda sig1, sig2: 0.5 * (sig1 + sig2)
        elif self.combining_rule == 'geometric':
            comb_sig = lambda sig1, sig2: sqrt(sig1 * sig2)
        pd = self.parm_data
        ntypes = self.pointers['NTYPES']
        fac = 2**(-1/6) * 2
        LJ_sigma = [x*fac for x in self.LJ_radius]
        fac = 2**(1/6)
        for i in range(ntypes):
            for j in range(i, ntypes):
                index = pd['NONBONDED_PARM_INDEX'][ntypes*i+j] - 1
                if index < 0:
                    # This indicates *either* a 10-12 potential for this pair
                    # _or_ it indicates a placeholder for HW atom type
                    # interactions and serves as a tag for fast water routines
                    # (or did in the past, anyway). So just skip to the next one
                    continue
                rij = comb_sig(LJ_sigma[i], LJ_sigma[j]) * fac
                wdij = sqrt(self.LJ_depth[i] * self.LJ_depth[j])
                pd["LENNARD_JONES_ACOEF"][index] = wdij * rij**12
                pd["LENNARD_JONES_BCOEF"][index] = 2 * wdij * rij**6

    #===================================================

    def has_NBFIX(self):
        """
        This routine determines whether there are any off-diagonal Lennard-Jones
        modifications (i.e., if any two atoms have a L-J pair interaction that
        does not match the combined L-J parameters for that pair).

        Returns
        -------
        nbfix : bool
            If True, off-diagonal elements in the combined Lennard-Jones matrix
            exist. If False, they do not.
        """
        assert self.combining_rule in ('lorentz', 'geometric'), "Unrecognized combining rule"
        if self.combining_rule == 'lorentz':
            comb_sig = lambda sig1, sig2: 0.5 * (sig1 + sig2)
        elif self.combining_rule == 'geometric':
            comb_sig = lambda sig1, sig2: sqrt(sig1 * sig2)
        fac = 2**(-1/6) * 2
        LJ_sigma = [x*fac for x in self.LJ_radius]
        pd = self.parm_data
        ntypes = self.parm_data['POINTERS'][PrmtopPointers.NTYPES]
        fac = 2**(1/6)
        for i in range(ntypes):
            for j in range(ntypes):
                idx = pd['NONBONDED_PARM_INDEX'][ntypes*i+j] - 1
                if idx < 0: continue
                rij = comb_sig(LJ_sigma[i], LJ_sigma[j]) * fac
                wdij = sqrt(self.LJ_depth[i] * self.LJ_depth[j])
                a = pd['LENNARD_JONES_ACOEF'][idx]
                b = pd['LENNARD_JONES_BCOEF'][idx]
                if a == 0 or b == 0:
                    if a != 0 or b != 0 or (wdij != 0 and rij != 0):
                        return True
                elif (abs((a - (wdij * rij**12)) / a) > 1e-6 or
                        abs((b - (2 * wdij * rij**6)) / b) > 1e-6):
                    return True
        return False

    #===================================================

    def has_1012(self):
        """
        This routine determines whether there are any defined 10-12
        Lennard-Jones interactions that are non-zero

        Returns
        -------
        has_10_12 : bool
            If True, 10-12 interactions *are* defined for this particular system
        """
        indices_with_1012 = []
        ntypes = self.parm_data['POINTERS'][PrmtopPointers.NTYPES]
        for i in range(ntypes):
            for j in range(ntypes):
                idx = self.parm_data['NONBONDED_PARM_INDEX'][i*ntypes+j] - 1
                if idx >= 0: continue
                # It was negative, so we should have ADDED 1 to adjust for
                # indexing from 0
                idx = -idx - 2
                a = self.parm_data['HBOND_ACOEF'][idx]
                b = self.parm_data['HBOND_BCOEF'][idx]
                if a == 0 and b == 0: continue
                indices_with_1012.append((i, j))
        if not indices_with_1012: return False
        # Now make sure that some of the atoms *have* those indices
        active_indices = set()
        for atom in self.atoms: active_indices.add(atom.nb_idx-1)
        for i, j in indices_with_1012:
            if i in active_indices and j in active_indices:
                return True
        return False

    #===================================================

    def load_rst7(self, rst7):
        """ Loads coordinates into the AmberParm class

        Parameters
        ----------
        rst7 : str or :class:`Rst7`
            The Amber coordinate file (either ASCII restart or NetCDF restart)
            object or filename to assign atomic coordinates from.
        """
        if not hasattr(rst7, 'coordinates'):
            rst7 = Rst7.open(rst7)
        self.coordinates = rst7.coordinates
        self.hasvels = rst7.hasvels
        self.box = _copy.copy(rst7.box)
        self.hasbox = self.box is not None
        if self.hasvels:
            self.velocities = rst7.vels

    #===================================================

    @needs_openmm
    def omm_nonbonded_force(self, nonbondedMethod=None,
                            nonbondedCutoff=8*u.angstroms,
                            switchDistance=0*u.angstroms,
                            ewaldErrorTolerance=0.0005,
                            reactionFieldDielectric=78.5):
        """
        Creates the OpenMM NonbondedForce (and CustomNonbondedForce if
        necessary) to define the nonbonded interatomic potential for this
        system. A CustomNonbondedForce is used for the r^-4 part of the 12-6-4
        Lennard-Jones potential as well as any modified off-diagonal (i.e.,
        NBFIX) terms

        Parameters
        ----------
        nonbondedMethod : cutoff method
            This is the cutoff method. It can be either the NoCutoff,
            CutoffNonPeriodic, CutoffPeriodic, PME, or Ewald objects from the
            openmm.app namespace
        nonbondedCutoff : float or distance Quantity
            The nonbonded cutoff must be either a floating point number
            (interpreted as nanometers) or a Quantity with attached units. This
            is ignored if nonbondedMethod is NoCutoff.
        switchDistance : float or distance Quantity
            The distance at which the switching function is turned on for van
            der Waals interactions. This is ignored when no cutoff is used, and
            no switch is used if switchDistance is 0, negative, or greater than
            the cutoff
        ewaldErrorTolerance : float=0.0005
            When using PME or Ewald, the Ewald parameters will be calculated
            from this value
        reactionFieldDielectric : float=78.5
            If the nonbondedMethod is CutoffPeriodic or CutoffNonPeriodic, the
            region beyond the cutoff is treated using a reaction field method
            with this dielectric constant. It should be set to 1 if another
            implicit solvent model is being used (e.g., GB)

        Returns
        -------
        NonbondedForce [, CustomNonbondedForce]
            If a CustomNonbondedForce is necessary, the return value is a
            2-element tuple of NonbondedForce, CustomNonbondedForce. If only a
            NonbondedForce is necessary, that is the return value
        """
        if not self.atoms: return None
        nonbfrc = super(AmberParm, self).omm_nonbonded_force(
            nonbondedMethod, nonbondedCutoff, switchDistance,
            ewaldErrorTolerance, reactionFieldDielectric
        )
        has1012 = self.has_1012()
        hasnbfix = self.has_NBFIX()
        has1264 = 'LENNARD_JONES_CCOEF' in self.flag_list
        if not hasnbfix and not has1264 and not has1012:
            if self.chamber:
                self._modify_nonb_exceptions(nonbfrc, None)
            return nonbfrc

        # If we have NBFIX, omm_nonbonded_force returned a tuple
        if hasnbfix:
            nonbfrc = nonbfrc[0]

        # If we have a 10-12 potential, toggle hasnbfix since the 12-6 terms
        # need to be handled via a lookup table (to avoid counting *both* 10-12
        # and 12-6 terms for the same pairs). Do this now, since nonbfrc is only
        # a tuple if hasnbfix is True *without* 10-12 detection.
        hasnbfix = hasnbfix or has1012

        # We need a CustomNonbondedForce... determine what it needs to calculate
        if hasnbfix and has1264:
            if has1012:
                force = mm.CustomNonbondedForce(
                    '(a/r6)^2-b/r6-c/r4+(ah/r6)^2-bh/(r6*r4); r6=r4*r2;'
                    'r4=r2^2; r2=r^2;'
                    'a=acoef(type1, type2);'
                    'b=bcoef(type1, type2);'
                    'c=ccoef(type1, type2);'
                    'ah=ahcoef(type1, type2);'
                    'bh=bhcoef(type1, type2);'
                )
            else:
                force = mm.CustomNonbondedForce(
                    '(a/r6)^2-b/r6-c/r4; r6=r4*r2;'
                    'r4=r2^2; r2=r^2;'
                    'a=acoef(type1, type2);'
                    'b=bcoef(type1, type2);'
                    'c=ccoef(type1, type2);'
                )
        elif hasnbfix:
            if has1012:
                force = mm.CustomNonbondedForce(
                    '(a/r6)^2-b/r6+(ah/r6)^2-bh/(r6*r4);'
                    'r6=r4*r2;r4=r2^2;r2=r^2;'
                    'a=acoef(type1, type2);'
                    'b=bcoef(type1, type2);'
                    'ah=ahcoef(type1, type2);'
                    'bh=bhcoef(type1, type2);'
                )
            else:
                force = mm.CustomNonbondedForce(
                    '(a/r6)^2-b/r6; r6=r2*r2*r2;'
                    'r2=r^2; a=acoef(type1, type2);'
                    'b=bcoef(type1, type2);'
                )
        elif has1264:
            force = mm.CustomNonbondedForce('-c/r^4;c=ccoef(type1, type2);')

        # Set up the force with all of the particles
        force.addPerParticleParameter('type')
        force.setForceGroup(self.NONBONDED_FORCE_GROUP)
        for atom in self.atoms:
            force.addParticle([atom.nb_idx-1])

        # Now construct the lookup tables
        ene_conv = u.kilocalories.conversion_factor_to(u.kilojoules)
        length_conv = u.angstroms.conversion_factor_to(u.nanometers)
        ntypes = self.parm_data['POINTERS'][PrmtopPointers.NTYPES]
        if hasnbfix:
            acoef = [0 for i in range(ntypes*ntypes)]
            parm_acoef = self.parm_data['LENNARD_JONES_ACOEF']
            bcoef = acoef[:]
            parm_bcoef = self.parm_data['LENNARD_JONES_BCOEF']
            if has1264:
                ccoef = acoef[:]
                parm_ccoef = self.parm_data['LENNARD_JONES_CCOEF']
            if has1012:
                ahcoef = acoef[:]
                bhcoef = acoef[:]
                parm_ahcoef = self.parm_data['HBOND_ACOEF']
                parm_bhcoef = self.parm_data['HBOND_BCOEF']
            afac = sqrt(ene_conv) * length_conv**6
            bfac = ene_conv * length_conv**6
            cfac = ene_conv * length_conv**4
            ahfac = sqrt(ene_conv) * length_conv**6
            bhfac = ene_conv * length_conv**10
            nbidx = self.parm_data['NONBONDED_PARM_INDEX']
            for i in range(ntypes):
                for j in range(ntypes):
                    idx = nbidx[ntypes*i+j] - 1
                    ii = i + ntypes * j
                    if idx < 0 and has1012:
                        idx = -idx - 2 # Fix adjust for 0 since it was negative
                        ahcoef[ii] = sqrt(parm_ahcoef[idx]) * ahfac
                        bhcoef[ii] = parm_bhcoef[idx] * bhfac
                        if has1264:
                            ccoef[ii] = parm_ccoef[idx] * cfac
                    elif idx >= 0:
                        acoef[ii] = sqrt(parm_acoef[idx]) * afac
                        bcoef[ii] = parm_bcoef[idx] * bfac
                        if has1264:
                            ccoef[ii] = parm_ccoef[idx] * cfac
            force.addTabulatedFunction('acoef', mm.Discrete2DFunction(ntypes, ntypes, acoef))
            force.addTabulatedFunction('bcoef', mm.Discrete2DFunction(ntypes, ntypes, bcoef))
            if has1264:
                force.addTabulatedFunction('ccoef', mm.Discrete2DFunction(ntypes, ntypes, ccoef))
            # Our CustomNonbondedForce is taking care of the LJ part of our
            # potential, so we need to go through nonbfrc and zero-out the
            # Lennard-Jones parameters, but keep the charge parameters in place.
            for i in range(nonbfrc.getNumParticles()):
                chg, sig, eps = nonbfrc.getParticleParameters(i)
                nonbfrc.setParticleParameters(i, chg, 0.5, 0.0)
            # We still let the NonbondedForce handle our nonbonded exceptions,
            # but we may have to modify the Lennard-Jones part with
            # off-diagonal modifications. Offload this to a private method so it
            # can be overridden in the ChamberParm class
            self._modify_nonb_exceptions(nonbfrc, force)
        elif has1264:
            # Here we have JUST the r^-4 or r^-10/r^-12 parts, since the
            # hasnbfix block above handled the "hasnbfix and has1264" case
            if has1264:
                ccoef = [0 for i in range(ntypes*ntypes)]
                parm_ccoef = self.parm_data['LENNARD_JONES_CCOEF']
            cfac = ene_conv * length_conv**4
            ahfac = ene_conv * length_conv**12
            bhfac = ene_conv * length_conv**10
            nbidx = self.parm_data['NONBONDED_PARM_INDEX']
            for i in range(ntypes):
                for j in range(ntypes):
                    idx = nbidx[ntypes*i+j] - 1
                    if idx < 0:
                        continue
                    ccoef[i+ntypes*j] = parm_ccoef[idx] * cfac
            force.addTabulatedFunction('ccoef', mm.Discrete2DFunction(ntypes, ntypes, ccoef))
            # Copy the exclusions
            for ii in range(nonbfrc.getNumExceptions()):
                i, j, qq, ss, ee = nonbfrc.getExceptionParameters(ii)
                force.addExclusion(i, j)
        if has1012:
            force.addTabulatedFunction('ahcoef', mm.Discrete2DFunction(ntypes, ntypes, ahcoef))
            force.addTabulatedFunction('bhcoef', mm.Discrete2DFunction(ntypes, ntypes, bhcoef))
        # Copy the switching function information to the CustomNonbondedForce
        if nonbfrc.getUseSwitchingFunction():
            force.setUseSwitchingFunction(True)
            force.setSwitchingDistance(nonbfrc.getSwitchingDistance())
        # Set the dispersion correction on (by default)
        force.setUseLongRangeCorrection(True)
        # Determine which nonbonded method we should use and transfer the
        # nonbonded cutoff
        assert nonbondedMethod in {app.NoCutoff, app.CutoffNonPeriodic, app.PME, app.LJPME, app.Ewald, app.CutoffPeriodic}, 'Bad nonbondedMethod'
        if nonbondedMethod is app.NoCutoff:
            force.setNonbondedMethod(mm.CustomNonbondedForce.NoCutoff)
        elif nonbondedMethod is app.CutoffNonPeriodic:
            force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic)
        elif nonbondedMethod in (app.PME, app.Ewald, app.CutoffPeriodic, app.LJPME):
            force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        force.setCutoffDistance(nonbfrc.getCutoffDistance())

        return nonbfrc, force

    #===================================================

    # Iterators for parameters with and without hydrogen
    @property
    def bonds_inc_h(self):
        """ All bonds including hydrogen """
        for bond in self.bonds:
            if bond.atom1.atomic_number == 1 or bond.atom2.atomic_number == 1:
                yield bond

    @property
    def bonds_without_h(self):
        """ All bonds without hydrogen """
        for bond in self.bonds:
            if bond.atom1.atomic_number == 1 or bond.atom2.atomic_number == 1:
                continue
            yield bond

    @property
    def angles_inc_h(self):
        """ All angles including hydrogen """
        for angle in self.angles:
            if (angle.atom1.atomic_number == 1 or angle.atom2.atomic_number == 1
                    or angle.atom3.atomic_number == 1):
                yield angle

    @property
    def angles_without_h(self):
        """ All angles including hydrogen """
        for angle in self.angles:
            if (angle.atom1.atomic_number == 1 or angle.atom2.atomic_number == 1
                    or angle.atom3.atomic_number == 1):
                continue
            yield angle

    @property
    def dihedrals_inc_h(self):
        """ All dihedrals including hydrogen """
        for dihed in self.dihedrals:
            if (dihed.atom1.atomic_number == 1
                    or dihed.atom2.atomic_number == 1
                    or dihed.atom3.atomic_number == 1
                    or dihed.atom4.atomic_number == 1):
                yield dihed

    @property
    def dihedrals_without_h(self):
        """ All dihedrals including hydrogen """
        for dihed in self.dihedrals:
            if (dihed.atom1.atomic_number == 1
                    or dihed.atom2.atomic_number == 1
                    or dihed.atom3.atomic_number == 1
                    or dihed.atom4.atomic_number == 1):
                continue
            yield dihed

    #===================================================

    @property
    def chamber(self):
        """ Whether this instance uses the CHARMM force field """
        return False

    @property
    def amoeba(self):
        """ Whether this instance uses the Amoeba force field """
        return False

    @property
    def has_cmap(self):
        """ Whether this instance has correction map terms or not """
        return len(self.cmaps) > 0 or (f'{self._cmap_prefix}CMAP_COUNT') in self.parm_data


    #===========  PRIVATE INSTANCE METHODS  ============

    def _truncate_array(self, section, length):
        """ Truncates an array to get the given length """
        self.parm_data[section] = self.parm_data[section][:length]

    #===================================================

    def _load_cmap_info(self):
        """ Loads the CHARMM CMAP types and array """
        if not self.has_cmap: return
        del self.cmaps[:]
        del self.cmap_types[:]
        resolution_key = f"{self._cmap_prefix}CMAP_RESOLUTION"
        parameter_key = f"{self._cmap_prefix}CMAP_PARAMETER_{{:02d}}"
        for i in range(self.pointers['CMAP_TYPES']):
            resolution = self.parm_data[resolution_key][i]
            grid = self.parm_data[parameter_key.format(i + 1)]
            cmts = self.parm_comments[parameter_key.format(i + 1)]
            self.cmap_types.append(CmapType(resolution, grid, cmts, list=self.cmap_types))
        it = iter(self.parm_data[self._cmap_prefix + 'CMAP_INDEX'])
        for i, j, k, l, m, n in zip(it, it, it, it, it, it):
            self.cmaps.append(
                Cmap(self.atoms[i-1], self.atoms[j-1], self.atoms[k-1],
                     self.atoms[l-1], self.atoms[m-1], self.cmap_types[n-1])
            )

    #===================================================

    def _check_section_lengths(self):
        """
        Checks that all of the raw sections have the appropriate length as
        specified by the POINTER section.

        Raises
        ------
        AmberError if any of the lengths are incorrect
        """
        def check_length(key, length, required=True):
            if not required and key not in self.parm_data:
                return
            if len(self.parm_data[key]) != length:
                raise AmberError(f'FLAG {key} has {len(self.parm_data[key])} elements; expected {length}')
        natom = self.ptr('NATOM')
        check_length('ATOM_NAME', natom)
        check_length('CHARGE', natom)
        check_length('MASS', natom)
        check_length('ATOM_TYPE_INDEX', natom)
        check_length('NUMBER_EXCLUDED_ATOMS', natom)
        check_length('JOIN_ARRAY', natom)
        check_length('IROTAT', natom)
        check_length('RADIUS', natom, False)
        check_length('SCREEN', natom, False)
        check_length('ATOMIC_NUMBER', natom, False)

        ntypes = self.ptr('NTYPES')
        check_length('NONBONDED_PARM_INDEX', ntypes*ntypes)
        check_length('LENNARD_JONES_ACOEF', ntypes*(ntypes+1)//2)
        check_length('LENNARD_JONES_BCOEF', ntypes*(ntypes+1)//2)
        check_length('LENNARD_JONES_CCOEF', ntypes*(ntypes+1)//2, False)

        nres = self.ptr('NRES')
        check_length('RESIDUE_LABEL', nres)
        check_length('RESIDUE_POINTER', nres)
        check_length('RESIDUE_CHAINID', nres, False)
        check_length('RESIDUE_ICODE', nres, False)
        check_length('RESIDUE_NUMBER', nres, False)

        check_length('BOND_FORCE_CONSTANT', self.ptr('NUMBND'))
        check_length('BOND_EQUIL_VALUE', self.ptr('NUMBND'))
        check_length('ANGLE_FORCE_CONSTANT', self.ptr('NUMANG'))
        check_length('ANGLE_EQUIL_VALUE', self.ptr('NUMANG'))
        check_length('DIHEDRAL_FORCE_CONSTANT', self.ptr('NPTRA'))
        check_length('DIHEDRAL_PERIODICITY', self.ptr('NPTRA'))
        check_length('DIHEDRAL_PHASE', self.ptr('NPTRA'))
        check_length('SCEE_SCALE_FACTOR', self.ptr('NPTRA'), False)
        check_length('SCNB_SCALE_FACTOR', self.ptr('NPTRA'), False)
        check_length('SOLTY', self.ptr('NATYP'))
        check_length('BONDS_INC_HYDROGEN', self.ptr('NBONH')*3)
        check_length('BONDS_WITHOUT_HYDROGEN', self.ptr('MBONA')*3)
        check_length('ANGLES_INC_HYDROGEN', self.ptr('NTHETH')*4)
        check_length('ANGLES_WITHOUT_HYDROGEN', self.ptr('NTHETA')*4)
        check_length('DIHEDRALS_INC_HYDROGEN', self.ptr('NPHIH')*5)
        check_length('DIHEDRALS_WITHOUT_HYDROGEN', self.ptr('NPHIA')*5)
        check_length('HBOND_ACOEF', self.ptr('NPHB'))
        check_length('HBOND_BCOEF', self.ptr('NPHB'))
        check_length('SOLVENT_POINTERS', 3, False)
        if 'SOLVENT_POINTERS' in self.parm_data:
            check_length('ATOMS_PER_MOLECULE', self.parm_data['SOLVENT_POINTERS'][1], False)
        if self.has_cmap:
            check_length(self._cmap_prefix + 'CMAP_COUNT', 2)
            check_length(self._cmap_prefix + 'CMAP_RESOLUTION', self.pointers['CMAP_TYPES'])
            resolution_key = self._cmap_prefix + 'CMAP_RESOLUTION'
            parameter_key = self._cmap_prefix + 'CMAP_PARAMETER_%02d'
            for i in range(self.pointers['CMAP_TYPES']):
                res = self.parm_data[resolution_key][i]
                check_length(parameter_key % (i+1), res*res)

    #===================================================

    def _load_atoms_and_residues(self):
        """
        Loads the atoms and residues (which are always done together) into the
        data structure
        """
        del self.residues[:]
        del self.atoms[:]
        # Figure out on which atoms the residues start and stop
        natom = self.parm_data['POINTERS'][PrmtopPointers.NATOM]
        res_ptr = self.parm_data['RESIDUE_POINTER'] + [natom+1]
        try:
            atnums = self.parm_data['ATOMIC_NUMBER']
        except KeyError:
            atnums = None
        try:
            res_icd = self.parm_data['RESIDUE_ICODE']
        except KeyError:
            res_icd = ['' for i in range(self.parm_data['POINTERS'][PrmtopPointers.NRES])]
        try:
            res_chn = self.parm_data['RESIDUE_CHAINID']
        except KeyError:
            res_chn = ['' for i in range(self.parm_data['POINTERS'][PrmtopPointers.NRES])]
        for i, resname in enumerate(self.parm_data['RESIDUE_LABEL']):
            resstart = res_ptr[i] - 1
            resend = res_ptr[i+1] - 1
            for j in range(resstart, resend):
                if atnums is None:
                    if self.parm_data['AMBER_ATOM_TYPE'][j] in ('EP', 'LP'):
                        atom = ExtraPoint()
                    else:
                        atom = Atom()
                elif atnums[j] == 0:
                    atom = ExtraPoint()
                else:
                    atom = Atom()
                self.add_atom(atom, resname, i, res_chn[i], res_icd[i])
        # Residue number may not be unique, since zeros may be assigned to all
        # solvent residues by addPDB. So we can't use the RESIDUE_NUMBER section
        # to fill in the residue sequence IDs in the add_atom call above. Add it
        # in as a post-hoc addition now if that information is present
        if 'RESIDUE_NUMBER' in self.parm_data:
            for res, num in zip(self.residues, self.parm_data['RESIDUE_NUMBER']):
                res.number = num

    #===================================================

    def _load_extra_exclusions(self):
        """
        Look through the exclusion list in the prmtop file and see if any
        _additional_ exclusions outside the basic ones defined for bonds,
        angles, and dihedrals are specified. If so, add those to the exclusion
        list.

        This also goes through all atoms and loads the proper bond, angle and
        dihedral partners into any extra points. The way extra points are
        handled is that they are considered to carry the same topological
        connectivity as they actual atom they are bonded to.
        """
        num_excluded = self.parm_data['NUMBER_EXCLUDED_ATOMS']
        excluded_list = self.parm_data['EXCLUDED_ATOMS_LIST']
        first_excl = 0
        for i, atom in enumerate(self.atoms):
            exclusions = set()
            bond_excl = set(atom._bond_partners + atom._angle_partners +
                            atom._dihedral_partners + atom._tortor_partners +
                            atom._exclusion_partners)
            nexcl = num_excluded[i]
            for j in range(first_excl, first_excl+nexcl):
                if excluded_list[j]: # Zeros are placeholders
                    exclusions.add(self.atoms[excluded_list[j]-1])
            for eatom in exclusions - bond_excl:
                atom.exclude(eatom)
            first_excl += nexcl

    #===================================================

    def _load_bond_info(self):
        """ Loads the bond types and bond arrays """
        del self.bond_types[:]
        del self.bonds[:]
        for k, req in zip(self.parm_data['BOND_FORCE_CONSTANT'], self.parm_data['BOND_EQUIL_VALUE']):
            self.bond_types.append(BondType(k, req, self.bond_types))
        it = iter(self.parm_data['BONDS_WITHOUT_HYDROGEN'])
        for i, j, k in zip(it, it, it):
            self.bonds.append(
                Bond(self.atoms[i//3], self.atoms[j//3], self.bond_types[k-1])
            )
        it = iter(self.parm_data['BONDS_INC_HYDROGEN'])
        for i, j, k in zip(it, it, it):
            self.bonds.append(
                Bond(self.atoms[i//3], self.atoms[j//3], self.bond_types[k-1])
            )

    #===================================================

    def _load_angle_info(self):
        """ Loads the angle types and angle arrays """
        del self.angle_types[:]
        del self.angles[:]
        for k, theteq in zip(
            self.parm_data['ANGLE_FORCE_CONSTANT'], self.parm_data['ANGLE_EQUIL_VALUE']
        ):
            theteq *= RAD_TO_DEG
            self.angle_types.append(AngleType(k, theteq, self.angle_types))
        it = iter(self.parm_data['ANGLES_WITHOUT_HYDROGEN'])
        for i, j, k, l in zip(it, it, it, it):
            self.angles.append(
                Angle(self.atoms[i//3], self.atoms[j//3], self.atoms[k//3], self.angle_types[l-1])
            )
        it = iter(self.parm_data['ANGLES_INC_HYDROGEN'])
        for i, j, k, l in zip(it, it, it, it):
            self.angles.append(
                Angle(self.atoms[i//3], self.atoms[j//3], self.atoms[k//3], self.angle_types[l-1])
            )

    #===================================================

    def _load_dihedral_info(self):
        """ Loads the dihedral types and dihedral arrays """
        del self.dihedral_types[:]
        del self.dihedrals[:]
        try:
            scee = self.parm_data['SCEE_SCALE_FACTOR']
        except KeyError:
            scee = [1.2 for i in self.parm_data['DIHEDRAL_FORCE_CONSTANT']]
        try:
            scnb = self.parm_data['SCNB_SCALE_FACTOR']
        except KeyError:
            scnb = [2.0 for i in self.parm_data['DIHEDRAL_FORCE_CONSTANT']]
        for k, per, ph, e, n in zip(
            self.parm_data['DIHEDRAL_FORCE_CONSTANT'],
            self.parm_data['DIHEDRAL_PERIODICITY'],
            self.parm_data['DIHEDRAL_PHASE'],
            scee,
            scnb,
        ):
            ph *= RAD_TO_DEG
            self.dihedral_types.append(
                DihedralType(k, per, ph, e, n, list=self.dihedral_types)
            )
        it = iter(self.parm_data['DIHEDRALS_WITHOUT_HYDROGEN'])
        for i, j, k, l, m in zip(it, it, it, it, it):
            ignore_end = k < 0
            improper = l < 0
            self.dihedrals.append(
                Dihedral(self.atoms[i//3], self.atoms[j//3],
                         self.atoms[abs(k)//3], self.atoms[abs(l)//3],
                         improper=improper, ignore_end=ignore_end,
                         type=self.dihedral_types[m-1])
            )
        it = iter(self.parm_data['DIHEDRALS_INC_HYDROGEN'])
        for i, j, k, l, m in zip(it, it, it, it, it):
            ignore_end = k < 0
            improper = l < 0
            self.dihedrals.append(
                Dihedral(
                    self.atoms[i//3], self.atoms[j//3],
                    self.atoms[abs(k)//3], self.atoms[abs(l)//3],
                    improper=improper, ignore_end=ignore_end,
                    type=self.dihedral_types[m-1])
            )

    #===================================================

    def _xfer_atom_info(self):
        """
        Sets the various topology file section data from the `atoms` list to the
        topology file data in `parm_data`
        """
        natom = len(self.atoms)
        data = self.parm_data
        data['POINTERS'][PrmtopPointers.NATOM] = natom
        self.pointers['NATOM'] = natom
        data['ATOM_NAME'] = [atom.name[:4] for atom in self.atoms]
        data['AMBER_ATOM_TYPE'] = [atom.type[:4] for atom in self.atoms]
        data['CHARGE'] = [atom.charge for atom in self.atoms]
        data['MASS'] = [atom.mass for atom in self.atoms]
        data['ATOM_TYPE_INDEX'] = [atom.nb_idx for atom in self.atoms]
        data['JOIN_ARRAY'] = [atom.join for atom in self.atoms]
        data['TREE_CHAIN_CLASSIFICATION'] = [atom.tree[:4] for atom in self.atoms]
        data['IROTAT'] = [atom.irotat for atom in self.atoms]
        data['NUMBER_EXCLUDED_ATOMS'] = [0 for atom in self.atoms]
        if 'RADII' in data:
            data['RADII'] = [atom.solvent_radius for atom in self.atoms]
        if 'SCREEN' in data:
            data['SCREEN'] = [atom.screen for atom in self.atoms]
        if 'ATOMIC_NUMBER' in data:
            data['ATOMIC_NUMBER'] = [atom.atomic_number for atom in self.atoms]
        # Do the non-bonded exclusions now
        data['EXCLUDED_ATOMS_LIST'] = []
        nextra = 0
        max_typ = 0
        for i, atom in enumerate(self.atoms):
            excl = atom.nonbonded_exclusions(index_from=1)
            if len(excl) == 0:
                excl = [0]
            data['EXCLUDED_ATOMS_LIST'] += excl
            data['NUMBER_EXCLUDED_ATOMS'][i] = len(excl)
            if atom.atomic_number == 0:
                nextra += 1
            max_typ = max(max_typ, atom.nb_idx)
        nnb = len(data['EXCLUDED_ATOMS_LIST'])
        data['POINTERS'][PrmtopPointers.NNB] = nnb
        self.pointers['NNB'] = self.pointers['NEXT'] = nnb
        data['POINTERS'][PrmtopPointers.NUMEXTRA] = nextra
        self.pointers['NUMEXTRA'] = nextra
        max_typ = max(data['POINTERS'][PrmtopPointers.NTYPES], max_typ)
        data['POINTERS'][PrmtopPointers.NTYPES] = max_typ
        self.pointers['NTYPES'] = max_typ

    #===================================================

    def _xfer_residue_info(self):
        """
        Sets the various topology file section data from the `residues` list to
        the topology file data in `parm_data`
        """
        data = self.parm_data
        nres = len(self.residues)
        data['POINTERS'][PrmtopPointers.NRES] = nres
        self.pointers['NRES'] = nres
        data['RESIDUE_LABEL'] = [r.name[:4] for r in self.residues]
        data['RESIDUE_POINTER'] = [r.atoms[0].idx+1 for r in self.residues]
        if 'RESIDUE_NUMBER' in data:
            data['RESIDUE_NUMBER'] = [r.number for r in self.residues]
        if 'RESIDUE_CHAINID' in data:
            data['RESIDUE_CHAINID'] = [res.chain for res in self.residues]
        if 'RESIDUE_ICODE' in data:
            data['RESIDUE_ICODE'] = [r.insertion_code for r in self.residues]
        nmxrs = max([len(res) for res in self.residues]) if self.residues else 0
        data['POINTERS'][PrmtopPointers.NMXRS] = nmxrs
        self.pointers['NMXRS'] = nmxrs

    #===================================================

    def _xfer_bond_info(self):
        """
        Sets the data for the various bond arrays in the raw data from the
        parameter lists
        """
        # First do the bond types
        data = self.parm_data
        for bond_type in self.bond_types:
            bond_type.used = False
        for bond in self.bonds:
            bond.type.used = True
        self.bond_types.prune_unused()
        data['BOND_FORCE_CONSTANT'] = [type.k for type in self.bond_types]
        data['BOND_EQUIL_VALUE'] = [type.req for type in self.bond_types]
        data['POINTERS'][PrmtopPointers.NUMBND] = len(self.bond_types)
        self.pointers['NUMBND'] = len(self.bond_types)
        # Now do the bond arrays
        data['BONDS_INC_HYDROGEN'] = bond_array = []
        bond_list = list(self.bonds_inc_h)
        for bond in bond_list:
            bond_array.extend([bond.atom1.idx*3, bond.atom2.idx*3, bond.type.idx+1])
        data['POINTERS'][PrmtopPointers.NBONH] = len(bond_list)
        self.pointers['NBONH'] = len(bond_list)
        data['BONDS_WITHOUT_HYDROGEN'] = bond_array = []
        bond_list = list(self.bonds_without_h)
        for bond in bond_list:
            bond_array.extend([bond.atom1.idx*3, bond.atom2.idx*3, bond.type.idx+1])
        data['POINTERS'][PrmtopPointers.MBONA] = len(bond_list)
        data['POINTERS'][PrmtopPointers.NBONA] = len(bond_list)
        self.pointers['MBONA'] = len(bond_list)
        self.pointers['NBONA'] = len(bond_list)

    #===================================================

    def _xfer_angle_info(self):
        """
        Sets the data for the various angle arrays in the raw data from the
        parameter lists
        """
        # First do the angle types
        data = self.parm_data
        for angle_type in self.angle_types:
            angle_type.used = False
        for angle in self.angles:
            angle.type.used = True
        self.angle_types.prune_unused()
        data['ANGLE_FORCE_CONSTANT'] = [type.k for type in self.angle_types]
        data['ANGLE_EQUIL_VALUE'] = [type.theteq*DEG_TO_RAD for type in self.angle_types]
        data['POINTERS'][PrmtopPointers.NUMANG] = len(self.angle_types)
        self.pointers['NUMANG'] = len(self.angle_types)
        # Now do the angle arrays
        data['ANGLES_INC_HYDROGEN'] = angle_array = []
        angle_list = list(self.angles_inc_h)
        for angle in angle_list:
            angle_array.extend(
                [angle.atom1.idx*3, angle.atom2.idx*3, angle.atom3.idx*3, angle.type.idx+1]
            )
        data['POINTERS'][PrmtopPointers.NTHETH] = len(angle_list)
        self.pointers['NTHETH'] = len(angle_list)
        data['ANGLES_WITHOUT_HYDROGEN'] = angle_array = []
        angle_list = list(self.angles_without_h)
        for angle in angle_list:
            angle_array.extend(
                [angle.atom1.idx*3, angle.atom2.idx*3, angle.atom3.idx*3, angle.type.idx+1]
            )
        data['POINTERS'][PrmtopPointers.NTHETA] = len(angle_list)
        data['POINTERS'][PrmtopPointers.MTHETA] = len(angle_list)
        self.pointers['NTHETA'] = len(angle_list)
        self.pointers['MTHETA'] = len(angle_list)

    #===================================================

    def _xfer_dihedral_info(self):
        """
        Sets the data for the various dihedral arrays in the raw data from the
        parameter lists
        """
        # First do the dihedral types
        data = self.parm_data
        for dihedral_type in self.dihedral_types:
            dihedral_type.used = False
        for dihed in self.dihedrals:
            dihed.type.used = True
        self.dihedral_types.prune_unused()
        data['DIHEDRAL_FORCE_CONSTANT'] = [type.phi_k for type in self.dihedral_types]
        data['DIHEDRAL_PERIODICITY'] = [type.per for type in self.dihedral_types]
        data['DIHEDRAL_PHASE'] = [type.phase*DEG_TO_RAD for type in self.dihedral_types]
        if 'SCEE_SCALE_FACTOR' in data:
            data['SCEE_SCALE_FACTOR'] = [type.scee for type in self.dihedral_types]
        if 'SCNB_SCALE_FACTOR' in data:
            data['SCNB_SCALE_FACTOR'] = [type.scnb for type in self.dihedral_types]
        data['POINTERS'][PrmtopPointers.NPTRA] = len(self.dihedral_types)
        self.pointers['NPTRA'] = len(self.dihedral_types)
        # Now do the dihedral arrays
        data['DIHEDRALS_INC_HYDROGEN'] = dihed_array = []
        dihed_list = list(self.dihedrals_inc_h)
        for dihed in dihed_list:
            imp_sign = -1 if dihed.improper else 1
            end_sign = -1 if dihed.ignore_end else 1
            if dihed.atom3.idx == 0 or dihed.atom4.idx == 0:
                dihed_array.extend(
                    [
                        dihed.atom4.idx *3,
                        dihed.atom3.idx *3,
                        dihed.atom2.idx *3 * end_sign,
                        dihed.atom1.idx *3 * imp_sign,
                        dihed.type.idx + 1,
                    ]
                )
            else:
                dihed_array.extend(
                    [
                        dihed.atom1.idx * 3,
                        dihed.atom2.idx * 3,
                        dihed.atom3.idx * 3 * end_sign,
                        dihed.atom4.idx * 3 * imp_sign,
                        dihed.type.idx + 1,
                    ]
                )
        data['POINTERS'][PrmtopPointers.NPHIH] = len(dihed_list)
        self.pointers['NPHIH'] = len(dihed_list)
        data['DIHEDRALS_WITHOUT_HYDROGEN'] = dihed_array = []
        dihed_list = list(self.dihedrals_without_h)
        for dihed in dihed_list:
            imp_sign = -1 if dihed.improper else 1
            end_sign = -1 if dihed.ignore_end else 1
            if dihed.atom3.idx == 0 or dihed.atom4.idx == 0:
                dihed_array.extend(
                    [
                        dihed.atom4.idx *3,
                        dihed.atom3.idx *3,
                        dihed.atom2.idx *3 * end_sign,
                        dihed.atom1.idx *3 * imp_sign,
                        dihed.type.idx + 1,
                    ]
                )
            else:
                dihed_array.extend(
                    [
                        dihed.atom1.idx * 3,
                        dihed.atom2.idx * 3,
                        dihed.atom3.idx * 3 * end_sign,
                        dihed.atom4.idx * 3 * imp_sign,
                        dihed.type.idx + 1,
                    ]
                )
        data['POINTERS'][PrmtopPointers.NPHIA] = len(dihed_list)
        data['POINTERS'][PrmtopPointers.MPHIA] = len(dihed_list)
        self.pointers['NPHIA'] = len(dihed_list)
        self.pointers['MPHIA'] = len(dihed_list)

    #===================================================

    def _xfer_cmap_properties(self):
        """ Sets the topology file section data from the cmap arrays """
        # If we have no cmaps, delete all remnants of the CMAP terms in the
        # prmtop and bail out
        if len(self.cmaps) == 0:
            # We have deleted all cmaps. Get rid of them from the parm file and
            # bail out. This is probably pretty unlikely, though...
            flag_prefix = f"{self._cmap_prefix}CMAP"
            flags_to_delete = [flag for flag in self.flag_list if flag.startswith(flag_prefix)]
            for flag in flags_to_delete:
                self.delete_flag(flag)
            if 'CMAP' in self.pointers:
                del self.pointers['CMAP']
            if 'CMAP_TYPES' in self.pointers:
                del self.pointers['CMAP_TYPES']
            return
        # Time to transfer our CMAP types
        data = self.parm_data
        for ct in self.cmap_types:
            ct.used = False
        for cmap in self.cmaps:
            cmap.type.used = True
        self.cmap_types.prune_unused()
        # All of our CMAP types are in different topology file sections. We need
        # to delete all of the CMAP_PARAMETER_XX sections and then
        # recreate them with the correct size and comments.  The comments have
        # been stored in the CMAP types themselves to prevent them from being
        # lost. We will also assume that the Fortran format we're going to use
        # is the same for all CMAP types, so just pull it from
        # CMAP_PARAMETER_01 (or fall back to 8(F9.5))
        parameter_key = self._cmap_prefix + 'CMAP_PARAMETER_%02d'
        try:
            fmt = str(self.formats[parameter_key % 1])
        except KeyError:
            fmt = '8(F9.5)'
        flags_to_delete = []
        for flag in self.flag_list:
            if 'CMAP_PARAMETER' in flag:
                flags_to_delete.append(flag)
        for flag in flags_to_delete:
            self.delete_flag(flag)
        # Now add them back
        after = self._cmap_prefix + 'CMAP_RESOLUTION'
        for i, ct in enumerate(self.cmap_types):
            newflag = self._cmap_prefix + 'CMAP_PARAMETER_%02d' % (i+1)
            self.add_flag(newflag, fmt, data=ct.grid, comments=ct.comments, after=after)
            after = newflag
        # Now do the CMAP_INDEX section
        data[self._cmap_prefix + 'CMAP_INDEX'] = cmap_array = []
        for cm in self.cmaps:
            cmap_array.extend([cm.atom1.idx+1, cm.atom2.idx+1, cm.atom3.idx+1, cm.atom4.idx+1,
                               cm.atom5.idx+1, cm.type.idx+1])
        data[self._cmap_prefix + 'CMAP_COUNT'] = [len(self.cmaps), len(self.cmap_types)]
        data[self._cmap_prefix + 'CMAP_RESOLUTION'] = [ct.resolution for ct in self.cmap_types]
        self.pointers['CMAP'] = len(self.cmaps)
        self.pointers['CMAP_TYPES'] = len(self.cmap_types)

    #===================================================

    def _add_standard_flags(self):
        """ Adds all of the standard flags to the parm_data array """
        self.set_version()
        self.add_flag('TITLE', '20a4', num_items=0)
        self.add_flag('POINTERS', '10I8', num_items=31)
        self.add_flag('ATOM_NAME', '20a4', num_items=0)
        self.add_flag('CHARGE', '5E16.8', num_items=0)
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
        self.add_flag('ANGLE_EQUIL_VALUE', '5E16.8', num_items=0)
        self.add_flag('DIHEDRAL_FORCE_CONSTANT', '5E16.8', num_items=0)
        self.add_flag('DIHEDRAL_PERIODICITY', '5E16.8', num_items=0)
        self.add_flag('DIHEDRAL_PHASE', '5E16.8', num_items=0)
        self.add_flag('SCEE_SCALE_FACTOR', '5E16.8', num_items=0)
        self.add_flag('SCNB_SCALE_FACTOR', '5E16.8', num_items=0)
        self.pointers['NATYP'] = 1
        self.parm_data['POINTERS'][PrmtopPointers.NATYP] = 1
        self.add_flag('SOLTY', '5E16.8', num_items=1)
        self.add_flag('LENNARD_JONES_ACOEF', '5E16.8', num_items=0)
        self.add_flag('LENNARD_JONES_BCOEF', '5E16.8', num_items=0)
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
        """
        Sets the tables of Lennard-Jones nonbonded interaction pairs
        """
        ntypes = self.parm_data['POINTERS'][PrmtopPointers.NTYPES]
        ntypes2 = ntypes * ntypes
        # Set up the index lookup tables (not a unique solution)
        self.parm_data['NONBONDED_PARM_INDEX'] = [0 for i in range(ntypes2)]
        holder = [0 for i in range(ntypes2)]
        idx = 0
        for i in range(ntypes):
            for j in range(i+1):
                idx += 1
                holder[ntypes*i+j] = holder[ntypes*j+i] = idx
        idx = 0
        for i in range(ntypes):
            for j in range(ntypes):
                self.parm_data['NONBONDED_PARM_INDEX'][idx] = holder[ntypes*i+j]
                idx += 1
        nttyp = ntypes * (ntypes + 1) // 2
        # Now build the Lennard-Jones arrays
        self.parm_data['LENNARD_JONES_ACOEF'] = [0 for i in range(nttyp)]
        self.parm_data['LENNARD_JONES_BCOEF'] = [0 for i in range(nttyp)]
        self.recalculate_LJ()
        # Now make any NBFIX modifications we had
        if nbfixes is not None:
            for i, fix in enumerate(nbfixes):
                for terms in fix:
                    j, rmin, eps, rmin14, eps14 = terms
                    i, j = min(i, j-1), max(i, j-1)
                    eps = abs(eps)
                    eps14 = abs(eps14)
                    idx = self.parm_data['NONBONDED_PARM_INDEX'][ntypes*i+j] - 1
                    self.parm_data['LENNARD_JONES_ACOEF'][idx] = eps * rmin**12
                    self.parm_data['LENNARD_JONES_BCOEF'][idx] = 2*eps * rmin**6

    #===================================================

    @needs_openmm
    def _modify_nonb_exceptions(self, nonbfrc, customforce):
        """
        Modifies the nonbonded force exceptions and the custom nonbonded force
        exclusions. The exceptions on the nonbonded force might need to be
        adjusted if off-diagonal modifications on the L-J matrix are present
        """
        # To get into this routine, either NBFIX is present OR this is a chamber
        # prmtop and we need to pull the 1-4 L-J parameters from the
        # LENNARD_JONES_14_A/BCOEF arrays
        length_conv = u.angstroms.conversion_factor_to(u.nanometers)
        ene_conv = u.kilocalories.conversion_factor_to(u.kilojoules)
        atoms = self.atoms
        try:
            acoef = self.parm_data['LENNARD_JONES_14_ACOEF']
            bcoef = self.parm_data['LENNARD_JONES_14_BCOEF']
        except KeyError:
            acoef = self.parm_data['LENNARD_JONES_ACOEF']
            bcoef = self.parm_data['LENNARD_JONES_BCOEF']
        nbidx = self.parm_data['NONBONDED_PARM_INDEX']
        ntypes = self.parm_data['POINTERS'][PrmtopPointers.NTYPES]
        sigma_scale = 2**(-1/6) * length_conv
        for ii in range(nonbfrc.getNumExceptions()):
            i, j, qq, ss, ee = nonbfrc.getExceptionParameters(ii)
            if qq.value_in_unit(u.elementary_charge**2) == 0 and (
                    ss.value_in_unit(u.angstroms) == 0 or
                    ee.value_in_unit(u.kilocalories_per_mole) == 0):
                # Copy this exclusion as-is... no need to modify the nonbfrc
                # exception parameters
                if customforce is not None:
                    customforce.addExclusion(i, j)
                continue
            # Figure out what the 1-4 scaling parameters were for this pair...
            unscaled_ee = sqrt(self.atoms[i].epsilon_14 * self.atoms[j].epsilon_14) * ene_conv
            try:
                one_scnb = ee.value_in_unit(u.kilojoules_per_mole) / unscaled_ee
            except ZeroDivisionError:
                one_scnb = 1
            id1 = atoms[i].nb_idx - 1
            id2 = atoms[j].nb_idx - 1
            idx = nbidx[ntypes*id1+id2] - 1
            if idx >= 0:
                a = acoef[idx]
                b = bcoef[idx]
            else:
                a = b = 0
            if b == 0:
                epsilon = 0.0
                sigma = 0.5
            elif idx >= 0:
                # b / a == 2 / r^6 --> (a / b * 2)^(1/6) = rmin
                rmin = (a / b * 2)**(1/6)
                epsilon = b / (2 * rmin**6) * ene_conv * one_scnb
                sigma = rmin * sigma_scale
            nonbfrc.setExceptionParameters(ii, i, j, qq, sigma, epsilon)
            if customforce is not None:
                customforce.addExclusion(i, j)

    #===================================================

    def _add_missing_13_14(self, ignore_inconsistent_vdw=False):
        """
        Uses the bond graph to fill in zero-parameter angles and dihedrals. The
        reason this is necessary is that Amber assumes that the list of angles
        and dihedrals encompasses *all* 1-3 and 1-4 pairs as determined by the
        bond graph, respectively. As a result, Amber programs use the angle and
        dihedral lists to set nonbonded exclusions and exceptions.

        Parameters
        ----------
        ignore_inconsistent_vdw : bool, optional
            If True, do not make inconsistent 1-4 vdW parameters fatal. For
            ChamberParm, the 1-4 specific vdW parameters can compensate. For
            AmberParm, the 1-4 scaling factor cannot represent arbitrary
            exceptions. Default is False (should only be True for ChamberParm)

        Returns
        -------
        n13, n14 : int, int
            The number of 1-3 and 1-4 pairs that needed to be added,
            respectively. Purely diagnostic
        """
        # We need to figure out what 1-4 scaling term to use if we don't have
        # explicit exceptions
        assert self.combining_rule in ('lorentz', 'geometric'), "Unrecognized combining rule"
        if not self.adjusts:
            scalings = defaultdict(int)
            for dih in self.dihedrals:
                if dih.ignore_end or dih.improper: continue
                scalings[(dih.type.scee, dih.type.scnb)] += 1
            if len(scalings) > 0:
                maxkey, maxval = next(iter(scalings.items()))
                for key, val in scalings.items():
                    if maxval < val:
                        maxkey, maxval = key, val
                scee, scnb = maxkey
            else:
                scee = scnb = 1e10
            zero_torsion = DihedralType(0, 1, 0, scee, scnb)
        else:
            # Turn list of exceptions into a dict so we can look it up quickly
            adjust_dict = dict()
            for pair in self.adjusts:
                adjust_dict[tuple(sorted([pair.atom1, pair.atom2]))] = pair
            ignored_torsion = None
            zero_torsion = None
            # Scan through existing dihedrals to make sure the exceptions match
            # the dihedral list
            if self.combining_rule == 'lorentz':
                comb_sig = lambda sig1, sig2: 0.5 * (sig1 + sig2)
            elif self.combining_rule == 'geometric':
                comb_sig = lambda sig1, sig2: sqrt(sig1 * sig2)
            fac = 2**(1/6)
            for dihedral in self.dihedrals:
                if dihedral.ignore_end: continue
                key = tuple(sorted([dihedral.atom1, dihedral.atom4]))
                eref = sqrt(dihedral.atom1.epsilon_14 * dihedral.atom4.epsilon_14)
                rref = comb_sig(dihedral.atom1.sigma_14, dihedral.atom4.sigma_14) * fac
                if key in adjust_dict:
                    pair = adjust_dict[key]
                    if pair.type.epsilon == 0:
                        scnb = 1e10
                    else:
                        scnb = eref / pair.type.epsilon
                    if pair.type.chgscale == 0:
                        scee = 1e10
                    else:
                        scee = 1 / pair.type.chgscale
                    if ignore_inconsistent_vdw:
                        scnb = 1.0
                    elif abs(rref - pair.type.rmin) > SMALL and pair.type.epsilon != 0:
                        raise TypeError('Cannot translate exceptions')
                    if (abs(scnb - dihedral.type.scnb) < SMALL and
                            abs(scee - dihedral.type.scee) < SMALL):
                        continue
                else:
                    scee = scnb = 1e10
                newtype = _copy.copy(dihedral.type)
                newtype.scee = scee
                newtype.scnb = scnb
                dihedral.type = newtype
                newtype.list = self.dihedral_types
                self.dihedral_types.append(newtype)

        zero_angle = AngleType(0, 0)

        n13 = n14 = 0
        if self.combining_rule == 'lorentz':
            comb_sig = lambda sig1, sig2: 0.5 * (sig1 + sig2)
        elif self.combining_rule == 'geometric':
            comb_sig = lambda sig1, sig2: sqrt(sig1 * sig2)
        fac = 2**(1/6)
        for atom in self.atoms:
            if isinstance(atom, ExtraPoint): continue
            for batom in atom.bond_partners:
                if isinstance(batom, ExtraPoint): continue
                for aatom in batom.bond_partners:
                    if isinstance(aatom, ExtraPoint) or aatom is atom: continue
                    for datom in aatom.bond_partners:
                        if isinstance(datom, ExtraPoint): continue
                        if (datom in atom.angle_partners + atom.bond_partners +
                                atom.dihedral_partners or datom is atom):
                            continue
                        # Add the missing dihedral
                        if not self.adjusts:
                            tortype = zero_torsion
                            if n14 == 0:
                                tortype.list = self.dihedral_types
                                self.dihedral_types.append(tortype)
                        else:
                            # Figure out what the scale factors must be
                            key = tuple(sorted([atom, datom]))
                            if key not in adjust_dict:
                                if ignored_torsion is None:
                                    ignored_torsion = DihedralType(0, 1, 0, 1e10, 1e10)
                                    self.dihedral_types.append(ignored_torsion)
                                    ignored_torsion.list = self.dihedral_types
                                tortype = ignored_torsion
                            elif 0 in (adjust_dict[key].type.epsilon, adjust_dict[key].type.rmin) \
                                 and adjust_dict[key].type.chgscale == 0:
                                if ignored_torsion is None:
                                    ignored_torsion = DihedralType(0, 1, 0, 1e10, 1e10,
                                                                   list=self.dihedral_types)
                                    self.dihedral_types.append(ignored_torsion)
                                tortype = ignored_torsion
                            else:
                                pair = adjust_dict[key]
                                epsilon = pair.type.epsilon
                                rmin = pair.type.rmin
                                # Compare it to the 1-4 parameters that are
                                # already present
                                eref = sqrt(pair.atom1.epsilon_14 * pair.atom2.epsilon_14)
                                if pair.type.epsilon == 0:
                                    scnb = 1e10
                                else:
                                    scnb = eref / epsilon
                                if pair.type.chgscale == 0:
                                    scee = 1e10
                                else:
                                    scee = 1 / pair.type.chgscale
                                rref = comb_sig(pair.atom1.sigma_14, pair.atom2.sigma_14) * fac
                                if abs(rmin - rref) > SMALL:
                                    if ignore_inconsistent_vdw:
                                        scnb = 1.0
                                    else:
                                        raise TypeError('Cannot translate exceptions')
                                tortype = DihedralType(0, 1, 0, scee, scnb, list=self.dihedral_types)
                                self.dihedral_types.append(tortype)
                        dihedral = Dihedral(atom, batom, aatom, datom, ignore_end=False,
                                            improper=False, type=tortype)
                        self.dihedrals.append(dihedral)
                        n14 += 1
                    if aatom in atom.angle_partners + atom.bond_partners:
                        continue
                    # Add the missing angle
                    self.angles.append(Angle(atom, batom, aatom, zero_angle))
                    n13 += 1

        if n13:
            self.angle_types.append(zero_angle)
            zero_angle.list = self.angle_types
        if n14: # See if there is some ambiguity here
            if not self.adjusts and len(scalings) > 1:
                warn('Multiple 1-4 scaling factors detected. Using the most-used values scee=%f '
                     'scnb=%f' % (scee, scnb), AmberWarning)
        return n13, n14

    #===================================================

    def _get_atom_collection_for_alternate_labels(self):
        atom_collection = [defaultdict(list) for r in self.residues]

        for adict, residue in zip(atom_collection, self.residues):
            for atom in residue.atoms:
                adict[atom.name].append(atom)
        return atom_collection

    def _label_alternates(self):
        atom_collection = self._get_atom_collection_for_alternate_labels()
        possible_labels = list('ABCDEFGHIJKLMNOPQRSTUVWXYZ')

        for _, adict in enumerate(atom_collection):
            for atom_name, atom_list in adict.items():
                if len(atom_list) > 1:
                    for i, atom in enumerate(atom_list):
                        label = possible_labels[i%len(possible_labels)]
                        atom.altloc = label

    #===================================================

    @property
    def box(self):
        if self._box is not None:
            return self._box[0]
        return None

    @box.setter
    def box(self, value):
        # Deleting the box is more complicated for AmberParm, so override it
        # here
        if value is None:
            # Delete all of the other box info in the prmtop
            self._box = None
            for flag in ('IPTRES', 'NSPM', 'NSPSOL'):
                if flag in self.pointers:
                    del self.pointers[flag]
            for flag in ('SOLVENT_POINTERS', 'ATOMS_PER_MOLECULE', 'BOX_DIMENSIONS'):
                self.delete_flag(flag)
            self.hasbox = False
        else:
            if isinstance(value, np.ndarray):
                box = value
            else:
                box = list(value)
                if len(box) != 6:
                    raise ValueError('Box information must be 6 floats')
                if u.is_quantity(box[0]):
                    box[0] = box[0].value_in_unit(u.angstroms)
                if u.is_quantity(box[1]):
                    box[1] = box[1].value_in_unit(u.angstroms)
                if u.is_quantity(box[2]):
                    box[2] = box[2].value_in_unit(u.angstroms)
                if u.is_quantity(box[3]):
                    box[3] = box[3].value_in_unit(u.degrees)
                if u.is_quantity(box[4]):
                    box[4] = box[4].value_in_unit(u.degrees)
                if u.is_quantity(box[5]):
                    box[5] = box[5].value_in_unit(u.degrees)
            box = np.asanyarray(box, dtype=np.float64).reshape((-1, 6))

            # We are adding a box for the first time, so make sure we add some flags
            if self._box is None:
                self._box = box
                # We need to add topology information
                if 'SOLVENT_POINTERS' not in self.flag_list:
                    self.add_flag('SOLVENT_POINTERS', '3I8', num_items=3, after='IROTAT')
                if 'ATOMS_PER_MOLECULE' not in self.flag_list:
                    self.add_flag('ATOMS_PER_MOLECULE', '10I8', data=[0], after='SOLVENT_POINTERS')
                if 'BOX_DIMENSIONS' not in self.flag_list:
                    self.add_flag('BOX_DIMENSIONS', '5E16.8', after='ATOMS_PER_MOLECULE',
                                  data=[box[0,3], box[0,0], box[0,1], box[0,2]])
                try:
                    self.rediscover_molecules(fix_broken=False)
                except MoleculeError:
                    # Do not reorder molecules here -- only do that when
                    # specifically requested. Otherwise we could get out-of-sync
                    # with coordinates.
                    pass
                self.load_pointers()
            else:
                self.parm_data['BOX_DIMENSIONS'] = [box[0,3], box[0,0], box[0,1], box[0,2]]
                self._box = box
        self._set_ifbox()

    def _set_ifbox(self):
        """ Sets the IFBOX pointers to 1 (ortho), 2 (octahedral) , or 3 (other) """
        if self.box is None:
            self.parm_data['POINTERS'][PrmtopPointers.IFBOX] = 0
            self.pointers['IFBOX'] = 0
        elif np.allclose(self.box[3:], 90):
            self.parm_data['POINTERS'][PrmtopPointers.IFBOX] = 1
            self.pointers['IFBOX'] = 1
        elif np.allclose(self.box[3:], TRUNCATED_OCTAHEDRON_ANGLE, atol=0.02):
            self.parm_data['POINTERS'][PrmtopPointers.IFBOX] = 2
            self.pointers['IFBOX'] = 2
        else:
            # General triclinic
            self.parm_data['POINTERS'][PrmtopPointers.IFBOX] = 3
            self.pointers['IFBOX'] = 3

    def _cleanup_dihedrals_with_periodicity_zero(self):
        """
        For torsions with only a single term and a periodicity set to 0, make sure pmemd still
        properly recognizes the necessary exception parameters. update_dihedral_exclusions will
        make sure that if a dihedral has a type pn0 *and* ignore_end is set to False (which means
        that it is required to specify exclusions), then it is the *only* torsion between those
        atoms in the system. This allows us to scan through our dihedrals, look for significant
        terms that have pn==0, and simply add another dihedral with pn=1 and k=0 to ensure that
        pmemd will always get that exception correct
        """
        new_dihedrals = []
        for dih in self.dihedrals:
            if dih.ignore_end or dih.type.per != 0:
                continue
            # If we got here, ignore_end must be False and out periodicity must be 0. So add
            # another dihedral
            dt = DihedralType(0, 1, 0, dih.type.scee, dih.type.scnb, list=self.dihedral_types)
            self.dihedral_types.append(dt)
            new_dihedrals.append(
                Dihedral(dih.atom1, dih.atom2, dih.atom3, dih.atom4, improper=dih.improper,
                         ignore_end=False, type=dt)
            )
            # Now that we added the above dihedral, we can start ignoring the end-group interactions
            # on this dihedral
            dih.ignore_end = True
        if new_dihedrals:
            self.dihedrals.extend(new_dihedrals)

    #===================================================

    _AMBERPARM_ATTRS = 'LJ_types LJ_radius LJ_depth parm_data pointers'.split()

    def __getstate__(self):
        d = Structure.__getstate__(self)
        d.update(AmberFormat.__getstate__(self))
        for attr in self._AMBERPARM_ATTRS:
            if getattr(self, attr, None) is not None:
                d[attr] = getattr(self, attr)
        return d

    def __setstate__(self, d):
        AmberFormat.__setstate__(self, d)
        Structure.__setstate__(self, d)
        for attr in self._AMBERPARM_ATTRS:
            if attr in d:
                setattr(self, attr, d[attr])


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class Rst7(object):
    """
    Amber input coordinate (or restart coordinate) file. Front-end for the
    readers and writers, supports both NetCDF and ASCII restarts.

    Parameters
    ----------
    filename : str, optional
        If a filename is provided, this file is parsed and the Rst7 data
        populated from that file. The format (ASCII or NetCDF) is autodetected
    natom : int, optional
        If no filename is provided, this value is required. If a filename is
        provided, this value is ignored (and instead set to the value of natom
        from the coordinate file). This is the number of atoms for which we have
        coordinates. If not provided for a new file, it *must* be set later.
    title : str, optional
        For a file that is to be written, this is the title that will be given
        to that file. Default is an empty string
    time : float, optional
        The time to write to the restart file. This is cosmetic. Default is 0
    """

    def __init__(self, filename=None, natom=None, title='', time=0.0):
        """
        Optionally takes a filename to read. This is deprecated, though, as the
        alternative constructor "open" should be used instead
        """
        self.coordinates = []
        self.vels = None
        self._box = None
        self.natom = natom
        self.title = title
        self.time = 0
        if filename is not None:
            self.filename = filename
            self._read(filename)

    @property
    def box(self):
        return self._box
    @box.setter
    def box(self, value):
        if value is None:
            self._box = None
        else:
            self._box = np.array(value).reshape((-1, 6))[0]

    @classmethod
    def open(cls, filename):
        """ Constructor that opens and parses an input coordinate file

        Parameters
        ----------
        filename : str
            Name of the file to parse
        """
        inst = cls()
        inst.filename = filename
        inst._read(filename)
        return inst

    def _read(self, filename):
        """
        Open and parse an input coordinate file in either ASCII or NetCDF format
        """
        try:
            f = AmberAsciiRestart(filename, 'r')
            self.natom = f.natom
        except ValueError:
            # Maybe it's a NetCDF file?
            try:
                f = NetCDFRestart.open_old(filename)
                self.natom = f.atom
            except (TypeError, RuntimeError):
                raise AmberError('Could not parse restart file %s' % filename)

        self.coordinates = f.coordinates
        if f.hasvels:
            self.vels = f.velocities
        if f.hasbox:
            self.box = f.box
        self.title = f.title
        self.time = f.time

    @classmethod
    def copy_from(cls, thing):
        """
        Copies the coordinates, velocities, and box information from another
        instance
        """
        inst = cls()
        inst.natom = thing.natom
        inst.title = thing.title
        inst.coordinates = thing.coordinates[:]
        if hasattr(thing, 'vels'): inst.vels = _copy.deepcopy(thing.vels)
        if hasattr(thing, 'box'): inst.box = _copy.deepcopy(thing.box)
        inst.time = thing.time

        return inst

    def __copy__(self):
        """ Copy constructor """
        return type(self).copy_from(self)

    def write(self, fname, netcdf=False):
        """ Writes the coordinates and/or velocities to a restart file """
        if netcdf:
            if self.natom is None:
                raise RuntimeError('Number of atoms must be set for NetCDF '
                                   'Restart files before write time')
            f = NetCDFRestart.open_new(fname, self.natom, self.box is not None,
                                       self.vels is not None, self.title)
        else:
            f = AmberAsciiRestart(fname, 'w', natom=self.natom,
                                  title=self.title)

        f.time = self.time
        # Now write the coordinates
        f.coordinates = self.coordinates
        if self.vels is not None:
            f.velocities = self.vels
        if self.box is not None:
            f.box = self.box
        f.close()

    @property
    def positions(self):
        """ Atomic coordinates with units """
        coordinates = self.coordinates.reshape(self.natom, 3)
        return [Vec3(*x) for x in coordinates] * u.angstroms

    @property
    def velocities(self):
        """ Atomic velocities in units of angstroms/picoseconds """
        return np.asarray(self.vels).reshape(self.natom, 3)

    @property
    def box_vectors(self):
        """ Unit cell vectors with units """
        if self.box is None: return None
        return box_lengths_and_angles_to_vectors(*self.box)

    @property
    def hasbox(self):
        """ Whether or not this Rst7 has unit cell information """
        return self.box is not None

    @property
    def hasvels(self):
        """ Whether or not this Rst7 has velocities """
        return self.vels is not None

def _zeros(length):
    """ Returns an array of zeros of the given length """
    return [0 for i in range(length)]
