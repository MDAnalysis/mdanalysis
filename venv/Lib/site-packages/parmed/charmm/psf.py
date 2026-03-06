"""
Provides a Python class for parsing a PSF file and setting up a system
structure for it
"""
import math
import re
from copy import copy as _copy
from functools import wraps
from ..constants import DEG_TO_RAD
from ..periodic_table import AtomicNum, element_by_mass
from ..topologyobjects import (Bond, Angle, Dihedral, Improper, AcceptorDonor, Group, Cmap,
                               UreyBradley, NoUreyBradley, Atom, DihedralType, ImproperType,
                               UnassignedAtomType, ExtraPoint, DrudeAtom, LocalCoordinatesFrame,
                               DrudeAnisotropy)
from ..exceptions import CharmmError, ParameterError
from ..structure import needs_openmm, Structure
from ..utils.io import genopen

def _catchindexerror(func):
    """
    Protects a function from raising an index error, and replace that exception
    with a CharmmError instead
    """
    @wraps(func)
    def newfunc(*args, **kwargs):
        """ Catch the index error """
        try:
            return func(*args, **kwargs)
        except IndexError as e:
            raise CharmmError(f'Array is too short: {e}') from e

    return newfunc

class _FileEOF(Exception):
    """ For control flow """

class _ZeroDict(dict):
    """
    Contains a dict that returns dummy (zero) arguments when a key is not
    present rather than raising a KeyError.  The return value for non-existent
    items is (0, []). It also special-case sections that have multiple pointers
    to avoid index errors if those are not present in the PSF file
    """
    def __getitem__(self, key):
        try:
            return dict.__getitem__(self, key)
        except KeyError:
            if key.startswith('NGRP'):
                for k in self:
                    if k.startswith('NGRP'):
                        return dict.__getitem__(self, k)
                return [0, 0], []
            elif key.startswith('NUMLP'):
                for k in self:
                    if k.startswith('NUMLP'):
                        return dict.__getitem__(self, k)
                return [0, 0], []
            return 0, []

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

_resre = re.compile(r'(-?\d+)([a-zA-Z]*)')

class CharmmPsfFile(Structure):
    """
    A chemical :class:`Structure` instantiated from CHARMM files.

    Parameters
    ----------
    psf_name : str, optional
        Name of the PSF file (it must exist)

    Raises
    ------
    IOError : If file ``psf_name`` does not exist
    CharmmPsfError : If any parsing errors are encountered
    """
    @staticmethod
    def _convert(string, type, message):
        """
        Converts a string to a specific type, making sure to raise
        CharmmError with the given message in the event of a failure.

        Parameters
        ----------
        string : str
            Input string to process
        type : type
            Type of data to convert to (e.g., ``int``)
        message : str
            Error message to put in exception if failed
        """
        try:
            return type(string)
        except ValueError:
            raise CharmmError(f'Could not convert {message} [{string}]')

    #===================================================

    @classmethod
    def _parse_psf_title_line(cls, line):
        words = line[:line.index('!')].split()
        title = line[line.index('!')+1:].strip().upper()
        # Strip out description
        if ':' in title:
            title = title[:title.index(':')]
        if len(words) == 1:
            pointers = cls._convert(words[0], int, 'pointer')
        else:
            pointers = tuple([cls._convert(w, int, 'pointer') for w in words])
        return title, pointers

    #===================================================

    @classmethod
    def _parse_psf_section(cls, psf):
        """
        This method parses a section of the PSF file

        Parameters
        ----------
        psf : file
            Open file that is pointing to the first line of the section that is
            to be parsed

        Returns
        -------
        title : str
            The label of the PSF section we are parsing
        pointers : (int/tuple of ints)
            If one pointer is set, pointers is simply the integer that is value
            of that pointer. Otherwise it is a tuple with every pointer value
            defined in the first line
        data : list
            A list of all data in the parsed section converted to integers
        """
        line = psf.readline()
        while not line.strip():
            if not line:
                raise _FileEOF('Unexpected EOF in PSF file')
            else:
                line = psf.readline()
        if '!' in line:
            title, pointers = cls._parse_psf_title_line(line)
        else:
            raise CharmmError('Could not determine section title') # pragma: no cover
        line = psf.readline().strip()
        if not line and title.startswith('NNB'):
            # This will correctly handle the NNB section (which has a spurious
            # blank line) as well as any sections that have 0 members.
            line = psf.readline().strip()
            # If this line has a title in it, then it's one of those weird PSF files that
            # has basically no NNB section. Just skip over NNB and take the next section
            # instead
            if '!' in line:
                title, pointers = cls._parse_psf_title_line(line)
                line = psf.readline().strip()
        data = []
        if title in {'NATOM', 'NTITLE', 'NUMLP NUMLPH', 'NUMANISO'}:
            # Store these two sections as strings to be parsed later.
            # The rest of the sections are integer pointers
            while line:
                data.append(line)
                line = psf.readline().strip()
        else:
            while line:
                words = line.split()
                data.extend([cls._convert(w, int, 'PSF data') for w in words])
                line = psf.readline().strip()
        return title, pointers, data

    #===================================================

    @_catchindexerror
    def __init__(self, psf_name=None):
        """
        Opens and parses a PSF file, then instantiates a CharmmPsfFile
        instance from the data.
        """
        from ..utils import tag_molecules
        super().__init__()
        # Bail out if we don't have a filename
        if psf_name is None:
            return
        # Open the PSF and read the first line. It must start with "PSF"
        if isinstance(psf_name, str):
            fileobj = genopen(psf_name, 'r')
            own_handle = True
        else:
            fileobj = psf_name
            own_handle = False
        try:
            self.name = psf_name if isinstance(psf_name, str) else ''
            line = fileobj.readline()
            if not line.startswith('PSF'):
                raise CharmmError(f'Unrecognized PSF file. First line is {line.strip()}')
            # Store the flags
            psf_flags = line.split()[1:]
            # Now get all of the sections and store them in a dict
            fileobj.readline()
            # Now get all of the sections
            psfsections = _ZeroDict()
            while True:
                try:
                    sec, ptr, data = self._parse_psf_section(fileobj)
                except _FileEOF:
                    break
                psfsections[sec] = (ptr, data)
            # store the title
            self.title = psfsections['NTITLE'][1]
            # Next is the number of atoms
            natom = self._convert(psfsections['NATOM'][0], int, 'natom')
            # Parse all of the atoms
            drude_alpha_thole = []
            is_drude = 'DRUDE' in psf_flags
            for i in range(natom):
                words = psfsections['NATOM'][1][i].split()
                atid = int(words[0])
                if atid != i + 1:
                    raise CharmmError('Nonsequential atoms detected!')
                segid = words[1]
                rematch = _resre.match(words[2])
                if not rematch:
                    raise CharmmError(f'Could not interpret residue number {words[2]}')
                resid, inscode = rematch.groups()
                resid = self._convert(resid, int, 'residue number')
                resname = words[3]
                name = words[4]
                attype = words[5]
                # Try to convert the atom type to an integer a la CHARMM
                try:
                    attype = int(attype)
                except ValueError:
                    pass
                charge = self._convert(words[6], float, 'partial charge')
                mass = self._convert(words[7], float, 'atomic mass')
                props = words[8:]
                if is_drude:
                    alpha = self._convert(words[9], float, 'alpha')
                    thole = self._convert(words[10], float, 'thole')
                    drude_alpha_thole.append((alpha, thole))
                if is_drude and i >= 1 and drude_alpha_thole[-2] != (0, 0):
                    # This assumes that the Drude atom is defined immediately after its parent atom.
                    # This always appears to be the case, but if it proves untrue at some point then
                    # this will need to be refactored to identify Drude atoms after identifying bonds.
                    my_alpha, my_thole = drude_alpha_thole[-2]
                    atom = DrudeAtom(name=name, type=attype, charge=charge, mass=mass, parent=self.atoms[-1],
                                     atomic_number=0, alpha=my_alpha, thole=my_thole, drude_type=attype)
                elif (name.startswith('LP') or (isinstance(attype, str) and attype.startswith("LP"))) and mass == 0:
                    atom = ExtraPoint(name=name, type=attype, charge=charge, mass=mass, atomic_number=0)
                else:
                    atom = Atom(name=name, type=attype, charge=charge, mass=mass,
                                atomic_number=AtomicNum[element_by_mass(mass)])
                atom.props = props
                self.add_atom(atom, resname, resid, chain=segid, inscode=inscode, segid=segid)
            # Now get the number of bonds
            nbond = self._convert(psfsections['NBOND'][0], int, 'number of bonds')
            if len(psfsections['NBOND'][1]) != nbond * 2:
                raise CharmmError(f"Got {len(psfsections['NBOND'][1])} indexes for {nbond} bonds")
            it = iter(psfsections['NBOND'][1])
            for i, j in zip(it, it):
                self.bonds.append(Bond(self.atoms[i-1], self.atoms[j-1]))
            # Now get the number of angles and the angle list
            ntheta = self._convert(psfsections['NTHETA'][0], int, 'number of angles')
            if len(psfsections['NTHETA'][1]) != ntheta * 3:
                raise CharmmError(f"Got {len(psfsections['NTHETA'][1])} indexes for {ntheta} angles")
            it = iter(psfsections['NTHETA'][1])
            for i, j, k in zip(it, it, it):
                self.angles.append(Angle(self.atoms[i-1], self.atoms[j-1], self.atoms[k-1]))
                self.angles[-1].funct = 5 # urey-bradley
            # Now get the number of torsions and the torsion list
            nphi = self._convert(psfsections['NPHI'][0], int, 'number of torsions')
            if len(psfsections['NPHI'][1]) != nphi * 4:
                raise CharmmError(f"Got {len(psfsections['NPHI'])} indexes for {nphi} torsions")
            it = iter(psfsections['NPHI'][1])
            for i, j, k, l in zip(it, it, it, it):
                self.dihedrals.append(
                    Dihedral(self.atoms[i-1], self.atoms[j-1], self.atoms[k-1], self.atoms[l-1])
                )
            self.dihedrals.split = False
            # Now get the number of improper torsions
            nimphi = self._convert(psfsections['NIMPHI'][0], int, 'number of impropers')
            if len(psfsections['NIMPHI'][1]) != nimphi * 4:
                raise CharmmError(f"Got {len(psfsections['NIMPHI'][1])} indexes for {nimphi} impropers")
            it = iter(psfsections['NIMPHI'][1])
            for i, j, k, l in zip(it, it, it, it):
                self.impropers.append(
                    Improper(self.atoms[i-1], self.atoms[j-1], self.atoms[k-1], self.atoms[l-1])
                )
            # Now handle the donors (what is this used for??)
            ndon = self._convert(psfsections['NDON'][0], int, 'number of donors')
            if len(psfsections['NDON'][1]) != ndon * 2:
                raise CharmmError(f"Got {len(psfsections['NDON'][1])} indexes for {ndon} donors")
            it = iter(psfsections['NDON'][1])
            for i, j in zip(it, it):
                self.donors.append(AcceptorDonor(self.atoms[i-1], self.atoms[j-1]))
            # Now handle the acceptors (what is this used for??)
            nacc = self._convert(psfsections['NACC'][0], int, 'number of acceptors')
            if len(psfsections['NACC'][1]) != nacc * 2:
                raise CharmmError(f"Got {len(psfsections['NACC'][1])} indexes for {nacc} acceptors")
            it = iter(psfsections['NACC'][1])
            for i, j in zip(it, it):
                self.acceptors.append(AcceptorDonor(self.atoms[i-1], self.atoms[j-1]))
            # Now get the group sections
            try:
                ngrp, nst2 = psfsections['NGRP NST2'][0]
            except ValueError: # pragma: no cover
                raise CharmmError('Could not unpack GROUP pointers') # pragma: no cover
            tmp = psfsections['NGRP NST2'][1]
            self.groups.nst2 = nst2
            # Now handle the groups
            if len(psfsections['NGRP NST2'][1]) != ngrp * 3:
                raise CharmmError(f"Got {len(tmp)} indexes for {ngrp} groups")
            it = iter(psfsections['NGRP NST2'][1])
            for i, j, k in zip(it, it, it):
                self.groups.append(Group(self.atoms[i], j, k))
            # Assign all of the atoms to molecules recursively
            tmp = psfsections['MOLNT'][1]
            tag_molecules(self)
            if len(psfsections["NUMLP NUMLPH"][1]) != 0:
                # We have a CHARMM PSF file; now do NUMLP/NUMLPH sections
                self._process_lonepair_section(psfsections["NUMLP NUMLPH"])
            # Now process the NUMANISO records if this is a drude PSF
            if is_drude:
                self._process_aniso_section(psfsections["NUMANISO"])
            # Now do the CMAPs
            ncrterm = self._convert(psfsections['NCRTERM'][0], int, 'Number of cross-terms')
            if len(psfsections['NCRTERM'][1]) != ncrterm * 8:
                raise CharmmError(f"Got {len(psfsections['NCRTERM'])} CMAP indexes for {ncrterm} cmap terms")
            it = iter(psfsections['NCRTERM'][1])
            for i, j, k, l, m, n, o, p in zip(it, it, it, it, it, it, it, it):
                self.cmaps.append(
                    Cmap.extended(
                        self.atoms[i - 1],
                        self.atoms[j - 1],
                        self.atoms[k - 1],
                        self.atoms[l - 1],
                        self.atoms[m - 1],
                        self.atoms[n - 1],
                        self.atoms[o - 1],
                        self.atoms[p - 1],
                    )
                )
            self.unchange()
            self.flags = psf_flags
        finally:
            if own_handle:
                fileobj.close()

    #===================================================

    def _process_lonepair_section(self, section):
        lp_distance_list = []
        lp_angle_list = []
        lp_dihedral_list = []
        lp_hostnum_list = []
        numlp, numlph = section[0]
        lp_lines = section[1]
        if len(lp_lines) != numlp + (numlph + 7) // 8:
            raise CharmmError("Unexpected number of lines in NUMLP/NUMLPH section")
        for i in range(numlp):
            lp_line = lp_lines[i].split()
            if len(lp_line) != 6 or lp_line[2] != 'F':
                raise CharmmError("lonepair format error")
            lp_hostnum_list.append(int(lp_line[0]))
            lp_distance_list.append(float(lp_line[3]))
            lp_angle_list.append(float(lp_line[4]))
            lp_dihedral_list.append(float(lp_line[5]))
        lp_cnt = 0
        for i in range(numlp):
            idall = []
            for j in range(lp_hostnum_list[i] + 1):
                words = lp_lines[int((lp_cnt + j) // 8) + numlp].split()
                icol = (lp_cnt + j) % 8
                idall.append(int(words[icol]) - 1)
            if len(idall) not in (3, 4):
                raise CharmmError(f"Unsupported lone pair configuration (must have 2 or 3 atoms in frame, not {len(idall) - 1})")
            if len(idall) == 4:
                frame = self._get_frame3(idall[1], idall[3], idall[2], lp_distance_list[i], lp_angle_list[i], lp_dihedral_list[i])
            else:
                frame = self._get_frame2(idall[1], idall[2], lp_distance_list[i])
            if not isinstance(self.atoms[idall[0]], ExtraPoint):
                raise CharmmError(f"Expected atom {idall[0]} to be an ExtraPoint but it is not")
            self.atoms[idall[0]].frame_type = frame
            lp_cnt += lp_hostnum_list[i] + 1

    #===================================================

    def _process_aniso_section(self, section):
        num_aniso, data = section
        if num_aniso == 0:
            return
        k_list = []
        pointers = []

        for i in range(num_aniso):
            k_list.append([float(x) for x in data[i].split()])
        for i in range((num_aniso + 1) // 2):
            pointers.extend([int(x) for x in data[i + num_aniso].split()])
        if len(pointers) != 4 * len(k_list):
            raise CharmmError(
                f"Bad NUMANISO section - mismatching numbers of atom indices and k ({len(pointers)} vs {len(k_list)})"
            )
        it = iter(pointers)
        parent_to_drude = {atom.parent: atom for atom in self.atoms if isinstance(atom, DrudeAtom)}
        for id1, id2, id3, id4, (k11, k22, k33) in zip(it, it, it, it, k_list):
            parent_atom: Atom = self.atoms[id1 - 1]
            drude_atom = parent_to_drude[parent_atom]
            # Calculate a11 and a12, per the code originally found in OpenMM
            a = k11 + k22 + (3 * k33)
            b = (2 * k11 * k22) + (4 * k11 * k33) + (6 * k33 * k33)
            c = 3 * k33 * (k11 + k33) * (k22 + k33)
            drude_k = (math.sqrt(b * b - 4 * a * c)) / (2 * a)
            a11 = round(drude_k / (k11 + k33 + drude_k), 5)
            a22 = round(drude_k / (k22 + k33 + drude_k), 5)
            at2, at3, at4 = self.atoms[id2 - 1], self.atoms[id3 - 1], self.atoms[id4 - 1]
            drude_atom.anisotropy = DrudeAnisotropy(
                parent_atom, at2, at3, at4, a11, a22, k11=k11, k22=k22, k33=k33
            )

    #===================================================

    def _get_frame3(
        self,
        a1idx: int,
        a2idx: int,
        a3idx: int,
        dist: float,
        ang: float,
        dihed: float,
    ) -> LocalCoordinatesFrame:

        if dist > 0:
            x_weights = [-1.0, 0.0, 1.0]
        elif dist < 0:
            x_weights = [-1.0, 0.5, 0.5]
        origin_weights = [1, 0, 0]
        y_weights = [0, -1, 1]
        return LocalCoordinatesFrame(
            self.atoms[a1idx], self.atoms[a2idx], self.atoms[a3idx],
            origin_weights, x_weights, y_weights, dist, ang, dihed, 3
        )

    def _get_frame2(
        self,
        a1idx: int,
        a2idx: int,
        dist: float,
    ) -> LocalCoordinatesFrame:

        # TODO - figure out if this can be done with another virtual site type. This comes from the OpenMM implementation
        # in the CHARMM parsers there
        a1 = self.atoms[a1idx]
        a2 = self.atoms[a2idx]
        for a3 in a1.bond_partners + a2.bond_partners:
            if a3 not in (a1, a2):
                break
        else:
            raise CharmmError("Could not find a third atom to define in the LocalCoordinatesFrame")
        origin_weights = [1, 0, 0]
        x_weights = [1, -1, 0]
        y_weights = [0, -1, 1]
        return LocalCoordinatesFrame(a1, a2, a3, origin_weights, x_weights, y_weights, dist, 0, 0, 2)

    #===================================================

    @classmethod
    def from_structure(cls, struct, copy=False):
        """
        Instantiates a CharmmPsfFile from an input Structure instance. This
        method makes sure all atom types have uppercase-only names

        Parameters
        ----------
        struct : :class:`parmed.structure.Structure`
            The input structure to convert to a CharmmPsfFile instance
        copy : bool, optional
            If True, a copy of all items are made. Otherwise, the resulting
            CharmmPsfFile is a shallow copy

        Returns
        -------
        psf : :class:`CharmmPsfFile`
            CHARMM PSF file

        Raises
        ------
        ValueError if the functional form is not recognized or cannot be
        implemented through the PSF and parameter/stream files

        Notes
        -----
        If copy is False, the original object may have its atom type names
        changed if any of them have lower-case letters
        """
        from .parameters import _typeconv as typeconv
        if (struct.rb_torsions or struct.trigonal_angles or
                struct.out_of_plane_bends or struct.pi_torsions or
                struct.stretch_bends or struct.torsion_torsions or
                struct.chiral_frames or struct.multipole_frames or
                struct.nrexcl != 3):
            raise ValueError('Unsupported functional form for CHARMM PSF')
        if copy:
            struct = _copy(struct)
        psf = cls()
        psf.atoms = struct.atoms
        psf.residues = struct.residues
        psf.bonds = struct.bonds
        psf.angles = struct.angles
        psf.urey_bradleys = struct.urey_bradleys
        psf.dihedrals = struct.dihedrals
        psf.impropers = struct.impropers
        psf.acceptors = struct.acceptors
        psf.donors = struct.donors
        psf.groups = struct.groups
        psf.cmaps = struct.cmaps

        psf.bond_types = struct.bond_types
        psf.angle_types = struct.angle_types
        psf.urey_bradley_types = struct.urey_bradley_types
        psf.dihedral_types = struct.dihedral_types
        psf.improper_types = struct.improper_types
        psf.cmap_types = struct.cmap_types

        for atom in psf.atoms:
            atom.type = typeconv(atom.type)
            if atom.atom_type is not UnassignedAtomType:
                atom.atom_type.name = typeconv(atom.atom_type.name)

        # If no groups are defined, make each residue its own group
        if not psf.groups:
            for residue in psf.residues:
                chg = sum(a.charge for a in residue)
                if chg < 1e-4:
                    psf.groups.append(Group(residue[0], 1, 0))
                else:
                    psf.groups.append(Group(residue[0], 2, 0))
            psf.groups.nst2 = 0

        return psf

    #===================================================

    def __str__(self):
        return self.name

    #===================================================

    @needs_openmm
    def createSystem(self, params=None, *args, **kwargs):
        """
        Creates an OpenMM System object from the CHARMM PSF file. This is a
        shortcut for calling `load_parameters` followed by
        Structure.createSystem. If params is not None, `load_parameters` will be
        called on that parameter set, and Structure.createSystem will be called
        with the remaining args and kwargs

        Parameters
        ----------
        params : CharmmParameterSet=None
            If not None, this parameter set will be loaded

        See Also
        --------
        :meth:`parmed.structure.Structure.createSystem`
            In addition to `params`, this method also takes all arguments for
            :meth:`parmed.structure.Structure.createSystem`
        """
        if params is not None:
            self.load_parameters(params)
        return super().createSystem(*args, **kwargs)

    #===================================================

    def load_parameters(self, parmset, copy_parameters=True):
        """
        Loads parameters from a parameter set that was loaded via CHARMM RTF,
        PAR, and STR files.

        Parameters
        ----------
        parmset : :class:`CharmmParameterSet`
            List of all parameters

        copy_parameters : bool, optional, default=True
            If False, parmset will not be copied.

            WARNING:
            -------
            Not copying parmset will cause ParameterSet and Structure to share
            references to types.  If you modify the original parameter set, the
            references in Structure list_types will be silently modified.
            However, if you change any reference in the parameter set, then that
            reference will no longer be shared with structure.

            Example where the reference in ParameterSet is changed. The
            following will NOT modify the parameters in the psf::

                psf.load_parameters(parmset, copy_parameters=False)
                parmset.angle_types[('a1', 'a2', a3')] = AngleType(1, 2)

            The following WILL change the parameter in the psf because the
            reference has not been changed in ``ParameterSet``::

                psf.load_parameters(parmset, copy_parameters=False)
                a = parmset.angle_types[('a1', 'a2', 'a3')]
                a.k = 10
                a.theteq = 100

            Extra care should be taken when trying this with dihedral_types.
            Since dihedral_type is a Fourier sequence, ParameterSet stores
            DihedralType for every term in DihedralTypeList. Therefore, the
            example below will STILL modify the type in the :class:`Structure`
            list_types::

                parmset.dihedral_types[('a', 'b', 'c', 'd')][0] = DihedralType(1, 2, 3)

            This assigns a new instance of DihedralType to an existing
            DihedralTypeList that ParameterSet and Structure are tracking and
            the shared reference is NOT changed.

            Use with caution!

        Notes
        -----
        - If any dihedral or improper parameters cannot be found, I will try
          inserting wildcards (at either end for dihedrals and as the two
          central atoms in impropers) and see if that matches.  Wild-cards will
          apply ONLY if specific parameters cannot be found.

        - This method will expand the dihedrals attribute by adding a separate
          Dihedral object for each term for types that have a multi-term
          expansion

        Raises
        ------
        ParameterError if any parameters cannot be found
        """
        if copy_parameters:
            parmset = _copy(parmset)
        self.combining_rule = parmset.combining_rule
        # First load the atom types
        for atom in self.atoms:
            try:
                if isinstance(atom.type, int):
                    atype = parmset.atom_types_int[atom.type]
                else:
                    atype = parmset.atom_types_str[atom.type.upper()]
            except KeyError:
                raise ParameterError(f"Could not find atom type for {atom.type}")
            atom.atom_type = atype
            # Change to string type to look up the rest of the parameters. Use
            # upper-case since all parameter sets were read in as upper-case
            atom.type = str(atom.atom_type).upper()
            atom.atomic_number = atype.atomic_number

        # Next load all of the bonds
        skipped_bonds = set()
        for bond in self.bonds:
            # Skip any bonds with drude atoms or virtual sites. They are not stored.
            # Depending on how Drude support is implemented in Amber (if that ever happens), we
            # may have to add dummy values here.
            if isinstance(bond.atom1, (DrudeAtom, ExtraPoint)) or isinstance(bond.atom2, (DrudeAtom, ExtraPoint)):
                skipped_bonds.add(bond)
                continue
            # Construct the key
            key = (min(bond.atom1.type, bond.atom2.type), max(bond.atom1.type, bond.atom2.type))
            try:
                bond.type = parmset.bond_types[key]
            except KeyError:
                raise ParameterError(f"Missing bond type for {bond}")
            bond.type.used = False
        # Build the bond_types list
        del self.bond_types[:]
        for bond in self.bonds:
            if bond in skipped_bonds:
                continue
            if bond.type.used:
                continue
            bond.type.used = True
            self.bond_types.append(bond.type)
            bond.type.list = self.bond_types
        # Next load all of the angles. If a Urey-Bradley term is defined for
        # this angle, also build the urey_bradley and urey_bradley_type lists
        del self.urey_bradleys[:]
        for ang in self.angles:
            # Construct the key
            key = (min(ang.atom1.type, ang.atom3.type),
                   ang.atom2.type,
                   max(ang.atom1.type, ang.atom3.type))
            try:
                ang.type = parmset.angle_types[key]
                ang.type.used = False
                ubt = parmset.urey_bradley_types[key]
                if ubt is not NoUreyBradley:
                    ub = UreyBradley(ang.atom1, ang.atom3, ubt)
                    self.urey_bradleys.append(ub)
                    ubt.used = False
            except KeyError:
                raise ParameterError(f"Missing angle type for {ang}")
        del self.urey_bradley_types[:]
        del self.angle_types[:]
        for ub in self.urey_bradleys:
            if ub.type.used:
                continue
            ub.type.used = True
            self.urey_bradley_types.append(ub.type)
            ub.type.list = self.urey_bradley_types
        for ang in self.angles:
            if ang.type.used:
                continue
            ang.type.used = True
            self.angle_types.append(ang.type)
            ang.type.list = self.angle_types
        # Next load all of the dihedrals.
        active_dih_list = set()
        for dih in self.dihedrals:
            # Store the atoms
            a1, a2, a3, a4 = dih.atom1, dih.atom2, dih.atom3, dih.atom4
            key = (a1.type, a2.type, a3.type, a4.type)
            # First see if the exact dihedral is specified
            if not key in parmset.dihedral_types:
                # Check for wild-cards
                key = ('X', a2.type, a3.type, 'X')
                if not key in parmset.dihedral_types:
                    raise ParameterError(f'No dihedral parameters found for {dih}')
            dih.type = parmset.dihedral_types[key]
            dih.type.used = False
            pair = (dih.atom1.idx, dih.atom4.idx) # To determine exclusions
            if (dih.atom1 in dih.atom4.bond_partners or
                dih.atom1 in dih.atom4.angle_partners):
                dih.ignore_end = True
            elif pair in active_dih_list:
                dih.ignore_end = True
            else:
                active_dih_list.add(pair)
                active_dih_list.add((dih.atom4.idx, dih.atom1.idx))
        del self.dihedral_types[:]
        for dihedral in self.dihedrals:
            if dihedral.type.used:
                continue
            dihedral.type.used = True
            self.dihedral_types.append(dihedral.type)
            dihedral.type.list = self.dihedral_types
        # Now do the impropers
        for imp in self.impropers:
            a1, a2, a3, a4 = imp.atom1.type, imp.atom2.type, imp.atom3.type, imp.atom4.type
            imp.type = parmset.match_improper_type(a1, a2, a3, a4)
            if imp.type is None:
                raise ParameterError(f"No improper type for {a1}, {a2}, {a3}, and {a4}")
            imp.type.used = False
        # prepare list of harmonic impropers present in system
        del self.improper_types[:]
        for improper in self.impropers:
            if improper.type.used:
                continue
            improper.type.used = True
            if isinstance(improper.type, ImproperType):
                self.improper_types.append(improper.type)
                improper.type.list = self.improper_types
            elif isinstance(improper.type, DihedralType):
                self.dihedral_types.append(improper.type)
                improper.type.list = self.dihedral_types
            else:
                assert False, 'Should not be here'
        # Look through the list of impropers -- if there are any periodic
        # impropers, move them over to the dihedrals list
        for i in reversed(range(len(self.impropers))):
            if isinstance(self.impropers[i].type, DihedralType):
                imp = self.impropers.pop(i)
                dih = Dihedral(imp.atom1, imp.atom2, imp.atom3, imp.atom4,
                               improper=True, ignore_end=True, type=imp.type)
                imp.delete()
                self.dihedrals.append(dih)
        # Now do the cmaps. These will not have wild-cards
        for cmap in self.cmaps:
            key = (cmap.atom1.type, cmap.atom2.type, cmap.atom3.type,
                   cmap.atom4.type, cmap.atom2.type, cmap.atom3.type,
                   cmap.atom4.type, cmap.atom5.type)
            try:
                cmap.type = parmset.cmap_types[key]
            except KeyError:
                raise ParameterError(f"No CMAP parameters found for {cmap}")
            cmap.type.used = False
        del self.cmap_types[:]
        for cmap in self.cmaps:
            if cmap.type.used: continue
            cmap.type.used = True
            self.cmap_types.append(cmap.type)
            cmap.type.list = self.cmap_types

    #===================================================

    def clear_cmap(self):
        " Clear the cmap list to prevent any CMAP parameters from being used "
        del self.cmaps[:]
