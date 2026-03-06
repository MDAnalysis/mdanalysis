"""
This module contains functionality relevant to loading a GROMACS topology file
and building a Structure from it
"""
from collections import OrderedDict, defaultdict
from contextlib import closing
import copy
from datetime import datetime
import math
import os
import re
from string import ascii_letters
import sys
import warnings

from ..constants import TINY, DEG_TO_RAD
from ..exceptions import GromacsError, GromacsWarning, ParameterError
from ..formats.registry import FileFormatType
from ..parameters import ParameterSet, _find_ureybrad_key
from ._gromacsfile import GromacsFile
from ..structure import Structure
from ..topologyobjects import (Atom, Bond, Angle, Dihedral, Improper,
            NonbondedException, ExtraPoint, BondType, Cmap, NoUreyBradley,
            AngleType, DihedralType, DihedralTypeList, ImproperType, CmapType,
            RBTorsionType, ThreeParticleExtraPointFrame, AtomType, UreyBradley,
            TwoParticleExtraPointFrame, OutOfPlaneExtraPointFrame,
            NonbondedExceptionType, UnassignedAtomType)
from ..periodic_table import element_by_mass, AtomicNum
from .. import unit as u
from ..utils.io import genopen

try:
    import pwd
    try:
        _username = pwd.getpwuid(os.getuid())[0]
    except KeyError:
        _username = 'username'
    _userid = os.getuid()
    _uname = os.uname()[1]
except ImportError:
    import getpass
    _username = getpass.getuser()   # pragma: no cover
    _userid = 0                     # pragma: no cover
    import platform                 # pragma: no cover
    _uname = platform.node()        # pragma: no cover



# Gromacs uses "funct" flags in its parameter files to indicate what kind of
# functional form is used for each of its different parameter types. This is
# taken from the topdirs.c source code file along with a table in the Gromacs
# user manual. The table below summarizes my findings, for reference:

# Bonds
# -----
#  1 - F_BONDS : simple harmonic potential
#  2 - F_G96BONDS : fourth-power potential
#  3 - F_MORSE : morse potential
#  4 - F_CUBICBONDS : cubic potential
#  5 - F_CONNBONDS : not even implemented in GROMACS
#  6 - F_HARMONIC : seems to be the same as (1) ??
#  7 - F_FENEBONDS : finietely-extensible-nonlinear-elastic (FENE) potential
#  8 - F_TABBONDS : bond function from tabulated function
#  9 - F_TABBONDSNC : bond function from tabulated function (no exclusions)
# 10 - F_RESTRBONDS : restraint bonds

# Angles
# ------
#  1 - F_ANGLES : simple harmonic potential
#  2 - F_G96ANGLES : cosine-based angle potential
#  3 - F_CROSS_BOND_BONDS : bond-bond cross term potential
#  4 - F_CROSS_BOND_ANGLES : bond-angle cross term potential
#  5 - F_UREY_BRADLEY : Urey-Bradley angle-bond potential
#  6 - F_QUARTIC_ANGLES : 4th-order polynomial potential
#  7 - F_TABANGLES : angle function from tabulated function
#  8 - F_LINEAR_ANGLES : angle function from tabulated function
#  9 - F_RESTRANGLES : restricted bending potential

# Dihedrals
# ---------
#  1 - F_PDIHS : periodic proper torsion potential [ k(1+cos(n*phi-phase)) ]
#  2 - F_IDIHS : harmonic improper torsion potential
#  3 - F_RBDIHS : Ryckaert-Bellemans torsion potential
#  4 - F_PIDIHS : periodic harmonic improper torsion potential (same as 1)
#  5 - F_FOURDIHS : Fourier dihedral torsion potential
#  8 - F_TABDIHS : dihedral potential from tabulated function
#  9 - F_PDIHS : Same as 1, but can be multi-term
# 10 - F_RESTRDIHS : Restricted torsion potential
# 11 - F_CBTDIHS : combined bending-torsion potential

_sectionre = re.compile(r'\[ (\w+) \]\s*$')

class _Defaults(object):
    """ Global properties of force fields as implemented in GROMACS """
    def __init__(self, nbfunc=1, comb_rule=2, gen_pairs='no', fudgeLJ=1.0, fudgeQQ=1.0):
        if int(nbfunc) not in (1, 2):
            raise ValueError('nbfunc must be 1 (L-J) or 2 (Buckingham)')
        if int(comb_rule) not in (1, 2, 3):
            raise ValueError('comb_rule must be 1, 2, or 3')
        if gen_pairs not in ('yes', 'no'):
            raise ValueError("gen_pairs must be 'yes' or 'no'")
        if float(fudgeLJ) < 0:
            raise ValueError('fudgeLJ must be non-negative')
        if float(fudgeQQ) < 0:
            raise ValueError('fudgeQQ must be non-negative')
        self.nbfunc = int(nbfunc)
        self.comb_rule = int(comb_rule)
        self.gen_pairs = gen_pairs
        self.fudgeLJ = float(fudgeLJ)
        self.fudgeQQ = float(fudgeQQ)

    def __repr__(self):
        return ('<_Defaults: nbfunc=%d, comb-rule=%d, gen-pairs="%s", '
                'fudgeLJ=%g, fudgeQQ=%g>' % (self.nbfunc, self.comb_rule,
                    self.gen_pairs, self.fudgeLJ, self.fudgeQQ))

    def __getitem__(self, idx):
        # Treat it like the array that it is in the topology file
        if idx < 0: idx += 5
        if idx == 0:
            return self.nbfunc
        if idx == 1:
            return self.comb_rule
        if idx == 2:
            return self.gen_pairs
        if idx == 3:
            return self.fudgeLJ
        if idx == 4:
            return self.fudgeQQ
        raise IndexError(f'Index {idx} out of range')

    def __eq__(self, other):
        return (self.nbfunc == other.nbfunc and
                self.comb_rule == other.comb_rule and
                self.gen_pairs == other.gen_pairs and
                self.fudgeLJ == other.fudgeLJ and
                self.fudgeQQ == other.fudgeQQ)

    def __setitem__(self, idx, value):
        if idx < 0: idx += 5
        if idx == 0:
            if int(value) not in (1, 2):
                raise ValueError('nbfunc must be 1 or 2')
            self.nbfunc = int(value)
        elif idx == 1:
            if int(value) not in (1, 2, 3):
                raise ValueError('comb_rule must be 1, 2, or 3')
            self.comb_rule = int(value)
        elif idx == 2:
            if value not in ('yes', 'no'):
                raise ValueError('gen_pairs must be "yes" or "no"')
            self.gen_pairs = value
        elif idx == 3:
            if float(value) < 0:
                raise ValueError('fudgeLJ must be non-negative')
            self.fudgeLJ = float(value)
        elif idx == 4:
            if float(value) < 0:
                raise ValueError('fudgeQQ must be non-negative')
            self.fudgeQQ = value
        else:
            raise IndexError(f'Index {idx} out of range')

class TopFromStructureMixin:

    @classmethod
    def from_structure(cls, struct, copy=False):
        """ Instantiates a GromacsTopologyFile instance from a Structure

        Parameters
        ----------
        struct : :class:`parmed.Structure`
            The input structure to generate from
        copy : bool, optional
            If True, assign from a *copy* of ``struct`` (this is a lot slower).
            Default is False

        Returns
        -------
        gmxtop : :class:`GromacsTopologyFile`
            The topology file defined by the given struct
        """
        from copy import copy as _copy
        clstop = cls()
        if copy:
            struct = _copy(struct)
            struct.join_dihedrals()
        clstop.atoms = struct.atoms
        clstop.residues = struct.residues
        clstop.bonds = struct.bonds
        clstop.angles = struct.angles
        clstop.dihedrals = struct.dihedrals
        clstop.impropers = struct.impropers
        clstop.cmaps = struct.cmaps
        clstop.rb_torsions = struct.rb_torsions
        clstop.urey_bradleys = struct.urey_bradleys
        clstop.adjusts = struct.adjusts
        clstop.bond_types = struct.bond_types
        clstop.angle_types = struct.angle_types
        clstop.dihedral_types = struct.dihedral_types
        clstop.improper_types = struct.improper_types
        clstop.cmap_types = struct.cmap_types
        clstop.rb_torsion_types = struct.rb_torsion_types
        clstop.urey_bradley_types = struct.urey_bradley_types
        clstop.adjust_types = struct.adjust_types
        clstop.combining_rule = struct.combining_rule
        clstop.box = struct.box
        clstop.nrexcl = struct.nrexcl
        if (struct.trigonal_angles or
                struct.out_of_plane_bends or
                struct.pi_torsions or
                struct.stretch_bends or
                struct.torsion_torsions or
                struct.chiral_frames or
                struct.multipole_frames):
            raise TypeError('GromacsTopologyFile does not support Amoeba FF')
        # Now check what the 1-4 scaling factors should be
        if hasattr(struct, 'defaults') and isinstance(struct.defaults, _Defaults):
            clstop.defaults = struct.defaults
        else:
            scee_values = set()
            scnb_values = set()
            if struct.adjusts:
                for adjust in struct.adjusts:
                    if adjust.type is None: continue
                    scee_values.add(1/adjust.type.chgscale)
                    # Do not add scnb_values, since we can just set explicit
                    # exception pair parameters in GROMACS (which this structure
                    # already has)
                # In order to specify specific pair parameters, we need to set
                # gen_pairs to 'no' so that the pair-specific L-J parameters are
                # printed to the topology file (rather than being auto-created)
                clstop.defaults.gen_pairs = 'no'
            else:
                for dihedral in struct.dihedrals:
                    if dihedral.type is None or dihedral.ignore_end: continue
                    if isinstance(dihedral.type, DihedralTypeList):
                        for dt in dihedral.type:
                            if dt.scee:
                                scee_values.add(dt.scee)
                            if dt.scnb:
                                scnb_values.add(dt.scnb)
                    else:
                        if dihedral.type.scee:
                            scee_values.add(dihedral.type.scee)
                        if dihedral.type.scnb:
                            scnb_values.add(dihedral.type.scnb)
            if len(set('%.5f' % x for x in scee_values)) > 1:
                raise GromacsError('Structure has mixed 1-4 scaling which is not supported by Gromacs')
            scee_values = list(scee_values)
            scnb_values = list(scnb_values)
            if len(set('%.5f' % x for x in scee_values)) == 1:
                clstop.defaults.fudgeQQ = 1/scee_values[0]
            else:
                clstop.defaults.fudgeQQ = 1.0
            if len(set('%.5f' % x for x in scnb_values)) == 1:
                clstop.defaults.fudgeLJ = 1/scnb_values[0]
            else:
                clstop.defaults.fudgeLJ = 1.0
        if clstop.combining_rule == 'geometric':
            clstop.defaults.comb_rule = 3

        clstop.parameterset = ParameterSet.from_structure(struct, allow_unequal_duplicates=True)
        return clstop

class GromacsTopologyFile(Structure, TopFromStructureMixin, metaclass=FileFormatType):
    """ Class providing a parser and writer for a GROMACS topology file

    Parameters
    ----------
    fname : str
        The name of the file to read
    defines : list of str=None
        If specified, this is the set of defines to use when parsing the
        topology file
    parametrized : bool, optional
        If True, parameters are assigned after parsing is done from the
        parametertypes sections. If False, only parameter types defined in the
        parameter sections themselves are loaded (i.e., on the same line as the
        parameter was defined). Default is True
    xyz : str or array, optional
        The source of atomic coordinates. It can be a string containing the name
        of a coordinate file from which to fill the coordinates (and optionally
        the unit cell information), or it can be an array with the coordinates.
        Default is None
    box : array, optional
        If provided, the unit cell information will be set from this variable.
        If provided, it must be a collection of 6 floats representing the unit
        cell dimensions a, b, c, alpha, beta, and gamma, respectively. Default
        is None.

    Notes
    -----
    If the ``xyz`` argument is a file name that contains the unit cell
    information, this unit cell information is set. However, the ``box``
    argument takes precedence and will override values given in the coordinate
    file unless it has its default value of ``None``.
    """

    #===================================================

    @staticmethod
    def id_format(filename):
        """ Identifies the file as a GROMACS topology file

        Parameters
        ----------
        filename : str
            Name of the file to check if it is a gromacs topology file

        Returns
        -------
        is_fmt : bool
            If it is identified as a gromacs topology, return True. False
            otherwise
        """
        with closing(genopen(filename)) as f:
            for line in f:
                if line.startswith(';'):
                    line = line[:line.index(';')]
                if not line.strip(): continue
                if line.startswith('#'):
                    if line.startswith('#if'): continue
                    if line.startswith('#define'): continue
                    if line.startswith('#include'): continue
                    if line.startswith('#undef'): continue
                    if line.startswith('#endif'): continue
                    return False
                rematch = _sectionre.match(line)
                if not rematch:
                    return False
                sec, = rematch.groups()
                return sec in {'atoms', 'atomtypes', 'defaults', 'moleculetype', 'system', 'bondtypes', 'angletypes',
                               'cmaptypes', 'dihedraltypes', 'bonds', 'angles', 'dihedrals', 'cmaps', 'molecules',
                               'exclusions', 'nonbond_params', 'position_restraints'}
            return False

    #===================================================

    def __init__(self, fname=None, defines=None, parametrize=True, xyz=None, box=None):
        from .. import load_file
        super(GromacsTopologyFile, self).__init__()
        self.parameterset = None
        self.defaults = _Defaults(gen_pairs='yes') # make ParmEd's default yes
        if fname is not None:
            self.read(fname, defines, parametrize)
            # Fill in coordinates and unit cell information if appropriate
            if xyz is not None:
                if isinstance(xyz, str):
                    f = load_file(xyz, skip_bonds=True)
                    if not hasattr(f, 'coordinates') or f.coordinates is None:
                        raise TypeError(f'File {xyz} does not have coordinates')
                    self.coordinates = f.coordinates
                    if box is None and hasattr(f, 'box'):
                        self.box = f.box
                else:
                    self.coordinates = xyz
            if box is not None:
                self.box = box
            self.unchange()
        elif xyz is not None or box is not None:
            raise ValueError('Cannot provide coordinates/box and NOT a top')

    #===================================================

    def read(self, fname, defines=None, parametrize=True):
        """ Reads the topology file into the current instance """
        from .. import gromacs as gmx
        params = self.parameterset = ParameterSet()
        molecules = self.molecules = dict()
        bond_types = dict()
        angle_types = dict()
        ub_types = dict()
        dihedral_types = dict()
        exc_types = dict()
        structure_contents = []
        molnames = []
        if defines is None:
            defines = OrderedDict(FLEXIBLE=1)
        proper_multiterm_dihedrals = dict()
        with closing(GromacsFile(fname, includes=[gmx.GROMACS_TOPDIR], defines=defines)) as f:
            current_section = None
            for line in f:
                line = line.strip()
                if not line: continue

                if line[0] == '[':
                    current_section = line[1:-1].strip()
                elif current_section == 'moleculetype':
                    molname, nrexcl = line.split()
                    nrexcl = int(nrexcl)
                    if molname in molecules:
                        raise GromacsError(f'Duplicate definition of molecule {molname}')
                    molecule = Structure()
                    molecules[molname] = (molecule, nrexcl)
                    molnames.append(molname)
                    molecule.nrexcl = nrexcl
                    bond_types = dict()
                    angle_types = dict()
                    ub_types = dict()
                    dihedral_types = dict()
                    exc_types = dict()
                elif current_section == 'atoms':
                    atom, resname, resnum, icode = self._parse_atoms(line, params)
                    molecule.add_atom(atom, resname, resnum, inscode=icode)
                elif current_section == 'bonds':
                    bond, bond_type = self._parse_bonds(line, bond_types, molecule.atoms)
                    molecule.bonds.append(bond)
                    if bond_type is not None:
                        molecule.bond_types.append(bond_type)
                        bond_type.list = molecule.bond_types
                elif current_section == 'pairs':
                    nbe, nbet = self._parse_pairs(line, exc_types, molecule.atoms)
                    molecule.adjusts.append(nbe)
                    if nbet is not None:
                        molecule.adjust_types.append(nbet)
                        nbet.list = molecule.adjust_types
                elif current_section == 'angles':
                    ang, ub, angt, ubt = self._parse_angles(line, angle_types, ub_types,
                                                            molecule.atoms)
                    molecule.angles.append(ang)
                    if ub is not None:
                        molecule.urey_bradleys.append(ub)
                    if angt is not None:
                        molecule.angle_types.append(angt)
                        angt.list = molecule.angle_types
                    if ubt is not None and ubt is not NoUreyBradley:
                        molecule.urey_bradley_types.append(ubt)
                        ubt.list = molecule.urey_bradley_types
                elif current_section == 'dihedrals':
                    self._parse_dihedrals(line, dihedral_types, proper_multiterm_dihedrals,
                                          molecule)
                elif current_section == 'cmap':
                    cmap = self._parse_cmaps(line, molecule.atoms)
                    molecule.cmaps.append(cmap)
                elif current_section == 'system':
                    self.title = line
                elif current_section == 'defaults':
                    words = line.split()
                    if len(words) < 2: # 3, 4, and 5 fields are optional
                        raise GromacsError('Too few fields in [ defaults ]')
                    if words[0] != '1':
                        warnings.warn('Unsupported nonbonded type; unknown functional',
                                      GromacsWarning)
                        self.unknown_functional = True
                    if words[1] in ('1', '3'):
                        self.combining_rule = 'geometric'
                    self.defaults = _Defaults(*words)
                elif current_section == 'molecules':
                    name, num = line.split()
                    num = int(num)
                    structure_contents.append((name, num))
                elif current_section == 'settles':
                    bnds, bndts = self._parse_settles(line, molecule.atoms)
                    molecule.bonds.extend(bnds)
                    molecule.bond_types.extend(bndts)
                    molecule.bond_types.claim()
                elif current_section in ('virtual_sites3', 'dummies3'):
                    try:
                        b, bt = self._parse_vsites3(line, molecule.atoms, params)
                    except KeyError:
                        raise GromacsError('Cannot determine vsite geometry '
                                           'without parameter types')
                    molecule.bonds.append(b)
                    molecule.bond_types.append(bt)
                    bt.list = molecule.bond_types
                elif current_section == 'exclusions':
                    atoms = [molecule.atoms[int(w)-1] for w in line.split()]
                    for a in atoms[1:]:
                        atoms[0].exclude(a)
                elif current_section == 'atomtypes':
                    attype, typ = self._parse_atomtypes(line)
                    params.atom_types[attype] = typ
                elif current_section == 'nonbond_params':
                    words = line.split()
                    a1, a2 = words[:2]
#                   func = int(words[2]) #... unused
                    sig, eps = (float(x) for x in words[3:5])
                    sig *= 10 # Convert to Angstroms
                    eps *= u.kilojoule.conversion_factor_to(u.kilocalorie)
                    params.nbfix_types[(a1, a2)] = (eps, sig*2**(1/6))
                    params.nbfix_types[(a2, a1)] = (eps, sig*2**(1/6))
                    params.atom_types[a1].add_nbfix(a2, sig*2**(1/6), eps)
                    params.atom_types[a2].add_nbfix(a1, sig*2**(1/6), eps)
                elif current_section == 'bondtypes':
                    a, b, t = self._parse_bondtypes(line)
                    params.bond_types[(a, b)] = t
                    params.bond_types[(b, a)] = t
                elif current_section == 'angletypes':
                    a, b, c, t, ut = self._parse_angletypes(line)
                    params.angle_types[(a, b, c)] = t
                    params.angle_types[(c, b, a)] = t
                    if ut is not None:
                        params.urey_bradley_types[(a, b, c)] = ut
                        params.urey_bradley_types[(c, b, a)] = ut
                elif current_section == 'dihedraltypes':
                    key, knd, t, replace = self._parse_dihedraltypes(line)
                    rkey = tuple(reversed(key))
                    if knd == 'normal':
                        if replace or key not in params.dihedral_types:
                            t = DihedralTypeList([t])
                            params.dihedral_types[key] = t
                            params.dihedral_types[rkey] = t
                        elif key in params.dihedral_types:
                            params.dihedral_types[key].append(t, override=True)
                    elif knd == 'improper':
                        params.improper_types[key] = t
                    elif knd == 'improper_periodic':
                        params.improper_periodic_types[key] = t
                        params.improper_periodic_types[rkey] = t
                    elif knd == 'rbtorsion':
                        params.rb_torsion_types[key] = t
                        params.rb_torsion_types[rkey] = t
                elif current_section == 'cmaptypes':
                    a1, a2, a3, a4, a5, t = self._parse_cmaptypes(line)
                    params.cmap_types[(a1, a2, a3, a4, a2, a3, a4, a5)] = t
                    params.cmap_types[(a5, a4, a3, a2, a4, a3, a2, a1)] = t
                elif current_section == 'pairtypes':
                    a, b, t = self._parse_pairtypes(line)
                    params.pair_types[(a, b)] = params.pair_types[(b, a)] = t
            itplist = f.included_files

        # If the file did not contain the molecules section, perhaps
        # because it was an itp-file. We assume that each molecule loaded
        # should be contained once in this structure
        if not structure_contents :
            for name in molnames :
                structure_contents.append((name, 1))

        # Combine first, then parametrize. That way, we don't have to create
        # copies of the ParameterType instances in self.parameterset
        for molname, num in structure_contents:
            if molname not in molecules:
                raise GromacsError('Structure contains %s molecules, but no '
                                   'template defined' % molname)
            molecule, nrexcl = molecules[molname]
            if nrexcl < 3 and _any_atoms_farther_than(molecule, nrexcl):
                warnings.warn('nrexcl %d not currently supported' % nrexcl,
                              GromacsWarning)
            elif nrexcl > 3 and _any_atoms_farther_than(molecule, 3):
                warnings.warn('nrexcl %d not currently supported' % nrexcl,
                              GromacsWarning)
            if num == 0:
                warnings.warn('Detected addition of 0 %s molecules in topology '
                              'file' % molname, GromacsWarning)
            if num == 1:
                self += molecules[molname][0]
            elif num > 1:
                self += molecules[molname][0] * num
            else:
                raise GromacsError("Can't add %d %s molecules" % (num, molname))
        self.itps = itplist
        if parametrize:
            self.parametrize()

    #===================================================

    # Private parsing helper functions

    def _parse_atoms(self, line, params):
        """ Parses an atom line. Returns an Atom, resname, resnum """
        words = line.split()
        try:
            attype = params.atom_types[words[1]]
        except KeyError:
            attype = None
        if len(words) < 8:
            if attype is not None:
                mass = attype.mass
                atomic_number = attype.atomic_number
            else:
                mass = -1
                atomic_number = -1
        else:
            mass = float(words[7])
            if attype is not None and attype.atomic_number >= 0:
                atomic_number = attype.atomic_number
            else:
                atomic_number = AtomicNum[element_by_mass(mass)]
        charge = float(words[6]) if len(words) > 6 else None
        if atomic_number == 0:
            atom = ExtraPoint(name=words[4], type=words[1], charge=charge)
        else:
            atom = Atom(atomic_number=atomic_number, name=words[4],
                        type=words[1], charge=charge, mass=mass)

        # check for insertion code and negative res number
        if words[2].isnumeric() or (words[2].startswith('-') and words[2][1:].isnumeric()):
            icode = ''
            resnum = int(words[2])
        else:
            icode = words[2][-1]
            resnum = int(words[2][:-1])
        return atom, words[3], resnum, icode

    def _parse_bonds(self, line, bond_types, atoms):
        """ Parses a bond line. Returns a Bond, BondType/None """
        words = line.split()
        i, j = int(words[0])-1, int(words[1])-1
        funct = int(words[2])
        if funct != 1:
            warnings.warn('bond funct != 1; unknown functional',
                          GromacsWarning)
            self.unknown_functional = True
        bond = Bond(atoms[i], atoms[j])
        bond.funct = funct
        bond_type = None
        if len(words) >= 5 and funct == 1:
            req, k = (float(x) for x in words[3:5])
            if (req, k) in bond_types:
                bond.type = bond_types[(req, k)]
            else:
                bond_type = BondType(
                        k*u.kilojoule_per_mole/u.nanometer**2/2,
                        req*u.nanometer
                )
                bond_types[(req, k)] = bond.type = bond_type
        return bond, bond_type

    def _parse_pairs(self, line, exc_types, atoms):
        """ Parses a pairs line. Returns NonbondedException, NEType/None """
        words = line.split()
        i, j = int(words[0])-1, int(words[1])-1
        funct = int(words[2])
        if funct != 1:
            # This is not even supported in Gromacs
            warnings.warn('pairs funct != 1; unknown functional',
                          GromacsWarning)
            self.unknown_functional = True
        nbe = NonbondedException(atoms[i], atoms[j])
        nbe.funct = funct
        nbet = None
        if funct == 1 and len(words) >= 5:
            sig = float(words[3]) * 2**(1/6)
            eps = float(words[4])
            if (sig, eps) in exc_types:
                nbe.type = exc_types[(sig, eps)]
            else:
                nbet = NonbondedExceptionType(sig*u.nanometers, eps*u.kilojoules_per_mole,
                                              self.defaults.fudgeQQ)
                exc_types[(sig, eps)] = nbe.type = nbet
        return nbe, nbet

    def _parse_angles(self, line, angle_types, ub_types, atoms):
        """ Parse an angles line, Returns Angle, UB/None, and types """
        words = line.split()
        i, j, k = [int(w)-1 for w in words[:3]]
        funct = int(words[3])
        if funct not in (1, 5):
            warnings.warn('angles funct != 1 or 5; unknown '
                          'functional', GromacsWarning)
            self.unknown_functional = True
        angt = ub = ubt = None
        ang = Angle(atoms[i], atoms[j], atoms[k])
        ang.funct = funct
        if funct == 5:
            ub = UreyBradley(atoms[i], atoms[k])
        if (funct == 1 and len(words) >= 6) or (funct == 5 and len(words) >= 8):
            theteq, k = (float(x) for x in words[4:6])
            if (theteq, k) in angle_types:
                ang.type = angle_types[(theteq, k)]
            else:
                angt = AngleType(k*u.kilojoule_per_mole/u.radian**2/2,
                                 theteq*u.degree)
                angle_types[(theteq, k)] = ang.type = angt
        if funct == 5 and len(words) >= 8:
            ubreq, ubk = (float(x) for x in words[6:8])
            if ubk > 0:
                if (ubreq, ubk) in ub_types:
                    ub.type = ub_types[(ubreq, ubk)]
                else:
                    ubt = BondType(
                        ubk*u.kilojoule_per_mole/u.nanometer**2/2,
                        ubreq*u.nanometer,
                    )
                    ub_types[(ubreq, ubk)] = ub.type = ubt
            else:
                ub.type = NoUreyBradley
        return ang, ub, angt, ubt

    def _parse_dihedrals(self, line, dihedral_types, PMD, molecule):
        """ Processes a dihedrals line, returns None """
        words = line.split()
        i, j, k, l = [int(x)-1 for x in words[:4]]
        funct = int(words[4])
        if funct in (1, 4) or (funct == 9 and len(words) < 8):
            dih, diht = self._process_normal_dihedral(words, molecule.atoms, i,
                                                      j, k, l, dihedral_types,
                                                      funct==4)
            molecule.dihedrals.append(dih)
            if diht is not None:
                molecule.dihedral_types.append(diht)
                diht.list = molecule.dihedral_types
        elif funct == 2:
            dih, impt = self._process_improper(words, i, j, k, l,
                                               molecule.atoms, dihedral_types)
            molecule.impropers.append(dih)
            if impt is not None:
                molecule.improper_types.append(impt)
                impt.list = molecule.improper_types
        elif funct == 3:
            dih, rbt = self._process_rbtorsion(words, i, j, k, l, molecule.atoms,
                                              dihedral_types)
            molecule.rb_torsions.append(dih)
            if rbt is not None:
                molecule.rb_torsion_types.append(rbt)
                rbt.list = molecule.rb_torsion_types
        elif funct == 9:
            # in-line parameters, since len(words) must be >= 8
            key = (molecule.atoms[i], molecule.atoms[j],
                   molecule.atoms[k], molecule.atoms[l])
            if key in PMD:
                diht = PMD[key]
                self._process_dihedral_series(words, diht)
                dih = None
            else:
                dih = Dihedral(*key)
                diht = self._process_dihedral_series(words)
                dih.type = PMD[key] = PMD[tuple(reversed(key))] = diht
                molecule.dihedrals.append(dih)
                molecule.dihedral_types.append(diht)
                diht.list = molecule.dihedral_types
        else:
            # ??? unknown funct
            warnings.warn('torsions funct != 1, 2, 3, 4, 9; unknown'
                          ' functional', GromacsWarning)
            dih = Dihedral(molecule.atoms[i], molecule.atoms[j],
                           molecule.atoms[k], molecule.atoms[l])
            molecule.dihedrals.append(dih)
            self.unknown_functional = True

        if dih is not None:
            dih.funct = funct

    def _parse_cmaps(self, line, atoms):
        """ Parses cmap terms, returns cmap """
        words = line.split()
        i, j, k, l, m = (int(w)-1 for w in words[:5])
        funct = int(words[5])
        if funct != 1:
            warnings.warn('cmap funct != 1; unknown functional',
                          GromacsWarning)
            self.unknown_functional = True
        cmap = Cmap(atoms[i], atoms[j], atoms[k], atoms[l], atoms[m])
        cmap.funct = funct
        return cmap

    def _parse_settles(self, line, atoms):
        """ Parses settles line; returns list of Bonds, list of BondTypes """
        # Instead of adding bonds that get constrained for waters (or other
        # 3-atom molecules), GROMACS uses a "settles" section to specify the
        # constraint geometry. We have to translate that into bonds.
        natoms = len([a for a in atoms if not isinstance(a, ExtraPoint)])
        if natoms != 3:
            raise GromacsError("Cannot SETTLE a %d-atom molecule" % natoms)
        try:
            oxy, = [atom for atom in atoms if atom.atomic_number == 8]
            hyd1, hyd2 = [atom for atom in atoms if atom.atomic_number == 1]
        except ValueError:
            raise GromacsError('Can only SETTLE water; wrong atoms')
        #TODO see if there's a bond_type entry in the parameter set
        #     that we can fill in? Wait until this is needed...
        try:
            i, funct, doh, dhh = line.split()
            doh, dhh = float(doh), float(dhh)
        except ValueError:
            raise GromacsError('Bad [ settles ] line')
        nm = u.nanometers
        bt_oh = BondType(5e5*u.kilojoules_per_mole/nm**2, doh*nm)
        bt_hh = BondType(5e5*u.kilojoules_per_mole/nm**2, dhh*nm)
        return [Bond(oxy, hyd1, bt_oh), Bond(oxy, hyd2, bt_oh),
                Bond(hyd1, hyd2, bt_hh)], [bt_oh, bt_hh]

    def _parse_vsites3(self, line, all_atoms, params):
        """ Parse vsites3/dummy3 line; returns Bond, BondType """
        words = line.split()
        vsite = all_atoms[int(words[0])-1]
        atoms = [all_atoms[int(i)-1] for i in words[1:4]]
        funct = int(words[4])
        if funct == 1:
            a, b = float(words[5]), float(words[6])
            if abs(a - b) > TINY:
                raise GromacsError("No vsite frames with different weights")
        else:
            raise GromacsError('Only 3-point vsite type 1 is supported')
        # We need to know the geometry of the frame in order to
        # determine the bond length between the virtual site and its
        # parent atom
        parent = atoms[0]
        if vsite in parent.bond_partners:
            raise GromacsError('Unexpected bond b/w vsite and its parent')
        kws = dict()
        for bond in parent.bonds:
            if atoms[1] in bond:
                key = (_gettype(parent), _gettype(atoms[1]))
                kws['dp1'] = (bond.type or params.bond_types[key]).req
            if atoms[2] in bond:
                key = (_gettype(bond.atom1), _gettype(bond.atom2))
                kws['dp2'] = (bond.type or params.bond_types[key]).req
        for angle in parent.angles:
            if parent is not angle.atom2: continue
            if atoms[0] not in angle or atoms[1] not in angle: continue
            key = (_gettype(angle.atom1), _gettype(angle.atom2),
                   _gettype(angle.atom3))
            kws['theteq'] = (angle.type or params.angle_types[key]).theteq
            break
        else: # Did not break, no theta found
            for bond in atoms[1].bonds:
                if atoms[2] in bond:
                    key = (_gettype(bond.atom1), _gettype(bond.atom2))
                    kws['d12'] = (bond.type or params.bond_types[key]).req
        bondlen = ThreeParticleExtraPointFrame.from_weights(parent, atoms[1],
                                                    atoms[2], a, b, **kws)
        bt_vs = BondType(0, bondlen*u.angstroms)
        return Bond(vsite, parent, bt_vs), bt_vs

    def _parse_atomtypes(self, line):
        """ Parses line from atomtypes section, returns str, AtomType """
        words = line.split()
        # Support the following spec, found in the Gromacs source code:
        # Field 0 (mandatory) : nonbonded type name (string)
        # Field 1 (optional)  : bonded type (string)
        # Field 2 (optional)  : atomic number (int)
        # Field 3 (mandatory) : mass (float)
        # Field 4 (mandatory) : charge (float)
        # Field 5 (mandatory) : particle type (single character)
        attype = words[0]
        if len(words[3]) == 1 and words[3] in ascii_letters:
            # Field 1 and Field 2 are both missing
            atnum = -1
            sigidx = 4
#           ptypeidx = 3 # ... unused
            massidx = 1
            bond_type = None
        elif len(words[5]) == 1 and words[5] in ascii_letters:
            # Both Field 1 and Field 2 are present
            sigidx = 6
#           ptypeidx = 5 # ... unused
            massidx = 3
            atnum = int(words[2])
            bond_type = words[1]
        else:
            # One of Field 1 or 2 are missing
#           ptypeidx = 4 # ... unused
            massidx = 2
            sigidx = 5
            try:
                atnum = int(words[1])
                bond_type = None
            except ValueError:
                # This must be a bonded type string
                bond_type = words[1]
                atnum = -1
        mass = float(words[massidx])
        if mass > 0 and atnum == -1:
            atnum = AtomicNum[element_by_mass(mass)]
        chg = float(words[massidx+1])
#       ptype = words[ptypeidx] # ... unused
        sig = float(words[sigidx]) * u.nanometers
        eps = float(words[sigidx+1]) * u.kilojoules_per_mole
        typ = AtomType(attype, None, mass, atnum, bond_type=bond_type, charge=chg)
        typ.set_lj_params(eps, sig*2**(1/6)/2)
        return attype, typ

    def _parse_bondtypes(self, line):
        """ Parse bondtypes line. Returns str, str, BondType """
        words = line.split()
        r = float(words[3]) * u.nanometers
        k = (float(words[4]) / 2) * (u.kilojoules_per_mole / u.nanometers**2)
        if words[2] != '1':
            warnings.warn('bondtypes funct != 1; unknown functional',
                          GromacsWarning)
            self.unknown_functional = True
        return words[0], words[1], BondType(k, r)

    def _parse_angletypes(self, line):
        """
        Parses angletypes line. Returns str, str, str, AngleType, BondType/None
        """
        words = line.split()
        theta = float(words[4]) * u.degrees
        k = (float(words[5]) / 2) * (u.kilojoules_per_mole / u.radians**2)
        if words[3] != '1' and words[3] != '5':
            warnings.warn('angletypes funct != 1 or 5; unknown functional',
                          GromacsWarning)
            self.unknown_functional = True
        ub = None
        if words[3] == '5':
            # Contains the angle with urey-bradley
            ub0 = float(words[6])
            cub = float(words[7]) / 2
            if cub == 0:
                ub = NoUreyBradley
            else:
                ub0 *= u.nanometers
                cub *= u.kilojoules_per_mole / u.nanometers**2
                ub = BondType(cub, ub0)
        return words[0], words[1], words[2], AngleType(k, theta), ub

    def _parse_dihedraltypes(self, line):
        """ Parse dihedraltypes, returns (str,str,str,str), str, Type, bool """
        words = line.split()
        replace = False
        dtype = 'normal'
        # Ugh. Gromacs allows only two atom types (the middle atom types) to be
        # specified. This signifies wild-cards
        if words[2] in ('1', '2', '3', '4', '5', '8', '9', '10', '11'):
            a1 = a4 = 'X'
            a2, a3 = words[:2]
            si = 2
        else:
            a1, a2, a3, a4 = words[:4]
            si = 4
        improper_periodic = False
        replace = words[si] in ('1', '2', '3', '4')
        improper_periodic = words[si] == '4'
        if words[si] == '2':
            dtype = 'improper'
        elif words[si] == '3':
            dtype = 'rbtorsion'
        elif words[si] not in ('1', '4', '9'):
            warnings.warn('dihedraltypes funct not supported', GromacsWarning)
            self.unknown_functional = True
        # Do the proper types
        if dtype == 'normal':
            phase = float(words[si+1]) * u.degrees
            phi_k = float(words[si+2]) * u.kilojoules_per_mole
            per = int(words[si+3])
            ptype = DihedralType(phi_k, per, phase,
                                 scee=1/self.defaults.fudgeQQ,
                                 scnb=1/self.defaults.fudgeLJ)
            if improper_periodic:
                # must do this here, since dtype has to be 'normal' above
                dtype = 'improper_periodic'
        elif dtype == 'improper':
            theta = float(words[si+1])*u.degrees
            k = float(words[si+2])*u.kilojoules_per_mole/u.radians**2/2
            a1, a2, a3, a4 = sorted([a1, a2, a3, a4])
            ptype = ImproperType(k, theta)
        elif dtype == 'rbtorsion':
            a1, a2, a3, a4 = words[:4]
            c0, c1, c2, c3, c4, c5 = (float(x)*u.kilojoules_per_mole
                                        for x in words[si+1:si+7])
            ptype = RBTorsionType(c0, c1, c2, c3, c4, c5,
                                  scee=1/self.defaults.fudgeQQ,
                                  scnb=1/self.defaults.fudgeLJ)
        return (a1, a2, a3, a4), dtype, ptype, replace

    def _parse_cmaptypes(self, line):
        words = line.split()
        a1, a2, a3, a4, a5 = words[:5]
#       funct = int(words[5]) # ... unused
        res1, res2 = int(words[6]), int(words[7])
        grid = [float(w) for w in words[8:]] * u.kilojoules_per_mole
        if len(grid) != res1 * res2:
            raise GromacsError('CMAP grid dimensions do not match resolution')
        if res1 != res2:
            raise GromacsError('Only square CMAPs are supported')
        return a1, a2, a3, a4, a5, CmapType(res1, grid)

    def _parse_pairtypes(self, line):
        words = line.split()
        a1, a2 = words[:2]
#       funct = int(words[2]) # ... unused
        cs6, cs12 = (float(x) for x in words[3:5])
        cs6 *= u.nanometers * 2**(1/6)
        cs12 *= u.kilojoules_per_mole
        return a1, a2, NonbondedExceptionType(cs6, cs12, self.defaults.fudgeQQ)

    #===================================================

    # Internal Dihedral processing routines for different kinds of dihedrals

    def _process_normal_dihedral(self, words, atoms, i, j, k, l,
                                 dihedral_types, imp):
        dih = Dihedral(atoms[i], atoms[j], atoms[k], atoms[l], improper=imp)
        diht = None
        if len(words) >= 8:
            phase, phi_k, per = (float(x) for x in words[5:8])
            if (phase, phi_k, per) in dihedral_types:
                dih.type = dihedral_types[(phase, phi_k, per)]
            else:
                diht = DihedralType(phi_k*u.kilojoule_per_mole,
                                    per, phase*u.degrees,
                                    scee=1/self.defaults.fudgeQQ,
                                    scnb=1/self.defaults.fudgeLJ)
                dihedral_types[(phase, phi_k, per)] = dih.type = diht
        return dih, diht

    def _process_dihedral_series(self, words, dihtype=None):
        phase, phi_k, per = (float(x) for x in words[5:8])
        dt = DihedralType(phi_k*u.kilojoule_per_mole,
                          per, phase*u.degrees,
                          scee=1/self.defaults.fudgeQQ,
                          scnb=1/self.defaults.fudgeLJ)
        if dihtype is not None:
            dihtype.append(dt)
            dtl = None
        else:
            dt = DihedralType(phi_k*u.kilojoule_per_mole,
                              per, phase*u.degrees,
                              scee=1/self.defaults.fudgeQQ,
                              scnb=1/self.defaults.fudgeLJ)
            dtl = DihedralTypeList()
            dtl.append(dt)
        return dtl

    def _process_improper(self, words, i, j, k, l, atoms, dihedral_types):
        """ Processes an improper, returns Improper, ImproperType """
        # Improper
        imp = Improper(atoms[i], atoms[j], atoms[k], atoms[l])
        impt = None
        if len(words) >= 7:
            psieq, k = (float(x) for x in words[5:7])
            if (psieq, k) in dihedral_types:
                imp.type = dihedral_types[(psieq, k)]
            else:
                impt = ImproperType(k*u.kilojoule_per_mole/u.radian**2/2,
                                    psieq*u.degree)
                imp.type = dihedral_types[(psieq, k)] = impt
        return imp, impt

    def _process_rbtorsion(self, words, i, j, k, l, atoms, dihedral_types):
        rb = Dihedral(atoms[i], atoms[j], atoms[k], atoms[l])
        rbt = None
        if len(words) >= 11:
            c0, c1, c2, c3, c4, c5 = (float(x) for x in words[5:11])
            if (c0, c1, c2, c3, c4, c5) in dihedral_types:
                rb.type = dihedral_types[(c0, c1, c2, c3, c4, c5)]
            else:
                kjpm = u.kilojoules_per_mole
                rbt = RBTorsionType(c0*kjpm, c1*kjpm, c2*kjpm,
                                    c3*kjpm, c4*kjpm, c5*kjpm,
                                    scee=1/self.defaults.fudgeQQ,
                                    scnb=1/self.defaults.fudgeLJ)
                dihedral_types[(c0, c1, c2, c3, c4, c5)] = rb.type = rbt
        return rb, rbt

    #===================================================

    def parametrize(self):
        """
        Assign parameters to the current structure. This should be called
        *after* `read`
        """
        if self.parameterset is None:
            raise RuntimeError('parametrize called before read')
        params = copy.copy(self.parameterset)
        def update_typelist_from(ptypes, types):
            added_types = set(id(typ) for typ in types)
            for k, typ in ptypes.items():
                if not typ.used: continue
                if id(typ) in added_types: continue
                added_types.add(id(typ))
                types.append(typ)
            types.claim()
        # Assign all of the parameters. If they've already been assigned (i.e.,
        # on the parameter line itself) keep the existing parameters
        for atom in self.atoms:
            atom.atom_type = params.atom_types[atom.type]
        # The list of ordered 2-tuples of atoms explicitly specified in [ pairs ].
        # Under most circumstances, this is the list of 1-4 pairs.
        gmx_pair = set()
        for pair in self.adjusts:
            if pair.atom1 > pair.atom2:
                gmx_pair.add((pair.atom2, pair.atom1))
            else:
                gmx_pair.add((pair.atom1, pair.atom2))
            if pair.type is not None: continue
            key = (_gettype(pair.atom1), _gettype(pair.atom2))
            if key in params.pair_types:
                pair.type = params.pair_types[key]
                pair.type.used = True
            elif self.defaults.gen_pairs == 'yes':
                assert self.combining_rule in ('geometric', 'lorentz'), \
                        'Unrecognized combining rule'
                if self.combining_rule == 'geometric':
                    eps = math.sqrt(pair.atom1.epsilon * pair.atom2.epsilon)
                    sig = math.sqrt(pair.atom1.sigma * pair.atom2.sigma)
                elif self.combining_rule == 'lorentz':
                    eps = math.sqrt(pair.atom1.epsilon * pair.atom2.epsilon)
                    sig = 0.5 * (pair.atom1.sigma + pair.atom2.sigma)
                eps *= self.defaults.fudgeLJ
                pairtype = NonbondedExceptionType(sig*2**(1/6), eps,
                            self.defaults.fudgeQQ, list=self.adjust_types)
                self.adjust_types.append(pairtype)
                pair.type = pairtype
                pair.type.used = True
            else:
                raise ParameterError('Not all pair parameters can be found')
        update_typelist_from(params.pair_types, self.adjust_types)
        # This is the list of 1-4 pairs determined from the bond graph.
        # If this is different from what's in [ pairs ], we print a warning
        # and make some adjustments (specifically, other programs assume
        # the 1-4 list is complete, so we zero out the parameters for
        # 1-4 pairs that aren't in [ pairs ].
        true_14 = set()
        for bond in self.bonds:
            for bpi in bond.atom1.bond_partners:
                for bpj in bond.atom2.bond_partners:
                    if len(set([bpi, bond.atom1, bond.atom2, bpj])) < 4:
                        continue
                    if bpi in bpj.bond_partners or bpi in bpj.angle_partners:
                        continue
                    if bpi > bpj:
                        true_14.add((bpj, bpi))
                    else:
                        true_14.add((bpi, bpj))
            if bond.type is not None: continue
            key = (_gettype(bond.atom1), _gettype(bond.atom2))
            if key in params.bond_types:
                bond.type = params.bond_types[key]
                bond.type.used = True
            else:
                raise ParameterError('Not all bond parameters found')
        if len(true_14 - gmx_pair) > 0:
            zero_pairtype = NonbondedExceptionType(0.0, 0.0, 0.0,
                                                   list=self.adjust_types)
            self.adjust_types.append(zero_pairtype)
            num_zero_14 = 0
            for a1, a2 in (true_14 - gmx_pair):
                self.adjusts.append(NonbondedException(a1, a2, zero_pairtype))
                num_zero_14 += 1
            warnings.warn('%i 1-4 pairs were missing from the [ pairs ] '
                          'section and were set to zero; make sure you '
                          'know what you\'re doing!' % num_zero_14,
                          GromacsWarning)
        if len(gmx_pair - true_14) > 0:
            warnings.warn('The [ pairs ] section contains %i exceptions that '
                          'aren\'t 1-4 pairs; make sure you know what '
                          'you\'re doing!' % (len(gmx_pair - true_14)),
                          GromacsWarning)
        update_typelist_from(params.bond_types, self.bond_types)
        for angle in self.angles:
            if angle.type is not None: continue
            key = (_gettype(angle.atom1), _gettype(angle.atom2),
                   _gettype(angle.atom3))
            if key in params.angle_types:
                angle.type = params.angle_types[key]
                angle.type.used = True
            else:
                raise ParameterError('Not all angle parameters found')
        update_typelist_from(params.angle_types, self.angle_types)
        for ub in self.urey_bradleys:
            if ub.type is not None: continue
            key = _find_ureybrad_key(ub)
            if key in params.urey_bradley_types:
                ub.type = params.urey_bradley_types[key]
                if ub.type is not NoUreyBradley:
                    ub.type.used = True
            else:
                raise ParameterError('Not all urey-bradley parameters found')
        # Now strip out all of the Urey-Bradley terms whose parameters are 0
        for i in reversed(range(len(self.urey_bradleys))):
            if self.urey_bradleys[i].type is NoUreyBradley:
                del self.urey_bradleys[i]
        update_typelist_from(params.urey_bradley_types, self.urey_bradley_types)
        for t in self.dihedrals:
            if t.type is not None: continue
            key = (_gettype(t.atom1), _gettype(t.atom2), _gettype(t.atom3),
                   _gettype(t.atom4))
            if not t.improper:
                wckey = ('X', _gettype(t.atom2), _gettype(t.atom3), 'X')
                wckey1 = (_gettype(t.atom1), _gettype(t.atom2),
                          _gettype(t.atom3), 'X')
                wckey2 = ('X', _gettype(t.atom2), _gettype(t.atom3),
                          _gettype(t.atom4))
                if key in params.dihedral_types:
                    t.type = params.dihedral_types[key]
                    t.type.used = True
                elif wckey1 in params.dihedral_types:
                    t.type = params.dihedral_types[wckey1]
                    t.type.used = True
                elif wckey2 in params.dihedral_types:
                    t.type = params.dihedral_types[wckey2]
                    t.type.used = True
                elif wckey in params.dihedral_types:
                    t.type = params.dihedral_types[wckey]
                    t.type.used = True
                else:
                    raise ParameterError('Not all torsion parameters found')
            else:
                if key in params.improper_periodic_types:
                    t.type = params.improper_periodic_types[key]
                    t.type.used = True
                else:
                    for wckey in [(key[0],key[1],key[2],'X'),
                                  ('X',key[1],key[2],key[3]),
                                  (key[0],key[1],'X','X'),
                                  ('X','X',key[2],key[3])]:
                        if wckey in params.improper_periodic_types:
                            t.type = params.improper_periodic_types[wckey]
                            t.type.used = True
                            break
                    else:
                        raise ParameterError('Not all improper torsion '
                                             'parameters found')
        update_typelist_from(params.dihedral_types, self.dihedral_types)
        update_typelist_from(params.improper_periodic_types, self.dihedral_types)
        for t in self.rb_torsions:
            if t.type is not None: continue
            key = (_gettype(t.atom1), _gettype(t.atom2), _gettype(t.atom3),
                   _gettype(t.atom4))
            wckey = ('X', _gettype(t.atom2), _gettype(t.atom3), 'X')
            wckey1 = (_gettype(t.atom1), _gettype(t.atom2),
                      _gettype(t.atom3), 'X')
            wckey2 = ('X', _gettype(t.atom2), _gettype(t.atom3),
                      _gettype(t.atom4))
            if key in params.rb_torsion_types:
                t.type = params.rb_torsion_types[key]
                t.type.used = True
            elif wckey1 in params.rb_torsion_types:
                t.type = params.rb_torsion_types[wckey1]
                t.type.used = True
            elif wckey2 in params.rb_torsion_types:
                t.type = params.rb_torsion_types[wckey2]
                t.type.used = True
            elif wckey in params.rb_torsion_types:
                t.type = params.rb_torsion_types[wckey]
                t.type.used = True
            else:
                raise ParameterError('Not all R-B torsion parameters found')
        update_typelist_from(params.rb_torsion_types, self.rb_torsion_types)
        self.update_dihedral_exclusions()
        for t in self.impropers:
            if t.type is not None: continue
            key = tuple(sorted([_gettype(t.atom1), _gettype(t.atom2),
                                _gettype(t.atom3), _gettype(t.atom4)]))
            if key in params.improper_types:
                t.type = params.improper_types[key]
                t.type.used = True
                continue
            # Now we will try to find a compatible wild-card... the first atom
            # is the central atom. So take each of the other three and plug that
            # one in
            for anchor in (_gettype(t.atom2), _gettype(t.atom3),
                           _gettype(t.atom4)):
                wckey = tuple(sorted([_gettype(t.atom1), anchor, 'X', 'X']))
                if wckey not in params.improper_types: continue
                t.type = params.improper_types[wckey]
                t.type.used = True
                break
            else:
                raise ParameterError('Not all improper parameters found')
        update_typelist_from(params.improper_types, self.improper_types)
        for c in self.cmaps:
            if c.type is not None: continue
            key = (_gettype(c.atom1), _gettype(c.atom2), _gettype(c.atom3),
                    _gettype(c.atom4), _gettype(c.atom5))
            key = (key[0],key[1],key[2],key[3],key[1],key[2],key[3],key[4])
            if key in params.cmap_types:
                c.type = params.cmap_types[key]
                c.type.used = True
            else:
                raise ParameterError('Not all cmap parameters found')
        update_typelist_from(params.cmap_types, self.cmap_types)

    #===================================================

    def copy(self, cls, split_dihedrals=False):
        """
        Makes a copy of the current structure as an instance of a specified
        subclass

        Parameters
        ----------
        cls : Structure subclass
            The returned object is a copy of this structure as a `cls` instance
        split_dihedrals : ``bool``
            If True, then the Dihedral entries will be split up so that each one
            is paired with a single DihedralType (rather than a
            DihedralTypeList)

        Returns
        -------
        *cls* instance
            The instance of the Structure subclass `cls` with a copy of the
            current Structure's topology information
        """
        c = super(GromacsTopologyFile, self).copy(cls, split_dihedrals)
        c.defaults = copy.copy(self.defaults)
        return c

    #===================================================

    def __getitem__(self, selection):
        """ See Structure.__getitem__ for documentation """
        # Make sure defaults is properly copied
        struct = super(GromacsTopologyFile, self).__getitem__(selection)
        if isinstance(struct, Atom):
            return struct
        struct.defaults = copy.copy(self.defaults)
        return struct

    #===================================================

    def write(self, dest, combine=None, parameters='inline', molfile=None, itp=False):
        """ Write a Gromacs Topology File from a Structure

        Parameters
        ----------
        dest : str or file-like
            The name of a file or a file object to write the Gromacs topology to
        combine : 'all', None, or list of iterables, optional
            If None, no molecules are combined into a single moleculetype. If
            'all', all molecules are combined into a single moleculetype.
            Otherwise, the list of molecule indices (start from 0) will control
            which atoms are combined into single moleculetype's. Default is None
        parameters : 'inline' or str or file-like object, optional
            This specifies where parameters should be printed. If 'inline'
            (default), the parameters are written on the same lines as the
            valence terms are defined on. Any other string is interpreted as a
            filename for an ITP that will be written to and then included at the
            top of `dest`. If it is a file-like object, parameters will be
            written there.  If parameters is the same as ``dest``, then the
            parameter types will be written to the same topologyfile.
        molfile : None or str of file-like object, optional
            If specified as other than None, the molecules will be written to a
            separate file that is included in the main topology file. The
            name of this file will be the provided srting. If None or
            the same as the ``dest'', the molecules will be written into the
            body of the topology file. If it is a file-like object,
            the molecules will be written there. Using this option can make
            it easier to combine multiple molecules into the same topology.
            This will change where the following topology sections are
            written: moleculetype, atoms, bonds, pairs, angles, dihedrals,
            cmap, settles, virtual_sites2, virtual_sites3 and exclusions.
        itp : bool, optional
            If True the following topology sections are not written:
            defaults, atomtypes, nonbond_params, bondtypes, pairtypes,
            angletypes, dihedraltypes, cmaptypes, system and molecules
            Thus only the individual molecules will be written in a stand-alone
            fashion, i.e. an itp-file.
            If True the molfile parameter will be set to None

        Raises
        ------
        ValueError if the same molecule number appears in multiple combine lists
        TypeError if the dest input cannot be parsed
        ValueError if the combine, parameters, or molfile input cannot be parsed
        """
        import parmed.gromacs as gmx
        from .. import __version__
        own_handle = False
        fname = ''
        params = ParameterSet.from_structure(self, allow_unequal_duplicates=True)
        if isinstance(dest, str):
            fname = '%s ' % dest
            dest = genopen(dest, 'w')
            own_handle = True
        elif not hasattr(dest, 'write'):
            raise TypeError('dest must be a file name or file-like object')

        # Determine where to write the parameters
        own_parfile_handle = False
        include_parfile = None
        if parameters == 'inline':
            parfile = dest
        elif isinstance(parameters, str):
            if parameters == fname.strip():
                parfile = dest
            else:
                own_parfile_handle = True
                parfile = genopen(parameters, 'w')
                include_parfile = parameters
        elif hasattr(parameters, 'write'):
            parfile = parameters
        else:
            raise ValueError('parameters must be "inline", a file name, or a file-like object')

        # Determine where to write the molecules
        if itp :
            molfile = None
        own_molfile_handle = False
        include_molfile = None
        if molfile is None:
            _molfile = dest
        elif isinstance(molfile, str):
            if molfile == fname.strip():
                _molfile = dest
            else:
                own_molfile_handle = True
                _molfile = genopen(molfile, 'w')
                include_molfile = molfile
        elif hasattr(molfile, 'write'):
            _molfile = molfile
            include_molfile = _molfile.name
            # I assume the file should still be included even if it's not passed
            # in as a file name. I'm not sure if all `write`-able objects have a
            # `name` property, though.
        else:
            raise ValueError('molfile must be "top", a file name, or a file-like object')

        # Error-checking for combine
        if combine is not None:
            if isinstance(combine, str):
                if combine.lower() != 'all':
                    raise ValueError('combine must be None, list of indices, or "all"')
            else:
                combine_lists = []
                for indices in combine:
                    indices = sorted(set(indices))
                    if any((indices[i+1] - indices[i]) != 1 for i in range(len(indices)-1)):
                        raise ValueError('Can only combine adjacent molecules')
                    combine_lists.append(indices)
        try:
            # Write the header
            now = datetime.now()
            dest.write(f'''\
;
;   File {fname} was generated
;   By user: {_username} ({_userid})
;   On host:{_uname} 
;   At date:{now.strftime('%a. %B  %w %X %Y')} 
;
;   This is a standalone topology file
;
;   Created by:
;   ParmEd:       {os.path.split(sys.argv[0])[1]}, VERSION{__version__} 
;   Executable:   {os.path.split(sys.argv[0])[1]}
;   Library dir:  {gmx.GROMACS_TOPDIR}
;   Command line:
;     {(' '.join(sys.argv)).encode('unicode_escape').decode('utf-8')}
;
''')
            if not itp :
                dest.write('\n[ defaults ]\n')
                dest.write('; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n')
                dest.write(
                    f'{self.defaults.nbfunc:<15d} {self.defaults.comb_rule:<15d} {self.defaults.gen_pairs:<15s} '
                    f'{self.defaults.fudgeLJ:<12.8g} {self.defaults.fudgeQQ:<12.8g}\n\n'
                )
            if include_parfile is not None:
                dest.write(f'#include "{include_parfile}"\n\n')
            # Print all atom types
            if not itp :
                parfile.write('[ atomtypes ]\n')
                if any(typ._bond_type is not None for key, typ in params.atom_types.items()):
                    print_bond_types = True
                else:
                    print_bond_types = False
                if all(typ.atomic_number != -1 for key, typ in params.atom_types.items()):
                    print_atnum = True
                else:
                    print_atnum = False
                parfile.write('; name    ')
                if print_bond_types:
                    parfile.write('bond_type ')
                if print_atnum:
                    parfile.write('at.num    ')
                parfile.write('mass    charge ptype  sigma      epsilon\n')
                econv = u.kilocalories.conversion_factor_to(u.kilojoules)
                for key, atom_type in params.atom_types.items():
                    parfile.write(f'{str(atom_type):<7s} ')
                    if print_bond_types:
                        parfile.write(f'{str(atom_type.bond_type):<8s} ')
                    if print_atnum:
                        parfile.write('%8d ' % atom_type.atomic_number)
                    parfile.write('%10.6f  %10.8f  A %14.8g %14.8g\n' % (
                                  atom_type.mass, atom_type.charge, atom_type.sigma/10,
                                  atom_type.epsilon*econv))
                parfile.write('\n')
            # Nonbonded parameters
            if not itp and self.has_NBFIX():
                typemap = dict(self.parameterset.nbfix_types)
                types_in_system = self.parameterset.atom_types.keys()
                dest.write('[ nonbond_params ]\n')
                eps_conversion = u.kilocalorie.conversion_factor_to(u.kilojoule)
                for key, val in typemap.items():
                    if key[0] in types_in_system and key[1] in types_in_system:
                        eps = val[0] # kcal
                        sig = val[1] # Angstrom
                        eps *= eps_conversion
                        sig *= 0.1
                        dest.write(f'{key[0]} {key[1]} 1 {sig / 2**(1/6)} {eps}\n')
            # Print all parameter types unless we asked for inline
            if not itp and parameters != 'inline':
                if params.bond_types:
                    parfile.write('[ bondtypes ]\n')
                    parfile.write('; i    j  func       b0          kb\n')
                    used_keys = set()
                    conv = (u.kilocalorie/u.angstrom**2).conversion_factor_to(
                                u.kilojoule/u.nanometer**2) * 2
                    for key, param in params.bond_types.items():
                        if key in used_keys: continue
                        used_keys.add(key)
                        used_keys.add(tuple(reversed(key)))
                        parfile.write(f'{key[0]:<5s} {key[1]:<5s}    1   {param.req / 10:.5f}   {param.k * conv:f}\n')
                    parfile.write('\n')
                if params.pair_types and self.defaults.gen_pairs == 'no':
                    parfile.write('[ pairtypes ]\n')
                    parfile.write('; i j   func    sigma1-4    epsilon1-4 ; THESE ARE 1-4 INTERACTIONS\n')
                    econv = u.kilocalorie.conversion_factor_to(u.kilojoule)
                    lconv = u.angstrom.conversion_factor_to(u.nanometer)
                    used_keys = set()
                    for key, param in params.pair_types.items():
                        if key in used_keys:
                            continue
                        used_keys.add(key)
                        used_keys.add(tuple(reversed(key)))
                        parfile.write(
                            f'{key[0]:<5s} {key[1]:<5s}  1  {param.sigma * lconv:.9f} {param.epsilon * econv:.9f}\n'
                        )
                    parfile.write('\n')
                if params.angle_types:
                    parfile.write('[ angletypes ]\n')
                    parfile.write(';  i    j    k  func       th0       cth    rub         kub\n')
                    used_keys = set()
                    conv = (u.kilocalorie/u.radian**2).conversion_factor_to(
                                u.kilojoule/u.radian**2) * 2
                    bconv = (u.kilocalorie/u.angstrom**2).conversion_factor_to(
                                u.kilojoule/u.nanometer**2) * 2
                    for key, param in params.angle_types.items():
                        if key in used_keys:
                            continue
                        used_keys.add(key)
                        used_keys.add(tuple(reversed(key)))
                        part = '%-5s %-5s %-5s    %%d   %12.7f   %12.7f' % (
                                key[0], key[1], key[2], param.theteq, param.k*conv)
                        if key in params.urey_bradley_types:
                            ub = params.urey_bradley_types[key]
                            parfile.write(part % 5)
                            parfile.write('  %12.7f  %12.7f\n' % (0.1*ub.req, ub.k*bconv))
                        else:
                            parfile.write(part % 1)
                            parfile.write('\n')
                    parfile.write('\n')
                if params.dihedral_types:
                    parfile.write('[ dihedraltypes ]\n')
                    parfile.write(';i  j   k  l  func      phase      kd      pn\n')
                    used_keys = set()
                    conv = u.kilocalories.conversion_factor_to(u.kilojoules)
                    fmt = '%-6s %-6s %-6s %-6s  %d   %.2f   %.6f   %d\n'
                    for key, param in params.dihedral_types.items():
                        if key in used_keys: continue
                        used_keys.add(key)
                        used_keys.add(tuple(reversed(key)))
                        for dt in param:
                            parfile.write(fmt % (key[0], key[1], key[2], key[3], 9, dt.phase,
                                          dt.phi_k*conv, int(dt.per)))
                    parfile.write('\n')
                if params.improper_periodic_types:
                    parfile.write('[ dihedraltypes ]\n')
                    parfile.write(';i  j   k  l  func      phase      kd      pn\n')
                    used_keys = set()
                    conv = u.kilojoules.conversion_factor_to(u.kilocalories)
                    fmt = '%-6s %-6s %-6s %-6s  %d   %.2f   %.6f   %d\n'
                    for key, param in params.improper_periodic_types.items():
                        if key in used_keys: continue
                        used_keys.add(key)
                        used_keys.add(tuple(reversed(key)))
                        parfile.write(fmt % (key[0], key[1], key[2], key[3],
                                      4, param.phase, param.phi_k*conv, int(param.per)))
                    parfile.write('\n')
                if params.improper_types:
                    # BUGBUG -- The ordering is borked here because that made it
                    # simpler for me to work with back when I wrote the CHARMM
                    # parsers. This needs to be fixed now and handled correctly.
                    parfile.write('[ dihedraltypes ]\n')
                    parfile.write('; i  j       k       l       func     q0    cq\n')
                    fmt = '%-6s %-6s %-6s %-6s    %d   %.6f   %.6f\n'
                    conv = u.kilocalories.conversion_factor_to(u.kilojoules)*2
                    for key, param in params.improper_types.items():
                        parfile.write(fmt % (key[0], key[1], key[2], key[3],
                                      2, param.psi_eq, param.psi_k*conv))
                    parfile.write('\n')
            # CMAP grids are never printed inline, so if we have them, we need
            # to write a dedicated section for them
            if not itp and params.cmap_types:
                    parfile.write('[ cmaptypes ]\n\n')
                    used_keys = set()
                    conv = u.kilocalories.conversion_factor_to(u.kilojoules)
                    for key, param in params.cmap_types.items():
                        if key in used_keys: continue
                        used_keys.add(key)
                        used_keys.add(tuple(reversed(key)))
                        parfile.write(
                            f'{key[0]:<6s} {key[1]:<6s} {key[2]:<6s} {key[3]:<6s} {key[7]:<6s}   '
                            f'1   {param.resolution:4d} {param.resolution:4d}'
                        )
                        res2 = param.resolution * param.resolution
                        for i in range(0, res2, 10):
                            parfile.write('\\\n')
                            end = min(i+10, res2)
                            parfile.write(' '.join(str(param.grid[j]*conv) for j in range(i, end)))
                        parfile.write('\n\n')
            if include_molfile is not None:
                dest.write(f'#include "{include_molfile}"\n\n')
            if combine is None:
                molecules = self.split()
                sysnum = 1
                names = []
                nameset = set()
                for molecule, num in molecules:
                    if len(molecule.residues) == 1:
                        title = molecule.residues[0].name if molecule.residues[0].name.strip() else "UNK"
                        if title in nameset:
                            orig = title
                            sfx = 2
                            while title in nameset:
                                title = f'{orig}{sfx}'
                                sfx += 1
                    else:
                        title = f'system{sysnum}'
                        sysnum += 1
                    names.append(title)
                    nameset.add(title)
                    GromacsTopologyFile._write_molecule(molecule, _molfile, title, params, parameters == 'inline')
                if not itp :
                    # System
                    dest.write('[ system ]\n; Name\n')
                    if self.title:
                        dest.write(self.title)
                    else:
                        dest.write('Generic title')
                    dest.write('\n\n')
                    # Molecules
                    dest.write('[ molecules ]\n; Compound       #mols\n')
                    total_mols = sum(len(m[1]) for m in molecules)
                    i = 0
                    while i < total_mols:
                        for j, (molecule, lst) in enumerate(molecules):
                            if i in lst:
                                break
                        else:
                            raise AssertionError(f'Could not find molecule {i:d} in list')
                        ii = i
                        while ii < total_mols and ii in lst:
                            ii += 1
                        dest.write(f'{names[j]:<15s} {ii - i:6d}\n')
                        i = ii
            elif isinstance(combine, str) and combine.lower() == 'all':
                GromacsTopologyFile._write_molecule(self, _molfile, 'system', params, parameters == 'inline')
                if not itp :
                    dest.write('[ system ]\n; Name\n')
                    if self.title:
                        dest.write(self.title)
                    else:
                        dest.write('Generic title') # pragma: no cover
                    dest.write('\n\n')
                    # Molecules
                    dest.write('[ molecules ]\n; Compound       #mols\n')
                    dest.write('%-15s %6d\n' % ('system', 1))
            else:
                molecules = self.split()
                nmols = sum(len(m[1]) for m in molecules)
                moleculedict = dict()
                # Hash our molecules by indices
                for m, num in molecules:
                    for i in num:
                        moleculedict[i] = m
                combined_molecules = []
                for cl in combine_lists:
                    counts = defaultdict(int)
                    mols_in_mol = []
                    for molid in cl:
                        try:
                            mol = moleculedict[molid]
                        except KeyError:
                            raise IndexError('Molecule ID out of range')
                        counts[id(moleculedict[molid])] += 1
                        if counts[id(moleculedict[molid])] == 1:
                            mols_in_mol.append(mol)
                    if counts[id(mols_in_mol[0])] > 1:
                        combmol = mols_in_mol[0] * counts[id(mols_in_mol[0])]
                    else:
                        combmol = copy.copy(mols_in_mol[0])
                    for i, mol in enumerate(mols_in_mol):
                        if i == 0: continue
                        assert id(mol) in counts and counts[id(mol)] > 0
                        if counts[id(mol)] > 1:
                            combmol += mol * counts[id(mol)]
                        else:
                            combmol += mol
                    combined_molecules.append((combmol, cl[0], len(cl)))
                    nmols -= (len(cl) - 1)
                # combined_molecules now contains a list of tuples, and that
                # tuple stores the combined molecule, first molecule index of
                # the pre-combined molecule, and how many molecules were
                # combined

                # Sort combined molecules by starting location
                combined_molecules.sort(key=lambda x: x[1])
                new_molecules = []
                counts = defaultdict(set)
                cmc = 0 # Combined Molecule Counter
                add = 0 # How many molecules to "skip" due to combining
                for i in range(nmols):
                    ii = i + add
                    if cmc < len(combined_molecules) and combined_molecules[cmc][1] == ii:
                        new_molecules.append([combined_molecules[cmc][0], set([i])])
                        add += combined_molecules[cmc][2] - 1
                        cmc += 1
                    elif len(counts[id(moleculedict[ii])]) == 0:
                        counts[id(moleculedict[ii])].add(i)
                        new_molecules.append([moleculedict[ii], counts[id(moleculedict[ii])]])
                    else:
                        counts[id(moleculedict[ii])].add(i)
                sysnum = 1
                names = []
                nameset = set()
                for molecule, num in new_molecules:
                    if len(molecule.residues) == 1:
                        title = molecule.residues[0].name if molecule.residues[0].name.strip() else "UNK"
                        if title in nameset:
                            orig = title
                            sfx = 2
                            while title in nameset:
                                title = f'{orig:s}{sfx:d}'
                                sfx += 1
                    else:
                        title = f'system{sysnum}'
                        sysnum += 1
                    names.append(title)
                    nameset.add(title)
                    GromacsTopologyFile._write_molecule(molecule, _molfile, title, params, parameters == 'inline')
                if not itp :
                    # System
                    dest.write('[ system ]\n; Name\n')
                    if self.title:
                        dest.write(self.title)
                    else:
                        dest.write('Generic title') # pragma: no cover
                    dest.write('\n\n')
                    # Molecules
                    dest.write('[ molecules ]\n; Compound       #mols\n')
                    total_mols = sum(len(m[1]) for m in new_molecules)
                    i = 0
                    while i < total_mols:
                        for j, (molecule, lst) in enumerate(new_molecules):
                            if i in lst:
                                break
                        else:
                            raise AssertionError(f'Could not find molecule {i} in list')
                        ii = i
                        while ii < total_mols and ii in lst:
                            ii += 1
                        dest.write(f'{names[j]:<15s} {ii - i:6d}\n')
                        i = ii
        finally:
            if own_handle:
                dest.close()
            if own_parfile_handle:
                parfile.close()
            if own_molfile_handle:
                _molfile.close()

    #===================================================

    @staticmethod
    def _write_molecule(struct, dest, title, params, writeparams):
        assert title.strip(), "title must not be whitespace only!"
        dest.write('\n[ moleculetype ]\n; Name            nrexcl\n')
        dest.write(f'{title}          {struct.nrexcl}\n\n')
        dest.write('[ atoms ]\n')
        dest.write(';   nr       type  resnr residue  atom   cgnr    '
                   'charge       mass  typeB    chargeB      massB\n')
        runchg = 0
        for residue in struct.residues:
            dest.write('; residue %4d %s rtp %s q %.1f\n' %
                       (residue.idx+1, residue.name, residue.name, sum(a.charge for a in residue)))
            for atom in residue:
                runchg += atom.charge
                dest.write('%5d %10s %6d %6s %6s %6d %10.8f %10.6f   ; '
                           'qtot %.6f\n' % (atom.idx+1, atom.type,
                            residue.idx+1, residue.name, atom.name,
                            atom.idx+1, atom.charge, atom.mass, runchg))
        dest.write('\n')
        # Do valence terms now
        EPs = [a for a in struct.atoms if isinstance(a, ExtraPoint)]
        settle = False
        if len(struct.atoms) - len(EPs) == 3:
            try:
                oxy, = (a for a in struct.atoms if a.atomic_number == 8)
                hyd1, hyd2 = (a for a in struct.atoms if a.atomic_number == 1)
                settle = True
            except ValueError:
                pass
        if struct.bonds:
            conv = (u.kilocalorie_per_mole/u.angstrom**2).conversion_factor_to(
                    u.kilojoule_per_mole/u.nanometer**2)*2
            if settle:
                dest.write('#ifdef FLEXIBLE\n\n')
            dest.write('[ bonds ]\n')
            dest.write(';%6s %6s %5s %10s %10s %10s %10s\n' % ('ai', 'aj',
                       'funct', 'c0', 'c1', 'c2', 'c3'))
            for bond in struct.bonds:
                if (isinstance(bond.atom1, ExtraPoint) or isinstance(bond.atom2, ExtraPoint)):
                    continue # pragma: no cover
                dest.write('%7d %6d %5d' % (bond.atom1.idx+1, bond.atom2.idx+1, bond.funct))
                if bond.type is None:
                    dest.write('\n')
                    continue # pragma: no cover
                key = (_gettype(bond.atom1), _gettype(bond.atom2))
                if writeparams or key not in params.bond_types or \
                        bond.type != params.bond_types[key]:
                    dest.write('   %.5f %f' % (bond.type.req/10, bond.type.k*conv))
                dest.write('\n')
            dest.write('\n')
        # Do the pair-exceptions
        if struct.adjusts:
            dest.write('[ pairs ]\n')
            dest.write(';%6s %6s %5s %10s %10s %10s %10s\n' % ('ai', 'aj',
                       'funct', 'c0', 'c1', 'c2', 'c3'))
            econv = u.kilocalories.conversion_factor_to(u.kilojoules)
            lconv = u.angstroms.conversion_factor_to(u.nanometer)
            for adjust in struct.adjusts:
                key = (_gettype(adjust.atom1), _gettype(adjust.atom2))
                dest.write('%7d %6d %5d' % (adjust.atom1.idx+1, adjust.atom2.idx+1, adjust.funct))
                if struct.defaults.gen_pairs == 'no' and (writeparams or
                        key not in params.pair_types or
                        adjust.type != params.pair_types[key]) and adjust.type is not None:
                    dest.write(' %.9f %.9f' % (adjust.type.sigma*lconv, adjust.type.epsilon*econv))
                dest.write('\n')
            dest.write('\n')
        elif struct.dihedrals:
            dest.write('[ pairs ]\n')
            dest.write(';%6s %6s %5s %10s %10s %10s %10s\n' % ('ai', 'aj',
                       'funct', 'c0', 'c1', 'c2', 'c3'))
            # Get the 1-4 pairs from the dihedral list
            struct.update_dihedral_exclusions()
            econv = u.kilocalories.conversion_factor_to(u.kilojoules)
            lconv = u.angstroms.conversion_factor_to(u.nanometer)
            for dihed in struct.dihedrals:
                if dihed.ignore_end or dihed.improper: continue
                a1, a2 = dihed.atom1, dihed.atom4
                if a1 in a2.bond_partners or a1 in a2.angle_partners:
                    continue # pragma: no cover
                dest.write('%7d %6d %5d' % (a1.idx+1, a2.idx+1, 1))
                if struct.defaults.gen_pairs == 'no':
                    dest.write('  %.9f  %.9f' % (0.5*(a1.sigma_14+a2.sigma_14)*lconv,
                                math.sqrt(a1.epsilon_14*a2.epsilon_14)*econv))
                dest.write('\n')
            dest.write('\n')
        # Angles
        if struct.angles:
            conv = (u.kilocalorie_per_mole/u.radian**2).conversion_factor_to(
                        u.kilojoule_per_mole/u.radian**2)*2
            conv2 = (u.kilocalorie_per_mole/u.angstrom**2).conversion_factor_to(
                    u.kilojoule_per_mole/u.nanometer**2)*2
            dest.write('[ angles ]\n')
            dest.write(';%6s %6s %6s %5s %10s %10s %10s %10s\n' %
                       ('ai', 'aj', 'ak', 'funct', 'c0', 'c1', 'c2', 'c3'))
            for angle in struct.angles:
                dest.write('%7d %6d %6d %5d' % (angle.atom1.idx+1, angle.atom2.idx+1,
                           angle.atom3.idx+1, angle.funct))
                if angle.type is None:
                    dest.write('\n')
                    continue
                key = (_gettype(angle.atom1), _gettype(angle.atom2), _gettype(angle.atom3))
                param_equal = params.angle_types.get(key) == angle.type
                if angle.funct == 5:
                    # Find the Urey-Bradley term, if it exists
                    for ub in struct.urey_bradleys:
                        if angle.atom1 in ub and angle.atom3 in ub:
                            ubtype = ub.type
                            break
                    else:
                        ubtype = NoUreyBradley
                    param_equal = param_equal and params.urey_bradley_types.get(key) == ubtype
                if writeparams or not param_equal:
                    dest.write('   %.7f %f' % (angle.type.theteq, angle.type.k*conv))
                    if angle.funct == 5:
                        dest.write(' %.7f %f' % (ubtype.req/10, ubtype.k*conv2))
                dest.write('\n')
            dest.write('\n')
        # Dihedrals
        if struct.dihedrals:
            dest.write('[ dihedrals ]\n')
            dest.write((';%6s %6s %6s %6s %5s'+' %10s'*6) % ('ai', 'aj',
                       'ak', 'al', 'funct', 'c0', 'c1', 'c2', 'c3', 'c4', 'c5'))
            dest.write('\n')
            conv = u.kilocalories.conversion_factor_to(u.kilojoules)
            for dihed in struct.dihedrals:
                dest.write('%7d %6d %6d %6d %5d' % (dihed.atom1.idx+1, dihed.atom2.idx+1,
                           dihed.atom3.idx+1, dihed.atom4.idx+1, dihed.funct))
                if dihed.type is None:
                    dest.write('\n')
                    continue
                if dihed.improper:
                    typedict = params.improper_periodic_types
                else:
                    typedict = params.dihedral_types
                key = (_gettype(dihed.atom1), _gettype(dihed.atom2),
                        _gettype(dihed.atom3), _gettype(dihed.atom4))
                if writeparams or key not in typedict or \
                        _diff_diheds(dihed.type, typedict[key]):
                    if isinstance(dihed.type, DihedralTypeList):
                        dest.write('  %.6f  %.6f  %d' % (dihed.type[0].phase,
                            dihed.type[0].phi_k*conv, int(dihed.type[0].per)))
                        for dt in dihed.type[1:]:
                            dest.write('\n%7d %6d %6d %6d %5d  %.5f  %.7f  %d' %
                                    (dihed.atom1.idx+1, dihed.atom2.idx+1,
                                     dihed.atom3.idx+1, dihed.atom4.idx+1,
                                     dihed.funct, dt.phase, dt.phi_k*conv, int(dt.per)))
                    else:
                        dest.write('  %.7f  %.7f  %d' % (dihed.type.phase,
                            dihed.type.phi_k*conv, int(dihed.type.per)))
                dest.write('\n')
            dest.write('\n')
        # RB-torsions
        if struct.rb_torsions:
            dest.write('[ dihedrals ]\n')
            dest.write((';%6s %6s %6s %6s %5s'+' %10s'*6) % ('ai', 'aj',
                       'ak', 'al', 'funct', 'c0', 'c1', 'c2', 'c3',
                       'c4', 'c5'))
            dest.write('\n')
            conv = u.kilocalories.conversion_factor_to(u.kilojoules)
            paramfmt = '  %12.7f  %12.7f  %12.7f  %12.7f  %12.7f  %12.7f'
            for dihed in struct.rb_torsions:
                dest.write('%7d %6d %6d %6d %5d' % (dihed.atom1.idx+1,
                           dihed.atom2.idx+1, dihed.atom3.idx+1,
                           dihed.atom4.idx+1, dihed.funct))
                if dihed.type is None:
                    dest.write('\n')
                    continue
                key = (_gettype(dihed.atom1), _gettype(dihed.atom2),
                        _gettype(dihed.atom3), _gettype(dihed.atom4))
                if writeparams or key not in params.rb_torsion_types or \
                        params.rb_torsion_types[key] != dihed.type:
                    dest.write(paramfmt % (
                        dihed.type.c0*conv, dihed.type.c1*conv, dihed.type.c2*conv,
                        dihed.type.c3*conv, dihed.type.c4*conv, dihed.type.c5*conv
                    ))
                    dest.write('\n')
            dest.write('\n')
        # Impropers
        if struct.impropers:
            dest.write('[ dihedrals ]\n')
            dest.write((';%6s %6s %6s %6s %5s'+' %10s'*4) % ('ai', 'aj',
                       'ak', 'al', 'funct', 'c0', 'c1', 'c2', 'c3'))
            dest.write('\n')
            conv = u.kilocalories.conversion_factor_to(u.kilojoules) * 2
            for dihed in struct.impropers:
                dest.write('%7d %6d %6d %6d %5d' % (dihed.atom1.idx+1,
                           dihed.atom2.idx+1, dihed.atom3.idx+1,
                           dihed.atom4.idx+1, dihed.funct))
                if dihed.type is None:
                    dest.write('\n')
                    continue
                # BUGBUG: We always write improper types since we don't
                # currently store the correct ordering of the types in the
                # improper section
                dest.write('  %12.7f  %12.7f\n' % (dihed.type.psi_eq, dihed.type.psi_k*conv))
            dest.write('\n')
        # Cmaps
        if struct.cmaps:
            dest.write('[ cmap ]\n')
            dest.write(';%6s %6s %6s %6s %6s %5s\n' % ('ai', 'aj', 'ak',
                       'al', 'am', 'funct'))
            for cmap in struct.cmaps:
                dest.write('%7d %6d %6d %6d %6d %5d\n' % (cmap.atom1.idx+1,
                           cmap.atom2.idx+1, cmap.atom3.idx+1,
                           cmap.atom4.idx+1, cmap.atom5.idx+1, cmap.funct))
        # See if this is a solvent molecule with 3 or fewer particles that can
        # be SETTLEd
        if settle:
            dest.write('\n#else\n\n')
            dest.write('[ settles ]\n')
            dest.write('; i     funct   doh     dhh\n')
            for b in oxy.bonds:
                if hyd1 in b:
                    # Use default values for TIPnP if no type exists
                    doh = 0.09572 if b.type is None else b.type.req / 10
                    break
            for b in hyd1.bonds:
                if hyd2 in b:
                    # Use default values for TIPnP if no type exists
                    dhh = 0.15139 if b.type is None else b.type.req / 10
                    break
            else:
                for a in oxy.angles:
                    if hyd1 in a and hyd2 in a:
                        theteq = a.type.theteq * DEG_TO_RAD
                        dhh = math.sqrt(2*doh*doh - 2*doh*doh*math.cos(theteq))
                        break
                else:
                    raise GromacsError('Cannot determine SETTLE geometry') # pragma: no cover
            dest.write('1     1   %.8f   %.8f\n\n#endif\n\n' % (doh, dhh))
        # Virtual sites
        if EPs:
            ftypes = set(type(a.frame_type) for a in EPs)
            for ftype in ftypes:
                if ftype is TwoParticleExtraPointFrame:
                    dest.write('[ virtual_sites2 ]\n')
                    dest.write('; Site  from    funct  a\n')
                    for EP in EPs:
                        if not isinstance(EP.frame_type, ftype): continue
                        a1, a2 = EP.frame_type.get_atoms()
                        dest.write('%-5d %-4d %-4d %-4d   %.6f\n' %
                                   (EP.idx+1, a1.idx+1, a2.idx+1, 1,
                                    EP.frame_type.get_weights()[0]))
                    dest.write('\n')
                elif ftype in (ThreeParticleExtraPointFrame,
                               OutOfPlaneExtraPointFrame):
                    dest.write('[ virtual_sites3 ]\n')
                    dest.write('; Site  from                   funct\n')
                    for EP in EPs:
                        if isinstance(EP.frame_type,
                                ThreeParticleExtraPointFrame):
                            a1, a2, a3 = EP.frame_type.get_atoms()
                            junk, w1, w2 = EP.frame_type.get_weights()
                            dest.write('%-5d %-4d %-4d %-4d %-4d   %.6f  %.6f\n'
                                       % (EP.idx+1, a1.idx+1, a2.idx+1, a3.idx+1, 1, w1, w2))
                        elif isinstance(EP.frame_type,
                                OutOfPlaneExtraPointFrame):
                            a1, a2, a3 = EP.frame_type.get_atoms()
                            w1, w2, w3 = EP.frame_type.get_weights()
                            dest.write('%-5d %-4d %-4d %-4d %-4d   %.6f  %.6f  %.6f\n' %
                                      (EP.idx+1, a1.idx+1, a2.idx+1, a3.idx+1, w1, w2, w3))
                    dest.write('\n')
        # Do we need to list exclusions for systems with EPs?
        if EPs or settle:
            dest.write('[ exclusions ]\n')
            for i, atom in enumerate(struct.atoms):
                dest.write('%d' % (i+1))
                for a in atom.bond_partners:
                    dest.write('  %d' % (a.idx+1))
                for a in atom.angle_partners:
                    dest.write('  %d' % (a.idx+1))
                dest.write('\n')
            dest.write('\n')

    #===================================================

    def __getstate__(self):
        d = Structure.__getstate__(self)
        d['parameterset'] = self.parameterset
        d['defaults'] = self.defaults
        return d

    def __setstate__(self, d):
        Structure.__setstate__(self, d)
        self.parameterset = d['parameterset']
        self.defaults = d['defaults']

def _any_atoms_farther_than(structure, limit=3):
    """
    This function checks to see if there are any atom pairs farther away in the
    bond graph than the desired limit

    Parameters
    ----------
    structure : :class:`Structure`
        The structure to search through
    limit : int, optional
        The most number of bonds away to check for. Default is 3

    Returns
    -------
    within : bool
        True if any atoms are *more* than ``limit`` bonds away from any other
        atom
    """
    import sys
    if len(structure.atoms) <= limit + 1: return False
    sys.setrecursionlimit(max(sys.getrecursionlimit(), limit+1))
    for atom in structure.atoms:
        for atom in structure.atoms: atom.marked = limit + 1
        _mark_graph(atom, 0)
        if any((atom.marked > limit for atom in structure.atoms)):
            return True
    return False

def _mark_graph(atom, num):
    """ Marks all atoms in the graph listing the minimum number of bonds each
    atom is away from the current atom

    Parameters
    ----------
    atom : :class:`Atom`
        The current atom to evaluate in the bond graph
    num : int
        The current depth in our search
    limit : int
        The maximum depth we want to search
    """
    atom.marked = num
    for a in atom.bond_partners:
        if a.marked <= num: continue
        _mark_graph(a, num+1)

def _diff_diheds(dt1, dt2):
    """ Determine if 2 dihedrals are *really* different. dt1 can either be a
    DihedralType or a DihedralTypeList or dt1 can be a DihedralType and dt2 can
    be a DihedralTypeList.  This returns True if dt1 == dt2 *or* dt1 is equal to
    the only element of dt2
    """
    if type(dt1) is type(dt2) and dt1 == dt2:
        return False
    if isinstance(dt2, DihedralTypeList) and isinstance(dt1, DihedralType):
        if len(dt2) == 1 and dt2[0] == dt1: return False
    return True

def _gettype(atom):
    if atom.atom_type not in (None, UnassignedAtomType):
        return atom.atom_type.bond_type
    return atom.type
