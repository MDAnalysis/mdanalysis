"""
List of all pointers and constants used in the Amber topology file.

Can be used like:
   from parmed.constants import *
"""
from enum import IntEnum

from math import pi as _pi, sqrt as _sqrt, log10 as _log10, acos as _acos

__all__ = ['AMBER_ELECTROSTATIC', 'AMBER_POINTERS', 'PrmtopPointers',
           'RAD_TO_DEG', 'DEG_TO_RAD', 'TRUNCATED_OCTAHEDRON_ANGLE']

AMBER_ELECTROSTATIC = 18.2223
CHARMM_ELECTROSTATIC = _sqrt(332.0716)

AMBER_POINTERS = """
NATOM  : total number of atoms
NTYPES : total number of distinct atom types
NBONH  : number of bonds containing hydrogen
MBONA  : number of bonds not containing hydrogen
NTHETH : number of angles containing hydrogen
MTHETA : number of angles not containing hydrogen
NPHIH  : number of dihedrals containing hydrogen
MPHIA  : number of dihedrals not containing hydrogen
NHPARM : currently not used
NPARM  : currently not used
NEXT   : number of excluded atoms
NRES   : number of residues
NBONA  : MBONA + number of constraint bonds
NTHETA : MTHETA + number of constraint angles
NPHIA  : MPHIA + number of constraint dihedrals
NUMBND : number of unique bond types
NUMANG : number of unique angle types
NPTRA  : number of unique dihedral types
NATYP  : number of atom types in parameter file, see SOLTY below
NPHB   : number of distinct 10-12 hydrogen bond pair types
IFPERT : set to 1 if perturbation info is to be read in
NBPER  : number of bonds to be perturbed
NGPER  : number of angles to be perturbed
NDPER  : number of dihedrals to be perturbed
MBPER  : number of bonds with atoms completely in perturbed group
MGPER  : number of angles with atoms completely in perturbed group
MDPER  : number of dihedrals with atoms completely in perturbed groups
IFBOX  : set to 1 if standard periodic box, 2 when truncated octahedral
NMXRS  : number of atoms in the largest residue
IFCAP  : set to 1 if the CAP option from edit was specified
NUMEXTRA: number of extra points
NCOPY  : Number of copies for advanded simulations
"""
# These global variables provide a more natural way of accessing
# the various pointers.  Most useful if they're loaded into the
# top-level namespace.
class PrmtopPointers(IntEnum):
    NATOM = 0
    NTYPES = 1
    NBONH = 2
    MBONA = 3
    NTHETH = 4
    MTHETA = 5
    NPHIH = 6
    MPHIA = 7
    NHPARM = 8
    NPARM = 9
    NEXT = 10
    NRES = 11
    NBONA = 12
    NTHETA = 13
    NPHIA = 14
    NUMBND = 15
    NUMANG = 16
    NPTRA = 17
    NATYP = 18
    NPHB = 19
    IFPERT = 20
    NBPER = 21
    NGPER = 22
    NDPER = 23
    MBPER = 24
    MGPER = 25
    MDPER = 26
    IFBOX = 27
    NMXRS = 28
    IFCAP = 29
    NUMEXTRA = 30
    NCOPY = 31
    NNB = NEXT

RAD_TO_DEG = 180.0 / _pi
DEG_TO_RAD = _pi / 180.0
TRUNCATED_OCTAHEDRON_ANGLE = _acos(-1/3) * 180 / _pi

# For use in floating point comparisons
TINY = 1.0e-8
SMALL = 1.0e-4
TINY_DIGITS = int(_log10(TINY) + 0.5)
SMALL_DIGITS = int(_log10(SMALL) + 0.5)

# For I/O
DEFAULT_ENCODING = 'UTF-8'
