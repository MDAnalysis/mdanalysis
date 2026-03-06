"""
The parmed package manipulates molecular structures and provides a way to go
between standard and amber file formats, manipulate structures, etc.
"""

# Version format should be "major.minor.patch". For beta releases, attach
# "-beta#" to the end. The beta number will be turned into another number in the
# version tuple
__author__ = 'Jason Swails'

__all__ = ['exceptions', 'periodic_table', 'residue', 'unit', 'utils', 'Structure', 'entos', 'dlpoly',
           'StructureView', 'amber', 'charmm', 'namd', 'gromacs', 'tinker', 'openmm', 'rosetta',
           'rdkit', 'formats', 'Vec3', 'ParameterSet', 'load_file', 'read_PDB', 'read_CIF',
           'load_rosetta', 'load_rdkit', 'download_PDB', 'download_CIF', 'tools', 'version']

from . import _version
__version__ = _version.get_versions()['version']

from . import exceptions, periodic_table, residue
from . import unit, utils
from .topologyobjects import *
from .structure import Structure, StructureView
from . import amber, charmm, gromacs, dlpoly, namd, openmm, rosetta, tinker, entos
from . import formats
from .vec3 import Vec3
from .parameters import ParameterSet
from .rdkit import load_rdkit
from . import rdkit
load_file = formats.load_file
read_PDB = formats.PDBFile.parse
read_CIF = formats.CIFFile.parse
read_pdb = read_PDB
read_cif = read_CIF
load_rosetta = rosetta.RosettaPose.load

download_PDB = formats.PDBFile.download
download_CIF = formats.CIFFile.download
download_pdb = download_PDB
download_cif = download_CIF

# The tools package depends on *everything*, so import this at the end to avoid
# circular imports.
from . import tools  # isort: skip

# Add all of the objects from parmed.topologyobjects to the top-level namespace
__all__ += topologyobjects.__all__

# Build a version tuple from __version__ for easy comparison
import re
from collections import namedtuple
class version(namedtuple('version', ['major', 'minor', 'patchlevel', 'commits_ahead'])):
    def __eq__(self, other):
        try:
            if len(other) == 3:
                return (other == self[:3] and self.commits_ahead == 0 and not self.dirty)
        except TypeError:
            return NotImplemented
        if tuple(other) != tuple(self):
            return False
        elif not hasattr(other, 'git_hash') or not hasattr(other, 'dirty'):
            return not self.dirty
        return self.git_hash == other.git_hash and self.dirty is other.dirty

    def __ne__(self, other):
        return not self == other

    def __le__(self, other):
        return self == other or self < other

    def __ge__(self, other):
        return self == other or self > other

# Build the version
_versionre = re.compile(r'(\d+)\.(\d+)\.(\d+)\+?(\d*)\.?g?([\dabcdefABCDEF]*)\.*(dirty)?')
if _versionre.match(__version__):
    versionlist = list(_versionre.match(__version__).groups())
    versionlist[3] = versionlist[3] or 0
    version = version(*[int(v) for v in versionlist[:4]])
    version.git_hash = versionlist[4]
    version.dirty = bool(versionlist[5])
else:
    versionlist = None
    version = version(0, 0, 0, 0)

# Clean up
del namedtuple, re, versionlist
