"""
A package dealing with different file formats and automatic detection of those
formats
"""

__all__ = ['load_file', 'PDBFile', 'CIFFile', 'Mol2File', 'PSFFile', 'PQRFile', 'SDFFile']

from .registry import load_file
from .mol2 import Mol2File
from .pdb import PDBFile, CIFFile
from .pqr import PQRFile
from .psf import PSFFile
from .sdf import SDFFile

# Now let's modify structure.Structure and add our write methods from our
# various formats
from ..structure import Structure
Structure.write_pdb = PDBFile.write
Structure.write_cif = CIFFile.write
Structure.write_psf = PSFFile.write
