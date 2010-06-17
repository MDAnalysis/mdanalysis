"""
Location of data files for the MDAnalysis unit tests
====================================================

Real MD simulation data are stored in the ``data/`` sub
directory. File names are in all caps. Use as ::

  from MDAnalysis.tests.datafiles import *

"""
__all__ = ["PSF", "DCD",                 # Charmm
           "PDB_small",                  # pdb
           "PDB", "XTC", "TRR", "TPR"]   # Gromacs

from pkg_resources import resource_filename

PSF = resource_filename(__name__, 'data/adk.psf')
DCD = resource_filename(__name__, 'data/adk_dims.dcd')
DCD_empty = resource_filename(__name__, 'data/empty.dcd')

PDB_small = resource_filename(__name__, 'data/adk_open.pdb')

PDB = resource_filename(__name__, 'data/adk_oplsaa.pdb')
XTC = resource_filename(__name__, 'data/adk_oplsaa.xtc')
TRR = resource_filename(__name__, 'data/adk_oplsaa.trr')
TPR = resource_filename(__name__, 'data/adk_oplsaa.tpr')
