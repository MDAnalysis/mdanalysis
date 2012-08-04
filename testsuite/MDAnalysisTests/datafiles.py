# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://mdanalysis.googlecode.com
# Copyright (c) 2006-2011 Naveen Michaud-Agrawal,
#               Elizabeth J. Denning, Oliver Beckstein,
#               and contributors (see website for details)
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
#     N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and
#     O. Beckstein. MDAnalysis: A Toolkit for the Analysis of
#     Molecular Dynamics Simulations. J. Comput. Chem. 32 (2011), 2319--2327,
#     doi:10.1002/jcc.21787
#

"""
Location of data files for the MDAnalysis unit tests
====================================================

Real MD simulation data are stored in the ``data/`` sub
directory. Use as ::

  from MDAnalysis.tests.datafiles import *

Note that the files are actually located in a separate package,
:mod:`MDAnalysisTestData` from where they are initially imported as ::

 from MDAnalysisTestData.datafiles import *
"""
__all__ = ["PSF", "DCD", "CRD",               # CHARMM
           "DCD_empty",
           "PSF_NAMD", "PDB_NAMD",            # NAMD
           "PDB_small","NUCL",                # PDB
           "PDB_closed",
           "PDB_multiframe",
           "PDB_helix",
           "PDB", "GRO", "XTC", "TRR", "TPR", "GRO_velocity",   # Gromacs
           "XYZ", "XYZ_psf", "XYZ_bz2",         # XYZ
           "PRM", "TRJ", "TRJ_bz2",             # Amber (no periodic box)
           "PRMpbc", "TRJpbc_bz2",              # Amber (periodic box)
           "PRM12", "TRJ12_bz2",                # Amber (v12 format, Issue 100)
           "PRMncdf", "TRJncdf", "NCDF",        # Amber (netcdf)
           "PQR",                               # PQR
           "PDBQT_input",                       # PDBQT
           "PDBQT_querypdb",
           "FASTA",                             # sequence alignment, Issue 112 + 113
           "PDB_HOLE",                          # gramicidin A
           ]

from pkg_resources import resource_filename

PSF = resource_filename(__name__, 'data/adk.psf')
DCD = resource_filename(__name__, 'data/adk_dims.dcd')
DCD_empty = resource_filename(__name__, 'data/empty.dcd')
CRD = resource_filename(__name__, 'data/adk_open.crd')

PSF_NAMD = resource_filename(__name__, 'data/namd_cgenff.psf')
PDB_NAMD = resource_filename(__name__, 'data/namd_cgenff.pdb')

PDB_small = resource_filename(__name__, 'data/adk_open.pdb')
PDB_closed = resource_filename(__name__, 'data/adk_closed.pdb')

NUCL =  resource_filename(__name__, 'data/1k5i.pdb')
PDB_multiframe = resource_filename(__name__, 'data/nmr_neopetrosiamide.pdb')
PDB_helix = resource_filename(__name__, 'data/A6PA6_alpha.pdb')

GRO = resource_filename(__name__, 'data/adk_oplsaa.gro')
GRO_velocity = resource_filename(__name__, 'data/sample_velocity_file.gro')
PDB = resource_filename(__name__, 'data/adk_oplsaa.pdb')
XTC = resource_filename(__name__, 'data/adk_oplsaa.xtc')
TRR = resource_filename(__name__, 'data/adk_oplsaa.trr')
TPR = resource_filename(__name__, 'data/adk_oplsaa.tpr')

XYZ_psf = resource_filename(__name__, 'data/2r9r-1b.psf')
XYZ_bz2 = resource_filename(__name__, 'data/2r9r-1b.xyz.bz2')
XYZ = resource_filename(__name__, 'data/2r9r-1b.xyz')

PRM = resource_filename(__name__, 'data/ache.prmtop')
TRJ = resource_filename(__name__, 'data/ache.mdcrd')
TRJ_bz2 = resource_filename(__name__, 'data/ache.mdcrd.bz2')

PRMpbc = resource_filename(__name__, 'data/capped-ala.prmtop')
TRJpbc_bz2 = resource_filename(__name__, 'data/capped-ala.mdcrd.bz2')

PRMncdf = resource_filename(__name__, 'data/bala.prmtop')
TRJncdf = resource_filename(__name__, 'data/bala.trj')
NCDF = resource_filename(__name__, 'data/bala.ncdf')

PRM12 = resource_filename(__name__, 'data/anti.top')
TRJ12_bz2 = resource_filename(__name__, 'data/anti_md1.mdcrd.bz2')

PQR = resource_filename(__name__, 'data/adk_open.pqr')

PDBQT_input = resource_filename(__name__, 'data/pdbqt_inputpdbqt.pdbqt')
PDBQT_querypdb = resource_filename(__name__, 'data/pdbqt_querypdb.pdb')

FASTA = resource_filename(__name__, 'data/test.fasta')

PDB_HOLE = resource_filename(__name__, 'data/1grm_single.pdb')
