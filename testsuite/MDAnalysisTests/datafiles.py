# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
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

__all__ = [
    "PSF", "DCD", "CRD",  # CHARMM (AdK example, DIMS trajectory from JMB 2009 paper)
    "PSF_notop", "PSF_BAD",  # Same as PSF but no bonds etc, malformed version of previous
    "DCD_empty",
    "PSF_TRICLINIC", "DCD_TRICLINIC",  # CHARMM c36 new unitcell, NPT 125 TIP3P (box vectors, see Issue 187 for details)
    "PSF_NAMD", "PDB_NAMD",  # NAMD
    "PSF_NAMD_TRICLINIC", "DCD_NAMD_TRICLINIC", # NAMD, triclinic unitcell (Issue 187)
    "PSF_nosegid",  # psf without a segid, Issue 121
    "PDB_small",  # PDB
    "PDB_closed",
    "PDB_multiframe",
    "PDB_helix",
    "XPDB_small",
    "PDB_full",   # PDB 4E43 (full HEADER, TITLE, COMPND, REMARK, altloc)
    "NUCL",  # nucleic acid (PDB)
    "INC_PDB",  # incomplete PDB file (Issue #396)
    "PDB", "GRO", "XTC", "TRR", "TPR", "GRO_velocity",  # Gromacs (AdK)
    "PDB_xvf", "TPR_xvf", "TRR_xvf",  # Gromacs coords/veloc/forces (cobrotoxin, OPLS-AA, Gromacs 4.5.5 tpr)
    "TPR400", "TPR402", "TPR403", "TPR404", "TPR405", "TPR406", "TPR407",
    "TPR450", "TPR451", "TPR452", "TPR453", "TPR454", "TPR455", "TPR455Double",
    "TPR460", "TPR461",
    "PDB_sub_sol", "PDB_sub_dry",  # TRRReader sub selection
    "TRR_sub_sol",
    "XYZ", "XYZ_psf", "XYZ_bz2",
    "XYZ_mini", "XYZ_five", # 3 and 5 atoms xyzs for an easy topology
    "PRM", "TRJ", "TRJ_bz2",  # Amber (no periodic box)
    "INPCRD",
    "PRMpbc", "TRJpbc_bz2",  # Amber (periodic box)
    "PRM12", "TRJ12_bz2",  # Amber (v12 format, Issue 100)
    "PRMncdf", "TRJncdf", "NCDF",  # Amber (netcdf)
    "PFncdf_Top", "PFncdf_Trj", # Amber ncdf with Positions and Forces
    "PQR",  # PQR
    "PDBQT_input",  # PDBQT
    "PDBQT_querypdb",
    "FASTA",  # sequence alignment, Issue 112 + 113
    "PDB_HOLE",  # gramicidin A
    "XTC_HOLE",  # gramicidin A, all frames identical, for Issue 129
    "DMS",
    "CONECT",  # HIV Reverse Transcriptase with inhibitor
    "TRZ", "TRZ_psf",
    "TRIC",
    "merge_protein", "merge_ligand", "merge_water",
    "mol2_molecules", "mol2_molecule", "mol2_broken_molecule",
    "capping_input", "capping_output", "capping_ace", "capping_nma",
    "LAMMPSdata", "trz4data", "LAMMPSdata_mini",
    "unordered_res",  # pdb file with resids non sequential
    "GMS_ASYMOPT",  # GAMESS C1  optimization
    "GMS_SYMOPT",   # GAMESS D4h optimization
    "GMS_ASYMSURF", # GAMESS C1  surface
    "two_water_gro", "two_water_gro_nonames",  # for bond guessing, 2 water molecules, one with weird names
    "DLP_CONFIG", "DLP_CONFIG_order", "DLP_CONFIG_minimal",  # dl_poly 4 config file
    "DLP_HISTORY", "DLP_HISTORY_order", "DLP_HISTORY_minimal",  # dl_poly 4 history file
    "waterPSF","waterDCD","rmsfArray",
    "HoomdXMLdata",
    "Make_Whole",  # for testing the function lib.mdamath.make_whole, has 9 atoms
]

from pkg_resources import resource_filename

PSF = resource_filename(__name__, 'data/adk.psf')
PSF_notop = resource_filename(__name__, 'data/adk_notop.psf')
PSF_BAD = resource_filename(__name__, 'data/adk_notop_BAD.psf')
DCD = resource_filename(__name__, 'data/adk_dims.dcd')
DCD_empty = resource_filename(__name__, 'data/empty.dcd')
CRD = resource_filename(__name__, 'data/adk_open.crd')
PSF_TRICLINIC = resource_filename(__name__, 'data/tip125_tric_C36.psf')
DCD_TRICLINIC = resource_filename(__name__, 'data/tip125_tric_C36.dcd')

PSF_NAMD = resource_filename(__name__, 'data/namd_cgenff.psf')
PDB_NAMD = resource_filename(__name__, 'data/namd_cgenff.pdb')
PSF_NAMD_TRICLINIC = resource_filename(__name__, 'data/SiN_tric_namd.psf')
DCD_NAMD_TRICLINIC = resource_filename(__name__, 'data/SiN_tric_namd.dcd')

PSF_nosegid = resource_filename(__name__, 'data/nosegid.psf')

PDB_small = resource_filename(__name__, 'data/adk_open.pdb')
PDB_closed = resource_filename(__name__, 'data/adk_closed.pdb')

NUCL = resource_filename(__name__, 'data/1k5i.pdb')
INC_PDB = resource_filename(__name__, 'data/incomplete.pdb')
PDB_multiframe = resource_filename(__name__, 'data/nmr_neopetrosiamide.pdb')
PDB_helix = resource_filename(__name__, 'data/A6PA6_alpha.pdb')

GRO = resource_filename(__name__, 'data/adk_oplsaa.gro')
GRO_velocity = resource_filename(__name__, 'data/sample_velocity_file.gro')
PDB = resource_filename(__name__, 'data/adk_oplsaa.pdb')
XTC = resource_filename(__name__, 'data/adk_oplsaa.xtc')
TRR = resource_filename(__name__, 'data/adk_oplsaa.trr')
TPR = resource_filename(__name__, 'data/adk_oplsaa.tpr')
PDB_sub_dry = resource_filename(__name__, 'data/cobrotoxin_dry_neutral_0.pdb')
TRR_sub_sol = resource_filename(__name__, 'data/cobrotoxin.trr')
PDB_sub_sol = resource_filename(__name__, 'data/cobrotoxin.pdb')

PDB_xvf = resource_filename(__name__, 'data/cobrotoxin.pdb')
TPR_xvf = resource_filename(__name__, 'data/cobrotoxin.tpr')
TRR_xvf = resource_filename(__name__, 'data/cobrotoxin.trr')

XPDB_small = resource_filename(__name__, 'data/5digitResid.pdb')
# number is the gromacs version
TPR400 = resource_filename(__name__, 'data/tprs/2lyz_gmx_4.0.tpr')
TPR402 = resource_filename(__name__, 'data/tprs/2lyz_gmx_4.0.2.tpr')
TPR403 = resource_filename(__name__, 'data/tprs/2lyz_gmx_4.0.3.tpr')
TPR404 = resource_filename(__name__, 'data/tprs/2lyz_gmx_4.0.4.tpr')
TPR405 = resource_filename(__name__, 'data/tprs/2lyz_gmx_4.0.5.tpr')
TPR406 = resource_filename(__name__, 'data/tprs/2lyz_gmx_4.0.6.tpr')
TPR407 = resource_filename(__name__, 'data/tprs/2lyz_gmx_4.0.7.tpr')
TPR450 = resource_filename(__name__, 'data/tprs/2lyz_gmx_4.5.tpr')
TPR451 = resource_filename(__name__, 'data/tprs/2lyz_gmx_4.5.1.tpr')
TPR452 = resource_filename(__name__, 'data/tprs/2lyz_gmx_4.5.2.tpr')
TPR453 = resource_filename(__name__, 'data/tprs/2lyz_gmx_4.5.3.tpr')
TPR454 = resource_filename(__name__, 'data/tprs/2lyz_gmx_4.5.4.tpr')
TPR455 = resource_filename(__name__, 'data/tprs/2lyz_gmx_4.5.5.tpr')
# double precision
TPR455Double = resource_filename(__name__, 'data/tprs/drew_gmx_4.5.5.double.tpr')
TPR460 = resource_filename(__name__, 'data/tprs/ab42_gmx_4.6.tpr')
TPR461 = resource_filename(__name__, 'data/tprs/ab42_gmx_4.6.1.tpr')

XYZ_psf = resource_filename(__name__, 'data/2r9r-1b.psf')
XYZ_bz2 = resource_filename(__name__, 'data/2r9r-1b.xyz.bz2')
XYZ = resource_filename(__name__, 'data/2r9r-1b.xyz')
XYZ_mini = resource_filename(__name__, 'data/mini.xyz')
XYZ_five = resource_filename(__name__, 'data/five.xyz')

PRM = resource_filename(__name__, 'data/ache.prmtop')
TRJ = resource_filename(__name__, 'data/ache.mdcrd')
INPCRD = resource_filename(__name__, 'data/test.inpcrd')
TRJ_bz2 = resource_filename(__name__, 'data/ache.mdcrd.bz2')
PFncdf_Top = resource_filename(__name__, 'data/posfor.top')
PFncdf_Trj = resource_filename(__name__, 'data/posfor.ncdf')

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
XTC_HOLE = resource_filename(__name__, 'data/gram_A_identical_frames.xtc')

DMS = resource_filename(__name__, 'data/adk_closed.dms')

CONECT = resource_filename(__name__, 'data/1hvr.pdb')

TRZ = resource_filename(__name__, 'data/trzfile.trz')
TRZ_psf = resource_filename(__name__, 'data/trz_psf.psf')

TRIC = resource_filename(__name__, 'data/dppc_vesicle_hg.gro')

PDB_full = resource_filename(__name__, "data/4E43.pdb")

merge_protein = resource_filename(__name__, "data/merge/2zmm/protein.pdb")
merge_ligand = resource_filename(__name__, "data/merge/2zmm/ligand.pdb")
merge_water = resource_filename(__name__, "data/merge/2zmm/water.pdb")

mol2_molecules = resource_filename(__name__, "data/mol2/Molecules.mol2")
mol2_molecule = resource_filename(__name__, "data/mol2/Molecule.mol2")
mol2_broken_molecule = resource_filename(__name__, "data/mol2/BrokenMolecule.mol2")

capping_input = resource_filename(__name__, "data/capping/aaqaa.gro")
capping_output = resource_filename(__name__, "data/capping/maestro_aaqaa_capped.pdb")
capping_ace = resource_filename(__name__, "data/capping/ace.pdb")
capping_nma = resource_filename(__name__, "data/capping/nma.pdb")

trz4data = resource_filename(__name__, "data/datatest.trz")
LAMMPSdata = resource_filename(__name__, "data/datatest.data")
LAMMPSdata_mini = resource_filename(__name__, "data/mini.data")

LAMMPSdata2 = resource_filename(__name__, "data/lammps2.data")
LAMMPSdcd2 = resource_filename(__name__, "data/lammps2.dcd")

unordered_res = resource_filename(__name__, "data/unordered_res.pdb")

GMS_ASYMOPT       = resource_filename(__name__, "data/gms/c1opt.gms.gz")
GMS_SYMOPT        = resource_filename(__name__, "data/gms/symopt.gms")
GMS_ASYMSURF      = resource_filename(__name__, "data/gms/surf2wat.gms")

two_water_gro = resource_filename(__name__, "data/two_water_gro.gro")
two_water_gro_nonames = resource_filename(__name__, "data/two_water_gro_nonames.gro")

DLP_CONFIG = resource_filename(__name__, "data/dlpoly/CONFIG")
DLP_CONFIG_order = resource_filename(__name__, "data/dlpoly/CONFIG_order")
DLP_CONFIG_minimal = resource_filename(__name__, "data/dlpoly/CONFIG_minimal")
DLP_HISTORY = resource_filename(__name__, "data/dlpoly/HISTORY")
DLP_HISTORY_order = resource_filename(__name__, "data/dlpoly/HISTORY_order")
DLP_HISTORY_minimal = resource_filename(__name__, "data/dlpoly/HISTORY_minimal")

waterPSF = resource_filename(__name__, 'data/watdyn.psf')
waterDCD = resource_filename(__name__, 'data/watdyn.dcd')

rmsfArray = resource_filename(__name__, 'data/adk_oplsaa_CA_rmsf.npy')

HoomdXMLdata = resource_filename(__name__, 'data/C12x64.xml.bz2')

Make_Whole = resource_filename(__name__, 'data/make_whole.gro')
