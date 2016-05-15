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
    "DCD2", # CHARMM (AdK example, DIMS trajectory from PLOS Comput Biol paper)
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
    "PDB_conect",
    "XPDB_small",
    "PDB_full",   # PDB 4E43 (full HEADER, TITLE, COMPND, REMARK, altloc)
    "ALIGN",  # Various way to align atom names in PDB files
    "NUCL",  # nucleic acid (PDB)
    "INC_PDB",  # incomplete PDB file (Issue #396)
    # for testing cryst before/after model headers
    "PDB_cm", "PDB_cm_bz2", "PDB_cm_gz",
    "PDB_mc", "PDB_mc_bz2", "PDB_mc_gz",
    "PDB", "GRO", "XTC", "TRR", "TPR", "GRO_velocity",  # Gromacs (AdK)
    "GRO_incomplete_vels",
    "GRO_large", #atom number truncation at > 100,000 particles, Issue 550
    "PDB_xvf", "TPR_xvf", "TRR_xvf",  # Gromacs coords/veloc/forces (cobrotoxin, OPLS-AA, Gromacs 4.5.5 tpr)
    "PDB_xlserial",
    "TPR400", "TPR402", "TPR403", "TPR404", "TPR405", "TPR406", "TPR407",
    "TPR450", "TPR451", "TPR452", "TPR453", "TPR454", "TPR455", "TPR455Double",
    "TPR460", "TPR461", "TPR502", "TPR504", "TPR505", "TPR510",
    "TPR510_bonded",
    "PDB_sub_sol", "PDB_sub_dry",  # TRRReader sub selection
    "TRR_sub_sol",
    "XTC_sub_sol",
    "XYZ", "XYZ_psf", "XYZ_bz2",
    "XYZ_mini", "XYZ_five", # 3 and 5 atoms xyzs for an easy topology
    "PRM", "TRJ", "TRJ_bz2",  # Amber (no periodic box)
    "INPCRD",
    "PRMpbc", "TRJpbc_bz2",  # Amber (periodic box)
    "PRM7", "NCDFtruncoct",  # Amber (cpptrj test trajectory, see Issue 488)
    "PRM12", "TRJ12_bz2",  # Amber (v12 format, Issue 100)
    "PRMncdf", "TRJncdf", "NCDF",  # Amber (netcdf)
    "PFncdf_Top", "PFncdf_Trj", # Amber ncdf with Positions and Forces
    "PQR",  # PQR
    "PDBQT_input",  # PDBQT
    "PDBQT_querypdb",
    "FASTA",  # sequence alignment, Issue 112 + 113
    "HELANAL_BENDING_MATRIX",  # HELANAL test (from PSF+DCD (AdK) helix 8)
    "PDB_HOLE",  # gramicidin A
    "XTC_HOLE",  # gramicidin A, all frames identical, for Issue 129
    "DMS",
    "CONECT",  # HIV Reverse Transcriptase with inhibitor
    "TRZ", "TRZ_psf",
    "TRIC",
    "XTC_single_frame",
    "XTC_multi_frame",
    "TRR_single_frame",
    "TRR_multi_frame",
    "merge_protein", "merge_ligand", "merge_water",
    "mol2_molecules", "mol2_molecule", "mol2_broken_molecule",
    "mol2_zinc",
    "capping_input", "capping_output", "capping_ace", "capping_nma",
    "contacts_villin_folded", "contacts_villin_unfolded", "contacts_file",
    "LAMMPSdata", "trz4data", "LAMMPSdata_mini",
    "LAMMPSdata2", "LAMMPSdcd2",
    "LAMMPScnt", "LAMMPScnt2",  # triclinic box
    "LAMMPShyd", "LAMMPShyd2",
    "unordered_res",  # pdb file with resids non sequential
    "GMS_ASYMOPT",  # GAMESS C1  optimization
    "GMS_SYMOPT",   # GAMESS D4h optimization
    "GMS_ASYMSURF", # GAMESS C1  surface
    "two_water_gro", "two_water_gro_nonames",  # for bond guessing, 2 water molecules, one with weird names
    "two_water_gro_widebox",  # Issue #548
    "DLP_CONFIG", "DLP_CONFIG_order", "DLP_CONFIG_minimal",  # dl_poly 4 config file
    "DLP_HISTORY", "DLP_HISTORY_order", "DLP_HISTORY_minimal",  # dl_poly 4 history file
    "waterPSF","waterDCD","rmsfArray",
    "HoomdXMLdata",
    "Make_Whole",  # for testing the function lib.mdamath.make_whole, has 9 atoms
    "Plength",
    "COORDINATES_XYZ",
    "COORDINATES_XYZ_BZ2",
    "Martini_membrane_gro", # for testing the leaflet finder
    "COORDINATES_XTC",
    "COORDINATES_TRR",
    "COORDINATES_TOPOLOGY",
    "NUCLsel",
    "GRO_empty_atom", "GRO_missing_atomname", # for testing GROParser exception raise
    "ENT" #for testing ENT file extension
]

from pkg_resources import resource_filename

ENT = resource_filename(__name__, 'data/testENT.ent')
GRO_missing_atomname = resource_filename(__name__, 'data/missing_atomname.gro')
GRO_empty_atom = resource_filename(__name__, 'data/empty_atom.gro')

COORDINATES_XYZ = resource_filename(__name__, 'data/coordinates/test.xyz')
COORDINATES_XYZ_BZ2 = resource_filename(
    __name__, 'data/coordinates/test.xyz.bz2')
COORDINATES_XTC = resource_filename(__name__, 'data/coordinates/test.xtc')
COORDINATES_TRR = resource_filename(__name__, 'data/coordinates/test.trr')
COORDINATES_TOPOLOGY = resource_filename(__name__, 'data/coordinates/test_topology.pdb')

PSF = resource_filename(__name__, 'data/adk.psf')
PSF_notop = resource_filename(__name__, 'data/adk_notop.psf')
PSF_BAD = resource_filename(__name__, 'data/adk_notop_BAD.psf')
DCD = resource_filename(__name__, 'data/adk_dims.dcd')
DCD_empty = resource_filename(__name__, 'data/empty.dcd')
CRD = resource_filename(__name__, 'data/adk_open.crd')
PSF_TRICLINIC = resource_filename(__name__, 'data/tip125_tric_C36.psf')
DCD_TRICLINIC = resource_filename(__name__, 'data/tip125_tric_C36.dcd')
DCD2 = resource_filename(__name__, 'data/adk_dims2.dcd')

PSF_NAMD = resource_filename(__name__, 'data/namd_cgenff.psf')
PDB_NAMD = resource_filename(__name__, 'data/namd_cgenff.pdb')
PSF_NAMD_TRICLINIC = resource_filename(__name__, 'data/SiN_tric_namd.psf')
DCD_NAMD_TRICLINIC = resource_filename(__name__, 'data/SiN_tric_namd.dcd')

PSF_nosegid = resource_filename(__name__, 'data/nosegid.psf')

PDB_small = resource_filename(__name__, 'data/adk_open.pdb')
PDB_closed = resource_filename(__name__, 'data/adk_closed.pdb')

ALIGN = resource_filename(__name__, 'data/align.pdb')
NUCL = resource_filename(__name__, 'data/1k5i.pdb')
INC_PDB = resource_filename(__name__, 'data/incomplete.pdb')
PDB_cm = resource_filename(__name__, 'data/cryst_then_model.pdb')
PDB_cm_gz = resource_filename(__name__, 'data/cryst_then_model.pdb.gz')
PDB_cm_bz2 = resource_filename(__name__, 'data/cryst_then_model.pdb.bz2')
PDB_mc = resource_filename(__name__, 'data/model_then_cryst.pdb')
PDB_mc_gz = resource_filename(__name__, 'data/model_then_cryst.pdb.gz')
PDB_mc_bz2 = resource_filename(__name__, 'data/model_then_cryst.pdb.bz2')
PDB_multiframe = resource_filename(__name__, 'data/nmr_neopetrosiamide.pdb')
PDB_helix = resource_filename(__name__, 'data/A6PA6_alpha.pdb')
PDB_conect = resource_filename(__name__, 'data/conect_parsing.pdb')

GRO = resource_filename(__name__, 'data/adk_oplsaa.gro')
GRO_velocity = resource_filename(__name__, 'data/sample_velocity_file.gro')
GRO_incomplete_vels = resource_filename(__name__, 'data/grovels.gro')
GRO_large = resource_filename(__name__, 'data/bigbox.gro.bz2')
PDB = resource_filename(__name__, 'data/adk_oplsaa.pdb')
XTC = resource_filename(__name__, 'data/adk_oplsaa.xtc')
TRR = resource_filename(__name__, 'data/adk_oplsaa.trr')
TPR = resource_filename(__name__, 'data/adk_oplsaa.tpr')
PDB_sub_dry = resource_filename(__name__, 'data/cobrotoxin_dry_neutral_0.pdb')
TRR_sub_sol = resource_filename(__name__, 'data/cobrotoxin.trr')
XTC_sub_sol = resource_filename(__name__, 'data/cobrotoxin.xtc')
PDB_sub_sol = resource_filename(__name__, 'data/cobrotoxin.pdb')
PDB_xlserial = resource_filename(__name__, 'data/xl_serial.pdb')
XTC_single_frame = resource_filename(
    __name__, 'data/xtc_test_only_single_frame_10_atoms.xtc')
XTC_multi_frame = resource_filename(
    __name__, 'data/xtc_test_only_10_frame_10_atoms.xtc'
)
TRR_single_frame = resource_filename(
    __name__, 'data/trr_test_only_single_frame_10_atoms.trr')
TRR_multi_frame = resource_filename(
    __name__, 'data/trr_test_only_10_frame_10_atoms.trr'
)

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
TPR502 = resource_filename(__name__, 'data/tprs/2lyz_gmx_5.0.2.tpr')
TPR504 = resource_filename(__name__, 'data/tprs/2lyz_gmx_5.0.4.tpr')
TPR505 = resource_filename(__name__, 'data/tprs/2lyz_gmx_5.0.5.tpr')
TPR510 = resource_filename(__name__, 'data/tprs/2lyz_gmx_5.1.tpr')
# double precision
TPR455Double = resource_filename(__name__, 'data/tprs/drew_gmx_4.5.5.double.tpr')
TPR460 = resource_filename(__name__, 'data/tprs/ab42_gmx_4.6.tpr')
TPR461 = resource_filename(__name__, 'data/tprs/ab42_gmx_4.6.1.tpr')
# all bonded interactions
TPR510_bonded = resource_filename(__name__, 'data/tprs/all_bonded/dummy_5.1.tpr')

XYZ_psf = resource_filename(__name__, 'data/2r9r-1b.psf')
XYZ_bz2 = resource_filename(__name__, 'data/2r9r-1b.xyz.bz2')
XYZ = resource_filename(__name__, 'data/2r9r-1b.xyz')
XYZ_mini = resource_filename(__name__, 'data/mini.xyz')
XYZ_five = resource_filename(__name__, 'data/five.xyz')

PRM = resource_filename(__name__, 'data/Amber/ache.prmtop')
TRJ = resource_filename(__name__, 'data/Amber/ache.mdcrd')
INPCRD = resource_filename(__name__, 'data/Amber/test.inpcrd')
TRJ_bz2 = resource_filename(__name__, 'data/Amber/ache.mdcrd.bz2')
PFncdf_Top = resource_filename(__name__, 'data/Amber/posfor.top')
PFncdf_Trj = resource_filename(__name__, 'data/Amber/posfor.ncdf')

PRMpbc = resource_filename(__name__, 'data/Amber/capped-ala.prmtop')
TRJpbc_bz2 = resource_filename(__name__, 'data/Amber/capped-ala.mdcrd.bz2')

PRMncdf = resource_filename(__name__, 'data/Amber/bala.prmtop')
TRJncdf = resource_filename(__name__, 'data/Amber/bala.trj')
NCDF = resource_filename(__name__, 'data/Amber/bala.ncdf')

PRM12 = resource_filename(__name__, 'data/Amber/anti.top')
TRJ12_bz2 = resource_filename(__name__, 'data/Amber/anti_md1.mdcrd.bz2')

PRM7 =  resource_filename(__name__, 'data/Amber/tz2.truncoct.parm7.bz2')
NCDFtruncoct =  resource_filename(__name__, 'data/Amber/tz2.truncoct.nc')

PQR = resource_filename(__name__, 'data/adk_open.pqr')

PDBQT_input = resource_filename(__name__, 'data/pdbqt_inputpdbqt.pdbqt')
PDBQT_querypdb = resource_filename(__name__, 'data/pdbqt_querypdb.pdb')

FASTA = resource_filename(__name__, 'data/test.fasta')
HELANAL_BENDING_MATRIX = resource_filename(__name__, 'data/helanal_bending_matrix_AdK_DIMS_H8.dat')


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
# MOL2 file without substructure field
mol2_zinc = resource_filename(__name__, "data/mol2/zinc_856218.mol2")

capping_input = resource_filename(__name__, "data/capping/aaqaa.gro")
capping_output = resource_filename(__name__, "data/capping/maestro_aaqaa_capped.pdb")
capping_ace = resource_filename(__name__, "data/capping/ace.pdb")
capping_nma = resource_filename(__name__, "data/capping/nma.pdb")

contacts_villin_folded = resource_filename(__name__, "data/contacts/villin_folded.gro.bz2")
contacts_villin_unfolded = resource_filename(__name__, "data/contacts/villin_unfolded.gro.bz2")
contacts_file = resource_filename(__name__, "data/contacts/2F4K_qlist5_remap.dat")

trz4data = resource_filename(__name__, "data/lammps/datatest.trz")
LAMMPSdata = resource_filename(__name__, "data/lammps/datatest.data")
LAMMPSdata_mini = resource_filename(__name__, "data/lammps/mini.data")
LAMMPSdata2 = resource_filename(__name__, "data/lammps/ifabp_apo_100mM.data.bz2")
LAMMPSdcd2 = resource_filename(__name__, "data/lammps/ifabp_apo_100mM.dcd")
LAMMPScnt = resource_filename(__name__, "data/lammps/cnt-hexagonal-class1.data")
LAMMPScnt2 = resource_filename(__name__, "data/lammps/cnt-hexagonal-class1.data2")
LAMMPShyd = resource_filename(__name__, "data/lammps/hydrogen-class1.data")
LAMMPShyd2 = resource_filename(__name__, "data/lammps/hydrogen-class1.data2")

unordered_res = resource_filename(__name__, "data/unordered_res.pdb")

GMS_ASYMOPT       = resource_filename(__name__, "data/gms/c1opt.gms.gz")
GMS_SYMOPT        = resource_filename(__name__, "data/gms/symopt.gms")
GMS_ASYMSURF      = resource_filename(__name__, "data/gms/surf2wat.gms")

two_water_gro = resource_filename(__name__, "data/two_water_gro.gro")
two_water_gro_nonames = resource_filename(__name__, "data/two_water_gro_nonames.gro")
two_water_gro_widebox = resource_filename(__name__, "data/two_water_gro_widebox.gro")

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

Plength = resource_filename(__name__, 'data/plength.gro')
Martini_membrane_gro = resource_filename(__name__, 'data/martini_dppc_chol_bilayer.gro')

# Contains one of each residue in 'nucleic' selections
NUCLsel = resource_filename(__name__, 'data/nucl_res.pdb')


# This should be the last line: clean up namespace
del resource_filename
