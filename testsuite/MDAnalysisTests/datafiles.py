# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

"""
Location of data files for the MDAnalysis unit tests
====================================================

Real MD simulation data are stored in the ```` sub
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
    "PSF_NAMD_GBIS", "DCD_NAMD_GBIS",  # NAMD, implicit solvent, 100 steps, #1819
    "PSF_nosegid",  # psf without a segid, Issue 121
    "PSF_cmap",  # ala3 PSF from ParmEd test files with cmap
    "PDB_small",  # PDB
    "PDB_closed",
    "PDB_multiframe",
    "PDB_helix",
    "PDB_conect",
    "PDB_conect2TER",  # Conect record to a TER entry (Issue 936)
    "PDB_singleconect",  # Conect record with one entry (Issue 937)
    "PDB_icodes",  # stripped down version of 1osm, has icodes!
    "XPDB_small",
    "PDB_full",   # PDB 4E43 (full HEADER, TITLE, COMPND, REMARK, altloc)
    "ALIGN",  # Various way to align atom names in PDB files
    "RNA_PSF", "RNA_PDB",  # nucleic acid (PDB 1K5I in CHARMM36m)
    "INC_PDB",  # incomplete PDB file (Issue #396)
    # for testing cryst before/after model headers
    "PDB_cm", "PDB_cm_bz2", "PDB_cm_gz",
    "PDB_mc", "PDB_mc_bz2", "PDB_mc_gz",
    "PDB_chainidnewres",  # Issue 1110
    "PDB_sameresid_diffresname", #Case where two residues share the same resid
    "PDB_chainidrepeat",  # Issue #1107
    "PDB", "GRO", "XTC", "TRR", "TPR", "GRO_velocity",  # Gromacs (AdK)
    "GRO_incomplete_vels",
    "COORDINATES_GRO_BZ2",
    "GRO_large", #atom number truncation at > 100,000 particles, Issue 550
    "GRO_residwrap",  # resids wrapping because of 5 digit field (Issue #728)
    "GRO_residwrap_0base",  # corner case of #728 with resid=0 for first atom
    "GRO_sameresid_diffresname", # Case where two residues share the same resid
    "PDB_xvf", "TPR_xvf", "TRR_xvf",  # Gromacs coords/veloc/forces (cobrotoxin, OPLS-AA, Gromacs 4.5.5 tpr)
    "H5MD_xvf",  # TPR_xvf + TRR_xvf converted to h5md format
    "XVG_BZ2",  # Compressed xvg file about cobrotoxin
    "PDB_xlserial",
    "TPR400", "TPR402", "TPR403", "TPR404", "TPR405", "TPR406", "TPR407",
    "TPR450", "TPR451", "TPR452", "TPR453", "TPR454", "TPR455", "TPR455Double",
    "TPR460", "TPR461", "TPR502", "TPR504", "TPR505", "TPR510", "TPR2016",
    "TPR2018", "TPR2019B3", "TPR2020B2", "TPR2020", "TPR2020Double",
    "TPR2021", "TPR2021Double", "TPR2022RC1", "TPR2023",
    "TPR510_bonded", "TPR2016_bonded", "TPR2018_bonded", "TPR2019B3_bonded",
    "TPR2020B2_bonded", "TPR2020_bonded", "TPR2020_double_bonded",
    "TPR2021_bonded", "TPR2021_double_bonded", "TPR2022RC1_bonded",
    "TPR334_bonded", "TPR2023_bonded",
    "TPR_EXTRA_2021", "TPR_EXTRA_2020", "TPR_EXTRA_2018",
    "TPR_EXTRA_2016", "TPR_EXTRA_407", "TPR_EXTRA_2022RC1",
    "TPR_EXTRA_2023",
    "PDB_sub_sol", "PDB_sub_dry",  # TRRReader sub selection
    "TRR_sub_sol",
    "XTC_sub_sol",
    "XYZ", "XYZ_psf", "XYZ_bz2",
    "XYZ_mini", "XYZ_five", # 3 and 5 atoms xyzs for an easy topology
    "TXYZ", "ARC", "ARC_PBC",       # Tinker files
    "PRM", "TRJ", "TRJ_bz2",  # Amber (no periodic box)
    "INPCRD",
    "PRMpbc", "TRJpbc_bz2",  # Amber (periodic box)
    "PRM7", "NCDFtruncoct",  # Amber (cpptrj test trajectory, see Issue 488)
    "PRM12", "TRJ12_bz2",  # Amber (v12 format, Issue 100)
    "PRMncdf", "TRJncdf", "NCDF",  # Amber (netcdf)
    "PFncdf_Top", "PFncdf_Trj", # Amber ncdf with Positions and Forces
    "PRMcs", # Amber (format, Issue 1331)
    "PRMNCRST", # Amber ncrst with positions/forces/velocities
    "PRM_NCBOX", "TRJ_NCBOX", # Amber parm7 + nc w/ pos/forces/vels/box
    "PRMNEGATIVE", # Amber negative ATOMIC_NUMBER (Issue 2306)
    "PRMErr1", "PRMErr2", "PRMErr3", # Amber TOP files to check raised errors
    "PRM_UreyBradley", # prmtop from ParmEd test files with Urey-Bradley angles
    "PRM7_ala2", "RST7_ala2",  # prmtop and rst files from ParmEd example files
    "PRM19SBOPC", #  prmtop w/ ff19SB CMAP terms and OPC water (Issue #2449)
    "PQR",  # PQR v1
    "PQR_icodes",  # PQR v2 with icodes
    "PDBQT_input",  # PDBQT
    "PDBQT_querypdb",
    "PDB_multipole",
    "FASTA",  # sequence alignment, Issue 112 + 113
    "HELANAL_BENDING_MATRIX",  # HELANAL test (from PSF+DCD (AdK) helix 8)
    "HELANAL_BENDING_MATRIX_SUBSET", # As above, slice of frames 10 to 79
    "PDB_HOLE",  # gramicidin A
    "MULTIPDB_HOLE", # gramicidin A, normal mode 7 from ElNemo
    "DMS",
    "DMS_DOMAINS",  # ADK closed with multiple segids
    "DMS_NO_SEGID",  # ADK closed with no segids or chains
    "CONECT",  # HIV Reverse Transcriptase with inhibitor
    "CONECT_ERROR",  # PDB file with corrupt CONECT
    "TRZ", "TRZ_psf",
    "TRIC",
    "XTC_multi_frame",
    "TRR_multi_frame",
    "TNG_traj",  # TNG trajectory from GROMACS physical validation testsuite, longish trajectory
    "TNG_traj_gro",  # topology for argon_npt_compressed_traj
    "TNG_traj_uneven_blocks",  # TNG trajectory with pos and vel deposited on different strides
    "TNG_traj_vels_forces",  # similar to above but with velocities and forces
    "merge_protein", "merge_ligand", "merge_water",
    "mol2_molecules", "mol2_molecule", "mol2_broken_molecule",
    "mol2_zinc", "mol2_comments_header", "mol2_ligand", "mol2_sodium_ion",
    "capping_input", "capping_output", "capping_ace", "capping_nma",
    "contacts_villin_folded", "contacts_villin_unfolded", "contacts_file",
    "LAMMPSdata", "trz4data", "LAMMPSdata_mini",
    "LAMMPSdata2", "LAMMPSdcd2",
    "LAMMPScnt", "LAMMPScnt2",  # triclinic box
    "LAMMPShyd", "LAMMPShyd2",
    "LAMMPSdata_deletedatoms",  # with deleted atoms
    "LAMMPSdata_triclinic", # lammpsdata file to test triclinic dimension parsing, albite with most atoms deleted
    "LAMMPSdata_PairIJ",  # lammps datafile with a PairIJ Coeffs section
    "LAMMPSDUMP",
    "LAMMPSDUMP_long",  # lammpsdump file with a few zeros sprinkled in the first column first frame
    "LAMMPSDUMP_allcoords",  # lammpsdump file with all coordinate conventions (x,xs,xu,xsu) present, from LAMMPS rdf example
    "LAMMPSDUMP_nocoords",  # lammpsdump file with no coordinates
    "LAMMPSDUMP_triclinic", # lammpsdump file to test triclinic dimension parsing, albite with most atoms deleted
    "LAMMPSDUMP_image_vf", # Lammps dump file with image flags, velocities, and forces.
    "LAMMPS_image_vf", # Lammps data file to go with LAMMPSDUMP_image_vf
    "unordered_res",  # pdb file with resids non sequential
    "GMS_ASYMOPT",  # GAMESS C1  optimization
    "GMS_SYMOPT",   # GAMESS D4h optimization
    "GMS_ASYMSURF", # GAMESS C1  surface
    "two_water_gro", "two_water_gro_nonames",  # for bond guessing, 2 water molecules, one with weird names
    "two_water_gro_multiframe",
    "two_water_gro_widebox",  # Issue #548
    "DLP_CONFIG", "DLP_CONFIG_order", "DLP_CONFIG_minimal",  # dl_poly 4 config file
    "DLP_HISTORY", "DLP_HISTORY_order", "DLP_HISTORY_minimal",  # dl_poly 4 history file
    "DLP_HISTORY_minimal_cell", # dl_poly 4 history file with cell parameters
    "DLP_HISTORY_classic",  # dl_poly classic history file
    "waterPSF","waterDCD","rmsfArray",
    "HoomdXMLdata",
    "Make_Whole",  # for testing the function lib.mdamath.make_whole, has 9 atoms
    "fullerene",  # for make_whole, a nice friendly C60 with bonds
    "Plength",
    "COORDINATES_XYZ",
    "COORDINATES_XYZ_BZ2",
    "COORDINATES_GRO",
    "COORDINATES_GRO_INCOMPLETE_VELOCITY",
    "Martini_membrane_gro", # for testing the leaflet finder
    "COORDINATES_XTC",
    "COORDINATES_TRR",
    "COORDINATES_TNG",
    "COORDINATES_H5MD",
    "COORDINATES_DCD",
    "COORDINATES_TOPOLOGY",
    "NUCLsel",
    "GRO_empty_atom", "GRO_missing_atomname", # for testing GROParser exception raise
    "ENT", #for testing ENT file extension
    "RANDOM_WALK",
    "RANDOM_WALK_TOPO", # garbage topology to go along with XTC positions above
    "AUX_XVG", "XVG_BAD_NCOL", #for testing .xvg auxiliary reader
    "AUX_XVG_LOWF", "AUX_XVG_HIGHF",
    "AUX_EDR", "AUX_EDR_TPR",
    "AUX_EDR_XTC", "AUX_EDR_RAW",
    "AUX_EDR_SINGLE_FRAME",  # for testing .edr auxiliary reader
    "MMTF", "MMTF_gz", 'MMTF_skinny',  # skinny - some optional fields stripped out
    "MMTF_skinny2",
    "ALIGN_BOUND",  # two component bound system
    "ALIGN_UNBOUND", # two component unbound system
    "legacy_DCD_ADK_coords", # frames 5 and 29 read in for adk_dims.dcd using legacy DCD reader
    "legacy_DCD_NAMD_coords", # frame 0 read in for SiN_tric_namd.dcd using legacy DCD reader
    "legacy_DCD_c36_coords", # frames 1 and 4 read in for tip125_tric_C36.dcd using legacy DCD reader
    "GSD", "GSD_bonds", "GSD_long",
    "GRO_MEMPROT", "XTC_MEMPROT", # YiiP transporter in POPE:POPG lipids with Na+, Cl-, Zn2+ dummy model without water
    "DihedralArray", "DihedralsArray", # time series of single dihedral
    "RamaArray", "GLYRamaArray", # time series of phi/psi angles
    "JaninArray", "LYSJaninArray", # time series of chi1/chi2 angles
    "PDB_rama", "PDB_janin", # for testing failures of Ramachandran and Janin classes
    "BATArray", # time series of bond-angle-torsion coordinates array from Molecule_comments_header.mol2
    # DOS line endings
    "WIN_PDB_multiframe", "WIN_DLP_HISTORY", "WIN_TRJ", "WIN_LAMMPSDUMP", "WIN_ARC",
    "GRO_huge_box", # for testing gro parser with hige box sizes
    "ITP", # for GROMACS generated itps
    "ITP_nomass", # for ATB generated itps
    "ITP_atomtypes",  # atom definitions to check atomtyes section parsing
    "ITP_charges", # atom definitions to test custom particle charge parsing.
    "NAMDBIN", # for NAMD generated binary file
    "ITP_edited", # to check different directives are read properly
    "ITP_tip5p", # tip5p water from opls-aa, edited with additional keywords
    "ITP_spce", # spce water from gromos54a7, edited with additional keywords,
    "GMX_TOP", # 2 ala10 chains + 3 spc water
    "GMX_DIR", # GROMACS directory
    "GMX_TOP_BAD", # file with an #include that doesn't exist
    "ITP_no_endif", # file missing an #endif
    "PDB_CRYOEM_BOX", # Issue 2599, Issue #2679, PR #2685
    "PDB_CHECK_RIGHTHAND_PA", # for testing right handedness of principal_axes
    "MMTF_NOCRYST", # File with meaningless CRYST1 record (Issue #2679, PR #2685)
    "FHIAIMS", # to test FHIAIMS coordinate files
    "SDF_molecule",  # MDL SDFile for rdkit
    "PDBX",  # PDBxfile
    "PDB_elements",  # PDB file with elements
    "PDB_charges",  # PDB file with formal charges
    "SURFACE_PDB",  # 111 FCC lattice topology for NSGrid bug #2345
    "SURFACE_TRR",  # full precision coordinates for NSGrid bug #2345
]

from importlib import resources
import MDAnalysisTests.data

data_ref = resources.files('MDAnalysisTests.data')

WIN_PDB_multiframe = data_ref / 'windows/WIN_nmr_neopetrosiamide.pdb'
WIN_DLP_HISTORY = data_ref / 'windows/WIN_HISTORY'
WIN_TRJ = data_ref / 'windows/WIN_ache.mdcrd'
WIN_ARC = data_ref / 'windows/WIN_test.arc'
WIN_LAMMPSDUMP = data_ref / 'windows/WIN_wat.lammpstrj'

legacy_DCD_NAMD_coords = data_ref / 'legacy_DCD_NAMD_coords.npy'
legacy_DCD_ADK_coords = data_ref / 'legacy_DCD_adk_coords.npy'
legacy_DCD_c36_coords = data_ref / 'legacy_DCD_c36_coords.npy'
AUX_XVG_LOWF = data_ref / 'test_lowf.xvg'
AUX_XVG_HIGHF = data_ref / 'test_highf.xvg'
XVG_BAD_NCOL = data_ref / 'bad_num_col.xvg'
AUX_XVG = data_ref / 'test.xvg'
AUX_EDR = data_ref / 'test.edr'
AUX_EDR_RAW = data_ref / 'aux_edr_raw.txt'
AUX_EDR_TPR = data_ref / 'aux_edr.tpr'
AUX_EDR_XTC = data_ref / 'aux_edr.xtc'
AUX_EDR_SINGLE_FRAME = data_ref / 'single_frame.edr'
ENT = data_ref / 'testENT.ent'
GRO_missing_atomname = data_ref / 'missing_atomname.gro'
GRO_empty_atom = data_ref / 'empty_atom.gro'
GRO_huge_box = data_ref / 'huge_box.gro'

COORDINATES_GRO = data_ref / 'coordinates/test.gro'
COORDINATES_GRO_INCOMPLETE_VELOCITY = data_ref / 'coordinates/test_incomplete_vel.gro'
COORDINATES_GRO_BZ2 = data_ref / 'coordinates/test.gro.bz2'
COORDINATES_XYZ = data_ref / 'coordinates/test.xyz'
COORDINATES_XYZ_BZ2 = data_ref /  'coordinates/test.xyz.bz2'
COORDINATES_XTC = data_ref / 'coordinates/test.xtc'
COORDINATES_TRR = data_ref / 'coordinates/test.trr'
COORDINATES_TNG = data_ref / 'coordinates/test.tng'
COORDINATES_H5MD = data_ref / 'coordinates/test.h5md'
COORDINATES_DCD = data_ref / 'coordinates/test.dcd'
COORDINATES_TOPOLOGY = data_ref / 'coordinates/test_topology.pdb'

PSF = data_ref / 'adk.psf'
PSF_notop = data_ref / 'adk_notop.psf'
PSF_BAD = data_ref / 'adk_notop_BAD.psf'
DCD = data_ref / 'adk_dims.dcd'
DCD_empty = data_ref / 'empty.dcd'
CRD = data_ref / 'adk_open.crd'
PSF_TRICLINIC = data_ref / 'tip125_tric_C36.psf'
DCD_TRICLINIC = data_ref / 'tip125_tric_C36.dcd'
DCD2 = data_ref / 'adk_dims2.dcd'

PSF_NAMD = data_ref / 'namd_cgenff.psf'
PDB_NAMD = data_ref / 'namd_cgenff.pdb'
PDB_multipole = data_ref / 'water_methane_acetic-acid_ammonia.pdb'
PSF_NAMD_TRICLINIC = data_ref / 'SiN_tric_namd.psf'
DCD_NAMD_TRICLINIC = data_ref / 'SiN_tric_namd.dcd'
PSF_NAMD_GBIS = data_ref / 'adk_closed_NAMD.psf'
DCD_NAMD_GBIS = data_ref / 'adk_gbis_tmd-fast1_NAMD.dcd'

PSF_nosegid = data_ref / 'nosegid.psf'

PSF_cmap = data_ref / 'parmed_ala3.psf'

PDB_small = data_ref / 'adk_open.pdb'
PDB_closed = data_ref / 'adk_closed.pdb'

ALIGN = data_ref / 'align.pdb'
RNA_PSF = data_ref / 'analysis/1k5i_c36.psf.gz'
RNA_PDB = data_ref / 'analysis/1k5i_c36.pdb.gz'
INC_PDB = data_ref / 'incomplete.pdb'
PDB_cm = data_ref / 'cryst_then_model.pdb'
PDB_cm_gz = data_ref / 'cryst_then_model.pdb.gz'
PDB_cm_bz2 = data_ref / 'cryst_then_model.pdb.bz2'
PDB_mc = data_ref / 'model_then_cryst.pdb'
PDB_mc_gz = data_ref / 'model_then_cryst.pdb.gz'
PDB_mc_bz2 = data_ref / 'model_then_cryst.pdb.bz2'
PDB_chainidnewres = data_ref / 'chainIDnewres.pdb.gz'
PDB_sameresid_diffresname = data_ref / 'sameresid_diffresname.pdb'
PDB_chainidrepeat = data_ref / 'chainIDrepeat.pdb.gz'
PDB_multiframe = data_ref / 'nmr_neopetrosiamide.pdb'
PDB_helix = data_ref / 'A6PA6_alpha.pdb'
PDB_conect = data_ref / 'conect_parsing.pdb'
PDB_conect2TER = data_ref / 'CONECT2TER.pdb'
PDB_singleconect = data_ref / 'SINGLECONECT.pdb'
PDB_icodes = data_ref / '1osm.pdb.gz'
PDB_CRYOEM_BOX = data_ref / '5a7u.pdb'
PDB_CHECK_RIGHTHAND_PA = data_ref / '6msm.pdb.bz2'
FHIAIMS = data_ref / 'fhiaims.in'

GRO = data_ref / 'adk_oplsaa.gro'
GRO_velocity = data_ref / 'sample_velocity_file.gro'
GRO_incomplete_vels = data_ref / 'grovels.gro'
GRO_large = data_ref / 'bigbox.gro.bz2'
GRO_residwrap = data_ref / 'residwrap.gro'
GRO_residwrap_0base = data_ref / 'residwrap_0base.gro'
GRO_sameresid_diffresname = data_ref / 'sameresid_diffresname.gro'
PDB = data_ref / 'adk_oplsaa.pdb'
XTC = data_ref / 'adk_oplsaa.xtc'
TRR = data_ref / 'adk_oplsaa.trr'
TPR = data_ref / 'adk_oplsaa.tpr'
PDB_sub_dry = data_ref / 'cobrotoxin_dry_neutral_0.pdb'
TRR_sub_sol = data_ref / 'cobrotoxin.trr'
XTC_sub_sol = data_ref / 'cobrotoxin.xtc'
PDB_sub_sol = data_ref / 'cobrotoxin.pdb'
PDB_xlserial = data_ref / 'xl_serial.pdb'
GRO_MEMPROT = data_ref / 'analysis/YiiP_lipids.gro.gz'
XTC_MEMPROT = data_ref / 'analysis/YiiP_lipids.xtc'
XTC_multi_frame = data_ref / 'xtc_test_only_10_frame_10_atoms.xtc'
TRR_multi_frame = data_ref / 'trr_test_only_10_frame_10_atoms.trr'
TNG_traj = data_ref / 'argon_npt_compressed.tng'
TNG_traj_gro = data_ref / 'argon_npt_compressed.gro.gz'
TNG_traj_uneven_blocks = data_ref / 'argon_npt_compressed_uneven.tng'
TNG_traj_vels_forces = data_ref / 'argon_npt_compressed_vels_forces.tng'
PDB_xvf = data_ref / 'cobrotoxin.pdb'
TPR_xvf = data_ref / 'cobrotoxin.tpr'
TRR_xvf = data_ref / 'cobrotoxin.trr'
H5MD_xvf = data_ref / 'cobrotoxin.h5md'
XVG_BZ2 = data_ref / 'cobrotoxin_protein_forces.xvg.bz2'

XPDB_small = data_ref / '5digitResid.pdb'
# number is the gromacs version
TPR400 = data_ref / 'tprs/2lyz_gmx_4.0.tpr'
TPR402 = data_ref / 'tprs/2lyz_gmx_4.0.2.tpr'
TPR403 = data_ref / 'tprs/2lyz_gmx_4.0.3.tpr'
TPR404 = data_ref / 'tprs/2lyz_gmx_4.0.4.tpr'
TPR405 = data_ref / 'tprs/2lyz_gmx_4.0.5.tpr'
TPR406 = data_ref / 'tprs/2lyz_gmx_4.0.6.tpr'
TPR407 = data_ref / 'tprs/2lyz_gmx_4.0.7.tpr'
TPR450 = data_ref / 'tprs/2lyz_gmx_4.5.tpr'
TPR451 = data_ref / 'tprs/2lyz_gmx_4.5.1.tpr'
TPR452 = data_ref / 'tprs/2lyz_gmx_4.5.2.tpr'
TPR453 = data_ref / 'tprs/2lyz_gmx_4.5.3.tpr'
TPR454 = data_ref / 'tprs/2lyz_gmx_4.5.4.tpr'
TPR455 = data_ref / 'tprs/2lyz_gmx_4.5.5.tpr'
TPR502 = data_ref / 'tprs/2lyz_gmx_5.0.2.tpr'
TPR504 = data_ref / 'tprs/2lyz_gmx_5.0.4.tpr'
TPR505 = data_ref / 'tprs/2lyz_gmx_5.0.5.tpr'
TPR510 = data_ref / 'tprs/2lyz_gmx_5.1.tpr'
TPR2016 = data_ref / 'tprs/2lyz_gmx_2016.tpr'
TPR2018 = data_ref / 'tprs/2lyz_gmx_2018.tpr'
TPR2019B3 = data_ref / 'tprs/2lyz_gmx_2019-beta3.tpr'
TPR2020B2 = data_ref / 'tprs/2lyz_gmx_2020-beta2.tpr'
TPR2020 = data_ref / 'tprs/2lyz_gmx_2020.tpr'
TPR2021 = data_ref / 'tprs/2lyz_gmx_2021.tpr'
TPR2022RC1 = data_ref / 'tprs/2lyz_gmx_2022-rc1.tpr'
TPR2023 = data_ref / 'tprs/2lyz_gmx_2023.tpr'
# double precision
TPR455Double = data_ref / 'tprs/drew_gmx_4.5.5.double.tpr'
TPR460 = data_ref / 'tprs/ab42_gmx_4.6.tpr'
TPR461 = data_ref / 'tprs/ab42_gmx_4.6.1.tpr'
TPR2020Double = data_ref / 'tprs/2lyz_gmx_2020_double.tpr'
TPR2021Double = data_ref / 'tprs/2lyz_gmx_2021_double.tpr'
# all bonded interactions
TPR334_bonded = data_ref / 'tprs/all_bonded/dummy_3.3.4.tpr'
TPR510_bonded = data_ref / 'tprs/all_bonded/dummy_5.1.tpr'
TPR2016_bonded = data_ref / 'tprs/all_bonded/dummy_2016.tpr'
TPR2018_bonded = data_ref / 'tprs/all_bonded/dummy_2018.tpr'
TPR2019B3_bonded = data_ref / 'tprs/all_bonded/dummy_2019-beta3.tpr'
TPR2020B2_bonded = data_ref / 'tprs/all_bonded/dummy_2020-beta2.tpr'
TPR2020_bonded = data_ref / 'tprs/all_bonded/dummy_2020.tpr'
TPR2020_double_bonded = data_ref / 'tprs/all_bonded/dummy_2020_double.tpr'
TPR2021_bonded = data_ref / 'tprs/all_bonded/dummy_2021.tpr'
TPR2021_double_bonded = data_ref / 'tprs/all_bonded/dummy_2021_double.tpr'
TPR2022RC1_bonded = data_ref / 'tprs/all_bonded/dummy_2022-rc1.tpr'
TPR2023_bonded = data_ref / 'tprs/all_bonded/dummy_2023.tpr'
# all interactions
TPR_EXTRA_2023 = data_ref / 'tprs/virtual_sites/extra-interactions-2023.tpr'
TPR_EXTRA_2022RC1 = data_ref / 'tprs/virtual_sites/extra-interactions-2022-rc1.tpr'
TPR_EXTRA_2021 = data_ref / 'tprs/virtual_sites/extra-interactions-2021.tpr'
TPR_EXTRA_2020 = data_ref / 'tprs/virtual_sites/extra-interactions-2020.tpr'
TPR_EXTRA_2018 = data_ref / 'tprs/virtual_sites/extra-interactions-2018.tpr'
TPR_EXTRA_2016 = data_ref / 'tprs/virtual_sites/extra-interactions-2016.3.tpr'
TPR_EXTRA_407 = data_ref / 'tprs/virtual_sites/extra-interactions-4.0.7.tpr'

XYZ_psf = data_ref / '2r9r-1b.psf'
XYZ_bz2 = data_ref / '2r9r-1b.xyz.bz2'
XYZ = data_ref / '2r9r-1b.xyz'
XYZ_mini = data_ref / 'mini.xyz'
XYZ_five = data_ref / 'five.xyz'
TXYZ = data_ref / 'coordinates/test.txyz'
ARC = data_ref / 'coordinates/test.arc'
ARC_PBC = data_ref / 'coordinates/new_hexane.arc'

PRM = data_ref / 'Amber/ache.prmtop'
TRJ = data_ref / 'Amber/ache.mdcrd'
INPCRD = data_ref / 'Amber/test.inpcrd'
TRJ_bz2 = data_ref / 'Amber/ache.mdcrd.bz2'
PFncdf_Top = data_ref / 'Amber/posfor.top'
PFncdf_Trj = data_ref / 'Amber/posfor.ncdf'

PRMpbc = data_ref / 'Amber/capped-ala.prmtop'
TRJpbc_bz2 = data_ref / 'Amber/capped-ala.mdcrd.bz2'

PRMncdf = data_ref / 'Amber/bala.prmtop'
TRJncdf = data_ref / 'Amber/bala.trj'
NCDF = data_ref / 'Amber/bala.ncdf'

PRM12 = data_ref / 'Amber/anti.top'
TRJ12_bz2 = data_ref / 'Amber/anti_md1.mdcrd.bz2'

PRM7 =  data_ref / 'Amber/tz2.truncoct.parm7.bz2'
NCDFtruncoct =  data_ref / 'Amber/tz2.truncoct.nc'

PRMcs = data_ref / 'Amber/chitosan.prmtop'

PRMNCRST = data_ref / 'Amber/ace_mbondi3.parm7'

PRM_NCBOX = data_ref / 'Amber/ace_tip3p.parm7'
TRJ_NCBOX = data_ref / 'Amber/ace_tip3p.nc'

PRMNEGATIVE = data_ref / 'Amber/ace_mbondi3.negative.parm7'

PRMErr1 = data_ref / 'Amber/ace_mbondi3.error1.parm7'
PRMErr2 = data_ref / 'Amber/ace_mbondi3.error2.parm7'
PRMErr3 = data_ref / 'Amber/ace_mbondi3.error3.parm7'

PRM_UreyBradley = data_ref / 'Amber/parmed_fad.prmtop'
PRM7_ala2 = data_ref / 'Amber/parmed_ala2_solv.parm7'
RST7_ala2 = data_ref / 'Amber/parmed_ala2_solv.rst7'

PRM19SBOPC = data_ref / 'Amber/ala.ff19SB.OPC.parm7.bz2'

PQR = data_ref / 'adk_open.pqr'
PQR_icodes = data_ref / '1A2C.pqr'

PDBQT_input = data_ref / 'pdbqt_inputpdbqt.pdbqt'
PDBQT_querypdb = data_ref / 'pdbqt_querypdb.pdb'

FASTA = data_ref / 'test.fasta'
HELANAL_BENDING_MATRIX = data_ref / 'helanal_bending_matrix_AdK_DIMS_H8.dat'
HELANAL_BENDING_MATRIX_SUBSET = data_ref / 'helanal_bending_matrix_AdK_DIMS_H8_frames10to79.dat'

PDB_HOLE = data_ref / '1grm_single.pdb'
MULTIPDB_HOLE = data_ref / '1grm_elNemo_mode7.pdb.bz2'

DMS = data_ref / 'adk_closed.dms'
DMS_DOMAINS = data_ref / 'adk_closed_domains.dms'
DMS_NO_SEGID = data_ref / 'adk_closed_no_segid.dms'

CONECT = data_ref / '1hvr.pdb'
CONECT_ERROR = data_ref / 'conect_error.pdb'

TRZ = data_ref / 'trzfile.trz'
TRZ_psf = data_ref / 'trz_psf.psf'

TRIC = data_ref / 'dppc_vesicle_hg.gro'

PDB_full = data_ref / "4E43.pdb"

merge_protein = data_ref / "merge/2zmm/protein.pdb"
merge_ligand = data_ref / "merge/2zmm/ligand.pdb"
merge_water = data_ref / "merge/2zmm/water.pdb"

mol2_molecules = data_ref / "mol2/Molecules.mol2"
mol2_molecule = data_ref / "mol2/Molecule.mol2"
mol2_ligand = data_ref / "mol2/Ligand.mol2"
mol2_broken_molecule = data_ref / "mol2/BrokenMolecule.mol2"
mol2_comments_header = data_ref / "mol2/Molecule_comments_header.mol2"
# MOL2 file without substructure field
mol2_zinc = data_ref / "mol2/zinc_856218.mol2"
# MOL2 file without bonds
mol2_sodium_ion = data_ref / "mol2/sodium_ion.mol2"

capping_input = data_ref / "capping/aaqaa.gro"
capping_output = data_ref / "capping/maestro_aaqaa_capped.pdb"
capping_ace = data_ref / "capping/ace.pdb"
capping_nma = data_ref / "capping/nma.pdb"

contacts_villin_folded = data_ref / "contacts/villin_folded.gro.bz2"
contacts_villin_unfolded = data_ref / "contacts/villin_unfolded.gro.bz2"
contacts_file = data_ref / "contacts/2F4K_qlist5_remap.dat"

trz4data = data_ref / "lammps/datatest.trz"
LAMMPSdata = data_ref / "lammps/datatest.data"
LAMMPSdata_mini = data_ref / "lammps/mini.data"
LAMMPSdata2 = data_ref / "lammps/ifabp_apo_100mM.data.bz2"
LAMMPSdcd2 = data_ref / "lammps/ifabp_apo_100mM.dcd"
LAMMPScnt = data_ref / "lammps/cnt-hexagonal-class1.data"
LAMMPScnt2 = data_ref / "lammps/cnt-hexagonal-class1.data2"
LAMMPShyd = data_ref / "lammps/hydrogen-class1.data"
LAMMPShyd2 = data_ref / "lammps/hydrogen-class1.data2"
LAMMPSdata_deletedatoms = data_ref / 'lammps/deletedatoms.data'
LAMMPSdata_triclinic = data_ref / "lammps/albite_triclinic.data"
LAMMPSdata_PairIJ = data_ref / "lammps/pairij_coeffs.data.bz2"
LAMMPSDUMP = data_ref / "lammps/wat.lammpstrj.bz2"
LAMMPSDUMP_long = data_ref / "lammps/wat.lammpstrj_long.bz2"
LAMMPSDUMP_allcoords = data_ref / "lammps/spce_all_coords.lammpstrj.bz2"
LAMMPSDUMP_nocoords = data_ref / "lammps/spce_no_coords.lammpstrj.bz2"
LAMMPSDUMP_triclinic = data_ref / "lammps/albite_triclinic.dump"
LAMMPSDUMP_image_vf = data_ref / "lammps/image_vf.lammpstrj"
LAMMPS_image_vf = data_ref / "lammps/image_vf.data"

unordered_res = data_ref / "unordered_res.pdb"

GMS_ASYMOPT       = data_ref / "gms/c1opt.gms.gz"
GMS_SYMOPT        = data_ref / "gms/symopt.gms"
GMS_ASYMSURF      = data_ref / "gms/surf2wat.gms"

two_water_gro = data_ref / "two_water_gro.gro"
two_water_gro_multiframe = data_ref / "two_water_gro_multiframe.gro"
two_water_gro_nonames = data_ref / "two_water_gro_nonames.gro"
two_water_gro_widebox = data_ref / "two_water_gro_widebox.gro"

DLP_CONFIG = data_ref / "dlpoly/CONFIG"
DLP_CONFIG_order = data_ref / "dlpoly/CONFIG_order"
DLP_CONFIG_minimal = data_ref / "dlpoly/CONFIG_minimal"
DLP_HISTORY = data_ref / "dlpoly/HISTORY"
DLP_HISTORY_order = data_ref / "dlpoly/HISTORY_order"
DLP_HISTORY_minimal = data_ref / "dlpoly/HISTORY_minimal"
DLP_HISTORY_minimal_cell = data_ref / "dlpoly/HISTORY_minimal_cell"
DLP_HISTORY_classic = data_ref / "dlpoly/HISTORY_classic"

waterPSF = data_ref / 'watdyn.psf'
waterDCD = data_ref / 'watdyn.dcd'

rmsfArray = data_ref / 'adk_oplsaa_CA_rmsf.npy'

HoomdXMLdata = data_ref / 'C12x64.xml.bz2'

Make_Whole = data_ref / 'make_whole.gro'
fullerene = data_ref / 'fullerene.pdb.gz'

Plength = data_ref / 'plength.gro'
Martini_membrane_gro = data_ref / 'martini_dppc_chol_bilayer.gro'

# Contains one of each residue in 'nucleic' selections
NUCLsel = data_ref / 'nucl_res.pdb'

RANDOM_WALK = data_ref / 'xyz_random_walk.xtc'
RANDOM_WALK_TOPO = data_ref / 'RANDOM_WALK_TOPO.pdb'

MMTF = data_ref / '173D.mmtf'
MMTF_gz = data_ref / '5KIH.mmtf.gz'
MMTF_skinny = data_ref / '1ubq-less-optional.mmtf'
MMTF_skinny2 = data_ref / '3NJW-onlyrequired.mmtf'
MMTF_NOCRYST = data_ref / "6QYR.mmtf.gz"

ALIGN_BOUND = data_ref / 'analysis/align_bound.pdb.gz'
ALIGN_UNBOUND = data_ref / 'analysis/align_unbound.pdb.gz'

GSD = data_ref / 'example.gsd'
GSD_bonds = data_ref / 'example_bonds.gsd'
GSD_long = data_ref / 'example_longer.gsd'

DihedralArray = data_ref / 'adk_oplsaa_dihedral.npy'
DihedralsArray = data_ref / 'adk_oplsaa_dihedral_list.npy'
RamaArray = data_ref / 'adk_oplsaa_rama.npy'
GLYRamaArray = data_ref / 'adk_oplsaa_GLY_rama.npy'
JaninArray = data_ref / 'adk_oplsaa_janin.npy'
LYSJaninArray = data_ref / 'adk_oplsaa_LYS_janin.npy'
PDB_rama = data_ref / '19hc.pdb.gz'
PDB_janin = data_ref / '1a28.pdb.gz'

BATArray = data_ref / 'mol2_comments_header_bat.npy'

ITP = data_ref / 'gromacs_ala10.itp'
ITP_nomass = data_ref / 'itp_nomass.itp'
ITP_atomtypes = data_ref / 'atomtypes.itp'
ITP_charges = data_ref / 'atomtypes_charge.itp'
ITP_edited = data_ref / 'edited_itp.itp'
ITP_tip5p = data_ref / "tip5p.itp"
ITP_spce = data_ref / 'spce.itp'

GMX_TOP = data_ref / 'gromacs_ala10.top'
GMX_DIR = data_ref / 'gromacs/'
GMX_TOP_BAD = data_ref / 'bad_top.top'
ITP_no_endif = data_ref / 'no_endif_spc.itp'

NAMDBIN = data_ref / 'adk_open.coor'

SDF_molecule = data_ref / 'molecule.sdf'

PDB_elements = data_ref / 'elements.pdb'
PDB_charges = data_ref / 'charges.pdb'

PDBX = data_ref / "4x8u.pdbx"

SURFACE_PDB = data_ref / 'surface.pdb.bz2'
SURFACE_TRR = data_ref / 'surface.trr'


# This should be the last line: clean up namespace
del resources
