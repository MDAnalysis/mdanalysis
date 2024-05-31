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

Real MD simulation data are stored in the ``data/`` sub
directory. Use as ::

  from MDAnalysis.tests.datafiles import *

Note that the files are actually located in a separate package,
:mod:`MDAnalysisTestData` from where they are initially imported as ::

 from MDAnalysisTestData.datafiles import *
"""

__all__ = [
    "PSF", "DCD", "CRD",  # CHARMM (AdK example, DIMS trajectory from JMB 2009 paper)
    "DCD2",  # CHARMM (AdK example, DIMS trajectory from PLOS Comput Biol paper)
    "PSF_notop", "PSF_BAD",  # Same as PSF but no bonds etc, malformed version of previous
    "DCD_empty",
    "PSF_TRICLINIC", "DCD_TRICLINIC",  # CHARMM c36 new unitcell, NPT 125 TIP3P (box vectors, see Issue 187 for details)
    "PSF_NAMD", "PDB_NAMD",  # NAMD
    "PSF_NAMD_TRICLINIC", "DCD_NAMD_TRICLINIC",  # NAMD, triclinic unitcell (Issue 187)
    "PSF_NAMD_GBIS", "DCD_NAMD_GBIS",  # NAMD, implicit solvent, 100 steps,  #1819
    "PSF_nosegid",  # psf without a segid, Issue 121
    "PSF_cmap",  # ala3 PSF from ParmEd test files with cmap
    "PSF_inscode",  # PSF file with insertion codes
    "PDB_small",  # PDB
    "PDB_closed",
    "PDB_multiframe",
    "PDB_helix",
    "PDB_conect",
    "PDB_conect2TER",  # Conect record to a TER entry (Issue 936)
    "PDB_singleconect",  # Conect record with one entry (Issue 937)
    "PDB_icodes",  # stripped down version of 1osm, has icodes!
    "PDB_varying",  # varying occupancies and tempfactors
    "XPDB_small",
    "PDB_full",   # PDB 4E43 (full HEADER, TITLE, COMPND, REMARK, altloc)
    "ALIGN",  # Various way to align atom names in PDB files
    "RNA_PSF", "RNA_PDB",  # nucleic acid (PDB 1K5I in CHARMM36m)
    "INC_PDB",  # incomplete PDB file (Issue #396)
    # for testing cryst before/after model headers
    "PDB_cm", "PDB_cm_bz2", "PDB_cm_gz",
    "PDB_mc", "PDB_mc_bz2", "PDB_mc_gz",
    "PDB_chainidnewres",  # Issue 1110
    "PDB_sameresid_diffresname",  #Case where two residues share the same resid
    "PDB_chainidrepeat",  # Issue #1107
    "PDB", "GRO", "XTC", "TRR", "TPR", "GRO_velocity",  # Gromacs (AdK)
    "GRO_incomplete_vels",
    "COORDINATES_GRO_BZ2",
    "GRO_large",  #atom number truncation at > 100,000 particles, Issue 550
    "GRO_residwrap",  # resids wrapping because of 5 digit field (Issue #728)
    "GRO_residwrap_0base",  # corner case of #728 with resid=0 for first atom
    "GRO_sameresid_diffresname",  # Case where two residues share the same resid
    "PDB_xvf", "TPR_xvf", "TRR_xvf",  # Gromacs coords/veloc/forces (cobrotoxin, OPLS-AA, Gromacs 4.5.5 tpr)
    "H5MD_xvf",  # TPR_xvf + TRR_xvf converted to h5md format
    "XVG_BZ2",  # Compressed xvg file about cobrotoxin
    "PDB_xlserial",
    "TPR400", "TPR402", "TPR403", "TPR404", "TPR405", "TPR406", "TPR407",
    "TPR450", "TPR451", "TPR452", "TPR453", "TPR454", "TPR455", "TPR455Double",
    "TPR460", "TPR461", "TPR502", "TPR504", "TPR505", "TPR510", "TPR2016",
    "TPR2018", "TPR2019B3", "TPR2020B2", "TPR2020", "TPR2020Double",
    "TPR2021", "TPR2021Double", "TPR2022RC1", "TPR2023", "TPR2024",
    "TPR510_bonded", "TPR2016_bonded", "TPR2018_bonded", "TPR2019B3_bonded",
    "TPR2020B2_bonded", "TPR2020_bonded", "TPR2020_double_bonded",
    "TPR2021_bonded", "TPR2021_double_bonded", "TPR2022RC1_bonded",
    "TPR334_bonded", "TPR2023_bonded", "TPR2024_bonded",
    "TPR_EXTRA_2021", "TPR_EXTRA_2020", "TPR_EXTRA_2018",
    "TPR_EXTRA_2016", "TPR_EXTRA_407", "TPR_EXTRA_2022RC1",
    "TPR_EXTRA_2023", "TPR_EXTRA_2024",
    "PDB_sub_sol", "PDB_sub_dry",  # TRRReader sub selection
    "TRR_sub_sol",
    "XTC_sub_sol",
    "XYZ", "XYZ_psf", "XYZ_bz2",
    "XYZ_mini", "XYZ_five",  # 3 and 5 atoms xyzs for an easy topology
    "TXYZ", "ARC", "ARC_PBC",       # Tinker files
    "PRM",
    "PRM_chainid_bz2",
    "TRJ",
    "TRJ_bz2",  # Amber (no periodic box)
    "INPCRD",
    "PRMpbc", "TRJpbc_bz2",  # Amber (periodic box)
    "PRM7", "NCDFtruncoct",  # Amber (cpptrj test trajectory, see Issue 488)
    "PRM12", "TRJ12_bz2",  # Amber (v12 format, Issue 100)
    "PRMncdf", "TRJncdf", "NCDF",  # Amber (netcdf)
    "PFncdf_Top", "PFncdf_Trj",  # Amber ncdf with Positions and Forces
    "CPPTRAJ_TRAJ_TOP", "CPPTRAJ_TRAJ",  # Amber ncdf extracted from CPPTRAJ without time variable
    "PRMcs",  # Amber (format, Issue 1331)
    "PRMNCRST",  # Amber ncrst with positions/forces/velocities
    "PRM_NCBOX", "TRJ_NCBOX",  # Amber parm7 + nc w/ pos/forces/vels/box
    "PRMNEGATIVE",  # Amber negative ATOMIC_NUMBER (Issue 2306)
    "PRMErr1",  # Amber TOP files to check raised errors
    "PRMErr2",
    "PRMErr3",
    "PRMErr4",
    "PRMErr5",
    "PRM_UreyBradley",  # prmtop from ParmEd test files with Urey-Bradley angles
    "PRM7_ala2", "RST7_ala2",  # prmtop and rst files from ParmEd example files
    "PRM19SBOPC",  #  prmtop w/ ff19SB CMAP terms and OPC water (Issue #2449)
    "PQR",  # PQR v1
    "PQR_icodes",  # PQR v2 with icodes
    "PDBQT_input",  # PDBQT
    "PDBQT_querypdb",
    "PDBQT_tyrosol",
    "PDB_multipole",
    "FASTA",  # sequence alignment, Issue 112 + 113
    "HELANAL_BENDING_MATRIX",  # HELANAL test (from PSF+DCD (AdK) helix 8)
    "HELANAL_BENDING_MATRIX_SUBSET",  # As above, slice of frames 10 to 79
    "PDB_HOLE",  # gramicidin A
    "MULTIPDB_HOLE",  # gramicidin A, normal mode 7 from ElNemo
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
    "LAMMPSdata_many_bonds",
    "LAMMPSdata_deletedatoms",  # with deleted atoms
    "LAMMPSdata_triclinic",  # lammpsdata file to test triclinic dimension parsing, albite with most atoms deleted
    "LAMMPSdata_PairIJ",  # lammps datafile with a PairIJ Coeffs section
    "LAMMPSdata_additional_columns",  # structure for the additional column lammpstrj
    "LAMMPSDUMP",
    "LAMMPSDUMP_long",  # lammpsdump file with a few zeros sprinkled in the first column first frame
    "LAMMPSDUMP_allcoords",  # lammpsdump file with all coordinate conventions (x,xs,xu,xsu) present, from LAMMPS rdf example
    "LAMMPSDUMP_nocoords",  # lammpsdump file with no coordinates
    "LAMMPSDUMP_triclinic",  # lammpsdump file to test triclinic dimension parsing, albite with most atoms deleted
    "LAMMPSDUMP_image_vf",  # Lammps dump file with image flags, velocities, and forces.
    "LAMMPS_image_vf",  # Lammps data file to go with LAMMPSDUMP_image_vf
    "LAMMPSDUMP_chain1", # Lammps dump file with chain reader
    "LAMMPSDUMP_chain2", # Lammps dump file with chain reader
    "LAMMPS_chain", # Lammps data file with chain reader
    "LAMMPSDUMP_additional_columns",  # lammpsdump file with additional data (an additional charge column)
    "unordered_res",  # pdb file with resids non sequential
    "GMS_ASYMOPT",  # GAMESS C1  optimization
    "GMS_SYMOPT",   # GAMESS D4h optimization
    "GMS_ASYMSURF",  # GAMESS C1  surface
    "two_water_gro", "two_water_gro_nonames",  # for bond guessing, 2 water molecules, one with weird names
    "two_water_gro_multiframe",
    "two_water_gro_widebox",  # Issue #548
    "DLP_CONFIG", "DLP_CONFIG_order", "DLP_CONFIG_minimal",  # dl_poly 4 config file
    "DLP_HISTORY", "DLP_HISTORY_order", "DLP_HISTORY_minimal",  # dl_poly 4 history file
    "DLP_HISTORY_minimal_cell",  # dl_poly 4 history file with cell parameters
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
    "Martini_membrane_gro",  # for testing the leaflet finder
    "COORDINATES_XTC",
    "COORDINATES_TRR",
    "COORDINATES_TNG",
    "COORDINATES_H5MD",
    "COORDINATES_DCD",
    "COORDINATES_TOPOLOGY",
    "NUCLsel",
    "GRO_empty_atom", "GRO_missing_atomname",  # for testing GROParser exception raise
    "ENT",  #for testing ENT file extension
    "RANDOM_WALK",
    "RANDOM_WALK_TOPO",  # garbage topology to go along with XTC positions above
    "AUX_XVG", "XVG_BAD_NCOL",  #for testing .xvg auxiliary reader
    "AUX_XVG_LOWF", "AUX_XVG_HIGHF",
    "AUX_EDR", "AUX_EDR_TPR",
    "AUX_EDR_XTC", "AUX_EDR_RAW",
    "AUX_EDR_SINGLE_FRAME",  # for testing .edr auxiliary reader
    "MMTF", "MMTF_gz", 'MMTF_skinny',  # skinny - some optional fields stripped out
    "MMTF_skinny2",
    "ALIGN_BOUND",  # two component bound system
    "ALIGN_UNBOUND",  # two component unbound system
    "legacy_DCD_ADK_coords",  # frames 5 and 29 read in for adk_dims.dcd using legacy DCD reader
    "legacy_DCD_NAMD_coords",  # frame 0 read in for SiN_tric_namd.dcd using legacy DCD reader
    "legacy_DCD_c36_coords",  # frames 1 and 4 read in for tip125_tric_C36.dcd using legacy DCD reader
    "GSD", "GSD_bonds", "GSD_long",
    "TRC_PDB_VAC", "TRC_TRAJ1_VAC", "TRC_TRAJ2_VAC",  # 2x 3 frames of vacuum trajectory from GROMOS11 tutorial
    "TRC_CLUSTER_VAC",  # three frames without TIMESTEP and GENBOX block but with unsupported POSITION block
    "TRC_TRICLINIC_SOLV", "TRC_TRUNCOCT_VAC",
    "TRC_GENBOX_ORIGIN", "TRC_GENBOX_EULER",
    "TRC_EMPTY",  # Empty file containing only one space
    "TRC_PDB_SOLV", "TRC_TRAJ_SOLV",  # 2 frames of solvated trajectory from GROMOS11 tutorial
    "GRO_MEMPROT", "XTC_MEMPROT",  # YiiP transporter in POPE:POPG lipids with Na+, Cl-, Zn2+ dummy model without water
    "DihedralArray", "DihedralsArray",  # time series of single dihedral
    "RamaArray", "GLYRamaArray",  # time series of phi/psi angles
    "JaninArray", "LYSJaninArray",  # time series of chi1/chi2 angles
    "PDB_rama", "PDB_janin",  # for testing failures of Ramachandran and Janin classes
    "BATArray",  # time series of bond-angle-torsion coordinates array from Molecule_comments_header.mol2
    # DOS line endings
    "WIN_PDB_multiframe", "WIN_DLP_HISTORY", "WIN_TRJ", "WIN_LAMMPSDUMP", "WIN_ARC",
    "GRO_huge_box",  # for testing gro parser with hige box sizes
    "ITP",  # for GROMACS generated itps
    "ITP_nomass",  # for ATB generated itps
    "ITP_atomtypes",  # atom definitions to check atomtyes section parsing
    "ITP_charges",  # atom definitions to test custom particle charge parsing.
    "NAMDBIN",  # for NAMD generated binary file
    "ITP_edited",  # to check different directives are read properly
    "ITP_tip5p",  # tip5p water from opls-aa, edited with additional keywords
    "ITP_spce",  # spce water from gromos54a7, edited with additional keywords,
    "GMX_TOP",  # 2 ala10 chains + 3 spc water
    "GMX_DIR",  # GROMACS directory
    "GMX_TOP_BAD",  # file with an #include that doesn't exist
    "ITP_no_endif",  # file missing an #endif
    "PDB_CRYOEM_BOX",  # Issue 2599, Issue #2679, PR #2685
    "PDB_CHECK_RIGHTHAND_PA",  # for testing right handedness of principal_axes
    "MMTF_NOCRYST",  # File with meaningless CRYST1 record (Issue #2679, PR #2685)
    "FHIAIMS",  # to test FHIAIMS coordinate files
    "SDF_molecule",  # MDL SDFile for rdkit
    "PDBX",  # PDBxfile
    "PDB_elements",  # PDB file with elements
    "PDB_charges",  # PDB file with formal charges
    "SURFACE_PDB",  # 111 FCC lattice topology for NSGrid bug #2345
    "SURFACE_TRR",  # full precision coordinates for NSGrid bug #2345
]

from importlib import resources
import MDAnalysisTests.data

_data_ref = resources.files('MDAnalysisTests.data')

WIN_PDB_multiframe = (_data_ref / 'windows/WIN_nmr_neopetrosiamide.pdb').as_posix()
WIN_DLP_HISTORY = (_data_ref / 'windows/WIN_HISTORY').as_posix()
WIN_TRJ = (_data_ref / 'windows/WIN_ache.mdcrd').as_posix()
WIN_ARC = (_data_ref / 'windows/WIN_test.arc').as_posix()
WIN_LAMMPSDUMP = (_data_ref / 'windows/WIN_wat.lammpstrj').as_posix()

legacy_DCD_NAMD_coords = (_data_ref / 'legacy_DCD_NAMD_coords.npy').as_posix()
legacy_DCD_ADK_coords = (_data_ref / 'legacy_DCD_adk_coords.npy').as_posix()
legacy_DCD_c36_coords = (_data_ref / 'legacy_DCD_c36_coords.npy').as_posix()
AUX_XVG_LOWF = (_data_ref / 'test_lowf.xvg').as_posix()
AUX_XVG_HIGHF = (_data_ref / 'test_highf.xvg').as_posix()
XVG_BAD_NCOL = (_data_ref / 'bad_num_col.xvg').as_posix()
AUX_XVG = (_data_ref / 'test.xvg').as_posix()
AUX_EDR = (_data_ref / 'test.edr').as_posix()
AUX_EDR_RAW = (_data_ref / 'aux_edr_raw.txt').as_posix()
AUX_EDR_TPR = (_data_ref / 'aux_edr.tpr').as_posix()
AUX_EDR_XTC = (_data_ref / 'aux_edr.xtc').as_posix()
AUX_EDR_SINGLE_FRAME = (_data_ref / 'single_frame.edr').as_posix()
ENT = (_data_ref / 'testENT.ent').as_posix()
GRO_missing_atomname = (_data_ref / 'missing_atomname.gro').as_posix()
GRO_empty_atom = (_data_ref / 'empty_atom.gro').as_posix()
GRO_huge_box = (_data_ref / 'huge_box.gro').as_posix()

COORDINATES_GRO = (_data_ref / 'coordinates/test.gro').as_posix()
COORDINATES_GRO_INCOMPLETE_VELOCITY = (_data_ref / 'coordinates/test_incomplete_vel.gro').as_posix()
COORDINATES_GRO_BZ2 = (_data_ref / 'coordinates/test.gro.bz2').as_posix()
COORDINATES_XYZ = (_data_ref / 'coordinates/test.xyz').as_posix()
COORDINATES_XYZ_BZ2 = (_data_ref /  'coordinates/test.xyz.bz2').as_posix()
COORDINATES_XTC = (_data_ref / 'coordinates/test.xtc').as_posix()
COORDINATES_TRR = (_data_ref / 'coordinates/test.trr').as_posix()
COORDINATES_TNG = (_data_ref / 'coordinates/test.tng').as_posix()
COORDINATES_H5MD = (_data_ref / 'coordinates/test.h5md').as_posix()
COORDINATES_DCD = (_data_ref / 'coordinates/test.dcd').as_posix()
COORDINATES_TOPOLOGY = (_data_ref / 'coordinates/test_topology.pdb').as_posix()

PSF = (_data_ref / 'adk.psf').as_posix()
PSF_notop = (_data_ref / 'adk_notop.psf').as_posix()
PSF_BAD = (_data_ref / 'adk_notop_BAD.psf').as_posix()
DCD = (_data_ref / 'adk_dims.dcd').as_posix()
DCD_empty = (_data_ref / 'empty.dcd').as_posix()
CRD = (_data_ref / 'adk_open.crd').as_posix()
PSF_TRICLINIC = (_data_ref / 'tip125_tric_C36.psf').as_posix()
DCD_TRICLINIC = (_data_ref / 'tip125_tric_C36.dcd').as_posix()
DCD2 = (_data_ref / 'adk_dims2.dcd').as_posix()

PSF_NAMD = (_data_ref / 'namd_cgenff.psf').as_posix()
PDB_NAMD = (_data_ref / 'namd_cgenff.pdb').as_posix()
PDB_multipole = (_data_ref / 'water_methane_acetic-acid_ammonia.pdb').as_posix()
PSF_NAMD_TRICLINIC = (_data_ref / 'SiN_tric_namd.psf').as_posix()
DCD_NAMD_TRICLINIC = (_data_ref / 'SiN_tric_namd.dcd').as_posix()
PSF_NAMD_GBIS = (_data_ref / 'adk_closed_NAMD.psf').as_posix()
DCD_NAMD_GBIS = (_data_ref / 'adk_gbis_tmd-fast1_NAMD.dcd').as_posix()

PSF_nosegid = (_data_ref / 'nosegid.psf').as_posix()

PSF_cmap = (_data_ref / 'parmed_ala3.psf').as_posix()

PSF_inscode = (_data_ref / '1a2c_ins_code.psf').as_posix()

PDB_varying = (_data_ref / 'varying_occ_tmp.pdb').as_posix()
PDB_small = (_data_ref / 'adk_open.pdb').as_posix()
PDB_closed = (_data_ref / 'adk_closed.pdb').as_posix()

ALIGN = (_data_ref / 'align.pdb').as_posix()
RNA_PSF = (_data_ref / 'analysis/1k5i_c36.psf.gz').as_posix()
RNA_PDB = (_data_ref / 'analysis/1k5i_c36.pdb.gz').as_posix()
INC_PDB = (_data_ref / 'incomplete.pdb').as_posix()
PDB_cm = (_data_ref / 'cryst_then_model.pdb').as_posix()
PDB_cm_gz = (_data_ref / 'cryst_then_model.pdb.gz').as_posix()
PDB_cm_bz2 = (_data_ref / 'cryst_then_model.pdb.bz2').as_posix()
PDB_mc = (_data_ref / 'model_then_cryst.pdb').as_posix()
PDB_mc_gz = (_data_ref / 'model_then_cryst.pdb.gz').as_posix()
PDB_mc_bz2 = (_data_ref / 'model_then_cryst.pdb.bz2').as_posix()
PDB_chainidnewres = (_data_ref / 'chainIDnewres.pdb.gz').as_posix()
PDB_sameresid_diffresname = (_data_ref / 'sameresid_diffresname.pdb').as_posix()
PDB_chainidrepeat = (_data_ref / 'chainIDrepeat.pdb.gz').as_posix()
PDB_multiframe = (_data_ref / 'nmr_neopetrosiamide.pdb').as_posix()
PDB_helix = (_data_ref / 'A6PA6_alpha.pdb').as_posix()
PDB_conect = (_data_ref / 'conect_parsing.pdb').as_posix()
PDB_conect2TER = (_data_ref / 'CONECT2TER.pdb').as_posix()
PDB_singleconect = (_data_ref / 'SINGLECONECT.pdb').as_posix()
PDB_icodes = (_data_ref / '1osm.pdb.gz').as_posix()
PDB_CRYOEM_BOX = (_data_ref / '5a7u.pdb').as_posix()
PDB_CHECK_RIGHTHAND_PA = (_data_ref / '6msm.pdb.bz2').as_posix()
FHIAIMS = (_data_ref / 'fhiaims.in').as_posix()

GRO = (_data_ref / 'adk_oplsaa.gro').as_posix()
GRO_velocity = (_data_ref / 'sample_velocity_file.gro').as_posix()
GRO_incomplete_vels = (_data_ref / 'grovels.gro').as_posix()
GRO_large = (_data_ref / 'bigbox.gro.bz2').as_posix()
GRO_residwrap = (_data_ref / 'residwrap.gro').as_posix()
GRO_residwrap_0base = (_data_ref / 'residwrap_0base.gro').as_posix()
GRO_sameresid_diffresname = (_data_ref / 'sameresid_diffresname.gro').as_posix()
PDB = (_data_ref / 'adk_oplsaa.pdb').as_posix()
XTC = (_data_ref / 'adk_oplsaa.xtc').as_posix()
TRR = (_data_ref / 'adk_oplsaa.trr').as_posix()
TPR = (_data_ref / 'adk_oplsaa.tpr').as_posix()
PDB_sub_dry = (_data_ref / 'cobrotoxin_dry_neutral_0.pdb').as_posix()
TRR_sub_sol = (_data_ref / 'cobrotoxin.trr').as_posix()
XTC_sub_sol = (_data_ref / 'cobrotoxin.xtc').as_posix()
PDB_sub_sol = (_data_ref / 'cobrotoxin.pdb').as_posix()
PDB_xlserial = (_data_ref / 'xl_serial.pdb').as_posix()
GRO_MEMPROT = (_data_ref / 'analysis/YiiP_lipids.gro.gz').as_posix()
XTC_MEMPROT = (_data_ref / 'analysis/YiiP_lipids.xtc').as_posix()
XTC_multi_frame = (_data_ref / 'xtc_test_only_10_frame_10_atoms.xtc').as_posix()
TRR_multi_frame = (_data_ref / 'trr_test_only_10_frame_10_atoms.trr').as_posix()
TNG_traj = (_data_ref / 'argon_npt_compressed.tng').as_posix()
TNG_traj_gro = (_data_ref / 'argon_npt_compressed.gro.gz').as_posix()
TNG_traj_uneven_blocks = (_data_ref / 'argon_npt_compressed_uneven.tng').as_posix()
TNG_traj_vels_forces = (_data_ref / 'argon_npt_compressed_vels_forces.tng').as_posix()
PDB_xvf = (_data_ref / 'cobrotoxin.pdb').as_posix()
TPR_xvf = (_data_ref / 'cobrotoxin.tpr').as_posix()
TRR_xvf = (_data_ref / 'cobrotoxin.trr').as_posix()
H5MD_xvf = (_data_ref / 'cobrotoxin.h5md').as_posix()
XVG_BZ2 = (_data_ref / 'cobrotoxin_protein_forces.xvg.bz2').as_posix()

XPDB_small = (_data_ref / '5digitResid.pdb').as_posix()
# number is the gromacs version
TPR400 = (_data_ref / 'tprs/2lyz_gmx_4.0.tpr').as_posix()
TPR402 = (_data_ref / 'tprs/2lyz_gmx_4.0.2.tpr').as_posix()
TPR403 = (_data_ref / 'tprs/2lyz_gmx_4.0.3.tpr').as_posix()
TPR404 = (_data_ref / 'tprs/2lyz_gmx_4.0.4.tpr').as_posix()
TPR405 = (_data_ref / 'tprs/2lyz_gmx_4.0.5.tpr').as_posix()
TPR406 = (_data_ref / 'tprs/2lyz_gmx_4.0.6.tpr').as_posix()
TPR407 = (_data_ref / 'tprs/2lyz_gmx_4.0.7.tpr').as_posix()
TPR450 = (_data_ref / 'tprs/2lyz_gmx_4.5.tpr').as_posix()
TPR451 = (_data_ref / 'tprs/2lyz_gmx_4.5.1.tpr').as_posix()
TPR452 = (_data_ref / 'tprs/2lyz_gmx_4.5.2.tpr').as_posix()
TPR453 = (_data_ref / 'tprs/2lyz_gmx_4.5.3.tpr').as_posix()
TPR454 = (_data_ref / 'tprs/2lyz_gmx_4.5.4.tpr').as_posix()
TPR455 = (_data_ref / 'tprs/2lyz_gmx_4.5.5.tpr').as_posix()
TPR502 = (_data_ref / 'tprs/2lyz_gmx_5.0.2.tpr').as_posix()
TPR504 = (_data_ref / 'tprs/2lyz_gmx_5.0.4.tpr').as_posix()
TPR505 = (_data_ref / 'tprs/2lyz_gmx_5.0.5.tpr').as_posix()
TPR510 = (_data_ref / 'tprs/2lyz_gmx_5.1.tpr').as_posix()
TPR2016 = (_data_ref / 'tprs/2lyz_gmx_2016.tpr').as_posix()
TPR2018 = (_data_ref / 'tprs/2lyz_gmx_2018.tpr').as_posix()
TPR2019B3 = (_data_ref / 'tprs/2lyz_gmx_2019-beta3.tpr').as_posix()
TPR2020B2 = (_data_ref / 'tprs/2lyz_gmx_2020-beta2.tpr').as_posix()
TPR2020 = (_data_ref / 'tprs/2lyz_gmx_2020.tpr').as_posix()
TPR2021 = (_data_ref / 'tprs/2lyz_gmx_2021.tpr').as_posix()
TPR2022RC1 = (_data_ref / 'tprs/2lyz_gmx_2022-rc1.tpr').as_posix()
TPR2023 = (_data_ref / 'tprs/2lyz_gmx_2023.tpr').as_posix()
TPR2024 = (_data_ref / 'tprs/2lyz_gmx_2024.tpr').as_posix()
# double precision
TPR455Double = (_data_ref / 'tprs/drew_gmx_4.5.5.double.tpr').as_posix()
TPR460 = (_data_ref / 'tprs/ab42_gmx_4.6.tpr').as_posix()
TPR461 = (_data_ref / 'tprs/ab42_gmx_4.6.1.tpr').as_posix()
TPR2020Double = (_data_ref / 'tprs/2lyz_gmx_2020_double.tpr').as_posix()
TPR2021Double = (_data_ref / 'tprs/2lyz_gmx_2021_double.tpr').as_posix()
# all bonded interactions
TPR334_bonded = (_data_ref / 'tprs/all_bonded/dummy_3.3.4.tpr').as_posix()
TPR510_bonded = (_data_ref / 'tprs/all_bonded/dummy_5.1.tpr').as_posix()
TPR2016_bonded = (_data_ref / 'tprs/all_bonded/dummy_2016.tpr').as_posix()
TPR2018_bonded = (_data_ref / 'tprs/all_bonded/dummy_2018.tpr').as_posix()
TPR2019B3_bonded = (_data_ref / 'tprs/all_bonded/dummy_2019-beta3.tpr').as_posix()
TPR2020B2_bonded = (_data_ref / 'tprs/all_bonded/dummy_2020-beta2.tpr').as_posix()
TPR2020_bonded = (_data_ref / 'tprs/all_bonded/dummy_2020.tpr').as_posix()
TPR2020_double_bonded = (_data_ref / 'tprs/all_bonded/dummy_2020_double.tpr').as_posix()
TPR2021_bonded = (_data_ref / 'tprs/all_bonded/dummy_2021.tpr').as_posix()
TPR2021_double_bonded = (_data_ref / 'tprs/all_bonded/dummy_2021_double.tpr').as_posix()
TPR2022RC1_bonded = (_data_ref / 'tprs/all_bonded/dummy_2022-rc1.tpr').as_posix()
TPR2023_bonded = (_data_ref / 'tprs/all_bonded/dummy_2023.tpr').as_posix()
TPR2024_bonded = (_data_ref / 'tprs/all_bonded/dummy_2024.tpr').as_posix()
# all interactions
TPR_EXTRA_2024 = (_data_ref / 'tprs/virtual_sites/extra-interactions-2024.tpr').as_posix()
TPR_EXTRA_2023 = (_data_ref / 'tprs/virtual_sites/extra-interactions-2023.tpr').as_posix()
TPR_EXTRA_2022RC1 = (_data_ref / 'tprs/virtual_sites/extra-interactions-2022-rc1.tpr').as_posix()
TPR_EXTRA_2021 = (_data_ref / 'tprs/virtual_sites/extra-interactions-2021.tpr').as_posix()
TPR_EXTRA_2020 = (_data_ref / 'tprs/virtual_sites/extra-interactions-2020.tpr').as_posix()
TPR_EXTRA_2018 = (_data_ref / 'tprs/virtual_sites/extra-interactions-2018.tpr').as_posix()
TPR_EXTRA_2016 = (_data_ref / 'tprs/virtual_sites/extra-interactions-2016.3.tpr').as_posix()
TPR_EXTRA_407 = (_data_ref / 'tprs/virtual_sites/extra-interactions-4.0.7.tpr').as_posix()

XYZ_psf = (_data_ref / '2r9r-1b.psf').as_posix()
XYZ_bz2 = (_data_ref / '2r9r-1b.xyz.bz2').as_posix()
XYZ = (_data_ref / '2r9r-1b.xyz').as_posix()
XYZ_mini = (_data_ref / 'mini.xyz').as_posix()
XYZ_five = (_data_ref / 'five.xyz').as_posix()
TXYZ = (_data_ref / 'coordinates/test.txyz').as_posix()
ARC = (_data_ref / 'coordinates/test.arc').as_posix()
ARC_PBC = (_data_ref / 'coordinates/new_hexane.arc').as_posix()

PRM = (_data_ref / 'Amber/ache.prmtop').as_posix()
TRJ = (_data_ref / 'Amber/ache.mdcrd').as_posix()
INPCRD = (_data_ref / 'Amber/test.inpcrd').as_posix()
TRJ_bz2 = (_data_ref / 'Amber/ache.mdcrd.bz2').as_posix()
PFncdf_Top = (_data_ref / 'Amber/posfor.top').as_posix()
PFncdf_Trj = (_data_ref / 'Amber/posfor.ncdf').as_posix()
PRM_chainid_bz2 = (_data_ref / "Amber/ache_chainid.prmtop.bz2").as_posix()

CPPTRAJ_TRAJ_TOP = (_data_ref / 'Amber/cpptraj_traj.prmtop').as_posix()
CPPTRAJ_TRAJ = (_data_ref / 'Amber/cpptraj_traj.nc').as_posix()

PRMpbc = (_data_ref / 'Amber/capped-ala.prmtop').as_posix()
TRJpbc_bz2 = (_data_ref / 'Amber/capped-ala.mdcrd.bz2').as_posix()

PRMncdf = (_data_ref / 'Amber/bala.prmtop').as_posix()
TRJncdf = (_data_ref / 'Amber/bala.trj').as_posix()
NCDF = (_data_ref / 'Amber/bala.ncdf').as_posix()

PRM12 = (_data_ref / 'Amber/anti.top').as_posix()
TRJ12_bz2 = (_data_ref / 'Amber/anti_md1.mdcrd.bz2').as_posix()

PRM7 =  (_data_ref / 'Amber/tz2.truncoct.parm7.bz2').as_posix()
NCDFtruncoct =  (_data_ref / 'Amber/tz2.truncoct.nc').as_posix()

PRMcs = (_data_ref / 'Amber/chitosan.prmtop').as_posix()

PRMNCRST = (_data_ref / 'Amber/ace_mbondi3.parm7').as_posix()

PRM_NCBOX = (_data_ref / 'Amber/ace_tip3p.parm7').as_posix()
TRJ_NCBOX = (_data_ref / 'Amber/ace_tip3p.nc').as_posix()

PRMNEGATIVE = (_data_ref / 'Amber/ace_mbondi3.negative.parm7').as_posix()

PRMErr1 = (_data_ref / "Amber/ace_mbondi3.error1.parm7").as_posix()
PRMErr2 = (_data_ref / "Amber/ace_mbondi3.error2.parm7").as_posix()
PRMErr3 = (_data_ref / "Amber/ace_mbondi3.error3.parm7").as_posix()
PRMErr4 = (_data_ref / "Amber/ace_mbondi3.error4.parm7").as_posix()
PRMErr5 = (_data_ref / "Amber/ache_chainid.error5.prmtop.bz2").as_posix()

PRM_UreyBradley = (_data_ref / 'Amber/parmed_fad.prmtop').as_posix()
PRM7_ala2 = (_data_ref / 'Amber/parmed_ala2_solv.parm7').as_posix()
RST7_ala2 = (_data_ref / 'Amber/parmed_ala2_solv.rst7').as_posix()

PRM19SBOPC = (_data_ref / 'Amber/ala.ff19SB.OPC.parm7.bz2').as_posix()

PQR = (_data_ref / 'adk_open.pqr').as_posix()
PQR_icodes = (_data_ref / '1A2C.pqr').as_posix()

PDBQT_input = (_data_ref / "pdbqt_inputpdbqt.pdbqt").as_posix()
PDBQT_querypdb = (_data_ref / "pdbqt_querypdb.pdb").as_posix()
PDBQT_tyrosol = (_data_ref / "tyrosol.pdbqt.bz2").as_posix()

FASTA = (_data_ref / 'test.fasta').as_posix()
HELANAL_BENDING_MATRIX = (_data_ref / 'helanal_bending_matrix_AdK_DIMS_H8.dat').as_posix()
HELANAL_BENDING_MATRIX_SUBSET = (_data_ref / 'helanal_bending_matrix_AdK_DIMS_H8_frames10to79.dat').as_posix()

PDB_HOLE = (_data_ref / '1grm_single.pdb').as_posix()
MULTIPDB_HOLE = (_data_ref / '1grm_elNemo_mode7.pdb.bz2').as_posix()

DMS = (_data_ref / 'adk_closed.dms').as_posix()
DMS_DOMAINS = (_data_ref / 'adk_closed_domains.dms').as_posix()
DMS_NO_SEGID = (_data_ref / 'adk_closed_no_segid.dms').as_posix()

CONECT = (_data_ref / '1hvr.pdb').as_posix()
CONECT_ERROR = (_data_ref / 'conect_error.pdb').as_posix()

TRZ = (_data_ref / 'trzfile.trz').as_posix()
TRZ_psf = (_data_ref / 'trz_psf.psf').as_posix()

TRIC = (_data_ref / 'dppc_vesicle_hg.gro').as_posix()

PDB_full = (_data_ref / "4E43.pdb").as_posix()

merge_protein = (_data_ref / "merge/2zmm/protein.pdb").as_posix()
merge_ligand = (_data_ref / "merge/2zmm/ligand.pdb").as_posix()
merge_water = (_data_ref / "merge/2zmm/water.pdb").as_posix()

mol2_molecules = (_data_ref / "mol2/Molecules.mol2").as_posix()
mol2_molecule = (_data_ref / "mol2/Molecule.mol2").as_posix()
mol2_ligand = (_data_ref / "mol2/Ligand.mol2").as_posix()
mol2_broken_molecule = (_data_ref / "mol2/BrokenMolecule.mol2").as_posix()
mol2_comments_header = (_data_ref / "mol2/Molecule_comments_header.mol2").as_posix()
# MOL2 file without substructure field
mol2_zinc = (_data_ref / "mol2/zinc_856218.mol2").as_posix()
# MOL2 file without bonds
mol2_sodium_ion = (_data_ref / "mol2/sodium_ion.mol2").as_posix()

capping_input = (_data_ref / "capping/aaqaa.gro").as_posix()
capping_output = (_data_ref / "capping/maestro_aaqaa_capped.pdb").as_posix()
capping_ace = (_data_ref / "capping/ace.pdb").as_posix()
capping_nma = (_data_ref / "capping/nma.pdb").as_posix()

contacts_villin_folded = (_data_ref / "contacts/villin_folded.gro.bz2").as_posix()
contacts_villin_unfolded = (_data_ref / "contacts/villin_unfolded.gro.bz2").as_posix()
contacts_file = (_data_ref / "contacts/2F4K_qlist5_remap.dat").as_posix()

trz4data = (_data_ref / "lammps/datatest.trz").as_posix()
LAMMPSdata = (_data_ref / "lammps/datatest.data").as_posix()
LAMMPSdata_mini = (_data_ref / "lammps/mini.data").as_posix()
LAMMPSdata2 = (_data_ref / "lammps/ifabp_apo_100mM.data.bz2").as_posix()
LAMMPSdcd2 = (_data_ref / "lammps/ifabp_apo_100mM.dcd").as_posix()
LAMMPScnt = (_data_ref / "lammps/cnt-hexagonal-class1.data").as_posix()
LAMMPScnt2 = (_data_ref / "lammps/cnt-hexagonal-class1.data2").as_posix()
LAMMPShyd = (_data_ref / "lammps/hydrogen-class1.data").as_posix()
LAMMPShyd2 = (_data_ref / "lammps/hydrogen-class1.data2").as_posix()
LAMMPSdata_deletedatoms = (_data_ref / 'lammps/deletedatoms.data').as_posix()
LAMMPSdata_triclinic = (_data_ref / "lammps/albite_triclinic.data").as_posix()
LAMMPSdata_PairIJ = (_data_ref / "lammps/pairij_coeffs.data.bz2").as_posix()
LAMMPSDUMP = (_data_ref / "lammps/wat.lammpstrj.bz2").as_posix()
LAMMPSDUMP_long = (_data_ref / "lammps/wat.lammpstrj_long.bz2").as_posix()
LAMMPSDUMP_allcoords = (_data_ref / "lammps/spce_all_coords.lammpstrj.bz2").as_posix()
LAMMPSDUMP_nocoords = (_data_ref / "lammps/spce_no_coords.lammpstrj.bz2").as_posix()
LAMMPSDUMP_triclinic = (_data_ref / "lammps/albite_triclinic.dump").as_posix()
LAMMPSDUMP_image_vf = (_data_ref / "lammps/image_vf.lammpstrj").as_posix()
LAMMPS_image_vf = (_data_ref / "lammps/image_vf.data").as_posix()
LAMMPSDUMP_chain1 = (_data_ref / "lammps/chain_dump_1.lammpstrj").as_posix()
LAMMPSDUMP_chain2 = (_data_ref / "lammps/chain_dump_2.lammpstrj").as_posix()
LAMMPS_chain = (_data_ref / "lammps/chain_initial.data").as_posix()
LAMMPSdata_many_bonds = (_data_ref / "lammps/a_lot_of_bond_types.data").as_posix()
LAMMPSdata_additional_columns = (_data_ref / "lammps/additional_columns.data").as_posix()
LAMMPSDUMP_additional_columns = (_data_ref / "lammps/additional_columns.lammpstrj").as_posix()

unordered_res = (_data_ref / "unordered_res.pdb").as_posix()

GMS_ASYMOPT       = (_data_ref / "gms/c1opt.gms.gz").as_posix()
GMS_SYMOPT        = (_data_ref / "gms/symopt.gms").as_posix()
GMS_ASYMSURF      = (_data_ref / "gms/surf2wat.gms").as_posix()

two_water_gro = (_data_ref / "two_water_gro.gro").as_posix()
two_water_gro_multiframe = (_data_ref / "two_water_gro_multiframe.gro").as_posix()
two_water_gro_nonames = (_data_ref / "two_water_gro_nonames.gro").as_posix()
two_water_gro_widebox = (_data_ref / "two_water_gro_widebox.gro").as_posix()

DLP_CONFIG = (_data_ref / "dlpoly/CONFIG").as_posix()
DLP_CONFIG_order = (_data_ref / "dlpoly/CONFIG_order").as_posix()
DLP_CONFIG_minimal = (_data_ref / "dlpoly/CONFIG_minimal").as_posix()
DLP_HISTORY = (_data_ref / "dlpoly/HISTORY").as_posix()
DLP_HISTORY_order = (_data_ref / "dlpoly/HISTORY_order").as_posix()
DLP_HISTORY_minimal = (_data_ref / "dlpoly/HISTORY_minimal").as_posix()
DLP_HISTORY_minimal_cell = (_data_ref / "dlpoly/HISTORY_minimal_cell").as_posix()
DLP_HISTORY_classic = (_data_ref / "dlpoly/HISTORY_classic").as_posix()

waterPSF = (_data_ref / 'watdyn.psf').as_posix()
waterDCD = (_data_ref / 'watdyn.dcd').as_posix()

rmsfArray = (_data_ref / 'adk_oplsaa_CA_rmsf.npy').as_posix()

HoomdXMLdata = (_data_ref / 'C12x64.xml.bz2').as_posix()

Make_Whole = (_data_ref / 'make_whole.gro').as_posix()
fullerene = (_data_ref / 'fullerene.pdb.gz').as_posix()

Plength = (_data_ref / 'plength.gro').as_posix()
Martini_membrane_gro = (_data_ref / 'martini_dppc_chol_bilayer.gro').as_posix()

# Contains one of each residue in 'nucleic' selections
NUCLsel = (_data_ref / 'nucl_res.pdb').as_posix()

RANDOM_WALK = (_data_ref / 'xyz_random_walk.xtc').as_posix()
RANDOM_WALK_TOPO = (_data_ref / 'RANDOM_WALK_TOPO.pdb').as_posix()

MMTF = (_data_ref / '173D.mmtf').as_posix()
MMTF_gz = (_data_ref / '5KIH.mmtf.gz').as_posix()
MMTF_skinny = (_data_ref / '1ubq-less-optional.mmtf').as_posix()
MMTF_skinny2 = (_data_ref / '3NJW-onlyrequired.mmtf').as_posix()
MMTF_NOCRYST = (_data_ref / "6QYR.mmtf.gz").as_posix()

ALIGN_BOUND = (_data_ref / 'analysis/align_bound.pdb.gz').as_posix()
ALIGN_UNBOUND = (_data_ref / 'analysis/align_unbound.pdb.gz').as_posix()

GSD = (_data_ref / 'example.gsd').as_posix()
GSD_bonds = (_data_ref / 'example_bonds.gsd').as_posix()
GSD_long = (_data_ref / 'example_longer.gsd').as_posix()

TRC_PDB_VAC = (_data_ref / 'gromos11/gromos11_traj_vac.pdb.gz').as_posix()
TRC_TRAJ1_VAC = (_data_ref / 'gromos11/gromos11_traj_vac_1.trc.gz').as_posix()
TRC_TRAJ2_VAC = (_data_ref / 'gromos11/gromos11_traj_vac_2.trc.gz').as_posix()
TRC_PDB_SOLV = (_data_ref / 'gromos11/gromos11_traj_solv.pdb.gz').as_posix()
TRC_TRAJ_SOLV = (_data_ref / 'gromos11/gromos11_traj_solv.trc.gz').as_posix()
TRC_CLUSTER_VAC = (_data_ref / 'gromos11/gromos11_cluster_vac.trj.gz').as_posix()
TRC_TRICLINIC_SOLV = (_data_ref / 'gromos11/gromos11_triclinic_solv.trc.gz').as_posix()
TRC_TRUNCOCT_VAC = (_data_ref / 'gromos11/gromos11_truncOcta_vac.trc.gz').as_posix()
TRC_GENBOX_ORIGIN = (_data_ref / 'gromos11/gromos11_genbox_origin.trc.gz').as_posix()
TRC_GENBOX_EULER = (_data_ref / 'gromos11/gromos11_genbox_euler.trc.gz').as_posix()
TRC_EMPTY = (_data_ref / 'gromos11/gromos11_empty.trc').as_posix()

DihedralArray = (_data_ref / 'adk_oplsaa_dihedral.npy').as_posix()
DihedralsArray = (_data_ref / 'adk_oplsaa_dihedral_list.npy').as_posix()
RamaArray = (_data_ref / 'adk_oplsaa_rama.npy').as_posix()
GLYRamaArray = (_data_ref / 'adk_oplsaa_GLY_rama.npy').as_posix()
JaninArray = (_data_ref / 'adk_oplsaa_janin.npy').as_posix()
LYSJaninArray = (_data_ref / 'adk_oplsaa_LYS_janin.npy').as_posix()
PDB_rama = (_data_ref / '19hc.pdb.gz').as_posix()
PDB_janin = (_data_ref / '1a28.pdb.gz').as_posix()

BATArray = (_data_ref / 'mol2_comments_header_bat.npy').as_posix()

ITP = (_data_ref / 'gromacs_ala10.itp').as_posix()
ITP_nomass = (_data_ref / 'itp_nomass.itp').as_posix()
ITP_atomtypes = (_data_ref / 'atomtypes.itp').as_posix()
ITP_charges = (_data_ref / 'atomtypes_charge.itp').as_posix()
ITP_edited = (_data_ref / 'edited_itp.itp').as_posix()
ITP_tip5p = (_data_ref / "tip5p.itp").as_posix()
ITP_spce = (_data_ref / 'spce.itp').as_posix()

GMX_TOP = (_data_ref / 'gromacs_ala10.top').as_posix()
GMX_DIR = (_data_ref / 'gromacs/').as_posix()
GMX_TOP_BAD = (_data_ref / 'bad_top.top').as_posix()
ITP_no_endif = (_data_ref / 'no_endif_spc.itp').as_posix()

NAMDBIN = (_data_ref / 'adk_open.coor').as_posix()

SDF_molecule = (_data_ref / 'molecule.sdf').as_posix()

PDB_elements = (_data_ref / 'elements.pdb').as_posix()
PDB_charges = (_data_ref / 'charges.pdb').as_posix()

PDBX = (_data_ref / "4x8u.pdbx").as_posix()

SURFACE_PDB = (_data_ref / 'surface.pdb.bz2').as_posix()
SURFACE_TRR = (_data_ref / 'surface.trr').as_posix()

# This should be the last line: clean up namespace
del resources
