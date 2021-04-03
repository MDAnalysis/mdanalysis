# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
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


# TPR parser and tpr support module
# Copyright (c) 2011 Zhuyi Xue
# Released under the  GNU Public Licence, v2

"""
TPRParser settings
==================

Definition of constants.

The currently read file format versions are defined in
:data:`SUPPORTED_VERSIONS`.

"""

#: Gromacs TPR file format versions that can be read by the TPRParser.
SUPPORTED_VERSIONS = (58, 73, 83, 100, 103, 110, 112, 116, 119, 122)

# Some constants
STRLEN = 4096
BIG_STRLEN = 1048576
DIM = 3
NR_RBDIHS = 6  # <gromacs-5.1-dir>/src/gromacs/topology/idef.h
NR_CBTDIHS = 6  # <gromacs-5.1-dir>/src/gromacs/topology/idef.h
NR_FOURDIHS = 4  # <gromacs-5.1-dir>/src/gromacs/topology/idef.h
egcNR = 10  # include/types/topolog.h
TPX_TAG_RELEASE = "release"  # <gromacs-5.1-dir>/src/gromacs/fileio/tpxio.c
tpx_version = 103    # <gromacs-5.1-dir>/src/gromacs/fileio/tpxio.c
tpx_generation = 27  # <gromacs-5.1-dir>/src/gromacs/fileio/tpxio.c
tpxv_RestrictedBendingAndCombinedAngleTorsionPotentials = 98
tpxv_GenericInternalParameters = 117
tpxv_VSite2FD = 118
tpxv_AddSizeField = 119
tpxv_VSite1 = 121


#: Function types from ``<gromacs_dir>/include/types/idef.h``
(
    F_BONDS, F_G96BONDS, F_MORSE, F_CUBICBONDS,
    F_CONNBONDS, F_HARMONIC, F_FENEBONDS, F_TABBONDS,
    F_TABBONDSNC, F_RESTRBONDS, F_ANGLES, F_G96ANGLES, F_RESTRANGLES,
    F_LINEAR_ANGLES, F_CROSS_BOND_BONDS, F_CROSS_BOND_ANGLES, F_UREY_BRADLEY,
    F_QUARTIC_ANGLES, F_TABANGLES, F_PDIHS, F_RBDIHS, F_RESTRDIHS, F_CBTDIHS,
    F_FOURDIHS, F_IDIHS, F_PIDIHS, F_TABDIHS,
    F_CMAP, F_GB12, F_GB13, F_GB14,
    F_GBPOL, F_NPSOLVATION, F_LJ14, F_COUL14,
    F_LJC14_Q, F_LJC_PAIRS_NB, F_LJ, F_BHAM,
    F_LJ_LR, F_BHAM_LR, F_DISPCORR, F_COUL_SR,
    F_COUL_LR, F_RF_EXCL, F_COUL_RECIP, F_LJ_RECIP, F_DPD,
    F_POLARIZATION, F_WATER_POL, F_THOLE_POL, F_ANHARM_POL,
    F_POSRES, F_FBPOSRES, F_DISRES, F_DISRESVIOL, F_ORIRES,
    F_ORIRESDEV, F_ANGRES, F_ANGRESZ, F_DIHRES,
    F_DIHRESVIOL, F_CONSTR, F_CONSTRNC, F_SETTLE, F_VSITE1,
    F_VSITE2, F_VSITE2FD, F_VSITE3, F_VSITE3FD, F_VSITE3FAD,
    F_VSITE3OUT, F_VSITE4FD, F_VSITE4FDN, F_VSITEN,
    F_COM_PULL, F_DENSITYFITTING, F_EQM, F_EPOT, F_EKIN,
    F_ETOT, F_ECONSERVED, F_TEMP, F_VTEMP_NOLONGERUSED,
    F_PDISPCORR, F_PRES, F_DHDL_CON, F_DVDL,
    F_DKDL, F_DVDL_COUL, F_DVDL_VDW, F_DVDL_BONDED,
    F_DVDL_RESTRAINT, F_DVDL_TEMPERATURE, F_NRE) = list(range(95))

#: Function types from ``<gromacs_dir>/src/gmxlib/tpxio.c``
ftupd = [
    (20, F_CUBICBONDS), (20, F_CONNBONDS), (20, F_HARMONIC), (34, F_FENEBONDS),
    (43, F_TABBONDS), (43, F_TABBONDSNC), (70, F_RESTRBONDS),
    (tpxv_RestrictedBendingAndCombinedAngleTorsionPotentials, F_RESTRANGLES),
    (76, F_LINEAR_ANGLES), (30, F_CROSS_BOND_BONDS), (30, F_CROSS_BOND_ANGLES),
    (30, F_UREY_BRADLEY), (34, F_QUARTIC_ANGLES), (43, F_TABANGLES),
    (tpxv_RestrictedBendingAndCombinedAngleTorsionPotentials, F_RESTRDIHS),
    (tpxv_RestrictedBendingAndCombinedAngleTorsionPotentials, F_CBTDIHS),
    (26, F_FOURDIHS), (26, F_PIDIHS), (43, F_TABDIHS), (65, F_CMAP),
    (60, F_GB12), (61, F_GB13), (61, F_GB14), (72, F_GBPOL),
    (72, F_NPSOLVATION), (41, F_LJC14_Q), (41, F_LJC_PAIRS_NB),
    (32, F_BHAM_LR), (32, F_RF_EXCL), (32, F_COUL_RECIP), (93, F_LJ_RECIP),
    (46, F_DPD), (30, F_POLARIZATION), (36, F_THOLE_POL), (90, F_FBPOSRES),
    (22, F_DISRESVIOL), (22, F_ORIRES), (22, F_ORIRESDEV),
    (26, F_DIHRES), (26, F_DIHRESVIOL), (49, F_VSITE4FDN),
    (50, F_VSITEN), (46, F_COM_PULL), (20, F_EQM),
    (46, F_ECONSERVED), (69, F_VTEMP_NOLONGERUSED), (66, F_PDISPCORR),
    (54, F_DHDL_CON), (76, F_ANHARM_POL), (79, F_DVDL_COUL),
    (79, F_DVDL_VDW,), (79, F_DVDL_BONDED,), (79, F_DVDL_RESTRAINT),
    (79, F_DVDL_TEMPERATURE),
    (tpxv_GenericInternalParameters, F_DENSITYFITTING),
    (tpxv_VSite1, F_VSITE1),
    (tpxv_VSite2FD, F_VSITE2FD),
]

#: Interaction types from ``<gromacs_dir>/gmxlib/ifunc.c``
interaction_types = [
    ("BONDS", "Bond", 2),
    ("G96BONDS", "G96Bond", 2),
    ("MORSE", "Morse", 2),
    ("CUBICBONDS", "Cubic Bonds", 2),
    ("CONNBONDS", "Connect Bonds", 2),
    ("HARMONIC", "Harmonic Pot.", 2),
    ("FENEBONDS", "FENE Bonds", 2),
    ("TABBONDS", "Tab. Bonds", 2),
    ("TABBONDSNC", "Tab. Bonds NC", 2),
    ("RESTRAINTPOT", "Restraint Pot.", 2),
    ("ANGLES", "Angle", 3),
    ("G96ANGLES", "G96Angle", 3),
    ("RESTRANGLES", "Restricted Angles", 3),
    ("LINEAR_ANGLES", "Lin. Angle", 3),
    ("CROSS_BOND_BOND", "Bond-Cross", 3),
    ("CROSS_BOND_ANGLE", "BA-Cross", 3),
    ("UREY_BRADLEY", "U-B", 3),
    ("QANGLES", "Quartic Angles", 3),
    ("TABANGLES", "Tab. Angles", 3),
    ("PDIHS", "Proper Dih.", 4),
    ("RBDIHS", "Ryckaert-Bell.", 4),
    ("RESTRDIHS", "Restricted Dih.", 4),
    ("CBTDIHS", "CBT Dih.", 4),
    ("FOURDIHS", "Fourier Dih.", 4),
    ("IDIHS", "Improper Dih.", 4),
    ("PIDIHS", "Improper Dih.", 4),
    ("TABDIHS", "Tab. Dih.", 4),
    ("CMAP", "CMAP Dih.", 5),
    ("GB12", "GB 1-2 Pol.", 2),
    ("GB13", "GB 1-3 Pol.", 2),
    ("GB14", "GB 1-4 Pol.", 2),
    ("GBPOL", "GB Polarization", None),
    ("NPSOLVATION", "Nonpolar Sol.", None),
    ("LJ14", "LJ-14", 2),
    ("COUL14", "Coulomb-14", None),
    ("LJC14_Q", "LJC-14 q", 2),
    ("LJC_NB", "LJC Pairs NB", 2),
    ("LJ_SR", "LJ (SR)", 2),
    ("BHAM", "Buck.ham (SR)", 2),
    ("LJ_LR", "LJ (LR)", None),
    ("BHAM_LR", "Buck.ham (LR)", None),
    ("DISPCORR", "Disper. corr.", None),
    ("COUL_SR", "Coulomb (SR)", None),
    ("COUL_LR", "Coulomb (LR)", None),
    ("RF_EXCL", "RF excl.", None),
    ("COUL_RECIP", "Coul. recip.", None),
    ("LJ_RECIP", "LJ recip.", None),
    ("DPD", "DPD", None),
    ("POLARIZATION", "Polarization", 2),
    ("WATERPOL", "Water Pol.", 5),
    ("THOLE", "Thole Pol.", 4),
    ("ANHARM_POL", "Anharm. Pol.", 2),
    ("POSRES", "Position Rest.", 1),
    ("FBPOSRES", "Flat-bottom posres", 1),
    ("DISRES", "Dis. Rest.", 2),
    ("DISRESVIOL", "D.R.Viol. (nm)", None),
    ("ORIRES", "Orient. Rest.", 2),
    ("ORDEV", "Ori. R. RMSD", None),
    ("ANGRES", "Angle Rest.", 4),
    ("ANGRESZ", "Angle Rest. Z", 2),
    ("DIHRES", "Dih. Rest.", 4),
    ("DIHRESVIOL", "Dih. Rest. Viol.", None),
    ("CONSTR", "Constraint", 2),
    ("CONSTRNC", "Constr. No Conn.", 2),
    ("SETTLE", "Settle", 3),
    ("VSITE1", "Virtual site 1", 2),
    ("VSITE2", "Virtual site 2", 3),
    ("VSITE2FD", "Virtual site 2fd", 3),
    ("VSITE3", "Virtual site 3", 4),
    ("VSITE3FD", "Virtual site 3fd", 4),
    ("VSITE3FAD", "Virtual site 3fad", 4),
    ("VSITE3OUT", "Virtual site 3out", 4),
    ("VSITE4FD", "Virtual site 4fd", 5),
    ("VSITE4FDN", "Virtual site 4fdn", 5),
    ("VSITEN", "Virtual site N", 2),
    ("COM_PULL", "COM Pull En.", None),
    ("DENSITYFIT", "Density fitting", None),
    ("EQM", "Quantum En.", None),
    ("EPOT", "Potential", None),
    ("EKIN", "Kinetic En.", None),
    ("ETOT", "Total Energy", None),
    ("ECONS", "Conserved En.", None),
    ("TEMP", "Temperature", None),
    ("VTEMP", "Vir. Temp. (not used)", None),
    ("PDISPCORR", "Pres. DC", None),
    ("PRES", "Pressure", None),
    ("DH/DL_CON", "dH/dl constr.", None),
    ("DV/DL", "dVremain/dl", None),
    ("DK/DL", "dEkin/dl", None),
    ("DVC/DL", "dVcoul/dl", None),
    ("DVV/DL", "dVvdw/dl", None),
    ("DVB/DL", "dVbonded/dl", None),
    ("DVR/DL", "dVrestraint/dl", None),
    ("DVT/DL", "dVtemperature/dl", None)
]
