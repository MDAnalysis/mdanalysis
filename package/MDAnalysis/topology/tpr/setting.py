# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://mdanalysis.googlecode.com
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


# TPR parser and tpr support module
# Copyright (c) 2011 Zhuyi Xue
# Released under the  GNU Public Licence, v2

"""
TPRParser settings
==================

Definition of constants and reading of the TPR header with
:func:`do_inputrec`.

The currently read file format versions are defined in
:data:`SUPPORTED_VERSIONS`.

"""

import utils as U

#: Gromacs TPR file format versions that can be read by the TPRParser.
SUPPORTED_VERSIONS = (58, 73, 83)

# Some constants
STRLEN = 4096
BIG_STRLEN = 1048576
DIM = 3
NR_RBDIHS = 6  # include/types/idef.h
egcNR = 10  # include/types/topolog.h
TPX_TAG_RELEASE = "release"  # <gromacs-4.6.1-dir>/src/gmxlib/tpxio.c
tpx_version = 83  # <gromacs-4.6.1-dir>/src/gmxlib/tpxio.c
tpx_generation = 24  # <gromacs-4.6.1-dir>/src/gmxlib/tpxio.c


#: Function types from ``<gromacs_dir>/include/types/idef.h``
(
    F_BONDS, F_G96BONDS, F_MORSE, F_CUBICBONDS,
    F_CONNBONDS, F_HARMONIC, F_FENEBONDS, F_TABBONDS,
    F_TABBONDSNC, F_RESTRBONDS, F_ANGLES, F_G96ANGLES,
    F_LINEAR_ANGLES, F_CROSS_BOND_BONDS, F_CROSS_BOND_ANGLES, F_UREY_BRADLEY,
    F_QUARTIC_ANGLES, F_TABANGLES, F_PDIHS, F_RBDIHS,
    F_FOURDIHS, F_IDIHS, F_PIDIHS, F_TABDIHS,
    F_CMAP, F_GB12, F_GB13, F_GB14,
    F_GBPOL, F_NPSOLVATION, F_LJ14, F_COUL14,
    F_LJC14_Q, F_LJC_PAIRS_NB, F_LJ, F_BHAM,
    F_LJ_LR, F_BHAM_LR, F_DISPCORR, F_COUL_SR,
    F_COUL_LR, F_RF_EXCL, F_COUL_RECIP, F_DPD,
    F_POLARIZATION, F_WATER_POL, F_THOLE_POL, F_ANHARM_POL,
    F_POSRES, F_DISRES, F_DISRESVIOL, F_ORIRES,
    F_ORIRESDEV, F_ANGRES, F_ANGRESZ, F_DIHRES,
    F_DIHRESVIOL, F_CONSTR, F_CONSTRNC, F_SETTLE,
    F_VSITE2, F_VSITE3, F_VSITE3FD, F_VSITE3FAD,
    F_VSITE3OUT, F_VSITE4FD, F_VSITE4FDN, F_VSITEN,
    F_COM_PULL, F_EQM, F_EPOT, F_EKIN,
    F_ETOT, F_ECONSERVED, F_TEMP, F_VTEMP_NOLONGERUSED,
    F_PDISPCORR, F_PRES, F_DHDL_CON, F_DVDL,
    F_DKDL, F_DVDL_COUL, F_DVDL_VDW, F_DVDL_BONDED,
    F_DVDL_RESTRAINT, F_DVDL_TEMPERATURE, F_NRE) = range(87)

#: Function types from ``<gromacs_dir>/src/gmxlib/tpxio.c``
ftupd = [
    (20, F_CUBICBONDS), (20, F_CONNBONDS), (20, F_HARMONIC),
    (34, F_FENEBONDS), (43, F_TABBONDS), (43, F_TABBONDSNC),
    (70, F_RESTRBONDS), (76, F_LINEAR_ANGLES), (30, F_CROSS_BOND_BONDS),
    (30, F_CROSS_BOND_ANGLES), (30, F_UREY_BRADLEY), (34, F_QUARTIC_ANGLES),
    (43, F_TABANGLES), (26, F_FOURDIHS), (26, F_PIDIHS),
    (43, F_TABDIHS), (65, F_CMAP), (60, F_GB12),
    (61, F_GB13), (61, F_GB14), (72, F_GBPOL),
    (72, F_NPSOLVATION), (41, F_LJC14_Q), (41, F_LJC_PAIRS_NB),
    (32, F_BHAM_LR), (32, F_RF_EXCL), (32, F_COUL_RECIP),
    (46, F_DPD), (30, F_POLARIZATION), (36, F_THOLE_POL),
    (22, F_DISRESVIOL), (22, F_ORIRES), (22, F_ORIRESDEV),
    (26, F_DIHRES), (26, F_DIHRESVIOL), (49, F_VSITE4FDN),
    (50, F_VSITEN), (46, F_COM_PULL), (20, F_EQM),
    (46, F_ECONSERVED), (69, F_VTEMP_NOLONGERUSED), (66, F_PDISPCORR),
    (54, F_DHDL_CON), (76, F_ANHARM_POL), (79, F_DVDL_COUL),
    (79, F_DVDL_VDW,), (79, F_DVDL_BONDED,), (79, F_DVDL_RESTRAINT),
    (79, F_DVDL_TEMPERATURE), (54, F_DHDL_CON)
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
    ("LINEAR_ANGLES", "Lin. Angle", 3),
    ("CROSS_BOND_BOND", "Bond-Cross", 3),
    ("CROSS_BOND_ANGLE", "BA-Cross", 3),
    ("UREY_BRADLEY", "U-B", 3),
    ("QANGLES", "Quartic Angles", 3),
    ("TABANGLES", "Tab. Angles", 3),
    ("PDIHS", "Proper Dih.", 4),
    ("RBDIHS", "Ryckaert-Bell.", 4),
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
    ("DPD", "DPD", None),
    ("POLARIZATION", "Polarization", 2),
    ("WATERPOL", "Water Pol.", 5),
    ("THOLE", "Thole Pol.", 4),
    ("ANHARM_POL", "Anharm. Pol.", 2),
    ("POSRES", "Position Rest.", 1),
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
    ("VSITE2", "Virtual site 2", 3),
    ("VSITE3", "Virtual site 3", 4),
    ("VSITE3FD", "Virtual site 3fd", 4),
    ("VSITE3FAD", "Virtual site 3fad", 4),
    ("VSITE3OUT", "Virtual site 3out", 4),
    ("VSITE4FD", "Virtual site 4fd", 5),
    ("VSITE4FDN", "Virtual site 4fdn", 5),
    ("VSITEN", "Virtual site N", 2),
    ("COM_PULL", "COM Pull En.", None),
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


def do_inputrec(data):
    data.unpack_int()  # ir_eI
    data.unpack_int()  # ir_nsteps = idum
    data.unpack_int()  # ir_init_step = idum =

    data.unpack_int()  # simulation_part
    # not relevant here
    # ir_nstcalcenergy = 1

    data.unpack_int()  # ir_ns_type
    data.unpack_int()  # ir_nslist
    data.unpack_int()  # ir_ndelta

    data.unpack_real()  # ir_rtpi
    data.unpack_int()  # ir_nstcomm
    abs(data.unpack_int())  # ir_comm_mode

    data.unpack_int()  # ir_nstcheckpoint
    data.unpack_int()  # ir_nstcgsteep
    data.unpack_int()  # ir_nbfgscorr

    data.unpack_int()  # ir_nstlog
    data.unpack_int()  # ir_nstxout
    data.unpack_int()  # ir_nstvout
    data.unpack_int()  # ir_nstfout
    data.unpack_int()  # ir_nstenergy
    data.unpack_int()  # ir_nstxtcout

    data.unpack_real()  # ir_init_t = rdum =
    data.unpack_real()  # ir_delta_t = rdum =

    data.unpack_real()  # ir_xtcprec
    ir_rlist = data.unpack_real()

    data.unpack_int()  # ir_coulombtype
    data.unpack_real()  # ir_rcoulomb_switch
    ir_rcoulomb = data.unpack_real()

    data.unpack_int()  # ir_rvdwtype
    data.unpack_real()  # ir_rvdw_switch
    ir_rvdw = data.unpack_real()

    max(ir_rlist, max(ir_rvdw, ir_rcoulomb))  # ir_rlistlong

    data.unpack_int()  # ir_eDispCorr
    data.unpack_real()  # ir_epsilon_r
    data.unpack_real()  # ir_epsilon_rf

    data.unpack_real()  # ir_tabext

    data.unpack_int()  # ir_gb_algorithm
    data.unpack_int()  # ir_nstgbradii
    data.unpack_real()  # ir_rgbradii
    data.unpack_real()  # ir_gb_saltconc
    data.unpack_int()  # ir_implicit_solvent

    data.unpack_real()  # ir_gb_epsilon_solvent
    data.unpack_real()  # ir_gb_obc_alpha
    data.unpack_real()  # ir_gb_obc_beta
    data.unpack_real()  # ir_gb_obc_gamma

    # not relevant here
    # ir_gb_dielectric_offset = 0.009
    # ir_sa_algorithm = 0                                     # esaAPPROX

    data.unpack_real()  # ir_sa_surface_tension

    data.unpack_int()  # ir_nkx
    data.unpack_int()  # ir_nky
    data.unpack_int()  # ir_nkz
    data.unpack_int()  # ir_pme_order
    data.unpack_real()  # ir_ewald_rtol
    data.unpack_int()  # ir_ewald_geometry

    data.unpack_real()  # ir_epsilon_surface

    data.unpack_bool()  # ir_bOptFFT
    data.unpack_bool()  # ir_bContinuation
    data.unpack_int()  # ir_etc

    # not relevant here
    # ir_nsttcouple = ir_nstcalcenergy

    data.unpack_int()  # ir_epcpressure coupling
    data.unpack_int()  # ir_epctepctype, e.g. isotropic

    # not relevant here
    # ir_nstpcouple = ir_nstcalcenergy
    data.unpack_real()  # tau_p

    data.unpack_farray(DIM, data.unpack_real)  # ir_ref_p_XX
    data.unpack_farray(DIM, data.unpack_real)  # ir_ref_p_YY
    data.unpack_farray(DIM, data.unpack_real)  # ir_ref_p_ZZ

    data.unpack_farray(DIM, data.unpack_real)  # ir_compress_XX
    data.unpack_farray(DIM, data.unpack_real)  # ir_compress_YY
    data.unpack_farray(DIM, data.unpack_real)  # ir_compress_ZZ

    data.unpack_int()  # ir_refcoord_scaling
    data.unpack_farray(DIM, data.unpack_real)  # ir_posres_com
    data.unpack_farray(DIM, data.unpack_real)  # ir_posres_comB

    data.unpack_int()  # ir_andersen_seed
    data.unpack_real()  # ir_shake_tol
    data.unpack_int()  # ir_efep

    data.unpack_real()  # ir_init_lambda = rdum =
    data.unpack_real()  # ir_delta_lambda = rdum =

    # Not relevant here
    # ir_n_flambda = 0
    # ir_flambda   = None

    data.unpack_real()  # ir_sc_alpha
    data.unpack_int()  # ir_sc_power
    data.unpack_real()  # ir_sc_sigma

    # not relevant here
    # ir_sc_sigma_min = 0
    # ir_nstdhdl = 1
    # ir_separate_dhdl_file  = 0                            # epdhdlfileYES;
    # ir_dhdl_derivatives = 0                               # dhdlderivativesYES
    # ir_dh_hist_size    = 0
    # ir_dh_hist_spacing = 0.1

    data.unpack_int()  # ir_eDisre
    data.unpack_int()  # ir_eDisre_weighting
    data.unpack_bool()  # ir_bDisreMixed

    data.unpack_real()  # ir_dr_fc
    data.unpack_real()  # ir_dr_tau
    data.unpack_int()  # ir_nstdisreout

    data.unpack_real()  # ir_orires_fc
    data.unpack_real()  # ir_orires_tau
    data.unpack_int()  # ir_nstorireout

    data.unpack_real()  # ir_dihre_fc

    data.unpack_real()  # ir_em_stepsize
    data.unpack_real()  # ir_em_tol

    data.unpack_bool()  # ir_bShakeSOR
    data.unpack_int()  # ir_niter

    data.unpack_real()  # ir_fc_stepsize

    data.unpack_int()  # ir_eConstrAlg
    data.unpack_int()  # ir_nProjOrder
    data.unpack_real()  # ir_LincsWarnAngle
    data.unpack_int()  # ir_nLincsIter

    data.unpack_real()  # ir_bd_fric
    data.unpack_int()  # ir_ld_seed

    U.ndo_rvec(data, DIM)  # ir_deform

    data.unpack_real()  # ir_cos_accel

    data.unpack_int()  # ir_userint1
    data.unpack_int()  # ir_userint2
    data.unpack_int()  # ir_userint3
    data.unpack_int()  # ir_userint4
    data.unpack_real()  # ir_userreal1
    data.unpack_real()  # ir_userreal2
    data.unpack_real()  # ir_userreal3
    data.unpack_real()  # ir_userreal4

    # pull_stuff
    data.unpack_int()  # ir_ePull

    # grpopts stuff
    ir_opts_ngtc = data.unpack_int()
    # not relevant here
    # ir_opts_nhchainlength = 1

    ir_opts_ngacc = data.unpack_int()
    ir_opts_ngfrz = data.unpack_int()
    ir_opts_ngener = data.unpack_int()

    U.ndo_real(data, ir_opts_ngtc)  # ir_nrdf
    U.ndo_real(data, ir_opts_ngtc)  # ir_ref_t
    U.ndo_real(data, ir_opts_ngtc)  # ir_tau_t

    if ir_opts_ngfrz > 0:
        U.ndo_ivec(data, ir_opts_ngfrz)  # ir_opts_nFreeze

    if ir_opts_ngacc > 0:
        U.ndo_rvec(data, ir_opts_ngacc)  # ir_opts_acc

    U.ndo_int(data, ir_opts_ngener ** 2)  # ir_opts_egp_flags
    U.ndo_int(data, ir_opts_ngtc)  # ir_opts_annealing
    ir_opts_anneal_npoints = U.ndo_int(data, ir_opts_ngtc)

    ir_opts_anneal_time = []
    ir_opts_anneal_temp = []
    for j in range(ir_opts_ngtc):
        k = ir_opts_anneal_npoints[j]
        ir_opts_anneal_time.append(U.ndo_int(data, k))
        ir_opts_anneal_temp.append(U.ndo_int(data, k))

    # Walls
    data.unpack_int()  # ir_nwall
    data.unpack_int()  # ir_nwall_type
    data.unpack_real()  # ir_wall_r_linpot

    # ir->wall_atomtype[0], ir->wall_atomtype[1]
    ir_wall_atomtype = []
    ir_wall_atomtype.append(data.unpack_int())
    ir_wall_atomtype.append(data.unpack_int())

    # ir->wall_density[0], ir->wall_density[1]
    ir_wall_density = []
    ir_wall_density.append(data.unpack_real())
    ir_wall_density.append(data.unpack_real())

    data.unpack_real()  # ir_wall_ewald_zfac

    # cosine stuff for electric fields
    ir_ex_n, ir_et_n, ir_ex_a, ir_ex_phi, ir_et_a, ir_et_phi = [], [], [], [], [], []
    for j in range(DIM):
        x = data.unpack_int()
        ir_ex_n.append(x)
        y = data.unpack_int()
        ir_et_n.append(y)

        ir_ex_a.append(U.ndo_real(data, x))
        ir_ex_phi.append(U.ndo_real(data, x))
        ir_et_a.append(U.ndo_real(data, y))
        ir_et_phi.append(U.ndo_real(data, y))

    # QMM stuff
    data.unpack_bool()  # ir_bQMMM
    data.unpack_int()  # ir_bQMMMscheme
    data.unpack_real()  # ir_scalefactor
    data.unpack_int()  # ir_opts_ngQM

    # if ir_opts_ngQM > 0:
    # do_something

    # indicating the parsing finishes properly
    data.done()
