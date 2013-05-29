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

# TPR parser and tpr support module
# Copyright (c) 2011 Zhuyi Xue
# Released under the  GNU Public Licence, v2

"""
Utilities for the TPRParser
===========================

Function calling order::

   (TPRParser.py call do_mtop)
   do_mtop -> do_symtab
           -> do_ffparams -> do_iparams
           -> do_moltype  -> do_atoms  -> do_atom
                                       -> do_resinfo
                          -> do_ilists
                          -> do_block
                          -> do_blocka
           -> do_molblock

Then compose the stuffs in the format :class:`MDAnalysis.Universe` reads in.

"""

from MDAnalysis.core.AtomGroup import Atom

import obj
import setting as S

def ndo_int(data, n):
    """mimic of gmx_fio_ndo_real in gromacs"""
    return [data.unpack_int() for i in xrange(n)]

def ndo_real(data, n):
    """mimic of gmx_fio_ndo_real in gromacs"""
    return [data.unpack_float() for i in xrange(n)]

def do_rvec(data):
    return data.unpack_farray(S.DIM, data.unpack_float)

def ndo_rvec(data, n):
    """mimic of gmx_fio_ndo_rvec in gromacs"""
    return [data.unpack_farray(S.DIM, data.unpack_float) for i in xrange(n)]

def ndo_ivec(data, n):
    """mimic of gmx_fio_ndo_rvec in gromacs"""
    return [data.unpack_farray(S.DIM, data.unpack_int) for i in xrange(n)]

def fver_err(fver):
    if fver not in S.SUPPORTED_VERSIONS:
        raise NotImplementedError(
            "Your tpx version is {0}, which this parser does not support, yet ".format(
                fver))

def read_tpxheader(data):
    """this function is now compatible with do_tpxheader in tpxio.c"""
    number = data.unpack_int()                     # ?
    ver_str = data.unpack_string()                 # version string e.g. VERSION 4.0.5
    precision = data.unpack_int()                  # e.g. 4
    fver = data.unpack_int()                       # version of tpx file
    fver_err(fver)

    fgen = data.unpack_int() if fver >= 26 else 0 # generation of tpx file, e.g. 17

    # Versions before 77 don't have the tag, set it to TPX_TAG_RELEASE file_tag
    # file_tag is used for comparing with tpx_tag. Only tpr files with a
    # tpx_tag from a lower or the same version of gromacs code can be parsed by
    # the tpxio.c

    if fver >= 80:
        data.unpack_int()       # the value is 8, but haven't found the
                                # corresponding code in the
                                # <gromacs-4.6.1-dir>/src/gmxlib/tpxio.c yet.
        file_tag = data.unpack_string()
    else:
        file_tag = S.TPX_TAG_RELEASE

    natoms = data.unpack_int()                    # total number of atoms
    ngtc = data.unpack_int() if fver >= 28 else 0 # number of groups for T-coupling

    if fver < 62:
        # not sure what these two are for.
        data.unpack_int()                             # idum
        data.unpack_float()                           # rdum

    fep_state = data.unpack_int() if fver > 79 else 0

    # actually, it's lambda, not sure what is it. us lamb because lambda is a
    # keywod in python
    lamb = data.unpack_float()
    bIr =  data.unpack_int()                          # has input record or not
    bTop =  data.unpack_int()                         # has topology or not
    bX =  data.unpack_int()                           # has coordinates or not
    bV =  data.unpack_int()                           # has velocity or not
    bF =  data.unpack_int()                           # has force or not
    bBox =  data.unpack_int()                         # has box or not

    th = obj.TpxHeader(number, ver_str, precision,
                   fver, fgen, file_tag, natoms, ngtc, fep_state, lamb,
                   bIr, bTop, bX, bV, bF, bBox)
    return th

def extract_box_info(data, fver):
    box = ndo_rvec(data, S.DIM)
    box_rel = ndo_rvec(data, S.DIM) if fver >= 51 else 0
    box_v = ndo_rvec(data, S.DIM) if fver >= 28 else None
    if (fver < 56):
        ndo_rvec(data, S.DIM)                         # mdum?

    return obj.Box(box, box_rel, box_v)

def do_mtop(data, fver):
    # mtop: the topology of the whole system
    symtab = do_symtab(data)
    do_symstr(data, symtab)                                 # system_name
    do_ffparams(data, fver)                                 # params

    nmoltype = data.unpack_int()
    moltypes = []                                     # non-gromacs
    for i in xrange(nmoltype):
        moltype = do_moltype(data, symtab, fver)
        moltypes.append(moltype)

    nmolblock = data.unpack_int()

    mtop = obj.Mtop(nmoltype, moltypes, nmolblock)


    ttop = obj.TPRTopology(*[[] for i in xrange(5)])

    atom_start_ndx = 0
    res_start_ndx = 0
    for i in xrange(mtop.nmolblock):
        # molb_type is just an index for moltypes/molecule_types
        mb = do_molblock(data)
        # segment is made to correspond to the molblock as in gromacs, the
        # naming is kind of arbitrary
        segid = "seg_{0}_{1}".format(i, mtop.moltypes[mb.molb_type].name)
        for j in xrange(mb.molb_nmol):
            mt = mtop.moltypes[mb.molb_type]                  # mt: molecule type
            for atomkind in mt.atomkinds:
                ttop.atoms.append(Atom(atomkind.id + atom_start_ndx,
                                       atomkind.name,
                                       atomkind.type,
                                       atomkind.resname,
                                       atomkind.resid + res_start_ndx,
                                       segid,
                                       atomkind.mass,
                                       atomkind.charge))
            # remap_ method returns [blah, blah, ..] or []
            ttop.bonds.extend(mt.remap_bonds(atom_start_ndx))
            ttop.angles.extend(mt.remap_angles(atom_start_ndx))
            ttop.dihe.extend(mt.remap_dihe(atom_start_ndx))
            ttop.impr.extend(mt.remap_impr(atom_start_ndx))

            atom_start_ndx += mt.number_of_atoms()
            res_start_ndx += mt.number_of_residues()

    # not useful here

    # data.unpack_int()                         # mtop_natoms
    # do_atomtypes(data)
    # mtop_ffparams_cmap_grid_ngrid        = 0
    # mtop_ffparams_cmap_grid_grid_spacing = 0.1
    # mtop_ffparams_cmap_grid_cmapdata     = 'NULL'
    # do_groups(data, symtab)
    return ttop

def do_symstr(data, symtab):
    #do_symstr: get a string based on index from the symtab
    ndx = data.unpack_int()
    return symtab[ndx]

def do_symtab(data):
    symtab_nr = data.unpack_int()                           # number of symbols
    symtab = []
    for i in xrange(symtab_nr):
        j = data.unpack_fstring(1)              # strings are separated by void
        j = data.unpack_string()
        symtab.append(j)
    return symtab

def do_ffparams(data, fver):
    atnr = data.unpack_int()
    if fver < 57:
        data.unpack_int()                                   # idum
    ntypes = data.unpack_int()
    functype = ndo_int(data, ntypes)
    reppow = data.unpack_double() if fver >= 66 else 12.0
    if fver >= 57:
        fudgeQQ = data.unpack_float()

    # mimicing the c code,
    # remapping the functype due to inconsistency in different versions
    for i in xrange(len(functype)):
        for k in S.ftupd:
            # j[0]: tpx_version, j[1] funtype
            if fver < k[0] and functype[i] >= k[1]:
                functype[i] += 1

    # parameters for different functions, None returned for now since not sure
    # what is iparams
    iparams = do_iparams(data, functype, fver)

    params = obj.Params(atnr, ntypes, functype, reppow, fudgeQQ, iparams)
    return params

def do_harm(data):
    data.unpack_float()                                     # rA
    data.unpack_float()                                     # krA
    data.unpack_float()                                     # rB
    data.unpack_float()                                     # krB

def do_iparams(data, functypes, fver):
    # Not all elif cases in this function has been used and tested
    for k, i in enumerate(functypes):
        if i in [S.F_ANGLES, S.F_G96ANGLES,
                   S.F_BONDS, S.F_G96BONDS,
                   S.F_HARMONIC, S.F_IDIHS]:
            do_harm(data)
        elif i in [S.F_FENEBONDS]:
            data.unpack_float()                             # fene.bm
            data.unpack_float()                             # fene.kb
        elif i in [S.F_RESTRBONDS]:
            data.unpack_float()                             # restraint.lowA
            data.unpack_float()                             # restraint.up1A
            data.unpack_float()                             # restraint.up2A
            data.unpack_float()                             # restraint.kA
            data.unpack_float()                             # restraint.lowB
            data.unpack_float()                             # restraint.up1B
            data.unpack_float()                             # restraint.up2B
            data.unpack_float()                             # restraint.kB
        elif i in [S.F_TABBONDS, S.F_TABBONDSNC, S.F_TABANGLES, S.F_TABDIHS]:
            data.unpack_float()                             # tab.kA
            data.unpack_int()                               # tab.table
            data.unpack_float()                             # tab.kB
        elif i in [S.F_CROSS_BOND_BONDS]:
            data.unpack_float()                             # cross_bb.r1e
            data.unpack_float()                             # cross_bb.r2e
            data.unpack_float()                             # cross_bb.krr
        elif i in [S.F_CROSS_BOND_ANGLES]:
            data.unpack_float()                             # cross_ba.r1e
            data.unpack_float()                             # cross_ba.r2e
            data.unpack_float()                             # cross_ba.r3e
            data.unpack_float()                             # cross_ba.krt
        elif i in [S.F_UREY_BRADLEY]:
            data.unpack_float()                             # u_b.theta
            data.unpack_float()                             # u_b.ktheta
            data.unpack_float()                             # u_b.r13
            data.unpack_float()                             # u_b.kUB
        elif i in [S.F_QUARTIC_ANGLES]:
            data.unpack_float()                             # qangle.theta
            ndo_real(data, 5);                              # qangle.c

        elif i in [S.F_BHAM]:
            data.unpack_float()                             # bham.a
            data.unpack_float()                             # bham.b
            data.unpack_float()                             # bham.c
        elif i in [S.F_MORSE]:
            data.unpack_float()                             # morse.b0
            data.unpack_float()                             # morse.cb
            data.unpack_float()                             # morse.beta
        elif i in [S.F_CUBICBONDS]:
            data.unpack_float()                             # cubic.b0
            data.unpack_float()                             # cubic.kb
            data.unpack_float()                             # cubic.kcub
        elif i in [S.F_CONNBONDS]:
            pass
        elif i in [S.F_POLARIZATION]:
            data.unpack_float()                             # polarize.alpha
        elif i in [S.F_WATER_POL]:
            if fver < 31:
                fver_err(fver)
            data.unpack_float()                             # wpol.al_x
            data.unpack_float()                             # wpol.al_y
            data.unpack_float()                             # wpol.al_z
            data.unpack_float()                             # wpol.rOH
            data.unpack_float()                             # wpol.rHH
            data.unpack_float()                             # wpol.rOD
        elif i in [S.F_THOLE_POL]:
            data.unpack_float()                             # thole.a
            data.unpack_float()                             # thole.alpha1
            data.unpack_float()                             # thole.alpha2
            data.unpack_float()                             # thole.rfac

        elif i in [S.F_LJ]:
            data.unpack_float()                             # lj_c6
            data.unpack_float()                             # lj_c9
        elif i in [S.F_LJ14]:
            data.unpack_float()                            # lj14_c6A
            data.unpack_float()                            # lj14_c12A
            data.unpack_float()                            # lj14_c6B
            data.unpack_float()                            # lj14_c12B
        elif i in [S.F_LJC14_Q]:
            data.unpack_float()                             # ljc14.fqq
            data.unpack_float()                             # ljc14.qi
            data.unpack_float()                             # ljc14.qj
            data.unpack_float()                             # ljc14.c6
            data.unpack_float()                             # ljc14.c12
        elif i in [S.F_LJC_PAIRS_NB]:
            data.unpack_float()                             # ljcnb.qi
            data.unpack_float()                             # ljcnb.qj
            data.unpack_float()                             # ljcnb.c6
            data.unpack_float()                             # ljcnb.c12

        elif i in [S.F_PIDIHS, S.F_ANGRES,
                   S.F_ANGRESZ, S.F_PDIHS]:
            data.unpack_float()                             # pdihs_phiA
            data.unpack_float()                             # pdihs_cpA
            if (i == S.F_ANGRES or i == S.F_ANGRESZ) and fver < 42:
                data.unpack_float()                         # harmonic.rB
                data.unpack_float()                         # harmonic.krB
            else:
                data.unpack_float()                         # pdihs_phiB
                data.unpack_float()                         # pdihs_cpB
                data.unpack_int()                           # pdihs_mult

        elif i in [S.F_DISRES]:
            data.unpack_int()                               # disres.label
            data.unpack_int()                               # disres.type
            data.unpack_float()                             # disres.low
            data.unpack_float()                             # disres.up1
            data.unpack_float()                             # disres.up2
            data.unpack_float()                             # disres.kfac

        elif i in [S.F_ORIRES]:
            data.unpack_int()                               # orires.ex
            data.unpack_int()                               # orires.label
            data.unpack_int()                               # orires.power
            data.unpack_float()                             # orires.c
            data.unpack_float()                             # orires.obs
            data.unpack_float()                             # orires.kfac

        elif i in [S.F_DIHRES]:
            data.unpack_int()                               # dihres.power
            data.unpack_int()                               # dihres.label
            data.unpack_float()                             # dihres.phi
            data.unpack_float()                             # dihres.dphi
            data.unpack_float()                             # dihres.kfac

        elif i in [S.F_POSRES]:
            do_rvec(data)                                   # posres.pos0A
            do_rvec(data)                                   # posres.fcA
            if fver < 27:
                fver_err(fver)
            else:
                do_rvec(data)                               # posres.pos0B
                do_rvec(data)                               # posres.fcB

        elif i in [S.F_RBDIHS]:
            ndo_real(data, S.NR_RBDIHS)                   # iparams_rbdihs_rbcA
            if fver >= 25:
                ndo_real(data, S.NR_RBDIHS)               # iparams_rbdihs_rbcB

        elif i in [S.F_FOURDIHS]:
            # Fourier dihedrals
            ndo_real(data, S.NR_RBDIHS);                    # rbdihs.rbcA
            ndo_real(data, S.NR_RBDIHS);                    # rbdihs.rbcB

        elif i in [S.F_CONSTR, S.F_CONSTRNC]:
            data.unpack_float()                             # dA
            data.unpack_float()                             # dB

        elif i in [S.F_SETTLE]:
            data.unpack_float()                             # settle.doh
            data.unpack_float()                             # settle.dhh

        elif i in [S.F_VSITE2]:
            data.unpack_float()                             # vsite.a

        elif i in [S.F_VSITE3, S.F_VSITE3FD, S.F_VSITE3FAD]:
            data.unpack_float()                             # vsite.a
            data.unpack_float()                             # vsite.b

        elif i in [S.F_VSITE3OUT, S.F_VSITE4FD, S.F_VSITE4FDN]:
            data.unpack_float()                             # vsite.a
            data.unpack_float()                             # vsite.b
            data.unpack_float()                             # vsite.c

        elif i in [S.F_VSITEN]:
            data.unpack_int()                               # vsiten.n
            data.unpack_float()                             # vsiten.a

        elif i in [S.F_GB12, S.F_GB13, S.F_GB14]:
            # /* We got rid of some parameters in version 68 */
            if fver < 68:
                data.unpack_float()                         # rdum
                data.unpack_float()                         # rdum
                data.unpack_float()                         # rdum
                data.unpack_float()                         # rdum
            data.unpack_float()                             # gb.sar
            data.unpack_float()                             # gb.st
            data.unpack_float()                             # gb.pi
            data.unpack_float()                             # gb.gbr
            data.unpack_float()                             # gb.bmlt

        elif i in [S.F_CMAP]:
            data.unpack_int()                               # cmap.cmapA
            data.unpack_int()                               # cmap.cmapB
        else:
            raise NotImplementedError("unknown functype: {0}".format(i))
    return

def do_moltype(data, symtab, fver):
    if fver >= 57:
        molname = do_symstr(data, symtab)

    # key info: about atoms
    atoms_obj = do_atoms(data, symtab, fver)

    #### start: MDAnalysis specific
    atomkinds = []
    for k, a in enumerate(atoms_obj.atoms):
        atomkinds.append(obj.AtomKind(
                k,
                atoms_obj.atomnames[k],
                atoms_obj.type[k],
                a.resind,
                atoms_obj.resnames[a.resind],
                a.m,
                a.q))
    #### end: MDAnalysis specific

    # key info: about bonds, angles, dih, improp dih.
    ilists = do_ilists(data, fver)

    #### start: MDAnalysis specific
    # these may not be available for certain molecules, e.g. tip4p
    bonds, angles, dihs, impr = None, None, None, None
    for ilist in ilists:
        if ilist.nr > 0:
            ik_obj = obj.InteractionKind(*ilist.ik)
            ias = ilist.iatoms

            # the following if..elif..else statement needs to be updated as new
            # type of interactions become interested
            if ik_obj.name in ['LJ14']:
                bonds = list(ik_obj.process(ias))
            elif ik_obj.name in ['ANGLES']:
                angles = list(ik_obj.process(ias))
            elif ik_obj.name in ['RBDIHS']:
                dihs = list(ik_obj.process(ias))
            elif ik_obj.name in ['PDIHS', 'VSITE3']:
                # this is possible a bug in gromacs, the so-named Proper Dih is
                # actually Improper Dih
                impr = list(ik_obj.process(ias))
            else:
                # other interaction types are not interested at the moment
                pass

    moltype = obj.MoleculeKind(molname, atomkinds, bonds, angles, dihs, impr)
    #### end: MDAnalysis specific

    # info in do_block and do_blocka is not interested, but has to be parsed
    # here so that all moltypes can be parsed properly
    do_block(data)
    do_blocka(data)
    return moltype

def do_atoms(data, symtab, fver):
    nr = data.unpack_int()     # number of atoms in a particular molecule
    nres = data.unpack_int() # number of residues in a particular molecule

    if fver < 57:
        fver_err(fver)

    atoms = []
    for i in xrange(nr):
        A = do_atom(data, fver)
        atoms.append(A)

    # do_strstr
    atomnames = [symtab[i] for i in ndo_int(data, nr)]

    if fver <= 20:
        fver_err(fver)
    else:
        type  = [symtab[i] for i in ndo_int(data, nr)]      # e.g. opls_111
        typeB = [symtab[i] for i in ndo_int(data, nr)]
    resnames = do_resinfo(data, symtab, fver, nres)

    if fver < 57:
        fver_err(fver)

    return obj.Atoms(atoms, nr, nres, type, typeB, atomnames, resnames)

def do_resinfo(data, symtab, fver, nres):
    if fver < 63:
        resnames = [symtab[i] for i in ndo_int(data, nres)]
    else:
        resnames = []
        for i in xrange(nres):
            resnames.append(symtab[data.unpack_int()])
            # assume the uchar in gmx is 8 byte, seems right
            data.unpack_fstring(8)
    return resnames

def do_atom(data, fver):
    m = data.unpack_float()                                 # mass
    q = data.unpack_float()                                 # charge
    mB = data.unpack_float()
    qB = data.unpack_float()
    tp = data.unpack_uint()                       # type is a keyword in python
    typeB = data.unpack_uint()
    ptype = data.unpack_int()       # regular atom, virtual site or others
    resind = data.unpack_int()      # index of residue

    if fver >= 52:
        atomnumber = data.unpack_int()       # index of atom type

    if fver < 23 or fver < 39  or  fver < 57:
        fver_err(fver)
    return obj.Atom(m, q, mB, qB, tp, typeB, ptype, resind, atomnumber)

def do_ilists(data, fver):
    nr = []                   # number of ilist
    iatoms = []               # atoms involved in a particular interaction type
    for j in xrange(S.F_NRE):  # total number of energies (i.e. interaction types)
        bClear = False
        for k in S.ftupd:
            if fver < k[0] and j == k[1]:
                bClear = True
        if bClear:
            nr.append(0)
            iatoms.append(None)
        else:
            if fver < 44:
                fver_err(fver)
            # do_ilist
            n = data.unpack_int()
            nr.append(n)
            l_ = []
            for i in xrange(n):
                l_.append(data.unpack_int())
            iatoms.append(l_)

    return [obj.Ilist(n, it, i) for n, it, i in
            zip(nr, S.interaction_types, iatoms)]

def do_molblock(data):
    molb_type = data.unpack_int()
    molb_nmol = data.unpack_int()            # number of moles in the molblock
    molb_natoms_mol = data.unpack_int()      # number of atoms in a molecule
    molb_nposres_xA = data.unpack_int()
    if molb_nposres_xA > 0:
        ndo_rvec(data, molb_nposres_xA)
    molb_nposres_xB = data.unpack_int() # The number of posres coords for top B
    if molb_nposres_xB > 0:
        ndo_rvec(data, molb_nposres_xB)

    return obj.Molblock(molb_type, molb_nmol, molb_natoms_mol,
                    molb_nposres_xA, molb_nposres_xB)

def do_block(data):
    block_nr = data.unpack_int()                       # for cgs: charge groups
    # starting or ending atom indices, based on which cgs are grouped
    ndo_int(data, block_nr + 1)
    return do_block

def do_blocka(data):
    block_nr = data.unpack_int()  # No. of atoms with excls
    block_nra = data.unpack_int() # total times fo appearance of atoms for excls
    ndo_int(data, block_nr + 1)
    ndo_int(data, block_nra)
    return block_nr, block_nra

##############UTILS FOR INFORMATION NOT INTERESTED AT THE MOMENT###############

def do_grps(data):
    grps_nr = []
    myngrps = ngrps = S.egcNR           # remind of version inconsistency
    for j in xrange(ngrps):
        if j < myngrps:
            v = data.unpack_int()
            grps_nr.append(v)
        ndo_int(data, v)
    return grps_nr

def do_groups(data, symtab):
    do_grps(data)

    ngrpname = data.unpack_int()
    # do_strstr, list of indices of group name: e.g. System, Protein,
    # Protein-H, etc. use symtab[i] to check
    ndo_int(data, ngrpname)

    ngrpnr = []
    grpnr = []
    for i in xrange(S.egcNR):
        x = data.unpack_int()
        ngrpnr.append(x)
        if x == 0:
            grpnr.append(None)
        else:
            l_ = []
            for i in xrange(x):
                l_.append(data.unpack_uint())
            grpnr.append(l_)
    # print ngrpnr
    # print [len(i) for i in grpnr if i]
    return

def do_atomtypes(data):
    at_nr = data.unpack_int()                               # at: atomtype
    at_radius = ndo_real(data, at_nr)
    at_vol = ndo_real(data, at_nr)
    at_surftens = ndo_real(data, at_nr)
    at_atomnumber = ndo_int(data, at_nr)
    return at_radius, at_vol, at_surftens, at_atomnumber
