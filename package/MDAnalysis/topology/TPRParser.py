#!/usr/bin/env python

__author__      = "Zhuyi Xue"
__copyright__   = "GNU Public Licence, v2"

"""
Interested: 
Atoms: number, name, type, resname, resid, segid, mass, charge, 
       [residue, segment, radius, bfactor, resnum]
Bonds:
Angels:
Dihedrals:
Impropers:

Potential Bug: in the result of gmxdump, the "Proper Dih.:" section is actually
a list of Improper Dih.

This tpr parser is written according to the following files
- <gromacs_dir>/src/kernel/gmxdump.c
- <gromacs_dir>/src/gmxlib/tpxio.c (the most important one)
- <gromacs_dir>/src/gmxlib/gmxfio_rw.c
- <gromacs_dir>/src/gmxlib/gmxfio_xdr.c
- <gromacs_dir>/include/gmxfiofio.h

function read_tpxheader is based on
http://code.google.com/p/mdanalysis/wiki/TPRReaderDevelopment

functions with name starting at read_ or do_ is trying to be similar to those
in gmxdump.c or tpxio.c, those with extract_ is new

"""

import xdrlib

from MDAnalysis.core.AtomGroup import Atom
from MDAnalysis.topology.core import guess_atom_type

import tpr_obj as obj
import tpr_setting as setting

class NotImplementedError(Exception):
    """This code only works for certain features in certain version of tpr files""" 
    pass

ndo_int = setting.ndo_int
ndo_real = setting.ndo_real
ndo_rvec = setting.ndo_rvec
ndo_ivec = setting.ndo_ivec

def do_symstr(data, symtab):
    #do_symstr: get a string based on index from the symtab 
    ndx = data.unpack_int()
    return symtab[ndx]

def read_tpxheader(data):
    number = data.unpack_int()                     # ?
    version_string = data.unpack_string()          # e.g. version VERSION 4.0.5
    precision = data.unpack_int()                  # e.g. 4
    fver = data.unpack_int()                       # version of tpx file
    if fver != 58:
        raise NotImplementedError(
            "The tpx version is {0}, this parser only works for version 58 currently".format(fver))

    fgen = data.unpack_int()                       # generation of tpx file, e.g. 17
    tpx_natoms = data.unpack_int()                 # total number of atoms
    tpx_ngtc = data.unpack_int()                   # number of groups for T-coupling
    idum = data.unpack_int()                       # ?
    rdum = data.unpack_float()                     # ?
    tpx_lambda = data.unpack_float()               # ?
    tpx_bIr =  data.unpack_int()                   # has input record or not
    tpx_bTop =  data.unpack_int()                  # has topology or not
    tpx_bX =  data.unpack_int()                    # has coordinates or not 
    tpx_bV =  data.unpack_int()                    # has velocity or not
    tpx_bF =  data.unpack_int()                    # has force or not
    tpx_bBox =  data.unpack_int()                  # has box or not

    return (number, version_string, precision, fver, 
            fgen, tpx_natoms, tpx_ngtc, idum, rdum, tpx_lambda, 
            tpx_bIr, tpx_bTop, tpx_bX, tpx_bV, tpx_bF, tpx_bBox)

def extract_box(data):
    box = ndo_rvec(data, setting.DIM)
    box_rel = ndo_rvec(data, setting.DIM)
    box_v = ndo_rvec(data, setting.DIM)
    return box, box_rel, box_v

def do_symtab(data):
    symtab_nr = data.unpack_int()                           # number of symbols
    symtab = []
    for i in range(symtab_nr):
        j = data.unpack_fstring(1)              # strings are separated by void
        j = data.unpack_string()
        symtab.append(j)
    return symtab

def do_ffparams(data, fver):
    data.unpack_int()                                       # ffparams_atnr
    ffparams_ntypes = data.unpack_int()
    functypes = ndo_int(data, ffparams_ntypes)
    data.unpack_float()                                     # fudgeQQ

    # mimicing the c code,
    # remapping the functypes due to inconsistency in different versions
    for k in range(len(functypes)):
        for j in setting.ftupd:
            if fver < j[0] and functypes[k] >= j[1]:
                functypes[k] += 1
    return functypes

def do_harm(data):
    data.unpack_float()                                     # rA
    data.unpack_float()                                     # krA
    data.unpack_float()                                     # rB
    data.unpack_float()                                     # krB
    
def do_iparams(data, functypes):
    for k, i in enumerate(functypes):
        if i in [setting.F_LJ]:
            data.unpack_float()                             # lj_c6
            data.unpack_float()                             # lj_c9
        elif i in [setting.F_ANGLES, setting.F_G96ANGLES, 
                   setting.F_BONDS, setting.F_G96BONDS, 
                   setting.F_HARMONIC, setting.F_IDIHS]:
            do_harm(data)
        elif i in [setting.F_PIDIHS, setting.F_ANGRES, 
                   setting.F_ANGRESZ, setting.F_PDIHS]:
            data.unpack_float()                             # pdihs_phiA
            data.unpack_float()                             # pdihs_cpA
            data.unpack_float()                             # pdihs_phiB
            data.unpack_float()                             # pdihs_cpB
            data.unpack_int()                               # pdihs_mult
        elif i in [setting.F_RBDIHS]:
            ndo_real(data, setting.NR_RBDIHS)             # iparams_rbdihs_rbcA
            ndo_real(data, setting.NR_RBDIHS)             # iparams_rbdihs_rbcB
        elif i in [setting.F_CONSTR, setting.F_CONSTRNC]:
            data.unpack_float()                             # dA
            data.unpack_float()                             # dB
        elif i in [setting.F_SETTLE]:
            data.unpack_float()                            # iparams_settle_doh
            data.unpack_float()                            # iparams_settle_dhh
        elif i in [setting.F_VSITE3, setting.F_VSITE3FD, setting.F_VSITE3FAD]:
            data.unpack_float()                             # vsite_a
            data.unpack_float()                             # vsite_b
        elif i in [setting.F_LJ14]:
             data.unpack_float()                            # lj14_c6A
             data.unpack_float()                            # lj14_c12A
             data.unpack_float()                            # lj14_c6B
             data.unpack_float()                            # lj14_c12B
        else:
            raise NotImplementedError("functype: {0} not implemented yet".format(i))

def do_moltype(data, symtab, fver):
    molname = do_symstr(data, symtab)
    # key info: about atoms
    atomkinds = do_atoms(data, symtab)
    # key info: about bonds, angles, dih, improp dih.
    zip_ilists = do_ilists(data, fver)
    # these may not be available for certain molecules, e.g. tip4p
    bonds, angles, dihs, impr = None, None, None, None
    for zip_ilist in zip_ilists:
        nr = zip_ilist[0]
        if nr:
            itrt_t, atoms_ndx = zip_ilist[1:]   # interaction_type, atoms index
            itrt_type = obj.InteractionType(*itrt_t)
            # the following if..elif..else statement needs to be updated as new
            # type of interactions become interested
            if itrt_type.name in ['LJ14']:
                bonds = list(itrt_type.process(atoms_ndx))
            elif itrt_type.name in ['ANGLES']:
                angles = list(itrt_type.process(atoms_ndx))
            elif itrt_type.name in ['RBDIHS']:
                dihs = list(itrt_type.process(atoms_ndx))
            elif itrt_type.name in ['PDIHS', 'VSITE3']:
                # this is possible a bug in gromacs, the so-named Proper Dih is
                # actually Improper Dih
                impr = list(itrt_type.process(atoms_ndx))
            else:
                # other interaction types are not interested at the moment
                pass
    molecule_type = obj.MoleculeType(molname, atomkinds, bonds, angles, dihs, impr)

    # info in do_block and do_blocka is not interested
    do_block(data)
    do_blocka(data)

    return molecule_type

def do_atoms(data, symtab):
    # do_atoms
    atoms_nr = data.unpack_int()
    atoms_nres = data.unpack_int()

    (m, q, mB, qB, type, typeB, ptype, resind, atomnumber) = (
        [], [], [], [], [], [], [], [], [])

    for i in range(atoms_nr):
        m_, q_, mB_, qB_, type_, typeB_, ptype_, resind_, atomnumber_ = do_atom(data)
        m.append(m_)
        q.append(q_)
        mB.append(mB_)
        qB.append(qB_)
        type.append(type_)
        typeB.append(typeB_)
        ptype.append(ptype_)            # regular atom, virtual site or others
        resind.append(resind_)          # index of residue
        atomnumber.append(atomnumber)   # index of atom type

    # do_strstr
    atomnames = [symtab[i] for i in ndo_int(data, atoms_nr)]
    ndo_int(data, atoms_nr) # list of ndx, type    : symtab[ndx], e.g. opls_111
    ndo_int(data, atoms_nr) # list of ndx, typeB   : symtab[ndx]
    resnames = [symtab[i] for i in ndo_int(data, atoms_nres)]

    atom_resnames = [resnames[i] for i in resind]
    types = [guess_atom_type(i) for i in atomnames]

    # since atomtype has already been used
    atomkinds = [obj.AtomKind(i, *v) for i, v in
                 enumerate(zip(atomnames, types, resind, atom_resnames, m, q))]
    return atomkinds

def do_atom(data):
    atom_m = data.unpack_float()
    atom_q = data.unpack_float()
    atom_mB = data.unpack_float()
    atom_qB = data.unpack_float()
    atom_type = data.unpack_uint()
    atom_typeB = data.unpack_uint()
    atom_ptype = data.unpack_int()
    atom_resind = data.unpack_int()
    atomnumber = data.unpack_int()
    return (atom_m, atom_q, atom_mB, atom_qB, atom_type, 
            atom_typeB, atom_ptype, atom_resind, atomnumber)

def do_ilists(data, fver):
    ilist_nr = []                                           # number of ilist
    ilist_iatoms = []
    for j in range(setting.F_NRE):                   # total number of energies
        bClear = False
        for k in setting.ftupd:
            if fver < k[0] and j == k[1]:
                bClear = True
        if bClear:
            ilist_nr.append(0)
            ilist_iatoms.append(None)
        else:
            # do_ilist
            n = data.unpack_int()
            ilist_nr.append(n)
            l_ = []
            for i in range(n):
                l_.append(data.unpack_int())
            ilist_iatoms.append(l_)

    return zip(ilist_nr, setting.interaction_types, ilist_iatoms)
    # print "lengths of ilist_iatoms: {0}".format([len(i) for i in ilist_iatoms if i])

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
    return molb_type, molb_nmol, molb_natoms_mol, molb_nposres_xA, molb_nposres_xB

def do_atomtypes(data):
    at_nr = data.unpack_int()                               # at: atomtype
    at_radius = ndo_real(data, at_nr)
    at_vol = ndo_real(data, at_nr)
    at_surftens = ndo_real(data, at_nr)
    at_atomnumber = ndo_int(data, at_nr)
    return at_radius, at_vol, at_surftens, at_atomnumber

def do_grps(data):
    grps_nr = []
    myngrps = ngrps = setting.egcNR           # remind of version inconsistency
    for j in range(ngrps):
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
    for i in range(setting.egcNR):
        x = data.unpack_int()
        ngrpnr.append(x)
        if x == 0:
            grpnr.append(None)
        else:
            l_ = []
            for i in range(x):
                l_.append(data.unpack_uint())
            grpnr.append(l_)
    # print ngrpnr
    # print [len(i) for i in grpnr if i]
    return

def do_mtop(data, fver):
    # mtop: the topology of the whole system
    symtab = do_symtab(data)
    do_symstr(data, symtab)                                 # system_name

    functypes = do_ffparams(data, fver)
    do_iparams(data, functypes)            # parameters for different functions

    nmoltype = data.unpack_int()
    molecule_types = []                                      # non-gromacs
    for i in range(nmoltype):
        moltype = do_moltype(data, symtab, fver)
        molecule_types.append(moltype)

    nmolblock = data.unpack_int()

    atoms, bonds, angles, dihe, impr = [], [], [], [], []
    atom_start_ndx = 0
    res_start_ndx = 0

    for i in range(nmolblock):
        # molb_type is just an index for moltypes/molecule_types
        (molb_type, molb_nmol, molb_natoms_mol, 
         molb_nposres_xA, molb_nposres_xB) = do_molblock(data)
        # segment is made to correspond to the molblock as in gromacs, the
        # naming is kind of arbitrary
        segid = "seg_{0}_{1}".format(i, molecule_types[molb_type].name)
        for j in xrange(molb_nmol):
            mt = molecule_types[molb_type]                  # mt: molecule type
            for atomkind in mt.atomkinds:
                atoms.append(Atom(atomkind.id + atom_start_ndx, 
                                  atomkind.name,
                                  atomkind.type,
                                  atomkind.resname,
                                  atomkind.resid + res_start_ndx,
                                  segid,
                                  atomkind.mass,
                                  atomkind.charge))
            # remap_ method returns [blah, blah, ..] or []
            bonds.extend(mt.remap_bonds(atom_start_ndx))
            angles.extend(mt.remap_angles(atom_start_ndx))
            dihe.extend(mt.remap_dihe(atom_start_ndx))
            impr.extend(mt.remap_impr(atom_start_ndx))

            atom_start_ndx += mt.number_of_atoms()
            res_start_ndx += mt.number_of_residues()

    data.unpack_int()                         # mtop_natoms
    do_atomtypes(data)

    # not useful here
    # mtop_ffparams_cmap_grid_ngrid        = 0
    # mtop_ffparams_cmap_grid_grid_spacing = 0.1
    # mtop_ffparams_cmap_grid_cmapdata     = 'NULL'

    do_groups(data, symtab)
    return atoms, bonds, angles, dihe, impr

def parse(tprfile):
    tprf = open(tprfile).read()
    data = xdrlib.Unpacker(tprf)
    (number, version_string, precision, fver, 
     fgen, tpx_natoms, tpx_ngtc, idum, rdum, tpx_lambda, 
     tpx_bIr, tpx_bTop, tpx_bX, tpx_bV, tpx_bF, tpx_bBox) = read_tpxheader(data)

    state_ngtc = tpx_ngtc         # done init_state() in src/gmxlib/tpxio.c
    if tpx_bBox:
        box, box_rel, box_v = extract_box(data)

    if state_ngtc > 0:
        # gmx_fio_ndo_real()
        if fver < 69:                   # redundancy due to  different versions
            ndo_real(data, state_ngtc)
        ndo_real(data, state_ngtc)        # relevant to Berendsen tcoupl_lambda

    if tpx_bTop:
        # In principle, the parsing may stop here since the following info are
        # not needed at the moment, but then it is difficult to check if tpr
        # has been parsed correctly
        atoms, bonds, angles, dihe, impr = do_mtop(data, fver)
        structure = {
            '_atoms': atoms,
            '_bonds': bonds,
            '_angles': angles,
            '_dihe': dihe,
            '_impr': impr
            }
    else:
        raise ValueError("No topology found in tpr file: {0}".formation(tprfile))

    if tpx_bX:
        ndo_rvec(data, tpx_natoms)

    if tpx_bV:
        ndo_rvec(data, tpx_natoms)

    if tpx_bF:
        ndo_rvec(data, tpx_natoms)

    ePBC = -1;
    bPeriodicMols = False
    if tpx_bIr:
        # update
        ePBC = data.unpack_int()
        bPeriodicMols = data.unpack_bool()
        # 17 < 23. and ir (ir is from the c code, seems not apply here
        if fgen < setting.tpx_generation:
            # a crazily long (670 lines) function in c, slightly better here
            # (240 lines), so put it in setting.py
            setting.do_inputrec(data)

    return structure

if __name__ == "__main__":
    parse("/Users/zyxue/Dropbox/help_printing/final.tpr")
    # parse("../../../../md.tpr")
