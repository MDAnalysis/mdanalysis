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

The module also contains the :func:`do_inputrec` to read the TPR header with.
"""
from __future__ import absolute_import

from six.moves import range
import numpy as np
import xdrlib
import struct

from . import obj
from . import setting as S
from ..base import squash_by
from ...core.topology import Topology
from ...core.topologyattrs import (
    Atomids,
    Atomnames,
    Atomtypes,
    Masses,
    Charges,
    Resids,
    Resnames,
    Moltypes,
    Molnums,
    Segids,
    Bonds,
    Angles,
    Dihedrals,
    Impropers
)


class TPXUnpacker(xdrlib.Unpacker):
    """
    Extend the standard XDR unpacker for the specificity of TPX files.
    """
    def __init__(self, data):
        # Using super here would be ideal, but it does not work on python2 as
        # TPXUnpacker is a `classobj` and not a `type`.
        # super().__init__(data)
        xdrlib.Unpacker.__init__(self, data)
        self._buf = self.get_buffer()

    # The parent class uses a dunder attribute to store the
    # cursor position. This property makes it easier to manipulate
    # this attribute that is otherwise "protected".
    # The use of this property works well in python3, but fails in python2
    # where direct access to the mangled attribute seems required.
    @property
    def _pos(self):
        return self.get_position()

    @_pos.setter
    def _pos(self, value):
        self.set_position(value)

    def _unpack_value(self, item_size, struct_template):
        # Ideally, we should use the _pos attribute, but it present some
        # unexpected behaviour on python2 where the position of the cursor is
        # not kept in sync between the method defined in TPXUnpacker and the
        # methods defined in the base class.
        # start_position = self._pos
        # end_position = self._pos = start_position + item_size
        start_position = self._Unpacker__pos  # pylint: disable=access-member-before-definition 
        end_position = self._pos = self._Unpacker__pos = start_position + item_size
        content = self._buf[start_position:end_position]
        if len(content) != item_size:
            raise EOFError
        return struct.unpack(struct_template, content)[0]

    def unpack_int64(self):
        return self._unpack_value(8, '>q')

    def unpack_uint64(self):
        return self._unpack_value(8, '>Q')

    def unpack_ushort(self):
        return self.unpack_uint()

    def unpack_uchar(self):
        # TPX files prior to gromacs 2020 (tpx version < 119) use unsigned ints
        # (4 bytes) instead of unsigned chars.
        return self._unpack_value(4, '>I')

    def do_string(self):
        """
        Emulate gmx_fio_do_string

        gmx_fio_do_string reads a string from a XDR file. On the contrary to the
        python unpack_string, gmx_fio_do_string reads the size as an unsigned
        integer before reading the actual string.

        See <gromacs-2016-src>/src/gromacs/fileio/gmx_system_xdr.c:454
        """
        self.unpack_int()
        return self.unpack_string()


class TPXUnpacker2020(TPXUnpacker):
    """
    Unpacker for TPX files from and later than gromacs 2020.

    A new implementation of the serializer (InMemorySerializer), introduced in
    gromacs 2020, changes le meaning of some types in the file body (the header
    keep using the previous implementation of the serializer).
    """
    @classmethod
    def from_unpacker(cls, unpacker):
        new_unpacker = cls(unpacker._buf)
        new_unpacker._pos = new_unpacker._Unpacker__pos = unpacker._Unpacker__pos
        if hasattr(unpacker, 'unpack_real'):
            if unpacker.unpack_real == unpacker.unpack_float:
                new_unpacker.unpack_real = new_unpacker.unpack_float
            elif unpacker.unpack_real == unpacker.unpack_double:
                new_unpacker.unpack_real = new_unpacker.unpack_double
            else:
                raise ValueError("Unrecognized precision")
        return new_unpacker

    def unpack_fstring(self, n):
        if n < 0:
            raise ValueError('Size of fstring cannot be negative.')
        start_position = self._Unpacker__pos # pylint: disable=access-member-before-definition
        end_position = self._pos = self._Unpacker__pos = start_position + n
        if end_position > len(self._buf):
            raise EOFError
        content = self._buf[start_position:end_position]
        return content

    def unpack_ushort(self):
        # The InMemorySerializer implements ushort according to the XDR standard
        # on the contrary to the IO serializer.
        return self._unpack_value(2, '>H')

    def unpack_uchar(self):
        # The InMemorySerializer implements uchar according to the XDR standard
        # on the contrary to the IO serializer.
        return self._unpack_value(1, '>B')

    def do_string(self):
        """
        Emulate gmx_fio_do_string
        """
        n = self.unpack_uint64()
        return self.unpack_fstring(n)


def ndo_int(data, n):
    """mimic of gmx_fio_ndo_real in gromacs"""
    return [data.unpack_int() for i in range(n)]


def ndo_real(data, n):
    """mimic of gmx_fio_ndo_real in gromacs"""
    return [data.unpack_real() for i in range(n)]


def do_rvec(data):
    return data.unpack_farray(S.DIM, data.unpack_real)


def ndo_rvec(data, n):
    """mimic of gmx_fio_ndo_rvec in gromacs"""
    return [data.unpack_farray(S.DIM, data.unpack_real) for i in range(n)]


def ndo_ivec(data, n):
    """mimic of gmx_fio_ndo_rvec in gromacs"""
    return [data.unpack_farray(S.DIM, data.unpack_int) for i in range(n)]


def fileVersion_err(fver):
    if fver not in S.SUPPORTED_VERSIONS:
        raise NotImplementedError(
            "Your tpx version is {0}, which this parser does not support, yet ".format(
                fver))


def define_unpack_real(prec, data):
    """Define an unpack_real method of data based on the float precision used"""
    if prec == 4:
        data.unpack_real = data.unpack_float
    elif prec == 8:
        data.unpack_real = data.unpack_double
    else:
        raise ValueError("unsupported precision: {0}".format(prec))


def read_tpxheader(data):
    """this function is now compatible with do_tpxheader in tpxio.cpp
    """
    # Last compatibility check with gromacs-2016
    ver_str = data.do_string()  # version string e.g. VERSION 4.0.5
    precision = data.unpack_int()  # e.g. 4
    define_unpack_real(precision, data)
    fileVersion = data.unpack_int()  # version of tpx file
    fileVersion_err(fileVersion)

    # This is for backward compatibility with development version 77-79 where
    # the tag was, mistakenly, placed before the generation.
    if 77 <= fileVersion <= 79:
        data.unpack_int()  # the value is 8, but haven't found the
        file_tag = data.do_string()

    fileGeneration = data.unpack_int() if fileVersion >= 26 else 0  # generation of tpx file, e.g. 17

    # Versions before 77 don't have the tag, set it to TPX_TAG_RELEASE file_tag
    # file_tag is used for comparing with tpx_tag. Only tpr files with a
    # tpx_tag from a lower or the same version of gromacs code can be parsed by
    # the tpxio.c

    file_tag = data.do_string() if fileVersion >= 81 else S.TPX_TAG_RELEASE

    natoms = data.unpack_int()  # total number of atoms
    ngtc = data.unpack_int() if fileVersion >= 28 else 0  # number of groups for T-coupling

    if fileVersion < 62:
        # not sure what these two are for.
        data.unpack_int()  # idum
        data.unpack_real()  # rdum

    fep_state = data.unpack_int() if fileVersion >= 79 else 0

    # actually, it's lambda, not sure what is it. us lamb because lambda is a
    # keyword in python
    lamb = data.unpack_real()
    bIr = data.unpack_int()  # has input record or not
    bTop = data.unpack_int()  # has topology or not
    bX = data.unpack_int()  # has coordinates or not
    bV = data.unpack_int()  # has velocity or not
    bF = data.unpack_int()  # has force or not
    bBox = data.unpack_int()  # has box or not

    sizeOfTprBody = None
    if fileVersion >= S.tpxv_AddSizeField and fileGeneration >= 27:
        sizeOfTprBody = data.unpack_int64()

    th = obj.TpxHeader(ver_str, precision, fileVersion, fileGeneration,
                       file_tag, natoms, ngtc, fep_state, lamb,
                       bIr, bTop, bX, bV, bF, bBox, sizeOfTprBody)
    return th


def extract_box_info(data, fileVersion):
    box = ndo_rvec(data, S.DIM)
    box_rel = ndo_rvec(data, S.DIM) if fileVersion >= 51 else 0
    box_v = ndo_rvec(data, S.DIM) if fileVersion >= 28 else None
    if (fileVersion < 56):
        ndo_rvec(data, S.DIM)  # mdum?

    return obj.Box(box, box_rel, box_v)


def do_mtop(data, fver):
    # mtop: the topology of the whole system
    symtab = do_symtab(data)
    do_symstr(data, symtab)  # system_name
    do_ffparams(data, fver)  # params

    nmoltype = data.unpack_int()
    moltypes = []  # non-gromacs
    for i in range(nmoltype):
        moltype = do_moltype(data, symtab, fver)
        moltypes.append(moltype)

    nmolblock = data.unpack_int()

    mtop = obj.Mtop(nmoltype, moltypes, nmolblock)

    bonds = []
    angles = []
    dihedrals = []
    impropers = []

    atomids = []
    segids = []
    resids = []
    resnames = []
    atomnames = []
    atomtypes = []
    moltypes = []
    molnums = []
    charges = []
    masses = []

    atom_start_ndx = 0
    res_start_ndx = 0
    molnum = 0
    for i in range(mtop.nmolblock):
        # molb_type is just an index for moltypes/molecule_types
        mb = do_molblock(data)
        # segment is made to correspond to the molblock as in gromacs, the
        # naming is kind of arbitrary
        molblock = mtop.moltypes[mb.molb_type].name.decode('utf-8')
        segid = "seg_{0}_{1}".format(i, molblock)
        for j in range(mb.molb_nmol):
            mt = mtop.moltypes[mb.molb_type]  # mt: molecule type
            for atomkind in mt.atomkinds:
                atomids.append(atomkind.id + atom_start_ndx)
                segids.append(segid)
                resids.append(atomkind.resid + res_start_ndx)
                resnames.append(atomkind.resname.decode())
                atomnames.append(atomkind.name.decode())
                atomtypes.append(atomkind.type.decode())
                moltypes.append(molblock)
                molnums.append(molnum)
                charges.append(atomkind.charge)
                masses.append(atomkind.mass)
            molnum += 1

            # remap_ method returns [blah, blah, ..] or []
            bonds.extend(mt.remap_bonds(atom_start_ndx))
            angles.extend(mt.remap_angles(atom_start_ndx))
            dihedrals.extend(mt.remap_dihe(atom_start_ndx))
            impropers.extend(mt.remap_impr(atom_start_ndx))

            atom_start_ndx += mt.number_of_atoms()
            res_start_ndx += mt.number_of_residues()

    # not useful here

    # data.unpack_int()                         # mtop_natoms
    # do_atomtypes(data)
    # mtop_ffparams_cmap_grid_ngrid        = 0
    # mtop_ffparams_cmap_grid_grid_spacing = 0.1
    # mtop_ffparams_cmap_grid_cmapdata     = 'NULL'
    # do_groups(data, symtab)

    atomids = Atomids(np.array(atomids, dtype=np.int32))
    atomnames = Atomnames(np.array(atomnames, dtype=object))
    atomtypes = Atomtypes(np.array(atomtypes, dtype=object))
    charges = Charges(np.array(charges, dtype=np.float32))
    masses = Masses(np.array(masses, dtype=np.float32))

    moltypes = np.array(moltypes, dtype=object)
    molnums = np.array(molnums, dtype=np.int32)
    segids = np.array(segids, dtype=object)
    resids = np.array(resids, dtype=np.int32)
    resnames = np.array(resnames, dtype=object)
    (residx, new_resids,
     (new_resnames,
      new_moltypes,
      new_molnums,
      perres_segids
      )
     ) = squash_by(resids,
                   resnames,
                   moltypes,
                   molnums,
                   segids)
    residueids = Resids(new_resids)
    residuenames = Resnames(new_resnames)
    residue_moltypes = Moltypes(new_moltypes)
    residue_molnums = Molnums(new_molnums)

    segidx, perseg_segids = squash_by(perres_segids)[:2]
    segids = Segids(perseg_segids)

    top = Topology(len(atomids), len(new_resids), len(perseg_segids),
                   attrs=[atomids, atomnames, atomtypes,
                          charges, masses,
                          residueids, residuenames,
                          residue_moltypes, residue_molnums,
                          segids],
                   atom_resindex=residx,
                   residue_segindex=segidx)
    top.add_TopologyAttr(Bonds([bond for bond in bonds if bond]))
    top.add_TopologyAttr(Angles([angle for angle in angles if angle]))
    top.add_TopologyAttr(Dihedrals([dihedral for dihedral in dihedrals
                                    if dihedral]))
    top.add_TopologyAttr(Impropers([improper for improper in impropers
                                    if improper]))

    return top


def do_symstr(data, symtab):
    #do_symstr: get a string based on index from the symtab
    ndx = data.unpack_int()
    return symtab[ndx]


def do_symtab(data):
    symtab_nr = data.unpack_int()  # number of symbols
    symtab = []
    for i in range(symtab_nr):
        j = data.do_string()
        symtab.append(j)
    return symtab


def do_ffparams(data, fver):
    atnr = data.unpack_int()
    if fver < 57:
        data.unpack_int()  # idum
    ntypes = data.unpack_int()
    functype = ndo_int(data, ntypes)
    reppow = data.unpack_double() if fver >= 66 else 12.0
    if fver >= 57:
        fudgeQQ = data.unpack_real()

    # mimicing the c code,
    # remapping the functype due to inconsistency in different versions
    for i in range(len(functype)):
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
    data.unpack_real()  # rA
    data.unpack_real()  # krA
    data.unpack_real()  # rB
    data.unpack_real()  # krB


def do_iparams(data, functypes, fver):
    # Not all elif cases in this function has been used and tested
    for k, i in enumerate(functypes):
        if i in [
            S.F_ANGLES, S.F_G96ANGLES,
            S.F_BONDS, S.F_G96BONDS,
            S.F_HARMONIC, S.F_IDIHS
        ]:
            do_harm(data)
        elif i in [S.F_RESTRANGLES]:
            data.unpack_real()  # harmonic.rA
            data.unpack_real()  # harmonic.krA
        elif i in [S.F_LINEAR_ANGLES]:
            data.unpack_real()  # linangle.klinA
            data.unpack_real()  # linangle.aA
            data.unpack()  # linangle.klinB
            data.unpack()  # linangle.aB);
        elif i in [S.F_FENEBONDS]:
            data.unpack_real()  # fene.bm
            data.unpack_real()  # fene.kb
        elif i in [S.F_RESTRBONDS]:
            data.unpack_real()  # restraint.lowA
            data.unpack_real()  # restraint.up1A
            data.unpack_real()  # restraint.up2A
            data.unpack_real()  # restraint.kA
            data.unpack_real()  # restraint.lowB
            data.unpack_real()  # restraint.up1B
            data.unpack_real()  # restraint.up2B
            data.unpack_real()  # restraint.kB
        elif i in [S.F_TABBONDS, S.F_TABBONDSNC, S.F_TABANGLES, S.F_TABDIHS]:
            data.unpack_real()  # tab.kA
            data.unpack_int()  # tab.table
            data.unpack_real()  # tab.kB
        elif i in [S.F_CROSS_BOND_BONDS]:
            data.unpack_real()  # cross_bb.r1e
            data.unpack_real()  # cross_bb.r2e
            data.unpack_real()  # cross_bb.krr
        elif i in [S.F_CROSS_BOND_ANGLES]:
            data.unpack_real()  # cross_ba.r1e
            data.unpack_real()  # cross_ba.r2e
            data.unpack_real()  # cross_ba.r3e
            data.unpack_real()  # cross_ba.krt
        elif i in [S.F_UREY_BRADLEY]:
            data.unpack_real()  # u_b.theta
            data.unpack_real()  # u_b.ktheta
            data.unpack_real()  # u_b.r13
            data.unpack_real()  # u_b.kUB
            if fver >= 79:
                data.unpack_real()  # u_b.thetaB
                data.unpack_real()  # u_b.kthetaB
                data.unpack_real()  # u_b.r13B
                data.unpack_real()  # u_b.kUBB
        elif i in [S.F_QUARTIC_ANGLES]:
            data.unpack_real()  # qangle.theta
            ndo_real(data, 5)   # qangle.c
        elif i in [S.F_BHAM]:
            data.unpack_real()  # bham.a
            data.unpack_real()  # bham.b
            data.unpack_real()  # bham.c
        elif i in [S.F_MORSE]:
            data.unpack_real()  # morse.b0
            data.unpack_real()  # morse.cb
            data.unpack_real()  # morse.beta
            if fver >= 79:
                data.unpack_real()  # morse.b0B
                data.unpack_real()  # morse.cbB
                data.unpack_real()  # morse.betaB
        elif i in [S.F_CUBICBONDS]:
            data.unpack_real()  # cubic.b0
            data.unpack_real()  # cubic.kb
            data.unpack_real()  # cubic.kcub
        elif i in [S.F_CONNBONDS]:
            pass
        elif i in [S.F_POLARIZATION]:
            data.unpack_real()  # polarize.alpha
        elif i in [S.F_ANHARM_POL]:
            data.unpack_real()  # anharm_polarize.alpha
            data.unpack_real()  # anharm_polarize.drcut
            data.unpack_real()  # anharm_polarize.khyp
        elif i in [S.F_WATER_POL]:
            if fver < 31:
                fver_err(fver)
            data.unpack_real()  # wpol.al_x
            data.unpack_real()  # wpol.al_y
            data.unpack_real()  # wpol.al_z
            data.unpack_real()  # wpol.rOH
            data.unpack_real()  # wpol.rHH
            data.unpack_real()  # wpol.rOD
        elif i in [S.F_THOLE_POL]:
            data.unpack_real()  # thole.a
            data.unpack_real()  # thole.alpha1
            data.unpack_real()  # thole.alpha2
            data.unpack_real()  # thole.rfac

        elif i in [S.F_LJ]:
            data.unpack_real()  # lj_c6
            data.unpack_real()  # lj_c9
        elif i in [S.F_LJ14]:
            data.unpack_real()  # lj14_c6A
            data.unpack_real()  # lj14_c12A
            data.unpack_real()  # lj14_c6B
            data.unpack_real()  # lj14_c12B
        elif i in [S.F_LJC14_Q]:
            data.unpack_real()  # ljc14.fqq
            data.unpack_real()  # ljc14.qi
            data.unpack_real()  # ljc14.qj
            data.unpack_real()  # ljc14.c6
            data.unpack_real()  # ljc14.c12
        elif i in [S.F_LJC_PAIRS_NB]:
            data.unpack_real()  # ljcnb.qi
            data.unpack_real()  # ljcnb.qj
            data.unpack_real()  # ljcnb.c6
            data.unpack_real()  # ljcnb.c12

        elif i in [S.F_PIDIHS, S.F_ANGRES, S.F_ANGRESZ, S.F_PDIHS]:
            data.unpack_real()  # pdihs_phiA
            data.unpack_real()  # pdihs_cpA
            if (i == S.F_ANGRES or i == S.F_ANGRESZ) and fver < 42:
                data.unpack_real()  # harmonic.rB
                data.unpack_real()  # harmonic.krB
            else:
                data.unpack_real()  # pdihs_phiB
                data.unpack_real()  # pdihs_cpB
                data.unpack_int()  # pdihs_mult

        elif i in [S.F_RESTRDIHS]:
            data.unpack_real()  # pdihs.phiA
            data.unpack_real()  # pdihs.cpA
        elif i in [S.F_DISRES]:
            data.unpack_int()  # disres.label
            data.unpack_int()  # disres.type
            data.unpack_real()  # disres.low
            data.unpack_real()  # disres.up1
            data.unpack_real()  # disres.up2
            data.unpack_real()  # disres.kfac

        elif i in [S.F_ORIRES]:
            data.unpack_int()  # orires.ex
            data.unpack_int()  # orires.label
            data.unpack_int()  # orires.power
            data.unpack_real()  # orires.c
            data.unpack_real()  # orires.obs
            data.unpack_real()  # orires.kfac

        elif i in [S.F_DIHRES]:
            if fver < 72:
                data.unpack_int()  # idum
                data.unpack_int()  # idum
            data.unpack_real()  # dihres.phiA
            data.unpack_real()  # dihres.dphiA
            data.unpack_real()  # dihres.kfacA
            if fver >= 72:
                data.unpack_real()  # dihres.phiB
                data.unpack_real()  # dihres.dphiB
                data.unpack_real()  # dihres.kfacB

        elif i in [S.F_POSRES]:
            do_rvec(data)  # posres.pos0A
            do_rvec(data)  # posres.fcA
            if fver < 27:
                fver_err(fver)
            else:
                do_rvec(data)  # posres.pos0B
                do_rvec(data)  # posres.fcB

        elif i in [S.F_FBPOSRES]:
            data.unpack_int()   # fbposres.geom
            do_rvec(data)       # fbposres.pos0
            data.unpack_real()  # fbposres.r
            data.unpack_real()  # fbposres.k

        elif i in [S.F_CBTDIHS]:
            ndo_real(data, S.NR_CBTDIHS)  # cbtdihs.cbtcA

        elif i in [S.F_RBDIHS]:
            ndo_real(data, S.NR_RBDIHS)  # iparams_rbdihs_rbcA
            if fver >= 25:
                ndo_real(data, S.NR_RBDIHS)  # iparams_rbdihs_rbcB

        elif i in [S.F_FOURDIHS]:
            # Fourier dihedrals
            ndo_real(data, S.NR_RBDIHS)  # rbdihs.rbcA
            ndo_real(data, S.NR_RBDIHS)  # rbdihs.rbcB

        elif i in [S.F_CONSTR, S.F_CONSTRNC]:
            data.unpack_real()  # dA
            data.unpack_real()  # dB

        elif i in [S.F_SETTLE]:
            data.unpack_real()  # settle.doh
            data.unpack_real()  # settle.dhh

        elif i in [S.F_VSITE2, S.F_VSITE2FD]:
            data.unpack_real()  # vsite.a

        elif i in [S.F_VSITE3, S.F_VSITE3FD, S.F_VSITE3FAD]:
            data.unpack_real()  # vsite.a
            data.unpack_real()  # vsite.b

        elif i in [S.F_VSITE3OUT, S.F_VSITE4FD, S.F_VSITE4FDN]:
            data.unpack_real()  # vsite.a
            data.unpack_real()  # vsite.b
            data.unpack_real()  # vsite.c

        elif i in [S.F_VSITEN]:
            data.unpack_int()  # vsiten.n
            data.unpack_real()  # vsiten.a

        elif i in [S.F_GB12, S.F_GB13, S.F_GB14]:
            # /* We got rid of some parameters in version 68 */
            if fver < 68:
                data.unpack_real()  # rdum
                data.unpack_real()  # rdum
                data.unpack_real()  # rdum
                data.unpack_real()  # rdum
            data.unpack_real()  # gb.sar
            data.unpack_real()  # gb.st
            data.unpack_real()  # gb.pi
            data.unpack_real()  # gb.gbr
            data.unpack_real()  # gb.bmlt

        elif i in [S.F_CMAP]:
            data.unpack_int()  # cmap.cmapA
            data.unpack_int()  # cmap.cmapB
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
    bonds, angles, dihs, impr = [], [], [], []
    for ilist in ilists:
        if ilist.nr > 0:
            ik_obj = obj.InteractionKind(*ilist.ik)
            ias = ilist.iatoms

            # the following if..elif..else statement needs to be updated as new
            # type of interactions become interested
            if ik_obj.name in ['BONDS', 'G96BONDS', 'MORSE', 'CUBICBONDS',
                               'CONNBONDS', 'HARMONIC', 'FENEBONDS',
                               'RESTRAINTPOT', 'CONSTR', 'CONSTRNC',
                               'TABBONDS', 'TABBONDSNC']:
                bonds += list(ik_obj.process(ias))
            elif ik_obj.name in ['ANGLES', 'G96ANGLES', 'CROSS_BOND_BOND',
                                 'CROSS_BOND_ANGLE', 'UREY_BRADLEY', 'QANGLES',
                                 'RESTRANGLES', 'TABANGLES']:
                angles += list(ik_obj.process(ias))
            elif ik_obj.name in ['PDIHS', 'RBDIHS', 'RESTRDIHS', 'CBTDIHS',
                                 'FOURDIHS', 'TABDIHS']:
                dihs += list(ik_obj.process(ias))
            elif ik_obj.name in ['IDIHS', 'PIDIHS']:
                impr += list(ik_obj.process(ias))
            elif ik_obj.name == 'SETTLE':
                # SETTLE interactions are optimized triangular constraints for
                # water molecules. They should be counted as a pair of bonds
                # between the oxigen and the hydrogens. In older versions of
                # the TPR format only specifies the index of the oxygen and
                # assumes that the next two atoms are the hydrogens.
                if len(ias) == 2:
                    # Old format. Only the first atom is specified.
                    base_atom = ias[1]
                    bonds += [
                        [base_atom, base_atom + 1],
                        [base_atom, base_atom + 2],
                    ]
                else:
                    all_settle = ik_obj.process(ias)
                    for settle in all_settle:
                        base_atom = settle[0]
                        bonds += [
                            [settle[0], settle[1]],
                            [settle[0], settle[2]],
                        ]
            else:
                # other interaction types are not interested at the moment
                pass

    bonds = bonds if bonds else None
    angles = angles if angles else None
    dihs = dihs if dihs else None
    impr = impr if impr else None
    moltype = obj.MoleculeKind(molname, atomkinds, bonds, angles, dihs, impr)
    #### end: MDAnalysis specific

    # info in do_block and do_blocka is not interested, but has to be parsed
    # here so that all moltypes can be parsed properly
    do_block(data)
    do_blocka(data)
    return moltype


def do_atoms(data, symtab, fver):
    nr = data.unpack_int()  # number of atoms in a particular molecule
    nres = data.unpack_int()  # number of residues in a particular molecule

    if fver < 57:
        fver_err(fver)

    atoms = []
    for i in range(nr):
        A = do_atom(data, fver)
        atoms.append(A)

    # do_strstr
    atomnames = [symtab[i] for i in ndo_int(data, nr)]

    if fver <= 20:
        fver_err(fver)
    else:
        type = [symtab[i] for i in ndo_int(data, nr)]  # e.g. opls_111
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
        for i in range(nres):
            resnames.append(symtab[data.unpack_int()])
            data.unpack_int()
            data.unpack_uchar()
    return resnames


def do_atom(data, fver):
    m = data.unpack_real()  # mass
    q = data.unpack_real()  # charge
    mB = data.unpack_real()
    qB = data.unpack_real()
    tp = data.unpack_ushort()  # type is a keyword in python
    typeB = data.unpack_ushort()
    ptype = data.unpack_int()  # regular atom, virtual site or others
    resind = data.unpack_int()  # index of residue

    if fver >= 52:
        atomnumber = data.unpack_int()  # index of atom type

    if fver < 23 or fver < 39 or fver < 57:
        fver_err(fver)
    return obj.Atom(m, q, mB, qB, tp, typeB, ptype, resind, atomnumber)


def do_ilists(data, fver):
    nr = []  # number of ilist
    iatoms = []  # atoms involved in a particular interaction type
    pos = []
    for j in range(S.F_NRE):  # total number of energies (i.e. interaction types)
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
            for i in range(n):
                l_.append(data.unpack_int())
            iatoms.append(l_)

    return [obj.Ilist(n, it, i) for n, it, i in zip(nr, S.interaction_types, iatoms)]


def do_molblock(data):
    molb_type = data.unpack_int()
    molb_nmol = data.unpack_int()  # number of moles in the molblock
    molb_natoms_mol = data.unpack_int()  # number of atoms in a molecule
    molb_nposres_xA = data.unpack_int()
    if molb_nposres_xA > 0:
        ndo_rvec(data, molb_nposres_xA)
    molb_nposres_xB = data.unpack_int()  # The number of posres coords for top B
    if molb_nposres_xB > 0:
        ndo_rvec(data, molb_nposres_xB)

    return obj.Molblock(molb_type, molb_nmol, molb_natoms_mol,
                        molb_nposres_xA, molb_nposres_xB)


def do_block(data):
    block_nr = data.unpack_int()  # for cgs: charge groups
    # starting or ending atom indices, based on which cgs are grouped
    ndo_int(data, block_nr + 1)
    return do_block


def do_blocka(data):
    block_nr = data.unpack_int()  # No. of atoms with excls
    block_nra = data.unpack_int()  # total times fo appearance of atoms for excls
    ndo_int(data, block_nr + 1)
    ndo_int(data, block_nra)
    return block_nr, block_nra


##############UTILS FOR INFORMATION NOT INTERESTED AT THE MOMENT###############

def do_grps(data):  # pragma: no cover
    grps_nr = []
    myngrps = ngrps = S.egcNR  # remind of version inconsistency
    for j in range(ngrps):
        if j < myngrps:
            v = data.unpack_int()
            grps_nr.append(v)
        ndo_int(data, v)
    return grps_nr


def do_groups(data, symtab):  # pragma: no cover
    do_grps(data)

    ngrpname = data.unpack_int()
    # do_strstr, list of indices of group name: e.g. System, Protein,
    # Protein-H, etc. use symtab[i] to check
    ndo_int(data, ngrpname)

    ngrpnr = []
    grpnr = []
    for i in range(S.egcNR):
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


def do_atomtypes(data):  # pragma: no cover
    at_nr = data.unpack_int()  # at: atomtype
    at_radius = ndo_real(data, at_nr)
    at_vol = ndo_real(data, at_nr)
    at_surftens = ndo_real(data, at_nr)
    at_atomnumber = ndo_int(data, at_nr)
    return at_radius, at_vol, at_surftens, at_atomnumber


def do_inputrec(data):  # pragma: no cover
    """Read through header information from TPR file *data* structure.

    Note that this function does not return any useful data itself. If
    your are interested in using the header information, use this
    functions as a starting point for your own code.
    """
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

    ndo_rvec(data, DIM)  # ir_deform

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

    ndo_real(data, ir_opts_ngtc)  # ir_nrdf
    ndo_real(data, ir_opts_ngtc)  # ir_ref_t
    ndo_real(data, ir_opts_ngtc)  # ir_tau_t

    if ir_opts_ngfrz > 0:
        ndo_ivec(data, ir_opts_ngfrz)  # ir_opts_nFreeze

    if ir_opts_ngacc > 0:
        ndo_rvec(data, ir_opts_ngacc)  # ir_opts_acc

    ndo_int(data, ir_opts_ngener ** 2)  # ir_opts_egp_flags
    ndo_int(data, ir_opts_ngtc)  # ir_opts_annealing
    ir_opts_anneal_npoints = ndo_int(data, ir_opts_ngtc)

    ir_opts_anneal_time = []
    ir_opts_anneal_temp = []
    for j in range(ir_opts_ngtc):
        k = ir_opts_anneal_npoints[j]
        ir_opts_anneal_time.append(ndo_int(data, k))
        ir_opts_anneal_temp.append(ndo_int(data, k))

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

        ir_ex_a.append(ndo_real(data, x))
        ir_ex_phi.append(ndo_real(data, x))
        ir_et_a.append(ndo_real(data, y))
        ir_et_phi.append(ndo_real(data, y))

    # QMM stuff
    data.unpack_bool()  # ir_bQMMM
    data.unpack_int()  # ir_bQMMMscheme
    data.unpack_real()  # ir_scalefactor
    data.unpack_int()  # ir_opts_ngQM

    # if ir_opts_ngQM > 0:
    # do_something

    # indicating the parsing finishes properly
    data.done()
