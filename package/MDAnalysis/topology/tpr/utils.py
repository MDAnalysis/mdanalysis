# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding: utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
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
# Released under the  Lesser GNU Public Licence, v2.1+

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

import numpy as np
from mda_xdrlib import xdrlib
import struct

from . import obj
from . import setting
from ..base import squash_by
from ...core.topology import Topology
from ...core.topologyattrs import (
    Atomids,
    Atomnames,
    Atomtypes,
    Masses,
    Charges,
    Elements,
    Resids,
    Resnames,
    Moltypes,
    Molnums,
    Segids,
    ChainIDs,
    Bonds,
    Angles,
    Dihedrals,
    Impropers,
)


class TPXUnpacker(xdrlib.Unpacker):
    """
    Extend the standard XDR unpacker for the specificity of TPX files.
    """
    def __init__(self, data):
        super().__init__(data)
        self._buf = self.get_buffer()

    # The parent class uses a dunder attribute to store the
    # cursor position. This property makes it easier to manipulate
    # this attribute that is otherwise "protected".
    @property
    def _pos(self):
        return self.get_position()

    @_pos.setter
    def _pos(self, value):
        self.set_position(value)

    def _unpack_value(self, item_size, struct_template):
        start_position = self._pos
        end_position = self._pos = start_position + item_size
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
        new_unpacker = cls(unpacker.get_buffer())
        new_unpacker._pos = unpacker.get_position()
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
        start_position = self._pos
        end_position = self._pos = start_position + n
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
    return data.unpack_farray(setting.DIM, data.unpack_real)


def ndo_rvec(data, n):
    """mimic of gmx_fio_ndo_rvec in gromacs"""
    return [data.unpack_farray(setting.DIM, data.unpack_real) for i in range(n)]


def ndo_ivec(data, n):
    """mimic of gmx_fio_ndo_rvec in gromacs"""
    return [data.unpack_farray(setting.DIM, data.unpack_int) for i in range(n)]


def fileVersion_err(fver):
    if fver not in setting.SUPPORTED_VERSIONS:
        raise NotImplementedError(
            f"Your tpx version is {fver}, which this parser does not support, yet "
        )


def define_unpack_real(prec, data):
    """Define an unpack_real method of data based on the float precision used"""
    if prec == 4:
        data.unpack_real = data.unpack_float
    elif prec == 8:
        data.unpack_real = data.unpack_double
    else:
        raise ValueError(f"unsupported precision: {prec}")


def read_tpxheader(data):
    """this function is now compatible with do_tpxheader in tpxio.cpp
    """
    # Last compatibility check with gromacs-2016
    ver_str = data.do_string()  # version string e.g. VERSION 4.0.5
    if not ver_str.startswith(b'VERSION'):
        raise ValueError('Input does not look like a TPR file.')
    precision = data.unpack_int()  # e.g. 4
    define_unpack_real(precision, data)
    fileVersion = data.unpack_int()  # version of tpx file
    fileVersion_err(fileVersion)

    # This is for backward compatibility with development version 77-79 where
    # the tag was, mistakenly, placed before the generation.
    # These versions are development versions between the 4.5 and 4.6 series.
    if 77 <= fileVersion <= 79:
        data.unpack_int()  # the value is 8, but haven't found the
        file_tag = data.do_string()

    fileGeneration = data.unpack_int() if fileVersion >= 26 else 0  # generation of tpx file, e.g. 17

    # Versions before 77 don't have the tag, set it to TPX_TAG_RELEASE file_tag
    # file_tag is used for comparing with tpx_tag. Only tpr files with a
    # tpx_tag from a lower or the same version of gromacs code can be parsed by
    # the tpxio.c

    file_tag = data.do_string() if fileVersion >= 81 else setting.TPX_TAG_RELEASE

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
    if fileVersion >= setting.tpxv_AddSizeField and fileGeneration >= 27:
        sizeOfTprBody = data.unpack_int64()

    th = obj.TpxHeader(ver_str, precision, fileVersion, fileGeneration,
                       file_tag, natoms, ngtc, fep_state, lamb,
                       bIr, bTop, bX, bV, bF, bBox, sizeOfTprBody)
    return th


def extract_box_info(data, fver):
    box = ndo_rvec(data, setting.DIM)
    box_rel = ndo_rvec(data, setting.DIM)
    box_v = ndo_rvec(data, setting.DIM)

    return obj.Box(box, box_rel, box_v)


def do_mtop(data, fver, tpr_resid_from_one=False, precision=4):
    # mtop: the topology of the whole system
    symtab = do_symtab(data)
    do_symstr(data, symtab)  # system_name
    ff_params = do_ffparams(data, fver)  # params

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
    chainIDs = []
    resids = []
    resnames = []
    atomnames = []
    atomtypes = []
    moltypes = []
    molnums = []
    charges = []
    masses = []
    elements = []

    atom_start_ndx = 0
    res_start_ndx = 0
    molnum = 0
    for i in range(mtop.nmolblock):
        # molb_type is just an index for moltypes/molecule_types
        mb = do_molblock(data)
        # segment is made to correspond to the molblock as in gromacs, the
        # naming is kind of arbitrary
        molblock = mtop.moltypes[mb.molb_type].name.decode('utf-8')
        segid = f"seg_{i}_{molblock}"
        chainID = molblock[14:] if molblock[:14] == "Protein_chain_" else molblock
        for j in range(mb.molb_nmol):
            mt = mtop.moltypes[mb.molb_type]  # mt: molecule type
            for atomkind in mt.atomkinds:
                atomids.append(atomkind.id + atom_start_ndx)
                segids.append(segid)
                chainIDs.append(chainID)
                resids.append(atomkind.resid + res_start_ndx)
                resnames.append(atomkind.resname.decode())
                atomnames.append(atomkind.name.decode())
                atomtypes.append(atomkind.type.decode())
                moltypes.append(molblock)
                molnums.append(molnum)
                charges.append(atomkind.charge)
                masses.append(atomkind.mass)
                elements.append(atomkind.element_symbol)
            molnum += 1

            # remap_ method returns [blah, blah, ..] or []
            bonds.extend(mt.remap_bonds(atom_start_ndx))
            angles.extend(mt.remap_angles(atom_start_ndx))
            dihedrals.extend(mt.remap_dihe(atom_start_ndx))
            impropers.extend(mt.remap_impr(atom_start_ndx))

            atom_start_ndx += mt.number_of_atoms()
            res_start_ndx += mt.number_of_residues()

    atomids = Atomids(np.array(atomids, dtype=np.int32))
    atomnames = Atomnames(np.array(atomnames, dtype=object))
    atomtypes = Atomtypes(np.array(atomtypes, dtype=object))
    charges = Charges(np.array(charges, dtype=np.float32))
    masses = Masses(np.array(masses, dtype=np.float32))

    moltypes = np.array(moltypes, dtype=object)
    molnums = np.array(molnums, dtype=np.int32)
    segids = np.array(segids, dtype=object)
    chainIDs = np.array(chainIDs, dtype=object)
    resids = np.array(resids, dtype=np.int32)
    if tpr_resid_from_one:
        resids += 1

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
    chainIDs = ChainIDs(chainIDs)

    top = Topology(
        len(atomids),
        len(new_resids),
        len(perseg_segids),
        attrs=[
            atomids,
            atomnames,
            atomtypes,
            charges,
            masses,
            residueids,
            residuenames,
            residue_moltypes,
            residue_molnums,
            segids,
            chainIDs,
        ],
        atom_resindex=residx,
        residue_segindex=segidx,
    )
    top.add_TopologyAttr(Bonds([bond for bond in bonds if bond]))
    top.add_TopologyAttr(Angles([angle for angle in angles if angle]))
    top.add_TopologyAttr(Dihedrals([dihedral for dihedral in dihedrals
                                    if dihedral]))
    top.add_TopologyAttr(Impropers([improper for improper in impropers
                                    if improper]))

    if any(elements):
        elements = Elements(np.array(elements, dtype=object))
        top.add_TopologyAttr(elements)

    # NOTE: the tpr striding code below serves the
    # purpose of placing us in a suitable "seek" position
    # in the binary file such that we are well placed
    # for reading coords and velocities if they are present
    # the order of operations is based on an analysis of
    # the C++ code in do_mtop() function in the GROMACS
    # source at:
    # src/gromacs/fileio/tpxio.cpp
    # TODO: expand tpx version support for striding to
    # the coordinates
    atnr = ff_params.atnr
    if fver >= 58:
        # TODO: the following value is important, and not sure
        # how to access programmatically yet...
        # from GMX source code:
        # api/legacy/include/gromacs/topology/topology_enums.h
        # worst case scenario we hard code it based on
        # tpx/GMX version?
        # for details, see gh-4873
        SimulationAtomGroupType_size = 10
        n_atoms = data.unpack_int()
        if fver < 116:
            for i in range(3 * atnr * int(precision / 4)):
                data.unpack_int()
        if fver < 129:
            # NOTE: speculative, Tyler did this by
            # inspecting binary data and relative file offsets
            # for older tpr files relative to more recent;
            # the stride skip appears related to `atnr` in GMX source, the
            # number of non-bonded atom types
            for i in range(atnr + 1):
                data.unpack_int()

        if 58 < fver <= 83:
            for i in range(atnr * int(precision / 4) * 2):
                data.unpack_int()
        if fver > 83:
            interm = data.unpack_uchar()
        if 83 < fver < 116:
            for i in range(2 * atnr):
                data.unpack_int()
        if fver > 58:
            ngrid = data.unpack_int()
            grid_spacing = data.unpack_int()
            n_elements = grid_spacing ** 2
            for i in range(ngrid):
                for j in range(n_elements):
                    ndo_real(data, 4)
        for i in range(SimulationAtomGroupType_size):
            group_size = data.unpack_int()
            ndo_int(data, group_size)
        n_group_names = data.unpack_int()
        for i in range(n_group_names):
            data.unpack_int()
        for i in range(SimulationAtomGroupType_size):
            n_grp_numbers = data.unpack_int()
            if n_grp_numbers != 0:
                for i in range(n_grp_numbers):
                    data.unpack_uchar()
        if fver > 119:
            im_excl_grp_size = data.unpack_int()
            ndo_int(data, im_excl_grp_size)
            # TODO: why is this needed?
            # NOTE: this `mystery_n` value is at least sometimes used
            # for tpx >= 137, but is often 0 otherwise
            mystery_n = data.unpack_int()
        if fver >= 137:
            # NOTE: Tyler observed that the first float32
            # positional coordinate was offset by `mystery_n`
            # 32-bit units in a sample tpx 137 tpr file
            # by empirical testing
            for i in range(mystery_n):
                data.unpack_int()

    return top


def do_symstr(data, symtab):
    # do_symstr: get a string based on index from the symtab
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
    ntypes = data.unpack_int()
    functype = ndo_int(data, ntypes)
    reppow = data.unpack_double() if fver >= 66 else 12.0
    fudgeQQ = data.unpack_real()

    # mimicing the c code,
    # remapping the functype due to inconsistency in different versions
    for i in range(len(functype)):
        for k in setting.ftupd:
            # j[0]: tpx_version, j[1] funtype
            if fver < k[0] and functype[i] >= k[1]:
                functype[i] += 1

    do_iparams(data, functype, fver)

    params = obj.Params(atnr, ntypes, functype, reppow, fudgeQQ)
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
            setting.F_ANGLES, setting.F_G96ANGLES,
            setting.F_BONDS, setting.F_G96BONDS,
            setting.F_HARMONIC, setting.F_IDIHS
        ]:
            do_harm(data)
        elif i in [setting.F_RESTRANGLES]:
            data.unpack_real()  # harmonic.rA
            data.unpack_real()  # harmonic.krA
            if fver >= 134:
                data.unpack_real()  # harmonic.rB
                data.unpack_real()  # harmonic.krB

        elif i in [setting.F_LINEAR_ANGLES]:
            data.unpack_real()  # linangle.klinA
            data.unpack_real()  # linangle.aA
            data.unpack_real()  # linangle.klinB
            data.unpack_real()  # linangle.aB);
        elif i in [setting.F_FENEBONDS]:
            data.unpack_real()  # fene.bm
            data.unpack_real()  # fene.kb
        elif i in [setting.F_RESTRBONDS]:
            data.unpack_real()  # restraint.lowA
            data.unpack_real()  # restraint.up1A
            data.unpack_real()  # restraint.up2A
            data.unpack_real()  # restraint.kA
            data.unpack_real()  # restraint.lowB
            data.unpack_real()  # restraint.up1B
            data.unpack_real()  # restraint.up2B
            data.unpack_real()  # restraint.kB
        elif i in [
            setting.F_TABBONDS, setting.F_TABBONDSNC,
            setting.F_TABANGLES, setting.F_TABDIHS,
        ]:
            data.unpack_real()  # tab.kA
            data.unpack_int()  # tab.table
            data.unpack_real()  # tab.kB
        elif i in [setting.F_CROSS_BOND_BONDS]:
            data.unpack_real()  # cross_bb.r1e
            data.unpack_real()  # cross_bb.r2e
            data.unpack_real()  # cross_bb.krr
        elif i in [setting.F_CROSS_BOND_ANGLES]:
            data.unpack_real()  # cross_ba.r1e
            data.unpack_real()  # cross_ba.r2e
            data.unpack_real()  # cross_ba.r3e
            data.unpack_real()  # cross_ba.krt
        elif i in [setting.F_UREY_BRADLEY]:
            data.unpack_real()  # u_b.theta
            data.unpack_real()  # u_b.ktheta
            data.unpack_real()  # u_b.r13
            data.unpack_real()  # u_b.kUB
            if fver >= 79:
                data.unpack_real()  # u_b.thetaB
                data.unpack_real()  # u_b.kthetaB
                data.unpack_real()  # u_b.r13B
                data.unpack_real()  # u_b.kUBB
        elif i in [setting.F_QUARTIC_ANGLES]:
            data.unpack_real()  # qangle.theta
            ndo_real(data, 5)   # qangle.c
        elif i in [setting.F_BHAM]:
            data.unpack_real()  # bham.a
            data.unpack_real()  # bham.b
            data.unpack_real()  # bham.c
        elif i in [setting.F_MORSE]:
            data.unpack_real()  # morse.b0
            data.unpack_real()  # morse.cb
            data.unpack_real()  # morse.beta
            if fver >= 79:
                data.unpack_real()  # morse.b0B
                data.unpack_real()  # morse.cbB
                data.unpack_real()  # morse.betaB
        elif i in [setting.F_CUBICBONDS]:
            data.unpack_real()  # cubic.b0g
            data.unpack_real()  # cubic.kb
            data.unpack_real()  # cubic.kcub
        elif i in [setting.F_CONNBONDS]:
            pass
        elif i in [setting.F_POLARIZATION]:
            data.unpack_real()  # polarize.alpha
        elif i in [setting.F_ANHARM_POL]:
            data.unpack_real()  # anharm_polarize.alpha
            data.unpack_real()  # anharm_polarize.drcut
            data.unpack_real()  # anharm_polarize.khyp
        elif i in [setting.F_WATER_POL]:
            data.unpack_real()  # wpol.al_x
            data.unpack_real()  # wpol.al_y
            data.unpack_real()  # wpol.al_z
            data.unpack_real()  # wpol.rOH
            data.unpack_real()  # wpol.rHH
            data.unpack_real()  # wpol.rOD
        elif i in [setting.F_THOLE_POL]:
            data.unpack_real()  # thole.a
            data.unpack_real()  # thole.alpha1
            data.unpack_real()  # thole.alpha2
            if fver < setting.tpxv_RemoveTholeRfac:
                data.unpack_real()  # thole.rfac

        elif i in [setting.F_LJ]:
            data.unpack_real()  # lj_c6
            data.unpack_real()  # lj_c9
        elif i in [setting.F_LJ14]:
            data.unpack_real()  # lj14_c6A
            data.unpack_real()  # lj14_c12A
            data.unpack_real()  # lj14_c6B
            data.unpack_real()  # lj14_c12B
        elif i in [setting.F_LJC14_Q]:
            data.unpack_real()  # ljc14.fqq
            data.unpack_real()  # ljc14.qi
            data.unpack_real()  # ljc14.qj
            data.unpack_real()  # ljc14.c6
            data.unpack_real()  # ljc14.c12
        elif i in [setting.F_LJC_PAIRS_NB]:
            data.unpack_real()  # ljcnb.qi
            data.unpack_real()  # ljcnb.qj
            data.unpack_real()  # ljcnb.c6
            data.unpack_real()  # ljcnb.c12

        elif i in [
            setting.F_PIDIHS, setting.F_ANGRES,
            setting.F_ANGRESZ, setting.F_PDIHS,
        ]:
            data.unpack_real()  # pdihs_phiA
            data.unpack_real()  # pdihs_cpA
            data.unpack_real()  # pdihs_phiB
            data.unpack_real()  # pdihs_cpB
            data.unpack_int()  # pdihs_mult

        elif i in [setting.F_RESTRDIHS]:
            data.unpack_real()  # pdihs.phiA
            data.unpack_real()  # pdihs.cpA
            if fver >= 134:
                data.unpack_real()  # pdihs.phiB
                data.unpack_real()  # pdihs.cpB
        elif i in [setting.F_DISRES]:
            data.unpack_int()  # disres.label
            data.unpack_int()  # disres.type
            data.unpack_real()  # disres.low
            data.unpack_real()  # disres.up1
            data.unpack_real()  # disres.up2
            data.unpack_real()  # disres.kfac

        elif i in [setting.F_ORIRES]:
            data.unpack_int()  # orires.ex
            data.unpack_int()  # orires.label
            data.unpack_int()  # orires.power
            data.unpack_real()  # orires.c
            data.unpack_real()  # orires.obs
            data.unpack_real()  # orires.kfac

        elif i in [setting.F_DIHRES]:
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

        elif i in [setting.F_POSRES]:
            do_rvec(data)  # posres.pos0A
            do_rvec(data)  # posres.fcA
            do_rvec(data)  # posres.pos0B
            do_rvec(data)  # posres.fcB

        elif i in [setting.F_FBPOSRES]:
            data.unpack_int()   # fbposres.geom
            do_rvec(data)       # fbposres.pos0
            data.unpack_real()  # fbposres.r
            data.unpack_real()  # fbposres.k

        elif i in [setting.F_CBTDIHS]:
            ndo_real(data, setting.NR_CBTDIHS)  # cbtdihs.cbtcA
            if fver >= 134:
                ndo_real(data, setting.NR_CBTDIHS)  # cbtdihs.cbtcB

        elif i in [setting.F_RBDIHS]:
            ndo_real(data, setting.NR_RBDIHS)  # iparams_rbdihs_rbcA
            ndo_real(data, setting.NR_RBDIHS)  # iparams_rbdihs_rbcB

        elif i in [setting.F_FOURDIHS]:
            # Fourier dihedrals
            ndo_real(data, setting.NR_RBDIHS)  # rbdihs.rbcA
            ndo_real(data, setting.NR_RBDIHS)  # rbdihs.rbcB

        elif i in [setting.F_CONSTR, setting.F_CONSTRNC]:
            data.unpack_real()  # dA
            data.unpack_real()  # dB

        elif i in [setting.F_SETTLE]:
            data.unpack_real()  # settle.doh
            data.unpack_real()  # settle.dhh

        elif i in [setting.F_VSITE1]:
            pass

        elif i in [setting.F_VSITE2, setting.F_VSITE2FD]:
            data.unpack_real()  # vsite.a

        elif i in [setting.F_VSITE3, setting.F_VSITE3FD, setting.F_VSITE3FAD]:
            data.unpack_real()  # vsite.a
            data.unpack_real()  # vsite.b

        elif i in [setting.F_VSITE3OUT, setting.F_VSITE4FD, setting.F_VSITE4FDN]:
            data.unpack_real()  # vsite.a
            data.unpack_real()  # vsite.b
            data.unpack_real()  # vsite.c

        elif i in [setting.F_VSITEN]:
            data.unpack_int()  # vsiten.n
            data.unpack_real()  # vsiten.a

        elif i in [setting.F_GB12, setting.F_GB13, setting.F_GB14]:
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

        elif i in [setting.F_CMAP]:
            data.unpack_int()  # cmap.cmapA
            data.unpack_int()  # cmap.cmapB
        elif i in [setting.F_ENNPOT]:
            # TODO: handle the new neural network potentials
            # supported from GROMACS 2025.0
            pass
        else:
            raise NotImplementedError(f"unknown functype: {i}")
    return


def do_moltype(data, symtab, fver):
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
            a.q,
            a.atomnumber,
        ))
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
                # between the oxygen and the hydrogens. In older versions of
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

    # info in do_block and do_blocka is not interesting, but has to be parsed
    # here so that all moltypes can be parsed properly
    do_block(data)
    do_blocka(data)
    return moltype


def do_atoms(data, symtab, fver):
    nr = data.unpack_int()  # number of atoms in a particular molecule
    nres = data.unpack_int()  # number of residues in a particular molecule

    atoms = []
    for i in range(nr):
        A = do_atom(data, fver)
        atoms.append(A)

    # do_strstr
    atomnames = [symtab[i] for i in ndo_int(data, nr)]

    type = [symtab[i] for i in ndo_int(data, nr)]  # e.g. opls_111
    typeB = [symtab[i] for i in ndo_int(data, nr)]
    resnames = do_resinfo(data, symtab, fver, nres)

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

    atomnumber = data.unpack_int()  # index of atom type

    return obj.Atom(m, q, mB, qB, tp, typeB, ptype, resind, atomnumber)


def do_ilists(data, fver):
    nr = []  # number of ilist
    iatoms = []  # atoms involved in a particular interaction type
    pos = []
    for j in range(setting.F_NRE):  # total number of energies (i.e. interaction types)
        bClear = False
        for k in setting.ftupd:
            if fver < k[0] and j == k[1]:
                bClear = True
        if bClear:
            nr.append(0)
            iatoms.append(None)
        else:
            # do_ilist
            n = data.unpack_int()
            nr.append(n)
            l_ = []
            for i in range(n):
                l_.append(data.unpack_int())
            iatoms.append(l_)

    return [
        obj.Ilist(n, it, i)
        for n, it, i in zip(nr, setting.interaction_types, iatoms)
    ]


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
