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

Wherever err(fver) is used in this and other tpr_*.py files, it means the tpx
version problem haven't be resolved for those other than 58 and 73 (or gromacs
version before 4.x)
"""

import sys
import xdrlib

import tpr_utils as U

ndo_int = U.ndo_int
ndo_real = U.ndo_real
ndo_rvec = U.ndo_rvec
ndo_ivec = U.ndo_ivec

def print_header(th):
    print "#" * 79
    print "Gromacs version   : {0}".format(th.ver_str)
    print "tpx     version   : {0}".format(th.fver)
    print "tpx     generation: {0}".format(th.fgen)
    print "#" * 79

def parse(tprfile):
    tprf = open(tprfile).read()
    data = xdrlib.Unpacker(tprf)
    try:
        th = U.read_tpxheader(data)                    # tpxheader
    except EOFError:
        print "{0}\nInvalid tpr file or cannot be recognized\n".format(tprfile)
        raise

    print_header(th)

    V = th.fver                                    # since it's used very often

    state_ngtc = th.ngtc         # done init_state() in src/gmxlib/tpxio.c
    if th.bBox:
        U.extract_box_info(data, V)

    if state_ngtc > 0 and V >= 28:
        if V < 69:                      # redundancy due to  different versions
            ndo_real(data, state_ngtc)
        ndo_real(data, state_ngtc)        # relevant to Berendsen tcoupl_lambda

    if V < 26:
        U.err(V)

    if th.bTop:
        tpr_top = U.do_mtop(data, V)
        structure = {
            '_atoms': tpr_top.atoms,
            '_bonds': tpr_top.bonds,
            '_angles': tpr_top.angles,
            '_dihe': tpr_top.dihe,
            '_impr': tpr_top.impr
            }
    else:
        raise ValueError("No topology found in tpr file: {0}".formation(tprfile))

    # THE FOLLOWING CODE IS WORKING FOR TPX VERSION 58, BUT SINCE THESE INFO IS
    # NOT INTERESTED, SO IT'S NOT COVERED IN ALL VERSIONS. PARSING STOPS HERE.

    # if th.bX:
    #     ndo_rvec(data, th.natoms)

    # if th.bV:
    #     ndo_rvec(data, th.natoms)

    # if th.bF:
    #     ndo_rvec(data, th.natoms)

    # not useful at the moment
    # ePBC = -1;
    # bPeriodicMols = False
    # if th.bIr:
    #     # update
    #     data.unpack_int()                                # ePBC
    #     data.unpack_bool()                               # bPeriodicMols
    #     # 17 < 23. and ir (ir is from the c code, seems not apply here
    #     if th.fgen < setting.tpx_generation:
    #         # a crazily long (670 lines) function in c, slightly better here
    #         # (240 lines), so put it in setting.py
    #         setting.do_inputrec(data)

    return structure

if __name__ == "__main__":
    try:
        parse(sys.argv[1])
    except EOFError:
        print "{0}\nInvalid tpr file or cannot be recognized\n".format(sys.argv[1])
    except IndexError:
        print "Please feed an tpr file or use this file as a module"

    # parse("../../../testsuite/MDAnalysisTests/data/adk_oplsaa.tpr")
    # parse("../../../testsuite/MDAnalysisTests/data/gmxv73.tpr")
