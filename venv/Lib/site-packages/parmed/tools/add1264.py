"""
This module contains the necessary machinery to add the LENNARD_JONES_CCOEF to a
topology file given a list of atomic polarizabilities
"""
from .exceptions import LJ12_6_4Error, DuplicateParamWarning
import warnings

WATER_POL = 1.444 # Polarizability of water

DEFAULT_C4_PARAMS = {
        'TIP3P' : {'Li1': 27.0, 'Na1': 0.0, 'K1': 8.0, 'Rb1': 0.0, 'Cs1': 2.0,
                   'Tl1': 50.0, 'Cu1': 7.0, 'Ag1': 83.0, 'F-1': -27.0,
                   'Cl-1': -38.0, 'Br-1': -39.0, 'I-1': -45.0, 'Be2': 186.5,
                   'Cu2': 290.9, 'Ni2': 212.8, 'Zn2': 231.6, 'Co2': 209.7,
                   'Cr2': 136.8, 'Fe2': 163.0, 'Mg2': 132.9, 'V2': 195.7,
                   'Mn2': 146.1, 'Hg2': 288.8, 'Cd2': 185.6, 'Ca2': 87.3,
                   'Sn2': 187.9, 'Sr2': 82.7, 'Ba2': 71.9, 'Al3': 399.0,
                   'Fe3': 428.0, 'Cr3': 258.0, 'In3': 347.0, 'Tl3': 456.0,
                   'Y3': 216.0, 'La3': 152.0, 'Ce3': 230.0, 'Pr3': 264.0,
                   'Nd3': 213.0, 'Sm3': 230.0, 'Eu3': 259.0, 'Gd3': 198.0,
                   'Tb3': 235.0, 'Dy3': 207.0, 'Er3': 251.0, 'Tm3': 282.0,
                   'Lu3': 249.0, 'Hf4': 827.0, 'Zr4': 761.0, 'Ce4': 706.0,
                   'U4': 1034.0, 'Pu4': 828.0, 'Th4': 512.0},
        'TIP4PEW' : {'Li1': 36.0, 'Na1': 9.0, 'K1': 24.0, 'Rb1': 13.0,
                     'Cs1': 16.0, 'Tl1': 65.0, 'Cu1': 21.0, 'Ag1': 94.0,
                     'F-1': -67.0, 'Cl-1': -66.0, 'Br-1': -68.0, 'I-1': -62.0,
                     'Be2': 228.5, 'Cu2': 339.2, 'Ni2': 259.2, 'Zn2': 272.3,
                     'Co2': 252.8, 'Cr2': 177.4, 'Fe2': 201.1, 'Mg2': 180.5,
                     'V2': 244.8, 'Mn2': 192.3, 'Hg2': 335.2, 'Cd2': 233.7,
                     'Ca2': 128.0, 'Sn2': 231.4, 'Sr2': 118.9, 'Ba2': 112.5,
                     'Al3': 488.0, 'Fe3': 519.0, 'Cr3': 322.0, 'In3': 425.0,
                     'Tl3': 535.0, 'Y3': 294.0, 'La3': 243.0, 'Ce3': 315.0,
                     'Pr3': 348.0, 'Nd3': 297.0, 'Sm3': 314.0, 'Eu3': 345.0,
                     'Gd3': 280.0, 'Tb3': 313.0, 'Dy3': 298.0, 'Er3': 328.0,
                     'Tm3': 356.0, 'Lu3': 331.0, 'Hf4': 956.0, 'Zr4': 895.0,
                     'Ce4': 835.0, 'U4': 1183.0, 'Pu4': 972.0, 'Th4': 625.0},
        'SPCE' :  {'Li1': 33.0, 'Na1': 6.0, 'K1': 19.0, 'Rb1': 7.0, 'Cs1': 12.0,
                   'Tl1': 61.0, 'Cu1': 9.0, 'Ag1': 92.0, 'F-1': -53.0,
                   'Cl-1': -55.0, 'Br-1': -51.0, 'I-1': -51.0, 'Be2': 188.1,
                   'Cu2': 304.4, 'Ni2': 205.2, 'Zn2': 231.2, 'Co2': 209.2,
                   'Cr2': 131.2, 'Fe2': 155.4, 'Mg2': 122.2, 'V2': 206.6,
                   'Mn2': 154.9, 'Hg2': 300.2, 'Cd2': 198.8, 'Ca2': 89.0,
                   'Sn2': 201.1, 'Sr2': 96.3, 'Ba2': 85.8, 'Al3': 406.0,
                   'Fe3': 442.0, 'Cr3': 254.0, 'In3': 349.0, 'Tl3': 455.0,
                   'Y3': 209.0, 'La3': 165.0, 'Ce3': 242.0, 'Pr3': 272.0,
                   'Nd3': 235.0, 'Sm3': 224.0, 'Eu3': 273.0, 'Gd3': 186.0,
                   'Tb3': 227.0, 'Dy3': 206.0, 'Er3': 247.0, 'Tm3': 262.0,
                   'Lu3': 247.0, 'Hf4': 810.0, 'Zr4': 760.0, 'Ce4': 694.0,
                   'U4': 1043.0, 'Pu4': 828.0, 'Th4': 513.0},
        'OPC3' :  {'Li1': 29.0, 'Na1': 2.0, 'K1': 16.0, 'Rb1': 8.0, 'Cs1': 6.0,
                   'Tl1': 63.0, 'Cu1': 12.0, 'Ag1': 83.0, 'F-1': -40.0,
                   'Cl-1': -47.0, 'Br-1': -43.0, 'I-1': -45.0,
                   'Be2': 186.0, 'Cu2': 269.0, 'Ni2': 207.0, 'Zn2': 199.0,
                   'Co2': 182.0, 'Cr2': 109.0, 'Fe2': 131.0, 'Mg2': 117.0,
                   'V2': 201.0, 'Mn2': 137.0, 'Hg2': 276.0, 'Cd2': 200.0,
                   'Ca2': 76.0, 'Sn2': 188.0, 'Sr2': 85.0, 'Ba2': 77.0,
                   'Al3': 363.0, 'Fe3': 429.0, 'Cr3': 209.0, 'In3': 330.0,
                   'Tl3': 437.0, 'Y3': 192.0, 'La3': 131.0, 'Ce3': 215.0,
                   'Pr3': 255.0, 'Nd3': 184.0, 'Sm3': 188.0, 'Eu3': 233.0,
                   'Gd3': 164.0, 'Tb3': 199.0, 'Dy3': 183.0, 'Er3': 228.0,
                   'Tm3': 246.0, 'Lu3': 222.0,
                   'Hf4': 718.0, 'Zr4': 707.0, 'Ce4': 653.0, 'U4': 980.0,
                   'Pu4': 817.0, 'Th4': 452.0},
        'OPC' :  {'Li1': 29.0, 'Na1': 0.0, 'K1': 20.0, 'Rb1': 6.0, 'Cs1': 13.0,
                   'Tl1': 60.0, 'Cu1': 16.0, 'Ag1': 83.0, 'F-1': -67.0,
                   'Cl-1': -69.0, 'Br-1': -60.0, 'I-1': -60.0,
                   'Be2': 214.0, 'Cu2': 291.0, 'Ni2': 212.0, 'Zn2': 225.0,
                   'Co2': 204.0, 'Cr2': 132.0, 'Fe2': 154.0, 'Mg2': 127.0,
                   'V2': 239.0, 'Mn2': 175.0, 'Hg2': 289.0, 'Cd2': 219.0,
                   'Ca2': 86.0, 'Sn2': 199.0, 'Sr2': 87.0, 'Ba2': 78.0,
                   'Al3': 399.0, 'Fe3': 531.0, 'Cr3': 243.0, 'In3': 413.0,
                   'Tl3': 479.0, 'Y3': 260.0, 'La3': 165.0, 'Ce3': 289.0,
                   'Pr3': 311.0, 'Nd3': 243.0, 'Sm3': 236.0, 'Eu3': 279.0,
                   'Gd3': 222.0, 'Tb3': 256.0, 'Dy3': 243.0, 'Er3': 298.0,
                   'Tm3': 314.0, 'Lu3': 289.0,
                   'Hf4': 847.0, 'Zr4': 804.0, 'Ce4': 789.0, 'U4': 1123.0,
                   'Pu4': 941.0, 'Th4': 598.0},
        'FB3' :  {'Li1': 30.0, 'Na1': 2.0, 'K1': 15.0, 'Rb1': 7.0, 'Cs1': 17.0,
                   'Tl1': 65.0, 'Cu1': 17.0, 'Ag1': 85.0, 'F-1': -45.0,
                   'Cl-1': -49.0, 'Br-1': -40.0, 'I-1': -52.0,
                   'Be2': 193.0, 'Cu2': 279.0, 'Ni2': 223.0, 'Zn2': 217.0,
                   'Co2': 192.0, 'Cr2': 138.0, 'Fe2': 157.0, 'Mg2': 128.0,
                   'V2': 212.0, 'Mn2': 149.0, 'Hg2': 289.0, 'Cd2': 201.0,
                   'Ca2': 92.0, 'Sn2': 205.0, 'Sr2': 91.0, 'Ba2': 78.0,
                   'Al3': 387.0, 'Fe3': 446.0, 'Cr3': 232.0, 'In3': 343.0,
                   'Tl3': 464.0, 'Y3': 218.0, 'La3': 155.0, 'Ce3': 251.0,
                   'Pr3': 291.0, 'Nd3': 211.0, 'Sm3': 218.0, 'Eu3': 261.0,
                   'Gd3': 191.0, 'Tb3': 226.0, 'Dy3': 219.0, 'Er3': 256.0,
                   'Tm3': 285.0, 'Lu3': 250.0,
                   'Hf4': 787.0, 'Zr4': 769.0, 'Ce4': 701.0, 'U4': 1044.0,
                   'Pu4': 817.0, 'Th4': 507.0},
        'FB4' :  {'Li1': 33.0, 'Na1': 8.0, 'K1': 25.0, 'Rb1': 9.0, 'Cs1': 13.0,
                   'Tl1': 68.0, 'Cu1': 25.0, 'Ag1': 90.0, 'F-1': -57.0,
                   'Cl-1': -55.0, 'Br-1': -51.0, 'I-1': -53.0,
                   'Be2': 227.0, 'Cu2': 313.0, 'Ni2': 218.0, 'Zn2': 239.0,
                   'Co2': 206.0, 'Cr2': 159.0, 'Fe2': 187.0, 'Mg2': 133.0,
                   'V2': 234.0, 'Mn2': 181.0, 'Hg2': 331.0, 'Cd2': 227.0,
                   'Ca2': 109.0, 'Sn2': 215.0, 'Sr2': 103.0, 'Ba2': 95.0,
                   'Al3': 427.0, 'Fe3': 502.0, 'Cr3': 286.0, 'In3': 403.0,
                   'Tl3': 514.0, 'Y3': 268.0, 'La3': 211.0, 'Ce3': 294.0,
                   'Pr3': 326.0, 'Nd3': 256.0, 'Sm3': 257.0, 'Eu3': 302.0,
                   'Gd3': 238.0, 'Tb3': 270.0, 'Dy3': 253.0, 'Er3': 304.0,
                   'Tm3': 320.0, 'Lu3': 303.0,
                   'Hf4': 837.0, 'Zr4': 845.0, 'Ce4': 771.0, 'U4': 1140.0,
                   'Pu4': 919.0, 'Th4': 601.0}
}

def params1264(parm, mask, c4file, watermodel, polfile, tunfactor):
   
    from .. import periodic_table as pt

    try:
        pollist = _get_params(polfile)
    except ValueError as err:
        raise LJ12_6_4Error(
            f'Bad polarizability file {polfile}. Expected a file with '
            '2 columns: <Amber Atom Type> <Polarizability>'
        ) from err
   
    if c4file is None:
        c4list = DEFAULT_C4_PARAMS[watermodel]
    else:
        try:
            c4list = _get_params(c4file)
        except ValueError as err:
            raise LJ12_6_4Error(
                f'Bad C4 parameter file {c4file}. Expected a file with 2 columns: '
                '<Atom Element> <C4 Parameter>'
            ) from err


    # Determine which atom type was treated as the center metal ion
    mettypdict = dict()
    for i in mask.Selected():
        mettypind = parm.parm_data['ATOM_TYPE_INDEX'][i]
        metchg = parm.parm_data['CHARGE'][i]
        if mettypind in mettypdict: continue
        mettypdict[mettypind] = (parm.atoms[i].atomic_number, int(metchg))
        print(f"The selected metal ion is {pt.Element[parm.atoms[i].atomic_number]}")
    mettypinds = sorted(mettypdict.keys())

    # 1. Get the dict between AMBER_ATOM_TYPE and ATOM_TYPE_INDEX for all 
    # the atoms in prmtop file

    ntypes = parm.pointers['NTYPES']
    typs = parm.parm_data['AMBER_ATOM_TYPE']
    typinds = parm.parm_data['ATOM_TYPE_INDEX']

    typdict = {}
    for ty, ind in zip(typs, typinds):
        if ind not in typdict:
            typdict[ind] = [ty]
        elif ty in typdict[ind]:
            continue
        else:
            typdict[ind].append(ty)

    for i in range(1, ntypes+1):
        if i not in typinds:
            typdict[i] = []

    # 2.Generate the C4 term for each atom type pair
    result = [0.0 for i in parm.parm_data['LENNARD_JONES_ACOEF']]

    for mettypind in mettypinds:
        # Obtain the C4 parameters
        c4 = c4list[pt.Element[mettypdict[mettypind][0]] + str(mettypdict[mettypind][1])]
        i = mettypind - 1
        for j in range(1, ntypes+1):
            jj = j - 1
            attypjs = typdict[j]
            if len(attypjs) >= 1:
                # Get polarizability
                for k, typjs in enumerate(attypjs):
                    if k == 0:
                        try:
                            pol = pollist[typjs]
                        except KeyError:
                            raise LJ12_6_4Error(f"Could not find parameters for ATOM_TYPE {typjs}")
                    else:
                        try:
                            anthpol = pollist[typjs]
                            if anthpol != pol:
                                raise LJ12_6_4Error(
                                    f'Polarizability parameter of AMBER_ATOM_TYPE {attypjs[0]} is '
                                    f'not the same as that of AMBER_ATOM_TYPE {typjs}, but '
                                    'their VDW parameters are the same. '
                                )
                        except KeyError:
                            raise LJ12_6_4Error(f"Could not find parameters for ATOM_TYPE {typjs}")
                # Get index
                if jj < i:
                    idx = parm.parm_data['NONBONDED_PARM_INDEX'][ntypes*jj+i]-1
                else:
                    idx = parm.parm_data['NONBONDED_PARM_INDEX'][ntypes*i+jj]-1

                # Caculate C4 terms
                if attypjs == ['OW']:
                    # There is only one C4 term exist between water and a
                    # certain ion
                    result[idx] = c4
                else:
                    # There are two C4 terms need to add together between two
                    # different ions
                    result[idx] += c4 / WATER_POL * pol * tunfactor
    return result

def _get_params(fname):
    params = dict()

    with open(fname, 'r') as f:
        for line in f:
            atomtype, param = line.split()[:2]
            param = float(param)
            if atomtype in params and abs(param - params[atomtype]) > 0.0001:
                warnings.warn(f'Atom type {atomtype} has multiple parameters in {fname}.',
                              DuplicateParamWarning)
            params[atomtype] = param

        return params
