# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
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

"""
MDAnalysis topology tables
==========================

The module contains static lookup tables for atom typing etc. The
tables are dictionaries that are indexed by the element.

.. autodata:: atomelements
.. autodata:: masses
.. autodata:: vdwradii

The original raw data are stored as multi-line strings that are
translated into dictionaries with :func:`kv2dict`. In the future,
these tables might be moved into external data files; see
:func:`kv2dict` for explanation of the file format.

.. autofunction:: kv2dict

The raw tables are stored in the strings

.. autodata:: TABLE_ATOMELEMENTS
.. autodata:: TABLE_MASSES
.. autodata:: TABLE_VDWRADII
"""


def kv2dict(s, convertor=str):
    """Primitive ad-hoc parser of a key-value record list.

    * The string *s* should contain each key-value pair on a separate
      line (separated by newline). The first white space after the key
      separates key and value.

    * Empty lines are allowed.

    * Comment lines (starting with #) are allowed.

    * Leading whitespace is ignored.

    The *convertor* is a function that converts its single argument to
    a valid Python type. The default is :func:`str` but other
    possibilities are :func:`int` (for integers) or :func:`float` for
    floating point numbers.
    """
    d = {}
    lines = s.splitlines()
    for line in lines:
        line = line.lstrip()
        values = line.split(None, 1)
        if len(values) == 0 or line.startswith("#"):
            continue
        d[values[0]] = convertor(values[1])
    return d


#: Table with hard-coded special atom names, used for guessing atom types
#: with :func:`MDAnalysis.topology.core.guess_atom_element`.
TABLE_ATOMELEMENTS = """
# translation of atomnames to types/element
# based on CHARMM and AMBER usage with a little bit of GROMOS (and PROPKA)
# NOTE: CL might be ambiguous and is interpreted as chloride!

# --------- ------------------
# atomname   element
# --------- ------------------

# Bromide
BR           BR

# Calcium
CAL          CA
C0           CA
CA2+         CA

# Cesium
CES          CS

# Chloride
CLA          CL
CLAL         CL
CL           CL
CL-          CL

# Iodide
IOD          I

# Iron
FE           FE
FE2          FE

# Lithium
LIT          LI
LI           LI
LI+          LI
QL           LI

# Magnesium
MG           MG
MG2+         MG

# Noble gases
## XXX collides with NE, HE in Arg  XXX
## XXX so we remove the noble gases XXX
##HE           HE
##NE           NE

# Potassium
K            K
POT          K
K+           K
QK           K

# Sodium
SOD          NA
NA           NA
NA+          NA
QN           NA

# Zink
ZN           ZN

# Copper
CU           CU

# Cesium
CS           CS
CS+          CS
CES          CS

# Cerium??
QC           CE

# Rubidium
RB           RB
QR           RB

# special carbons (Amber?)
BC           C
AC           C

# dummy atom types
MW           DUMMY

# other types are guessed from the name; see
# topology.core.guess_atom_elements()
"""

#: Dictionary with hard-coded special atom names, used for guessing atom types
#: with :func:`MDAnalysis.topology.core.guess_atom_type`.
atomelements = kv2dict(TABLE_ATOMELEMENTS)

elements = ['H',
            'LI', 'BE', 'B', 'C', 'N', 'O', 'F',
            'NA', 'MG', 'AL', 'P', 'SI', 'S', 'CL',
            'K']

#: Plain-text table with atomic masses in u.
TABLE_MASSES = """
# masses for elements in atomic units (u)
# (taken from CHARMM and Gromacs atommass.dat)

#------------ -----------
# atomtype    mass
#------------ -----------
Ac    227.028
Al    26.981539
Am    243
Sb    121.757
Ar    39.948
As    74.92159
At    210
Ba    137.327
Bk    247
Be    9.012182
Bi    208.98037
Bh    264
B     10.811
BR    79.90400
Cd    112.411
CA    40.08000
Cf    251
C     12.01100
Ce    140.11600
CS    132.90000
CL    35.45000
Cr    51.9961
Co    58.9332
CU    63.54600
Cm    247
Db    262
Dy    162.5
Es    252
Er    167.26
Eu    151.965
Fm    257
F     18.99800
Fr    223
Gd    157.25
Ga    69.723
Ge    72.61
Au    196.96654
Hf    178.49
Hs    265
HE    4.00260
Ho    164.93032
H     1.00800
In    114.82
I     126.90450
Ir    192.22
FE    55.84700
Kr    83.8
La    138.9055
Lr    262
Pb    207.2
Li    6.941
Lu    174.967
MG    24.30500
Mn    54.93805
Mt    266
Md    258
Hg    200.59
Mo    95.94
N     14.00700
NA    22.98977
Nd    144.24
NE    20.17970
Np    237.048
Ni    58.6934
Nb    92.90638
No    259
Os    190.2
O     15.99900
Pd    106.42
P     30.97400
Pt    195.08
Pu    244
Po    209
K     39.10200
Pr    140.90765
Pm    145
Pa    231.0359
Ra    226.025
Rn    222
Re    186.207
Rh    102.9055
RB    85.46780
Ru    101.07
Rf    261
Sm    150.36
Sc    44.95591
Sg    263
Se    78.96
Si    28.0855
Ag    107.8682
Na    22.989768
Sr    87.62
S     32.06000
Ta    180.9479
Tc    98
Te    127.6
Tb    158.92534
Tl    204.3833
Th    232.0381
Tm    168.93421
Sn    118.71
Ti    47.88
W     183.85
U     238.0289
V     50.9415
Xe    131.29
Yb    173.04
Y     88.90585
ZN    65.37000
Zr    91.224
DUMMY 0.0
"""

#: Dictionary table with atomic masses in u, indexed by the element from
#: :data:`atomelements`.
masses = kv2dict(TABLE_MASSES, convertor=float)

#: Plain-text table with vdw radii.
TABLE_VDWRADII = r"""
# Van der Waals radii taken from
# [1] Bondi, A. (1964). "Van der Waals Volumes and Radii".
#     J. Phys. Chem. 68 (3): 441-451. doi:10.1021/j100785a001.
# [2] Rowland and Taylor (1996). "Intermolecular Nonbonded Contact Distances in Organic Crystal Structures:
#                                 Comparison with Distances Expected from van der Waals Radii".
#     J. Phys. Chem., 1996, 100 (18), 7384.7391. doi:10.1021/jp953141+.
# [3] Mantina, et al. (2009). "Consistent van der Waals Radii for the Whole Main Group".
#     J. Phys. Chem. A, 2009, 113 (19), 5806-5812. doi:10.1021/jp8111556.
#------------ -----------
# atomtype    r_vdw
#------------ -----------
H    1.10
HE   1.40
LI   1.82
BE   1.53
B    1.92
C    1.70
N    1.55
O    1.52
F    1.47
NE   1.54
NA   2.27
MG   1.73
AL   1.84
SI   2.10
P    1.80
S    1.80
CL   1.75
AR   1.88
K    2.75
CA   2.31
NI   1.63
CU   1.40
ZN   1.39
GA   1.87
GE   2.11
AA   1.85
SE   1.90
BR   1.85
KR   2.02
RR   3.03
SR   2.49
PD   1.63
AG   1.72
CD   1.58
IN   1.93
SN   2.17
SB   2.06
TE   2.06
I    1.98
XE   2.16
CS   3.43
BA   2.68
PT   1.75
AU   1.66
HH   1.55
TL   1.96
PB   2.02
BI   2.07
PO   1.97
AT   2.02
RN   2.20
FR   3.48
RA   2.83
U    1.86
"""

#: Dictionary table with vdw radii, indexed by the element from
#: :data:`atomelements`.
#: .. SeeAlso:: :func:`MDAnalysis.topology.core.guess_bonds`
vdwradii = kv2dict(TABLE_VDWRADII, convertor=float)

Z2SYMB = {1: 'H',                                                                2: 'He', 
          3: 'Li',   4: 'Be',  5: 'B',   6: 'C',   7: 'N',   8: 'O',   9: 'F',  10: 'Ne', 
          11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P',  16: 'S',  17: 'Cl', 18: 'Ar', 
          19: 'K',  20: 'Ca', 21: 'Sc', 22: 'Ti', 23: 'V',  24: 'Cr', 25: 'Mn', 26: 'Fe', 
          27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 
          35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y',  40: 'Zr', 41: 'Nb', 42: 'Mo', 
          43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn', 
          51: 'Sb', 52: 'Te', 53: 'I',  54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 
          59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 
          67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 
          75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg', 81: 'Tl', 82: 'Pb', 
          83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th', 
          91: 'Pa', 92: 'U',  93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk', 98: 'Cf', 
          99: 'Es', 100: 'Fm', 101: 'Md', 102: 'No', 103: 'Lr', 104: 'Rf', 105: 'Db', 
          106: 'Sg', 107: 'Bh', 108: 'Hs', 109: 'Mt', 110: 'Ds', 111: 'Rg', 112: 'Cn', 
          113: 'Nh', 114: 'Fl', 115: 'Mc', 116: 'Lv', 117: 'Ts', 118: 'Og'}

SYMB2Z = {v:k for k, v in Z2SYMB.items()}

# Conversion between SYBYL atom types and corresponding elements
# Tripos MOL2 file format:
#   https://web.archive.org/web/*/http://chemyang.ccnu.edu.cn/ccb/server/AIMMS/mol2.pdf
SYBYL2SYMB = {
    "H": "H", "H.spc": "H", "H.t3p": "H",
    "C.3": "C", "C.2": "C", "C.1": "C", "C.ar": "C", "C.cat": "C",
    "N.3": "N", "N.2": "N", "N.1": "N", "N.ar": "N",
    "N.am": "N", "N.pl3": "N", "N.4": "N",
    "O.3": "O", "O.2": "O", "O.co2": "O", "O.spc": "O", "O.t3p": "O",
    "S.3": "S", "S.2": "S", "S.O": "S", "S.O2": "S",
    "S.o": "S", "S.o2": "S",  # Non-standard but often found in the wild...
    "P.3": "P",
    "F": "F",
    "Li": "Li",
    "Na": "Na",
    "Mg": "Mg",
    "Al": "Al",
    "Si": "Si",
    "K": "K",
    "Ca": "Ca",
    "Cr.th": "Cr",
    "Cr.oh": "Cr",
    "Mn": "Mn",
    "Fe": "Fe",
    "Co.oh": "Co",
    "Cu": "Cu",
    "Cl": "Cl",
    "Br": "Br",
    "I": "I",
    "Zn": "Zn",
    "Se": "Se",
    "Mo": "Mo",
    "Sn": "Sn",
}