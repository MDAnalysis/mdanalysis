# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
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

The original raw data are stored as multi-line strings that are
translated into dictionaries with :func:`kv2dict`. In the future,
these tables might be moved into external data files; see
:func:`kv2dict` for explanation of the file format.

.. autofunction:: kv2dict

The raw tables are stored in the strings

.. autodata:: TABLE_ATOMELEMENTS
.. autodata:: TABLE_MASSES
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

# other types are guessed from the name; see
# topology.core.guess_atom_elements()
"""

#: Dictionary with hard-coded special atom names, used for guessing atom types
#: with :func:`MDAnalysis.topology.core.guess_atom_type`.
atomelements = kv2dict(TABLE_ATOMELEMENTS)

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
Bh    262
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
"""

#: Dictionary table with atomic masses in u, indexed by the element from
#: :data:`atomelements`.
masses = kv2dict(TABLE_MASSES, convertor=float)

#: Van der Waals radii (taken from GROMACS_, ``/usr/share/gromacs/top/vdwradii.dat``)
#: and converted to ångström.
#: .. _GROMACS: http://www.gromacs.org
#:
#: .. SeeAlso:: :func:`MDAnalysis.topology.core.guess_bonds`
vdwradii = {
    "C":     1.5,
    "F":     1.2,
    "H":     0.4,
    "N":     1.10,
    "O":     1.05,
    "S":     1.6,
    "P":     1.6,
}
