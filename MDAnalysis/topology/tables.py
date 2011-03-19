# MDAnalysis
# 
"""
MDAnalysis topology tables
==========================

The module contains static lookup tables for atom typing etc. 

.. autodata:: atomtypes
.. autodata:: masses

"""

def kv2dict(s, convertor=str):
    """Primitive ad-hoc parser of a key-value record list.

    The string s should contain each key-value pair on a separate
    line. The first white space after the key separates key and value.

    Empty lines are allowed.

    Comment lines (starting with #) are allowed.
    """
    d = {}
    lines = s.splitlines()
    for line in lines:
        line = line.lstrip()
        values = line.split(None,1)
        if len(values) == 0 or line.startswith("#"):
            continue
        d[values[0]] = convertor(values[1])
    return d

#: Table with hard-coded special atom names, used for guessing atom types
#: with :func:`MDAnalysis.topology.core.guess_atom_type`.
TABLE_ATOMTYPES = """
# translation of atomnames to types/element
# based on CHARMM and AMBER usage with a little bit of GROMOS
# NOTE: CL might be ambiguous and is interpreted as chloride!

# --------- ------------------
# atomname   element
# --------- ------------------

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

# Iron
FE           FE

# Lithium
LIT          LI
LI           LI
LI+          LI
QL           LI

# Magnesium
MG           MG
MG2+         MG

# Noble gases
HE           HE
NE           NE

# Potassium
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
# topology.core.guess_atom_type()
"""

#: Dictionary with hard-coded special atom names, used for guessing atom types
#: with :func:`MDAnalysis.topology.core.guess_atom_type`.
atomtypes = kv2dict(TABLE_ATOMTYPES)

#: Plain-text table with atomic masses in u.
TABLE_MASSES = """
# masses for elements/atomtypes in atomic units (u)
# (taken from CHARMM and Gromacs atommass.dat)

#------------ -----------
# atomtype    mass
#------------ -----------
H               1.00800
C              12.01100
N              14.00700
O              15.99900
S              32.06000
P              30.97400
HE              4.00260
NE             20.17970
FE             55.84700
NA             22.98977
K              39.10200
RB             85.46780  
CS            132.90000
MG             24.30500
CA             40.08000
CU             63.54600
ZN             65.37000
Ce            140.11600
F              18.99800
CL             35.45000
BR             79.90400
I             126.90450    
"""

#: Dictionary table with atomic masses in u, indexed by the atom types from
#: :data:`atomtypes`.
masses = kv2dict(TABLE_MASSES, convertor=float)
