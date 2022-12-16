#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This tests reading a LAMMPS data file with a
"PairIJ Coeffs" section, issue #3336. 

"""

import MDAnalysis as mda


PAIRIJ_COEFFS_DATA = "PR3959_test.data"

u = mda.Universe(PAIRIJ_COEFFS_DATA)

print(str(u))
