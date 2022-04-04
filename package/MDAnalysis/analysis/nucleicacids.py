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
Nucleic acid analysis
=====================

:Author: Alia Lescoulie
:Year: 2022
:copyright: GNU Public Licence v3

"""

from typing import List, Tuple

import numpy as np
import pandas as pd
from math import pi, sin, cos, atan2, sqrt, pow

from MDAnalysis.lib import mdamath
from MDAnalysis.analysis.base import AnalysisBase
from package import MDAnalysis
from package.MDAnalysis.exceptions import SelectionError



class WCBase(AnalysisBase):
    """Watson-Crick basepair distance for selected residues"""

    _sel: List[Tuple[str, str]]
    _unv: MDAnalysis.Universe
    _s1: MDAnalysis.AtomGroup
    _s2: MDAnalysis.AtomGroup

    def __init__(self, universe: MDAnalysis.Universe, selections: List[Tuple[str, str]], N1_name: str ='N1', N3_name: str ="N3", **kwargs) -> None:
        super(WCBase, self).__init__(universe.trajectory, **kwargs)
        self._sel = selections

        self._s1 = MDAnalysis.AtomGroup()
        self._s2 = MDAnalysis.AtomGroup()

        for s in self._sel:
            # Checking number of selections and building atomgroups
            if (len(s) != 2):
                raise ValueError("Each base pair selection given as a tuple of lenght 2")
            try:
                self._s1 += self._unv.select_atoms(f'({s[0]} and name {N1_name})')

            except SelectionError:
                raise SelectionError(f"Invalid selection in {s[0]}")

            try: 
                self._s2 += self._unv.select_atoms(f'({s[0]} and name {N3_name})')
            except SelectionError:
                raise SelectionError(f"Invalid selection in {s[1]}")
        
    def _prepare(self) -> None:
        self._res_dict = {k: [] for k in range(len(self._sel))}
    
    def _single_frame(self):
        wc = mdamath.norm(np.subtract(self._s1.positions, self._s2.positions))

        for k in range(len(self._sel)):
            self._res_dict[k].append(wc[k])
    
    def _conclude(self) -> None:
        self.results: pd.DataFrame = pd.DataFrame.from_dict(self._res_dict)