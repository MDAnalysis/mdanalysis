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

from typing import List, Tuple, Dict, Union, NamedTuple

import numpy as np
import pandas as pd
from math import pi, sin, cos, atan2, sqrt, pow

import MDAnalysis
from MDAnalysis.lib import mdamath
from MDAnalysis.analysis.base import AnalysisBase
from package import MDAnalysis
from package.MDAnalysis.exceptions import SelectionError


class BaseSelect(NamedTuple):
    """A named tuple for storing the selection information of a specific base pair of a nucleic acid.
    This makes passing in selections more organized.

    :Arguments:
    *segid*
        The name of the segment the base is on

    *base*
        The number of the base pair in the segment

    .. rubric:: Example

    Usage::

        U115Sel = BaseSelect('RNAA', 1)

    """
    segid: str
    resid: int


class NucPairDist(AnalysisBase):
    """An class for developing analyses that run on a pair of base pairs """

    _sel: List[Tuple[str, str]]
    _unv: MDAnalysis.Universe
    _s1: MDAnalysis.AtomGroup
    _s2: MDAnalysis.AtomGroup
    _res_dict: Dict[str, List[Union[float, str]]]
    _names: List[str]
    results: pd.DataFrame

    def __init__(self, universe: MDAnalysis.Universe, selections: List[Tuple[str, str]], **kwargs) -> None:
        super(NucPairDist, self).__init__(universe.trajectory, **kwargs)
        self._unv = universe
        self._sel = selections
        self._s1 = MDAnalysis.AtomGroup()
        self._s2 = MDAnalysis.AtomGroup()

    def _check_selections(self) -> None:
        for s in self._sel:
            # Checking number of selections and building AtomGroups
            if len(s) != 2:
                raise ValueError("Each base pair selection given as a tuple of lenght 2")
            try:
                self._s1 += self._unv.select_atoms(s[0])

            except SelectionError:
                raise SelectionError(f"Invalid selection in {s[0]}")

            try:
                self._s2 += self._unv.select_atoms(s[1])
            except SelectionError:
                raise SelectionError(f"Invalid selection in {s[1]}")

            if len(self._s1) != len(self._sel) or len(self._s2) != len(self._sel):
                raise ValueError("must only select 1 atom per selection")

    def _prepare(self) -> None:
        self._col = ('selection', 'time', 'distance')
        self._res_dict = {k: [] for k in range(len(self._col))}
        self.results: pd.DataFrame = pd.DataFrame(columns=self._col)
        self._names = [s[0] + s[1] for s in self._sel]

    def _single_frame(self) -> None:
        wc: np.ndarray = mdamath.pnorm(np.subtract(self._s1.positions, self._s2.positions))

        for i, n in enumerate(self._names):
            self._res_dict['selection'].append(n)
            self._res_dict['time'].append(self._ts.time)
            self._res_dict['distance'].append(wc[i])

    def _conclude(self) -> None:
        for k in self._col:
            self.results[k] = self._res_dict[k]


class WCDist(NucPairDist):
    """Watson-Crick basepair distance for selected residues"""

    def __init__(self, universe: MDAnalysis.Universe, selections: List[Tuple[BaseSelect, BaseSelect]],
                 n1_name: str = 'N1', n3_name: str = "N3", **kwargs) -> None:
        sel_str = [(f'segid {s[0].segid} and resid {s[0].resid} and {n1_name}',
                    f'segid {s[1].segid} and resid {s[1].resid} and {n3_name}') for s in selections]
        super(WCDist, self).__init__(universe.trajectory, sel_str, **kwargs)
