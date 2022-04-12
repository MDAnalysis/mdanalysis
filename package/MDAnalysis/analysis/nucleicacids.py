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

import MDAnalysis as mda
from ..lib import mdamath
from ..analysis.base import AnalysisBase, Results

from ..exceptions import SelectionError


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
    """A class for developing analyses that run on a pair of base pairs """

    _sel: List[Tuple[str, str]]
    _unv: mda.Universe
    _s1: mda.AtomGroup
    _s2: mda.AtomGroup
    _res_dict: Dict[str, List[Union[float, str]]]
    _names: List[str]

    def __init__(self, universe: mda.Universe, selections: List[Tuple[str, str]], **kwargs) -> None:
        super(NucPairDist, self).__init__(universe.trajectory, **kwargs)
        self._unv = universe
        self._sel = selections
        self._check_selections()
        self.results = Results()

    def _check_selections(self) -> None:
        ag1: List[mda.AtomGroup] = []
        ag2: List[mda.AtomGroup] = []

        for s in self._sel:
            # Checking number of selections and building mda.AtomGroups
            if len(s) != 2:
                raise ValueError("Each base pair selection given as a tuple of length 2")
            try:
                ag1.append(self._unv.select_atoms(s[0]))

            except SelectionError:
                raise SelectionError(f"Invalid selection in {s[0]}")

            try:
                ag2.append(self._unv.select_atoms(s[1]))
            except SelectionError:
                raise SelectionError(f"Invalid selection in {s[1]}")

        self._s1 = mda.AtomGroup([ag[0] for ag in ag1])
        self._s2 = mda.AtomGroup([ag[0] for ag in ag2])

    def _prepare(self) -> None:
        self._names = [s[0] + s[1] for s in self._sel]
        self._res_dict = {k: [] for k in self._names}
        self._res_dict['time'] = []

    def _single_frame(self) -> None:
        dist = self._s1.positions - self._s2.positions
        wc: np.ndarray = mdamath.pnorm(dist)

        for i, n in enumerate(self._names):
            self._res_dict[n].append(wc[i])
            self._res_dict['time'].append(self._ts.time)

    def _conclude(self) -> None:
        self.results['times'] = np.array(self._res_dict['time'])
        self.results['selections'] = self._names
        for i, n in enumerate(self._names):
            self.results[i] = np.array(self._res_dict[n])


class WCDist(NucPairDist):
    """Watson-Crick basepair distance for selected residues"""

    def __init__(self, universe: mda.Universe, selections: List[Tuple[BaseSelect, BaseSelect]],
                 n1_name: str = 'N1', n3_name: str = "N3", **kwargs) -> None:
        sel_str: List[Tuple[str, str]] = []

        for s in selections:
            if universe.select_atoms(f" resid {s[0].resid} ").resnames[0] in ["DC", "DT", "U", "C", "T", "CYT", "THY", "URA"]:
                a1, a2 = n3_name, n1_name
            elif universe.select_atoms(f" resid {s[0].resid} ").resnames[0] in ["DG", "DA", "A", "G", "ADE", "GUA"]:
                a1, a2 = n1_name, n3_name
            else:
                raise SelectionError(f"Resid {s[0].resid} is not a Nucleic acid")

            sel_str.append((f'segid {s[0].segid} and resid {s[0].resid} and name {a1}',
                            f'segid {s[1].segid} and resid {s[1].resid} and name {a2}'))
        super(WCDist, self).__init__(universe, sel_str, **kwargs)


class MinorDist(NucPairDist):
    """"""

    def __init__(self, universe: mda.Universe, selections: List[Tuple[BaseSelect, BaseSelect]],
                 c2_name: str = "C2", o2_name: str = "O2", **kwargs) -> None:
        sel_str: List[Tuple[str, str]] = []

        for s in selections:
            if universe.select_atoms(f" resid {s[0].resid} ").resnames[0] in ["DC", "DT", "U", "C", "T", "CYT", "THY",
                                                                              "URA"]:
                a1, a2 = o2_name, c2_name
            elif universe.select_atoms(f" resid {s[0].resid} ").resnames[0] in ["DG", "DA", "A", "G", "ADE", "GUA"]:
                a1, a2 = c2_name, o2_name
            else:
                raise SelectionError(f"Resid {s[0].resid} is not a Nucleic acid")

            sel_str.append((f'segid {s[0].segid} and resid {s[0].resid} and name {a1}',
                            f'segid {s[1].segid} and resid {s[1].resid} and name {a2}'))
        super(MinorDist, self).__init__(universe, sel_str, **kwargs)


class MajorDist(NucPairDist):
    """"""

    def __init__(self, universe: mda.Universe, selections: List[Tuple[BaseSelect, BaseSelect]],
                 n4_name: str = "N4", o6_name: str = "O6", **kwargs) -> None:
        sel_str: List[Tuple[str, str]] = []

        for s in selections:
            if universe.select_atoms(f" resid {s[0].resid} ").resnames[0] in ["DC", "DT", "U", "C", "T", "CYT", "THY",
                                                                              "URA"]:
                a1, a2 = n4_name, o6_name
            elif universe.select_atoms(f" resid {s[0].resid} ").resnames[0] in ["DG", "DA", "A", "G", "ADE", "GUA"]:
                a1, a2 = o6_name, n4_name
            else:
                raise SelectionError(f"Resid {s[0].resid} is not a Nucleic acid")

            sel_str.append((f'segid {s[0].segid} and resid {s[0].resid} and name {a1}',
                            f'segid {s[1].segid} and resid {s[1].resid} and name {a2}'))
        super(MajorDist, self).__init__(universe, sel_str, **kwargs)
