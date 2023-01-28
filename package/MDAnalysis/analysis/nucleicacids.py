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

r"""
Updated nucleic acid analysis --- :mod:`MDAnalysis.analysis.nucleicacids`
=========================================================================

:Author: Alia Lescoulie
:Year: 2022
:copyright: GNU Public Licence v3

The module provides classes for analyzing nucleic acids structures.
This is an updated, higher performance version of previous nucleic acid tools.
For applications see :cite:p:`b-Denning2011,b-Denning2012`.

.. rubric:: References

.. bibliography::
    :filter: False
    :style: MDA
    :keyprefix: b-
    :labelprefix: ᵇ

    Denning2011
    Denning2012

Distances
---------

.. autoclass:: NucPairDist
    :members:
    :inherited-members:

.. autoclass:: WatsonCrickDist
    :members:
    :inherited-members:

.. versionadded 2.2.0

"""

from typing import List, Dict
import warnings

import numpy as np

import MDAnalysis as mda
from .distances import calc_bonds
from .base import AnalysisBase, Results
from MDAnalysis.core.groups import Residue


class NucPairDist(AnalysisBase):
    r"""Atom Pair distance calculation base class.

    Takes two lists of :class:`~MDAnalysis.core.groups.AtomGroup` and
    computes the distances between them over a trajectory. Used as a
    superclass for the other nucleic acid distances classes. The distance
    will be measured between atoms sharing an index in the two lists of
    :class:`~MDAnalysis.core.groups.AtomGroup`.

    Parameters
    ----------
    selection1: List[AtomGroup]
        List of :class:`~MDAnalysis.core.groups.AtomGroup` containing an atom
        of each nucleic acid being analyzed.
    selection2: List[AtomGroup]
        List of :class:`~MDAnalysis.core.groups.AtomGroup` containing an atom
        of each nucleic acid being analyzed.
    kwargs: dict
        Arguments for :class:`~MDAnalysis.analysis.base.AnalysisBase`

    Attributes
    ----------
    times: numpy.ndarray
        Simulation times for analysis.
    results.pair_distances: numpy.ndarray
        2D array of pair distances. First dimension is simulation time, second
        dimension contains the pair distances for each each entry pair in
        selection1 and selection2.

        .. versionadded:: 2.4.0


    Raises
    ------

    ValueError
        If the selections given are not the same length


    .. versionchanged:: 2.5.0
       The ability to access by passing selection indices to :attr:`results` is
       is now removed as of MDAnalysis version 2.5.0. Please use
       :attr:`results.pair_distances` instead.
       The :attr:`results.times` was deprecated and is now removed as of
       MDAnalysis 2.5.0. Please use the class attribute :attr:`times` instead.
    """

    _s1: mda.AtomGroup
    _s2: mda.AtomGroup
    _n_sel: int
    _res_dict: Dict[int, List[float]]
    _times = List[float]

    def __init__(self, selection1: List[mda.AtomGroup],
                 selection2: List[mda.AtomGroup],
                 **kwargs) -> None:
        super(NucPairDist, self).__init__(selection1[0].universe.trajectory, **kwargs)

        if len(selection1) != len(selection2):
            raise ValueError("Selections must be same length")

        self._n_sel: int = len(selection1)

        self._s1 = selection1[0]
        self._s2 = selection2[0]

        for i in range(1, self._n_sel):
            self._s1 += selection1[i]
            self._s2 += selection2[i]

    def _prepare(self) -> None:
        self._res_array: np.ndarray = np.zeros([self.n_frames, self._n_sel])

    def _single_frame(self) -> None:
        dist: np.ndarray = calc_bonds(self._s1.positions, self._s2.positions)
        self._res_array[self._frame_index, :] = dist

    def _conclude(self) -> None:
        self.results['pair_distances'] = self._res_array


class WatsonCrickDist(NucPairDist):
    r"""Watson-Crick basepair distance for selected residues over a trajectory.

    Takes two lists of :class:`~MDAnalysis.core.groups.Residue` objects and calculates
    the Watson-Crick distance between them over the trajectory. Bases are matched by
    their index in the lists given as arguments.

    Parameters
    ----------
    strand1: List[Residue]
        First list of bases
    strand2: List[Residue]
        Second list of bases
    n1_name: str (optional)
        Name of Nitrogen 1 of nucleic acids, by default assigned to N1
    n3_name: str (optional)
        Name of Nitrogen 3 of nucleic acids, by default assigned to N3
    g_name: str (optional)
        Name of Guanine in topology, by default assigned to G
    a_name: str (optional)
        Name of Adenine in topology, by default assigned to A
    u_name: str (optional)
        Name of Uracil in topology, by default assigned to U
    t_name: str (optional)
        Name of Thymine in topology, by default assigned to T
    c_name: str (optional)
        Name of Cytosine in topology, by default assigned to C
    **kwargs: dict
        arguments for :class:`~MDAnalysis.analysis.base.AnalysisBase`

    Attributes
    ----------
    times: numpy.ndarray
        Simulation times for analysis.
    results.pair_distances: numpy.ndarray
        2D array of Watson-Crick basepair distances. First dimension is
        simulation time, second dimension contains the pair distances for
        each each entry pair in strand1 and strand2.

        .. versionadded:: 2.4.0


    Raises
    ------
    ValueError
        If the residues given are not amino acids
    ValueError
        If the selections given are not the same length


    .. versionchanged:: 2.5.0
       Accessing results by passing strand indices to :attr:`results` is
       was deprecated and is now removed as of MDAnalysis version 2.5.0. Please
       use :attr:`results.pair_distances` instead.
       The :attr:`results.times` was deprecated and is now removed as of
       MDAnalysis 2.5.0. Please use the class attribute :attr:`times` instead.
    """

    def __init__(self, strand1: List[Residue], strand2: List[Residue],
                 n1_name: str = 'N1', n3_name: str = "N3",
                 g_name: str = 'G', a_name: str = 'A', u_name: str = 'U',
                 t_name: str = 'T', c_name: str = 'C',
                 **kwargs) -> None:
        sel1: List[mda.AtomGroup] = []
        sel2: List[mda.AtomGroup] = []
        strand = zip(strand1, strand2)

        for s in strand:
            if s[0].resname[0] in [c_name, t_name, u_name]:
                a1, a2 = n3_name, n1_name
            elif s[0].resname[0] in [a_name, g_name]:
                a1, a2 = n1_name, n3_name
            else:
                raise ValueError(f"{s} are not valid nucleic acids")

            sel1.append(s[0].atoms.select_atoms(f'name {a1}'))
            sel2.append(s[1].atoms.select_atoms(f'name {a2}'))

        super(WatsonCrickDist, self).__init__(sel1, sel2, **kwargs)
