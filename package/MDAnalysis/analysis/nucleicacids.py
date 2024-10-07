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
:Year: 2022-2023
:copyright: LGPLv2.1

The module provides classes for analyzing nucleic acids structures.
This is an updated, higher performance version of previous nucleic acid tools.
For applications see :footcite:p:`Denning2011,Denning2012`.

.. rubric:: References

.. footbibliography::

Distances
---------

.. autoclass:: NucPairDist
    :members:
    :inherited-members:

.. autoclass:: WatsonCrickDist
    :members:
    :exclude-members: select_strand_atoms
    :inherited-members:

.. autoclass:: MinorPairDist
    :members:
    :exclude-members: select_strand_atoms
    :inherited-members:

.. autoclass:: MajorPairDist
    :members:
    :exclude-members: select_strand_atoms
    :inherited-members:

.. versionadded 2.2.0

"""

from typing import List, Dict, Tuple, Union
import warnings

import numpy as np

import MDAnalysis as mda
from .distances import calc_bonds
from .base import AnalysisBase, Results
from MDAnalysis.core.groups import Residue, ResidueGroup


# Deprecation: In 3.0.0 change type to just
# ResidueClass = ResidueGroup
ResidueClass = Union[List[Residue], ResidueGroup]
r"""A type alias for :code:`Union[List[Residue], ResidueGroup]`

Used as an alias for methods where either class is acceptable.
"""


class NucPairDist(AnalysisBase):
    r"""Atom pair distance calculation base class.

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
    results.pair_distances: numpy.ndarray
        2D array of pair distances. First dimension is simulation time,
        second dimension contains the pair distances for each each entry
        pair in selection1 and selection2.

        .. versionadded:: 2.4.0

        .. note::
            `results.pair_distances` is slated for deprecation in
            version 3.0.0, use `results.distances` instead.
        .. deprecated:: 2.7.0
            `results.pair_distances` will be removed in
            version 3.0.0, use :attr:`results.distances` instead.

    results.distances: numpy.ndarray
        stored in a 2d numpy array with first index selecting the
        Residue pair, and the second index selecting the frame number
        Distances are stored in a 2d numpy array with axis 0 (first index)
        indexing the trajectory frame and axis 1 (second index) selecting the
        Residue pair.

        .. versionadded:: 2.7.0

    times: numpy.ndarray
        Simulation times for analysis.


    Raises
    ------
    ValueError
        If the selections given are not the same length
    ValueError
        An :class:`~MDAnalysis.core.groups.AtomGroup` in one of the
        strands not a valid nucleic acid
    ValueError
        If a given residue pair from the provided strands returns an empty
        :class:`~MDAnalysis.core.groups.AtomGroup` when selecting the atom
        pairs used in the distance calculations


    *Version Info*

    .. versionchanged:: 2.5.0
       The ability to access by passing selection indices to :attr:`results`
       is now removed as of MDAnalysis version 2.5.0. Please use
       :attr:`results.pair_distances` instead.
       The :attr:`results.times` was deprecated and is now removed as of
       MDAnalysis 2.5.0.
       Please use the class attribute :attr:`times` instead.

    .. versionchanged:: 2.7.0
        Added static method :attr:`select_strand_atoms` as a
        helper for selecting atom pairs for distance analysis.
    """

    _analysis_algorithm_is_parallelizable = True

    @classmethod
    def get_supported_backends(cls):
        return ('serial', 'multiprocessing', 'dask',)
    
    _s1: mda.AtomGroup
    _s2: mda.AtomGroup
    _n_sel: int
    _res_dict: Dict[int, List[float]]

    def __init__(self, selection1: List[mda.AtomGroup],
                 selection2: List[mda.AtomGroup],
                 **kwargs) -> None:
        super(
            NucPairDist,
            self).__init__(
            selection1[0].universe.trajectory,
            **kwargs)

        if len(selection1) != len(selection2):
            raise ValueError("Selections must be same length")

        self._n_sel: int = len(selection1)

        self._s1 = selection1[0]
        self._s2 = selection2[0]

        for i in range(1, self._n_sel):
            self._s1 += selection1[i]
            self._s2 += selection2[i]

    @staticmethod
    def select_strand_atoms(
        strand1: ResidueGroup, strand2: ResidueGroup,
        a1_name: str, a2_name: str, g_name: str = 'G',
        a_name: str = 'A', u_name: str = 'U',
        t_name: str = 'T', c_name: str = 'C'
    ) -> Tuple[List[mda.AtomGroup], List[mda.AtomGroup]]:
        r"""
        A helper method for nucleic acid pair distance analyses.
        Used for selecting specific atoms from two strands of nucleic acids.


        Parameters
        ----------
        strand1: List[Residue]
            The first nucleic acid strand
        strand2: List[Residue]
            The second nucleic acid strand
        a1_name: str
            The selection for the purine base of the strand pair
        a2_name: str
            the selection for the pyrimidine base of the strand pair
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

        Returns
        -------
        Tuple[List[AtomGroup], List[AtomGroup]]
            returns a tuple containing two lists of
            :class:`~MDAnalysis.core.groups.AtomGroup`\s
            corresponding to the provided selections from each strand.

        Raises
        ------
        ValueError:
            An :class:`~MDAnalysis.core.groups.AtomGroup`
            in one of the strands not a valid nucleic acid
        ValueError:
            An :class:`~MDAnalysis.core.groups.Residue` returns an empty
            :class:`~MDAnalysis.core.groups.AtomGroup`
            with the provided selection


        .. versionadded:: 2.7.0
        """
        pyrimidines: List[str] = [c_name, t_name, u_name]
        purines: List[str] = [a_name, g_name]

        sel1: List[mda.AtomGroup] = []
        sel2: List[mda.AtomGroup] = []

        for pair in zip(strand1.residues, strand2.residues):
            if pair[0].resname[0] in pyrimidines:
                a1, a2 = a2_name, a1_name
            elif pair[0].resname[0] in purines:
                a1, a2 = a1_name, a2_name
            else:
                raise ValueError(
                    f"AtomGroup in {pair} is not a valid nucleic acid"
                )

            ag1 = pair[0].atoms.select_atoms(f'name {a1}')
            ag2 = pair[1].atoms.select_atoms(f'name {a2}')

            if not all(len(ag) > 0 for ag in [ag1, ag2]):
                err_info: Tuple[Residue, str] = (pair[0], a1) \
                    if len(ag1) == 0 else (pair[1], a2)

                raise ValueError(
                    (
                        f"{err_info[0]} returns an empty AtomGroup"
                        "with selection string \"name {a2}\""
                    )
                )

            sel1.append(ag1)
            sel2.append(ag2)

        return (sel1, sel2)

    def _prepare(self) -> None:
        self._res_array: np.ndarray = np.zeros(
            [self.n_frames, self._n_sel]
        )

    def _single_frame(self) -> None:
        dist: np.ndarray = calc_bonds(
            self._s1.positions, self._s2.positions
        )

        self._res_array[self._frame_index, :] = dist

    def _conclude(self) -> None:
        self.results['distances'] = self._res_array
        self.results['pair_distances'] = self.results['distances']
        # TODO: remove pair_distances in 3.0.0

    def _get_aggregator(self):
        return ResultsGroup(lookup={
        'distances': ResultsGroup.ndarray_vstack,
        'pair_distances': ResultsGroup.ndarray_vstack,}
        )

class WatsonCrickDist(NucPairDist):
    r"""
    Watson-Crick base pair distance for selected
    residues over a trajectory.

    Takes two :class:`~MDAnalysis.core.groups.ResidueGroup`
    objects or two lists of :class:`~MDAnalysis.core.groups.Residue`
    and calculates the distance between the nitrogen atoms in the
    Watson-Crick hydrogen bond over the trajectory. Bases are matched
    either by their index in the two
    :class:`~MDAnalysis.core.groups.ResidueGroup` provided as arguments,
    or based on the indices of the provided lists of
    :class:`~MDAnalysis.core.groups.Residue` objects depending
    on which is provided.

    .. note::
        Support for :class:`~MDAnalysis.core.groups.Residue` is slated for
        deprecation and will raise a warning when used. It still works but
        :class:`~MDAnalysis.core.groups.ResidueGroup` is preferred.

    Parameters
    ----------
    strand1: ResidueClass
        First list of bases

        .. deprecated:: 2.7.0
           Using a list of :class:`~MDAnalysis.core.groups.Residue` will
           be removed in 3.0.0. Use a
           :class:`~MDAnalysis.core.groups.ResidueGroup`.

    strand2: ResidueClass
        Second list of bases

        .. deprecated:: 2.7.0
           Using a list of :class:`~MDAnalysis.core.groups.Residue` will
           be removed in 3.0.0. Use a
           :class:`~MDAnalysis.core.groups.ResidueGroup`.

    n1_name: str (optional)
        Name of Nitrogen 1 of nucleic acids, by default assigned to "N1"
    n3_name: str (optional)
        Name of Nitrogen 3 of nucleic acids, by default assigned to "N3"
    g_name: str (optional)
        Name of Guanine in topology, by default assigned to "G"
    a_name: str (optional)
        Name of Adenine in topology, by default assigned to "A"
    u_name: str (optional)
        Name of Uracil in topology, by default assigned to "U"
    t_name: str (optional)
        Name of Thymine in topology, by default assigned to "T"
    c_name: str (optional)
        Name of Cytosine in topology, by default assigned to C
    **kwargs: dict
        Key word arguments for
        :class:`~MDAnalysis.analysis.base.AnalysisBase`

    Attributes
    ----------
    results.distances: numpy.ndarray
        Distances are stored in a 2d numpy array with axis 0 (first index)
        indexing the trajectory frame and axis 1 (second index) selecting the
        Residue pair.

        .. versionadded:: 2.7.0

    results.pair_distances: numpy.ndarray
        2D array of pair distances. First dimension is
        simulation time, second dimension contains the
        pair distances for each each entry pair in
        selection1 and selection2.

        .. versionadded:: 2.4.0

        .. deprecated:: 2.7.0
            `results.pair_distances` will be removed in version 3.0.0,
            use :attr:`results.distances` instead.

    times: numpy.ndarray
        Simulation times for analysis.

    Raises
    ------
    TypeError
        If the provided list of :class:`~MDAnalysis.core.Residue` contains
        non-Residue elements

        .. deprecated:: 2.7.0
           Starting with version 3.0.0, this exception will no longer
           be raised because only
           :class:`~MDAnalysis.core.groups.ResidueGroup` will be allowed.

    ValueError
        If `strand1` and `strand2` are not the same length
    ValueError:
        An :class:`~MDAnalysis.core.groups.AtomGroup`
        in one of the strands not a valid nucleic acid
    ValueError
        If a given residue pair from the provided strands returns an empty
        :class:`~MDAnalysis.core.groups.AtomGroup` when selecting the atom
        pairs used in the distance calculations


    *Version Info*

    .. versionchanged:: 2.5.0
       Accessing results by passing strand indices to :attr:`results`
       was deprecated and is now removed as of MDAnalysis version 2.5.0.
       Please use :attr:`results.pair_distances` instead.
       The :attr:`results.times` was deprecated and is now removed as of
       MDAnalysis 2.5.0. Please use the class attribute
       :attr:`times` instead.

    .. versionchanged:: 2.7.0
        `strand1` and `strand2` now also accept a
        :class:`~MDAnalysis.core.groups.ResidueGroup` as input.
        The previous input type, ``List[Residue]`` is still supported,
        but it is **deprecated** and will be removed in release 3.0.0.
    """

    def __init__(self, strand1: ResidueClass, strand2: ResidueClass,
                 n1_name: str = 'N1', n3_name: str = "N3",
                 g_name: str = 'G', a_name: str = 'A', u_name: str = 'U',
                 t_name: str = 'T', c_name: str = 'C',
                 **kwargs) -> None:

        def verify_strand(strand: ResidueClass) -> ResidueGroup:
            # Helper method to verify the strands

            if isinstance(strand, list):  # Checking if a list is given
                # verify list is only Residues
                if not all(isinstance(resid, Residue) for resid in strand):
                    raise TypeError(f"{strand} contains non-Residue elements")

                warnings.warn(
                    DeprecationWarning(
                        (
                            f"ResidueGroup should be used for {strand} instead"
                            "of giving a Residue list"
                        )
                    )
                )
                # Convert to a ResidueGroup
                strand: ResidueGroup = ResidueGroup(strand)

            return strand

        strand1: ResidueGroup = verify_strand(strand1)
        strand2: ResidueGroup = verify_strand(strand2)

        strand_atomgroups: Tuple[List[mda.AtomGroup], List[mda.AtomGroup]] = \
            self.select_strand_atoms(
                strand1, strand2, n1_name, n3_name,
                g_name=g_name, a_name=a_name,
                t_name=t_name, u_name=u_name, c_name=c_name
        )

        super(WatsonCrickDist, self).__init__(
            strand_atomgroups[0], strand_atomgroups[1], **kwargs
        )


class MinorPairDist(NucPairDist):
    r"""Minor-Pair basepair distance for selected residues over a trajectory.

    Takes two :class:`~MDAnalysis.core.groups.ResidueGroup` objects and
    calculates the Minor-groove hydrogen bond length between the
    nitrogen and oxygen atoms over the trajectory. Bases are
    matched by their index in the two
    :class:`~MDAnalysis.core.groups.ResidueGroup` provided as arguments.

    Parameters
    ----------
    strand1: List[Residue]
        First list of bases
    strand2: List[Residue]
        Second list of bases
    o2_name: str (optional)
        Name of Oxygen 2 of nucleic acids;
        by default assigned to "O2";
    c2_name: str (optional)
        Name of Carbon 2 of nucleic acids;
        by default assigned to "C2";
    g_name: str (optional)
        Name of Guanine in topology;
        by default assigned to "G";
    a_name: str (optional)
        Name of Adenine in topology
        by default assigned to "A";
    u_name: str (optional)
        Name of Uracil in topology;
        by default assigned to "U";
    t_name: str (optional)
        Name of Thymine in topology;
        by default assigned to "T";
    c_name: str (optional)
        Name of Cytosine in topology;
        by default assigned to "C";
    **kwargs:
        keyword arguments for
        :class:`~MDAnalysis.analysis.base.AnalysisBase`

    Attributes
    ----------
    results.distances: numpy.ndarray
        stored in a 2d numpy array with first index selecting
        the Residue pair, and the second index selecting the frame number
    times: numpy.ndarray
        Simulation times for analysis.

    Raises
    ------
    ValueError
        If the selections given are not the same length
        A :class:`~MDAnalysis.core.Residue` in
        one of the strands not a valid nucleic acid
    ValueError
        If a given residue pair from the provided strands returns an empty
        :class:`~MDAnalysis.core.groups.AtomGroup` when selecting the atom
        pairs used in the distance calculations


    .. versionadded:: 2.7.0
    """

    def __init__(self, strand1: ResidueGroup, strand2: ResidueGroup,
                 o2_name: str = 'O2', c2_name: str = "C2",
                 g_name: str = 'G', a_name: str = 'A', u_name: str = 'U',
                 t_name: str = 'T', c_name: str = 'C',
                 **kwargs) -> None:

        selections: Tuple[List[mda.AtomGroup], List[mda.AtomGroup]] = \
            self.select_strand_atoms(
                strand1, strand2, c2_name, o2_name,
                g_name=g_name, a_name=a_name,
                t_name=t_name, u_name=u_name, c_name=c_name
        )

        super(MinorPairDist, self).__init__(
            selections[0], selections[1], **kwargs
        )


class MajorPairDist(NucPairDist):
    r"""Minor-Pair base pair distance for
    selected residues over a trajectory.

    Takes two :class:`~MDAnalysis.core.groups.ResidueGroup` objects and
    calculates the Major-groove hydrogen bond length between the nitrogen
    and oxygen atoms over the trajectory. Bases are matched by their index
    in the two :class:`~MDAnalysis.core.groups.ResidueGroup`
    provided as arguments.

    Parameters
    ----------
    strand1: List[Residue]
        First list of bases
    strand2: List[Residue]
        Second list of bases
    o6_name: str (optional)
        Name of Oxygen 6 of nucleic acids;
        by default assigned to "O6"
    n4_name: str (optional)
        Name of Nitrogen 4 of nucleic acids;
        by default assigned to "N4"
    g_name: str (optional)
        Name of Guanine in topology;
        by default assigned to "G"
    a_name: str (optional)
        Name of Adenine in topology;
        by default assigned to "A"
    u_name: str (optional)
        Name of Uracil in topology;
        by default assigned to "U"
    t_name: str (optional)
        Name of Thymine in topology;
        by default assigned to "T"
    c_name: str (optional)
        Name of Cytosine in topology;
        by default assigned to "C"
    **kwargs:
        arguments for :class:`~MDAnalysis.analysis.base.AnalysisBase`

    Attributes
    ----------
    results.distances: numpy.ndarray
        Distances are stored in a 2d numpy array with axis 0 (first index)
        indexing the trajectory frame and axis 1 (second index) selecting the
        Residue pair.
    times: numpy.ndarray
        Simulation times for analysis.

    Raises
    ------
    ValueError
        A :class:`~MDAnalysis.core.Residue`
        in one of the strands not a valid nucleic acid
    ValueError
        If a given residue pair from the provided strands returns an empty
        :class:`~MDAnalysis.core.groups.AtomGroup` when selecting the atom
        pairs used in the distance calculations
    ValueError
        if the selections given are not the same length


    .. versionadded:: 2.7.0
    """

    def __init__(self, strand1: ResidueGroup, strand2: ResidueGroup,
                 n4_name: str = 'N4', o6_name: str = "O6",
                 g_name: str = 'G', a_name: str = 'A', u_name: str = 'U',
                 t_name: str = 'T', c_name: str = 'C',
                 **kwargs) -> None:

        selections: Tuple[List[mda.AtomGroup], List[mda.AtomGroup]] = \
            self.select_strand_atoms(
                strand1, strand2, o6_name, n4_name, g_name=g_name,
                a_name=a_name, t_name=t_name, u_name=u_name,
                c_name=c_name
        )

        super(MajorPairDist, self).__init__(
            selections[0], selections[1], **kwargs
        )
