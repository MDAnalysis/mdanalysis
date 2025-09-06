# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
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

:Author: Alia Lescoulie, with contributions from the MDAnalysis team
:Year: 2022-2025
:copyright: LGPLv2.1

The module provides classes for analyzing nucleic acids structures, including
both double-stranded and single-stranded nucleic acids, as well as support for
common RNA modifications. This is an updated, higher performance version of 
previous nucleic acid tools. For applications see :footcite:p:`Denning2011,Denning2012`.

.. versionchanged:: 3.0.0
   - Added support for single-stranded nucleic acids
   - Added support for common RNA modifications
   - Improved error handling and validation
   - Added visualization integration examples

.. rubric:: Examples

Analyzing double-stranded DNA::

    import MDAnalysis as mda
    from MDAnalysis.analysis import nucleicacids
    
    # Load a DNA structure
    u = mda.Universe("dna.pdb")
    
    # Select the two strands (assuming they're named 'strand1' and 'strand2')
    strand1 = u.select_atoms("segid A")
    strand2 = u.select_atoms("segid B")
    
    # Calculate Watson-Crick distances
    wc_dist = nucleicacids.WatsonCrickDist(
        strand1.residues, 
        strand2.residues
    )
    wc_dist.run()
    
    # Calculate minor groove widths
    minor_dist = nucleicacids.MinorPairDist(
        strand1.residues,
        strand2.residues
    )
    minor_dist.run()

Analyzing single-stranded RNA with modifications::

    # Load an RNA structure with modifications
    u = mda.Universe("rna_with_mods.pdb")
    
    # Select the single strand (may contain modified residues like 7MG, PSU, etc.)
    rna = u.select_atoms("segid A")
    
    # Analyze as single-stranded
    wc_dist = nucleicacids.WatsonCrickDist(
        rna.residues,
        None  # Indicates single-stranded analysis
    )
    wc_dist.run()
    
    # The distances will be between consecutive residues
    print("Distances between consecutive residues:", wc_dist.results.distances)

Visualization with NGLView::

    import nglview as nv
    
    # Create a view of the structure
    view = nv.show_mdanalysis(u)
    
    # Add distance measurements for visualization
    for i, (res1, res2) in enumerate(zip(strand1.residues, strand2.residues)):
        # Get the atom pairs used in the analysis
        n1 = res1.atoms.select_atoms("name N1")[0]
        n3 = res2.atoms.select_atoms("name N3")[0]
        
        # Add distance measurement
        view.add_distance(
            atom_pair=[[n1.ix, n3.ix]],
            label_type="label+length",
            color="red",
            label_color="white"
        )
    
    # Display the view
    view

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

from typing import List, Dict, Tuple, Union, Optional, Callable
import warnings
import numpy as np

import MDAnalysis as mda
from .distances import calc_bonds
from .base import AnalysisBase, ResultsGroup
from MDAnalysis.core.groups import Residue, ResidueGroup, AtomGroup
from functools import partial


# Deprecation: In 3.0.0 change type to just
# ResidueClass = ResidueGroup
ResidueClass = Union[List[Residue], ResidueGroup]
r"""A type alias for :code:`Union[List[Residue], ResidueGroup]`

Used as an alias for methods where either class is acceptable.
"""


class NucPairDist(AnalysisBase):
    r"""Atom pair distance calculation base class.

    Takes two lists of :class:`~MDAnalysis.core.groups.AtomGroup` and
    computes the distances between them over a trajectory. Can also handle
    single-stranded nucleic acids when strand2 is None. Used as a superclass
    for the other nucleic acid distances classes. The distance will be measured
    between atoms sharing an index in the two lists of
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

    .. versionchanged:: 2.9.0
       Enabled **parallel execution** with the ``multiprocessing`` and ``dask``
       backends; use the new method :meth:`get_supported_backends` to see all
       supported backends.
    """

    _analysis_algorithm_is_parallelizable = True

    @classmethod
    def get_supported_backends(cls):
        return ("serial", "multiprocessing", "dask")

    _s1: mda.AtomGroup
    _s2: mda.AtomGroup
    _n_sel: int

    def __init__(
        self,
        selection1: List[mda.AtomGroup],
        selection2: Optional[List[mda.AtomGroup]] = None,
        **kwargs,
    ) -> None:
        super(NucPairDist, self).__init__(
            selection1[0].universe.trajectory, **kwargs
        )

        self._single_strand = selection2 is None
        
        if not self._single_strand and len(selection1) != len(selection2):
            raise ValueError("When using two strands, selections must be same length")

        self._n_sel: int = len(selection1)

        self._s1 = selection1[0]
        self._s2 = selection1[0]  # Default to self if single strand
        
        if not self._single_strand:
            self._s2 = selection2[0]

        for i in range(1, self._n_sel):
            self._s1 += selection1[i]
            if not self._single_strand:
                self._s2 += selection2[i]

    @staticmethod
    def select_strand_atoms(
        strand1: ResidueGroup,
        strand2: Optional[ResidueGroup] = None,
        a1_name: str = "",
        a2_name: str = "",
        g_name: str = "G",
        a_name: str = "A",
        u_name: str = "U",
        t_name: str = "T",
        c_name: str = "C",
        # Common RNA modifications
        m7g_name: str = "7MG",  # 7-methylguanosine
        psu_name: str = "PSU",  # Pseudouridine
        m5c_name: str = "5MC",  # 5-methylcytidine
        m6a_name: str = "M6A",  # N6-methyladenosine
        i_name: str = "I",      # Inosine
        omc_name: str = "OMC",  # 5-methylcytidine (alternative)
        um_name: str = "5MU",   # 5-methyluridine
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
        pyrimidines: List[str] = [c_name, t_name, u_name, m5c_name, omc_name, um_name, psu_name]
        purines: List[str] = [a_name, g_name, m7g_name, m6a_name, i_name]
        
        # Map of common RNA modifications to their canonical bases
        mod_to_canonical = {
            m7g_name: g_name,
            psu_name: u_name,
            m5c_name: c_name,
            m6a_name: a_name,
            i_name: g_name,  # Inosine pairs with C, similar to G
            omc_name: c_name,
            um_name: u_name
        }

        sel1: List[mda.AtomGroup] = []
        sel2: List[mda.AtomGroup] = []
        
        # Handle single strand case
        if strand2 is None:
            for res in strand1.residues:
                # For single strand, just select the first atom for distance calculation
                ag = res.atoms[0:1]  # Take first atom of each residue
                sel1.append(ag)
            return (sel1, [])
            
        # Handle double strand case
        for res1, res2 in zip(strand1.residues, strand2.residues):
            # Check for modifications and map to canonical bases
            res1_name = res1.resname.strip()
            res2_name = res2.resname.strip()
            
            # Get canonical base name if modified
            base1 = mod_to_canonical.get(res1_name, res1_name[0] if res1_name else '')
            base2 = mod_to_canonical.get(res2_name, res2_name[0] if res2_name else '')
            
            # Determine atom names based on base pairing
            if base1[0] in pyrimidines and base2[0] in purines:
                a1, a2 = a2_name, a1_name
            elif base1[0] in purines and base2[0] in pyrimidines:
                a1, a2 = a1_name, a2_name
            else:
                raise ValueError(
                    f"Invalid base pair: {res1_name} and {res2_name} are not complementary"
                )

            ag1 = res1.atoms.select_atoms(f"name {a1}")
            ag2 = res2.atoms.select_atoms(f"name {a2}")

            if not ag1 or not ag2:
                # Try alternative atom names for modified bases
                alt_ag1 = res1.atoms[0:1]  # Fallback to first atom
                alt_ag2 = res2.atoms[0:1]
                
                if not ag1 and not alt_ag1:
                    raise ValueError(
                        f"Could not select atom {a1} from residue {res1.resname} (resid {res1.resid})"
                    )
                if not ag2 and not alt_ag2:
                    raise ValueError(
                        f"Could not select atom {a2} from residue {res2.resname} (resid {res2.resid})"
                    )
                
                ag1 = ag1 if ag1 else alt_ag1
                ag2 = ag2 if ag2 else alt_ag2

            sel1.append(ag1)
            sel2.append(ag2)

        return (sel1, sel2)

    def _prepare(self) -> None:
        self.results.distances: np.ndarray = np.zeros(
            [self.n_frames, self._n_sel]
        )

    def _single_frame(self) -> None:
        if self._single_strand:
            # For single strand, calculate distances between consecutive residues
            dist = np.array([
                np.linalg.norm(p1 - p2) 
                for p1, p2 in zip(self._s1.positions[:-1], self._s1.positions[1:])
            ])
        else:
            # For double strand, calculate distances between pairs
            dist = calc_bonds(self._s1.positions, self._s2.positions)

        self.results.distances[self._frame_index, :len(dist)] = dist

    def _conclude(self) -> None:
        self.results["pair_distances"] = self.results["distances"]
        # TODO: remove pair_distances in 3.0.0

    def _get_aggregator(self):
        return ResultsGroup(
            lookup={
                "distances": ResultsGroup.ndarray_vstack,
            }
        )


class WatsonCrickDist(NucPairDist):
    r"""
    Watson-Crick base pair distance for selected
    residues over a trajectory.

    Takes one or two :class:`~MDAnalysis.core.groups.ResidueGroup`
    objects and calculates the distance between the nitrogen atoms in the
    Watson-Crick hydrogen bond over the trajectory. 
    
    For double-stranded nucleic acids, bases are matched by their index in the two
    :class:`~MDAnalysis.core.groups.ResidueGroup` provided as arguments.
    For single-stranded nucleic acids, pass only `strand1` and set `strand2=None`.

    .. note::
        Support for :class:`~MDAnalysis.core.groups.Residue` is slated for
        deprecation and will raise a warning when used. It still works but
        :class:`~MDAnalysis.core.groups.ResidueGroup` is preferred.

    Parameters
    ----------
    strand1: ResidueClass
        First list of bases (required)

        .. deprecated:: 2.7.0
           Using a list of :class:`~MDAnalysis.core.groups.Residue` will
           be removed in 3.0.0. Use a
           :class:`~MDAnalysis.core.groups.ResidueGroup`.

    strand2: ResidueClass, optional
        Second list of bases. If None, treats as single-stranded nucleic acid.
        Default is None.

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

    def __init__(
        self,
        strand1: ResidueClass,
        strand2: Optional[ResidueClass] = None,
        n1_name: str = "N1",
        n3_name: str = "N3",
        g_name: str = "G",
        a_name: str = "A",
        u_name: str = "U",
        t_name: str = "T",
        c_name: str = "C",
        **kwargs,
    ) -> None:
        def verify_strand(strand: ResidueClass) -> Optional[ResidueGroup]:
            # Helper method to verify the strands
            if strand is None:
                return None
                
            if isinstance(strand, list):  # Checking if a list is given
                # verify list is only Residues
                if not all(isinstance(resid, Residue) for resid in strand):
                    raise TypeError(f"{strand} contains non-Residue elements")

                warnings.warn(
                    DeprecationWarning(
                        (
                            f"ResidueGroup should be used for {strand} instead "
                            "of giving a Residue list"
                        )
                    )
                )
                # Convert to a ResidueGroup
                strand = ResidueGroup(strand)

            return strand

        strand1_group: ResidueGroup = verify_strand(strand1)
        strand2_group: Optional[ResidueGroup] = verify_strand(strand2)

        # For single strand, we'll use the base class's single-strand handling
        if strand2_group is None:
            # For single strand, we'll just use the base class's single-strand handling
            # by passing the same group as both arguments
            strand_atomgroups: Tuple[List[mda.AtomGroup], List[mda.AtomGroup]] = (
                self.select_strand_atoms(
                    strand1_group,
                    None,  # Indicates single strand
                    n1_name,
                    n3_name,
                    g_name=g_name,
                    a_name=a_name,
                    t_name=t_name,
                    u_name=u_name,
                    c_name=c_name,
                )
            )
            # For single strand, we pass the same group as both arguments
            # The base class will handle it appropriately
            super(WatsonCrickDist, self).__init__(
                strand_atomgroups[0], None, **kwargs
            )
        else:
            # For double strand, proceed as before
            strand_atomgroups = self.select_strand_atoms(
                strand1_group,
                strand2_group,
                n1_name,
                n3_name,
                g_name=g_name,
                a_name=a_name,
                t_name=t_name,
                u_name=u_name,
                c_name=c_name,
            )
            super(WatsonCrickDist, self).__init__(
                strand_atomgroups[0], strand_atomgroups[1], **kwargs
            )


class MinorPairDist(NucPairDist):
    r"""Minor-Pair basepair distance for selected residues over a trajectory.

    Takes one or two :class:`~MDAnalysis.core.groups.ResidueGroup` objects and
    calculates the Minor-groove hydrogen bond length between the
    nitrogen and oxygen atoms over the trajectory. 
    
    For double-stranded nucleic acids, bases are matched by their index in the two
    :class:`~MDAnalysis.core.groups.ResidueGroup` provided as arguments.
    For single-stranded nucleic acids, pass only `strand1` and set `strand2=None`.

    Parameters
    ----------
    strand1: ResidueClass
        First list of bases (required)
    strand2: ResidueClass, optional
        Second list of bases. If None, treats as single-stranded nucleic acid.
        Default is None.
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

    def __init__(
        self,
        strand1: ResidueClass,
        strand2: Optional[ResidueClass] = None,
        o2_name: str = "O2",
        c2_name: str = "C2",
        g_name: str = "G",
        a_name: str = "A",
        u_name: str = "U",
        t_name: str = "T",
        c_name: str = "C",
        **kwargs,
    ) -> None:
        def verify_strand(strand: ResidueClass) -> Optional[ResidueGroup]:
            # Helper method to verify the strands
            if strand is None:
                return None
                
            if isinstance(strand, list):  # Checking if a list is given
                # verify list is only Residues
                if not all(isinstance(resid, Residue) for resid in strand):
                    raise TypeError(f"{strand} contains non-Residue elements")

                warnings.warn(
                    DeprecationWarning(
                        (
                            f"ResidueGroup should be used for {strand} instead "
                            "of giving a Residue list"
                        )
                    )
                )
                # Convert to a ResidueGroup
                strand = ResidueGroup(strand)

            return strand

        strand1_group: ResidueGroup = verify_strand(strand1)
        strand2_group: Optional[ResidueGroup] = verify_strand(strand2)

        if strand2_group is None:
            # For single strand, use the base class's single-strand handling
            strand_atomgroups: Tuple[List[mda.AtomGroup], List[mda.AtomGroup]] = (
                self.select_strand_atoms(
                    strand1_group,
                    None,  # Indicates single strand
                    c2_name,
                    o2_name,
                    g_name=g_name,
                    a_name=a_name,
                    t_name=t_name,
                    u_name=u_name,
                    c_name=c_name,
                )
            )
            super(MinorPairDist, self).__init__(
                strand_atomgroups[0], None, **kwargs
            )
        else:
            # For double strand, proceed as before
            strand_atomgroups = self.select_strand_atoms(
                strand1_group,
                strand2_group,
                c2_name,
                o2_name,
                g_name=g_name,
                a_name=a_name,
                t_name=t_name,
                u_name=u_name,
                c_name=c_name,
            )
            super(MinorPairDist, self).__init__(
                strand_atomgroups[0], strand_atomgroups[1], **kwargs
            )


class MajorPairDist(NucPairDist):
    r"""Major-Pair base pair distance for
    selected residues over a trajectory.

    Takes one or two :class:`~MDAnalysis.core.groups.ResidueGroup` objects and
    calculates the Major-groove hydrogen bond length between the nitrogen
    and oxygen atoms over the trajectory. 
    
    For double-stranded nucleic acids, bases are matched by their index in the two
    :class:`~MDAnalysis.core.groups.ResidueGroup` provided as arguments.
    For single-stranded nucleic acids, pass only `strand1` and set `strand2=None`.

    Parameters
    ----------
    strand1: ResidueClass
        First list of bases (required)
    strand2: ResidueClass, optional
        Second list of bases. If None, treats as single-stranded nucleic acid.
        Default is None.
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

    def __init__(
        self,
        strand1: ResidueClass,
        strand2: Optional[ResidueClass] = None,
        n4_name: str = "N4",
        o6_name: str = "O6",
        g_name: str = "G",
        a_name: str = "A",
        u_name: str = "U",
        t_name: str = "T",
        c_name: str = "C",
        **kwargs,
    ) -> None:
        def verify_strand(strand: ResidueClass) -> Optional[ResidueGroup]:
            # Helper method to verify the strands
            if strand is None:
                return None
                
            if isinstance(strand, list):  # Checking if a list is given
                # verify list is only Residues
                if not all(isinstance(resid, Residue) for resid in strand):
                    raise TypeError(f"{strand} contains non-Residue elements")

                warnings.warn(
                    DeprecationWarning(
                        (
                            f"ResidueGroup should be used for {strand} instead "
                            "of giving a Residue list"
                        )
                    )
                )
                # Convert to a ResidueGroup
                strand = ResidueGroup(strand)

            return strand

        strand1_group: ResidueGroup = verify_strand(strand1)
        strand2_group: Optional[ResidueGroup] = verify_strand(strand2)

        if strand2_group is None:
            # For single strand, use the base class's single-strand handling
            strand_atomgroups: Tuple[List[mda.AtomGroup], List[mda.AtomGroup]] = (
                self.select_strand_atoms(
                    strand1_group,
                    None,  # Indicates single strand
                    o6_name,
                    n4_name,
                    g_name=g_name,
                    a_name=a_name,
                    t_name=t_name,
                    u_name=u_name,
                    c_name=c_name,
                )
            )
            super(MajorPairDist, self).__init__(
                strand_atomgroups[0], None, **kwargs
            )
        else:
            # For double strand, proceed as before
            strand_atomgroups = self.select_strand_atoms(
                strand1_group,
                strand2_group,
                o6_name,
                n4_name,
                g_name=g_name,
                a_name=a_name,
                t_name=t_name,
                u_name=u_name,
                c_name=c_name,
            )
            super(MajorPairDist, self).__init__(
                strand_atomgroups[0], strand_atomgroups[1], **kwargs
            )
