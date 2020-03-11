# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2020 The MDAnalysis Development Team and contributors
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

"""Secondary structure analysis --- :mod:`MDAnalysis.analysis.secondary_structure`
==================================================================================

:Authors: Lily Wang
:Year: 2020
:Copyright: GNU Public License v3

.. versionadded:: 1.0.0


"""

from __future__ import absolute_import, division

import numpy as np

from ...core.topologyattrs import ResidueAttr
from ..base import AnalysisBase


class SecondaryStructure(ResidueAttr):
    """Single letter code for secondary structure"""
    attrname = 'secondary_structures'
    singular = 'secondary_structure'
    dtype = '<U1'

    @staticmethod
    def _gen_initial_values(na, nr, ns):
        return np.empty(nr, dtype='<U1')


class SecondaryStructureBase(AnalysisBase):
    """Base class for implementing secondary structure analysis.

    Subclasses should implement ``_compute_dssp``.

    Parameters
    ----------
    universe: Universe or AtomGroup
        The Universe or AtomGroup to apply the analysis to. As secondary
        structure is a residue property, the analysis is applied to every
        residue in your chosen atoms.
    select: string, optional
        The selection string for selecting atoms from ``universe``. The
        analysis is applied to the residues in this subset of atoms.
    add_topology_attr: bool, optional
        Whether to add the most common secondary structure as a topology
        attribute ``secondary_structure`` to your residues.
    verbose: bool, optional
        Turn on more logging and debugging.


    Attributes
    ----------
    residues: :class:`~MDAnalysis.core.groups.ResidueGroup`
        The residues to which the analysis is applied.
    ss_codes: :class:`numpy.ndarray` of shape (n_frames, n_residues)
        Single-letter secondary structure codes
    ss_names: :class:`numpy.ndarray` of shape (n_frames, n_residues)
        Secondary structure names
    ss_simple: :class:`numpy.ndarray` of shape (n_frames, n_residues)
        Simplified secondary structure names
    ss_counts: dict of {code: :class:`numpy.ndarray` of shape (n_frames,)}
        Dictionary that counts number of residues with each secondary
        structure code, for each frame
    simple_counts: dict of {name: :class:`numpy.ndarray` of shape (n_frames,)}
        Dictionary that counts number of residues with each simplified
        secondary structure, for each frame
    ss_mode: :class:`numpy.ndarray` of shape (n_residues,)
        The most common secondary structure for each residue
    ss_codes_to_names: dict of {code: name}
        Dictionary converting each single-letter code to the full name of
        the secondary structure
    ss_codes_to_simple: dict of {code: name}
        Dictionary converting each single-letter code to simplified
        secondary structures

    """

    ss_codes_to_names = {
        'G': 'Helix (3-10)',
        'H': 'Helix (alpha)',
        'I': 'Helix (pi)',
        'E': 'Strand',
        'B': 'Bridge',
        'T': 'Turn',
        'S': 'Bend',
        'C': 'Coil',
        '': '',
    }

    ss_codes_to_simple = {
        'G': 'Helix',
        'H': 'Helix',
        'I': 'Helix',
        'E': 'Strand',
        'B': 'Strand',
        'T': 'Coil',
        'S': 'Coil',
        'C': 'Coil',
        '': '',
    }

    def __init__(self, universe, select='backbone', verbose=False,
                 add_topology_attr=False):
        super(SecondaryStructureBase, self).__init__(universe.universe.trajectory,
                                                     verbose=verbose)
        self._universe = universe.universe
        self.atomgroup = universe.select_atoms(select)
        self.residues = self.atomgroup.residues
        self.n_residues = len(self.residues)
        self._add_topology_attr = add_topology_attr

    def _prepare(self):
        self.frames = np.zeros(self.n_frames, dtype=int)
        self.ss_codes = np.empty((self.n_frames, self.n_residues),
                                 dtype='<U1')
        self.ss_counts = dict.fromkeys(self.ss_codes_to_simple)
        self.ss_counts.pop('')
        for k in self.ss_counts:
            self.ss_counts[k] = np.zeros(self.n_frames, dtype=int)

        simple_ids = set(self.ss_codes_to_simple.values())
        self.simple_counts = dict.fromkeys(simple_ids)
        self.simple_counts.pop('')
        for k in self.simple_counts:
            self.simple_counts[k] = np.zeros(self.n_frames, dtype=int)

    def _compute_dssp(self):
        raise NotImplementedError

    def _single_frame(self):
        self.frames[self._frame_index] = self._ts.frame
        self._compute_dssp()
        row = self.ss_codes[self._frame_index]
        codes, counts = np.unique(row[np.nonzero(row)[0]], return_counts=True)
        for c, n in zip(codes, counts):
            self.ss_counts[c][self._frame_index] = n

    def _conclude(self):
        # convert to full names and simple codes
        codes, idx = np.unique(self.ss_codes, return_inverse=True)
        names = np.array([self.ss_codes_to_names[c] for c in codes])
        self.ss_names = names[idx].reshape(self.ss_codes.shape)
        simp = np.array([self.ss_codes_to_simple[c] for c in codes])
        self.ss_simple = simp[idx].reshape(self.ss_codes.shape)

        # count number of simple structures
        for code, counts in self.ss_counts.items():
            if code:
                simple = self.ss_codes_to_simple[code]
                self.simple_counts[simple] += counts

        # get most common secondary structure
        self.ss_mode = np.empty(len(self.residues), dtype='<U1')
        for i, col in enumerate(self.ss_codes.T):
            code, counts = np.unique(
                col[np.nonzero(col)[0]], return_counts=True)
            try:
                self.ss_mode[i] = code[np.argmax(counts)]
            except ValueError:  # counts is empty
                pass
        simple_mode = [self.ss_codes_to_simple[s] for s in self.ss_mode]
        self.simple_mode = np.array(simple_mode)

        # possibly add this as attribute
        if self._add_topology_attr:
            attr = SecondaryStructure.attrname
            self._universe.add_TopologyAttr(attr)
            setattr(self.residues, attr, self.ss_mode)

    def plot_content(self, ax=None, simple=False, figsize=(8, 6),
                     n_xticks=10, kind='bar'):
        """Plot the counts of secondary structures over frames.

        Parameters
        ----------
        ax: :class: `matplotlib.axes.Axes`, optional
            If no `ax` is supplied or set to ``None`` then the plot will
            be created on new axes.
        simple: bool, optional
            If ``True``, plots the counts of the simple secondary
            structures. If ``False``, plots the counts of each distinct
            category.

        Returns
        -------
        ax: :class: `matplotlib.axes.Axes`

        """
        try:
            import pandas as pd
        except ImportError:
            raise ImportError('pandas is required for this function.')
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError('matplotlib is required for this function.')

        if not simple:
            data = self.ss_counts
        else:
            data = self.simple_counts

        df = pd.DataFrame(data, index=self.frames)
        df = 100 * df/(df.sum(axis=1).max())
        if kind == 'bar':
            ax = df.plot(kind='bar', stacked=True, ax=ax, figsize=figsize)
        else:
            ax = df.plot(kind='line', ax=ax, figsize=figsize)

        # float for python 2.7
        step = int(np.ceil(float(self.n_frames)/n_xticks))
        step = max([1, step])
        idx = np.arange(0, self.n_frames, step, dtype=int)
        labels = self.frames[idx]
        plt.xticks(idx, labels)
        plt.ylabel('% secondary structure')
        plt.xlabel('Frame')
        return ax
