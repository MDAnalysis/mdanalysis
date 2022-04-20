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
Linear Density --- :mod:`MDAnalysis.analysis.lineardensity`
===========================================================

A tool to compute mass and charge density profiles along the three
cartesian axes [xyz] of the simulation cell. Works only for orthorombic,
fixed volume cells (thus for simulations in canonical NVT ensemble).
"""

import numpy as np
from MDAnalysis.analysis.base import AnalysisBase, Results


class LinearDensity(AnalysisBase):
    """Linear density profile

    Parameters
    ----------
    select : AtomGroup
          any atomgroup
    grouping : str {'atoms', 'residues', 'segments', 'fragments'}
          Density profiles will be computed either on the atom positions (in
          the case of 'atoms') or on the center of mass of the specified
          grouping unit ('residues', 'segments', or 'fragments').
    binsize : float
          Bin width in Angstrom used to build linear density
          histograms. Defines the resolution of the resulting density
          profile (smaller --> higher resolution)

    Attributes
    ----------
    results.x.dim : int
           index of the [xyz] axes
    results.x.pos : numpy.ndarray
           mass density in [xyz] direction
    results.x.pos_std : numpy.ndarray
           standard deviation of the mass density in [xyz] direction
    results.x.char : numpy.ndarray
           charge density in [xyz] direction
    results.x.char_std : numpy.ndarray
           standard deviation of the charge density in [xyz] direction
    results.x.slice_volume : float
           volume of bin in [xyz] direction

    Example
    -------
    First create a :class:`LinearDensity` object by supplying a selection,
    then use the :meth:`run` method. Finally access the results
    stored in results, i.e. the mass density in the x direction.

    .. code-block:: python

       ldens = LinearDensity(selection)
       ldens.run()
       print(ldens.results.x.pos)


    Alternatively, other types of grouping can be selected using the
    ``grouping`` keyword. For example to calculate the density based on
    a grouping of the :class:`~MDAnalysis.core.groups.ResidueGroup`
    of the input :class:`~MDAnalysis.core.groups.AtomGroup`.

    .. code-block:: python

       ldens = LinearDensity(selection, grouping='residues', binsize=1.0)
       ldens.run()



    .. versionadded:: 0.14.0

    .. versionchanged:: 1.0.0
       Support for the ``start``, ``stop``, and ``step`` keywords has been
       removed. These should instead be passed to :meth:`LinearDensity.run`.
       The ``save()`` method was also removed, you can use ``np.savetxt()`` or
       ``np.save()`` on the :attr:`LinearDensity.results` dictionary contents
       instead.

    .. versionchanged:: 1.0.0
       Changed `selection` keyword to `select`

    .. versionchanged:: 2.0.0
       Results are now instances of
       :class:`~MDAnalysis.core.analysis.Results` allowing access
       via key and attribute.

    .. versionchanged:: 2.2.0
       Fixed a bug that caused LinearDensity to fail if grouping="residues"
       or grouping="segments" were set.
       Residues, segments, and fragments will be analysed based on their centre
       of mass, not centre of geometry as previously stated.
    """

    def __init__(self, select, grouping='atoms', binsize=0.25, **kwargs):
        super(LinearDensity, self).__init__(select.universe.trajectory,
                                            **kwargs)
        # allows use of run(parallel=True)
        self._ags = [select]
        self._universe = select.universe

        self.binsize = binsize

        # group of atoms on which to compute the COM (same as used in
        # AtomGroup.wrap())
        self.grouping = grouping

        # Initiate result instances
        self.results["x"] = Results(dim=0)
        self.results["y"] = Results(dim=1)
        self.results["z"] = Results(dim=2)
        # Box sides
        self.dimensions = self._universe.dimensions[:3]
        self.volume = np.prod(self.dimensions)
        # number of bins
        bins = (self.dimensions // self.binsize).astype(int)

        # Here we choose a number of bins of the largest cell side so that
        # x, y and z values can use the same "coord" column in the output file
        self.nbins = bins.max()
        slices_vol = self.volume / bins

        self.keys = ['pos', 'pos_std', 'char', 'char_std']

        # Initialize results array with zeros
        for dim in self.results:
            idx = self.results[dim]['dim']
            self.results[dim]['slice_volume'] = slices_vol[idx]
            for key in self.keys:
                self.results[dim][key] = np.zeros(self.nbins)

        # Variables later defined in _prepare() method
        self.masses = None
        self.charges = None
        self.totalmass = None

    def _prepare(self):
        # Get masses and charges for the selection
        if self.grouping == "atoms":
            self.masses = self._ags[0].masses
            self.charges = self._ags[0].charges

        elif self.grouping in ["residues", "segments", "fragments"]:
            self.masses = self._ags[0].total_mass(compound=self.grouping)
            self.charges = self._ags[0].total_charge(compound=self.grouping)

        else:
            raise AttributeError(
                f"{self.grouping} is not a valid value for grouping.")

        self.totalmass = np.sum(self.masses)

    def _single_frame(self):
        self.group = getattr(self._ags[0], self.grouping)
        self._ags[0].wrap(compound=self.grouping)

        # Find position of atom/group of atoms
        if self.grouping == 'atoms':
            positions = self._ags[0].positions  # faster for atoms
        else:
            # Centre of mass for residues, segments, fragments
            positions = self._ags[0].center_of_mass(compound=self.grouping)

        for dim in ['x', 'y', 'z']:
            idx = self.results[dim]['dim']

            key = 'pos'
            key_std = 'pos_std'
            # histogram for positions weighted on masses
            hist, _ = np.histogram(positions[:, idx],
                                   weights=self.masses,
                                   bins=self.nbins,
                                   range=(0.0, max(self.dimensions)))

            self.results[dim][key] += hist
            self.results[dim][key_std] += np.square(hist)

            key = 'char'
            key_std = 'char_std'
            # histogram for positions weighted on charges
            hist, _ = np.histogram(positions[:, idx],
                                   weights=self.charges,
                                   bins=self.nbins,
                                   range=(0.0, max(self.dimensions)))

            self.results[dim][key] += hist
            self.results[dim][key_std] += np.square(hist)

    def _conclude(self):
        k = 6.022e-1  # divide by avodagro and convert from A3 to cm3

        # Average results over the number of configurations
        for dim in ['x', 'y', 'z']:
            for key in ['pos', 'pos_std', 'char', 'char_std']:
                self.results[dim][key] /= self.n_frames
            # Compute standard deviation for the error
            # For certain tests in testsuite, floating point imprecision
            # can lead to negative radicands of tiny magnitude (yielding nan).
            # radicand_pos and radicand_char are therefore calculated first and
            # and negative values set to 0 before the square root
            # is calculated.
            radicand_pos = self.results[dim][
                'pos_std'] - np.square(self.results[dim]['pos'])
            radicand_pos[radicand_pos < 0] = 0
            self.results[dim]['pos_std'] = np.sqrt(radicand_pos)

            radicand_char = self.results[dim][
                'char_std'] - np.square(self.results[dim]['char'])
            radicand_char[radicand_char < 0] = 0
            self.results[dim]['char_std'] = np.sqrt(radicand_char)

        for dim in ['x', 'y', 'z']:
            norm = k * self.results[dim]['slice_volume']
            for key in self.keys:
                self.results[dim][key] /= norm

    def _add_other_results(self, other):
        # For parallel analysis
        results = self.results
        for dim in ['x', 'y', 'z']:
            key = 'pos'
            key_std = 'pos_std'
            results[dim][key] += other[dim][key]
            results[dim][key_std] += other[dim][key_std]

            key = 'char'
            key_std = 'char_std'
            results[dim][key] += other[dim][key]
            results[dim][key_std] += other[dim][key_std]
