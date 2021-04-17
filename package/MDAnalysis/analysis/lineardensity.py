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
cartesian axes of the simulation cell. Works only for orthorombic,
fixed volume cells (thus for simulations in canonical NVT ensemble).
"""
import os.path as path

import numpy as np

from MDAnalysis.analysis.base import AnalysisBase


class LinearDensity(AnalysisBase):
    """Linear density profile

    Parameters
    ----------
    select : AtomGroup
          any atomgroup
    grouping : str {'atoms', 'residues', 'segments', 'fragments'}
          Density profiles will be computed on the center of geometry
          of a selected group of atoms
    binsize : float
          Bin width in Angstrom used to build linear density
          histograms. Defines the resolution of the resulting density
          profile (smaller --> higher resolution)
    verbose : bool (optional)
          Show detailed progress of the calculation if set to ``True``; the
          default is ``False``.

    Attributes
    ----------
    results.x_pos : numpy.ndarray
           mass density in x direction
    results.x_pos_std : numpy.ndarray
           standard deviation of the mass density in x direction
    results.x_char : numpy.ndarray
           charge density in x direction
    results.x_char_std : numpy.ndarray
           standard deviation of the charge density in x direction
    results.x_slice_volume : float
           volume of bin in x direction
    results.y_pos : numpy.ndarray
           mass density in y direction
    results.y_pos_std : numpy.ndarray
           standard deviation of the mass density in y direction
    results.y_char : numpy.ndarray
           charge density in y direction
    results.y_char_std : numpy.ndarray
           standard deviation of the charge density in y direction
    results.y_slice_volume : float
           volume of bin in y direction
    results.z_pos : numpy.ndarray
           mass density in z direction
    results.z_pos_std : numpy.ndarray
           standard deviation of the mass density in z direction
    results.z_char : numpy.ndarray
           charge density in z direction
    results.z_char_std : numpy.ndarray
           standard deviation of the charge density in z direction
    results.z_slice_volume : float
           volume of bin in z direction


    Example
    -------
    First create a LinearDensity object by supplying a selection,
    then use the :meth:`run` method. Finally access the results
    stored in results, i.e. the mass density in the x direction.

    .. code-block:: python

        ldens = LinearDensity(selection)
        ldens.run()
        print(ldens.results.x_pos)


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
        Changed structure of the the `results` dictionary to an object.
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

        # Box sides
        self.dimensions = self._universe.dimensions[:3]
        self.volume = np.prod(self.dimensions)
        # number of bins
        bins = (self.dimensions // self.binsize).astype(int)

        # Here we choose a number of bins of the largest cell side so that
        # x, y and z values can use the same "coord" column in the output file
        self.nbins = bins.max()
        slices_vol = self.volume / bins

        self._keys = ['pos', 'pos_std', 'char', 'char_std']
        # Create an alias for the results object to save writing
        self._results = self.results.__dict__

        # Initialize results array with zeros
        for idx in [0, 1, 2]:
            dim = "xyz"[idx]
            self._results[f'{dim}_slice volume'] = slices_vol[idx]
            for attr in self._keys:
                self._results[f"{dim}_{attr}"] =  np.zeros(self.nbins)

        # Variables later defined in _prepare() method
        self.masses = None
        self.charges = None
        self.totalmass = None

    def _prepare(self):
        # group must be a local variable, otherwise there will be
        # issues with parallelization
        group = getattr(self._ags[0], self.grouping)

        # Get masses and charges for the selection
        try:  # in case it's not an atom
            self.masses = np.array([elem.total_mass() for elem in group])
            self.charges = np.array([elem.total_charge() for elem in group])
        except AttributeError:  # much much faster for atoms
            self.masses = self._ags[0].masses
            self.charges = self._ags[0].charges

        self.totalmass = np.sum(self.masses)

    def _single_frame(self):
        self.group = getattr(self._ags[0], self.grouping)
        self._ags[0].wrap(compound=self.grouping)

        # Find position of atom/group of atoms
        if self.grouping == 'atoms':
            positions = self._ags[0].positions  # faster for atoms
        else:
            # COM for res/frag/etc
            positions = np.array([elem.centroid() for elem in self.group])

        for idx in [0, 1, 2]:
            dim = 'xyz'[idx]

            key = 'pos'
            key_std = 'pos_std'
            # histogram for positions weighted on masses
            hist, _ = np.histogram(positions[:, idx],
                                   weights=self.masses,
                                   bins=self.nbins,
                                   range=(0.0, max(self.dimensions)))

            self._results[f"{dim}_{key}"] += hist
            self._results[f"{dim}_{key_std}"] += np.square(hist)

            key = 'char'
            key_std = 'char_std'
            # histogram for positions weighted on charges
            hist, _ = np.histogram(positions[:, idx],
                                   weights=self.charges,
                                   bins=self.nbins,
                                   range=(0.0, max(self.dimensions)))

            self._results[f"{dim}_{key}"] += hist
            self._results[f"{dim}_{key_std}"] += np.square(hist)

    def _conclude(self):
        k = 6.022e-1  # divide by avodagro and convert from A3 to cm3

        # Average results over the  number of configurations
        for dim in ['x', 'y', 'z']:
            for key in self._keys:
                self._results[f"{dim}_{key}"] /= self.n_frames
            # Compute standard deviation for the error
            self._results[f'{dim}_pos_std'] = np.sqrt(
                self._results[f'{dim}_pos_std'] -
                np.square(self._results[f'{dim}_pos']))
            self._results[f'{dim}_char_std'] = np.sqrt(
                self._results[f'{dim}_char_std'] -
                np.square(self._results[f'{dim}_char']))

        for dim in ['x', 'y', 'z']:
            norm = k * self._results[f'{dim}_slice volume']
            for attr in self._keys:
                self._results[f'{dim}_{attr}'] /= norm

    def _add_other_results(self, other):
        # For parallel analysis
        for dim in ['x', 'y', 'z']:
            key = 'pos'
            key_std = 'pos_std'
            self._results[f"{dim}_{key}"] += other[f"{dim}_{key}"]
            self._results[f"{dim}_{key_std}"] += other[f"{dim}_{key_std}"]

            key = 'char'
            key_std = 'char_std'
            self._results[f"{dim}_{key}"] += other[f"{dim}_{key}"]
            self._results[f"{dim}_{key_std}"] += other[f"{dim}_{key}"]
