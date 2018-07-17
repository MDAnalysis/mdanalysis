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
from __future__ import division, absolute_import

import os.path as path

import numpy as np

from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.lib.util import deprecate

class LinearDensity(AnalysisBase):
    """Linear density profile

    Parameters
    ----------
    selection : AtomGroup
          any atomgroup
    grouping : str {'atoms', 'residues', 'segments', 'fragments'}
          Density profiles will be computed on the center of geometry
          of a selected group of atoms ['atoms']
    binsize : float
          Bin width in Angstrom used to build linear density
          histograms. Defines the resolution of the resulting density
          profile (smaller --> higher resolution) [0.25]
    start : int
          The frame to start at [0]
    stop : int
          The frame to end at [-1]
    step : int
          The step size through the trajectory in frames [0]
    verbose : bool (optional)
          Show detailed progress of the calculation if set to ``True``; the
          default is ``False``.

    Attributes
    ----------
    results : dict
          Keys 'x', 'y', and 'z' for the three directions. Under these
          keys, find 'pos', 'pos_std' (mass-weighted density and
          standard deviation), 'char', 'char_std' (charge density and
          its standard deviation), 'slice_volume' (volume of bin).

    Example
    -------
    First create a LinearDensity object by supplying a selection,
    then use the :meth:`run` method::

      ldens = LinearDensity(selection)
      ldens.run()


    .. versionadded:: 0.14.0

    """

    def __init__(self, selection, grouping='atoms', binsize=0.25, **kwargs):
        super(LinearDensity, self).__init__(selection.universe.trajectory,
                                            **kwargs)
        # allows use of run(parallel=True)
        self._ags = [selection]
        self._universe = selection.universe

        self.binsize = binsize

        # group of atoms on which to compute the COM (same as used in
        # AtomGroup.wrap())
        self.grouping = grouping

        # Dictionary containing results
        self.results = {'x': {'dim': 0}, 'y': {'dim': 1}, 'z': {'dim': 2}}
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
            self.results[dim].update({'slice volume': slices_vol[idx]})
            for key in self.keys:
                self.results[dim].update({key: np.zeros(self.nbins)})

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

        # Average results over the  number of configurations
        for dim in ['x', 'y', 'z']:
            for key in ['pos', 'pos_std', 'char', 'char_std']:
                self.results[dim][key] /= self.n_frames
            # Compute standard deviation for the error
            self.results[dim]['pos_std'] = np.sqrt(self.results[dim][
                'pos_std'] - np.square(self.results[dim]['pos']))
            self.results[dim]['char_std'] = np.sqrt(self.results[dim][
                'char_std'] - np.square(self.results[dim]['char']))

        for dim in ['x', 'y', 'z']:
            self.results[dim]['pos'] /= self.results[dim]['slice volume'] * k
            self.results[dim]['char'] /= self.results[dim]['slice volume'] * k
            self.results[dim]['pos_std'] /= self.results[dim]['slice volume'] * k
            self.results[dim]['char_std'] /= self.results[dim]['slice volume'] * k

    @deprecate(release="0.19.0", remove="1.0.0",
               message="Instead save the :attr:`results` dictionary directly in "
               "your favorite format (pickle, json, hdf5, ...).")
    def save(self, description='', form='txt'):
        """Save density profile to file

        Allows to save the density profile to either a ASCII txt file or a
        binary numpy npz file. Output file has extension 'ldens' and begins
        with the name of the trajectory file.

        Parameters
        ----------
        description : str
          An arbitrary description added to the output filename. Can be useful
        form : str {'txt', 'npz'}
          Format of the output. 'txt' will generate an ASCII text file while
          'npz' will produce a numpy binary file.

        Example
        -------
        After initializing and running a `LinearDensity` object, results can be
        written to file as follows::

          ldens.save(description='mydensprof', form='txt')

        which will output the linear density profiles in a file named
        `<trajectory_filename>.mydensprof_<grouping>.ldens`.


        """
        # Take root of trajectory filename for output file naming
        trajname = path.splitext(path.basename(
            self._universe.trajectory.filename))[0]
        # additional string for naming the output file
        description = description + "_" + str(self.grouping)
        filename = trajname + "." + description + ".ldens"

        if form is 'txt':
            self._savetxt(filename)
        elif form is 'npz':
            self._savez(filename)
        else:
            raise ValueError('form argument must be either txt or npz')

    def _savetxt(self, filename):
        # DEPRECATED: remove in 1.0.0
        bins = np.linspace(0.0, max(self.dimensions), num=self.nbins)

        # Create list of results which will be output
        output = [bins]

        for dim in ['x', 'y', 'z']:
            output.append(self.results[dim]['pos'])
            output.append(self.results[dim]['pos_std'])

        for dim in ['x', 'y', 'z']:
            output.append(self.results[dim]['char'])
            output.append(self.results[dim]['char_std'])

        density = self.totalmass / self.volume
        header = ("1 coord [Ang] 2-7 mass density (x,sx,y,sz,z,sz) [g/cm^3]"
                  "8-13 charge density (x,sx,y,sz,z,sz) [e/A^3]\n Average "
                  "density: {} g/cm3".format(density))
        np.savetxt(filename,
                   np.column_stack(output),
                   fmt='%10.5f',
                   header=header)

    def _savez(self, filename):
        # DEPRECATED: remove in 1.0.0
        bins = np.linspace(0.0, max(self.dimensions), num=self.nbins)
        dictionary = {'bins': bins}

        for dim in self.results:
            self.results[dim].pop('dim')
            self.results[dim].pop('slice volume')
            for key in self.results[dim]:
                dictionary[dim + "_" + key] = self.results[dim][key]

        np.savez(filename, **dictionary)

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
