# -*- coding: utf-8 -*-
#     pylint: disable=E1101
from __future__ import division
import os.path as path
from MDAnalysis.analysis.base import AnalysisBase
import numpy as np

class LinearDensity(AnalysisBase):
    def __init__(self, universe, selection, description="", grouping='atoms', binsize=0.25,
                 start=None, stop=None, step=None):
        self._ags = [selection] # allows use of run(parallel=True)
        self._setup_frames(universe, start, stop, step)

        self.binsize = binsize
        
        # group of atoms on which to compute the COM (same as used in
        # AtomGroup.wrap())
        self.grouping = grouping

        # Take root of trajectory filename for output file naming
        self.trajname = path.splitext(path.basename(universe.trajectory.dcdfilename))[0]
        # additional string for naming the output file
        self.description = description+"_"+str(grouping)

        # Dictionary containing results
        self.results = {'x': {'dim': 0}, 'y': {'dim': 1}, 'z': {'dim': 2}}
        # Box sides
        self.dimensions = [side for side in self._universe.dimensions[0:3]]
        self.volume = np.prod(self.dimensions)
        bins = [int(side//self.binsize) for side in self.dimensions] # number of bins

        # Here we choose a number of bins of the largest cell side so that
        # x, y and z values can use the same "coord" column in the output file
        self.nbins = max(bins)
        slices_vol = [self.volume/n for n in bins]

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
        try: # in case it's not an atom
            self.masses = [elem.total_mass() for elem in group]
            self.charges = [elem.total_charge() for elem in group]
        except AttributeError: # much much faster for atoms
            self.masses = self._ags[0].masses
            self.charges = self._ags[0].charges

        self.totalmass = np.sum(self.masses)

    #@profile
    def _single_frame(self, timestep):
        self.group = getattr(self._ags[0], self.grouping)
        self._ags[0].wrap(compound=self.grouping)

        # Find position of atom/group of atoms
        if self.grouping == 'atoms':
            positions = self._ags[0].positions # faster for atoms
        else:
            positions = np.array([elem.centroid() for elem in self.group]) # COM for res/frag/etc

        for dim in ['x', 'y', 'z']:
            idx = self.results[dim]['dim']

            key = 'pos'
            key_std = 'pos_std'
            # histogram for positions weighted on masses
            hist, _ = np.histogram(positions[:, idx], weights=self.masses,
                                   bins=self.nbins, range=(0.0, max(self.dimensions)))

            self.results[dim][key] = np.sum([self.results[dim][key], hist], axis=0)
            self.results[dim][key_std] = np.sum([self.results[dim][key_std],
                                                 np.square(hist)], axis=0)

            key = 'char'
            key_std = 'char_std'
            # histogram for positions weighted on charges
            hist, _ = np.histogram(positions[:, idx], weights=self.charges,
                                   bins=self.nbins, range=(0.0, max(self.dimensions)))

            self.results[dim][key] = np.sum([self.results[dim][key], hist], axis=0)
            self.results[dim][key_std] = np.sum([self.results[dim][key_std], np.square(hist)],
                                                axis=0)

    def _conclude(self):
        k = 6.022e-1 # divide by avodagro and convert from A3 to cm3
        bins = np.linspace(0.0, max(self.dimensions), num=self.nbins)

        # Average results over the  number of configurations
        for dim in ['x', 'y', 'z']:
            for key in ['pos', 'pos_std', 'char', 'char_std']:
                self.results[dim][key] /= self.nframes
            # Computed standard deviation for the error
            self.results[dim]['pos_std'] = np.sqrt(self.results[dim]['pos_std']
                                                   - np.square(self.results[dim]['pos']))
            self.results[dim]['char_std'] = np.sqrt(self.results[dim]['char_std']
                                                    - np.square(self.results[dim]['char']))


        # Create list of results which will be output
        output = [bins]

        for dim in ['x', 'y', 'z']:
            output.append(self.results[dim]['pos']/(self.results[dim]['slice volume']*k))
            output.append(self.results[dim]['pos_std']/(self.results[dim]['slice volume']*k))

        for dim in ['x', 'y', 'z']:
            output.append(self.results[dim]['char']/(self.results[dim]['slice volume']*k))
            output.append(self.results[dim]['char_std']/(self.results[dim]['slice volume']*k))

        # Define filename and write to output
        filename = self.trajname+"."+self.description+".ldens"
        density = self.totalmass/self.volume
        header = "1 coord [Ang] 2-7 mass density (x,sx,y,sz,z,sz) [g/cm^3]" + \
            "8-13 charge density (x,sx,y,sz,z,sz) [e/A^3]\n Average density: "\
            +str(density)+" g/cm3"
        np.savetxt(filename, np.column_stack(output), header=header)

    def _add_other_results(self, other):
        # For parallel analysis
        results = self.results
        for dim in ['x', 'y', 'z']:
            key = 'pos'
            key_std = 'pos_std'
            results[dim][key] = np.sum([results[dim][key], other[dim][key]], axis=0)
            results[dim][key_std] = np.sum([results[dim][key_std], other[dim][key_std]], axis=0)

            key = 'char'
            key_std = 'char_std'
            results[dim][key] = np.sum([results[dim][key], other[dim][key]], axis=0)
            results[dim][key_std] = np.sum([results[dim][key_std], other[dim][key_std]], axis=0)
