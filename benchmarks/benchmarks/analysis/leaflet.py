from __future__ import division, absolute_import, print_function

import MDAnalysis

# use a lipid bilayer system for leaflet testing
# since this is a common use case
try:
    from MDAnalysisTests.datafiles import Martini_membrane_gro
except:
    pass

try:
    from MDAnalysis.analysis import leaflet
except:
    pass

class LeafletBench(object):
    """Benchmarks for MDAnalysis.analysis.leaflet
    functions.
    """

    params = ([7.0, 15.0, 23.0],
              [None, True, False],
              [True, False])
    param_names = ['cutoff', 'sparse', 'pbc']

    def setup(self, cutoff, sparse, pbc):
        self.u = MDAnalysis.Universe(Martini_membrane_gro)
        self.headgroup_sel = 'name PO4'

    def time_leafletfinder(self, cutoff, sparse, pbc):
        """Benchmark LeafletFinder for test lipid
        membrane system.
        """
        leaflet.LeafletFinder(self.u,
                              self.headgroup_sel,
                              cutoff,
                              pbc,
                              sparse)

    def time_optimize_cutoff(self, cutoff, sparse, pbc):
        """Benchmark optimize_cutoff for test lipid
        membrane system using default network distance
        range.
        """
        leaflet.optimize_cutoff(self.u,
                                self.headgroup_sel,
                                pbc=pbc,
                                sparse=sparse)
