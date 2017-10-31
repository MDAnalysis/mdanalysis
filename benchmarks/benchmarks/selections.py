from __future__ import division, absolute_import, print_function

import MDAnalysis

try:
    from MDAnalysisTests.datafiles import GRO
except:
    pass

class SimpleSelectionBench(object):
    """Benchmarks for the various MDAnalysis
    simple selection strings.
    """
    params = ('protein',
              'backbone',
              'nucleic',
              'nucleicbackbone',
              'resid 1:10',
              'resnum 1:10',
              'resname LYS',
              'name CA',
              'bynum 0:10')

    param_names = ['selection_string']

    def setup(self, selection_string):
        self.u = MDAnalysis.Universe(GRO)

    def time_simple_selections(self, selection_string):
        """Benchmark simple selections on the protein-based
        standard test GRO file.
        """
        if hasattr(MDAnalysis.Universe, 'select_atoms'):
            self.u.select_atoms(selection_string)
        else:
            self.u.selectAtoms(selection_string)
