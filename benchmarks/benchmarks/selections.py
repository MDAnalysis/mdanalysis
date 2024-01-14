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

class GeoSelectionBench(object):
    """Benchmarks for the various MDAnalysis
    geometric selection strings.
    """

    # all selection strings verified
    # to produce non-zero atom groups
    # with GRO test file
    params = (['around 5.0 resid 1',
               'sphlayer 2.4 6.0 (protein)',
               'sphzone 6.0 (protein)',
               'cylayer 5 10 10 -8 protein',
               'cyzone 15 4 -8 protein',
               'point 5.0 5.0 5.0 3.5',
               'prop z >= 5.0',
               'prop abs z <= 5.0'],
              [True, False], # updating flags
              [[False, True], [True, False]]) # periodic flags


    # benchmarks should include static &
    # dynamic selections & periodic
    # vs non-periodic
    param_names = ['selection_string',
                   'dynamic_selection',
                   'periodic_selection']


    def setup(self,
              selection_string,
              dynamic_selection,
              periodic_selection):
        self.u = MDAnalysis.Universe(GRO)

    def time_geometric_selections(self,
                                  selection_string,
                                  dynamic_selection,
                                  periodic_selection):

        # TODO: Do we need a kwarg similar to old `use_KDTree_routines`
        # flag? We used to benchmark that.
        self.u.select_atoms(selection_string,
                            updating=dynamic_selection,
                            periodic=periodic_selection[0],
                            )
