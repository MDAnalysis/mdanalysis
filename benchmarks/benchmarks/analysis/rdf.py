import MDAnalysis

try:
    from MDAnalysisTests.datafiles import TPR, XTC
except:
    pass

try:
    from MDAnalysis.analysis.rdf import InterRDF
except:
    pass

class SimpleRdfBench(object):
    """Benchmarks for MDAnalysis.analysis.rdf
    """

    params = ([20,75,200],
              [[0,5], [0,15], [0,20]],
              [1, 100, 1000, 10000])

    param_names = ['nbins',
                   'range_val',
                   'natoms']

    def setup(self, nbins, range_val, natoms):
        
        self.sel_str = 'name OW'

        self.u = MDAnalysis.Universe(TPR, XTC)

        try:
            self.sel = self.u.select_atoms(self.sel_str)[:natoms]
        except AttributeError:
            self.sel = self.u.selectAtoms(self.sel_str)[:natoms]

        # do not include initialization of the
        # InterRDF object in the benchmark itself

        self.rdf = InterRDF(g1=self.sel,
                            g2=self.sel,
                            nbins=nbins,
                            range=range_val)

    def time_interrdf(self, nbins, range_val, natoms):
        """Benchmark a full trajectory parse
        by MDAnalysis.analysis.rdf.InterRDF
        """
        self.rdf.run()
