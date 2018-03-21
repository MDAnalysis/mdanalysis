from __future__ import division, absolute_import, print_function

import MDAnalysis

try:
    from MDAnalysisTests.datafiles import PSF, DCD
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
              [[0,5], [0,15], [0,20]])

    param_names = ['nbins',
                   'range_val']

    def setup(self, nbins, range_val):
        
        # use explicit OR in selections
        # for increased backward compatibility
        self.g1_sel_str = ('(resname ALA) or '
                           '(resname ARG) or '
                           '(resname ASN)')

        self.g2_sel_str = ('(resname THR) or '
                           '(resname TYR) or '
                           '(resname VAL)')

        self.u = MDAnalysis.Universe(PSF, DCD)

        try:
            self.sel1 = self.u.select_atoms(self.g1_sel_str)
            self.sel2 = self.u.select_atoms(self.g2_sel_str)
        except AttributeError:
            self.sel1 = self.u.selectAtoms(self.g1_sel_str)
            self.sel2 = self.u.selectAtoms(self.g2_sel_str)

        # do not include initialization of the
        # InterRDF object in the benchmark itself

        self.rdf = InterRDF(g1=self.sel1,
                            g2=self.sel2,
                            nbins=nbins,
                            range=range_val)

    def time_interrdf(self, nbins, range_val):
        """Benchmark a full trajectory parse
        by MDAnalysis.analysis.rdf.InterRDF
        """
        self.rdf.run()
