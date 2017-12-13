from __future__ import division, absolute_import, print_function

import MDAnalysis

try:
    from MDAnalysisTests.datafiles import PSF, DCD
except:
    pass

try:
    from MDAnalysis.analysis import rms
except:
    pass

class SimpleRmsBench(object):
    """Benchmarks for MDAnalysis.analysis.rms.rmsd
    """

    params = ([100, 500, 2000],
              [None, [1.0, 0.5]],
              [False, True],
              [False, True])
    param_names = ['num_atoms',
                   'weights',
                   'center',
                   'superposition']

    def setup(self, num_atoms, weights, center, superposition):
        # mimic rmsd docstring example code
        self.u = MDAnalysis.Universe(PSF, DCD)
        # ag.positions is the new syntax
        # but older commit hashes will need to use
        # ag.coordinates()
        try:
            self.A = self.u.atoms.positions.copy()[:num_atoms]
            self.u.trajectory[-1]
            self.B = self.u.atoms.positions.copy()[:num_atoms]
        except:
            self.A = self.u.atoms.coordinates().copy()[:num_atoms]
            self.u.trajectory[-1]
            self.B = self.u.atoms.coordinates().copy()[:num_atoms]

    def time_rmsd(self, num_atoms, weights, center, superposition):
        """Benchmark rmsd function using a setup similar to
        its docstring example code along with several possible
        permutations of parameters.
        """
        rms.rmsd(a=self.A,
                 b=self.B,
                 weights=weights,
                 center=center,
                 superposition=superposition)


