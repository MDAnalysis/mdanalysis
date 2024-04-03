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
            self.A = self.u.atoms.positions.copy()[:num_atoms]
            self.u.trajectory[-1]
            self.B = self.u.atoms.positions.copy()[:num_atoms]

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

class RmsdTrajBench(object):
    """Benchmarks for MDAnalysis.analysis.rms.RMSD
    """

    # TODO: RMSD has many parameters / options,
    # some of which are apparently still considered
    # experimental -- we'll eventually want to
    # benchmark more of these

    params = (['all', 'backbone'],
              [None, 'mass'])

    param_names = ['select',
                   'weights']

    def setup(self, select, weights):
        self.u = MDAnalysis.Universe(PSF, DCD)
        self.RMSD_inst = rms.RMSD(atomgroup=self.u,
                                  reference=None,
                                  select=select,
                                  weights=weights)

    def time_RMSD(self, select, weights):
        """Benchmark RMSD.run() method, which parses
        over the entire trajectory.
        """
        self.RMSD_inst.run()


class RmsfTrajBench(object):
    """Benchmarks for MDAnalysis.analysis.rms.RMSF
    """

    params = ([100,500,2000],
              [None, 3],
              [None,'mass'])

    param_names = ['n_atoms',
                   'step',
                   'weights']

    def setup(self, n_atoms, step, weights):
        self.u = MDAnalysis.Universe(PSF, DCD)
        self.ag = self.u.atoms[:n_atoms]
        self.RMSF_inst = rms.RMSF(atomgroup=self.ag,
                                  weights=weights)

    def time_RMSF(self, n_atoms, step, weights):
        """Benchmark RMSF.run() method, which parses
        over the entire trajectory.
        """
        self.RMSF_inst.run(step=step)
