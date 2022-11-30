import MDAnalysis
import numpy as np

try:
    from MDAnalysis.analysis import psa
except:
    pass

class PSA_sqnormBench(object):
        """Benchmarks for MDAnalysis.analysis.psa.
        sqnorm
        """

        params = ([2,3,4],
                  [100,1000,10000],
                  [None, 0, 1, -1])

        # num_cols is equivalent to dimensions
        # num_rows is equivalent to i.e., num atoms
        param_names = ['num_cols',
                       'num_rows',
                       'axis']

        def setup(self, num_cols, num_rows, axis):
            np.random.seed(170089)
            self.v = np.random.rand(num_rows, num_cols)

        def time_sqnorm(self, num_cols, num_rows, axis):
            """Benchmark sqnorm in psa module
            """
            psa.sqnorm(v=self.v, axis=axis)

class PSA_get_msd_matrixBench(object):
        """Benchmarks for MDAnalysis.analysis.psa.
        get_msd_matrix
        """

        params = ([10,100,1000],
                  [5,25,50])

        # since the function is defined to work with
        # 3N dimension data sets, we will restrict
        # benchmarks to that dimensionality
        param_names = ['time_steps',
                       'n_atoms']

        def setup(self, time_steps, n_atoms):
            np.random.seed(170089)
            self.P = np.random.rand(time_steps,
                                    n_atoms,
                                    3)
            np.random.seed(971132)
            self.Q = np.random.rand(time_steps,
                                    n_atoms,
                                    3)

        def time_get_msd_matrix(self, time_steps, n_atoms):
            """Benchmark for get_msd_matrix in psa module
            """
            # only default argument for axis is benchmarked
            psa.get_msd_matrix(P=self.P,
                               Q=self.Q,
                               axis=None)

class PSA_get_coord_axesBench(object):
    """Benchmarks for MDAnalysis.analysis.psa.
    get_coord_axes
    """

    params = ([10,100,1000],
              [5, 25, 50])

    param_names = ['time_steps',
                   'n_atoms']

    def setup(self, time_steps, n_atoms):
        np.random.seed(170089)
        # only using condensed path input
        # data structure for now
        self.path = np.random.rand(time_steps,
                                   n_atoms * 3)

    def time_get_coord_axes(self, time_steps, n_atoms):
        """Benchmark get_coord_axes in psa module
        """
        psa.get_coord_axes(path=self.path)

class PSA_get_path_metric_funcBench(object):
    """Benchmark for MDAnalysis.analysis.psa.
    get_path_metric_func
    """

    params = (['hausdorff',
               'weighted_average_hausdorff',
               'average_hausdorff',
               'hausdorff_neighbors',
               'discrete_frechet'])

    param_names = ['path_metric']

    def time_get_path_metric_func(self, path_metric):
        """Benchmark for get_path_metric_func in psa
        module
        """
        psa.get_path_metric_func(name=path_metric)

class PSA_metricBench(object):
    """Benchmarks for the various path metric calculations
    in the psa module.
    """

    params = ([10,100,200],
              [5,25,50])

    param_names = ['time_steps',
                   'n_atoms']

    def setup(self, time_steps, n_atoms):
        np.random.seed(170089)
        self.P = np.random.rand(time_steps,
                                n_atoms,
                                3)
        np.random.seed(971132)
        self.Q = np.random.rand(time_steps,
                                n_atoms,
                                3)

    def time_hausdorff(self, time_steps, n_atoms):
        """Benchmark for hausdorff() in psa module.
        """
        psa.hausdorff(P=self.P,
                      Q=self.Q)

    def time_hausdorff_wavg(self, time_steps, n_atoms):
        """Benchmark for hausdorff_wavg() in psa module.
        """
        psa.hausdorff_wavg(P=self.P,
                           Q=self.Q)

    def time_hausdorff_avg(self, time_steps, n_atoms):
        """Benchmark for hausdorff_avg() in psa module.
        """
        psa.hausdorff_avg(P=self.P,
                          Q=self.Q)


    def time_hausdorff_neighbors(self, time_steps, n_atoms):
        """Benchmark for hausdorff_neighbors() in psa module.
        """
        psa.hausdorff_neighbors(P=self.P,
                                Q=self.Q)

    def time_discrete_frechet(self, time_steps, n_atoms):
        """Benchmark for discrete_frechet() in psa module.
        """
        psa.discrete_frechet(P=self.P,
                             Q=self.Q)
