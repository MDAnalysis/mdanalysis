from .results import Results, ResultsGroup
import dask.array as da
import numpy as np


class DaskTimeSeriesAnalysisBase:

    def __init__(self, dask_timeseries, verbose=False, **kwargs):
        self._dts = dask_timeseries
        self._verbose = verbose
        self.results = Results()

    def _prepare(self):
        pass  # pylint: disable=unnecessary-pass

    def _compute(self):
        pass

    def _conclude(self):
        pass  # pylint: disable=unnecessary-pass

    def run(self, verbose):
        self._prepare()
        self._compute()
        self._conclude()
        return self


class DaskRMSF(DaskTimeSeriesAnalysisBase):
    def __init__(self, dask_timeseries, verbose=False, **kwargs):
        super().__init__(dask_timeseries, verbose=verbose, **kwargs)
        self._kwargs = kwargs

    def _prepare(self):
        n_atoms = len(self._dts[0])
        self.results["rmsf"] = np.zeros((n_atoms, 3))

    def _compute(self):
        positions = self._dts
        mean_positions = positions.mean(axis=0)
        subtracted_positions = positions - mean_positions
        squared_deviations = subtracted_positions**2
        avg_squared_deviations = squared_deviations.mean(axis=0)
        sqrt_avg_squared_deviations = da.sqrt(avg_squared_deviations)
        self.results.rmsf = da.sqrt(
            (sqrt_avg_squared_deviations**2).sum(axis=1)
        )

    def _conclude(self):
        pass
