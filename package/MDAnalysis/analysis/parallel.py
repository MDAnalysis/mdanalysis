from collections import UserDict
import numpy as np
from importlib import import_module
import warnings
from typing import Callable, Iterable, Sequence, Union

from MDAnalysis.lib.util import is_installed


def is_valid_distributed_client(client: object):
    try:
        import dask.distributed as distributed
    except ImportError:
        return False
    return isinstance(client, distributed.Client)


class Results(UserDict):
    r"""Container object for storing results.

    :class:`Results` are dictionaries that provide two ways by which values
    can be accessed: by dictionary key ``results["value_key"]`` or by object
    attribute, ``results.value_key``. :class:`Results` stores all results
    obtained from an analysis after calling :meth:`~AnalysisBase.run()`.

    The implementation is similar to the :class:`sklearn.utils.Bunch`
    class in `scikit-learn`_.

    .. _`scikit-learn`: https://scikit-learn.org/

    Raises
    ------
    AttributeError
        If an assigned attribute has the same name as a default attribute.

    ValueError
        If a key is not of type ``str`` and therefore is not able to be
        accessed by attribute.

    Examples
    --------
    >>> from MDAnalysis.analysis.base import Results
    >>> results = Results(a=1, b=2)
    >>> results['b']
    2
    >>> results.b
    2
    >>> results.a = 3
    >>> results['a']
    3
    >>> results.c = [1, 2, 3, 4]
    >>> results['c']
    [1, 2, 3, 4]


    .. versionadded:: 2.0.0
    """

    def _validate_key(self, key):
        if key in dir(self):
            raise AttributeError(f"'{key}' is a protected dictionary " "attribute")
        elif isinstance(key, str) and not key.isidentifier():
            raise ValueError(f"'{key}' is not a valid attribute")

    def __init__(self, *args, **kwargs):
        kwargs = dict(*args, **kwargs)
        if "data" in kwargs.keys():
            raise AttributeError(f"'data' is a protected dictionary attribute")
        self.__dict__["data"] = {}
        self.update(kwargs)

    def __setitem__(self, key, item):
        self._validate_key(key)
        super().__setitem__(key, item)

    def __setattr__(self, attr, val):
        if attr == "data":
            super().__setattr__(attr, val)
        else:
            self.__setitem__(attr, val)

    def __getattr__(self, attr):
        try:
            return self[attr]
        except KeyError as err:
            raise AttributeError("'Results' object has no " f"attribute '{attr}'") from err

    def __delattr__(self, attr):
        try:
            del self[attr]
        except KeyError as err:
            raise AttributeError("'Results' object has no " f"attribute '{attr}'") from err

    def __getstate__(self):
        return self.data

    def __setstate__(self, state):
        self.data = state


class BackendBase:
    def __init__(self, n_workers: int):
        self.n_workers = n_workers
        self._validate()

    def _get_checks(self):
        return {
            isinstance(self.n_workers, int)
            and self.n_workers > 0: f"n_workers should be positive integer, got {self.n_workers=}",
        }

    def _get_warnings(self):
        return dict()

    def _validate(self):
        for check, msg in self._get_checks().items():
            if not check:
                raise ValueError(msg)
        for check, msg in self._get_warnings().items():
            if not check:
                warnings.warn(msg)

    def apply(self, func: Callable, computations: list) -> list:
        raise NotImplementedError("Should be re-implemented in subclasses")


class BackendSerial(BackendBase):
    def _get_warnigns(self):
        return {self.n_workers > 1, "n_workers > 1 will be ignored while executing with backend='serial'"}

    def apply(self, func: Callable, computations: list) -> list:
        return [func(task) for task in computations]


class BackendMultiprocessing(BackendBase):
    def apply(self, func: Callable, computations: list) -> list:
        from multiprocessing import Pool

        with Pool(processes=self.n_workers) as pool:
            results = pool.map(func, computations)
        return results


class BackendDask(BackendBase):
    def apply(self, func: Callable, computations: list) -> list:
        from dask.delayed import delayed
        import dask

        computations = [delayed(func)(task) for task in computations]
        results = dask.compute(computations, scheduler="processes", chunksize=1, n_workers=self.n_workers)[0]
        return results

    def _get_checks(self):
        base_checks = super()._get_checks()
        checks = {is_installed("dask"): "module 'dask' should be installed: run 'python3 -m pip install dask"}
        return base_checks | checks


class BackendDaskDistributed(BackendBase):
    def __init__(self, client, n_workers: int = None):
        self.client = client
        self.n_workers = n_workers
        self._validate()
        self.n_workers = sum(client.nthreads().values())

    def _get_checks(self):
        from dask.distributed import Client

        checks = {
            is_installed("dask.distributed"): "module 'dask' should be installed: run 'python3 -m pip install dask",
            isinstance(
                self.client, Client
            ): f"self.client should be instance of dask.distributed.Client, got {type(self.client)=}",
        }
        return checks

    def _get_warnings(self):
        return {
            self.n_workers
            is not None: "n_workers will be ignored with 'dask.distributed' backend and deduced from self.client directly",
        }

    def apply(self, func: Callable, computations: list) -> list:
        from dask.delayed import delayed

        computations = [delayed(func)(task) for task in computations]
        return [obj.result() for obj in self.client.compute(computations)]


class ParallelExecutor:
    r"""Executor object that abstracts away running parallel computations
    with various backends. Mainly is used for AnalysisBase, but can be used for other purposes as well."""

    @classmethod
    @property
    def available_backends(cls):
        """Declares backends that are implemented for parallel execution.
        Importantly, backend names are the same as packages you must import
        in order to run them -- e.g. to run 'dask.distributed', you need to run 'import dask.distributed'

        Returns
        -------
        list[str]
            list of available backends as strings
        """
        return ("local", "multiprocessing", "dask", "dask.distributed")


class ResultsGroup:
    """Simple class that does aggregation of the results from independent remote workers in AnalysisBase"""

    def __init__(self, lookup: dict[str, Callable] = None):
        """Initializes the class, saving aggregation functions in _lookup attribute.

        Parameters
        ----------
        lookup : dict[str, Callable], optional
            aggregation functions lookup dict, by default None
        """
        self._lookup = lookup

    def merge(self, objects: Sequence[Results], do_raise: bool = True) -> Results:
        """Merge results into a single object.
        If objects contain single element, returns it while ignoring _lookup attribute.

        Parameters
        ----------
        do_raise : bool, optional
            if you want to raise an exception when no aggregation function for a particular argument is found, by default True

        Returns
        -------
        Results
            merged Result object

        Raises
        ------
        ValueError
            if no aggregation function for a key is found and do_raise = True
        """
        if len(objects) == 1:
            rv = objects[0]
        else:
            rv = Results()
            for key in objects[0].keys():
                agg_function = self._lookup.get(key, None)
                if agg_function is None and do_raise:
                    raise ValueError(f"No aggregation function for {key=}")
                results_of_t = [obj[key] for obj in objects]
                rv[key] = agg_function(results_of_t)
        return rv

    @staticmethod
    def flatten_sequence(arrs: list[list]):
        return [item for sublist in arrs for item in sublist]

    @staticmethod
    def ndarray_sum(arrs: np.ndarray):
        return np.array(arrs).sum(axis=0)

    @staticmethod
    def ndarray_mean(arrs: np.ndarray):
        return np.array(arrs).sum(axis=0)

    @staticmethod
    def float_mean(floats: list[float]):
        return sum(floats) / len(floats)
