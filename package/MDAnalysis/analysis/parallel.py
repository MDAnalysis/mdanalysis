from collections import UserDict
import numpy as np
from importlib import import_module
import warnings
from typing import Callable, Iterable, Sequence

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

    def __init__(self, backend: str, n_workers: int, client: object = None):
        """Generates an instance of a ParallelExecutor class, setting up correct
        computation parameters upon initialization.

        Parameters
        ----------
        backend : str
            type of backend (one of `self.available_backends`)
        n_workers : int
            number of workers to split the workload between
        client : object, optional
            pre-configured client for running computations, by default None

        Raises
        ------
        ValueError
            Configuration is correct if n_workers > 0 and integer, 'backend' is in `available_backends`
            and respective module is installed, or there is a pre-configured `client`
            which is an instance of `dask.distributed.Client`. Any other configuration gets error.
        """
        errors = {
            f"Should have only one of 'client' and 'backend' as non-None value, got {client=} and {backend=}": sum(
                (val is None for val in (client, backend))
            )
            != 1,
            f"'client' must be dask.distributed.Client instance, got {type(client)=}": client is not None
            and not is_valid_distributed_client(client),
            "if using dask.distributed backend, should set 'backend=None' and provide client argument": client is None
            and backend == "dask.distributed",
            f"'backend' should be one of the available backends {self.available_backends=}, got {backend=}": backend
            is not None
            and backend not in self.available_backends,
            f"'n_workers' must be a positive integer, got {n_workers=}": not (
                isinstance(n_workers, int) and n_workers > 0
            ),
            f"backend {backend} is not installed, please run 'python3 -m pip install {backend} to fix this": backend
            not in ["local", None]
            and not is_installed(backend),
        }

        for message, failing_condition in errors.items():
            if failing_condition:
                raise ValueError(message)

        warns = {
            f"Can be only 1 worker with backend='local', got {n_workers=}": backend == "local" and n_workers > 1,
        }
        for message, condition in warns.items():
            if condition:
                warnings.warn(message)

        self.backend = backend
        self.client = client
        self.n_workers = n_workers

    def apply(self, func: Callable, computations: list) -> list:
        """Applies function to the list of computations using pre-configured backend or client parameters

        Parameters
        ----------
        func : Callable
            function to be applied
        computations : list
            list of computations (single argument only) to apply the function to

        Returns
        -------
        list
            list of results of each computation

        Raises
        ------
        ValueError
            if none of the conditions for applying each computing method is met
        """
        options = {
            self._compute_with_client: self.client is not None,
            self._compute_with_dask: self.backend == "dask",
            self._compute_with_multiprocessing: self.backend == "multiprocessing",
            self._compute_with_local: self.backend == "local",
        }
        for map_, condition in options.items():
            if condition:
                return map_(func, computations)
        raise ValueError(f"Backend or client is not set properly: {self.backend=}, {self.client=}")

    def _compute_with_local(self, func: Callable, computations: list) -> list:
        """Performs all computations within current process, without using any parallel approach.

        Parameters
        ----------
        func : Callable
            function to be executed
        computations : list
            list of arguments to the function (without unpacking)

        Returns
        -------
        list
            list of results
        """
        return [func(task) for task in computations]

    def _compute_with_multiprocessing(self, func: Callable, computations: list) -> list:
        """Performs all computations using multiprocessing.Pool parallelization.

        Parameters
        ----------
        func : Callable
            function to be executed
        computations : list
            list of arguments to the function (without unpacking)

        Returns
        -------
        list
            list of results
        """
        from multiprocessing import Pool

        with Pool(processes=self.n_workers) as pool:
            results = pool.map(func, computations)
        return results

    def _compute_with_dask(self, func: Callable, computations: list) -> list:
        """Performs all computations using dask.compute(), which under the hood uses multiprocessing.Pool
        but a different method for serialization (instead of default `pickle`), hence allowing for
        potentially more functions to be executed in parallel.

        Parameters
        ----------
        func : Callable
            function to be executed
        computations : list
            list of arguments to the function (without unpacking)

        Returns
        -------
        list
            list of results
        """
        from dask.delayed import delayed
        import dask

        computations = [delayed(func)(task) for task in computations]
        results = dask.compute(computations, scheduler="processes", chunksize=1, n_workers=self.n_workers)[0]
        return results

    def _compute_with_client(self, func: Callable, computations: list) -> list:
        """Performs all computations using dask.distributed.Client.compute(),
        which sends data to potentially remote workers communicating via tcp,
        and uses the same serialization method as `dask`.

        Parameters
        ----------
        func : Callable
            function to be executed
        computations : list
            list of arguments to the function (without unpacking)

        Returns
        -------
        list
            list of results
        """
        from dask.delayed import delayed

        computations = [delayed(func)(task) for task in computations]
        return [obj.result() for obj in self.client.compute(computations)]


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
