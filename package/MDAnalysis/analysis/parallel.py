from collections import UserDict
import numpy as np
from importlib import import_module
import warnings
from typing import Callable, Iterable, Sequence, Union

from MDAnalysis.lib.util import is_installed


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
    """Base class for backend implementation. Initializes an instance and performs checks for its validity, such as n_workers and possibly other ones.

    Parameters
    ----------
    n_workers : int
        positive integer with number of workers (usually, processes) to split the work between

    Examples
    --------
    >>> # implement a thread-based backend
    >>> from MDAnalysis.analysis.parallel import BackendBase
    >>> class ThreadsBackend(BackendBase):
            def apply(self, func, computations):
                from multiprocessing.dummy import Pool

                with Pool(processes=self.n_workers) as pool:
                    results = pool.map(func, computations)
                return results
    >>> from MDAnalysis.analysis.rms import RMSD
    >>> R = RMSD(...) # setup the run
    >>> n_workers = 2
    >>> backend = ThreadsBackend(n_workers=n_workers)
    >>> R.run(backend=backend)

    .. versionadded:: 2.7.0
    """

    def __init__(self, n_workers: int):
        self.n_workers = n_workers
        self._validate()

    def _get_checks(self):
        """Get dictionary with `condition: error_message` pairs that ensure the validity of the backend instance

        Returns
        -------
        dict
            dictionary with `condition: error_message` pairs that will get checked during _validate() run

        .. versionadded: 2.7.0
        """
        return {
            isinstance(self.n_workers, int)
            and self.n_workers > 0: f"n_workers should be positive integer, got {self.n_workers=}",
        }

    def _get_warnings(self):
        """Get dictionary with `condition: warning_message` pairs that ensure the good usage of the backend instance

        Returns
        -------
        dict
            dictionary with `condition: warning_message` pairs that will get checked during _validate() run

        .. versionadded: 2.7.0
        """
        return dict()

    def _validate(self):
        """Check correctness (e.g. `dask` is installed if using `backend='dask'`)
        and good usage (e.g. `n_workers=1` if backend is serial) of the backend

        Raises
        ------
        ValueError
            if one of the conditions in :meth:`self._get_checks()` is True

        .. versionadded: 2.7.0
        """
        for check, msg in self._get_checks().items():
            if not check:
                raise ValueError(msg)
        for check, msg in self._get_warnings().items():
            if not check:
                warnings.warn(msg)

    def apply(self, func: Callable, computations: list) -> list:
        """Main function that will get called when using an instance of an object, mapping function to all tasks
        in the `computations` list. Should effectively be equivalent to running [func(item) for item in computations]
        while using the parallel backend capabilities.

        Parameters
        ----------
        func : Callable
            function to be called on each of the tasks in computations list
        computations : list
            computation tasks to apply function to

        Returns
        -------
        list
            list of results of the function

        .. versionadded: 2.7.0
        """
        raise NotImplementedError("Should be re-implemented in subclasses")


class BackendSerial(BackendBase):
    """A built-in backend that does serial execution of the function, without any parallelization

    .. versionadded: 2.7.0
    """

    def _get_warnigns(self):
        return {self.n_workers > 1, "n_workers > 1 will be ignored while executing with backend='serial'"}

    def apply(self, func: Callable, computations: list) -> list:
        return [func(task) for task in computations]


class BackendMultiprocessing(BackendBase):
    """A built-in backend that executes a given function using multiprocessing.Pool.map method

    .. versionadded: 2.7.0
    """

    def apply(self, func: Callable, computations: list) -> list:
        from multiprocessing import Pool

        with Pool(processes=self.n_workers) as pool:
            results = pool.map(func, computations)
        return results


class BackendDask(BackendBase):
    """A built-in backend that executes a given function using dask.delayed.compute method with `scheduler='processes'`
    and `chunksize=1`. Requires `dask` module to be installed.

    .. versionadded: 2.7.0
    """

    def apply(self, func: Callable, computations: list) -> list:
        from dask.delayed import delayed
        import dask

        computations = [delayed(func)(task) for task in computations]
        results = dask.compute(computations, scheduler="processes", chunksize=1, n_workers=self.n_workers)[0]
        return results

    def _get_checks(self):
        base_checks = super()._get_checks()
        checks = {is_installed("dask"): "module 'dask' should be installed: run 'python3 -m pip install dask'"}
        return base_checks | checks


class ResultsGroup:
    """
    Management and aggregation of results stored in :class:`Result` instances.

    A :class:`ResultsGroup` is an optional description for :class:`Result` "dictionaries"
    that are used in analysis classes based on :class:`AnalysisBase`. For each *key* in a
    :class:`Result` it describes how multiple pieces of the data held under the key are
    to be aggregated. This approach is necessary when parts of a trajectory are analyzed
    independently (e.g., in parallel) and then need to me merged (with :meth:`merge`) to
    obtain a complete data set.

    Parameters
    ----------
    lookup : dict[str, Callable], optional
        aggregation functions lookup dict, by default None

    Examples
    --------
    >>> from MDAnalysis.analysis.parallel import ResultsGroup, Results
    >>> group = ResultsGroup(lookup={'mass': ResultsGroup.float_mean})
    >>> obj1 = Results(mass=1)
    >>> obj2 = Results(mass=3)
    >>> group.merge([obj1, obj2])
    {'mass': 2.0}

    >>> # you can also set `lookup[attribute]=None` to those attributes that you want to skip
    >>> lookup = {'mass': ResultsGroup.float_mean, 'trajectory': None}
    >>> group = ResultsGroup(lookup)
    >>> objects = [Results(mass=1, skip=None), Results(mass=3, skip=object)]
    >>> group.merge(objects, require_all_aggregators=False)
    {'mass': 2.0}

    .. versionadded: 2.7.0
    """

    def __init__(self, lookup: dict[str, Callable] = None):
        self._lookup = lookup

    def merge(self, objects: Sequence[Results], require_all_aggregators: bool = True) -> Results:
        """Merge results into a single object. If objects contain single element, returns it while ignoring _lookup attribute.

        Parameters
        ----------
        require_all_aggregators : bool, optional
            if you want to raise an exception when no aggregation function for a particular argument is found, by default True.
            Allows to skip aggregation for the parameters that aren't needed in the final object:

        Returns
        -------
        Results
            merged Results object

        Raises
        ------
        ValueError
            if no aggregation function for a key is found and `require_all_aggregators=True`

        .. versionadded: 2.7.0
        """
        if len(objects) == 1:
            rv = objects[0]
        else:
            rv = Results()
            for key in objects[0].keys():
                agg_function = self._lookup.get(key, None)
                if agg_function is not None:
                    results_of_t = [obj[key] for obj in objects]
                    rv[key] = agg_function(results_of_t)
                else:
                    if require_all_aggregators:
                        raise ValueError(f"No aggregation function for {key=}")
        return rv

    @staticmethod
    def flatten_sequence(arrs: list[list]):
        """Flatten a list of lists into a list

        Parameters
        ----------
        arrs : list[list]
            list of lists

        Returns
        -------
        list
            flattened list
        """
        return [item for sublist in arrs for item in sublist]

    @staticmethod
    def ndarray_sum(arrs: list[np.ndarray]):
        """sums an ndarray along `axis=0`

        Parameters
        ----------
        arrs : list[np.ndarray]
            list of input arrays. Must have the same shape.

        Returns
        -------
        np.ndarray
            sum of input arrays
        """
        return np.array(arrs).sum(axis=0)

    @staticmethod
    def ndarray_mean(arrs: list[np.ndarray]):
        """calculates mean of input ndarrays along `axis=0`

        Parameters
        ----------
        arrs : list[np.ndarray]
            list of input arrays. Must have the same shape.

        Returns
        -------
        np.ndarray
            mean of input arrays
        """
        return np.array(arrs).mean(axis=0)

    @staticmethod
    def float_mean(floats: list[float]):
        """calculates mean of input float values

        Parameters
        ----------
        floats : list[float]
            list of float values

        Returns
        -------
        float
            mean value
        """
        return np.array(floats).mean()
