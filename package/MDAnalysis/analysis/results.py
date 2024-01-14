"""Analysis results and their aggregation --- :mod:`MDAnalysis.analysis.results`
==============================================================
"""
from collections import UserDict
import numpy as np
from typing import Callable, Sequence


class Results(UserDict):
    r"""Container object for storing results.

    :class:`Results` are dictionaries that provide two ways by which values
    can be accessed: by dictionary key ``results["value_key"]`` or by object
    attribute, ``results.value_key``. :class:`Results` stores all results
    obtained from an analysis after calling :meth:`~AnalysisBase.run()`.

    The implementation is similar to the :class:`sklearn.utils.Bunch`
    class in `scikit-learn`_.

    .. _`scikit-learn`: https://scikit-learn.org/
    .. _`sklearn.utils.Bunch`: https://scikit-learn.org/stable/modules/generated/sklearn.utils.Bunch.html

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

    .. versionchanged:: 2.8.0
        Moved :class:`Results` to `MDAnalysis.analysis.results`
    """

    def _validate_key(self, key):
        if key in dir(self):
            raise AttributeError(f"'{key}' is a protected dictionary attribute")
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
            raise AttributeError(f"'Results' object has no attribute '{attr}'") from err

    def __delattr__(self, attr):
        try:
            del self[attr]
        except KeyError as err:
            raise AttributeError(f"'Results' object has no attribute '{attr}'") from err

    def __getstate__(self):
        return self.data

    def __setstate__(self, state):
        self.data = state


class ResultsGroup:
    """
    Management and aggregation of results stored in :class:`Results` instances.

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
    >>> from MDAnalysis.analysis.results import ResultsGroup, Results
    >>> group = ResultsGroup(lookup={'mass': ResultsGroup.float_mean})
    >>> obj1 = Results(mass=1)
    >>> obj2 = Results(mass=3)
    >>> group.merge([obj1, obj2])
    {'mass': 2.0}

    >>> # you can also set `None` for those attributes that you want to skip
    >>> lookup = {'mass': ResultsGroup.float_mean, 'trajectory': None}
    >>> group = ResultsGroup(lookup)
    >>> objects = [Results(mass=1, skip=None), Results(mass=3, skip=object)]
    >>> group.merge(objects, require_all_aggregators=False)
    {'mass': 2.0}

    .. versionadded:: 2.8.0
    """

    def __init__(self, lookup: dict[str, Callable] = None):
        self._lookup = lookup

    def merge(self, objects: Sequence[Results], require_all_aggregators: bool = True) -> Results:
        """Merge results into a single object. If objects contain single element, returns it while ignoring _lookup attribute.

        Parameters
        ----------
        require_all_aggregators : bool, optional
            if True, raise an exception when no aggregation function for a particular argument is found.
            Allows to skip aggregation for the parameters that aren't needed in the final object -- see :class:`ResultsGroup`.

        Returns
        -------
        Results
            merged :class:`Results`

        Raises
        ------
        ValueError
            if no aggregation function for a key is found and `require_all_aggregators=True`

        .. versionadded:: 2.8.0
        """
        if len(objects) == 1:
            merged_results = objects[0]
        else:
            merged_results = Results()
            for key in objects[0].keys():
                agg_function = self._lookup.get(key, None)
                if agg_function is not None:
                    results_of_t = [obj[key] for obj in objects]
                    merged_results[key] = agg_function(results_of_t)
                else:
                    if require_all_aggregators:
                        raise ValueError(f"No aggregation function for {key=}")
        return merged_results

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

    @staticmethod
    def ndarray_hstack(arrs: list[np.ndarray]):
        """Performs horizontal stack of input arrays

        Parameters
        ----------
        arrs : list[np.ndarray]
            input numpy arrays

        Returns
        -------
        np.ndarray
            result of stacking
        """
        return np.hstack(arrs)

    @staticmethod
    def ndarray_vstack(arrs: list[np.ndarray]):
        """Performs vertical stack of input arrays

        Parameters
        ----------
        arrs : list[np.ndarray]
            input numpy arrays

        Returns
        -------
        np.ndarray
            result of stacking
        """
        return np.vstack(arrs)
