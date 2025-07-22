"""Analysis results and their aggregation --- :mod:`MDAnalysis.analysis.results`
================================================================================

Module introduces two classes, :class:`Results` and :class:`ResultsGroup`,
used for storing and aggregating data in
:meth:`MDAnalysis.analysis.base.AnalysisBase.run()`, respectively.


Classes
-------

The :class:`Results` class is an extension of a built-in dictionary
type, that holds all assigned attributes in :attr:`self.data` and 
allows for access either via dict-like syntax, or via class-like syntax:

.. code-block:: python

    from MDAnalysis.analysis.results import Results
    r = Results()
    r.array = [1, 2, 3, 4]
    assert r['array'] == r.array == [1, 2, 3, 4]


The :class:`ResultsGroup` can merge multiple :class:`Results` objects.
It is mainly used by :class:`MDAnalysis.analysis.base.AnalysisBase` class, 
that uses :meth:`ResultsGroup.merge()` method to aggregate results from
multiple workers, initialized during a parallel run:

.. code-block:: python

    from MDAnalysis.analysis.results import Results, ResultsGroup
    import numpy as np
    
    r1, r2 = Results(), Results()
    r1.masses = [1, 2, 3, 4, 5]
    r2.masses = [0, 0, 0, 0]
    r1.vectors = np.arange(10).reshape(5, 2)
    r2.vectors = np.arange(8).reshape(4, 2)

    group = ResultsGroup(
        lookup = {
            'masses': ResultsGroup.flatten_sequence,
            'vectors': ResultsGroup.ndarray_vstack
            }
        )

    r = group.merge([r1, r2])
    assert r.masses == list((*r1.masses, *r2.masses))
    assert (r.vectors == np.vstack([r1.vectors, r2.vectors])).all()
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
        Moved :class:`Results` to :mod:`MDAnalysis.analysis.results`
    """

    def _validate_key(self, key):
        if key in dir(self):
            raise AttributeError(
                f"'{key}' is a protected dictionary attribute"
            )
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
            raise AttributeError(
                f"'Results' object has no attribute '{attr}'"
            ) from err

    def __delattr__(self, attr):
        try:
            del self[attr]
        except KeyError as err:
            raise AttributeError(
                f"'Results' object has no attribute '{attr}'"
            ) from err

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

    .. code-block:: python

        from MDAnalysis.analysis.results import ResultsGroup, Results
        group = ResultsGroup(lookup={'mass': ResultsGroup.float_mean})
        obj1 = Results(mass=1)
        obj2 = Results(mass=3)
        assert {'mass': 2.0} == group.merge([obj1, obj2])


    .. code-block:: python

        # you can also set `None` for those attributes that you want to skip
        lookup = {'mass': ResultsGroup.float_mean, 'trajectory': None}
        group = ResultsGroup(lookup)
        objects = [Results(mass=1, skip=None), Results(mass=3, skip=object)]
        assert group.merge(objects, require_all_aggregators=False) == {'mass': 2.0}

    .. versionadded:: 2.8.0
    """

    def __init__(self, lookup: dict[str, Callable] = None):
        self._lookup = lookup

    def merge(
        self, objects: Sequence[Results], require_all_aggregators: bool = True
    ) -> Results:
        """Merge multiple Results into a single Results instance.

        Merge multiple :class:`Results` instances into a single one, using the
        `lookup` dictionary to determine the appropriate aggregator functions for
        each named results attribute. If the resulting object only contains a single
        element, it just returns it without using any aggregators.

        Parameters
        ----------
        objects : Sequence[Results]
            Multiple :class:`Results` instances with the same data attributes.
        require_all_aggregators : bool, optional
            if True, raise an exception when no aggregation function for a
            particular argument is found. Allows to skip aggregation for the
            parameters that aren't needed in the final object --
            see :class:`ResultsGroup`.

        Returns
        -------
        Results
            merged :class:`Results`

        Raises
        ------
        ValueError
            if no aggregation function for a key is found and ``require_all_aggregators=True``
        """
        if len(objects) == 1:
            merged_results = objects[0]
            return merged_results

        merged_results = Results()
        for key in objects[0].keys():
            agg_function = self._lookup.get(key, None)
            if agg_function is not None:
                results_of_t = [obj[key] for obj in objects]
                merged_results[key] = agg_function(results_of_t)
            elif require_all_aggregators:
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
        """sums an ndarray along ``axis=0``

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
        """calculates mean of input ndarrays along ``axis=0``

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
