"""Analysis backends --- :mod:`MDAnalysis.analysis.backends`
============================================================

.. versionadded:: 2.8.0


The :mod:`backends` module provides :class:`BackendBase` base class to
implement custom backends for
:meth:`MDAnalysis.analysis.base.AnalysisBase.run` and its
subclasses.

.. SeeAlso:: :ref:`parallel-analysis`

.. _backends:

Backends
--------

Three built-in backend classes classes are provided:

* *serial*: :class:`BackendSerial`, that is equal to using no
  parallelization and is default

* *multiprocessing*: :class:`BackendMultiprocessing` that supports
  parallelization via standard Python :mod:`multiprocessing` module
  and uses default :mod:`pickle` serialization

* *dask*: :class:`BackendDask`, that uses the same process-based
  parallelization as :class:`BackendMultiprocessing`, but different
  serialization algorithm via `dask <https://dask.org/>`_ (see `dask
  serialization algorithms
  <https://distributed.dask.org/en/latest/serialization.html>`_ for details)

Classes
-------

"""
import warnings
from typing import Callable
from MDAnalysis.lib.util import is_installed


class BackendBase:
    """Base class for backend implementation.

    Initializes an instance and performs checks for its validity, such as
    ``n_workers`` and possibly other ones.

    Parameters
    ----------
    n_workers : int
        number of workers (usually, processes) over which the work is split

    Examples
    --------
    .. code-block:: python

        from MDAnalysis.analysis.backends import BackendBase

        class ThreadsBackend(BackendBase):
            def apply(self, func, computations):
                from multiprocessing.dummy import Pool

                with Pool(processes=self.n_workers) as pool:
                    results = pool.map(func, computations)
                return results

        import MDAnalysis as mda
        from MDAnalysis.tests.datafiles import PSF, DCD
        from MDAnalysis.analysis.rms import RMSD

        u = mda.Universe(PSF, DCD)
        ref = mda.Universe(PSF, DCD)

        R = RMSD(u, ref)

        n_workers = 2
        backend = ThreadsBackend(n_workers=n_workers)
        R.run(backend=backend, unsupported_backend=True)

    .. warning::
        Using `ThreadsBackend` above will lead to erroneous results, since it
        is an educational example. Do not use it for real analysis.


    .. versionadded:: 2.8.0

    """

    def __init__(self, n_workers: int):
        self.n_workers = n_workers
        self._validate()

    def _get_checks(self):
        """Get dictionary with ``condition: error_message`` pairs that ensure the
        validity of the backend instance

        Returns
        -------
        dict
            dictionary with ``condition: error_message`` pairs that will get
            checked during ``_validate()`` run
        """
        return {
            isinstance(self.n_workers, int) and self.n_workers > 0:
            f"n_workers should be positive integer, got {self.n_workers=}",
        }

    def _get_warnings(self):
        """Get dictionary with ``condition: warning_message`` pairs that ensure
        the good usage of the backend instance

        Returns
        -------
        dict
            dictionary with ``condition: warning_message`` pairs that will get
            checked during ``_validate()`` run
        """
        return dict()

    def _validate(self):
        """Check correctness (e.g. ``dask`` is installed if using ``backend='dask'``)
        and good usage (e.g. ``n_workers=1`` if backend is serial) of the backend

        Raises
        ------
        ValueError
            if one of the conditions in :meth:`_get_checks` is ``True``
        """
        for check, msg in self._get_checks().items():
            if not check:
                raise ValueError(msg)
        for check, msg in self._get_warnings().items():
            if not check:
                warnings.warn(msg)

    def apply(self, func: Callable, computations: list) -> list:
        """map function `func` to all tasks in the `computations` list

        Main method that will get called when using an instance of
        ``BackendBase``.  It is equivalent to running ``[func(item) for item in
        computations]`` while using the parallel backend capabilities.

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

        """
        raise NotImplementedError


class BackendSerial(BackendBase):
    """A built-in backend that does serial execution of the function, without any
    parallelization.

    Parameters
    ----------
    n_workers : int
        Is ignored in this class, and if ``n_workers`` > 1, a warning will be
        given.


    .. versionadded:: 2.8.0
    """

    def _get_warnings(self):
        """Get dictionary with ``condition: warning_message`` pairs that ensure
        the good usage of the backend instance. Here, it checks if the number
        of workers is not 1, otherwise gives warning.

        Returns
        -------
        dict
            dictionary with ``condition: warning_message`` pairs that will get
            checked during ``_validate()`` run
        """
        return {
            self.n_workers == 1:
            "n_workers is ignored when executing with backend='serial'"
        }

    def apply(self, func: Callable, computations: list) -> list:
        """
        Serially applies `func` to each task object in ``computations``.

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
        """
        return [func(task) for task in computations]


class BackendMultiprocessing(BackendBase):
    """A built-in backend that executes a given function using the
    :meth:`multiprocessing.Pool.map <multiprocessing.pool.Pool.map>` method.

    Parameters
    ----------
    n_workers : int
        number of processes in :class:`multiprocessing.Pool
        <multiprocessing.pool.Pool>` to distribute the workload
        between. Must be a positive integer.

    Examples
    --------

    .. code-block:: python

        from MDAnalysis.analysis.backends import BackendMultiprocessing
        import multiprocessing as mp

        backend_obj = BackendMultiprocessing(n_workers=mp.cpu_count())


    .. versionadded:: 2.8.0

    """

    def apply(self, func: Callable, computations: list) -> list:
        """Applies `func` to each object in ``computations`` using `multiprocessing`'s `Pool.map`.

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
        """
        from multiprocessing import Pool

        with Pool(processes=self.n_workers) as pool:
            results = pool.map(func, computations)
        return results


class BackendDask(BackendBase):
    """A built-in backend that executes a given function with *dask*.

    Execution is performed with the :func:`dask.compute` function of
    :class:`dask.delayed.Delayed` object (created with
    :func:`dask.delayed.delayed`) with ``scheduler='processes'`` and
    ``chunksize=1`` (this ensures uniform distribution of tasks among
    processes). Requires the `dask package <https://docs.dask.org/en/stable/>`_
    to be `installed <https://docs.dask.org/en/stable/install.html>`_.

    Parameters
    ----------
    n_workers : int
        number of processes in to distribute the workload
        between. Must be a positive integer. Workers are actually
        :class:`multiprocessing.pool.Pool` processes, but they use a different and
        more flexible `serialization protocol
        <https://docs.dask.org/en/stable/phases-of-computation.html#graph-serialization>`_.

    Examples
    --------

    .. code-block:: python

        from MDAnalysis.analysis.backends import BackendDask
        import multiprocessing as mp

        backend_obj = BackendDask(n_workers=mp.cpu_count())


    .. versionadded:: 2.8.0

    """

    def apply(self, func: Callable, computations: list) -> list:
        """Applies `func` to each object in ``computations``.

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
        """
        from dask.delayed import delayed
        import dask

        computations = [delayed(func)(task) for task in computations]
        results = dask.compute(computations,
                               scheduler="processes",
                               chunksize=1,
                               num_workers=self.n_workers)[0]
        return results

    def _get_checks(self):
        """Get dictionary with ``condition: error_message`` pairs that ensure the
        validity of the backend instance. Here checks if ``dask`` module is
        installed in the environment.

        Returns
        -------
        dict
            dictionary with ``condition: error_message`` pairs that will get
            checked during ``_validate()`` run
        """
        base_checks = super()._get_checks()
        checks = {
            is_installed("dask"):
            ("module 'dask' is missing. Please install 'dask': "
             "https://docs.dask.org/en/stable/install.html")
        }
        return base_checks | checks
