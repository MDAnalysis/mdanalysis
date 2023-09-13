# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
"""Analysis building blocks --- :mod:`MDAnalysis.analysis.base`
============================================================

MDAnalysis provides building blocks for creating analysis classes. One can
think of each analysis class as a "tool" that performs a specific analysis over
the trajectory frames and stores the results in the tool.

Analysis classes are derived from :class:`AnalysisBase` by subclassing. This
inheritance provides a common workflow and API for users and makes many
additional features automatically available (such as frame selections and a
verbose progressbar). The important points for analysis classes are:

#. Analysis tools are Python classes derived from :class:`AnalysisBase`.
#. When instantiating an analysis, the :class:`Universe` or :class:`AtomGroup`
   that the analysis operates on is provided together with any other parameters
   that are kept fixed for the specific analysis.
#. The analysis is performed with :meth:`~AnalysisBase.run` method. It has a
   common set of arguments such as being able to select the frames the analysis
   is performed on. The `verbose` keyword argument enables additional output. A
   progressbar is shown by default that also shows an estimate for the
   remaining time until the end of the analysis.
#. Results are always stored in the attribute :attr:`AnalysisBase.results`,
   which is an instance of :class:`Results`, a kind of dictionary that allows
   allows item access via attributes. Each analysis class decides what and how
   to store in :class:`Results` and needs to document it. For time series, the
   :attr:`AnalysisBase.times` contains the time stamps of the analyzed frames.


Example of using a standard analysis tool
-----------------------------------------

For example, the :class:`MDAnalysis.analysis.rms.RMSD` performs a
root-mean-square distance analysis in the following way:

.. code-block:: python

   import MDAnalysis as mda
   from MDAnalysisTests.datafiles import TPR, XTC

   from MDAnalysis.analysis import rms

   u = mda.Universe(TPR, XTC)

   # (2) instantiate analysis
   rmsd = rms.RMSD(u, select='name CA')

   # (3) the run() method can select frames in different ways
   # run on all frames (with progressbar)
   rmsd.run(verbose=True)

   # or start, stop, and step can be used
   rmsd.run(start=2, stop=8, step=2)

   # a list of frames to run the analysis on can be passed
   rmsd.run(frames=[0,2,3,6,9])

   # a list of booleans the same length of the trajectory can be used
   rmsd.run(frames=[True, False, True, True, False, False, True, False,
                    False, True])

   # (4) analyze the results, e.g., plot
   t = rmsd.times
   y = rmsd.results.rmsd[:, 2]   # RMSD at column index 2, see docs

   import matplotlib.pyplot as plt
   plt.plot(t, y)
   plt.xlabel("time (ps)")
   plt.ylabel("RMSD (Å)")


Writing new analysis tools
--------------------------

In order to write new analysis tools, derive a class from :class:`AnalysisBase`
and define at least the :meth:`_single_frame` method, as described in
:class:`AnalysisBase`.

.. SeeAlso::

   The chapter `Writing your own trajectory analysis`_ in the *User Guide*
   contains a step-by-step example for writing analysis tools with
   :class:`AnalysisBase`.


.. _`Writing your own trajectory analysis`:
   https://userguide.mdanalysis.org/stable/examples/analysis/custom_trajectory_analysis.html


Classes
-------

The :class:`Results` and :class:`AnalysisBase` classes are the essential
building blocks for almost all MDAnalysis tools in the
:mod:`MDAnalysis.analysis` module. They aim to be easily useable and
extendable.

:class:`AnalysisFromFunction` and the :func:`analysis_class` functions are
simple wrappers that make it even easier to create fully-featured analysis
tools if only the single-frame analysis function needs to be written.

"""
import inspect
import itertools
import logging
import warnings
from functools import partial
from typing import Callable, Iterable, Sequence, Union

import numpy as np
from MDAnalysis import coordinates
from MDAnalysis.core.groups import AtomGroup
from MDAnalysis.lib.log import ProgressBar

from .parallel import Results, ResultsGroup, BackendDask, BackendMultiprocessing, BackendSerial, BackendBase

logger = logging.getLogger(__name__)


class AnalysisBase(object):
    r"""Base class for defining multi-frame analysis

    The class is designed as a template for creating multi-frame analyses.
    This class will automatically take care of setting up the trajectory
    reader for iterating, and it offers to show a progress meter.
    Computed results are stored inside the :attr:`results` attribute.

    To define a new Analysis, :class:`AnalysisBase` needs to be subclassed
    and :meth:`_single_frame` must be defined. It is also possible to define
    :meth:`_prepare` and :meth:`_conclude` for pre- and post-processing.
    All results should be stored as attributes of the :class:`Results`
    container.

    Parameters
    ----------
    trajectory : MDAnalysis.coordinates.base.ReaderBase
        A trajectory Reader
    verbose : bool, optional
        Turn on more logging and debugging

    Attributes
    ----------
    times: numpy.ndarray
        array of Timestep times. Only exists after calling
        :meth:`AnalysisBase.run`
    frames: numpy.ndarray
        array of Timestep frame indices. Only exists after calling
        :meth:`AnalysisBase.run`
    results: :class:`Results`
        results of calculation are stored after call
        to :meth:`AnalysisBase.run`


    Example
    -------
    .. code-block:: python

       from MDAnalysis.analysis.base import AnalysisBase

       class NewAnalysis(AnalysisBase):
           def __init__(self, atomgroup, parameter, **kwargs):
               super(NewAnalysis, self).__init__(atomgroup.universe.trajectory,
                                                 **kwargs)
               self._parameter = parameter
               self._ag = atomgroup

           def _prepare(self):
               # OPTIONAL
               # Called before iteration on the trajectory has begun.
               # Data structures can be set up at this time
               self.results.example_result = []

           def _single_frame(self):
               # REQUIRED
               # Called after the trajectory is moved onto each new frame.
               # store an example_result of `some_function` for a single frame
               self.results.example_result.append(some_function(self._ag,
                                                                self._parameter))

           def _conclude(self):
               # OPTIONAL
               # Called once iteration on the trajectory is finished.
               # Apply normalisation and averaging to results here.
               self.results.example_result = np.asarray(self.example_result)
               self.results.example_result /=  np.sum(self.result)

    Afterwards the new analysis can be run like this

    .. code-block:: python

       import MDAnalysis as mda
       from MDAnalysisTests.datafiles import PSF, DCD

       u = mda.Universe(PSF, DCD)

       na = NewAnalysis(u.select_atoms('name CA'), 35)
       na.run(start=10, stop=20)
       print(na.results.example_result)
       # results can also be accessed by key
       print(na.results["example_result"])


    .. versionchanged:: 1.0.0
        Support for setting `start`, `stop`, and `step` has been removed. These
        should now be directly passed to :meth:`AnalysisBase.run`.

    .. versionchanged:: 2.0.0
        Added :attr:`results`

    .. versionadded:: 2.7.0
        Added ability to run analysis in parallel using either a built-in backend (`multiprocessing` or `dask`)
        of a custom `parallel.BackendBase` instance with an implemented `apply` method that is used to run
        the computations.

    """

    @classmethod
    @property
    def available_backends(cls):
        """Tuple with backends supported by the core library for a given class.
        User can pass either one of these values as `backend=...` to :meth:`run()` method,
        or a custom object that has `apply` method (see documentation for :meth:`run()`):
         - 'serial': no parallelization
         - 'multiprocessing': parallelization using `multiprocessing.Pool`
         - 'dask': parallelization using `dask.delayed.compute()`. Requires installation of `mdanalysis[parallel]`

        If you want to add your own backend to an existing class, either pass a `parallel.BackendBase` subclass
        (see its documentation to learn how to implement it properly).

        Returns
        -------
        tuple
            names of built-in backends that can be used in :meth:`.run(backend=...)`
        """
        return ("serial",)

    @classmethod
    @property
    def _is_parallelizable(cls):
        """Boolean mark showing that a given class can be parallelizable with split-apply-combine procedure.
        Namely, if we can safely distribute :meth:`_single_frame` to multiple workers and then combine them
        with a proper :meth:`_conclude` call.
        If set to `False`, no backends except for `serial` are supported.

        Returns
        -------
        bool
            if a given `AnalysisBase` subclass is parallelizable with split-apply-combine, or not
        """
        return False

    def __init__(self, trajectory, verbose=False, **kwargs):
        self._trajectory = trajectory
        self._verbose = verbose
        self.results = Results()

    def _setup_frames(self, trajectory, start=None, stop=None, step=None, frames=None):
        """Pass a Reader object and define the desired iteration pattern
        through the trajectory

        Parameters
        ----------
        trajectory : mda.Reader
            A trajectory Reader
        start : int, optional
            start frame of analysis
        stop : int, optional
            stop frame of analysis
        step : int, optional
            number of frames to skip between each analysed frame

        .. versionadded:: 2.2.0
        frames : array_like, optional
            array of integers or booleans to slice trajectory; cannot be
            combined with `start`, `stop`, `step`


        Raises
        ------
        ValueError
            if *both* `frames` and at least one of `start`, `stop`, or `frames`
            is provided (i.e., set to another value than ``None``)


        .. versionchanged:: 1.0.0
            Added .frames and .times arrays as attributes

        .. versionchanged:: 2.2.0
            Added ability to iterate through trajectory by passing a list of
            frame indices in the `frames` keyword argument
        """
        self._trajectory = trajectory
        if frames is not None:
            if not all(opt is None for opt in [start, stop, step]):
                raise ValueError("start/stop/step cannot be combined with frames")
            slicer = frames
        else:
            start, stop, step = trajectory.check_slice_indices(start, stop, step)
            slicer = slice(start, stop, step)
        self._sliced_trajectory = trajectory[slicer]
        self.start = start
        self.stop = stop
        self.step = step
        self.n_frames = len(self._sliced_trajectory)
        self.frames = np.zeros(self.n_frames, dtype=int)
        self.times = np.zeros(self.n_frames)

    def _single_frame(self):
        """Calculate data from a single frame of trajectory

        Don't worry about normalising, just deal with a single frame.
        """
        raise NotImplementedError("Only implemented in child classes")

    def _prepare(self):
        """Set things up before the analysis loop begins"""
        pass  # pylint: disable=unnecessary-pass

    def _conclude(self):
        """Finalize the results you've gathered.

        Called at the end of the :meth:`run` method to finish everything up.
        """
        pass  # pylint: disable=unnecessary-pass

    def _compute(self, indexed_frames: np.ndarray, verbose=None, *, progressbar_kwargs={}) -> "AnalysisBase":
        """Perform the calculation on a balanced slice of frames
        that have been setup prior to that using _setup_computation_groups()

        Parameters
        ----------
        indexed_frames : np.ndarray
            np.ndarray of (n, 2) shape,
            where first column is frame iteration indices
            and second is frame numbers

        verbose : bool, optional
            Turn on verbosity

        progressbar_kwargs : dict, optional
            ProgressBar keywords with custom parameters regarding progress bar position, etc;
            see :class:`MDAnalysis.lib.log.ProgressBar` for full list.

        .. versionadded:: 2.7.0
        """
        logger.info("Choosing frames to analyze")
        # if verbose unchanged, use class default
        verbose = getattr(self, "_verbose", False) if verbose is None else verbose

        logger.info("Starting preparation")
        logger.info("Starting analysis loop over %d trajectory frames", self.n_frames)

        if len(indexed_frames) == 0:  # if `frames` were empty in `run` or `stop=0`
            return self

        frame_indices, frames = (
            indexed_frames[:, 0],
            indexed_frames[:, 1],
        )
        trajectory = self._trajectory[frames]
        for idx, ts in enumerate(ProgressBar(trajectory, verbose=verbose, **progressbar_kwargs)):
            i = frame_indices[idx]
            self._frame_index = i
            self._ts = ts
            self.frames[i] = ts.frame
            self.times[i] = ts.time
            self._single_frame()
        logger.info("Finishing up")
        return self

    def _setup_computation_groups(
        self, n_parts: int, start=None, stop=None, step=None, frames=None
    ) -> list[np.ndarray]:
        """
        Splits the trajectory frames, defined by `start/stop/step` or `frames`, into
        `n_parts` even groups, preserving their indices.

        Parameters
        ----------
        n_parts : int
            number of parts to split the workload into
        start : int, optional
            start frame
        stop : int, optional
            stop frame
        step : int, optional
            number of frames to skip between each analysed frame
        frames : array_like, optional
            array of integers or booleans to slice trajectory; `frames` can
            only be used *instead* of `start`, `stop`, and `step`. Setting
            *both* `frames` and at least one of `start`, `stop`, `step` to a
            non-default value will raise a :exc:`ValueError`.

        Raises
        ------
        ValueError
            if *both* `frames` and at least one of `start`, `stop`, or `frames`
            is provided (i.e., set to another value than ``None``)

        Returns
        -------
        computation_groups : list[np.ndarray]
            list of (n, 2) shaped np.ndarrays with frame indices and numbers

        .. versionadded:: 2.7.0
        """
        if frames is None:
            start, stop, step = self._trajectory.check_slice_indices(start, stop, step)
            used_frames = np.arange(start, stop, step)
        elif not all(opt is None for opt in [start, stop, step]):
            raise ValueError("start/stop/step cannot be combined with frames")
        else:
            used_frames = frames

        if all((isinstance(obj, bool) for obj in used_frames)):
            arange = np.arange(len(used_frames))
            used_frames = arange[used_frames]

        # this numpy thing is similar to list(enumerate(frames))
        enumerated_frames = np.vstack([np.arange(len(used_frames)), used_frames]).T

        return np.array_split(enumerated_frames, n_parts)

    def _configure_backend(self, backend: Union[str, BackendBase], n_workers: int, unsafe: bool = False):
        """Matches a passed backend string value with class attributes :meth:`_is_parallelizable` and :meth:`available_backends`
        to check if downstream calculations can be performed.

        Parameters
        ----------
        backend : Union[str, BackendBase]
            backend to be used:
               - `str` is matched to a builtin backend (one of `serial`, `multiprocessing` and `dask`)
               - `BackendBase` subclass is checked for the presence of `apply` method
        n_workers : int
            positive integer with number of workers (processes, in case of built-in backends) to split the work between
        unsafe : bool, optional
            if you want to run your custom backend on a parallelizable class that has not been tested by developers, by default False

        Returns
        -------
        BackendBase
            instance of a `BackendBase` class that will be used for computations

        Raises
        ------
        ValueError
            if :meth:`_is_parallelizable` is set to `False` but backend is not `serial`
        ValueError
            if `_is_parallelizable` and you're using custom backend instance without specifying `unsafe=True`
        ValueError
            if your trajectory has associated parallelizable transformations but backend is not serial
        ValueError
            if your backend object instance doesn't have an `apply` method

        .. versionadded:: 2.7.0
        """
        builtin_backends = {"serial": BackendSerial, "multiprocessing": BackendMultiprocessing, "dask": BackendDask}

        backend_class = builtin_backends.get(backend, backend)
        available_backend_classes = [builtin_backends.get(b) for b in self.available_backends]

        # check for serial-only classes
        if not self._is_parallelizable and backend_class is not BackendSerial:
            raise ValueError(f"Can not parallelize class {self.__class__}")

        # make sure user enabled 'unsafe=True' for custom classes
        if not unsafe and self._is_parallelizable and backend_class not in available_backend_classes:
            raise ValueError(
                f"Must specify 'unsafe=True' if you want to use a custom {backend_class=} for {self.__class__}"
            )

        # check for the presence of parallelizable transformations
        if backend_class is not BackendSerial and any((t.parallelizable for t in self._trajectory.transformations)):
            raise ValueError("Trajectory should not have associated parallelizable transformations")

        # conclude mapping from string to backend class if it's a builtin backend
        if isinstance(backend, str):
            return backend_class(n_workers=n_workers)

        # or pass along an instance of the class itself after ensuring it has apply method
        if not isinstance(backend, BackendBase) or not hasattr(backend, "apply"):
            raise ValueError(
                f"{backend=} is invalid: should have 'apply' method and be instance of MDAnalysis.analysis.parallel.BackendBase"
            )
        return backend

    def run(
        self,
        start: int = None,
        stop: int = None,
        step: int = None,
        frames: Iterable = None,
        verbose: bool = None,
        n_workers: int = None,
        n_parts: int = None,
        backend: Union[str, BackendBase] = None,
        *,
        unsafe: bool = False,
        progressbar_kwargs={},
    ):
        """Perform the calculation

        Parameters
        ----------
        start : int, optional
            start frame of analysis
        stop : int, optional
            stop frame of analysis
        step : int, optional
            number of frames to skip between each analysed frame
        frames : array_like, optional
            array of integers or booleans to slice trajectory; `frames` can
            only be used *instead* of `start`, `stop`, and `step`. Setting
            *both* `frames` and at least one of `start`, `stop`, `step` to a
            non-default value will raise a :exc:`ValueError`.

        .. versionadded:: 2.2.0
        verbose : bool, optional
            Turn on verbosity

        .. versionadded:: 2.5.0
        progressbar_kwargs : dict, optional
            ProgressBar keywords with custom parameters regarding progress bar position, etc;
            see :class:`MDAnalysis.lib.log.ProgressBar` for full list.
            Available only for `backend='serial'`
        backend : Union[str, BackendBase], optional
            By default, performs calculations in a serial fashion.
            Otherwise, user can choose a backend:
               - `str` is matched to a builtin backend (one of `serial`, `multiprocessing` and `dask`)
               - `BackendBase` subclass is checked for the presence of `apply` method
        n_workers : int
            positive integer with number of workers (processes, in case of built-in backends) to split the work between
        n_parts : int, optional
            number of parts to split computations across. Can be more than number of workers,
            especially if you want to get a progressbar with `dask.distributed` backend.
        unsafe : bool, optional
            if you want to run your custom backend on a parallelizable class that has not been tested by developers, by default False

        .. versionchanged:: 2.2.0
            Added ability to analyze arbitrary frames by passing a list of
            frame indices in the `frames` keyword argument.

        .. versionchanged:: 2.5.0
            Add `progressbar_kwargs` parameter,
            allowing to modify description, position etc of tqdm progressbars
        """
        # default to serial execution
        backend = "serial" if backend is None else backend

        if progressbar_kwargs and backend != "serial":
            raise ValueError("Can not display progressbar with non-serial backend")

        n_workers = 1 if n_workers is None else n_workers

        # set n_parts and check that is has a reasonable value
        n_parts = n_workers if n_parts is None else n_parts

        # do this as early as possible to check client parameters before any computations occur
        executor = self._configure_backend(backend=backend, n_workers=n_workers, unsafe=unsafe)
        if (
            hasattr(executor, "n_workers") and n_parts < executor.n_workers
        ):  # using executor's value here for non-default executors
            warnings.warn(f"likely running not at full capacity: {executor.n_workers=} is greater than {n_parts=}")

        # start preparing the run
        worker_func = partial(self._compute, progressbar_kwargs=progressbar_kwargs, verbose=verbose)
        self._setup_frames(trajectory=self._trajectory, start=start, stop=stop, step=step, frames=frames)
        self._prepare()
        computation_groups = self._setup_computation_groups(
            start=start, stop=stop, step=step, frames=frames, n_parts=n_parts
        )

        # gather all remote objects
        remote_objects: list["AnalysisBase"] = executor.apply(worker_func, computation_groups)
        self.frames = np.array([obj.frames for obj in remote_objects]).sum(axis=0)
        self.times = np.array([obj.times for obj in remote_objects]).sum(axis=0)

        # aggregate results
        remote_results = [obj.results for obj in remote_objects]
        results_aggregator = self._get_aggregator()
        self.results = results_aggregator.merge(remote_results)

        self._conclude()
        return self

    def _get_aggregator(self):
        """Returns a default aggregator that takes entire results if there is a single object, 
        and raises ValueError otherwise

        Returns
        -------
        ResultsGroup
            aggregating object

        .. versionadded:: 2.7.0
        """
        return ResultsGroup(lookup=None)


class AnalysisFromFunction(AnalysisBase):
    r"""Create an :class:`AnalysisBase` from a function working on AtomGroups

    Parameters
    ----------
    function : callable
        function to evaluate at each frame
    trajectory : MDAnalysis.coordinates.Reader, optional
        trajectory to iterate over. If ``None`` the first AtomGroup found in
        args and kwargs is used as a source for the trajectory.
    *args : list
        arguments for `function`
    **kwargs : dict
        arguments for `function` and :class:`AnalysisBase`

    Attributes
    ----------
    results.frames : numpy.ndarray
            simulation frames used in analysis
    results.times : numpy.ndarray
            simulation times used in analysis
    results.timeseries : numpy.ndarray
            Results for each frame of the wrapped function,
            stored after call to :meth:`AnalysisFromFunction.run`.

    Raises
    ------
    ValueError
        if `function` has the same `kwargs` as :class:`AnalysisBase`

    Example
    -------
    .. code-block:: python

        def rotation_matrix(mobile, ref):
            return mda.analysis.align.rotation_matrix(mobile, ref)[0]

        rot = AnalysisFromFunction(rotation_matrix, trajectory,
                                    mobile, ref).run()
        print(rot.results.timeseries)


    .. versionchanged:: 1.0.0
        Support for directly passing the `start`, `stop`, and `step` arguments
        has been removed. These should instead be passed to
        :meth:`AnalysisFromFunction.run`.

    .. versionchanged:: 2.0.0
        Former :attr:`results` are now stored as :attr:`results.timeseries`
    """

    @classmethod
    @property
    def available_backends(cls):
        # multiprocessing won't work because self.function won't get pickled
        return ("serial", "dask")

    @classmethod
    @property
    def _is_parallelizable(cls):
        return True

    def __init__(self, function, trajectory=None, *args, **kwargs):
        if (trajectory is not None) and (not isinstance(trajectory, coordinates.base.ProtoReader)):
            args = (trajectory,) + args
            trajectory = None

        if trajectory is None:
            # all possible places to find trajectory
            for arg in itertools.chain(args, kwargs.values()):
                if isinstance(arg, AtomGroup):
                    trajectory = arg.universe.trajectory
                    break

        if trajectory is None:
            raise ValueError("Couldn't find a trajectory")

        self.function = function
        self.args = args

        self.kwargs = kwargs

        super(AnalysisFromFunction, self).__init__(trajectory)

    def _prepare(self):
        self.results.timeseries = []

    def _get_aggregator(self):
        return ResultsGroup({"timeseries": ResultsGroup.flatten_sequence})

    def _single_frame(self):
        self.results.timeseries.append(self.function(*self.args, **self.kwargs))

    def _conclude(self):
        self.results.frames = self.frames
        self.results.times = self.times
        self.results.timeseries = np.asarray(self.results.timeseries)


def analysis_class(function):
    r"""Transform a function operating on a single frame to an
    :class:`AnalysisBase` class.

    Parameters
    ----------
    function : callable
        function to evaluate at each frame

    Attributes
    ----------
    results.frames : numpy.ndarray
            simulation frames used in analysis
    results.times : numpy.ndarray
            simulation times used in analysis
    results.timeseries : numpy.ndarray
            Results for each frame of the wrapped function,
            stored after call to :meth:`AnalysisFromFunction.run`.

    Raises
    ------
    ValueError
        if `function` has the same `kwargs` as :class:`AnalysisBase`

    Examples
    --------

    For use in a library, we recommend the following style

    .. code-block:: python

        def rotation_matrix(mobile, ref):
            return mda.analysis.align.rotation_matrix(mobile, ref)[0]
        RotationMatrix = analysis_class(rotation_matrix)

    It can also be used as a decorator

    .. code-block:: python

        @analysis_class
        def RotationMatrix(mobile, ref):
            return mda.analysis.align.rotation_matrix(mobile, ref)[0]

        rot = RotationMatrix(u.trajectory, mobile, ref).run(step=2)
        print(rot.results.timeseries)


    .. versionchanged:: 2.0.0
        Former :attr:`results` are now stored as :attr:`results.timeseries`
    """

    class WrapperClass(AnalysisFromFunction):
        def __init__(self, trajectory=None, *args, **kwargs):
            super(WrapperClass, self).__init__(function, trajectory, *args, **kwargs)

    return WrapperClass


def _filter_baseanalysis_kwargs(function, kwargs):
    """
    Create two dictionaries with `kwargs` separated for `function` and
    :class:`AnalysisBase`

    Parameters
    ----------
    function : callable
        function to be called
    kwargs : dict
        keyword argument dictionary

    Returns
    -------
    base_args : dict
        dictionary of AnalysisBase kwargs
    kwargs : dict
        kwargs without AnalysisBase kwargs

    Raises
    ------
    ValueError
        if `function` has the same `kwargs` as :class:`AnalysisBase`

    """
    try:
        # pylint: disable=deprecated-method
        base_argspec = inspect.getfullargspec(AnalysisBase.__init__)
    except AttributeError:
        # pylint: disable=deprecated-method
        base_argspec = inspect.getargspec(AnalysisBase.__init__)

    n_base_defaults = len(base_argspec.defaults)
    base_kwargs = {name: val for name, val in zip(base_argspec.args[-n_base_defaults:], base_argspec.defaults)}

    try:
        # pylint: disable=deprecated-method
        argspec = inspect.getfullargspec(function)
    except AttributeError:
        # pylint: disable=deprecated-method
        argspec = inspect.getargspec(function)

    for base_kw in base_kwargs.keys():
        if base_kw in argspec.args:
            raise ValueError(
                "argument name '{}' clashes with AnalysisBase argument."
                "Now allowed are: {}".format(base_kw, base_kwargs.keys())
            )

    base_args = {}
    for argname, default in base_kwargs.items():
        base_args[argname] = kwargs.pop(argname, default)

    return base_args, kwargs
