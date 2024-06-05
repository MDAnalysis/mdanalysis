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

If you want to run two or more different analyses on the same trajectory you
can also efficently combine them using the
:class:`MDAnalysis.analysis.base.AnalysisCollection` class.


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
from collections import UserDict
import inspect
import logging
import itertools

import numpy as np
from MDAnalysis import coordinates
from MDAnalysis.core.groups import AtomGroup
from MDAnalysis.lib.log import ProgressBar

logger = logging.getLogger(__name__)


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
            raise AttributeError(f"'{key}' is a protected dictionary "
                                 "attribute")
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
        if attr == 'data':
            super().__setattr__(attr, val)
        else:
            self.__setitem__(attr, val)

    def __getattr__(self, attr):
        try:
            return self[attr]
        except KeyError as err:
            raise AttributeError("'Results' object has no "
                                 f"attribute '{attr}'") from err

    def __delattr__(self, attr):
        try:
            del self[attr]
        except KeyError as err:
            raise AttributeError("'Results' object has no "
                                 f"attribute '{attr}'") from err

    def __getstate__(self):
        return self.data

    def __setstate__(self, state):
        self.data = state


class AnalysisCollection:
    """
    Class for running a collection of analysis classes on a single trajectory.

    Running a collection of analyses with ``AnalysisCollection`` can result in
    a speedup compared to running the individual analyses since the trajectory
    loop ins only performed once.

    The class assumes that each analysis is a child of
    :class:`MDAnalysis.analysis.base.AnalysisBase`. Additionally, the
    trajectory of all `analysis_instances` must be the same.

    By default, it is ensured that all analysis instances use the
    *same original* timestep and not an altered one from a previous analysis
    object.

    Parameters
    ----------
    *analysis_instances : tuple
        List of analysis classes to run on the same trajectory.

    Raises
    ------
    AttributeError
        If all the provided ``analysis_instances`` do not work on the same
        trajectory.
    AttributeError
        If an ``analysis_object`` is not a child of
        :class:`MDAnalysis.analysis.base.AnalysisBase`.

    Example
    -------
    .. code-block:: python

        import MDAnalysis as mda
        from MDAnalysis.analysis.rdf import InterRDF
        from MDAnalysis.analysis.base import AnalysisCollection
        from MDAnalysisTests.datafiles import TPR, XTC

        u = mda.Universe(TPR, XTC)

        # Select atoms
        ag_O = u.select_atoms("name O")
        ag_H = u.select_atoms("name H")

        # Create individual analysis instances
        rdf_OO = InterRDF(ag_O, ag_O)
        rdf_OH = InterRDF(ag_H, ag_H)

        # Create a collection for common trajectory
        collection = AnalysisCollection(rdf_OO, rdf_OH)

        # Run the collected analysis
        collection.run(start=0, stop=100, step=10)

        # Results are stored in the individual instances
        print(rdf_OO.results)
        print(rdf_OH.results)


    .. versionadded:: 2.5.0

    """

    def __init__(self, *analysis_instances):
        for analysis_object in analysis_instances:
            if analysis_instances[0]._trajectory != analysis_object._trajectory:
                raise ValueError(
                    "`analysis_instances` do not have the same trajectory."
                )
            if not isinstance(analysis_object, AnalysisBase):
                raise AttributeError(
                    f"Analysis object {analysis_object} is "
                    "not a child of `AnalysisBase`."
                )

        self._analysis_instances = analysis_instances

    def run(
        self,
        start=None,
        stop=None,
        step=None,
        frames=None,
        verbose=None,
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
        verbose : bool, optional
            Turn on verbosity
        """

        # Ensure compatibility with API of version 0.15.0
        if not hasattr(self, "_analysis_instances"):
            self._analysis_instances = (self,)

        logger.info("Choosing frames to analyze")
        # if verbose unchanged, use class default
        verbose = getattr(self, "_verbose", False) if verbose is None else verbose

        logger.info("Starting preparation")

        for analysis_object in self._analysis_instances:
            analysis_object._setup_frames(
                analysis_object._trajectory,
                start=start,
                stop=stop,
                step=step,
                frames=frames,
            )
            analysis_object._prepare()

        logger.info(
            f"Starting analysis loop over {self._analysis_instances[0].n_frames} "
            "trajectory frames."
        )

        for i, ts in enumerate(
            ProgressBar(self._analysis_instances[0]._sliced_trajectory, verbose=verbose)
        ):
            ts_original = ts.copy()

            for analysis_object in self._analysis_instances:
                # Set attributes before calling `_single_frame()`. Setting
                # these attributes explicitly is mandatory so that each
                # instance can access the information of the current timestep.
                analysis_object._frame_index = i
                analysis_object._ts = ts
                analysis_object.frames[i] = ts.frame
                analysis_object.times[i] = ts.time

                # Call the actual analysis of each instance.
                analysis_object._single_frame()

                ts = ts_original

        logger.info("Finishing up")

        for analysis_object in self._analysis_instances:
            analysis_object._conclude()

        return self


class AnalysisBase(AnalysisCollection):
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

    """

    def __init__(self, trajectory, verbose=False, **kwargs):
        self._trajectory = trajectory
        self._verbose = verbose
        self.results = Results()
        super().__init__(self)

    def _setup_frames(self, trajectory, start=None, stop=None, step=None,
                      frames=None):
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
        frames : array_like, optional
            array of integers or booleans to slice trajectory; `frames` can
            only be used *instead* of `start`, `stop`, and `step`. Setting
            *both* `frames` and at least one of `start`, `stop`, `step` to a
            non-default value will raise a :exc:`ValueError`.

            .. versionadded:: 2.2.0

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
                raise ValueError("start/stop/step cannot be combined with "
                                 "frames")
            slicer = frames
        else:
            start, stop, step = trajectory.check_slice_indices(start, stop,
                                                               step)
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

    def run(self, start=None, stop=None, step=None, frames=None,
            verbose=None, *, progressbar_kwargs=None):
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

        progressbar_kwargs : dict, optional
            ProgressBar keywords with custom parameters regarding progress bar
            position, etc; see :class:`MDAnalysis.lib.log.ProgressBar` for full
            list.


        .. versionchanged:: 2.2.0
            Added ability to analyze arbitrary frames by passing a list of
            frame indices in the `frames` keyword argument.

        .. versionchanged:: 2.5.0
            Add `progressbar_kwargs` parameter,
            allowing to modify description, position etc of tqdm progressbars
        """
        logger.info("Choosing frames to analyze")
        # if verbose unchanged, use class default
        verbose = getattr(self, '_verbose',
                          False) if verbose is None else verbose

        self._setup_frames(self._trajectory, start=start, stop=stop,
                           step=step, frames=frames)
        logger.info("Starting preparation")
        self._prepare()
        logger.info("Starting analysis loop over %d trajectory frames",
                    self.n_frames)
        if progressbar_kwargs is None:
            progressbar_kwargs = {}

        for i, ts in enumerate(ProgressBar(
                self._sliced_trajectory,
                verbose=verbose,
                **progressbar_kwargs)):
            self._frame_index = i
            self._ts = ts
            self.frames[i] = ts.frame
            self.times[i] = ts.time
            self._single_frame()
        logger.info("Finishing up")
        self._conclude()
        return self


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

    def __init__(self, function, trajectory=None, *args, **kwargs):
        if (trajectory is not None) and (not isinstance(
                trajectory, coordinates.base.ProtoReader)):
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

        super().__init__(trajectory)

    def _prepare(self):
        self.results.timeseries = []

    def _single_frame(self):
        self.results.timeseries.append(self.function(*self.args,
                                                     **self.kwargs))

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
            super().__init__(function, trajectory, *args, **kwargs)

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
    base_kwargs = {name: val
                   for name, val in zip(base_argspec.args[-n_base_defaults:],
                                        base_argspec.defaults)}

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
                "Now allowed are: {}".format(base_kw, base_kwargs.keys()))

    base_args = {}
    for argname, default in base_kwargs.items():
        base_args[argname] = kwargs.pop(argname, default)

    return base_args, kwargs
