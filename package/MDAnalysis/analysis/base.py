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
"""
Analysis building blocks --- :mod:`MDAnalysis.analysis.base`
============================================================

A collection of useful building blocks for creating Analysis
classes.

"""
import inspect
import logging
import itertools

import numpy as np
from MDAnalysis import coordinates
from MDAnalysis.core.groups import AtomGroup
from MDAnalysis.lib.log import ProgressBar

logger = logging.getLogger(__name__)


class Results(dict):
    r"""Container object for storing results.

    Results are dictionaries that provide two ways in which can values be
    accessed: by dict key ``results["value_key"]`` or by object attribute,
    ``results.value_key``. They store all results obatined from an analysis
    after the ``run`` call.

    The current implementation is similar to the :class:`sklearn.utils.Bunch`
    class in `scikit-learn`_.
  
    .. _`scikit-learn`: https://scikit-learn.org/
     

    Examples
    --------
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
    """

    def __init__(self, **kwargs):
        super().__init__(kwargs)

    def __setattr__(self, key, value):
        self[key] = value

    def __dir__(self):
        return self.keys()

    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError as err:
            raise AttributeError("'Results' object has no "
                                 f"attribute '{key}'") from err


class AnalysisBase(object):
    r"""Base class for defining multi frame analysis

    The class it is designed as a template for creating multiframe analyses.
    This class will automatically take care of setting up the trajectory
    reader for iterating, and it offers to show a progress meter.
    Computed results are stored inside the :attr:`results` attribute.

    To define a new Analysis, `AnalysisBase` needs to be subclassed
    :meth:`_single_frame` must be defined. It is also possible to define
    :meth:`_prepare` and :meth:`_conclude` for pre and post processing.
    All results should be stored as attributes of the :class:`Results`
    container. See the example below.

    .. code-block:: python

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
               # store a example_result of `some_function` for a single frame
               self.results.example_result.append(some_function(self._ag,
                                                              self._parameter))

           def _conclude(self):
               # OPTIONAL
               # Called once iteration on the trajectory is finished.
               # Apply normalisation and averaging to results here.
               self.results.example_result = np.asarray(self.example_result)
               self.results.example_result /=  np.sum(self.result)

    Afterwards the new analysis can be run like this.

    .. code-block:: python

       na = NewAnalysis(u.select_atoms('name CA'), 35).run(start=10, stop=20)
       print(na.results.example_result)
       # results can also accessed by key
       print(na.results["example_result"])

    Attributes
    ----------
    times: np.ndarray
        array of Timestep times. Only exists after calling run()
    frames: np.ndarray
        array of Timestep frame indices. Only exists after calling run()
    results: :class:`Results`
        results of calculation are stored after call to ``run``

    """

    def __init__(self, trajectory, verbose=False, **kwargs):
        """
        Parameters
        ----------
        trajectory : mda.Reader
            A trajectory Reader
        verbose : bool, optional
           Turn on more logging and debugging, default ``False``


        .. versionchanged:: 1.0.0
           Support for setting ``start``, ``stop``, and ``step`` has been
           removed. These should now be directly passed to
           :meth:`AnalysisBase.run`.
        """
        self._trajectory = trajectory
        self._verbose = verbose
        self.results = Results()

    def _setup_frames(self, trajectory, start=None, stop=None, step=None):
        """
        Pass a Reader object and define the desired iteration pattern
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


        .. versionchanged:: 1.0.0
            Added .frames and .times arrays as attributes

        """
        self._trajectory = trajectory
        start, stop, step = trajectory.check_slice_indices(start, stop, step)
        self.start = start
        self.stop = stop
        self.step = step
        self.n_frames = len(range(start, stop, step))
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
        """Finalise the results you've gathered.

        Called at the end of the run() method to finish everything up.
        """
        pass  # pylint: disable=unnecessary-pass

    def run(self, start=None, stop=None, step=None, verbose=None):
        """Perform the calculation

        Parameters
        ----------
        start : int, optional
            start frame of analysis
        stop : int, optional
            stop frame of analysis
        step : int, optional
            number of frames to skip between each analysed frame
        verbose : bool, optional
            Turn on verbosity
        """
        logger.info("Choosing frames to analyze")
        # if verbose unchanged, use class default
        verbose = getattr(self, '_verbose',
                          False) if verbose is None else verbose

        self._setup_frames(self._trajectory, start, stop, step)
        logger.info("Starting preparation")
        self._prepare()
        for i, ts in enumerate(ProgressBar(
                self._trajectory[self.start:self.stop:self.step],
                verbose=verbose)):
            self._frame_index = i
            self._ts = ts
            self.frames[i] = ts.frame
            self.times[i] = ts.time
            # logger.info("--> Doing frame {} of {}".format(i+1, self.n_frames))
            self._single_frame()
        logger.info("Finishing up")
        self._conclude()
        return self


class AnalysisFromFunction(AnalysisBase):
    r"""Create an :class:`AnalysisBase` from a function working on AtomGroups

    .. code-block:: python

        def rotation_matrix(mobile, ref):
            return mda.analysis.align.rotation_matrix(mobile, ref)[0]

        rot = AnalysisFromFunction(rotation_matrix, trajectory,
                                   mobile, ref).run()
        print(rot.results.timeseries)

    Attributes
    ----------
    results.frames : numpy.ndarray
            simulatiom frames taken for evaluation
    results.times : numpy.ndarray
            simulatiom times taken for evaluation
    results.timeseries : asarray
        Results for each frame of the underlaying function
        stored after call to ``run``.

    Raises
    ------
    ValueError : if ``function`` has the same kwargs as ``BaseAnalysis``

    .. versionchanged:: 2.0.0
        Former :attr:`results` are now stored as :attr:`results.timeseries`
    """

    def __init__(self, function, trajectory=None, *args, **kwargs):
        """Parameters
        ----------
        function : callable
            function to evaluate at each frame
        trajectory : mda.coordinates.Reader (optional)
            trajectory to iterate over. If ``None`` the first AtomGroup found in
            args and kwargs is used as a source for the trajectory.
        *args : list
            arguments for ``function``
        **kwargs : dict
            arguments for ``function`` and ``AnalysisBase``

        .. versionchanged:: 1.0.0
           Support for directly passing the ``start``, ``stop``, and ``step``
           arguments has been removed. These should instead be passed
           to :meth:`AnalysisFromFunction.run`.

        """
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

        super(AnalysisFromFunction, self).__init__(trajectory)

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

    For an usage in a library we recommend the following style

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

    Attributes
    ----------
    results.frames : numpy.ndarray
            simulatiom frames taken for evaluation
    results.times : numpy.ndarray
            simulatiom times taken for evaluation
    results.timeseries : asarray
        Results for each frame of the underlaying function
        stored after call to ``run``.

    Raises
    ------
    ValueError : if ``function`` has the same kwargs as ``BaseAnalysis``

    .. versionchanged:: 2.0.0
        Former `results` are now stored as `results.timeseries`
    """

    class WrapperClass(AnalysisFromFunction):
        def __init__(self, trajectory=None, *args, **kwargs):
            super(WrapperClass, self).__init__(function, trajectory,
                                               *args, **kwargs)

    return WrapperClass


def _filter_baseanalysis_kwargs(function, kwargs):
    """
    create two dictionaries with kwargs separated for function and AnalysisBase

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
    ValueError : if ``function`` has the same kwargs as ``BaseAnalysis``
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
