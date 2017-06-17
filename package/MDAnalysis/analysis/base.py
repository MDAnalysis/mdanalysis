# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.mdanalysis.org
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
from __future__ import absolute_import
import six
from six.moves import range, zip
import inspect
import logging

import numpy as np
from MDAnalysis import coordinates
from MDAnalysis.core.groups import AtomGroup
from MDAnalysis.lib.log import ProgressMeter, _set_verbose

logger = logging.getLogger(__name__)


class AnalysisBase(object):
    """Base class for defining multi frame analysis

    The class it is designed as a template for creating multiframe analyses.
    This class will automatically take care of setting up the trajectory
    reader for iterating, and it offers to show a progress meter.

    To define a new Analysis, `AnalysisBase` needs to be subclassed
    `_single_frame` must be defined. It is also possible to define
    `_prepare` and `_conclude` for pre and post processing. See the example
    below.

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
               self.result = []

           def _single_frame(self):
               # REQUIRED
               # Called after the trajectory is moved onto each new frame.
               # store result of `some_function` for a single frame
               self.result.append(some_function(self._ag, self._parameter))

           def _conclude(self):
               # OPTIONAL
               # Called once iteration on the trajectory is finished.
               # Apply normalisation and averaging to results here.
               self.result = np.asarray(self.result) / np.sum(self.result)

    Afterwards the new analysis can be run like this.

    .. code-block:: python

       na = NewAnalysis(u.select_atoms('name CA'), 35).run()
       print(na.result)

    """

    def __init__(self, trajectory, start=None,
                 stop=None, step=None, verbose=None, quiet=None):
        """
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
        verbose : bool, optional
            Turn on verbosity
        """
        self._verbose = _set_verbose(verbose, quiet, default=False)
        self._quiet = not self._verbose
        self._setup_frames(trajectory, start, stop, step)

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
        """
        self._trajectory = trajectory
        start, stop, step = trajectory.check_slice_indices(start, stop, step)
        self.start = start
        self.stop = stop
        self.step = step
        self.n_frames = len(range(start, stop, step))
        interval = int(self.n_frames // 100)
        if interval == 0:
            interval = 1

        # ensure _verbose is set when __init__ wasn't called, this is to not
        # break pre 0.16.0 API usage of AnalysisBase
        if not hasattr(self, '_verbose'):
            if hasattr(self, '_quiet'):
                # Here, we are in the odd case where a children class defined
                # self._quiet without going through AnalysisBase.__init__.
                warnings.warn("The *_quiet* attribute of analyses is "
                              "deprecated (from 0.16)use *_verbose* instead.",
                              DeprecationWarning)
                self._verbose = not self._quiet
            else:
                self._verbose = True
                self._quiet = not self._verbose
        self._pm = ProgressMeter(self.n_frames if self.n_frames else 1,
                                 interval=interval, verbose=self._verbose)

    def _single_frame(self):
        """Calculate data from a single frame of trajectory

        Don't worry about normalising, just deal with a single frame.
        """
        raise NotImplementedError("Only implemented in child classes")

    def _prepare(self):
        """Set things up before the analysis loop begins"""
        pass

    def _conclude(self):
        """Finalise the results you've gathered.

        Called at the end of the run() method to finish everything up.
        """
        pass

    def run(self):
        """Perform the calculation"""
        logger.info("Starting preparation")
        self._prepare()
        for i, ts in enumerate(
                self._trajectory[self.start:self.stop:self.step]):
            self._frame_index = i
            self._ts = ts
            # logger.info("--> Doing frame {} of {}".format(i+1, self.n_frames))
            self._single_frame()
            self._pm.echo(self._frame_index)
        logger.info("Finishing up")
        self._conclude()
        return self


class AnalysisFromFunction(AnalysisBase):
    """
    Create an analysis from a function working on AtomGroups

    Attributes
    ----------
    results : ndarray
        results of calculation are stored after call to ``run``

    Example
    -------
    >>> def rotation_matrix(mobile, ref):
    >>>     return mda.analysis.align.rotation_matrix(mobile, ref)[0]

    >>> rot = AnalysisFromFunction(rotation_matrix, trajectory, mobile, ref).run()
    >>> print(rot.results)

    Raises
    ------
    ValueError : if ``function`` has the same kwargs as ``BaseAnalysis``
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
           arugments for ``function`` and ``AnalysisBase``

        """
        if (trajectory is not None) and (not isinstance(
                trajectory, coordinates.base.ProtoReader)):
            args = args + (trajectory,)
            trajectory = None

        if trajectory is None:
            for arg in args:
                if isinstance(arg, AtomGroup):
                    trajectory = arg.universe.trajectory
            # when we still didn't find anything
            if trajectory is None:
                for arg in six.itervalues(kwargs):
                    if isinstance(arg, AtomGroup):
                        trajectory = arg.universe.trajectory

        if trajectory is None:
            raise ValueError("Couldn't find a trajectory")

        self.function = function
        self.args = args
        base_kwargs, self.kwargs = _filter_baseanalysis_kwargs(self.function,
                                                               kwargs)

        super(AnalysisFromFunction, self).__init__(trajectory, **base_kwargs)

    def _prepare(self):
        self.results = []

    def _single_frame(self):
        self.results.append(self.function(*self.args, **self.kwargs))

    def _conclude(self):
        self.results = np.asarray(self.results)


def analysis_class(function):
    """
    Transform a function operating on a single frame to an analysis class

    For an usage in a library we recommend the following style:

    >>> def rotation_matrix(mobile, ref):
    >>>     return mda.analysis.align.rotation_matrix(mobile, ref)[0]
    >>> RotationMatrix = analysis_class(rotation_matrix)

    It can also be used as a decorator:

    >>> @analysis_class
    >>> def RotationMatrix(mobile, ref):
    >>>     return mda.analysis.align.rotation_matrix(mobile, ref)[0]

    >>> rot = RotationMatrix(u.trajectory, mobile, ref, step=2).run()
    >>> print(rot.results)
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
    base_argspec = inspect.getargspec(AnalysisBase.__init__)
    n_base_defaults = len(base_argspec.defaults)
    base_kwargs = {name: val
                   for name, val in zip(base_argspec.args[-n_base_defaults:],
                                        base_argspec.defaults)}

    argspec = inspect.getargspec(function)
    for base_kw in six.iterkeys(base_kwargs):
        if base_kw in argspec.args:
            raise ValueError(
                "argument name '{}' clashes with AnalysisBase argument."
                "Now allowed are: {}".format(base_kw, base_kwargs.keys()))

    base_args = {}
    for argname, default in six.iteritems(base_kwargs):
        base_args[argname] = kwargs.pop(argname, default)

    return base_args, kwargs
