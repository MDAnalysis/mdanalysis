# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
# Beckstein and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
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
from six.moves import range

import logging
from MDAnalysis.lib.log import ProgressMeter

logger = logging.getLogger(__name__)


class AnalysisBase(object):
    """Base class for defining multi frame analysis, it is designed as a
    template for creating multiframe analysis. This class will automatically
    take care of setting up the trajectory reader for iterating and offers to
    show a progress meter.

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
                 stop=None, step=None, quiet=True):
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
        quiet : bool, optional
            Turn off verbosity
        """
        self._quiet = quiet
        self._setup_frames(trajectory, start, stop, step)

    def _setup_frames(self, trajectory, start=None,
                      stop=None, step=None):
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
        start, stop, step = trajectory.check_slice_indices(
            start, stop, step)
        self.start = start
        self.stop = stop
        self.step = step
        self.n_frames = len(range(start, stop, step))
        interval = int(self.n_frames // 100)
        if interval == 0:
            interval = 1

        # ensure _quiet is set when __init__ wasn't called, this is to not
        # break pre 0.16.0 API usage of AnalysisBase
        if not hasattr(self, '_quiet'):
            self._quiet = True
        self._pm = ProgressMeter(self.n_frames if self.n_frames else 1,
                                 interval=interval, quiet=self._quiet)

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
