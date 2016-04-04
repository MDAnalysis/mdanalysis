# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
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

import numpy as np
import logging


logger = logging.getLogger(__name__)


class AnalysisBase(object):
    """Base class for defining multi frame analysis

    The analysis base class is designed as a template for creating
    multiframe analysis.  

    The class implements the following methods:

    _setup_frames(trajectory, start=None, stop=None, step=None)
      Pass a Reader object and define the desired iteration pattern
      through the trajectory
      
    run
      The user facing run method.  Calls the analysis methods
      defined below

    Your analysis can implement the following methods, which are
    called from run:

    _prepare
      Called before iteration on the trajectory has begun.
      Data structures can be set up at this time, however most
      error checking should be done in the __init__

    _single_frame
      Called after the trajectory is moved onto each new frame.

    _conclude
      Called once iteration on the trajectory is finished.
      Apply normalisation and averaging to results here.

    """
    def _setup_frames(self, trajectory, start=None,
                      stop=None, step=None):
        self._trajectory = trajectory
        start, stop, step = trajectory.check_slice_indices(
            start, stop, step)
        self.start = start
        self.stop = stop
        self.step = step
        self.nframes = len(range(start, stop, step))

    def _single_frame(self):
        """Calculate data from a single frame of trajectory
 
        Don't worry about normalising, just deal with a single frame.
        """
        pass

    def _prepare(self):
        """Set things up before the analysis loop begins"""
        pass

    def _conclude(self):
        """Finalise the results you've gathered.

        Called at the end of the run() method to finish everything up.
        """
        pass

    def run(self, **kwargs):
        """Perform the calculation"""
        logger.info("Starting preparation")
        self._prepare()
        for i, ts in enumerate(
                self._trajectory[self.start:self.stop:self.step]):
            self._ts = ts
            #logger.info("--> Doing frame {} of {}".format(i+1, self.nframes))
            self._single_frame()
        logger.info("Finishing up")
        self._conclude()


