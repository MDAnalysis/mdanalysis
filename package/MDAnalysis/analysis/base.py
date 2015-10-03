# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
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
import numpy as np
import logging


logger = logging.getLogger(__name__)


class AnalysisBase(object):
    """Base class for defining multi frame analysis"""
    def _setup_frames(self, trajectory, start=None,
                      stop=None, step=None):
        self._trajectory = trajectory
        start, stop, step = trajectory._check_slice_indices(
            start, stop, step)
        self.start = start
        self.stop = stop
        self.step = step
        self.nframes = len(xrange(start, stop, step))

    def _single_frame(self):
        """Calculate data from a single frame of trajectory
 
        Don't worry about normalising, just deal with a single frame.
        """
        pass

    def _normalise(self):
        """Finalise the results you've gathered.

        Called at the end of the run() method to finish everything up.
        """
        pass

    def run(self):
        """Perform the calculation"""
        for i, ts in enumerate(
                self._trajectory[self.start:self.stop:self.step]):
            self._ts = ts
            logger.info("--> Doing frame {} of {}".format(i+1, self.nframes))
            self._single_frame()
        logger.info("Applying normalisation")
        self._normalise()


