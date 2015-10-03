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
    """Base class for defining multi frame analysis

    Defines common functions for setting up frames to analyse
    """
    def _setup_frames(self, trajectory, start=0, stop=-1, skip=1):
        """Look for the keywords which control trajectory iteration
        and build the list of frames to analyse.

        Returns an array of the identified frames
        """
        # We already have some nice frame checking in base.Reader
        # we should really pull that into lib so everyone can enjoy

        self._trajectory = trajectory

        nframes = len(self._trajectory)

        if stop < 0:
            stop += nframes

        logger.debug("_setup_frames")
        logger.debug(" * settings: start {} stop {} skip {}".format(
            start, stop, skip))

        frames = xrange(start, stop+1, skip)

        logger.debug(" * identified frames:\n"
                     " * {}".format(frames))
        self.frames = frames

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
        for frame in self.frames:
            logger.info("Seeking frame {}".format(frame))
            self._ts = self._trajectory[frame]
            print("Doing frame {} of {}".format(frame, self.frames[-1]))
            logger.info("--> Doing single frame")
            self._single_frame()
        logger.info("Applying normalisation")
        self._normalise()


