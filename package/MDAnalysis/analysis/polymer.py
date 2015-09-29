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
Polymer analysis --- :mod:`MDAnalysis.analysis.polymer`
=======================================================


:Author: Richard J. Gowers
:Year: 2015
:Copyright: GNU Public License v3

This module contains various commonly used tools in analysing polymers.

"""
import numpy as np
import logging

from ..lib.util import blocks_of
from ..lib.distances import calc_bonds

logger = logging.getLogger(__name__)


class _AnalysisBase(object):
    """Base class for defining multi frame analysis

    Defines common functions for setting up frames to analyse
    """
    def _setup_frames(self, **kwargs):
        """Look for the keywords which control trajectory iteration
        and build the list of frames to analyse.

        Returns an array of the identified frames
        """
        self._trajectory = kwargs.pop('traj', None)
        if self._trajectory is None:
            raise ValueError("Must supply the 'traj' keyword")

        nframes = len(self._trajectory)

        start = kwargs.pop('start', 0)
        stop = kwargs.pop('stop', nframes)
        skip = kwargs.pop('skip', 1)

        logger.debug("_setup_frames")
        logger.debug(" * settings: start {} stop {} skip {}".format(
            start, stop, skip))

        frames = np.arange(start, stop, skip)

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


class PersistenceLength(_AnalysisBase):
    """Calculate the persistence length for polymer chains

    Notes
    -----
    This analysis requires that the trajectory supports indexing

    .. versionadded:: 0.11.1
    """
    def __init__(self, atomgroups, **kwargs):
        """Calculate the persistence length for polymer chains

        Parameters
        ----------
        atomgroups : list
            List of atomgroups.  Each atomgroup should represent a single
            polymer chain, ordered in the correct order.
        start : int, optional
            First frame of trajectory to analyse, Default: 0
        stop : int, optional
            Last frame of trajectory to analyse, Default: -1
        skip : int, optional
            Step between frames to analyse, Default: 1

        """
        self._atomgroups = atomgroups

        # Check that all chains are the same length
        lens = [len(ag) for ag in atomgroups]
        chainlength = len(atomgroups[0])
        if not all([l == chainlength for l in lens]):
            raise ValueError("Not all AtomGroups were the same size")

        kwargs.update({'traj': atomgroups[0].universe.trajectory})
        self._setup_frames(**kwargs)

        self._results = np.zeros(chainlength - 1, dtype=np.float32)

    def _single_frame(self):
        # could optimise this by writing a "self dot array"
        # we're only using the upper triangle of np.inner
        # function would accept a bunch of coordinates and spit out the
        # decorrel for that 
        n = len(self._atomgroups[0])

        for chain in self._atomgroups:
            # Vector from each atom to next
            vecs = chain.positions[1:] - chain.positions[:-1]
            # Normalised to unit vectors
            vecs /= np.sqrt((vecs * vecs).sum(axis=1))[:, None]

            inner_pr = np.inner(vecs, vecs)
            for i in xrange(n-1):
                self._results[:(n-1)-i] += inner_pr[i, i:]

    def _normalise(self):
        n = len(self._atomgroups[0])

        norm = np.linspace(n - 1, 1, n - 1)
        norm *= len(self._atomgroups) * len(self.frames)

        self._results /= norm

    def _calc_bond_length(self):
        # calculate average bond length
        bs = []
        for ag in self._atomgroups:
            pos = ag.positions
            b = calc_bonds(pos[:-1], pos[1:]).mean()
            bs.append(b)
        self._lb = np.mean(bs)
