# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
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
from __future__ import division

from six.moves import range

from numpy.testing import (
    assert_, dec, raises
)

import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase

from MDAnalysisTests.datafiles import PSF, DCD
from MDAnalysisTests import parser_not_found


class FrameAnalysis(AnalysisBase):
    """Just grabs frame numbers of frames it goes over"""
    def __init__(self, reader, **kwargs):
        super(FrameAnalysis, self).__init__(reader, **kwargs)
        self.traj = reader
        self.frames = []

    def _single_frame(self):
        self.frames.append(self._ts.frame)


class IncompleteAnalysis(AnalysisBase):
    def __init__(self, reader, **kwargs):
        super(IncompleteAnalysis, self).__init__(reader, **kwargs)


class TestAnalysisBase(object):
    @dec.skipif(parser_not_found('DCD'),
                'DCD parser not available. Are you using python 3?')
    def setUp(self):
        # has 98 frames
        self.u = mda.Universe(PSF, DCD)

    def tearDown(self):
        del self.u

    def test_default(self):
        an = FrameAnalysis(self.u.trajectory).run()
        assert_(an.n_frames == len(self.u.trajectory))
        assert_(an.frames == list(range(len(self.u.trajectory))))

    def test_start(self):
        an = FrameAnalysis(self.u.trajectory, start=20).run()
        assert_(an.n_frames == len(self.u.trajectory) - 20)
        assert_(an.frames == list(range(20, len(self.u.trajectory))))

    def test_stop(self):
        an = FrameAnalysis(self.u.trajectory, stop=20).run()
        assert_(an.n_frames == 20)
        assert_(an.frames == list(range(20)))

    def test_step(self):
        an = FrameAnalysis(self.u.trajectory, step=20).run()
        assert_(an.n_frames == 5)
        assert_(an.frames == list(range(98))[::20])

    def test_quiet(self):
        a = FrameAnalysis(self.u.trajectory, quiet=False)
        assert_(not a._quiet)

    @raises(NotImplementedError)
    def test_incomplete_defined_analysis(self):
        IncompleteAnalysis(self.u.trajectory).run()
