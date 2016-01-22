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
from __future__ import division

from numpy.testing import (
    assert_,
)
from multiprocessing import cpu_count

import MDAnalysis as mda
from MDAnalysis.analysis.base import AnalysisBase

from MDAnalysisTests.datafiles import PSF, DCD


class FrameAnalysis(AnalysisBase):
    """Just grabs frame numbers of frames it goes over"""
    def __init__(self, universe, start=None, stop=None, step=None):
        self._setup_frames(universe,
                           start=start,
                           stop=stop,
                           step=step)

        self.results = {}
        self.results['frames'] = []

    def _single_frame(self, timestep):
        self.results['frames'].append(timestep.frame)

    def _add_other_results(self, other_result):
        self.results['frames'] += other_result['frames']

class TestAnalysisBase(object):
    def setUp(self):
        # has 98 frames
        self.u = mda.Universe(PSF, DCD)
        self.frames = len(self.u.trajectory)

    def tearDown(self):
        del self.u

    def test_default(self):
        an = FrameAnalysis(self.u)
        assert_(an.nframes == self.frames)
        an.run()
        assert_(an.results['frames'] == range(self.frames))

        for cores in range(1, cpu_count() + 1):
                an_par = FrameAnalysis(self.u)
                assert_(an_par.nframes == self.frames)
                an_par.run(parallel=True, nthreads=cores)
                assert_(an_par.results['frames'] == range(self.frames))

    def test_start(self):
        an = FrameAnalysis(self.u, start=20)
        assert_(an.nframes == self.frames - 20)

        an.run()
        assert_(an.results['frames'] == range(20, self.frames))

        for cores in range(1, cpu_count() + 1):
                an_par = FrameAnalysis(self.u, start=20)
                assert_(an_par.nframes == self.frames - 20)
                an_par.run(parallel=True, nthreads=cores)
                assert_(an_par.results['frames'] == range(20, self.frames))        

    def test_stop(self):
        an = FrameAnalysis(self.u, stop=20)
        assert_(an.nframes == 20)

        an.run()
        assert_(an.results['frames'] == range(20))

        for cores in range(1, cpu_count() + 1):
                an_par = FrameAnalysis(self.u, stop=20)
                assert_(an_par.nframes == 20)
                an_par.run(parallel=True, nthreads=cores)
                assert_(an_par.results['frames'] == range(20)) 
        
    def test_step(self):
        an = FrameAnalysis(self.u, step=20)
        assert_(an.nframes == 5)

        an.run()
        assert_(an.results['frames'] == range(98)[::20])

        for cores in range(1, cpu_count() + 1):
                an_par = FrameAnalysis(self.u, step=20)
                assert_(an_par.nframes == 5)
                an_par.run(parallel=True, nthreads=cores)
                assert_(an_par.results['frames'] == range(98)[::20])        
