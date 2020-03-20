# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
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
from __future__ import division, absolute_import

import pytest
from six.moves import range

import numpy as np

from numpy.testing import assert_equal, assert_almost_equal

import MDAnalysis as mda
from MDAnalysis.analysis import base

from MDAnalysisTests.datafiles import PSF, DCD, TPR, XTC
from MDAnalysisTests.util import no_deprecated_call


class FrameAnalysis(base.AnalysisBase):
    """Just grabs frame numbers of frames it goes over"""
    def __init__(self, reader, **kwargs):
        super(FrameAnalysis, self).__init__(reader, **kwargs)
        self.traj = reader
        self.frames = []

    def _single_frame(self):
        self.frames.append(self._ts.frame)


class IncompleteAnalysis(base.AnalysisBase):
    def __init__(self, reader, **kwargs):
        super(IncompleteAnalysis, self).__init__(reader, **kwargs)


class OldAPIAnalysis(base.AnalysisBase):
    """for version 0.15.0"""
    def __init__(self, reader, **kwargs):
        self._setup_frames(reader, **kwargs)

    def _single_frame(self):
        pass


@pytest.fixture()
def u():
    return mda.Universe(PSF, DCD)


def test_default(u):
    an = FrameAnalysis(u.trajectory).run()
    assert an.n_frames == len(u.trajectory)
    assert_equal(an.frames, list(range(len(u.trajectory))))


def test_start(u):
    an = FrameAnalysis(u.trajectory).run(start=20)
    assert an.n_frames == len(u.trajectory) - 20
    assert_equal(an.frames, list(range(20, len(u.trajectory))))


def test_stop(u):
    an = FrameAnalysis(u.trajectory).run(stop=20)
    assert an.n_frames == 20
    assert_equal(an.frames, list(range(20)))


def test_step(u):
    an = FrameAnalysis(u.trajectory).run(step=20)
    assert an.n_frames == 5
    assert_equal(an.frames, list(range(98))[::20])


def test_verbose(u):
    a = FrameAnalysis(u.trajectory, verbose=True)
    assert a._verbose


def test_verbose_run(u):
    a = FrameAnalysis(u.trajectory, verbose=False).run(verbose = True)
    assert not a._verbose
    assert a._pm.verbose
    

def test_incomplete_defined_analysis(u):
    with pytest.raises(NotImplementedError):
        IncompleteAnalysis(u.trajectory).run()


def test_old_api(u):
    OldAPIAnalysis(u.trajectory).run()


def test_filter_baseanalysis_kwargs_VE():
    def bad_f(mobile, verbose=2):
        pass

    kwargs = {'step': 3, 'foo': None}

    with pytest.raises(ValueError):
        base._filter_baseanalysis_kwargs(bad_f, kwargs)


def test_filter_baseanalysis_kwargs():
    def good_f(mobile, ref):
        pass

    kwargs = {'step': 3, 'foo': None}

    base_kwargs, kwargs = base._filter_baseanalysis_kwargs(good_f, kwargs)

    assert 2 == len(kwargs)
    assert kwargs['foo'] == None

    assert len(base_kwargs) == 1
    assert base_kwargs['verbose'] is False


def simple_function(mobile):
    return mobile.center_of_geometry()


@pytest.mark.parametrize('start, stop, step, nframes', [
        (None, None, 2, 49),
        (None, 50, 2, 25),
        (20, 50, 2, 15),
        (20, 50, None, 30)
    ])
def test_AnalysisFromFunction(u, start, stop, step, nframes):
    ana1 = base.AnalysisFromFunction(simple_function, mobile=u.atoms)
    ana1.run(start=start, stop=stop, step=step)

    ana2 = base.AnalysisFromFunction(simple_function, u.atoms)
    ana2.run(start=start, stop=stop, step=step)

    ana3 = base.AnalysisFromFunction(simple_function, u.trajectory, u.atoms)
    ana3.run(start=start, stop=stop, step=step)

    results = []

    for ts in u.trajectory[start:stop:step]:
        results.append(simple_function(u.atoms))

    results = np.asarray(results)

    assert np.size(results, 0) == nframes

    for ana in (ana1, ana2, ana3):
        assert_equal(results, ana.results)


def mass_xyz(atomgroup1, atomgroup2, masses):
        return atomgroup1.positions * masses

def test_AnalysisFromFunction_args_content(u):
    protein = u.select_atoms('protein')
    masses = protein.masses.reshape(-1, 1)
    another = mda.Universe(TPR, XTC).select_atoms("protein")
    ans = base.AnalysisFromFunction(mass_xyz, protein, another, masses)
    assert len(ans.args) == 3
    result = np.sum(ans.run().results)
    assert_almost_equal(result, -317054.67757345125, decimal=6)
    assert (ans.args[0] is protein) and (ans.args[1] is another)
    assert  ans._trajectory is protein.universe.trajectory


def test_analysis_class():
    ana_class = base.analysis_class(simple_function)
    assert issubclass(ana_class, base.AnalysisBase)
    assert issubclass(ana_class, base.AnalysisFromFunction)

    u = mda.Universe(PSF, DCD)
    step = 2
    ana = ana_class(u.atoms).run(step=step)

    results = []
    for ts in u.trajectory[::step]:
        results.append(simple_function(u.atoms))
    results = np.asarray(results)

    assert_equal(results, ana.results)
    with pytest.raises(ValueError):
        ana_class(2)


def test_analysis_class_decorator():
    # Issue #1511
    # analysis_class should not raise
    # a DeprecationWarning
    u = mda.Universe(PSF, DCD)

    def distance(a, b):
        return np.linalg.norm((a.centroid() - b.centroid()))

    Distances = base.analysis_class(distance)

    with no_deprecated_call():
        d = Distances(u.atoms[:10], u.atoms[10:20]).run()
