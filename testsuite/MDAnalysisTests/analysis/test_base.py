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
from collections import UserDict
import pickle

import pytest

import numpy as np

from numpy.testing import assert_equal, assert_almost_equal

import MDAnalysis as mda
from MDAnalysis.analysis import base

from MDAnalysisTests.datafiles import PSF, DCD, TPR, XTC
from MDAnalysisTests.util import no_deprecated_call


class Test_Results:

    @pytest.fixture
    def results(self):
        return base.Results(a=1, b=2)

    def test_get(self, results):
        assert results.a == results["a"] == 1

    def test_no_attr(self, results):
        msg = "'Results' object has no attribute 'c'"
        with pytest.raises(AttributeError, match=msg):
            results.c

    def test_set_attr(self, results):
        value = [1, 2, 3, 4]
        results.c = value
        assert results.c is results["c"] is value

    def test_set_key(self, results):
        value = [1, 2, 3, 4]
        results["c"] = value
        assert results.c is results["c"] is value

    @pytest.mark.parametrize('key', dir(UserDict) + ["data"])
    def test_existing_dict_attr(self, results, key):
        msg = f"'{key}' is a protected dictionary attribute"
        with pytest.raises(AttributeError, match=msg):
            results[key] = None

    @pytest.mark.parametrize('key', dir(UserDict) + ["data"])
    def test_wrong_init_type(self, key):
        msg = f"'{key}' is a protected dictionary attribute"
        with pytest.raises(AttributeError, match=msg):
            base.Results(**{key: None})

    @pytest.mark.parametrize('key', ("0123", "0j", "1.1", "{}", "a b"))
    def test_weird_key(self, results, key):
        msg = f"'{key}' is not a valid attribute"
        with pytest.raises(ValueError, match=msg):
            results[key] = None

    def test_setattr_modify_item(self, results):
        mylist = [1, 2]
        mylist2 = [3, 4]
        results.myattr = mylist
        assert results.myattr is mylist
        results["myattr"] = mylist2
        assert results.myattr is mylist2
        mylist2.pop(0)
        assert len(results.myattr) == 1
        assert results.myattr is mylist2

    def test_setitem_modify_item(self, results):
        mylist = [1, 2]
        mylist2 = [3, 4]
        results["myattr"] = mylist
        assert results.myattr is mylist
        results.myattr = mylist2
        assert results.myattr is mylist2
        mylist2.pop(0)
        assert len(results["myattr"]) == 1
        assert results["myattr"] is mylist2

    def test_delattr(self, results):
        assert hasattr(results, "a")
        delattr(results, "a")
        assert not hasattr(results, "a")

    def test_missing_delattr(self, results):
        assert not hasattr(results, "d")
        msg = "'Results' object has no attribute 'd'"
        with pytest.raises(AttributeError, match=msg):
            delattr(results, "d")

    def test_pop(self, results):
        assert hasattr(results, "a")
        results.pop("a")
        assert not hasattr(results, "a")

    def test_update(self, results):
        assert not hasattr(results, "spudda")
        results.update({"spudda": "fett"})
        assert results.spudda == "fett"

    def test_update_data_fail(self, results):
        msg = f"'data' is a protected dictionary attribute"
        with pytest.raises(AttributeError, match=msg):
            results.update({"data": 0})

    def test_pickle(self, results):
        results_p = pickle.dumps(results)
        results_new = pickle.loads(results_p)

    @pytest.mark.parametrize("args, kwargs, length", [
        (({"darth": "tater"},), {}, 1),
        ([], {"darth": "tater"}, 1),
        (({"darth": "tater"},), {"yam": "solo"}, 2),
        (({"darth": "tater"},), {"darth": "vader"}, 1),
    ])
    def test_initialize_arguments(self, args, kwargs, length):
        results = base.Results(*args, **kwargs)
        ref = dict(*args, **kwargs)
        assert ref == results
        assert len(results) == length

    def test_different_instances(self, results):
        new_results = base.Results(darth="tater")
        assert new_results.data is not results.data


class FrameAnalysis(base.AnalysisBase):
    """Just grabs frame numbers of frames it goes over"""

    def __init__(self, reader, **kwargs):
        super(FrameAnalysis, self).__init__(reader, **kwargs)
        self.traj = reader
        self.found_frames = []

    def _single_frame(self):
        self.found_frames.append(self._ts.frame)


class IncompleteAnalysis(base.AnalysisBase):
    def __init__(self, reader, **kwargs):
        super(IncompleteAnalysis, self).__init__(reader, **kwargs)


class OldAPIAnalysis(base.AnalysisBase):
    """for version 0.15.0"""

    def __init__(self, reader, **kwargs):
        self._setup_frames(reader, **kwargs)

    def _single_frame(self):
        pass


@pytest.fixture(scope='module')
def u():
    return mda.Universe(PSF, DCD)


@pytest.fixture(scope='module')
def u_xtc():
    return mda.Universe(TPR, XTC)  # dt = 100


FRAMES_ERR = 'AnalysisBase.frames is incorrect'
TIMES_ERR = 'AnalysisBase.times is incorrect'


@pytest.mark.parametrize('run_kwargs,frames', [
    ({}, np.arange(98)),
    ({'start': 20}, np.arange(20, 98)),
    ({'stop': 30}, np.arange(30)),
    ({'step': 10}, np.arange(0, 98, 10))
])
def test_start_stop_step(u, run_kwargs, frames):
    an = FrameAnalysis(u.trajectory).run(**run_kwargs)
    assert an.n_frames == len(frames)
    assert_equal(an.found_frames, frames)
    assert_equal(an.frames, frames, err_msg=FRAMES_ERR)
    assert_almost_equal(an.times, frames+1, decimal=4, err_msg=TIMES_ERR)


@pytest.mark.parametrize('run_kwargs, frames', [
    ({'frames': [4, 5, 6, 7, 8, 9]}, np.arange(4, 10)),
    ({'frames': [0, 2, 4, 6, 8]}, np.arange(0, 10, 2)),
    ({'frames': [4, 6, 8]}, np.arange(4, 10, 2)),
    ({'frames': [0, 3, 4, 3, 5]}, [0, 3, 4, 3, 5]),
    ({'frames': [True, True, False, True, False, True, True, False, True,
                 False]}, (0, 1, 3, 5, 6, 8)),
])
def test_frame_slice(u_xtc, run_kwargs, frames):
    an = FrameAnalysis(u_xtc.trajectory).run(**run_kwargs)
    assert an.n_frames == len(frames)
    assert_equal(an.found_frames, frames)
    assert_equal(an.frames, frames, err_msg=FRAMES_ERR)


@pytest.mark.parametrize('run_kwargs', [
    ({'start': 4, 'frames': [4, 5, 6, 7, 8, 9]}),
    ({'stop': 6, 'frames': [0, 1, 2, 3, 4, 5]}),
    ({'step': 2, 'frames': [0, 2, 4, 6, 8]}),
    ({'start': 4, 'stop': 7, 'frames': [4, 5, 6]}),
    ({'stop': 6, 'step': 2, 'frames': [0, 2, 4, 6]}),
    ({'start': 4, 'step': 2, 'frames': [4, 6, 8]}),
    ({'start': 0, 'stop': 0, 'step': 0, 'frames': [4, 6, 8]}),
])
def test_frame_fail(u, run_kwargs):
    an = FrameAnalysis(u.trajectory)
    msg = 'start/stop/step cannot be combined with frames'
    with pytest.raises(ValueError, match=msg):
        an.run(**run_kwargs)


def test_frame_bool_fail(u_xtc):
    an = FrameAnalysis(u_xtc.trajectory)
    frames = [True, True, False]
    msg = 'boolean index did not match indexed array along (axis|dimension) 0'
    with pytest.raises(IndexError, match=msg):
        an.run(frames=frames)


def test_rewind(u_xtc):
    FrameAnalysis(u_xtc.trajectory).run(frames=[0, 2, 3, 5, 9])
    assert_equal(u_xtc.trajectory.ts.frame, 0)


def test_frames_times(u_xtc):
    an = FrameAnalysis(u_xtc.trajectory).run(start=1, stop=8, step=2)
    frames = np.array([1, 3, 5, 7])
    assert an.n_frames == len(frames)
    assert_equal(an.found_frames, frames)
    assert_equal(an.frames, frames, err_msg=FRAMES_ERR)
    assert_almost_equal(an.times, frames*100, decimal=4, err_msg=TIMES_ERR)


def test_verbose(u):
    a = FrameAnalysis(u.trajectory, verbose=True)
    assert a._verbose


def test_verbose_progressbar(u, capsys):
    FrameAnalysis(u.trajectory).run()
    _, err = capsys.readouterr()
    expected = ''
    actual = err.strip().split('\r')[-1]
    assert actual == expected


def test_verbose_progressbar_run(u, capsys):
    FrameAnalysis(u.trajectory).run(verbose=True)
    _, err = capsys.readouterr()
    expected = u'100%|██████████| 98/98 [00:00<00:00, 8799.49it/s]'
    actual = err.strip().split('\r')[-1]
    assert actual[:24] == expected[:24]

def test_verbose_progressbar_run_with_kwargs(u, capsys):
    FrameAnalysis(u.trajectory).run(
        verbose=True, progressbar_kwargs={'desc': 'custom'})
    _, err = capsys.readouterr()
    expected = u'custom: 100%|██████████| 98/98 [00:00<00:00, 8799.49it/s]'
    actual = err.strip().split('\r')[-1]
    assert actual[:30] == expected[:30]

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


def test_results_type(u):
    an = FrameAnalysis(u.trajectory)
    assert type(an.results) == base.Results


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

    frames = []
    times = []
    timeseries = []

    for ts in u.trajectory[start:stop:step]:
        frames.append(ts.frame)
        times.append(ts.time)
        timeseries.append(simple_function(u.atoms))

    frames = np.asarray(frames)
    times = np.asarray(times)
    timeseries = np.asarray(timeseries)

    assert np.size(timeseries, 0) == nframes

    for ana in (ana1, ana2, ana3):
        assert_equal(frames, ana.results.frames)
        assert_equal(times, ana.results.times)
        assert_equal(timeseries, ana.results.timeseries)


def mass_xyz(atomgroup1, atomgroup2, masses):
    return atomgroup1.positions * masses


def test_AnalysisFromFunction_args_content(u):
    protein = u.select_atoms('protein')
    masses = protein.masses.reshape(-1, 1)
    another = mda.Universe(TPR, XTC).select_atoms("protein")
    ans = base.AnalysisFromFunction(mass_xyz, protein, another, masses)
    assert len(ans.args) == 3
    result = np.sum(ans.run().results.timeseries)
    assert_almost_equal(result, -317054.67757345125, decimal=6)
    assert (ans.args[0] is protein) and (ans.args[1] is another)
    assert ans._trajectory is protein.universe.trajectory


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

    assert_equal(results, ana.results.timeseries)
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
