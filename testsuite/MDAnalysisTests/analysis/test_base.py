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

from numpy.testing import assert_equal, assert_allclose

import MDAnalysis as mda
import numpy as np
import pytest
from MDAnalysis.analysis import base, backends
from MDAnalysisTests.datafiles import DCD, PSF, TPR, XTC
from MDAnalysisTests.util import no_deprecated_call
from numpy.testing import assert_almost_equal, assert_equal


class FrameAnalysis(base.AnalysisBase):
    """Just grabs frame numbers of frames it goes over"""

    @classmethod
    def get_supported_backends(cls): return ('serial', 'dask', 'multiprocessing')

    _analysis_algorithm_is_parallelizable = True

    def __init__(self, reader, **kwargs):
        super(FrameAnalysis, self).__init__(reader, **kwargs)
        self.traj = reader

    def _prepare(self):
        self.results.found_frames = []

    def _single_frame(self):
        self.results.found_frames.append(self._ts.frame)

    def _conclude(self):
        self.found_frames = list(self.results.found_frames)

    def _get_aggregator(self):
        return base.ResultsGroup({'found_frames': base.ResultsGroup.ndarray_hstack})

class IncompleteAnalysis(base.AnalysisBase):
    def __init__(self, reader, **kwargs):
        super(IncompleteAnalysis, self).__init__(reader, **kwargs)


class OldAPIAnalysis(base.AnalysisBase):
    """for version 0.15.0"""

    def __init__(self, reader, **kwargs):
        self._setup_frames(reader, **kwargs)

    def _single_frame(self):
        pass

    def _prepare(self):
        self.results = base.Results()


@pytest.fixture(scope='module')
def u():
    return mda.Universe(PSF, DCD)


@pytest.fixture(scope='module')
def u_xtc():
    return mda.Universe(TPR, XTC)  # dt = 100


FRAMES_ERR = 'AnalysisBase.frames is incorrect'
TIMES_ERR = 'AnalysisBase.times is incorrect'

class Parallelizable(base.AnalysisBase):
    _analysis_algorithm_is_parallelizable = True
    @classmethod
    def get_supported_backends(cls): return ('multiprocessing', 'dask')
    def _single_frame(self): pass

class SerialOnly(base.AnalysisBase):
    def _single_frame(self): pass

class ParallelizableWithDaskOnly(base.AnalysisBase):
    _analysis_algorithm_is_parallelizable = True
    @classmethod
    def get_supported_backends(cls): return ('dask',)
    def _single_frame(self): pass

class CustomSerialBackend(backends.BackendBase):
    def apply(self, func, computations):
        return [func(task) for task in computations]

class ManyWorkersBackend(backends.BackendBase):
    def apply(self, func, computations):
        return [func(task) for task in computations]

def test_incompatible_n_workers(u):
    backend = ManyWorkersBackend(n_workers=2)
    with pytest.raises(ValueError):
        FrameAnalysis(u).run(backend=backend, n_workers=3)

@pytest.mark.parametrize('run_class,backend,n_workers', [
    (Parallelizable, 'not-existing-backend', 2),
    (Parallelizable, 'not-existing-backend', None),
    (SerialOnly, 'not-existing-backend', 2),
    (SerialOnly, 'not-existing-backend', None),
    (SerialOnly, 'multiprocessing', 2),
    (SerialOnly, 'dask', None),
    (ParallelizableWithDaskOnly, 'multiprocessing', None),
    (ParallelizableWithDaskOnly, 'multiprocessing', 2),
])
def test_backend_configuration_fails(u, run_class, backend, n_workers):
    u = mda.Universe(TPR, XTC)  # dt = 100
    with pytest.raises(ValueError):
        _ = run_class(u.trajectory).run(backend=backend, n_workers=n_workers, stop=0)

@pytest.mark.parametrize('run_class,backend,n_workers', [
    (Parallelizable, CustomSerialBackend, 2),
    (ParallelizableWithDaskOnly, CustomSerialBackend, 2),
])
def test_backend_configuration_works_when_unsupported_backend(u, run_class, backend, n_workers):
    u = mda.Universe(TPR, XTC)  # dt = 100
    backend_instance = backend(n_workers=n_workers)
    _ = run_class(u.trajectory).run(backend=backend_instance, n_workers=n_workers, stop=0, unsupported_backend=True)

@pytest.mark.parametrize('run_class,backend,n_workers', [
    (Parallelizable, CustomSerialBackend, 1),
    (ParallelizableWithDaskOnly, CustomSerialBackend, 1),
])
def test_custom_backend_works(u, run_class, backend, n_workers):
    backend_instance = backend(n_workers=n_workers)
    u = mda.Universe(TPR, XTC)  # dt = 100
    _ = run_class(u.trajectory).run(backend=backend_instance, n_workers=n_workers, unsupported_backend=True)

@pytest.mark.parametrize('run_class,backend_instance,n_workers', [
    (Parallelizable, map, 1),
    (SerialOnly, list, 1),
    (ParallelizableWithDaskOnly, object, 1),
])
def test_fails_incorrect_custom_backend(u, run_class, backend_instance, n_workers):
    u = mda.Universe(TPR, XTC)  # dt = 100
    with pytest.raises(ValueError):
        _ = run_class(u.trajectory).run(backend=backend_instance, n_workers=n_workers, unsupported_backend=True)

    with pytest.raises(ValueError):
        _ = run_class(u.trajectory).run(backend=backend_instance, n_workers=n_workers)

@pytest.mark.parametrize('run_class,backend,n_workers', [
    (SerialOnly, CustomSerialBackend, 1),
    (SerialOnly, 'multiprocessing', 1),
    (SerialOnly, 'dask', 1),
])
def test_fails_for_unparallelizable(u, run_class, backend, n_workers):
    u = mda.Universe(TPR, XTC)  # dt = 100
    with pytest.raises(ValueError):
        if not isinstance(backend, str):
            backend_instance = backend(n_workers=n_workers)
            _ = run_class(u.trajectory).run(backend=backend_instance, n_workers=n_workers, unsupported_backend=True)
        else:
            _ = run_class(u.trajectory).run(backend=backend, n_workers=n_workers, unsupported_backend=True)

@pytest.mark.parametrize('run_kwargs,frames', [
    ({}, np.arange(98)),
    ({'start': 20}, np.arange(20, 98)),
    ({'stop': 30}, np.arange(30)),
    ({'step': 10}, np.arange(0, 98, 10))
])
def test_start_stop_step_parallel(u, run_kwargs, frames, client_FrameAnalysis):
    # client_FrameAnalysis is defined [here](testsuite/MDAnalysisTests/analysis/conftest.py),
    # and determines a set of parameters ('backend', 'n_workers'), taking only backends
    # that are implemented for a given subclass, to run the test against.
    an = FrameAnalysis(u.trajectory).run(**run_kwargs, **client_FrameAnalysis)
    assert an.n_frames == len(frames)
    assert_equal(an.found_frames, frames)
    assert_equal(an.frames, frames, err_msg=FRAMES_ERR)
    assert_almost_equal(an.times, frames+1, decimal=4, err_msg=TIMES_ERR)


def test_reset_n_parts_to_n_frames(u):
    # Issue #4685
    a = FrameAnalysis(u.trajectory)
    with pytest.warns(UserWarning, match='Set to'):
        a.run(backend='multiprocessing',
              start=0,
              stop=1,
              n_workers=2,
              n_parts=2)


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
    assert_allclose(an.times, frames+1, rtol=0, atol=1.5e-4, err_msg=TIMES_ERR)


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


@pytest.mark.parametrize('run_kwargs, frames', [
    ({'frames': [4, 5, 6, 7, 8, 9]}, np.arange(4, 10)),
    ({'frames': [0, 2, 4, 6, 8]}, np.arange(0, 10, 2)),
    ({'frames': [4, 6, 8]}, np.arange(4, 10, 2)),
    ({'frames': [0, 3, 4, 3, 5]}, [0, 3, 4, 3, 5]),
    ({'frames': [True, True, False, True, False, True, True, False, True,
                 False]}, (0, 1, 3, 5, 6, 8)),
])
def test_frame_slice_parallel(run_kwargs, frames, client_FrameAnalysis):
    u = mda.Universe(TPR, XTC)  # dt = 100
    an = FrameAnalysis(u.trajectory).run(**run_kwargs, **client_FrameAnalysis)
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
def test_frame_fail(u, run_kwargs, client_FrameAnalysis):
    an = FrameAnalysis(u.trajectory)
    msg = 'start/stop/step cannot be combined with frames'
    with pytest.raises(ValueError, match=msg):
        an.run(**client_FrameAnalysis, **run_kwargs)

def test_parallelizable_transformations():
    # pick any transformation that would allow 
    # for parallelizable attribute
    from MDAnalysis.transformations import NoJump 
    u = mda.Universe(XTC)
    u.trajectory.add_transformations(NoJump())

    # test that serial works
    FrameAnalysis(u.trajectory).run()

    # test that parallel fails
    with pytest.raises(ValueError):
        FrameAnalysis(u.trajectory).run(backend='multiprocessing')

def test_frame_bool_fail(client_FrameAnalysis):
    u = mda.Universe(TPR, XTC)  # dt = 100
    an = FrameAnalysis(u.trajectory)
    frames = [True, True, False]
    msg = 'boolean index did not match indexed array along (axis|dimension) 0'
    with pytest.raises(IndexError, match=msg):
        an.run(**client_FrameAnalysis, frames=frames)


def test_rewind(client_FrameAnalysis):
    u = mda.Universe(TPR, XTC)  # dt = 100
    an = FrameAnalysis(u.trajectory).run(**client_FrameAnalysis, frames=[0, 2, 3, 5, 9])
    assert_equal(u.trajectory.ts.frame, 0)


def test_frames_times(client_FrameAnalysis):
    u = mda.Universe(TPR, XTC)  # dt = 100
    an = FrameAnalysis(u.trajectory).run(start=1, stop=8, step=2, **client_FrameAnalysis)
    frames = np.array([1, 3, 5, 7])
    assert an.n_frames == len(frames)
    assert_equal(an.found_frames, frames)
    assert_equal(an.frames, frames, err_msg=FRAMES_ERR)
    assert_allclose(an.times, frames*100, rtol=0, atol=1.5e-4, err_msg=TIMES_ERR)


def test_verbose(u):
    a = FrameAnalysis(u.trajectory, verbose=True)
    assert a._verbose


def test_warn_nparts_nworkers(u):
    a = FrameAnalysis(u.trajectory)
    with pytest.warns(UserWarning):
        a.run(backend='multiprocessing', n_workers=3, n_parts=2)


@pytest.mark.parametrize(
    "classname,is_parallelizable",
    [
        (base.AnalysisBase, False),
        (base.AnalysisFromFunction, True),
        (FrameAnalysis, True)
    ]
)
def test_not_parallelizable(u, classname, is_parallelizable):
    assert classname._analysis_algorithm_is_parallelizable == is_parallelizable


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


def test_progressbar_multiprocessing(u):
    with pytest.raises(ValueError):
        FrameAnalysis(u.trajectory).run(backend='multiprocessing', verbose=True)


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
def test_AnalysisFromFunction(u, start, stop, step, nframes, client_AnalysisFromFunction):
    # client_AnalysisFromFunction is defined [here](testsuite/MDAnalysisTests/analysis/conftest.py),
    # and determines a set of parameters ('backend', 'n_workers'), taking only backends
    # that are implemented for a given subclass, to run the test against.
    ana1 = base.AnalysisFromFunction(simple_function, mobile=u.atoms)
    ana1.run(start=start, stop=stop, step=step, **client_AnalysisFromFunction)

    ana2 = base.AnalysisFromFunction(simple_function, u.atoms)
    ana2.run(start=start, stop=stop, step=step, **client_AnalysisFromFunction)

    ana3 = base.AnalysisFromFunction(simple_function, u.trajectory, u.atoms)
    ana3.run(start=start, stop=stop, step=step, **client_AnalysisFromFunction)

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


def test_AnalysisFromFunction_args_content(u, client_AnalysisFromFunction):
    protein = u.select_atoms('protein')
    masses = protein.masses.reshape(-1, 1)
    another = mda.Universe(TPR, XTC).select_atoms("protein")
    ans = base.AnalysisFromFunction(mass_xyz, protein, another, masses)
    assert len(ans.args) == 3
    result = np.sum(ans.run(**client_AnalysisFromFunction).results.timeseries)
    assert_allclose(result, -317054.67757345125, rtol=0, atol=1.5e-6)
    assert_almost_equal(result, -317054.67757345125, decimal=6)
    assert (ans.args[0] is protein) and (ans.args[1] is another)
    assert ans._trajectory is protein.universe.trajectory


def test_analysis_class(client_AnalysisFromFunctionAnalysisClass):
    ana_class = base.analysis_class(simple_function)
    assert issubclass(ana_class, base.AnalysisBase)
    assert issubclass(ana_class, base.AnalysisFromFunction)

    u = mda.Universe(PSF, DCD)
    step = 2
    ana = ana_class(u.atoms).run(step=step, **client_AnalysisFromFunctionAnalysisClass)

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
