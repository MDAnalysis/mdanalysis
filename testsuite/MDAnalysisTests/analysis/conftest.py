import pytest

from MDAnalysis.analysis.base import AnalysisBase, AnalysisFromFunction
from MDAnalysisTests.analysis.test_base import (
    FrameAnalysis,
    IncompleteAnalysis,
    OldAPIAnalysis,
)
from MDAnalysis.analysis.rms import RMSD, RMSF
from MDAnalysis.lib.util import is_installed
from MDAnalysis.analysis.gnm import GNMAnalysis


def params_for_cls(cls, exclude: list[str] = None):
    """
    This part contains fixtures for simultaneous testing
    of all available (=installed & supported) backends
    for analysis subclasses.

    If for some reason you want to limit backends,
    simply pass "exclude: list[str]" to the function
    that parametrizes fixture.

    Parameters
    ----------
    exclude : list[str], optional
        list of backends to exclude from parametrization, by default None

    Returns
    -------
    dict
        dictionary with all tested keyword combinations for the run
    """
    exclude = [] if exclude is None else exclude
    possible_backends = cls.get_supported_backends()
    installed_backends = [
        b for b in possible_backends if is_installed(b) and b not in exclude
    ]

    params = [
        pytest.param({
            "backend": backend,
            "n_workers": nproc
        }, ) for backend in installed_backends for nproc in (2, )
        if backend != "serial"
    ]
    params.extend([{"backend": "serial"}])
    return params


@pytest.fixture(scope='module', params=params_for_cls(FrameAnalysis))
def client_FrameAnalysis(request):
    return request.param


@pytest.fixture(scope='module', params=params_for_cls(AnalysisBase))
def client_AnalysisBase(request):
    return request.param


@pytest.fixture(scope='module', params=params_for_cls(AnalysisFromFunction))
def client_AnalysisFromFunction(request):
    return request.param


@pytest.fixture(scope='module',
                params=params_for_cls(AnalysisFromFunction,
                                      exclude=['multiprocessing']))
def client_AnalysisFromFunctionAnalysisClass(request):
    return request.param


@pytest.fixture(scope='module', params=params_for_cls(IncompleteAnalysis))
def client_IncompleteAnalysis(request):
    return request.param


@pytest.fixture(scope='module', params=params_for_cls(OldAPIAnalysis))
def client_OldAPIAnalysis(request):
    return request.param


@pytest.fixture(scope='module', params=params_for_cls(RMSD))
def client_RMSD(request):
    return request.param


@pytest.fixture(scope='module', params=params_for_cls(RMSF))
def client_RMSF(request):
    return request.param


@pytest.fixture(scope='module', params=params_for_cls(GNMAnalysis))
def client_GNMAnalysis(request):
    return request.param
