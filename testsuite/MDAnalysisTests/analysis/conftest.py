import pytest
import functools
import importlib
import multiprocessing
from MDAnalysis.analysis.base import Client, AnalysisBase, AnalysisFromFunction
from MDAnalysis.analysis.align import AverageStructure
from MDAnalysis.analysis.atomicdistances import AtomicDistances

from MDAnalysisTests.analysis.test_base import (
    FrameAnalysis,
    IncompleteAnalysis,
    OldAPIAnalysis,
)

import importlib
from MDAnalysis.analysis.base import Client


def is_installed(modulename):
    """Check if modulename is present"""
    if modulename == "local":
        return True
    try:
        importlib.import_module(modulename)
        return True
    except ImportError:
        return False


def create_fixture_params_for(cls: type):
    possible_backends = cls.available_backends
    installed_backends = [b for b in possible_backends if is_installed(b)]

    params = [
        pytest.param(
            {"backend": backend, "n_workers": nproc},
        )
        for backend in installed_backends
        for nproc in range(1, 3)
        if backend != "dask.distributed"
    ]

    if is_installed("dask.distributed"):
        params.extend([pytest.param({"client": get_running_dask_client()})])

    return params


def get_running_dask_client():
    # the solution found here: https://stackoverflow.com/questions/59070260/dask-client-detect-local-default-cluster-already-running
    # note: the used API is non-public and can change any time
    import dask

    return dask.distributed.client._get_global_client()


def dask_client():
    from dask import distributed

    lc = distributed.LocalCluster(n_workers=2, processes=True)
    client = distributed.Client(lc)
    return client


@pytest.fixture(scope="session", params=(1, 2, multiprocessing.cpu_count()))
def setup_dask_distributed(tmpdir_factory, request):
    import dask

    client = dask.distributed.client._get_global_client()
    if client is None:
        with tmpdir_factory.mktemp("dask_cluster").as_cwd():
            lc = dask.distributed.LocalCluster(n_workers=request.param, processes=True)
            client = dask.distributed.Client(lc)
    yield client
    client.close()
    lc.close()


@pytest.fixture(params=create_fixture_params_for(FrameAnalysis))
def client_FrameAnalysis(request):
    return request.param


@pytest.fixture(params=create_fixture_params_for(IncompleteAnalysis))
def client_IncompleteAnalysis(request):
    return request.param


@pytest.fixture(params=create_fixture_params_for(AnalysisFromFunction))
def client_AnalysisFromFunction(request):
    return request.param

@pytest.fixture(params=create_fixture_params_for(AverageStructure))
def client_AverageStructure(request):
    return request.param
    
@pytest.fixture(params=create_fixture_params_for(AtomicDistances))
def client_AtomicDistances(request):
    return request.param