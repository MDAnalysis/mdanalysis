import pytest
import functools
import importlib
import multiprocessing

from MDAnalysis.analysis.base import AnalysisBase, AnalysisFromFunction
from MDAnalysis.analysis.rms import RMSD, RMSF
from MDAnalysisTests.analysis.test_base import (
    FrameAnalysis,
    IncompleteAnalysis,
    OldAPIAnalysis,
)

from MDAnalysis.analysis.base import ParallelExecutor
from MDAnalysisTests.util import get_running_dask_client
from MDAnalysis.lib.util import is_installed


def create_fixture_params_for(cls: type):
    possible_backends = cls.available_backends
    installed_backends = [b for b in possible_backends if is_installed(b)]

    params = [
        pytest.param(
            {"backend": backend, "n_workers": nproc},
        )
        for backend in installed_backends
        for nproc in range(1, 3)
        if backend not in ("dask.distributed", "local")
    ]
    params.extend([{"backend": "local"}, {"backend": "local", "n_workers": 1}])

    if is_installed("dask.distributed"):
        params.extend([pytest.param({"client": get_running_dask_client()})])

    return params


def dask_client():
    from dask import distributed

    lc = distributed.LocalCluster(n_workers=2, processes=True)
    client = distributed.Client(lc)
    return client


@pytest.fixture(scope="session", params=(1, 2, min(3, multiprocessing.cpu_count())), autouse=True)
def setup_dask_distributed(tmpdir_factory, request):
    if is_installed("dask"):
        import dask

        client = dask.distributed.client._get_global_client()
        if client is None:  # meaning, in current process there is no active Client
            from dask.distributed import Client

            try:
                client = Client("localhost:8787", timeout=2)  # try connecting to client of different process
            except OSError:
                with tmpdir_factory.mktemp("dask_cluster").as_cwd():  # set up your own client
                    lc = dask.distributed.LocalCluster(n_workers=request.param, processes=True)
                    client = dask.distributed.Client(lc)
        yield client
        client.close()
        lc.close()
    else:
        yield None


@pytest.fixture(params=create_fixture_params_for(FrameAnalysis))
def client_FrameAnalysis(request):
    return request.param


@pytest.fixture(params=create_fixture_params_for(IncompleteAnalysis))
def client_IncompleteAnalysis(request):
    return request.param


@pytest.fixture(params=create_fixture_params_for(AnalysisFromFunction))
def client_AnalysisFromFunction(request):
    return request.param


@pytest.fixture(params=create_fixture_params_for(RMSF))
def client_RMSF(request):
    return request.param


@pytest.fixture(params=create_fixture_params_for(RMSD))
def client_RMSD(request):
    return request.param
