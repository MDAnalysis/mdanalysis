import pytest
import functools
import importlib
import multiprocessing
from contextlib import contextmanager

from MDAnalysis.analysis.base import AnalysisBase, AnalysisFromFunction
from MDAnalysis.analysis.rms import RMSD, RMSF
from MDAnalysisTests.analysis.test_base import (
    FrameAnalysis,
    IncompleteAnalysis,
    OldAPIAnalysis,
)

from MDAnalysis.analysis.base import ParallelExecutor
from MDAnalysis.lib.util import is_installed
from MDAnalysisTests.util import get_running_dask_client


def create_fixture_params_for(cls: type):
    possible_backends = cls.available_backends
    installed_backends = [b for b in possible_backends if is_installed(b)]

    params = [
        pytest.param(
            {"backend": backend, "n_workers": nproc},
        )
        for backend in installed_backends
        for nproc in (1, 2)
        if backend not in ("dask.distributed", "local")
    ]
    params.extend(
        [
            pytest.param(
                {"backend": "local"},
            )
        ]
    )

    if is_installed("dask.distributed"):
        params.extend([pytest.param({"client": get_running_dask_client()}) for n in (1, 2)])

    return params


@contextmanager
def set_tmpdir_as_cwd():
    import os
    from pathlib import Path
    import tempfile

    origin = Path().absolute()
    with tempfile.TemporaryDirectory() as tmpdirname:
        try:
            os.chdir(tmpdirname)
            yield
        finally:
            os.chdir(origin)


@contextmanager
def LocalClient(*args, **kwargs):
    from dask.distributed import LocalCluster, Client

    cluster = LocalCluster(*args, **kwargs)
    client = Client(cluster)
    yield client
    client.close()
    cluster.close()


@pytest.fixture(scope="package", autouse=True)
def setup_client():
    if is_installed("dask") and is_installed("dask.distributed"):
        from multiprocessing import cpu_count
        import warnings

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with set_tmpdir_as_cwd():
                with LocalClient(n_workers=cpu_count(), processes=True, threads_per_worker=1) as client:
                    yield client
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
