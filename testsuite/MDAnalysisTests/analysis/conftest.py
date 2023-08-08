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
        params.extend([pytest.param({"client": dask_client(n_workers=n)}) for n in (1, 2)])

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


def dask_client(n_workers: int = None):
    import socket

    from dask.distributed import Client
    import dask

    n_workers = 1 if n_workers is None else n_workers

    with socket.socket() as s:
        s.bind(("", 0))
        port = s.getsockname()[1]
        with set_tmpdir_as_cwd():
            lc = dask.distributed.LocalCluster(n_workers=n_workers, processes=True, port=port)
            client = dask.distributed.Client(lc)
            return client


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
