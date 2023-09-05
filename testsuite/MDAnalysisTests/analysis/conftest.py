import pytest
import inspect
import functools
import importlib
import multiprocessing
import warnings

from MDAnalysis.analysis.base import AnalysisBase, AnalysisFromFunction
from MDAnalysisTests.analysis.test_base import (
    FrameAnalysis,
    IncompleteAnalysis,
    OldAPIAnalysis,
)
from MDAnalysis.analysis.rms import RMSD, RMSF

from MDAnalysis.analysis.parallel import BackendDaskDistributed
from MDAnalysis.lib.util import is_installed


@pytest.fixture(scope="module")
def dask_client_1():
    if not is_installed("dask"):
        yield "NoClient"
    else:
        import dask

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            lc = dask.distributed.LocalCluster(
                n_workers=1, processes=True, threads_per_worker=1, dashboard_address=None
            )
            client = dask.distributed.Client(lc)
            yield client
            lc.close()
            client.close()


@pytest.fixture(scope="module")
def dask_client_2():
    if not is_installed("dask"):
        yield "NoClient"
    else:
        import dask

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)

            lc = dask.distributed.LocalCluster(
                n_workers=1, processes=True, threads_per_worker=1, dashboard_address=None
            )
            client = dask.distributed.Client(lc)
            yield client
            lc.close()
            client.close()


def generate_client_fixture(cls):
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
    params.extend([{"backend": "local"}, {"backend": "local", "n_workers": 1}])
    if is_installed("dask.distributed") and "dask.distributed" in possible_backends:
        params.extend(["dask_client_1", "dask_client_2"])

    @pytest.fixture(scope="module", params=params)
    def generated_fixture(request, dask_client_1, dask_client_2):
        if request.param == "dask_client_1":
            request.param = {"backend": BackendDaskDistributed(client=dask_client_1)}
        elif request.param == "dask_client_2":
            request.param = {"backend": BackendDaskDistributed(client=dask_client_2)}
        return request.param

    return generated_fixture


def inject_testing_fixture(fixture_name: str, class_name: type):
    """Dynamically inject a fixture at runtime"""
    # we need the caller's global scope for this hack to work hence the use of the inspect module
    caller_globals = inspect.stack()[1][0].f_globals
    # for an explanation of this trick and why it works go here: https://github.com/pytest-dev/pytest/issues/2424
    caller_globals[fixture_name] = generate_client_fixture(class_name)


classes = [AnalysisBase, AnalysisFromFunction, FrameAnalysis, IncompleteAnalysis, OldAPIAnalysis, RMSD, RMSF]
for cls in classes:
    name = cls.__name__
    fixture_name = f"client_{name}"
    inject_testing_fixture(fixture_name=fixture_name, class_name=cls)
