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
from MDAnalysis.lib.util import is_installed


def generate_client_fixture(cls):
    possible_backends = cls.available_backends
    installed_backends = [b for b in possible_backends if is_installed(b)]

    params = [
        pytest.param(
            {"backend": backend, "n_workers": nproc},
        )
        for backend in installed_backends
        for nproc in (1, 2)
        if backend != "serial"
    ]
    params.extend([{"backend": "serial"}])

    @pytest.fixture(scope="module", params=params)
    def generated_fixture(request):
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
