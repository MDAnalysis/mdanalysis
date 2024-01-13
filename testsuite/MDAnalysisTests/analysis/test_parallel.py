import pytest
from MDAnalysis.analysis import parallel
from MDAnalysis.lib.util import is_installed


def square(x: int):
    return x**2


def noop(x):
    return x


def upper(s):
    return s.upper()


class Test_Backends:
    @pytest.mark.parametrize(
        "backend_cls,n_workers",
        [
            (parallel.BackendBase, -1),
            (parallel.BackendSerial, None),
            (parallel.BackendMultiprocessing, "string"),
            (parallel.BackendDask, ()),
        ],
    )
    def test_fails_incorrect_n_workers(self, backend_cls, n_workers):
        with pytest.raises(ValueError):
            _ = backend_cls(n_workers=n_workers)

    @pytest.mark.parametrize(
        "func,iterable,answer",
        [
            (square, (1, 2, 3), [1, 4, 9]),
            (square, (), []),
            (noop, list(range(10)), list(range(10))),
            (upper, "asdf", list("ASDF")),
        ],
    )
    def test_all_backends_give_correct_results(self, func, iterable, answer):
        backend_instances = [
            parallel.BackendMultiprocessing(n_workers=2),
            parallel.BackendSerial(n_workers=1),
        ]
        if is_installed("dask"):
            backend_instances.append(parallel.BackendDask(n_workers=2))

        backends = {b: b.apply(func, iterable) for b in backend_instances}
        for answ in backends.values():
            assert answ == answer