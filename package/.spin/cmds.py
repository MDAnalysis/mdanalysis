import spin


@spin.util.extend_command(spin.cmds.meson.test, doc="")
def test(*, parent_callback, pytest_args, tests, coverage,
         **kwargs):
    pytest_args = ('MDAnalysisTests',) + pytest_args
    parent_callback(**{"pytest_args": pytest_args, "tests": tests,
                    "coverage": coverage, **kwargs})
