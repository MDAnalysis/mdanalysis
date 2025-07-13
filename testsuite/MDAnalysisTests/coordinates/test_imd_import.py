"""Test for MDAnalysis trajectory reader expectations
"""

import sys
import importlib
import pytest
from types import ModuleType

from MDAnalysis.coordinates.IMD import HAS_IMDCLIENT

if HAS_IMDCLIENT:
    import imdclient
    from imdclient.tests.utils import (
        get_free_port,
        create_default_imdsinfo_v3,
    )
    from imdclient.tests.server import InThreadIMDServer

from MDAnalysis.coordinates.IMD import IMDReader


def test_IMDCLIENT_import(monkeypatch):
    backup = sys.modules.copy()

    try:
        module_name = "imdclient"

        # Create mock modules
        mocked_module = ModuleType(module_name)
        IMDClient_module = ModuleType(f"{module_name}.IMDClient")

        class MockIMDClient:
            pass

        IMDClient_module.IMDClient = MockIMDClient
        mocked_module.IMDClient = IMDClient_module
        mocked_module.__version__ = "0.1.4"

        utils_module = ModuleType(f"{module_name}.utils")
        utils_module.parse_host_port = lambda x: ("localhost", 12345)
        mocked_module.utils = utils_module

        monkeypatch.setitem(sys.modules, module_name, mocked_module)
        monkeypatch.setitem(
            sys.modules, f"{module_name}.IMDClient", IMDClient_module
        )
        monkeypatch.setitem(sys.modules, f"{module_name}.utils", utils_module)

        sys.modules.pop("MDAnalysis.coordinates.IMD", None)

        # check if imdclient is new enough
        import MDAnalysis.coordinates.IMD

        importlib.reload(MDAnalysis.coordinates.IMD)
        from MDAnalysis.coordinates.IMD import HAS_IMDCLIENT

        assert HAS_IMDCLIENT

        # check if imdclient version is too old
        mocked_module.__version__ = "0.0.0"
        importlib.reload(MDAnalysis.coordinates.IMD)
        from MDAnalysis.coordinates.IMD import HAS_IMDCLIENT
        from MDAnalysis.coordinates.IMD import IMDReader as IMDReader_NOClient

        assert not HAS_IMDCLIENT

        # test initialization error
        with pytest.raises(
            ImportError, match="IMDReader requires the imdclient"
        ):
            IMDReader_NOClient("imd://localhost:12345", n_atoms=5)
    finally:
        # Restore sys.modules to avoid side effects on other tests
        sys.modules.clear()
        sys.modules.update(backup)
