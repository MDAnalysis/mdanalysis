import sys
if sys.version_info >= (3, 10):
    from importlib.metadata import entry_points, EntryPoint
else:
    from importlib_metadata import entry_points, EntryPoint
from typing import Type, Any, Dict, List

_runtime_plugins: Dict[str, Dict[str, Type]] = {
    'sphinxcontrib.bibtex.style.referencing': {}}


# wrapper to work around missing type annotations for entry_points function
def _entry_points(group: str, name: str) -> List[EntryPoint]:
    return entry_points(group=group, name=name)  # type: ignore


def find_plugin(group: str, name: str) -> Type[Any]:
    """Load a sphinxcontrib-bibtex plugin, either from the runtime store,
    or from the entry points.
    """
    global _runtime_plugins
    if group not in _runtime_plugins:
        raise ImportError(f"plugin group {group} not found")
    try:
        return _runtime_plugins[group][name]
    except KeyError:
        for entry_point in _entry_points(group=group, name=name):
            return entry_point.load()
    raise ImportError(f"plugin {group}.{name} not found")


def register_plugin(group: str, name: str, klass: Type[Any],
                    force: bool = False) -> bool:
    """Register a sphinxcontrib-bibtex plugin into the runtime store."""
    global _runtime_plugins
    if group not in _runtime_plugins:
        raise ImportError(f"plugin group {group} not found")
    eps: List[Any]
    try:
        eps = [_runtime_plugins[group][name]]
    except KeyError:
        eps = _entry_points(group=group, name=name)
    if not eps or force:
        _runtime_plugins[group][name] = klass
        return True
    else:
        return False
