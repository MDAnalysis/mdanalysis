from __future__ import unicode_literals
from itertools import chain
# Copyright (c) 2006-2021  Andrey Golovizin
# Copyright (c) 2014  Matthias C. M. Troffaes
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


import os.path  # splitext
import sys

from pybtex.exceptions import PybtexError

if sys.version_info < (3, 10):
    from importlib_metadata import entry_points
else:
    from importlib.metadata import entry_points


class Plugin(object):
    pass


#: default pybtex plugins
_DEFAULT_PLUGINS = {
    "pybtex.database.input": "bibtex",
    "pybtex.database.output": "bibtex",
    "pybtex.backends": "latex",
    "pybtex.style.labels": "number",
    "pybtex.style.names": "plain",
    "pybtex.style.sorting": "none",
    "pybtex.style.formatting": "unsrt",
}

_RUNTIME_PLUGINS = {}


class PluginGroupNotFound(PybtexError):

    def __init__(self, group_name):
        message = u'plugin group {group_name} not found'.format(
            group_name=group_name,
        )
        super(PluginGroupNotFound, self).__init__(message)


class PluginNotFound(PybtexError):

    def __init__(self, plugin_group, name):
        if not name.startswith('.'):
            message = u'plugin {plugin_group}.{name} not found'.format(
                plugin_group=plugin_group,
                name=name,
            )
        else:
            assert plugin_group.endswith('.suffixes')
            message = (
                u'plugin {plugin_group} for suffix {suffix} not found'.format(
                    plugin_group=plugin_group,
                    suffix=name,
                )
            )

        super(PluginNotFound, self).__init__(message)


def _load_entry_point(group, name, use_aliases=False):
    groups = [group, group + '.aliases'] if use_aliases else [group]
    for search_group in groups:
        # first check in runtime plugins registered via register_plugin
        klass = _RUNTIME_PLUGINS.get(search_group, {}).get(name)
        if klass is not None:
            return klass

        # then check installed entry-points
        for entry_point in entry_points(group=search_group, name=name):
            return entry_point.load()
    raise PluginNotFound(group, name)


def find_plugin(plugin_group, name=None, filename=None):
    """Find a :class:`Plugin` class within *plugin_group* which
    matches *name*, or *filename* if *name* is not specified, or
    the default plugin if neither *name* nor *filename* is
    specified.

    If *name* is specified, return the :class:`Plugin` class
    registered under *name*. If *filename* is specified, look at
    its suffix (i.e. extension) and return the :class:`Plugin`
    class registered for this suffix.

    If *name* is not a string, but a plugin class, just return it back.
    (Used to make functions like :function:`make_bibliography` accept already
    loaded plugins as well as plugin names.

    """
    if isinstance(name, type) and issubclass(name, Plugin):
        return name

    if plugin_group not in _DEFAULT_PLUGINS:
        raise PluginGroupNotFound(plugin_group)
    if name:
        return _load_entry_point(plugin_group, name, use_aliases=True)
    elif filename:
        suffix = os.path.splitext(filename)[1]
        return _load_entry_point(plugin_group + '.suffixes', suffix)
    else:
        return _load_entry_point(plugin_group, _DEFAULT_PLUGINS[plugin_group])


def enumerate_plugin_names(plugin_group):
    """Enumerate all plugin names for the given *plugin_group*."""
    runtime_plugins = (name for name in _RUNTIME_PLUGINS.get(plugin_group, {}))
    ep_plugins = (ep.name for ep in entry_points(group=plugin_group))
    return chain(runtime_plugins, ep_plugins)


def register_plugin(plugin_group, name, klass, force=False):
    """Register a plugin on the fly.

    This works by adding *klass* as a pybtex entry point, under entry
    point group *plugin_group* and entry point name *name*.

    To register a suffix, append ".suffixes" to the plugin group, and
    *name* is then simply the suffix, which should start with a
    period.

    To register an alias plugin name, append ".aliases" to the plugin group.
    Aliases work just as real names, but are not returned by
    :function:`enumerate_plugin_names` and not advertised in the command line
    help.

    If *force* is ``False``, then existing entry points are not
    overwritten. If an entry point with the given group and name
    already exists, then returns ``False``, otherwise returns
    ``True``.

    If *force* is ``True``, then existing entry points are
    overwritten, and the function always returns ``True``.

    """
    if plugin_group.endswith(".suffixes"):
        base_group, _ = plugin_group.rsplit(".", 1)
        if not name.startswith('.'):
            raise ValueError("a suffix must start with a period")
    elif plugin_group.endswith(".aliases"):
        base_group, _ = plugin_group.rsplit(".", 1)
    else:
        base_group = plugin_group
    if base_group not in _DEFAULT_PLUGINS:
        raise PluginGroupNotFound(base_group)

    if len(entry_points(group=plugin_group, name=name)) > 0 and not force:
        return False

    plugins = _RUNTIME_PLUGINS.setdefault(plugin_group, {})

    if force or name not in plugins:
        plugins[name] = klass
        return True

    return plugins[name] is klass
