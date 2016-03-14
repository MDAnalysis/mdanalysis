# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver Beckstein
# and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
# Plugin management written by Manuel Nuno Melo, 2015

"""
===========================
Nose plugins for MDAnalysis
===========================

:Author: Manuel Nuno Melo
:Year: 2015

.. versionadded:: 0.11.0

Writing plugins
===============

Nose plugins can be written to expose several hooks along the nose execution
chain. This allows very fine control of testing.

The `nose plugin guidelines`_, and examples therein, are a good starting point
on how to accomplish this, with more details available in the `API specification`_.
You can also base your plugins in the ones here or in the `nose` suite.

This module will take care of registering all plugins with nose. If you write a
plugin all you need to do in `__init__.py` is to import it and add its name
to `__all__`.
You must also include in its module a `plugin_class` variable that points to
the plugin class, or list of plugin classes, that you are implementing (it must
come up after the class definition, otherwise the name will be undefined). Again,
see the existing modules for examples.

Beware that command-line invocation of nose via the `nosetests` script will
bypass the loading of these plugins. If you have any code that depends on that,
have it use the `:func:_check_plugins_loaded` function to act accordingly.

Finally, avoid any imports from `MDAnalysis` in plugin code, unless said imports
happen only at test runtime. See Issue 344 for details.

.. _nose plugin guidelines:
   http://nose.readthedocs.org/en/latest/plugins/writing.html
.. _API specification:
   http://nose.readthedocs.org/en/latest/plugins/interface.html
"""

# Do NOT import MDAnalysis at this level. Tests should do it themselves.
# If MDAnalysis is imported here coverage accounting might fail because all the import
#  code won't be run again under coverage's watch.

# Don't forget to also add your plugin to the import further ahead
__all__ = ['memleak', 'capture_err', 'knownfailure', 'cleanup']

import distutils.version
try:
    import nose
    if distutils.version.LooseVersion(nose.__version__) < distutils.version.LooseVersion('1.3.7'):
        raise ImportError
except ImportError:
    raise ImportError("""nose >= 1.3.7 is required to run the test suite. Please install it first. """
                      """(For example, try "pip install 'nose>=1.3.7'").""")

import nose.plugins.multiprocess

def _nose_config():
    """Function that exposes nose's configuration via a hack in one of our plugins

    The external plugins managed by this module are scanned for the :attr:`config` attribute,
    which at least one should implement upon configuration. Plugins need only to be loaded
    for this to work, not necessarily enabled.
    """
    for plugin in loaded_plugins.values():
        if hasattr(plugin, "config"):
            return plugin.config
    return None

def _check_plugins_loaded():
    """Function that checks whether external plugins were loaded.

    It can be used to ascertain whether nose tests were launched using `nosetests` from the
    command-line, and hence that external plugins aren't available.

    It checks whether one of our plugins was configured.
    """
    return _nose_config() is not None

# This dictionary holds the external plugin instances that are loaded into nose.
# Beware that when multiprocessing this dict won't accurately reflect
# the plugin objects in actual use since they're reinstantiated for
# each child process.
loaded_plugins = dict()

# ADD HERE your plugin import
from . import memleak, capture_err, knownfailure, cleanup

plugin_classes = []
for plugin in __all__:
    cls = globals()[plugin].plugin_class
    if issubclass(cls, nose.plugins.Plugin):
        plugin_classes.append(cls)
    else:
        plugin_classes.extend(cls)

for p_class in plugin_classes:
    if p_class.name is None: # some plugins might not implement a name
        loaded_plugins[p_class.__name__] = p_class()
    else:
        loaded_plugins[p_class.name] = p_class()
    # Add it to multiprocess' list of plugins to instantiate
    try:
        nose.plugins.multiprocess._instantiate_plugins.append(p_class)
    except AttributeError:
        nose.plugins.multiprocess._instantiate_plugins = [p_class]

