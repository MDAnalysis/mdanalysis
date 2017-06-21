# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- http://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
# This module was based on code from the importing module from the PEAK
# package (see http://peak.telecommunity.com/DevCenter/FrontPage). The PEAK
# package is released under the following license, reproduced here:
#
#  Copyright (C) 1996-2004 by Phillip J. Eby and Tyler C. Sarna.
#  All rights reserved.  This software may be used under the same terms
#  as Zope or Python.  THERE ARE ABSOLUTELY NO WARRANTIES OF ANY KIND.
#  Code quality varies between modules, from "beta" to "experimental
#  pre-alpha".  :)
#
# The following list summarizes the modifications to the importing code:
#  - a replacement of lazyModule (import_module, which defers most work to
#    _import_module) is implemented that uses an alternative LazyModule class;
#  - a different LazyModule class is created per instance, so that reverting
#    the __getattribute__ behavior can be done safely;
#  - a function to lazily import module functions was added.


"""
Lazy module loading --- :mod:`MDAnalysis.lib.lazy`
====================================================

Functions and classes for lazy module loading that also delay import errors.
Heavily borrowed from the `importing`_ module.

.. versionadded:: 0.16.2
.. _`importing`: http://peak.telecommunity.com/DevCenter/Importing

Files and directories
---------------------

.. autofunction:: import_module
.. autofunction:: import_function

"""

__all__ = ['import_module', 'import_function']

from types import ModuleType
import sys
try:
    # imp is deprecated since python 3.4 but there's no clear alternative to
    # the lock mechanism, other than to import directly from _imp.
    from imp import acquire_lock, release_lock 
except ImportError:
    from _imp import acquire_lock, release_lock 

import six
from six.moves import reload_module

_MSG = ("{0} attempted to use a functionality that requires module {1}, but "
        "it couldn't be loaded. Please install {2} and retry.")

_MSG_FN = ("{0} attempted to use a functionality that requires function {1} "
           "of module {2}, but it couldn't be found in that module. Please "
           "install a version of {2} that has {1} and retry.")


class LazyModule(ModuleType):
    """Class for lazily-loaded modules that triggers proper loading on access

    Instantiation should be made from a subclass of
    :class:`MDAnalysis.lib.lazy.LazyModule`, with one subclass per instantiated
    module. Regular attribute set/access can then be recovered by setting the
    subclass's :meth:`__getattribute__` and :meth:`__setattribute__` to those
    of :class:`types.ModuleType`.
    """
    # peak.util.imports sets __slots__ to (), but it seems pointless because
    # the base ModuleType doesn't itself set __slots__.
    def __init__(self, modname):
        super(ModuleType, self).__setattr__('__name__', modname)

    def __getattribute__(self, attr):
        # IPython tries to be too clever and constantly inspects, asking for
        #  modules' attrs, which causes premature module loading and unesthetic
        #  internal errors if the lazily-loaded module doesn't exist. Returning
        #  Nones seems to satisfy those needs:
        caller_base = _caller_name().partition('.')[0]
        if run_from_ipython and caller_base in ('inspect', 'IPython'):
            return None
        _load_module(self)
        return ModuleType.__getattribute__(self, attr)

    def __setattr__(self, attr, value):
        _load_module(self)
        return ModuleType.__setattr__(self, attr, value)


def _caller_name(depth=2):
    """Returns the name of the calling namespace

    """
    # the presence of sys._getframe might be implementation-dependent.
    # It isn't that serious if we can't get the caller's name.
    try:
        return sys._getframe(depth).f_globals['__name__']
    except AttributeError:
        return 'MDAnalysis'


def run_from_ipython():
    # Taken from https://stackoverflow.com/questions/5376837
    try:
        __IPYTHON__
        return True
    except NameError:
        return False


def import_module(modname, level='leaf'):
    """Function allowing lazy importing of a module into the namespace

    Parameters
    ----------
    modname : str
         The module to import.
    level : str, optional
         Which submodule reference to return. Either a reference to the 'leaf'
         module (the default) or to the 'base' module. For 'base'::

            MDAnalysis = import_module("MDAnalysis.analysis.distances",
                                       level='base')
            # 'MDAnalysis' becomes defined in the current namespace, with
            #  (sub)attributes 'MDAnalysis.analysis' and
            #  'MDAnalysis.analysis.distances'.
            # Equivalent to:
            import MDAnalysis.analysis.distances

        For 'leaf'::

            distances = import_module("MDAnalysis.analysis.distances",
                                      level='leaf')
            # Only 'distances' becomes set in the current namespace.
            # Equivalent to:
            from MDAnalysis.analysis import distances

    Returns
    -------
    module
        The module specified by *modname*, or its base, depending on *level*.
        The module isn't immediately imported. Instead, a
        :class:`MDAnalysis.lib.lazy.LazyModule` instance is returned. Upon
        access to any of its attributes, the module is finally loaded.

    .. versionadded:: 0.16.2

    """
    mod = _import_module(modname, _caller_name())
    if level == 'base':
        return sys.modules[modname.split('.')[0]]
    elif level == 'leaf':
        return mod
    else:
        raise ValueError("Parameter 'level' must be one of ('base', 'leaf')")


def _import_module(modname, caller_name):
    acquire_lock()
    try:
        fullmodname = modname
        fullsubmodname = None
        # ensure parent module/package is in sys.modules
        # and parent.modname=module, as soon as the parent is imported   
        while modname:
            try:
                mod = sys.modules[modname]
                # We reached a (base) module that's already loaded. Let's stop
                # the cycle.
                modname = ''
            except KeyError:
                class _LazyModule(LazyModule):
                    _mda_lazy_caller_name = caller_name
                mod = sys.modules[modname] = _LazyModule(modname)
            if fullsubmodname:
                ModuleType.__setattr__(mod, submodname,
                                       sys.modules[fullsubmodname])
            fullsubmodname = modname
            modname, _, submodname = modname.rpartition('.')
        return sys.modules[fullmodname]
    finally:
        release_lock()


def import_function(modname, *funcnames):
    """Performs lazy importing of one or more functions into the namespace

    Parameters
    ----------
    modname : str
         The base module from where to import the function(s) in *funcnames*,
         or a full 'module_name.function_name' string.
    funcnames : str (optional)
         The function name(s) to import from the module specified by *modname*.
         If left empty *modname* is assumed to also include the function name
         to import.

    Returns
    -------
    function or list of functions
        If *funcnames* is passed, returns a list of imported functions, one for
        each element in *funcnames*.
        If only *modname* is passed it is assumed to be a full
        'module_name.function_name' string, in which case the imported function
        is returned directly, and not in a list.
        The module specified by *modname* is always imported lazily, via
        :func:`MDAnalysis.lib.lazy.import_module`.
        
    See Also
    --------
    :func:`MDAnalysis.lib.lazy.import_module`

    .. versionadded:: 0.16.2

    """
    if not funcnames:
        # We allow passing a single string as 'modname.funcname',
        # in which case the function is returned directly and not as a list.
        modname, _, funcname = modname.rpartition(".")
        return _import_function(modname, funcname, _caller_name())
    else:
        return [_import_function(modname, fn, _caller_name()) for fn in funcnames]


def _import_function(modname, funcname, caller_name):
    module = _import_module(modname, caller_name)

    def retfun(*args, **kwargs):
        try:
            return getattr(module, funcname)(*args, **kwargs)
        except AttributeError:
            raise AttributeError(_MSG_FN.format(caller_name, funcname, modname))
    return retfun


def _load_module(module):
    """Ensures that a module, and its parents, are properly loaded

    """
    modclass = type(module)
    # We only take care of our own LazyModule instances
    if not issubclass(modclass, LazyModule):
        return
    acquire_lock()
    try:
        try:
            # Alreay-loaded _LazyModule classes lose their
            # _mda_lazy_caller_name attr. No need to redo
            # those cases.
            caller_name = modclass._mda_lazy_caller_name
        except AttributeError:
            return
        modclass.__getattribute__ = ModuleType.__getattribute__
        modclass.__setattr__ = ModuleType.__setattr__
        del modclass._mda_lazy_caller_name

        # First, ensure the parent is loaded
        # (using recursion; negligible chance we'll ever hit a stack limit
        #  in this case).
        parent, _, modname = module.__name__.rpartition('.')
        try:
            if parent:
                _load_module(sys.modules[parent])
                setattr(sys.modules[parent], modname, module)
            # Get Python to do the real import!
            reload_module(module)           
        except:
            # We reset our state
            del modclass.__getattribute__
            del modclass.__setattr__
            modclass._mda_lazy_caller_name = caller_name
            raise
    except (AttributeError, ImportError) as err:
        # Under Python 3 reloading our dummy LazyModule instances causes an
        #  AttributeError if the module can't be found. Would be preferrable if
        #  we could always rely on an ImportError. As it is we vet the
        #  AttributeError as thoroughly as possible.
        if (six.PY3 and isinstance(err, AttributeError) and
                err.args[0] != "'NoneType' object has no attribute 'name'"):
            # Not the AttributeError we were looking for.
            raise
        modname = ModuleType.__getattribute__(module, '__name__')
        base_modname = modname.split(".")[0]
        # Way to silence context tracebacks in Python 3 but with a syntax
        # compatible with Python 2. This would normally be:
        #  raise ImportError(...) from None
        exc = ImportError(_MSG.format(caller_name, modname, base_modname))
        exc.__suppress_context__ = True
        raise exc
    finally:
        release_lock()

