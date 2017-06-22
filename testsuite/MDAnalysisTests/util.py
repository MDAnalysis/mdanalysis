# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
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
"""
Useful functions for running tests

"""
from __future__ import absolute_import

try:
    import __builtin__
    builtins_name = '__builtin__'
    importer = __builtin__.__import__
except ImportError:
    import builtins
    builtins_name = 'builtins'
    importer = builtins.__import__

from contextlib import contextmanager
from functools import wraps
import importlib
import mock
import os

from numpy.testing import assert_warns

def block_import(package):
    """Block import of a given package

    eg:

    @block_import('numpy')
    def try_and_do_something():
        import numpy as np  # this will fail!

    Will also block imports of subpackages ie block_import('numpy') should
    block 'import numpy.matrix'

    Shadows the builtin import method, sniffs import requests
    and blocks the designated package.
    """
    def blocker_wrapper(func):
        @wraps(func)
        def func_wrapper(*args, **kwargs):
            with mock.patch('{}.__import__'.format(builtins_name),
                            wraps=importer) as mbi:
                def blocker(*args, **kwargs):
                    if package in args[0]:
                        raise ImportError("Blocked by block_import")
                    else:
                        # returning DEFAULT allows the real function to continue
                        return mock.DEFAULT
                mbi.side_effect = blocker
                func(*args, **kwargs)
        return func_wrapper
    return blocker_wrapper


def executable_not_found(*args):
    """Return ``True`` if none of the executables in args can be found.

    ``False`` otherwise (i.e. at least one was found).

    To be used as the argument of::

    @dec.skipif(executable_not_found("binary_name"), msg="skip test because binary_name not available")
    """
    # This must come here so that MDAnalysis isn't imported prematurely,
    #  which spoils coverage accounting (see Issue 344).
    import MDAnalysis.lib.util
    for name in args:
        if MDAnalysis.lib.util.which(name) is not None:
            return False
    return True


def module_not_found(module):
    try:
        importlib.import_module(module)
    except ImportError:
        return True
    else:
        return False


def parser_not_found(parser_name):
    """Return ``True`` if the parser of the given name cannot be found.

    This allows to skip a test when the parser is unavailable (e.g in python 3
    when the parser is not compatible).

    To be used as the argument of::

    @dec.skipif(parser_not_found('DCD'),
                'DCD parser is not available. Are you on python 3?')
    """
    import MDAnalysis.coordinates
    try:
        getattr(MDAnalysis.coordinates, parser_name)
    except AttributeError:
        return True
    else:
        return False


@contextmanager
def in_dir(dirname):
    """Context manager for safely changing directories.

    Arguments
    ---------
    dirname : string
        directory to change into

    Example
    -------
    Change into a temporary directory and always change back to the
    current one::

      with in_dir("/tmp") as tmpdir:
          # do stuff

    SeeAlso
    -------
    The :mod:`tmpdir` module provides functionality such as :func:`tmpdir.in_tmpdir`
    to create temporary directories that are automatically deleted once they are no
    longer used.
    """

    old_path = os.getcwd()
    os.chdir(dirname)
    yield dirname
    os.chdir(old_path)


def assert_nowarns(warning_class, *args, **kwargs):
    """Fail if the given callable throws the specified warning.

    A warning of class warning_class should NOT be thrown by the callable when
    invoked with arguments args and keyword arguments kwargs.
    If a different type of warning is thrown, it will not be caught.

    Parameters
    ----------
    warning_class : class
        The class defining the warning that `func` is expected to throw.
    func : callable
        The callable to test.
    \*args : Arguments
        Arguments passed to `func`.
    \*\*kwargs : Kwargs
        Keyword arguments passed to `func`.

    Returns
    -------
    True
         if no `AssertionError` is raised

    Note
    ----
    numpy.testing.assert_warn returns the value returned by `func`; we would
    need a second func evaluation so in order to avoid it, only True is
    returned if no assertion is raised.

    SeeAlso
    -------
    numpy.testing.assert_warn

    """
    func = args[0]
    args = args[1:]
    try:
        value = assert_warns(DeprecationWarning, func, *args, **kwargs)
    except AssertionError:
        # a warning was NOT emitted: all good
        return True
    else:
        # There was a warning even though we do not want to see one.
        raise AssertionError("function {0} raises warning of class {1}".format(
            func.__name__, warning_class.__name__))
