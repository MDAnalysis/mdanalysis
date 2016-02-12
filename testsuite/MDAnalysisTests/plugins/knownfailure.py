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

"""
Wrapper decorator that bridges numpy's :class:`KnownFailure` plugin and provides known failure
exceptions to account for tests that are expected to fail.

Enhances numpy's own plugin by falling back to :class:`SkipTest` exceptions when the plugin 
isn't loaded (typically, when using nose's command-line `nosetests` script, which
does not allow for runtime loading of external plugins).

Beware that the decorator must be used as a function call: `@knownfailure()`, with parentheses 
and, optionally, arguments.
"""

from MDAnalysisTests.plugins import loaded_plugins, _check_plugins_loaded
import nose

try: # numpy < 1.11 had different plugin/exception names
    from numpy.testing.noseclasses import KnownFailurePlugin, KnownFailureException
except ImportError:
    from numpy.testing.noseclasses import KnownFailure as KnownFailurePlugin
    from numpy.testing.noseclasses import KnownFailureTest as KnownFailureException

plugin_class = KnownFailurePlugin
plugin_class.name = "knownfailure"

def knownfailure(msg="Test skipped due to expected failure", exc_type=AssertionError, mightpass=False):
    """If decorated function raises exception *exc_type* skip test, else raise AssertionError."""
    def knownfailure_decorator(f):
        def inner(*args, **kwargs):
            try:
                f(*args, **kwargs)
            except exc_type:
                # We have to allow for a number of cases where KnownFailureException won't be properly caught
                #  (running from the command-line nosetests or using too old a multiprocess plugin)
                if _check_plugins_loaded() and loaded_plugins["knownfailure"].enabled:
                    raise KnownFailureException(msg)
                else: #Fallback if run from command-line and the plugin isn't loaded 
                    raise nose.SkipTest(msg)
            else:
                if not mightpass:
                    raise AssertionError('Failure expected')
        return nose.tools.make_decorator(f)(inner)
    return knownfailure_decorator


