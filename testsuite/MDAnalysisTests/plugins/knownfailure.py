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
# Since we're already doing import checks in MDAnalysisTests and .plugins no need to clutter here
from numpy.testing.noseclasses import KnownFailure, KnownFailureTest
import nose

plugin_class = KnownFailure

def knownfailure(msg="Test skipped due to expected failure", exc_type=AssertionError, mightpass=False):
    """If decorated function raises exception *exc_type* skip test, else raise AssertionError."""
    def knownfailure_decorator(f):
        def inner(*args, **kwargs):
            try:
                f(*args, **kwargs)
            except exc_type:
                # We have to allow for a number of cases where KnownFailureTest won't be properly caught
                #  (running from the command-line nosetests or using too old a multiprocess plugin)
                if _check_plugins_loaded() and loaded_plugins["KnownFailure"].enabled:
                    raise KnownFailureTest(msg)
                else: #Fallback if run from command-line and the plugin isn't loaded 
                    raise nose.SkipTest(msg)
            else:
                if not mightpass:
                    raise AssertionError('Failure expected')
        return nose.tools.make_decorator(f)(inner)
    return knownfailure_decorator


