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
Wrapper decorator that bridges numpy's :class:`KnownFailure` plugin and provides known failure
exceptions to account for tests that are expected to fail.

Enhances numpy's own plugin by falling back to :class:`SkipTest` exceptions when the plugin 
isn't loaded (typically, when using nose's command-line `nosetests` script, which
does not allow for runtime loading of external plugins).

The decorator can be used without parentheses in case of default arguments as well as 
a function call: `@knownfailure()`, with parentheses in case of optional arguments.
"""
from __future__ import absolute_import

from MDAnalysisTests.plugins import loaded_plugins, _check_plugins_loaded
import nose

try: # numpy < 1.11 had different plugin/exception names
    from numpy.testing.noseclasses import KnownFailurePlugin, KnownFailureException
except ImportError:
    from numpy.testing.noseclasses import KnownFailure as KnownFailurePlugin
    from numpy.testing.noseclasses import KnownFailureTest as KnownFailureException

plugin_class = KnownFailurePlugin
plugin_class.name = "knownfailure"

def knownfailure(args=None, msg="Test skipped due to expected failure", exc_type=AssertionError, mightpass=False):
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
    if callable(args):
        return knownfailure_decorator(args)
    else:      
        return knownfailure_decorator


