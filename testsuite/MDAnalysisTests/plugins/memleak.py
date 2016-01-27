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

"""Test plugin that reports memleaks on a test-by-test basis.

The plugin works by clearing a test's namespace after it has run, and then
forcing a garbage collection round and checking for uncollectable objects.

Implementation uses the startTest hook to register our memory leak check
as a cleanup to the test.
"""

import gc
import nose
from nose.plugins.errorclass import ErrorClass, ErrorClassPlugin

_leakedobjs = set() # We must keep track of seen leaks to avoid raising multiple 
                    # errors for the same ones.

def memleak_check(test, pretest_attrs):
    """Memory leak check, to be registered as a :class:`unittest.TestCase` cleanup method.

    Registration can be done using :meth:`unittest.TestCase.addCleanup`. The
    test itself and its `__dict__.keys()` must be registered as arguments.

    This function works by clearing a test's namespace after it has run, and
    then forcing a garbage collection round and checking for uncollectable
    objects.
    """
    attrs = []
    for attrname in pretest_attrs:
        try:
            attrs.append((attrname, getattr(test, attrname)))
        except AttributeError:
            pass
    test.__dict__.clear()
    gc.collect()
    latest_leaks = [ obj for obj in gc.garbage if obj not in _leakedobjs ]
    _leakedobjs.update(gc.garbage)
    # Restore the pre-test stuff
    for name, val in attrs:
        setattr(test, name, val)
    if latest_leaks:
        raise MemleakError("GC failed to collect the following: {0}".format(latest_leaks))

class MemleakError(Exception):
    """Exception raised internally to mark the test as memleaking."""
    pass

class Memleak(ErrorClassPlugin):
    """Report memleaks on a test-by-test basis."""
    name = "memleak"
    enabled = False
    memleak = ErrorClass(MemleakError,
                         label='MEMLEAK',
                         isfailure=True)

    def configure(self, options, conf):
        super(Memleak, self).configure(options, conf)
        self.config = conf # This will let other tests know about config settings.

    def beforeTest(self, test):
        # We register our check function and record with it the names that
        # the test has at birth so that later we only wipe stuff created
        # during testing.
        # We really don't want a smart dict_keys object here, as it'd be
        # changed during testing.
        test.test.addCleanup(memleak_check, test.test, list(test.test.__dict__.keys()))

plugin_class = Memleak

