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

"""Test plugin that reports memleaks (objects deemed uncollectable by python's garbage collector) on a test-by-test basis.

The plugin works by clearing a test's namespace after it has run, and then
forcing a garbage collection round and checking for uncollectable objects.

Implementation uses the startTest and stopTest hooks. To prevent premature
success declaration (since the memleak check only runs after the test) the test's
:class:`nose.proxy.ResultProxy` is monkey-patched with a silent :meth:`addSuccess` method. This
might seem hackish, but it is how other `:class:nose.plugins.errorclass.ErrorClassPlugins`
are implemented.
"""

import gc
import nose
from nose.plugins.errorclass import ErrorClass, ErrorClassPlugin
from nose.pyversion import make_instancemethod

_leakedobjs = set() # We must keep track of seen leaks to avoid raising multiple 
                    # errors for the same ones.

class MemleakError(Exception):
    """Exception raised internally to mark the test as memleaking."""
    pass

class _MemleakResult(nose.proxy.ResultProxy):
    """Placeholder class where a replacement :meth:`addSuccess` function is defined.

    The :meth:`addSuccess` method will be later monkey-patched onto result objects, saving
    success declarations for later.
    """
    def addSuccess(self, test):
        self._ml_is_success = True

class Memleak(ErrorClassPlugin):
    """
    Report memleaks (objects deemed uncollectable by python's
    garbage collector) on a test-by-test basis.
    """
    name = "memleak"
    enabled = False
    memleak = ErrorClass(MemleakError,
                         label='MEMLEAK',
                         isfailure=True)

    def configure(self, options, conf):
        super(Memleak, self).configure(options, conf)
        self.config = conf # This will let other tests know about config settings.

    def startTest(self, test):
        # we want to change the ResultProxy object, not just the result,
        #  therefore we do it here, not in prepareTestResult
        rp = test.test._resultForDoCleanups
        # we store the original func
        rp._ml_add_success = rp.addSuccess
        # and replace it
        rp.addSuccess = make_instancemethod(_MemleakResult.addSuccess, rp)
        # also add a default failure flag
        rp._ml_is_success = False
        # Finally, we store with the test the names that it has at birth,
        #  so that we only wipe stuff created during the test.
        test._mda_pretest_attrs = test.test.__dict__.keys()

    def stopTest(self, test):
        rp = test.test._resultForDoCleanups
        # cleanup of the patched attributes as soon as possible, to minimize potential clashes
        rp.addSuccess = rp._ml_add_success
        del rp._ml_add_success
        attrs = []
        for attrname in test._mda_pretest_attrs:
            try:
                attrs.append((attrname, getattr(test.test, attrname)))
            except AttributeError:
                pass
        test.test.__dict__.clear()
        gc.collect()
        latest_leaks = [ obj for obj in gc.garbage if obj not in _leakedobjs ]
        _leakedobjs.update(gc.garbage)
        # Restore the pre-test stuff
        for name, val in attrs:
            setattr(test.test, name, val)
        if latest_leaks:
            raise MemleakError("GC failed to collect the following: {}".format(latest_leaks))
        elif rp._ml_is_success:
            rp.addSuccess(test)
            del rp._ml_is_success

plugin_class = Memleak

