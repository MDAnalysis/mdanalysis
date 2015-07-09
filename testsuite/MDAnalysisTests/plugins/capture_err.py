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
# This plugin was adapted from nose's capture.py

"""
This plugin captures stderr during test execution. If the test fails
or raises an error, the captured output will be appended to the error
or failure output. It is enabled by default but can be disabled with
the option `--no-errorcapture`.
"""

import os
import sys
from nose.plugins.base import Plugin
from nose.pyversion import exc_to_unicode, force_unicode
from nose.util import ln
from StringIO import StringIO

class CaptureStdErr(Plugin):
    """
    Output capture plugin. Enabled by default. Disable with
    ``--no-errorcapture``. This plugin captures stderr during test execution,
    appending any output captured to the error or failure output,
    should the test fail or raise an error.
    """
    enabled = True
    env_opt = 'NOSE_NO_ERRORCAPTURE'
    name = 'capture-error'
    score = 1600

    def __init__(self):
        self.stderr = []
        self._buf = None

    def options(self, parser, env):
        """Register commandline options
        """
        parser.add_option(
            "--no-errorcapture", action="store_false",
            default=not env.get(self.env_opt), dest="capture_error",
            help="Don't capture stderr (any stderr output "
            "will be printed immediately) [NOSE_NO_ERRORCAPTURE]")

    def configure(self, options, conf):
        """Configure plugin. Plugin is enabled by default.
        """
        self.config = conf
        try:
            self.enabled = options.capture_error
        except AttributeError:
            self.enabled = False

    def afterTest(self, test):
        """Clear capture buffer.
        """
        self.end()
        self._buf = None

    def begin(self):
        """Replace sys.stderr with capture buffer.
        """
        self.start() # get an early handle on sys.stdout

    def beforeTest(self, test):
        """Flush capture buffer.
        """
        self.start()

    def formatError(self, test, err):
        """Add captured output to error report.
        """
        test.capturedOutput = output = self.buffer
        self._buf = None
        if not output:
            # Don't return None as that will prevent other
            # formatters from formatting and remove earlier formatters
            # formats, instead return the err we got
            return err
        ec, ev, tb = err
        return (ec, self.addCaptureToErr(ev, output), tb)

    def formatFailure(self, test, err):
        """Add captured output to failure report.
        """
        return self.formatError(test, err)

    def addCaptureToErr(self, ev, output):
        ev = exc_to_unicode(ev)
        output = force_unicode(output)
        return u'\n'.join([ev, ln(u'>> begin captured stderr <<'),
                           output, ln(u'>> end captured stderr <<')])

    def start(self):
        self.stderr.append(sys.stderr)
        self._buf = StringIO()
        sys.stderr = self._buf

    def end(self):
        if self.stderr:
            sys.stderr = self.stderr.pop()

    def finalize(self, result):
        """Restore stderr.
        """
        while self.stderr:
            self.end()

    def _get_buffer(self):
        if self._buf is not None:
            return self._buf.getvalue()

    buffer = property(_get_buffer, None, None,
                      """Captured stderr output.""")

plugin_class = CaptureStdErr
