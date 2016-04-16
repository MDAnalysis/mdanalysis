# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- http://www.MDAnalysis.org
# Copyright (c) 2006-2015 Naveen Michaud-Agrawal, Elizabeth J. Denning, Oliver
# Beckstein and contributors (see AUTHORS for the full list)
#
# Released under the GNU Public Licence, v2 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

# initial simple tests for logging module
from six.moves import StringIO

import sys
import os
import logging

from numpy.testing import TestCase, assert_

from six.moves import range

import MDAnalysis
import MDAnalysis.lib.log

from MDAnalysisTests import tempdir


class TestLogging(TestCase):
    name = "MDAnalysis"

    def setUp(self):
        self.tempdir = tempdir.TempDir()
        self.outfile = os.path.join(self.tempdir.name, 'test.log')

    def tearDown(self):
        del self.tempdir

    def test_start_stop_logging(self):
        try:
            MDAnalysis.log.start_logging()
            logger = logging.getLogger(self.name)
            logger.info("Using the MDAnalysis logger works")
        except Exception as err:
            raise AssertionError("Problem with logger: {0}".format(err))
        finally:
            MDAnalysis.log.stop_logging()


class RedirectedStderr(object):
    """Temporarily replaces sys.stderr with *stream*.

    Deals with cached stderr, see
    http://stackoverflow.com/questions/6796492/temporarily-redirect-stdout-stderr
    """

    def __init__(self, stream=None):
        self._stderr = sys.stderr
        self.stream = stream or sys.stdout

    def __enter__(self):
        self.old_stderr = sys.stderr
        self.old_stderr.flush()
        sys.stderr = self.stream

    def __exit__(self, exc_type, exc_value, traceback):
        self._stderr.flush()
        sys.stderr = self.old_stderr

class TestProgressMeter(TestCase):
    def setUp(self):
        self.buf = StringIO()

    def tearDown(self):
        del self.buf

    def _assert_in(self, output, string):
        assert_(string in output,
                "Output '{0}' does not match required format '{1}'.".format(
                output.replace('\r', '\\r'), string))


    def test_default_ProgressMeter(self, n=101, interval=10):
        format = "Step %(step)5d/%(numsteps)d [%(percentage)5.1f%%]\r"
        with RedirectedStderr(self.buf):
            pm = MDAnalysis.lib.log.ProgressMeter(n, interval=interval)
            for frame in range(n):
                pm.echo(frame)
        self.buf.seek(0)
        output = "".join(self.buf.readlines())
        self._assert_in(output, format % {'step': 1, 'numsteps': n, 'percentage': 100./n})
        # last line always has \n instead of \r!
        self._assert_in(output, format.replace('\r', '\n') %
                        {'step': n, 'numsteps': n, 'percentage': 100.})

    def test_custom_ProgressMeter(self, n=51, interval=7):
        format = "RMSD %(rmsd)5.2f at %(step)03d/%(numsteps)4d [%(percentage)5.1f%%]\r"
        with RedirectedStderr(self.buf):
            pm = MDAnalysis.lib.log.ProgressMeter(n, interval=interval,
                                                  format=format, offset=1)
            for frame in range(n):
                rmsd = 0.02 * frame * (n+1)/float(n)  # n+1/n correction for 0-based frame vs 1-based counting
                pm.echo(frame, rmsd=rmsd)
        self.buf.seek(0)
        output = "".join(self.buf.readlines())
        self._assert_in(output, format %
                        {'rmsd': 0.0, 'step': 1, 'numsteps': n, 'percentage': 100./n})
        # last line always has \n instead of \r!
        self._assert_in(output, format.replace('\r', '\n') %
                        {'rmsd': 0.02*n, 'step': n, 'numsteps': n, 'percentage': 100.0})
