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
from __future__ import division, absolute_import
# initial simple tests for logging module
from six.moves import StringIO

import sys
import os
import logging
import warnings

from numpy.testing import (TestCase, assert_, assert_equal,
                           assert_raises, assert_warns)

from six.moves import range

import MDAnalysis
import MDAnalysis.lib.log
from MDAnalysis.lib.log import _set_verbose

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
                output.replace('\r', '\\r'), string.replace('\r', '\\r')))

    def test_default_ProgressMeter(self, n=101, interval=10):
        format = "Step {step:5d}/{numsteps} [{percentage:5.1f}%]"
        with RedirectedStderr(self.buf):
            pm = MDAnalysis.lib.log.ProgressMeter(n, interval=interval)
            for frame in range(n):
                pm.echo(frame)
        self.buf.seek(0)
        output = "".join(self.buf.readlines())
        self._assert_in(output, ('\r' + format).format(**{'step': 1, 'numsteps': n, 'percentage': 100./n}))
        # last line always ends with \n!
        self._assert_in(output,
                        ('\r' + format + '\n').format(**{'step': n, 'numsteps': n,
                                                  'percentage': 100.}))

    def test_custom_ProgressMeter(self, n=51, interval=7):
        format = "RMSD {rmsd:5.2f} at {step:03d}/{numsteps:4d} [{percentage:5.1f}%]"
        with RedirectedStderr(self.buf):
            pm = MDAnalysis.lib.log.ProgressMeter(n, interval=interval,
                                                  format=format, offset=1)
            for frame in range(n):
                rmsd = 0.02 * frame * (n+1)/float(n)  # n+1/n correction for 0-based frame vs 1-based counting
                pm.echo(frame, rmsd=rmsd)
        self.buf.seek(0)
        output = "".join(self.buf.readlines())
        self._assert_in(output,
                        ('\r' + format).format(**{'rmsd': 0.0, 'step': 1,
                                         'numsteps': n, 'percentage': 100./n}))
        # last line always ends with \n!
        self._assert_in(output,
                        ('\r' + format + '\n').format(
                            **{'rmsd': 0.02*n, 'step': n,
                               'numsteps': n, 'percentage': 100.0}))

    def test_legacy_ProgressMeter(self, n=51, interval=7):
        format = "RMSD %(rmsd)5.2f at %(step)03d/%(numsteps)4d [%(percentage)5.1f%%]"
        with RedirectedStderr(self.buf):
            pm = MDAnalysis.lib.log.ProgressMeter(n, interval=interval,
                                                  format=format, offset=1)
            for frame in range(n):
                rmsd = 0.02 * frame * (n+1)/float(n)  # n+1/n correction for 0-based frame vs 1-based counting
                pm.echo(frame, rmsd=rmsd)
        self.buf.seek(0)
        output = "".join(self.buf.readlines())
        self._assert_in(output,
                        ('\r' + format) % {'rmsd': 0.0, 'step': 1,
                                           'numsteps': n, 'percentage': 100./n})
        # last line always ends with \n!
        self._assert_in(output,
                        ('\r' + format + '\n') % {'rmsd': 0.02*n, 'step': n,
                                                  'numsteps': n, 'percentage': 100.0})

    def test_not_dynamic_ProgressMeter(self, n=51, interval=10):
        format = "Step {step:5d}/{numsteps} [{percentage:5.1f}%]"
        with RedirectedStderr(self.buf):
            pm = MDAnalysis.lib.log.ProgressMeter(n, interval=interval,
                                                  dynamic=False)
            for frame in range(n):
                pm.echo(frame)
        self.buf.seek(0)
        output = "".join(self.buf.readlines())
        self._assert_in(output, (format + '\n').format(**{'step': 1, 'numsteps': n, 'percentage': 100./n}))
        self._assert_in(output, (format + '\n').format(**{'step': n, 'numsteps': n, 'percentage': 100.}))


def test__set_verbose():
    # Everything agrees verbose should be True
    assert_equal(_set_verbose(verbose=True, quiet=False, default=True), True)
    # Everything agrees verbose should be False
    assert_equal(_set_verbose(verbose=False, quiet=True, default=False), False)
    # Make sure the default does not overwrite the user choice
    assert_equal(_set_verbose(verbose=True, quiet=False, default=False), True)
    assert_equal(_set_verbose(verbose=False, quiet=True, default=True), False)
    # Quiet is not provided
    assert_equal(_set_verbose(verbose=True, quiet=None, default=False), True)
    assert_equal(_set_verbose(verbose=False, quiet=None, default=False), False)
    # Verbose is not provided
    assert_equal(_set_verbose(verbose=None, quiet=True, default=False), False)
    assert_equal(_set_verbose(verbose=None, quiet=False, default=False), True)
    # Nothing is provided
    assert_equal(_set_verbose(verbose=None, quiet=None, default=True), True)
    assert_equal(_set_verbose(verbose=None, quiet=None, default=False), False)
    # quiet and verbose contradict each other
    assert_raises(ValueError, _set_verbose, verbose=True, quiet=True)
    assert_raises(ValueError, _set_verbose, verbose=False, quiet=False)
    # A deprecation warning is issued when quiet is set

    # The following tests are commented out because they fail only when the file `test_log.py`
    # is run individually. Initially seen in #1370

    # assert_warns(DeprecationWarning, _set_verbose, verbose=None, quiet=True)
    # assert_warns(DeprecationWarning, _set_verbose, verbose=False, quiet=True)
