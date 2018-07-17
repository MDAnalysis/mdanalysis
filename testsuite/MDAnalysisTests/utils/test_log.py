# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- https://www.mdanalysis.org
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

from six.moves import StringIO, range

import logging
import sys

import MDAnalysis
import MDAnalysis.lib.log
import pytest
from MDAnalysis.lib.log import _set_verbose


def test_start_stop_logging():
    try:
        MDAnalysis.log.start_logging()
        logger = logging.getLogger("MDAnalysis")
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


@pytest.fixture()
def buffer():
    return StringIO()


def _assert_in(output, string):
    assert string in output, "Output '{0}' does not match required format '{1}'.".format(output.replace('\r', '\\r'), string.replace('\r', '\\r'))


def test_default_ProgressMeter(buffer, n=101, interval=10):
    template = "Step {step:5d}/{numsteps} [{percentage:5.1f}%]"
    with RedirectedStderr(buffer):
        pm = MDAnalysis.lib.log.ProgressMeter(n, interval=interval)
        for frame in range(n):
            pm.echo(frame)
    buffer.seek(0)
    output = "".join(buffer.readlines())
    _assert_in(output, ('\r' + template).format(**{'step': 1, 'numsteps': n, 'percentage': 100./n}))
    # last line always ends with \n!
    _assert_in(output,
                    ('\r' + template + '\n').format(**{'step': n, 'numsteps': n,
                                              'percentage': 100.}))


def test_custom_ProgressMeter(buffer, n=51, interval=7):
    template = "RMSD {rmsd:5.2f} at {step:03d}/{numsteps:4d} [{percentage:5.1f}%]"
    with RedirectedStderr(buffer):
        pm = MDAnalysis.lib.log.ProgressMeter(n, interval=interval,
                                              format=template, offset=1)
        for frame in range(n):
            # n+1/n correction for 0-based frame vs 1-based counting
            rmsd = 0.02 * frame * (n+1)/ n
            pm.echo(frame, rmsd=rmsd)
    buffer.seek(0)
    output = "".join(buffer.readlines())
    _assert_in(output,
                    ('\r' + template).format(**{'rmsd': 0.0, 'step': 1,
                                     'numsteps': n, 'percentage': 100./n}))
    # last line always ends with \n!
    _assert_in(output,
                    ('\r' + template + '\n').format(
                        **{'rmsd': 0.02*n, 'step': n,
                           'numsteps': n, 'percentage': 100.0}))


def test_legacy_ProgressMeter(buffer, n=51, interval=7):
    template = "RMSD %(rmsd)5.2f at %(step)03d/%(numsteps)4d [%(percentage)5.1f%%]"
    with RedirectedStderr(buffer):
        pm = MDAnalysis.lib.log.ProgressMeter(n, interval=interval,
                                              format=template, offset=1)
        for frame in range(n):
            # n+1/n correction for 0-based frame vs 1-based counting
            rmsd = 0.02 * frame * (n+1)/ n
            pm.echo(frame, rmsd=rmsd)
    buffer.seek(0)
    output = "".join(buffer.readlines())
    _assert_in(output,
                    ('\r' + template) % {'rmsd': 0.0, 'step': 1,
                                       'numsteps': n, 'percentage': 100./n})
    # last line always ends with \n!
    _assert_in(output,
                    ('\r' + template + '\n') % {'rmsd': 0.02*n, 'step': n,
                                              'numsteps': n, 'percentage': 100.0})


@pytest.mark.parametrize('step, percentage', [
    (1, 100./51),
    (51, 100.)
])
def test_not_dynamic_ProgressMeter(buffer, step, percentage, n=51, interval=10):
    template = "Step {step:5d}/{numsteps} [{percentage:5.1f}%]"
    with RedirectedStderr(buffer):
        pm = MDAnalysis.lib.log.ProgressMeter(n, interval=interval,
                                              dynamic=False)
        for frame in range(n):
            pm.echo(frame)
    buffer.seek(0)
    output = "".join(buffer.readlines())
    _assert_in(output, (template + '\n').format(**{'step': step, 'numsteps': n, 'percentage': percentage}))


class TestSetVerbose(object):

    @pytest.mark.parametrize('verbose, quiet, default, result', [
        (True, False, True, True),  # Everything agrees verbose should be True
        (False, True, False, False),# Everything agrees verbose should be False
        (True, False, False, True), # Make sure the default does not overwrite the user choice
        (False, True, True, False), # Make sure the default does not overwrite the user choice
        (None, True, False, False), # Verbose is not provided
        (None, False, False, True), # Verbose is not provided
    ])
    def test__set_verbose_deprecated(self, verbose, quiet, default, result):
        with pytest.deprecated_call():
            assert _set_verbose(verbose=verbose, quiet=quiet, default=default) == result

    @pytest.mark.parametrize('verbose, quiet, default, result', [
        (True, None, False, True),  # Quiet is not provided
        (False, None, False, False),# Quiet is not provided
        (None, None, True, True),   # Nothing is provided
        (None, None, False, False), # Nothing is provided
    ])
    def test__set_verbose(self, verbose, quiet, default, result):
        assert _set_verbose(verbose=verbose, quiet=quiet, default=default) == result

    @pytest.mark.parametrize('verbose', (True, False))
    def test__set_verbose_invalid_args(self, verbose):
        # can't combine the two context managers
        with pytest.deprecated_call():
            with pytest.raises(ValueError):
                # setting quiet==verbose is a contradiction
                _set_verbose(verbose=verbose, quiet=verbose, default=None)


    @pytest.mark.parametrize('verbose, quiet', [
        (None, True),
        (False, True)
    ])
    def test_warnings__set_verbose(self, verbose, quiet):
        pytest.deprecated_call(_set_verbose, verbose=verbose, quiet=quiet)
