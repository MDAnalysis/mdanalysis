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
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from io import StringIO

import logging
import sys

import MDAnalysis
import MDAnalysis.lib.log
import pytest


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

