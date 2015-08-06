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

# initial simple tests for logging module

import os
import tempfile
import logging

from numpy.testing import TestCase

import MDAnalysis
import MDAnalysis.lib.log


class TestLogging(TestCase):
    name = "MDAnalysis"

    def setUp(self):
        fd, self.logfile = tempfile.mkstemp(suffix=".log")

    def tearDown(self):
        try:
            os.unlink(self.logfile)
        except OSError:
            pass

    def test_start_stop_logging(self):
        try:
            MDAnalysis.log.start_logging()
            logger = logging.getLogger(self.name)
            logger.info("Using the MDAnalysis logger works")
        except Exception as err:
            raise AssertionError("Problem with logger: {0}".format(err))
        finally:
            MDAnalysis.log.stop_logging()
