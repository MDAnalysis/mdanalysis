# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4 fileencoding=utf-8
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
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

import pytest
import sys
import logging
import MDAnalysis as mda

from os.path import basename
from MDAnalysis.lib.log import ProgressBar


class TestConvenienceFunctions:
    def test_start_logging(self, tmp_path):
        mda.start_logging(tmp_path / "MDAnalysis.log")
        logger = logging.getLogger("MDAnalysis")

        # Test expected handlers' presence and behavior
        assert any(isinstance(h, logging.NullHandler) for h in logger.handlers)
        any(
            isinstance(h, logging.FileHandler)
            and basename(h.stream.name) == "MDAnalysis.log"
            and h.level == logging.DEBUG
            for h in logger.handlers
        )
        assert any(
            isinstance(h, logging.StreamHandler)
            and h.stream is sys.stdout
            and h.level == logging.INFO
            for h in logger.handlers
        )

    def test_stop_logging(self, tmp_path):
        mda.lib.log.start_logging(tmp_path / "MDAnalysis.log")
        logger = logging.getLogger("MDAnalysis")
        mda.lib.log.stop_logging()

        assert len(logger.handlers) == 0

def test_message_console(tmp_path):
    pass

def test_message_file(tmp_path):
    pass

# TODO need to make a fixture that can clear all handlers per test
class TestCreateBehaviors:

    def test_input_path(self, tmp_path):
        mda.lib.log.create(stream=tmp_path / "foo.log")

        assert (tmp_path / "foo.log").exists()

    def test_input_string(self, tmp_path):
        mda.lib.log.create(stream=str(tmp_path / "foo.log"))

        assert (tmp_path / "foo.log").exists()

    # NOTE Assert state could be cleaned up after clear_handlers() fixture is implemented
    def test_input_stream(self):
        mda.lib.log.create(stream=sys.stdout)

        logger = logging.getLogger("MDAnalysis")
        assert any(
            isinstance(h, logging.StreamHandler) and h.stream is sys.stdout
            for h in logger.handlers
        )

    def test_exception(tmp_path):
        with pytest.raises(
            TypeError,
            match="Input Stream is neither a string, PathLike object or a stream",
        ):
            mda.lib.log.create(stream=2)

def test_level_parameter():
    pass

def test_fmt_parameter():
    pass

def test_mode_parameter():
    pass

class TestProgressBar(object):

    def test_output(self, capsys):
        for i in ProgressBar(list(range(10))):
            pass
        out, err = capsys.readouterr()
        expected = "100%|██████████"
        actual = err.strip().split("\r")[-1]
        assert actual[:15] == expected

    def test_disable(self, capsys):
        for i in ProgressBar(list(range(10)), disable=True):
            pass
        out, err = capsys.readouterr()
        expected = ""
        actual = err.strip().split("\r")[-1]
        assert actual == expected

    def test_verbose_disable(self, capsys):
        for i in ProgressBar(list(range(10)), verbose=False):
            pass
        out, err = capsys.readouterr()
        expected = ""
        actual = err.strip().split("\r")[-1]
        assert actual == expected
