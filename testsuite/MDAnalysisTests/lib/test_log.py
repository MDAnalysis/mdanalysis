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

from __future__ import absolute_import
import warnings
import pytest

from MDAnalysis.lib.log import ProgressMeter, ProgressBar, echo


class TestProgressMeter(object):

    def test_deprecated(self, capsys):
        with warnings.catch_warnings(record=True) as w:
            # Cause all warnings to always be triggered.
            warnings.simplefilter("always")
            # Trigger a warning.
            pm = ProgressMeter(10)
            # Verify the warning
            assert len(w) == 1
            assert issubclass(w[-1].category, DeprecationWarning)
            assert "MDAnalysis.lib.log.ProgressBar" in str(w[-1].message)

    def test_output(self, capsys):
        pm = ProgressMeter(10, interval=1)
        for i in range(10):
            pm.echo(i)
        out, err = capsys.readouterr()
        expected = 'Step    10/10 [100.0%]'
        actual = err.strip().split('\r')[-1]
        assert actual == expected


def test_echo_deprecated():
    with pytest.deprecated_call():
        echo("F=ma -- Isaac Newton")

class TestProgressBar(object):

    def test_output(self, capsys):
        for i in ProgressBar(list(range(10))):
            pass
        out, err = capsys.readouterr()
        expected = u'100%|██████████| 10/10 [00:00<00:00, 583.67it/s]'
        actual = err.strip().split('\r')[-1]
        assert actual[:24] == expected[:24]

    def test_disable(self, capsys):
        for i in ProgressBar(list(range(10)), disable=True):
            pass
        out, err = capsys.readouterr()
        expected = ''
        actual = err.strip().split('\r')[-1]
        assert actual == expected

    def test_verbose_disable(self, capsys):
        for i in ProgressBar(list(range(10)), verbose=False):
            pass
        out, err = capsys.readouterr()
        expected = ''
        actual = err.strip().split('\r')[-1]
        assert actual == expected
