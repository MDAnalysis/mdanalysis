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
from MDAnalysis.lib.log import ProgressBar


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
