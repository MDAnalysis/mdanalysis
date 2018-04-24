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
from __future__ import absolute_import
import os
import pytest

import importlib

# environment variable DUECREDIT_ENABLE is set to yes in MDAnalysisTests/__init__.py
# (if set to 'no', the tests will be SKIPPED; has to be yes, true, or 1 for duecredit
# to work; duecredit must also be installed)
import MDAnalysis as mda


@pytest.mark.skipif((os.environ.get('DUECREDIT_ENABLE', 'yes').lower()
                     in ('no', '0', 'false')),
                    reason=
                    "duecredit is explicitly disabled with DUECREDIT_ENABLE=no")
class TestDuecredits(object):

    def test_duecredit_active(self):
        assert mda.due.active == True

    @pytest.mark.parametrize("module,path,citekey", [
        ("MDAnalysis", "MDAnalysis", "gowers2016"),
        ("MDAnalysis", "MDAnalysis", "10.1002/jcc.21787"),
    ])
    def test_duecredit_collector_primary(self, module, path, citekey):
        assert mda.due.citations[(path, citekey)].cites_module == True

    @pytest.mark.parametrize("module,path,citekey", [
        ("MDAnalysis.analysis.psa",
         "MDAnalysis.analysis.psa",
         "10.1371/journal.pcbi.1004568"),
        ("MDAnalysis.analysis.hbonds.hbond_autocorrel",
         "MDAnalysis.analysis.hbonds.hbond_autocorrel",
         "10.1063/1.4922445"),
        ])
    def test_duecredit_collector_analysis(self, module, path, citekey):
        importlib.import_module(module)
        assert mda.due.citations[(path, citekey)].cites_module == True
