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
import os
import pytest

import importlib

# environment variable DUECREDIT_ENABLE is set to yes in MDAnalysisTests/__init__.py
# (if set to 'no', the tests will be SKIPPED; has to be yes, true, or 1 for duecredit
# to work; duecredit must also be installed)
import MDAnalysis as mda
from MDAnalysisTests.datafiles import MMTF

# duecredit itself is not needed in the name space but this is a
# convenient way to skip all tests if duecredit is not installed
# (see https://github.com/MDAnalysis/mdanalysis/issues/1906)
pytest.importorskip('duecredit')

@pytest.mark.skipif((os.environ.get('DUECREDIT_ENABLE', 'yes').lower()
                     in ('no', '0', 'false')),
                    reason=
                    "duecredit is explicitly disabled with DUECREDIT_ENABLE=no")
class TestDuecredit(object):

    def test_duecredit_active(self):
        assert mda.due.active == True

    @pytest.mark.parametrize("module,path,citekey", [
        ("MDAnalysis", "MDAnalysis", "10.25080/majora-629e541a-00e"),
        ("MDAnalysis", "MDAnalysis", "10.1002/jcc.21787"),
    ])
    def test_duecredit_collector_primary(self, module, path, citekey):
        assert mda.due.citations[(path, citekey)].cites_module == True

    # note: citekeys are *all lower case*
    @pytest.mark.parametrize("module,path,citekey", [
        ("MDAnalysis.analysis.psa",
         "MDAnalysis.analysis.psa",
         "10.1371/journal.pcbi.1004568"),
        ("MDAnalysis.analysis.hbonds.hbond_autocorrel",
         "MDAnalysis.analysis.hbonds.hbond_autocorrel",
         "10.1063/1.4922445"),
        ("MDAnalysis.analysis.leaflet",
         "MDAnalysis.analysis.leaflet",
         "10.1002/jcc.21787"),
        ("MDAnalysis.analysis.hole2",
         "MDAnalysis.analysis.hole2",
         "10.1016/s0006-3495(93)81293-1"),
        ("MDAnalysis.analysis.hole2",
         "MDAnalysis.analysis.hole2",
         "10.1016/s0263-7855(97)00009-x"),
        ("MDAnalysis.analysis.hole2",
         "MDAnalysis.analysis.hole2",
         "10.1016/j.jmb.2013.10.024"),
        ("MDAnalysis.lib.qcprot",
         "MDAnalysis.lib.qcprot",
         "10.1107/s0108767305015266"),
        ("MDAnalysis.lib.qcprot",
         "MDAnalysis.lib.qcprot",
         "qcprot2"),
        ("MDAnalysis.analysis.encore",
         "MDAnalysis.analysis.encore",
         "10.1371/journal.pcbi.1004415"),
        ])
    def test_duecredit_collector_analysis_modules(self, module, path, citekey):
        importlib.import_module(module)
        assert mda.due.citations[(path, citekey)].cites_module == True

    def test_duecredit_mmtf(self):
        # doesn't trigger on import but on use of either parser or reader
        u = mda.Universe(MMTF)

        assert mda.due.citations[('MDAnalysis.coordinates.MMTF',
                                  '10.1371/journal.pcbi.1005575')].cites_module
        assert mda.due.citations[('MDAnalysis.topology.MMTFParser',
                                  '10.1371/journal.pcbi.1005575')].cites_module
