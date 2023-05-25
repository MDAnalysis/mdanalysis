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
import pytest
import sys
from numpy.testing import assert_equal


def test_failed_import(monkeypatch):
    # Putting this test first to avoid datafiles already being loaded
    errmsg = "MDAnalysisTests package not installed."

    monkeypatch.setitem(sys.modules, 'MDAnalysisTests.datafiles', None)

    if 'MDAnalysis.tests.datafiles' in sys.modules:
        monkeypatch.delitem(sys.modules, 'MDAnalysis.tests.datafiles')

    with pytest.raises(ImportError, match=errmsg):
        import MDAnalysis.tests.datafiles


def test_import():
    try:
        import MDAnalysis.tests.datafiles
    except ImportError:
        pytest.fail("Failed to 'import MDAnalysis.tests.datafiles --- install MDAnalysisTests")


def test_all_exports():
    import MDAnalysisTests.datafiles
    missing = [name for name in dir(MDAnalysisTests.datafiles)
               if
               not name.startswith('_') and name not in MDAnalysisTests.datafiles.__all__ and name != 'MDAnalysisTests']
    assert_equal(missing, [], err_msg="Variables need to be added to __all__.")


def test_export_variables():
    import MDAnalysisTests.datafiles
    import MDAnalysis.tests.datafiles
    missing = [name for name in MDAnalysisTests.datafiles.__all__
               if name not in dir(MDAnalysis.tests.datafiles)]
    assert_equal(missing, [], err_msg="Variables not exported to MDAnalysis.tests.datafiles")
