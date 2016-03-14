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

from numpy.testing import assert_array_equal


class TestDatafiles(object):
    def test_import(self):
        try:
            import MDAnalysis.tests.datafiles
        except ImportError:
            raise AssertionError("Failed to 'import MDAnalysis.tests.datafiles --- install MDAnalysisTests")

    def test_all_exports(self):
        import MDAnalysisTests.datafiles
        missing = [name for name in dir(MDAnalysisTests.datafiles)
                   if not name.startswith('_') and name not in MDAnalysisTests.datafiles.__all__]
        assert_array_equal(missing, [], err_msg="Variables need to be added to __all__.")

    def test_export_variables(self):
        import MDAnalysisTests.datafiles
        import MDAnalysis.tests.datafiles
        missing = [name for name in MDAnalysisTests.datafiles.__all__
                   if name not in dir(MDAnalysis.tests.datafiles)]
        assert_array_equal(missing, [], err_msg="Variables not exported to MDAnalysis.tests.datafiles")
