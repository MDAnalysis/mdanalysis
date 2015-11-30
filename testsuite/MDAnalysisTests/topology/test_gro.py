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
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#
from numpy.testing import (
    assert_,
)

import MDAnalysis as mda

from MDAnalysisTests.datafiles import (
    two_water_gro_widebox,
)



class TestGROWideBox(object):
    """Tests for Issue #548"""
    def test_atoms(self):
        parser = mda.topology.GROParser.GROParser
        with parser(two_water_gro_widebox) as p:
            s = p.parse()
        assert_(len(s['atoms']) == 6)
