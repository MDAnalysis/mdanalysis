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

from numpy.testing import (
    assert_,
)

import MDAnalysis as mda
from MDAnalysis.analysis.rdf import InterRDF

from MDAnalysisTests.datafiles import two_water_gro


class TestInterRDF(object):
    def setUp(self):
        self.u = mda.Universe(two_water_gro)

    def tearDown(self):
        del self.u

    def _linear_water(self):
        """Set the positions of the two water molecules to be linear

        y =   | 1   |  2  | 3   |
        ------+-----------------|
        x = 1 | HW1 - OW - HW2  |
        ------+-----------------|
        x = 2 | HW2 - OW - HW2  |
        ------+-----------------|

        """
        # 2, 1, 3 because OW is first
        for at, (x, y) in zip(self.u.atoms,
                              zip([1] * 3 + [2] * 3, [2, 1, 3] * 2)):
            at.position = x, y, 0.0

    def _get_sels(self):
        s1 = self.u.atoms.OW
        s2 = self.u.atoms.HW1 + self.u.atoms.HW2

        return s1, s2

    def test_nbins(self):
        s1 = self.u.atoms[:3]
        s2 = self.u.atoms[3:]

        rdf = InterRDF(s1, s2, nbins=412)

        assert_(len(rdf.bins) == 412)

    def test_range(self):
        s1 = self.u.atoms[:3]
        s2 = self.u.atoms[3:]

        rmin, rmax = 1.0, 13.0

        rdf = InterRDF(s1, s2, range=(rmin, rmax))

        assert_(rdf.edges[0] == rmin)
        assert_(rdf.edges[-1] == rmax)

    def test_count_sum(self):
        # OW vs HW
        # should see 8 comparisons in count
        self._linear_water()

        s1, s2 = self._get_sels()

        rdf = InterRDF(s1, s2)
        rdf.run()

        assert_(rdf.count.sum() == 8)

    def test_count(self):
        # should see two distances with 4 counts each
        self._linear_water()

        s1, s2 = self._get_sels()

        rdf = InterRDF(s1, s2)
        rdf.run()

        assert_(len(rdf.count[rdf.count == 4]) == 2)
        
    def test_exclusion(self):
        # should see two distances with 4 counts each
        self._linear_water()

        s1, s2 = self._get_sels()

        rdf = InterRDF(s1, s2, exclusion_block=(1, 2))
        rdf.run()

        assert_(rdf.count.sum() == 4)
