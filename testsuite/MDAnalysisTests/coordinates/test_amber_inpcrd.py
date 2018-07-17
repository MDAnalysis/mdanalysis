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

import numpy as np
from six.moves import zip

from numpy.testing import (assert_allclose, assert_equal)

import MDAnalysis as mda
from MDAnalysisTests.datafiles import (INPCRD, XYZ_five)


class TestINPCRDReader(object):
    """Test reading Amber restart coordinate files"""

    @staticmethod
    def _check_ts(ts):
        # Check a ts has the right values in
        ref_pos = np.array([[6.6528795, 6.6711416, -8.5963255],
                            [7.3133773, 5.8359736, -8.8294175],
                            [8.3254058, 6.2227613, -8.7098593],
                            [7.0833200, 5.5038197, -9.8417650],
                            [7.1129439, 4.6170351, -7.9729560]])
        for ref, val in zip(ref_pos, ts._pos):
            assert_allclose(ref, val)

    def test_reader(self):
        from MDAnalysis.coordinates.INPCRD import INPReader
        r = INPReader(INPCRD)
        assert_equal(r.n_atoms, 5)
        self._check_ts(r.ts)

    def test_universe_inpcrd(self):
        u = mda.Universe(XYZ_five, INPCRD)
        self._check_ts(u.trajectory.ts)

    def test_universe_restrt(self):
        u = mda.Universe(XYZ_five, INPCRD, format='RESTRT')
        self._check_ts(u.trajectory.ts)
