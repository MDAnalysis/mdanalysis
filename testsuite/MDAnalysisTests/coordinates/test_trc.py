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
from numpy.testing import assert_almost_equal

import MDAnalysis as mda
from MDAnalysis.coordinates.TRC import TRCReader
from MDAnalysisTests.datafiles import TRC_PDB, TRC_TRAJ1, TRC_TRAJ2


class TestTRCReader:
    @pytest.fixture(scope='class')
    def TRC_U(self):
        return mda.Universe(TRC_PDB, [TRC_TRAJ1, TRC_TRAJ2], 
                            format="TRC", continuous=True)
    
    def test_trc_positions(self, TRC_U):
        # first frame first particle
        TRC_U.trajectory[0]
        assert_almost_equal(TRC_U.atoms.positions[0],
                            [2.19782507, 24.65064345, 29.39783426])
        # second frame first particle
        TRC_U.trajectory[5]
        assert_almost_equal(TRC_U.atoms.positions[0],
                            [0.37026654, 22.78805010, 3.69695262])
    
    def test_trc_dimensions(self, TRC_U):
        ts = TRC_U.trajectory[0]
        assert_almost_equal(
            ts.dimensions,
            [30.70196350, 30.70196350, 30.70196350, 90., 90., 90.]
        )

    def test_trc_n_frames(self, TRC_U):
        assert len(TRC_U.trajectory) == 6

    def test_trc_frame(self, TRC_U):
        assert TRC_U.trajectory[0].frame == 0
        assert TRC_U.trajectory[4].frame == 4

    def test_trc_time(self, TRC_U):
        assert TRC_U.trajectory[0].time == 0
        assert TRC_U.trajectory[4].time == 80

    def test_trc_data_step(self, TRC_U):
        assert TRC_U.trajectory[0].data['step'] == 0
        assert TRC_U.trajectory[4].data['step'] == 10000
