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

import numpy as np
from numpy.testing import assert_equal, assert_almost_equal
import pytest

import MDAnalysis as mda
from MDAnalysisTests.datafiles import PSF, DCD
from MDAnalysis.analysis.bat import BAT


class TestBAT(object):

    @pytest.fixture()
    def selected_residues(self):
        u = mda.Universe(PSF, DCD)
        ag = u.select_atoms("resid 5-10")
        return ag


    def test_bat(self, selected_residues):
        R = BAT(selected_residues)
        R.run()

        assert_equal(len(R.bat), 98,
                     err_msg="error: list is not length of trajectory")

        bat = R.bat[0]
        XYZ = R.Cartesian(bat)

        assert_almost_equal(XYZ, selected_residues.positions, 5,
          err_msg="error: Reconstructed Cartesian coordinates " + \
                  "don't match original")
