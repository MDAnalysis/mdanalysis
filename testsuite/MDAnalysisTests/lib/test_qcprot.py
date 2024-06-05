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
import numpy as np
from numpy.testing import assert_allclose
import pytest

import MDAnalysis as mda
from MDAnalysis.lib import qcprot

from MDAnalysisTests.datafiles import PSF, DCD


@pytest.mark.parametrize('dtype', [np.float64, np.float32])
class TestQCProt:
    def test_dummy(self, dtype):
        a = np.array([[1.0, 1.0, 2.0]], dtype=dtype)
        b = np.array([[1.0, 1.0, 1.0]], dtype=dtype)
        rot = np.zeros(9, dtype=dtype)

        r = qcprot.CalcRMSDRotationalMatrix(a, b, 1, rot, None)

        err = 0.0001 if dtype is np.float32 else 0.000001

        assert r == pytest.approx(0.7174389197644676, abs=err)
        assert rot.dtype == dtype
        assert_allclose(rot.reshape(3, 3), np.eye(3))

    def test_regression(self, dtype):
        u = mda.Universe(PSF, DCD)
        prot = u.select_atoms('protein')
        weights = prot.masses.astype(dtype)
        weights /= np.mean(weights)
        p1 = prot.positions.astype(dtype)
        u.trajectory[1]
        p2 = prot.positions.astype(dtype)
        rot = np.zeros(9, dtype=dtype)

        r = qcprot.CalcRMSDRotationalMatrix(p1, p2, len(prot), rot, weights)

        rot_ref = np.array([0.999998, 0.001696, 0.001004,
                            -0.001698,  0.999998,  0.001373,
                            -0.001002, -0.001375,  0.999999],
                           dtype=dtype)

        err = 0.001 if dtype is np.float32 else 0.000001
        assert r == pytest.approx(0.6057544485785074, abs=err)
        assert_allclose(rot, rot_ref, atol=err)
