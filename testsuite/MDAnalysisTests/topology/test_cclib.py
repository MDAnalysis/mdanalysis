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

import MDAnalysis as mda
import pytest
import numpy as np
from numpy.testing import assert_allclose

from MDAnalysisTests.datafiles import (
    QCHEM,
)


class TestQChem(object):
    @pytest.fixture(scope='class')
    def u(self):
        return mda.Universe(QCHEM, format='QCHEM')

    @pytest.mark.parametrize('attr', ['mulliken_charges', 'masses',
                                      'core_electrons', 'names'])
    def test_expected_attrs(self, u, attr):
        assert hasattr(u.atoms, attr)

    def test_size(self, u):
        assert len(u.atoms) == 20
        assert len(u.residues) == 1
        assert len(u.segments) == 1

    def test_names(self, u):
        refnames = ['C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'C', 'C',
                   'H', 'H', 'C', 'C', 'H', 'H', 'H', 'H']
        for atom, refname in zip(u.atoms, refnames):
            assert atom.name == refname

    def test_charges(self, u):
        refcharges = np.array([-0.004404, -0.004405, -0.076752, -0.076752,
                               -0.077083, -0.077083, 0.077995, 0.077995,
                               0.079316, 0.079316, -0.076366, -0.076366,
                               0.076693, 0.076693, -0.154814, -0.154814,
                               0.076117, 0.076117, 0.079299, 0.079299])
        for atom, refcharge in zip(u.atoms, refcharges):
            assert atom.mulliken_charge == refcharge

    def test_trajectory(self, u):
        assert len(u.trajectory) == 1

    def test_positions(self, u):
        refpos = np.array([[ 1.41308741,  0.24316249,  0.        ],
                           [-1.41308741, -0.24316249,  0.        ],
                           [ 0.48308439,  1.31893657,  0.        ],
                           [-0.48308439, -1.31893657,  0.        ],
                           [-0.89935459,  1.08276512,  0.        ],
                           [ 0.89935459, -1.08276512,  0.        ],
                           [ 0.85806792,  2.35168462,  0.        ],
                           [-0.85806792, -2.35168462,  0.        ],
                           [-1.59563597,  1.93048601,  0.        ],
                           [ 1.59563597, -1.93048601,  0.        ],
                           [-2.87969547, -0.54746309,  0.        ],
                           [ 2.87969547,  0.54746309,  0.        ],
                           [-3.12563689, -1.62062971,  0.        ],
                           [ 3.12563689,  1.62062971,  0.        ],
                           [-3.88337819,  0.34673302,  0.        ],
                           [ 3.88337819, -0.34673302,  0.        ],
                           [-3.70790402,  1.42856941,  0.        ],
                           [ 3.70790402, -1.42856941,  0.        ],
                           [-4.93208953,  0.02597472,  0.        ],
                           [ 4.93208953, -0.02597472,  0.        ]])
        assert_allclose(u.atoms.positions, refpos)
