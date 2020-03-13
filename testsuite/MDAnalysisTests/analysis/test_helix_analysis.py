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
import re

import numpy as np
import pytest
from numpy.testing import assert_equal, assert_almost_equal

import MDAnalysis as mda
from MDAnalysis.analysis import helix_analysis as hel
from MDAnalysisTests.datafiles import (GRO, XTC, PSF, DCD, PDB_small,
                                       HELANAL_BENDING_MATRIX,
                                       HELANAL_BENDING_MATRIX_SUBSET,
                                       XYZ)

# reference data from a single PDB file:
#   data = MDAnalysis.analysis.helanal.helanal_main(PDB_small,
#                    select="name CA and resnum 161-187")
# keys renamed
HELANAL_SINGLE_DATA = {
    'global_tilts': np.rad2deg(1.3309656332019535),
    'local_heights summary': np.array([1.5286051,  0.19648294,  0.11384312],
                                      dtype=np.float32),
    'local_bends':
        np.array([3.44526005,   4.85425806,   4.69548464,   2.39473653,
                  3.56172442,   3.97128344,   3.41563916,   1.86140978,
                  5.22997046,   5.41381264,  27.49601364,  39.69839478,
                  35.05921936,  21.78928566,   9.8632431,   8.80066967,
                  5.5344553,   6.14356709,  10.15450764,  11.07686138,
                  9.23541832], dtype=np.float32),
    'residues_per_turn summary': np.array([3.64864163,  0.152694,  0.1131402]),
    'local_screw':
        np.array([87.80540079, -171.86019984,  -75.2341296,   24.61695962,
                  121.77104796, -134.94786976,  -35.07857424,   58.9621866,
                  159.40210233, -104.83368122,   -7.54816243,   87.40202629,
                  -176.13071955,  -89.13196878,   17.17321345,  105.26627814,
                  -147.00075298,  -49.36850775,   54.24557615,  156.06486532,
                  -110.82698327,   -5.72138626,   85.36050546, -167.28218858,
                  -68.23076936]),
    'local_twists summary': np.array([98.83011627,   4.08171701,   3.07253003],
                                     dtype=np.float32),
    'local_twists':
        np.array([97.23709869,   99.09676361,   97.25350952,  101.76019287,
                  100.42689514,   97.08784485,   97.07430267,   98.33553314,
                  97.86578369,   95.45792389,   97.10089111,   95.26415253,
                  87.93136597,  108.38458252,   95.27167511,  104.01931763,
                  100.7199707,  101.48034668,   99.64170074,   94.78946686,
                  102.50147247,   97.25154877,  104.54204559,  101.42829895],
                 dtype=np.float32),
}


def read_bending_matrix(fn):
    """Read helanal_bending_matrix.dat into dict of numpy arrays.

    This is a quick and dirty parser with no error checking.

    Format::

      Mean
         0.0    5.7   10.9 ....
         .

      SD
        ....
        .

      ABDEV
       ....
       .

    Returns
    -------
       data : dict
           dictionary with keys 'Mean', 'SD', 'ABDEV' and NxN matrices.
    """
    data = {}
    current = None
    with open(fn) as bendingmatrix:
        for line in bendingmatrix:
            line = line.strip()
            if not line:
                continue
            if re.match("\D+", line):
                label = line.split()[0]
                current = data[label] = []
            else:
                current.append([float(x) for x in line.split()])
    for k, v in data.items():
        data[k] = np.array(v)
    return data


class TestHELANAL(object):

    @pytest.fixture()
    def universe(self):
        return mda.Universe(PSF, DCD)

    @pytest.fixture()
    def helanal(self, universe):
        ha = hel.HELANAL(universe, select='name CA and resnum 161-187',
                         flatten_single_helix=True)
        return ha.run(start=10, stop=80)

    def test_summary(self, helanal):
        bends = helanal.summary['all_bends']
        old_helanal = read_bending_matrix(HELANAL_BENDING_MATRIX_SUBSET)
        assert_almost_equal(np.triu(bends['mean'], 1), old_helanal['Mean'],
                            decimal=1)
        assert_almost_equal(np.triu(bends['sample_sd'], 1), old_helanal['SD'],
                            decimal=1)
        assert_almost_equal(np.triu(bends['abs_dev'], 1), old_helanal['ABDEV'],
                            decimal=1)

    def test_correct_values(self):
        u = mda.Universe(PDB_small)
        ha = hel.HELANAL(u, select='name CA and resnum 161-187',
                         flatten_single_helix=True)
        ha.run()

        for key, value in HELANAL_SINGLE_DATA.items():
            if 'summary' in key:
                data = getattr(ha, key.split()[0])
                calculated = [data.mean(), data.std(ddof=1),
                              np.fabs(data-data.mean()).mean()]
            else:
                calculated = getattr(ha, key)[0]

            msg = 'Mismatch between calculated and reference {}'
            assert_almost_equal(calculated, value,
                                decimal=4,
                                err_msg=msg.format(key))

    def test_multiple_selections(self, universe):
        ha = hel.HELANAL(universe, flatten_single_helix=True,
                         select=('name CA and resnum 30-40',
                                 'name CA and resnum 60-80'))
        ha.run()
        n_frames = len(universe.trajectory)
        assert len(ha.summary) == 2
        assert len(ha.all_bends) == 2
        assert ha.all_bends[0].shape == (n_frames, 8, 8)
        assert ha.all_bends[1].shape == (n_frames, 18, 18)

    def test_universe_from_origins(self, helanal):
        u = helanal.universe_from_origins()
        assert isinstance(u, mda.Universe)
        assert len(u.atoms) == len(helanal.atomgroups[0])-2
        assert len(u.trajectory) == 70

    def test_multiple_atoms_per_residues(self):
        u = mda.Universe(XYZ)
        with pytest.warns(UserWarning) as rec:
            ha = hel.HELANAL(u, select='name H')
        assert len(rec) == 1
        assert 'multiple atoms' in rec[0].message.args[0]


def test_pdot():
    arr = np.random.rand(4, 3)
    matrix_dot = hel.pdot(arr, arr)
    list_dot = [np.dot(a, a) for a in arr]
    assert_almost_equal(matrix_dot, list_dot)


def test_pnorm():
    arr = np.random.rand(4, 3)
    matrix_norm = hel.pnorm(arr)
    list_norm = [np.linalg.norm(a) for a in arr]
    assert_almost_equal(matrix_norm, list_norm)


def test_vector_of_best_fit():
    line = np.random.rand(3)
    unit = line / np.linalg.norm(line)
    points = line*np.arange(1000)[:, np.newaxis]
    noise = np.random.normal(size=(1000, 3))
    data = points + noise

    vector = hel.vector_of_best_fit(data)
    cos = np.dot(vector, unit)
    assert_almost_equal(abs(cos), 1.0, decimal=5)
