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
from __future__ import absolute_import, division
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
    'local_nres_per_turn summary': np.array([3.64864163,  0.152694,  0.1131402]),
    # changed this implementation
    # 'local_screw_angles':
    #     np.array([87.80540079, -171.86019984,  -75.2341296,   24.61695962,
    #               121.77104796, -134.94786976,  -35.07857424,   58.9621866,
    #               159.40210233, -104.83368122,   -7.54816243,   87.40202629,
    #               -176.13071955,  -89.13196878,   17.17321345,  105.26627814,
    #               -147.00075298,  -49.36850775,   54.24557615,  156.06486532,
    #               -110.82698327,   -5.72138626,   85.36050546, -167.28218858,
    #               -68.23076936]),
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


def test_helix_analysis_zigzag():
    #      x    x    x    x    x
    #     / \  / \  / \  / \  /
    #   x    x    x    x    x

    n_atoms = 100
    u = mda.Universe.empty(100, trajectory=True)
    xyz = np.array(list(zip([1, -1]*(n_atoms//2),  # x \in {0, 1}
                            [0]*n_atoms,  # y == 0
                            range(n_atoms))))  # z rises continuously
    u.load_new(xyz)
    properties = hel.helix_analysis(u.atoms.positions, ref_axis=[0, 0, 1])
    assert_almost_equal(properties['local_twists'], 180, decimal=4)
    assert_almost_equal(properties['local_nres_per_turn'], 2, decimal=4)
    assert_almost_equal(properties['global_axis'],
                        [0, 0, -1], decimal=4)
    # all 0 vectors
    assert_almost_equal(properties['local_axes'], 0, decimal=4)
    assert_almost_equal(properties['local_bends'], 0, decimal=4)
    assert_almost_equal(properties['all_bends'], 0, decimal=4)
    assert_almost_equal(properties['local_heights'], 0, decimal=4)
    assert_almost_equal(properties['local_helix_directions'][0::2] - [-1, 0, 0],
                        0, decimal=4)
    assert_almost_equal(properties['local_helix_directions'][1::2] - [1, 0, 0],
                        0, decimal=4)
    origins = xyz[1:-1]
    origins[:, 0] = 0
    assert_almost_equal(properties['local_origins'], origins, decimal=4)
    assert_almost_equal(properties['local_screw_angles'], [
                        180, 0]*49, decimal=4)


def test_helix_analysis_square_oct():
    # square-octagon-square-octagon
    n_rep = 10
    sq2 = 0.5 ** 0.5
    square = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    octagon = [(1, 0), (sq2, sq2), (0, 1), (-sq2, sq2),
               (-1, 0), (-sq2, -sq2), (0, -1), (sq2, -sq2)]
    shapes = (square+octagon)*n_rep
    xyz = np.array(list(zip(np.arange(len(shapes)),
                            *zip(*shapes))))
    n_atoms = len(xyz)
    u = mda.Universe.empty(n_atoms, trajectory=True)
    u.load_new(xyz)

    properties = hel.helix_analysis(u.atoms.positions, ref_axis=[0, 0, 1])
    twist_trans = [102.76438, 32.235607]
    twists = ([90]*2 + twist_trans + [45]*6 + twist_trans[::-1]) * n_rep
    assert_almost_equal(properties['local_twists'], twists[:n_atoms-3],
                        decimal=4)
    res_trans = [3.503159, 11.167775]
    res = ([4]*2 + res_trans + [8]*6 + res_trans[::-1]) * n_rep
    assert_almost_equal(properties['local_nres_per_turn'], res[:n_atoms-3],
                        decimal=4)
    assert_almost_equal(properties['global_axis'],
                        [-1, 0, 0], decimal=3)
    assert_almost_equal(properties['local_axes']-[1, 0, 0], 0, decimal=4)
    assert_almost_equal(properties['local_bends'], 0, decimal=4)
    assert_almost_equal(properties['all_bends'], 0, decimal=4)
    assert_almost_equal(properties['local_heights'], 1, decimal=4)

    loc_rot = [[0.,  0.,  1.],
               [0., -1.,  0.],
               [0.,  0., -1.],
               [0.,  0.97528684,  0.2209424],
               [0.,  0.70710677,  0.70710677],  # 0.5 ** 0.5
               [0.,  0.,  1.],
               [0., -0.70710677,  0.70710677],
               [0., -1.,  0.],
               [0., -0.70710677, -0.70710677],
               [0.,  0., -1.],
               [0.,  0.70710677, -0.70710677],
               [0.,  0.97528684, -0.2209424]] * n_rep
    assert_almost_equal(properties['local_helix_directions'],
                        loc_rot[:n_atoms-2], decimal=4)

    origins = xyz[1:-1]
    origins[:, 1:] = ([[0.,   0.],  # square
                       [0.,   0.],
                       [0.,  -0.33318555],  # square -> octagon
                       [-1.7878988,  -0.6315732],  # square -> octagon
                       [0.,   0.],  # octagon
                       [0.,   0.],
                       [0.,   0.],
                       [0.,   0.],
                       [0.,   0.],
                       [0.,   0.],
                       [-1.3141878,   1.3141878],  # octagon -> square
                       [0.34966463,   0.14732757]]*n_rep)[:len(origins)]
    assert_almost_equal(properties['local_origins'], origins,
                        decimal=4)

    # a bit off-90
    screw = [0, 89.99, 180, 77.2356, 45, 0,
             45, 89.99, 135, 180, 135, 102.7644] * n_rep
    assert_almost_equal(properties['local_screw_angles'], screw[:-2],
                        decimal=2)


class TestHELANAL(object):

    @pytest.fixture()
    def universe(self):
        return mda.Universe(PSF, DCD)

    @pytest.fixture()
    def helanal(self, universe):
        ha = hel.HELANAL(universe, select='name CA and resnum 161-187',
                         flatten_single_helix=True)
        return ha.run(start=10, stop=80)

    def test_regression_summary(self, helanal):
        bends = helanal.summary['all_bends']
        old_helanal = read_bending_matrix(HELANAL_BENDING_MATRIX_SUBSET)
        assert_almost_equal(np.triu(bends['mean'], 1), old_helanal['Mean'],
                            decimal=1)
        assert_almost_equal(np.triu(bends['sample_sd'], 1), old_helanal['SD'],
                            decimal=1)
        assert_almost_equal(np.triu(bends['abs_dev'], 1), old_helanal['ABDEV'],
                            decimal=1)

    def test_regression_values(self):
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


def test_vector_of_best_fit():
    line = np.random.rand(3)
    unit = line / np.linalg.norm(line)
    points = line*np.arange(1000)[:, np.newaxis]
    noise = np.random.normal(size=(1000, 3))
    data = points + noise

    vector = hel.vector_of_best_fit(data)
    cos = np.dot(vector, unit)
    assert_almost_equal(abs(cos), 1.0, decimal=5)
