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
import numpy as np
from numpy.testing import assert_almost_equal, assert_equal

import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PSF, DCD, PDB, TPR, XTC
from MDAnalysis.coordinates.base import Timestep
from MDAnalysis.analysis import align, diffusionmap

@pytest.fixture()
def u1():
    # 3341 atoms, 98 frames
    return mda.Universe(PSF, DCD)

@pytest.fixture()
def u2():
    # 47681 atoms, 1 frame
    return mda.Universe(PDB)

@pytest.fixture()
def u3():
    # 47681 atoms, 10 frames
    return mda.Universe(TPR, XTC)

@pytest.fixture()
def ens(u1, u2, u3):
    return mda.Ensemble([u1, u2, u3])
class TestEnsemble(object):

    def test_analysis_multiple_same(self, u3):
        ens = mda.Ensemble([u3, u3]).select_atoms("name CA")
        dist_mat = diffusionmap.DistanceMatrix(ens).run()
        dm = dist_mat.dist_matrix
        assert dm[0, 0] == 0
        assert dm[0, 9] != 0
        assert dm[0, 19] != 0
        assert dm[19, 10] != 0
        assert dm[0, 10] == 0
        assert dm[10, 0] == 0
        assert dm[-1, -2] != 0
        assert dm[1, 2] == dm[1, 12]


    def test_analysis_multiple_same_in_memory(self, u3):
        ens = mda.Ensemble([u3, u3])
        ens.transfer_to_memory()
        dist_mat = diffusionmap.DistanceMatrix(ens).run()
        dm = dist_mat.dist_matrix
        assert dm[0, 0] == 0
        assert dm[0, 9] != 0
        assert dm[0, 19] != 0
        assert dm[19, 10] != 0
        assert dm[0, 10] == 0
        assert dm[10, 0] == 0
        assert dm[-1, -2] != 0
        assert dm[1, 2] == dm[1, 12]

    def test_analysis_multiple_in_memory(self, u1, u3):
        ens = mda.Ensemble([u1, u3])
        ens.transfer_to_memory()
        dist_mat = diffusionmap.DistanceMatrix(ens).run()
        dm = dist_mat.dist_matrix
        assert dm[-1, -2] != 0


    def test_from_universe(self, u1):
        ens = mda.Ensemble(u1)
        assert len(ens.universes) == ens.n_universes == 1
        assert len(ens) == ens.n_frames == len(u1.trajectory)
        assert ens.labels == [DCD]
        assert ens._universe_labels == {DCD: 0}
        assert ens.trajectory is ens
        assert ens._ts_u == 0
        assert ens.n_atoms == 3341
        assert ens._ags == [u1.atoms]
        assert_equal(ens._ag_indices, [u1.atoms.ix])
        assert_equal(ens.frames, np.arange(98))
        assert_equal(ens.traj_frames, (98,))
        assert_equal(ens._frame_edges, [0, 98])

    def test_from_universes(self, u1, u2, u3):
        ens = mda.Ensemble([u1, u2, u3])
        assert len(ens.universes) == ens.n_universes == 3
        assert len(ens) == ens.n_frames == 109
        assert ens.labels == [DCD, PDB, XTC]
        assert ens._universe_labels == {DCD: 0, PDB: 1, XTC: 2}
        assert ens.trajectory is ens
        assert ens._ts_u == 0
        assert ens.n_atoms == 3341
        assert ens._ags == [u1.atoms, u2.atoms[:3341], u3.atoms[:3341]]
        assert_equal(ens._ag_indices, [u1.atoms.ix]*3)
        assert_equal(ens.frames, np.arange(109))
        assert_equal(ens.traj_frames, (98, 1, 10))
        assert_equal(ens._frame_edges, [0, 98, 99, 109])

    def test_from_atomgroups(self, u1, u2, u3):
        prot = [u.select_atoms('name CA') for u in [u3, u1, u2]]
        ens = mda.Ensemble(prot)
        assert len(ens.universes) == 3
        assert ens.labels == [XTC, DCD, PDB]
        assert ens._universe_labels == {DCD: 1, PDB: 2, XTC: 0}
        assert ens.n_atoms == 214
        assert ens._ags == prot
        assert_equal(ens._frame_edges, [0, 10, 108, 109])
    
    @pytest.mark.parametrize('start,stop,step', [
        (0, 10, None),  # default
        (10, 109, None),  # slice within universe
        (97, 109, None),  # slice across universe
        (90, 109, 2)  # slice with step
    ])
    def test_from_ensemble_slice(self, ens, start, stop, step):
        ens2 = ens[start:stop:step]
        frames = np.arange(start=start, stop=stop, step=step)
        assert len(ens2) == len(frames), 'incorrect number of frames'
        assert_equal(ens2.frames, frames, err_msg='incorrect frame indices')
        assert ens.ts is ens2.ts
        assert ens2.ts.frame == 0
        assert ens2._ts_u == 0, '_ts_u not set to 0 on creation'

    @pytest.mark.parametrize('ens_frames,traj_frames', [
        ([0, 1, 2], [0, 1, 2]),
        ([97, 98, 99, 100], [97, 0, 0, 1]),
        ([98, 1, 3, -2, 1, 1], [0, 1, 3, 8, 1, 1])
    ])
    def test_from_frame_list(self, ens, ens_frames, traj_frames):
        ens2 = ens[ens_frames]
        assert len(ens2) == len(traj_frames)
        for ts, i in zip(ens2, traj_frames):
            assert ts.frame == i, 'frame mismatch'

    def test_from_bool_arr(self, ens):
        mask = np.zeros(109, dtype=bool)
        mask[95:105:2] = True
        ens2 = ens[mask]
        assert len(ens2) == 5
        assert_equal(ens2.frames, np.arange(95, 105, 2))

    @pytest.mark.parametrize('index,n_universe,frame', [
        (0, 0, 0),  # first
        (-109, 0, 0),  # first
        (-1, 2, 9),  # last
        (108, 2, 9),  # last
        (97, 0, 97),
        (98, 1, 0),
        (99, 2, 0),
        (100, 2, 1),
        (4, 0, 4),
        (102, 2, 3),
    ])
    def test_ensemble_index(self, ens, index, n_universe, frame):
        ts = ens[index]
        assert isinstance(ts, Timestep)
        assert ts is ens.ts
        assert ts.frame == frame
        assert ens._ts_u == n_universe
        
    @pytest.mark.parametrize('start,stop,step,n_universe,frame4', [
        (0, 10, None, 0, 3),  # default
        (10, 109, None, 0, 13),  # slice within universe
        (97, 109, None, 2, 1),  # slice across universe
        (92, 109, 2, 1, 0)  # slice with step
    ])
    def test_sliced_ensemble_index(self, ens, start, stop, step, n_universe,
                                   frame4):
        ens2 = ens[start:stop:step]
        ts4 = ens2[3]
        assert ens2._ts_u == n_universe, 'incorrect Universe selected'
        assert ens2.ts is ts4, 'incorrect frame selected'
        assert ens2.ts.frame == frame4, 'incorrect frame selected'
        assert ens._ts_u == 0, 'initial Ensemble has persisting changes'

    def test_select_atoms_creation(self, u1, u2):
        sel = 'name CA'
        ens = mda.Ensemble([u1, u2], select=sel)
        assert len(ens) == 99
        assert len(ens.atoms) == 214
        ix = [u1.select_atoms(sel).ix, u2.select_atoms(sel).ix]
        assert_equal(ens._ag_indices, ix)

    def test_select_atoms_from_ens(self, u2, u3):
        sel = 'name CA'
        ens = mda.Ensemble([u2, u3])
        assert len(ens) == 11
        assert len(ens.atoms) == 47681
        ensca = ens.select_atoms(sel)
        assert len(ensca.atoms) == 214
        ix = [u2.select_atoms(sel).ix, u3.select_atoms(sel).ix]
        assert_equal(ensca._ag_indices, ix)

    def test_get_universe_from_file_label(self, ens, u1, u2, u3):
        assert ens[DCD] is u1.trajectory.ts
        assert ens[PDB] is u2.trajectory.ts
        assert ens[XTC] is u3.trajectory.ts

    def test_get_universe_from_given_label(self, u1, u2, u3):
        ens = mda.Ensemble([u1, u2, u3], labels=list('ABC'))
        assert ens['A'] is u1.trajectory.ts
        assert ens['B'] is u2.trajectory.ts
        assert ens['C'] is u3.trajectory.ts

    @pytest.mark.parametrize('start,stop,step,frames,n_ags', [
        (None, None, None, None, [98, 1, 10]),
        (10, 105, 2, None, [44, 1, 5]),
        (10, 105, 2, [97, 98, 99, 100], [1, 1, 2]),
    ])
    def test_iter_over_ag(self, ens, start, stop, step, frames, n_ags):
        ags = [a for ag, n in zip(ens._ags, n_ags) for a in [ag]*n]
        for ag, iag in zip(ags, ens.iterate_over_atomgroups(start=start,
                                                            stop=stop,
                                                            step=step,
                                                            frames=frames)):
            assert_almost_equal(ag.positions, iag.positions)

